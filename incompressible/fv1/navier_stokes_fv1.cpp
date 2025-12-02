/*
 * Copyright (c) 2010-2014:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "navier_stokes_fv1.h"

#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesFV1<TDomain>::NavierStokesFV1(const char* functions,
                                          const char* subsets)
: IncompressibleNavierStokesBase<TDomain>(functions, subsets)
{
    init();
};

template<typename TDomain>
NavierStokesFV1<TDomain>::NavierStokesFV1(const std::vector<std::string>& vFct,
                                          const std::vector<std::string>& vSubset)
: IncompressibleNavierStokesBase<TDomain>(vFct, vSubset)
{
    init();
};


template<typename TDomain>
void NavierStokesFV1<TDomain>::init()
{
    //	check number of functions
    if(this->num_fct() != dim+1)
        UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
                 " needs exactly "<<dim+1<<" symbolic function.");
    
    
    
    
    m_imDensitySCVF.set_comp_lin_defect(false);
    m_imDensitySCV.set_comp_lin_defect(false);
    
    m_imSourceSCVF.set_comp_lin_defect(false);
    m_imSourceSCV.set_comp_lin_defect(false);
    
    m_imKinViscosity.set_comp_lin_defect(false);
    m_imKinViscositySCV.set_comp_lin_defect(false);
    
    m_imVolumeFraction.set_comp_lin_defect(false);
    m_imSurfaceNormal.set_comp_lin_defect(false);
    
    
    //	register imports
    this->register_import(m_imSourceSCV);
    this->register_import(m_imSourceSCVF);
    this->register_import(m_imKinViscosity);
    this->register_import(m_imKinViscositySCV);
    this->register_import(m_imDensitySCVF);
    this->register_import(m_imDensitySCV);
    this->register_import(m_imSurfaceNormal);
    this->register_import(m_imVolumeFraction);
    
    
    m_imSourceSCV.set_rhs_part();
    m_imSourceSCVF.set_rhs_part();
    m_imDensitySCV.set_mass_part();
    
    
    
    //	default value for density
    base_type::set_density(1.0);
    
    //	update assemble functions
    register_all_funcs(false);
    
}

template<typename TDomain>
void NavierStokesFV1<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
    if(bNonRegularGrid)
        UG_THROW("NavierStokes: only regular grid implemented.");
    
    //	check number
    if(vLfeID.size() != dim+1)
        UG_THROW("NavierStokes: Need exactly "<<dim+1<<" functions");
    
    for(int d = 0; d <= dim; ++d)
        if(vLfeID[d].type() != LFEID::LAGRANGE || vLfeID[d].order() != 1)
            UG_THROW("NavierStokes: 'fv1' expects Lagrange P1 trial space "
                     "for velocity and pressure.");
    
    //	update assemble functions
    register_all_funcs(false);
}

template<typename TDomain>
void NavierStokesFV1<TDomain>::
set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data)
{
    m_imKinViscosity.set_data(data);
    m_imKinViscositySCV.set_data(data);
}

template<typename TDomain>
void NavierStokesFV1<TDomain>::
set_density(SmartPtr<CplUserData<number, dim> > data)
{
    m_imDensitySCVF.set_data(data);
    m_imDensitySCV.set_data(data);
}


template<typename TDomain>
void NavierStokesFV1<TDomain>::
set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
{
    m_imSourceSCV.set_data(data);
    m_imSourceSCVF.set_data(data);
}

template<typename TDomain>
void NavierStokesFV1<TDomain>::
set_interface_normal(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
{
    m_imSurfaceNormal.set_data(data);
}

template<typename TDomain>
void NavierStokesFV1<TDomain>::
set_vol_fraction(SmartPtr<CplUserData<number, dim> > data)
{
    m_imVolumeFraction.set_data(data);
}

template<typename TDomain>
void NavierStokesFV1<TDomain>::
set_relative_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
{
    m_imRelativeVelocity.set_data(data);
}

////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
    // 	Only first order implementation
    if(!(TFVGeom::order == 1))
        UG_THROW("Only first order implementation, but other Finite Volume"
                 " Geometry set.");
    
    //	check, that stabilization has been set
    if(m_spStab.invalid())
        UG_THROW("Stabilization has not been set.");
    
    //	init stabilization for element type
    m_spStab->template set_geometry_type<TFVGeom >();
    
    if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
    {
        //	check, that convective upwinding has been set
        if(m_spConvUpwind.invalid())
            UG_THROW("Upwinding for convective Term in Momentum eq. not set.");
        
        //	init convection stabilization for element type
        if(m_spConvStab.valid())
            m_spConvStab->template set_geometry_type<TFVGeom >();
        
        //	init convection stabilization for element type
        if(m_spConvUpwind.valid())
            m_spConvUpwind->template set_geometry_type<TFVGeom >();
    }
    
    //	check, that kinematic Viscosity has been set
    if(!m_imKinViscosity.data_given())
        UG_THROW("NavierStokes::prep_elem_loop:"
                 " Kinematic Viscosity has not been set, but is required.");
    
    //	check, that Density has been set
    if(!m_imDensitySCVF.data_given())
        UG_THROW("NavierStokes::prep_elem_loop:"
                 " Density has not been set, but is required.");
    
    //	check, that Density has been set
    if(!m_imDensitySCV.data_given())
        UG_THROW("NavierStokes::prep_elem_loop:"
                 " Density has not been set, but is required.");
    
    
    
    //	set local positions for imports
    if(!TFVGeom::usesHangingNodes)
    {
        static const int refDim = TElem::dim;
        TFVGeom& geo = GeomProvider<TFVGeom>::get();
        const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
        const size_t numSCVFip = geo.num_scvf_ips();
        const MathVector<refDim>* vSCVip = geo.scv_local_ips();
        const size_t numSCVip = geo.num_scv_ips();
        m_imKinViscosity.template set_local_ips<refDim>(vSCVFip,numSCVFip,true);
        m_imKinViscositySCV.template set_local_ips<refDim>(vSCVip,numSCVip,true);
        m_imDensitySCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip,true);
        m_imDensitySCV.template set_local_ips<refDim>(vSCVip,numSCVip,true);
        m_imSourceSCV.template set_local_ips<refDim>(vSCVip,numSCVip);
        m_imSourceSCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip);
        m_imSurfaceNormal.template set_local_ips<refDim>(vSCVip,numSCVip);
        m_imVolumeFraction.template set_local_ips<refDim>(vSCVip,numSCVip,true);
        
    }
    
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
    // 	Update Geometry for this element
    TFVGeom& geo = GeomProvider<TFVGeom>::get();
    try{
        geo.update(elem, vCornerCoords, &(this->subset_handler()));
    }
    UG_CATCH_THROW("NavierStokes::prep_elem:"
                   " Cannot update Finite Volume Geometry.");
    
    //	set local positions for imports
    if(TFVGeom::usesHangingNodes)
    {
        //	request ip series
        static const int refDim = TElem::dim;
        const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
        const size_t numSCVFip = geo.num_scvf_ips();
        const MathVector<refDim>* vSCVip = geo.scv_local_ips();
        const size_t numSCVip = geo.num_scv_ips();
        m_imKinViscosity.template set_local_ips<refDim>(vSCVFip,numSCVFip,true);
        m_imKinViscositySCV.template set_local_ips<refDim>(vSCVip,numSCVip,true);
        m_imDensitySCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip,true);
        m_imDensitySCV.template set_local_ips<refDim>(vSCVip,numSCVip,true);
        m_imSourceSCV.template set_local_ips<refDim>(vSCVip,numSCVip);
        m_imSourceSCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip);
        m_imSurfaceNormal.template set_local_ips<refDim>(vSCVip,numSCVip);
        m_imVolumeFraction.template set_local_ips<refDim>(vSCVip,numSCVip,true);
        
    }
    
    //	set global positions for imports
    const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
    const size_t numSCVFip = geo.num_scvf_ips();
    const MathVector<dim>* vSCVip = geo.scv_global_ips();
    const size_t numSCVip = geo.num_scv_ips();
    m_imKinViscosity.set_global_ips(vSCVFip, numSCVFip);
    m_imKinViscositySCV.set_global_ips(vSCVip, numSCVip);
    m_imDensitySCVF.set_global_ips(vSCVFip, numSCVFip);
    m_imDensitySCV.set_global_ips(vSCVip, numSCVip);
    m_imSourceSCV.set_global_ips(vSCVip, numSCVip);
    m_imSourceSCVF.set_global_ips(vSCVFip, numSCVFip);
    m_imSurfaceNormal.set_global_ips(vSCVip, numSCVip);
    m_imVolumeFraction.set_global_ips(vSCVip, numSCVip);
    
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
    static const size_t numSCVF = TFVGeom::numSCVF;
    //  number of shape functions
    static const size_t numSh = reference_element_traits<TElem>::reference_element_type::numCorners;
    
    
    // 	Only first order implementation
    UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");
    
    // 	get finite volume geometry
    static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
        
    
    //	check for solutions to pass to stabilization in time-dependent case
    const LocalVector *pSol = &u, *pOldSol = NULL;
    number dt = 0.0;
    if(this->is_time_dependent())
    {
        //	get and check current and old solution
        const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
        if(vLocSol->size() != 2)
            UG_THROW("NavierStokes::add_jac_A_elem: "
                     " Stabilization needs exactly two time points.");
        
        //	remember local solutions
        pSol = &vLocSol->solution(0);
        pOldSol = &vLocSol->solution(1);
        dt = vLocSol->time(0) - vLocSol->time(1);
    }
    
    //	interpolate velocity at ip with standard lagrange interpolation
    MathVector<dim> StdVel[numSCVF];
    MathVector<dim> Vel_ip[numSCVF];
    number Rho_up[numSCVF];
    number Rho_do[numSCVF];
   
    
    std_vel(  u,  geo,  StdVel, m_imDensitySCV,Rho_up,Rho_do);
    
    
	m_spStab->update(&geo, *pSol, StdVel, m_bStokes, m_imKinViscosity, m_imKinViscositySCV, m_imDensitySCVF, NULL, m_imDensitySCV, NULL, m_imSourceSCVF, m_imSourceSCV, pOldSol, dt);
    
    if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
    {
        //	compute stabilized velocities and shapes for convection upwind
        if(m_spConvStab.valid())
            if(m_spConvStab != m_spStab)
                m_spConvStab->update(&geo, *pSol, StdVel, false, m_imKinViscosity, m_imKinViscositySCV, m_imDensitySCVF, NULL, m_imDensitySCV, NULL, m_imSourceSCVF, m_imSourceSCV, pOldSol, dt);
    }
    //    get a const (!!) reference to the stabilization
    const INavierStokesFV1Stabilization<dim>& stab = *m_spStab;
    // Estimation of the veloctity at ip for upwind shape and continuity equations
    for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
        Vel_ip[ip] = stab.stab_vel(ip);
    
    if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
    {
	//	compute upwind shapes
		if(m_spConvUpwind.valid())
			if(m_spStab->upwind() != m_spConvUpwind)
            {
                m_spConvUpwind->update(&geo, StdVel);
                m_spConvUpwind->update_downwind(&geo, StdVel);
            }
	}


	const INavierStokesFV1Stabilization<dim>& convStab = *m_spConvStab;
	const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;
    

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	// 	loop shape functions
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
        {
            ////////////////////////////////////////////////////
            ////////////////////////////////////////////////////
            // Momentum Equation (conservation of momentum)
            ////////////////////////////////////////////////////
            ////////////////////////////////////////////////////
            
            ////////////////////////////////////////////////////
            // Diffusive Term (Momentum Equation)
            ////////////////////////////////////////////////////
            
            //     Compute flux derivative at IP
            const number flux_sh =  -1.0 * m_imKinViscosity[ip] * m_imDensitySCVF[ip] *
             VecDot(scvf.global_grad(sh), scvf.normal())   ;
            
            //     Add flux derivative  to local matrix
            for(int d1 = 0; d1 < dim; ++d1)
            {
                J(d1, scvf.from(), d1, sh) += flux_sh;
                J(d1, scvf.to()  , d1, sh) -= flux_sh;
            }
            
            if(!m_bLaplace)
            {
                for(int d1 = 0; d1 < dim; ++d1)
                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        const number flux2_sh = -1.0 * m_imKinViscosity[ip] * m_imDensitySCVF[ip] *
                        ( scvf.global_grad(sh)[d1] * scvf.normal()[d2] );
                        J(d1, scvf.from(), d2, sh) += flux2_sh;
                        J(d1, scvf.to()  , d2, sh) -= flux2_sh;
                    }
            }
            
            
            ////////////////////////////////////////////////////
            // Pressure Term (Momentum Equation)
            ////////////////////////////////////////////////////
	
			//	Add flux derivative for local matrix
			for(int d1 = 0; d1 < dim; ++d1)
			{
				const number flux_sh = scvf.shape(sh) * scvf.normal()[d1];
				J(d1, scvf.from(), _P_, sh) += flux_sh;
				J(d1, scvf.to()  , _P_, sh) -= flux_sh;
			}


            
            ////////////////////////////////////////////////////
            // Convective Term (Momentum Equation)
            ////////////////////////////////////////////////////
            
            // note: StdVel will be the advecting velocity (not upwinded)
            //       UpwindVel/StabVel will be the transported velocity
            
            if (! m_bStokes ) // no convective terms in the Stokes equation
            {
                //	compute upwind velocity
                MathVector<dim> UpwindMomentum;
                
                //	switch PAC
                if(m_spConvUpwind.valid())  UpwindMomentum = upwind.upwind_momentum(ip, u, m_imDensitySCV.values());
                else if (m_spConvStab.valid()) VecScale(UpwindMomentum, convStab.stab_vel(ip),m_imDensitySCVF[ip]);
                else UG_THROW("Cannot find upwind for convective term.");
                
                //	peclet blend
                number w = 1.0;
                if(m_bPecletBlend)
                    w = peclet_blend(UpwindMomentum, geo, ip, StdVel[ip], m_imKinViscosity[ip],m_imDensitySCVF[ip]);
                
                //	compute product of stabilized vel and normal
                const number prod = VecProd(StdVel[ip], scvf.normal());// * m_imDensitySCVF[ip];
                
                ///////////////////////////////////
                //	Add fixpoint linearization
                ///////////////////////////////////
                
                //	Stabilization used as upwind
                if(m_spConvStab.valid())
                {
                    //	velocity derivatives
                    if(stab.vel_comp_connected())
                        for(int d1 = 0; d1 < dim; ++d1)
                            for(int d2 = 0; d2 < dim; ++d2)
                            {
                                const number convFlux_vel = prod * w * convStab.stab_shape_vel(ip, d1, d2, sh) * m_imDensitySCVF[ip];
                                J(d1, scvf.from(), d2, sh) += convFlux_vel;
                                J(d1, scvf.to()  , d2, sh) -= convFlux_vel;
                            }
                    else
                        for(int d1 = 0; d1 < dim; ++d1)
                        {
                            const number convFlux_vel = prod * w * convStab.stab_shape_vel(ip, d1, d1, sh) * m_imDensitySCVF[ip];
                            J(d1, scvf.from(), d1, sh) += convFlux_vel;
                            J(d1, scvf.to()  , d1, sh) -= convFlux_vel;
                        }
                    
                    //	pressure derivative
                    for(int d1 = 0; d1 < dim; ++d1)
                    {
                        const number convFlux_p = prod * w * convStab.stab_shape_p(ip, d1, sh) * m_imDensitySCVF[ip];
                        
                        J(d1, scvf.from(), _P_, sh) += convFlux_p;
                        J(d1, scvf.to()  , _P_, sh) -= convFlux_p;
                    }
                }
                
                //	Upwind used as upwind
                if(m_spConvUpwind.valid())
                {
                    number convFlux_vel =  m_imDensitySCV[sh] * upwind.upwind_shape_sh(ip, sh);
                    
                    //	in some cases (e.g. PositiveUpwind, RegularUpwind) the upwind
                    //	velocity in an ip depends also on the upwind velocity in
                    //	other ips. This is reflected by the fact, that the ip
                    //	shapes are non-zero. In that case, we can interpolate an
                    //	approximate upwind only from the corner velocities by using
                    //	u_up = \sum shape_co U_co + \sum shape_ip \tilde{u}_ip
                    //	     = \sum shape_co U_co + \sum \sum shape_ip norm_shape_co|_ip * U_co
                    if(m_spConvUpwind->non_zero_shape_ip())
                    {
                        UG_THROW("Momentum trasport not implemented for u (Momentum upwind is missing).");
                        for(size_t ip2 = 0; ip2 < geo.num_scvf(); ++ip2)
                        {
                            const typename TFVGeom::SCVF& scvf2 = geo.scvf(ip2);
                            convFlux_vel += scvf2.shape(sh) * upwind.upwind_shape_ip(ip, ip2);
                        }
                    }
                    
                    convFlux_vel *= prod * w;
                    
                    for(int d1 = 0; d1 < dim; ++d1)
                    {
                        J(d1, scvf.from(), d1, sh) += convFlux_vel;
                        J(d1, scvf.to()  , d1, sh) -= convFlux_vel;
                    }
                }
                
                //	derivative due to peclet blending
                if(m_bPecletBlend)
                {
                    const number convFluxPe = prod * (1.0-w) * m_imDensitySCV[sh] * scvf.shape(sh);// * (m_imDensitySCVF[ip] / m_imDensitySCVF[ip]);
                    for(int d1 = 0; d1 < dim; ++d1)
                    {
                        J(d1, scvf.from(), d1, sh) += convFluxPe;
                        J(d1, scvf.to()  , d1, sh) -= convFluxPe;
                    }
                }
                
                /////////////////////////////////////////
                //	Add full jacobian (remaining part)
                /////////////////////////////////////////
                
                //	Add remaining term for exact jacobian
                if(m_boolFullNewton)
                {
                    // Jacobian for std Vel as convecting velocity
                    for(int d = 0; d < dim; ++d)
                    {
                        for(int d1 = 0; d1 < dim; ++d1)
                        {
                            /*number VelShape = (Rho_up[ip] * upwind.downwind_conv_length(ip) * upwind.upwind_shape_sh(ip, sh)  + Rho_do[ip] * upwind.upwind_conv_length(ip) * upwind.downwind_shape_sh(ip, sh) ) / ( Rho_up[ip] * upwind.downwind_conv_length(ip) + Rho_do[ip] * upwind.upwind_conv_length(ip));
                            number contFlux_vel =  m_bFullNewtonFactor * VelShape * scvf.normal()[d1];*/
                            number contFlux_vel =  m_bFullNewtonFactor * scvf.shape(sh) * scvf.normal()[d1];
                            J(d, scvf.from(), d1, sh) += contFlux_vel * UpwindMomentum[d];
                            J(d, scvf.to()  , d1, sh) -= contFlux_vel * UpwindMomentum[d];
                            
                        }
                        

                    }
                    
                } // end exact jacobian part
                
                
            } // end of if (! m_bStokes) for the convective terms
            
            ////////////////////////////////////////////////////
            ////////////////////////////////////////////////////
            // Continuity Equation (conservation of mass)
            ////////////////////////////////////////////////////
            ////////////////////////////////////////////////////
            
            //	Add derivative of stabilized flux w.r.t velocity comp to local matrix
            if(stab.vel_comp_connected())
            {
                for(int d1 = 0; d1 < dim; ++d1)
                {
                    number contFlux_vel = 0.0;
                    for(int d2 = 0; d2 < dim; ++d2)
                        contFlux_vel += stab.stab_shape_vel(ip, d2, d1, sh)
                        * scvf.normal()[d2]; //* m_imDensitySCVF[ip];
                    
                    J(_P_, scvf.from(), d1, sh) += contFlux_vel;
                    J(_P_, scvf.to()  , d1, sh) -= contFlux_vel;
                }
            }
            else
            {
                for(int d1 = 0; d1 < dim; ++d1)
                {
                    const number contFlux_vel = stab.stab_shape_vel(ip, d1, d1, sh)
                    * scvf.normal()[d1]; //* m_imDensitySCVF[ip];
                    
                    J(_P_, scvf.from(), d1, sh) += contFlux_vel;
                    J(_P_, scvf.to()  , d1, sh) -= contFlux_vel;
                }
            }
            
        
            //	Add derivative of stabilized flux w.r.t pressure to local matrix
            number contFlux_p = 0.0;
            for(int d1 = 0; d1 < dim; ++d1)
                contFlux_p += stab.stab_shape_p(ip, d1, sh) * scvf.normal()[d1]; //* m_imDensitySCVF[ip];
            
            J(_P_, scvf.from(), _P_, sh) += contFlux_p;
            J(_P_, scvf.to()  , _P_, sh) -= contFlux_p;
            
        }
	}
    
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	
// 	Only first order implemented
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    static const size_t numSCVF = TFVGeom::numSCVF;
    static const size_t numSh = reference_element_traits<TElem>::reference_element_type::numCorners;
   

//	check for solutions to pass to stabilization in time-dependent case
	const LocalVector *pSol = &u, *pOldSol = NULL;
	number dt = 0.0;
	if(this->is_time_dependent())
	{
	//	get and check current and old solution 
		const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
		if(vLocSol->size() != 2)
			UG_THROW("NavierStokes::add_def_A_elem: "
							" Stabilization needs exactly two time points.");

	//	remember local solutions
		pSol = &vLocSol->solution(0);
		pOldSol = &vLocSol->solution(1);
		dt = vLocSol->time(0) - vLocSol->time(1);
	}

//	interpolate velocity at ip with standard lagrange interpolation
    
	MathVector<dim> StdVel[numSCVF];
    MathVector<dim> Vel_ip[numSCVF];
    number Rho_up[numSCVF];
    number Rho_d[numSCVF];
   
    
    std_vel(  u,  geo,  StdVel, m_imDensitySCV,Rho_up,Rho_d);
    
    
//	compute stabilized velocities and shapes for continuity equation
	// \todo: (optional) Here we can skip the computation of shapes, implement?

	m_spStab->update(&geo, *pSol, StdVel, m_bStokes, m_imKinViscosity, m_imKinViscositySCV, m_imDensitySCVF, NULL, m_imDensitySCV, NULL, m_imSourceSCVF, m_imSourceSCV, pOldSol, dt);
    
    if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
    {
        //	compute stabilized velocities and shapes for convection upwind
        if(m_spConvStab.valid())
            if(m_spConvStab != m_spStab)
                m_spConvStab->update(&geo, *pSol, StdVel, false, m_imKinViscosity, m_imKinViscositySCV, m_imDensitySCVF, NULL, m_imDensitySCV, NULL, m_imSourceSCVF, m_imSourceSCV, pOldSol, dt);
    }
    
    //    get a const (!!) reference to the stabilization
    const INavierStokesFV1Stabilization<dim>& stab = *m_spStab;
    
    
    // Estimation of the veloctity at ip for upwind shape and continuity equations
    for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
        Vel_ip[ip] = stab.stab_vel(ip);

    
    if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
    {
	//	compute upwind shapes
		if(m_spConvUpwind.valid())
			if(m_spStab->upwind() != m_spConvUpwind)
            {
                m_spConvUpwind->update(&geo, StdVel);
                m_spConvUpwind->update_downwind(&geo, StdVel);
            }
	}


	const INavierStokesFV1Stabilization<dim>& convStab = *m_spConvStab;
	const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;


// 	loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Momentum Equation (conservation of momentum)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////

	//	\todo: we could also add all fluxes at once in order to save several
	//			accesses to the local defect, implement and loose clarity?

		////////////////////////////////////////////////////
		// Diffusive Term (Momentum Equation)
		////////////////////////////////////////////////////

	// 	1. Interpolate Functional Matrix of velocity at ip
		MathMatrix<dim, dim> grad;
		for(int d1 = 0; d1 < dim; ++d1)
			for(int d2 = 0; d2 <dim; ++d2)
			{
			//	sum up contributions of each shape
                grad(d1, d2) = 0.0;
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                    grad(d1, d2) += scvf.global_grad(sh)[d2] * u(d1, sh);
			}

	//	2. Compute flux
		MathVector<dim> diffFlux;

	//	Add (\nabla u) \cdot \vec{n}
		MatVecMult(diffFlux, grad, scvf.normal());

	//	Add (\nabla u)^T \cdot \vec{n}
		if(!m_bLaplace)
			TransposedMatVecMultAdd(diffFlux, grad, scvf.normal());

	//	scale by viscosity
		VecScale(diffFlux, diffFlux, (-1.0) * m_imDensitySCVF[ip] * m_imKinViscosity[ip]);

	//	3. Add flux to local defect
		for(int d1 = 0; d1 < dim; ++d1)
		{
			d(d1, scvf.from()) += diffFlux[d1];
			d(d1, scvf.to()  ) -= diffFlux[d1];
		}
        /*
    //     1. Interpolate Functional Matrix of velocity at ip
        MathMatrix<dim, dim> gradMomentum;
        for(int d1 = 0; d1 < dim; ++d1)
            for(int d2 = 0; d2 <dim; ++d2)
            {
            //    sum up contributions of each shape
                gradMomentum(d1, d2) = 0.0;
                for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                    gradMomentum(d1, d2) += scvf.global_grad(sh)[d2] * m_imDensitySCV[sh] * u(d1, sh);
                gradMomentum(d1, d2) -= StdVel[ip][d1] * RhoGrad[ip][d2];
            }

    //    2. Compute flux
        MathVector<dim> diffFlux;

    //    Add (\nabla u) \cdot \vec{n}
        MatVecMult(diffFlux, gradMomentum, scvf.normal());

    //    Add (\nabla u)^T \cdot \vec{n}
        if(!m_bLaplace)
            TransposedMatVecMultAdd(diffFlux, gradMomentum, scvf.normal());

    //    scale by viscosity
        VecScale(diffFlux, diffFlux, (-1.0) * m_imKinViscosity[ip]);

    //    3. Add flux to local defect
        for(int d1 = 0; d1 < dim; ++d1)
        {
            d(d1, scvf.from()) += diffFlux[d1];
            d(d1, scvf.to()  ) -= diffFlux[d1];
        }*/
		////////////////////////////////////////////////////
		// Convective Term (Momentum Equation)
		////////////////////////////////////////////////////

		// note: StdVel will be the advecting velocity (not upwinded)
		//       UpwindVel/StabVel will be the transported velocity

		if (! m_bStokes ) // no convective terms in the Stokes equation
		{
		//	find the upwind velocity at ip
			MathVector<dim> UpwindMomentum;
	
		//	switch PAC
            if(m_spConvUpwind.valid())  UpwindMomentum = upwind.upwind_momentum(ip, u, m_imDensitySCV.values());
            else if (m_spConvStab.valid())  VecScale(UpwindMomentum, convStab.stab_vel(ip),m_imDensitySCVF[ip]);
			else UG_THROW("Cannot find upwind for convective term.");
	
		//	Peclet Blend
			if(m_bPecletBlend)
                peclet_blend(UpwindMomentum, geo, ip, StdVel[ip], m_imKinViscosity[ip],m_imDensitySCVF[ip]);
            
            
	
		//	compute product of standard velocity and normal
            const number prod = VecProd(StdVel[ip], scvf.normal()) ;
	
		//	Add contributions to local velocity components
			for(int d1 = 0; d1 < dim; ++d1)
			{
				d(d1, scvf.from()) += UpwindMomentum[d1] * prod;
				d(d1, scvf.to()  ) -= UpwindMomentum[d1] * prod;
			}
            
            
            
		}

		////////////////////////////////////////////////////
		// Pressure Term (Momentum Equation)
		////////////////////////////////////////////////////
		//	1. Interpolate pressure at ip
		number pressure = 0.0;
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			pressure += scvf.shape(sh) * u(_P_, sh);
		
		//	2. Add contributions to local defect
		for(int d1 = 0; d1 < dim; ++d1)
		{
			d(d1, scvf.from()) += pressure * scvf.normal()[d1];
			d(d1, scvf.to()  ) -= pressure * scvf.normal()[d1];
		}
        

		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Continuity Equation (conservation of mass)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////

	//	compute flux at ip
        const number contFlux = VecProd(Vel_ip[ip], scvf.normal());// * m_imDensitySCVF[ip];

	//	Add contributions to local defect
		d(_P_, scvf.from()) += contFlux;
		d(_P_, scvf.to()  ) -= contFlux;
        
	}
    
    
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    static const size_t numSCV = TFVGeom::numSCV;
    
    /*number Vol = 0.0;
    number fac = 1.0;
    
    for(size_t ip = 0; ip < geo.num_scv(); ++ip)
    {
    //     get current SCV
        const typename TFVGeom::SCV& scv = geo.scv(ip);
        Vol += scv.volume();
    }

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int sh = scv.node_id();

	// 	loop velocity components
		for(int d1 = 0; d1 < dim; ++d1)
		{
		// 	Add to local matrix
            //J(d1, sh, d1, sh) += scv.volume() * m_imDensitySCV[ip];//(0.5 * m_imDensitySCV[ip] + 0.5 * Rho);
            J(d1, sh, d1, sh) += scv.volume() * (fac * m_imDensitySCV[ip]);//printf("Change the lin_def");
            //J(d1, sh, d1, sh) += scv.volume() * (0.0 * m_imDensitySCV[ip] + 1.0 * Rho);
            for(size_t ip2 = 0; ip2 < geo.num_scv(); ++ip2)
            {
                const typename TFVGeom::SCV& scv2 = geo.scv(ip2);
                const int sh2 = scv2.node_id();
                J(d1, sh, d1, sh2) += (1.0-fac) * scv.volume() * (m_imDensitySCV[ip] * scv2.volume()/Vol);
                
            }
		}
    }*/
    /*for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
    {
        //     get current SCV
        const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
        
    //    get associated SubControlVolumes
        const typename TFVGeom::SCV& scvFrom = geo.scv(scvf.from());
        const typename TFVGeom::SCV& scvTo   = geo.scv(scvf.to());
        
        const int from = scvFrom.node_id();
        const int to = scvTo.node_id();
        
        const number RhoDiff = (scvFrom.volume()*m_imDensitySCV[from]+scvTo.volume()*m_imDensitySCV[to]) / (scvFrom.volume()+scvTo.volume()) ;
        
        
        //     Add to local rhs
        for(int d1 = 0; d1 < dim; ++d1){
            J(d1, scvf.from(), d1, scvf.from()) +=  RhoDiff * scvFrom.volume() / 2.0;
            J(d1, scvf.to()  , d1, scvf.to()  ) +=  RhoDiff * scvTo.volume() / 2.0;
        
            //d(d1, scvf.from()) +=  RhoDiff * VecProd(m_imSourceSCVF[ip],scvf.normal()) * scvFrom.volume() * scvf.normal()[d1]/ (4.0*VecProd(scvf.normal(),scvf.normal()));
            //d(d1, scvf.to()  ) +=  RhoDiff * VecProd(m_imSourceSCVF[ip],scvf.normal()) * scvTo.volume() * scvf.normal()[d1]/ (4.0*VecProd(scvf.normal(),scvf.normal()));
        }
    }*/
    number RhoSCV[numSCV];
    number Count[numSCV];
    
    for(size_t ip = 0; ip < geo.num_scv(); ++ip)
        Count[ip]=0.0;
    
    for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
    {
        //     get current SCV
        const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
        
    //    get associated SubControlVolumes
        const typename TFVGeom::SCV& scvFrom = geo.scv(scvf.from());
        const typename TFVGeom::SCV& scvTo   = geo.scv(scvf.to());
        
        const int from = scvFrom.node_id();
        const int to = scvTo.node_id();

        RhoSCV[from] += m_imDensitySCVF[ip];
        RhoSCV[to] += m_imDensitySCVF[ip];
        
        Count[from] += 1.0;
        Count[to] += 1.0;
        
    }

    for(size_t ip = 0; ip < geo.num_scv(); ++ip)
    {
        
    //     get current SCV
        const typename TFVGeom::SCV& scv = geo.scv(ip);

    //     get associated node
        const int sh = scv.node_id();
        
        
        //     Add to local rhs
        for(int d1 = 0; d1 < dim; ++d1){
            J(d1, sh, d1, sh) += scv.volume() * RhoSCV[sh]/ Count[sh];//printf("Change the lin_def");
            
            

        }
    }
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    static const size_t numSCV = TFVGeom::numSCV;
    

    /*MathVector<dim> VRho; VecSet(VRho,0.0);
    number Vol = 0.0;
    number fac = 1.0;
    for(size_t ip = 0; ip < geo.num_scv(); ++ip)
    {
    //     get current SCV
        const typename TFVGeom::SCV& scv = geo.scv(ip);
        const int sh = scv.node_id();
        Vol += scv.volume();
        for(int d1 = 0; d1 < dim; ++d1)
            VRho[d1] += u(d1, sh) * m_imDensitySCV[ip] * scv.volume();
        
    }
    VRho *= 1.0/Vol;

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int sh = scv.node_id();

	// 	loop velocity components
		for(int d1 = 0; d1 < dim; ++d1)
		{
		// 	Add to local matrix
            //d(d1, sh) += u(d1, sh) * scv.volume() * m_imDensitySCV[ip];//(0.5 * m_imDensitySCV[ip] + 0.5 * Rho);
            d(d1, sh) +=  scv.volume() * (fac * u(d1, sh) * m_imDensitySCV[ip] + (1-fac) * VRho[d1]); //printf("Change the lin_def");
            //d(d1, sh) += u(d1, sh) * scv.volume() * (0.0 * m_imDensitySCV[ip] + 1.0 * Rho); //printf("Change the lin_def");
		}
        
    }*/
    /*for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
    {
        //     get current SCV
        const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
        
    //    get associated SubControlVolumes
        const typename TFVGeom::SCV& scvFrom = geo.scv(scvf.from());
        const typename TFVGeom::SCV& scvTo   = geo.scv(scvf.to());
        
        const int from = scvFrom.node_id();
        const int to = scvTo.node_id();
        
        const number RhoDiff = (scvFrom.volume()*m_imDensitySCV[from]+scvTo.volume()*m_imDensitySCV[to]) / (scvFrom.volume()+scvTo.volume()) ;
        
        
        //     Add to local rhs
        for(int d1 = 0; d1 < dim; ++d1){
            d(d1, scvf.from()) +=  u(d1, scvf.from())*RhoDiff * scvFrom.volume() / 2.0;
            d(d1, scvf.to()  ) +=  u(d1, scvf.to())*RhoDiff * scvTo.volume() / 2.0;
            
            //d(d1, scvf.from()) +=  RhoDiff * VecProd(m_imSourceSCVF[ip],scvf.normal()) * scvFrom.volume() * scvf.normal()[d1]/ (4.0*VecProd(scvf.normal(),scvf.normal()));
            //d(d1, scvf.to()  ) +=  RhoDiff * VecProd(m_imSourceSCVF[ip],scvf.normal()) * scvTo.volume() * scvf.normal()[d1]/ (4.0*VecProd(scvf.normal(),scvf.normal()));
        }
    }*/
    
    number RhoSCV[numSCV];
    number Count[numSCV];
    
    for(size_t ip = 0; ip < geo.num_scv(); ++ip)
        Count[ip]=0.0;
    
    for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
    {
        //     get current SCV
        const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
        
    //    get associated SubControlVolumes
        const typename TFVGeom::SCV& scvFrom = geo.scv(scvf.from());
        const typename TFVGeom::SCV& scvTo   = geo.scv(scvf.to());
        
        const int from = scvFrom.node_id();
        const int to = scvTo.node_id();

        RhoSCV[from] += m_imDensitySCVF[ip];
        RhoSCV[to] += m_imDensitySCVF[ip];
        
        Count[from] += 1.0;
        Count[to] += 1.0;
        
    }

    for(size_t ip = 0; ip < geo.num_scv(); ++ip)
    {
        
    //     get current SCV
        const typename TFVGeom::SCV& scv = geo.scv(ip);

    //     get associated node
        const int sh = scv.node_id();
        
        
        //     Add to local rhs
        for(int d1 = 0; d1 < dim; ++d1){
            d(d1, sh) +=  u(d1, sh) * RhoSCV[sh] * scv.volume()/ Count[sh];
            
            

        }
    }
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

//	if zero data given, return
    if(!m_imSourceSCV.data_given()) return;

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    /*if(m_imSourceSCV.data_given())
    {
        //     loop Sub Control Volumes (SCV)

        
        for(size_t ip = 0; ip < geo.num_scv(); ++ip)
        {
            //     get current SCV
            const typename TFVGeom::SCV& scv = geo.scv(ip);
            const int sh = scv.node_id();
            
            //     Add to local rhs
            for(int d1 = 0; d1 < dim; ++d1){
                d(d1, sh) += m_imSourceSCV[ip][d1]*scv.volume();

            }
        }

    }*/
    if(m_imSourceSCVF.data_given() )
    {

        for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
        {
            //     get current SCV
            const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
            
            //    get associated SubControlVolumes
            const typename TFVGeom::SCV& scvFrom = geo.scv(scvf.from());
            const typename TFVGeom::SCV& scvTo   = geo.scv(scvf.to());
            
            //const int from = scvFrom.node_id();
            //const int to = scvTo.node_id();
            
            if(m_imSourceSCVF.data_given())
            {
                //     Add to local rhs
                for(int d1 = 0; d1 < dim; ++d1)
                {
                    d(d1, scvf.from()) += m_imSourceSCVF[ip][d1] * scvFrom.volume() / 2.0;
                    d(d1, scvf.to()  ) += m_imSourceSCVF[ip][d1] * scvTo.volume() / 2.0;
                    
                    //d(d1, scvf.from()) +=  RhoDiff * VecProd(m_imSourceSCVF[ip],scvf.normal()) * scvFrom.volume() * scvf.normal()[d1]/ (4.0*VecProd(scvf.normal(),scvf.normal()));
                    //d(d1, scvf.to()  ) +=  RhoDiff * VecProd(m_imSourceSCVF[ip],scvf.normal()) * scvTo.volume() * scvf.normal()[d1]/ (4.0*VecProd(scvf.normal(),scvf.normal()));
                }
            }//end m_imSourceSCVF
        }// end of loop ips
    

    }
    /*if(m_imSourceSCVF.data_given())
    {
        //     loop Sub Control Volumes (SCV)

        MathVector<dim> RhoGrad;
        VecSet(RhoGrad, 0.0);
        for(size_t ip = 0; ip < geo.num_scv(); ++ip)
        {
            const typename TFVGeom::SCV& scv = geo.scv(ip);
            const int sh = scv.node_id();
            
            VecScaleAppend(RhoGrad, m_imDensitySCV[ip], scv.global_grad(sh));
        }
        

        for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
        {
            //     get current SCV
            const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
            
        //    get associated SubControlVolumes
            const typename TFVGeom::SCV& scvFrom = geo.scv(scvf.from());
            const typename TFVGeom::SCV& scvTo   = geo.scv(scvf.to());
            
            
            //     Add to local rhs
            for(int d1 = 0; d1 < dim; ++d1){
                d(d1, scvf.from()) -= 0.5*m_imSourceSCVF[ip][d1]*RhoGrad[dim-1]*geo.scvf_global_ips()[ip][dim-1]*scvFrom.volume();
                d(d1, scvf.to()  ) -= 0.5*m_imSourceSCVF[ip][d1]*RhoGrad[dim-1]*geo.scvf_global_ips()[ip][dim-1]*scvTo.volume();
            }
        }

    }*/
    /*if(m_imSourceSCVF.data_given())
    {
        //     loop Sub Control Volumes (SCV)

        MathVector<dim> RhoGrad;
        VecSet(RhoGrad, 0.0);
        for(size_t ip = 0; ip < geo.num_scv(); ++ip)
        {
            const typename TFVGeom::SCV& scv = geo.scv(ip);
            const int sh = scv.node_id();
            
            VecScaleAppend(RhoGrad, m_imDensitySCV[ip], scv.global_grad(sh));
        }
        

        for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
        {
            //     get current SCV
            const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
            
        //    get associated SubControlVolumes
            //const typename TFVGeom::SCV& scvFrom = geo.scv(scvf.from());
            //const typename TFVGeom::SCV& scvTo   = geo.scv(scvf.to());
            
            
            //     Add to local rhs
            for(int d1 = 0; d1 < dim; ++d1){
                d(d1, scvf.from()) += -m_imDensitySCVF[ip] * m_imSourceSCVF[ip][dim-1] * geo.scvf_global_ips()[ip][dim-1] * scvf.normal()[d1];
                d(d1, scvf.to()  ) += -m_imDensitySCVF[ip] * m_imSourceSCVF[ip][dim-1] * geo.scvf_global_ips()[ip][dim-1] * scvf.normal()[d1];
            }
        }

    }*/
    /*if(m_imSourceSCVF.data_given())
    {
        //     loop Sub Control Volumes (SCV)
        
        for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
        {
            //     get current SCV
            const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
            
        //    get associated SubControlVolumes
            const typename TFVGeom::SCV& scvFrom = geo.scv(scvf.from());
            const typename TFVGeom::SCV& scvTo   = geo.scv(scvf.to());
            
            const int from = scvFrom.node_id();
            const int to = scvTo.node_id();
            
            const number RhoDiff = (scvFrom.volume()*m_imDensitySCV[from]+scvTo.volume()*m_imDensitySCV[to]) / (scvFrom.volume()+scvTo.volume()) ;
            
            
            //     Add to local rhs
            for(int d1 = 0; d1 < dim; ++d1){
                d(d1, scvf.from()) +=  RhoDiff * m_imSourceSCVF[ip][d1] * scvFrom.volume() / 2.0;
                d(d1, scvf.to()  ) +=  RhoDiff * m_imSourceSCVF[ip][d1] * scvTo.volume() / 2.0;
                
                //d(d1, scvf.from()) +=  RhoDiff * VecProd(m_imSourceSCVF[ip],scvf.normal()) * scvFrom.volume() * scvf.normal()[d1]/ (4.0*VecProd(scvf.normal(),scvf.normal()));
                //d(d1, scvf.to()  ) +=  RhoDiff * VecProd(m_imSourceSCVF[ip],scvf.normal()) * scvTo.volume() * scvf.normal()[d1]/ (4.0*VecProd(scvf.normal(),scvf.normal()));
            }
        }

    }*/
    
}

template<typename TDomain>
template<typename TFVGeom>
inline
number
NavierStokesFV1<TDomain>::
peclet_blend(MathVector<dim>& UpwindMomentum, const TFVGeom& geo, size_t ip, 
             const MathVector<dim>& StdVel, number kinVisco, number densitySCVF)
{
	const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
//	compute peclet number
	number Pe = VecProd(StdVel, scvf.normal())/VecLength(scvf.normal())
	 * VecDistance(geo.corners() [scvf.to()], geo.corners() [scvf.from()]) / kinVisco;

//	compute weight
	const number Pe2 = Pe * Pe;
	const number w = Pe2 / (5.0 + Pe2);

//	compute upwind vel
	VecScaleAdd(UpwindMomentum, w, UpwindMomentum, (1.0-w) * densitySCVF, StdVel);

	return w;
}
template<typename TDomain>
template<typename TFVGeom>
inline
void
NavierStokesFV1<TDomain>::
std_vel( const LocalVector& u, const TFVGeom& geo, MathVector<dim>* StdVel, const DataImport<number, dim>& densitySCV, number* Rho_up, number* Rho_do)
{
    /*for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
    {
        const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
        
        VecSet(StdVel[ip], 0.0);
        for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
        {
            for(int d1 = 0; d1 < dim; ++d1)
            {
                StdVel[ip][d1] += scvf.shape(sh) * u(d1, sh) ;
            }
        }
    }*/
    for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
    {
        const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
        size_t to = scvf.to();
        size_t from = scvf.from();
        
        VecSet(StdVel[ip], 0.0);

        for(int d1 = 0; d1 < dim; ++d1)
        {
            StdVel[ip][d1] += (densitySCV[to] * u(d1, to) + densitySCV[from] * u(d1, from)) / ( densitySCV[to] + densitySCV[from]) ;
        }
    }
    if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
    {
    //    compute upwind shapes
        if(m_spConvUpwind.valid())
            {
                m_spConvUpwind->update(&geo, StdVel);
                m_spConvUpwind->update_downwind(&geo, StdVel);
                
                const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;

                for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
                {
                    const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
                    Rho_up[ip] = 0.0;
                    Rho_do[ip] = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                    {
                        Rho_up[ip] += upwind.upwind_shape_sh(ip, sh) * densitySCV[sh];
                        Rho_do[ip] += upwind.downwind_shape_sh(ip, sh) * densitySCV[sh];
                    }
                    VecSet(StdVel[ip], 0.0);
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        for(int d1 = 0; d1 < dim; ++d1)
                        {
                            StdVel[ip][d1] +=  (Rho_up[ip] * upwind.downwind_conv_length(ip) * upwind.upwind_shape_sh(ip, sh)  + Rho_do[ip] * upwind.upwind_conv_length(ip) * upwind.downwind_shape_sh(ip, sh) ) * u(d1, sh)/ ( Rho_up[ip] * upwind.downwind_conv_length(ip) + Rho_do[ip] * upwind.upwind_conv_length(ip));
                        }
                }
            }
    }
    
}

//    computes the linearized defect w.r.t to the density SCV
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
lin_def_densitySCV(const LocalVector& u,
                     std::vector<std::vector<number> > vvvLinDef[],
                     const size_t nip)
{
    //    request geometry
    const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    //    check for source term to pass to the stabilization
    
    //    reset the values for the linearized defect
    for(size_t ip = 0; ip < nip; ++ip)
        for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
            for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
                vvvLinDef[ip][c][sh] = 0.0;
    number Vol = 0.0;
    if(this->is_time_dependent())
        for(size_t ip = 0; ip < geo.num_scv(); ++ip)
        {
            const typename TFVGeom::SCV& scv = geo.scv(ip);
            Vol += scv.volume();
            
        }
    for(size_t ip = 0; ip < geo.num_scv(); ++ip)
    {
        
        //     get current SCV
        const typename TFVGeom::SCV& scv = geo.scv(ip);
        
        //     get associated node
        const int co = scv.node_id();
        
        //     set lin defect
        if(this->is_time_dependent())
        {
            for(size_t c = 0; c < dim; ++c)
            {
                vvvLinDef[ip][c][co] += u(c, co) * scv.volume() * 0.5;
                for(size_t ip = 0; ip < geo.num_scv(); ++ip)
                {
                    const int co2 = scv.node_id();
                    vvvLinDef[ip][c][co2] += u(c, co) *  0.5 * scv.volume()/Vol;
                }
            }
        }
        
        /*if(m_imSourceSCV.data_given())
        {
            for(size_t c = 0; c < dim; ++c)
            {
                vvvLinDef[ip][c][co] += m_imSourceSCV[ip][c] * scv.volume();
            }
        }*/
    }
    
    
    //     loop Sub Control Volume Faces (SCVF)
    /*if (! m_bStokes) // no convective terms in the Stokes equation
    {
        
    //    interpolate velocity at ip with standard lagrange interpolation
        
        MathVector<dim> StdVel[numSCVF];
        std_vel(  u,  geo,  StdVel, m_imDensitySCV);
         


        //    compute upwind shapes
        if(m_spConvUpwind.valid())
            if(m_spStab->upwind() != m_spConvUpwind)
                m_spConvUpwind->update(&geo, StdVel);

        
        //const INavierStokesFV1Stabilization<dim>& convStab = *m_spConvStab;
        const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;
        for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
        {
            //     get current SCVF
            const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
            for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
            {
                
                ////////////////////////////////////////////////////
                ////////////////////////////////////////////////////
                // Momentum Equation (conservation of momentum)
                ////////////////////////////////////////////////////
                ////////////////////////////////////////////////////
                
                ////////////////////////////////////////////////////
                // Convective Term (Momentum Equation)
                ////////////////////////////////////////////////////
                
                number convFlux_vel;
                //    Add contributions to local velocity components
                for(int d1 = 0; d1 < dim; ++d1)
                {
                    //    find the upwind velocity at ip
                    convFlux_vel = u(d1,sh) * upwind.upwind_shape_sh(ip, sh);
                    //    switch PAC
                    
                    
                    
                    //    compute product of standard velocity and normal
                    const number prod = VecProd(StdVel[ip], scvf.normal()) ;
                

                    
                    vvvLinDef[ip][d1][scvf.from()] += convFlux_vel * prod;
                    vvvLinDef[ip][d1][scvf.to()  ] -= convFlux_vel * prod;
                }
            }
            
        }
    }*/
    
}
//    computes the linearized defect w.r.t to the density SCVF
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
lin_def_densitySCVF(const LocalVector& u,
                     std::vector<std::vector<number> > vvvLinDef[],
                     const size_t nip)
{
//    request geometry
    const TFVGeom& geo = GeomProvider<TFVGeom>::get();

    static const size_t numSCVF = TFVGeom::numSCVF;
    //  number of shape functions
    static const size_t numSh = reference_element_traits<TElem>::reference_element_type::numCorners;
    
    //    reset the values for the linearized defect
    for(size_t ip = 0; ip < nip; ++ip)
        for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
            for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
                vvvLinDef[ip][c][sh] = 0.0;

    const LocalVector *pSol = &u, *pOldSol = NULL;
    number dt = 0.0;
    if(this->is_time_dependent())
    {
    //    get and check current and old solution
        const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
        if(vLocSol->size() != 2)
            UG_THROW("NavierStokes::add_def_A_elem: "
                            " Stabilization needs exactly two time points.");

    //    remember local solutions
        pSol = &vLocSol->solution(0);
        pOldSol = &vLocSol->solution(1);
        dt = vLocSol->time(0) - vLocSol->time(1);
    }
    
    MathVector<dim> StdVel[numSCVF];
    number Rho_up[numSCVF];
    number Rho_d[numSCVF];
   
    
    std_vel(  u,  geo,  StdVel, m_imDensitySCV,Rho_up,Rho_d);
    
    
    bool Phase2[numSCVF];
    number mu_l, rho_l;
    number mu_g, rho_g;
    MathVector<dim> Source_1;
    MathVector<dim> Source_g;
    

//     loop Sub Control Volume Faces (SCVF)
    for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
    {
    //     get current SCVF
        const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

        ////////////////////////////////////////////////////
        ////////////////////////////////////////////////////
        // Momentum Equation (conservation of momentum)
        ////////////////////////////////////////////////////
        ////////////////////////////////////////////////////


        ////////////////////////////////////////////////////
        // Diffusive Term (Momentum Equation)
        ////////////////////////////////////////////////////

    //     1. Interpolate Functional Matrix of velocity at ip
        MathMatrix<dim, dim> gradVel;
        for(int d1 = 0; d1 < dim; ++d1)
            for(int d2 = 0; d2 <dim; ++d2)
            {
            //    sum up contributions of each shape
                gradVel(d1, d2) = 0.0;
                for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                {
                    gradVel(d1, d2) += scvf.global_grad(sh)[d2] * u(d1, sh);
                }
            }

    //    2. Compute flux
        MathVector<dim> diffFlux;

    //    Add (\nabla u) \cdot \vec{n}
        MatVecMult(diffFlux, gradVel, scvf.normal());

    //    Add (\nabla u)^T \cdot \vec{n}
        if(!m_bLaplace)
            TransposedMatVecMultAdd(diffFlux, gradVel, scvf.normal());

    //    scale by viscosity
        VecScale(diffFlux, diffFlux, (-1.0) * m_imKinViscosity[ip]);

    //    3. Add flux to local defect
        for(int c = 0; c < dim; ++c)
        {
            vvvLinDef[ip][c][scvf.from()] += diffFlux[c];
            vvvLinDef[ip][c][scvf.to()  ] -= diffFlux[c];
        }
            

    }
        
}

//    computes the linearized defect w.r.t to the viscosity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
lin_def_viscosity(const LocalVector& u,
                     std::vector<std::vector<number> > vvvLinDef[],
                     const size_t nip)
{
//    request geometry
    const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    
    static const size_t numSCVF = TFVGeom::numSCVF;
    //  number of shape functions
    static const size_t numSh = reference_element_traits<TElem>::reference_element_type::numCorners;

    
    //    reset the values for the linearized defect
    for(size_t ip = 0; ip < nip; ++ip)
        for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
            for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
                vvvLinDef[ip][c][sh] = 0.0;
    //UG_LOG("Anfang add_def_A_elem");
    
//     Only first order implemented
    UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

    const LocalVector *pSol = &u, *pOldSol = NULL;
    number dt = 0.0;
    if(this->is_time_dependent())
    {
    //    get and check current and old solution
        const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
        if(vLocSol->size() != 2)
            UG_THROW("NavierStokes::add_def_A_elem: "
                            " Stabilization needs exactly two time points.");

    //    remember local solutions
        pSol = &vLocSol->solution(0);
        pOldSol = &vLocSol->solution(1);
        dt = vLocSol->time(0) - vLocSol->time(1);
    }
    
    MathVector<dim> StdVel[numSCVF];
    number Rho_up[numSCVF];
    number Rho_d[numSCVF];
   
    
    std_vel(  u,  geo,  StdVel, m_imDensitySCV,Rho_up,Rho_d);
    
    bool Phase2[numSCVF];
    number mu_l, rho_l;
    number mu_g, rho_g;
    MathVector<dim> Source_1;
    MathVector<dim> Source_g;
    

//     loop Sub Control Volume Faces (SCVF)
    for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
    {
    //     get current SCVF
        const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

        ////////////////////////////////////////////////////
        // Diffusive Term (Momentum Equation)
        ////////////////////////////////////////////////////

    //     1. Interpolate Functional Matrix of velocity at ip
        MathMatrix<dim, dim> gradVel;
        for(int d1 = 0; d1 < dim; ++d1)
            for(int d2 = 0; d2 <dim; ++d2)
            {
            //    sum up contributions of each shape
                gradVel(d1, d2) = 0.0;
                for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                {
                    gradVel(d1, d2) += scvf.global_grad(sh)[d2] * u(d1, sh);
                    
                }
            }

    //    2. Compute flux
        MathVector<dim> diffFlux;

    //    Add (\nabla u) \cdot \vec{n}
        MatVecMult(diffFlux, gradVel, scvf.normal());

    //    Add (\nabla u)^T \cdot \vec{n}
        if(!m_bLaplace)
            TransposedMatVecMultAdd(diffFlux, gradVel, scvf.normal());

    //    scale by viscosity
        VecScale(diffFlux, diffFlux, (-1.0) * m_imDensitySCVF[ip]);

    //    3. Add flux to local defect
        for(int c = 0; c < dim; ++c)
        {
            vvvLinDef[ip][c][scvf.from()] += diffFlux[c];
            vvvLinDef[ip][c][scvf.to()  ] -= diffFlux[c];
        }
    }
        
}

//    computes the linearized defect w.r.t to the source
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
lin_def_sourceSCV(const LocalVector& u,
               std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
               const size_t nip)
{
    // get finite volume geometry
    static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//    reset the values for the linearized defect
    for(size_t ip = 0; ip < nip; ++ip)
        for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
            for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
                vvvLinDef[ip][c][sh] = 0.0;

    // loop Sub Control Volumes Faces (SCVF)
    for ( size_t ip = 0; ip < geo.num_scv(); ++ip ) {
        // get current SCVF
        const typename TFVGeom::SCV& scv = geo.scv( ip );
        //     get associated node
        const int co = scv.node_id();
        
        //     set lin defect
        for(size_t c = 0; c < dim; ++c)
        {
            vvvLinDef[ip][c][co] += scv.volume() * m_imDensitySCV[ip];
        }
    }
}


//	prepares the nodal velocities for the export parameter
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
ex_nodal_velocity(MathVector<dim> vValue[],
        const MathVector<dim> vGlobIP[],
        number time, int si,
        const LocalVector& u,
        GridObject* elem,
        const MathVector<dim> vCornerCoords[],
        const MathVector<TFVGeom::dim> vLocIP[],
        const size_t nip,
        bool bDeriv,
        std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
    if(bDeriv)
    {
        for(size_t ip = 0; ip < nip; ++ip)
            for(size_t c = 0; c < vvvDeriv[ip].size(); ++c)
                for(size_t sh = 0; sh < vvvDeriv[ip][c].size(); ++sh)
                    for(size_t d = 0; d < vvvDeriv[ip][c][sh].size(); ++d)
                        vvvDeriv[ip][c][sh][d] = 0.0;

    }
// 	Get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

//  reference dimension
	static const int refDim = ref_elem_type::dim;

//  number of shape functions
	static const size_t numSH = ref_elem_type::numCorners;
    //static const size_t numSCVF = TFVGeom::numSCVF;
    
//	FV1 SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		//  Loop dimensions
			for(int d = 0; d < dim; ++d)
			{
                vValue[ip][d] = 0.0;
			//	Loop the shape functions
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
				//	Inerpolate the value
					vValue[ip][d] += u(d, sh) * scvf.shape(sh);
                        
					if(bDeriv)
                    {
                        vvvDeriv[ip][d][sh][d] += scvf.shape(sh);
                                                
                    }
				}
			}
			if(bDeriv)
			{
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					VecSet(vvvDeriv[ip][_P_][sh],0.0);
				}
			}
		}
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		number vLocShape[numSH];

	//	Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		ReferenceMapping<ref_elem_type, dim> mapping(vCornerCoords);

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.shapes(vLocShape, vLocIP[ip]);

		//  Loop dimensions
			for(int d = 0; d < dim; ++d)
			{
			//	Loop the shape functions
				for(size_t sh = 0; sh < numSH; ++sh)
				{
				//	Inerpolate the value
					vValue[ip][d] += u(d, sh) * vLocShape[sh];
					if(bDeriv)
						vvvDeriv[ip][d][sh][d] = +vLocShape[sh];
				}
			}
			if(bDeriv)
			{
				for(size_t sh = 0; sh < numSH; ++sh)
				{
					VecSet(vvvDeriv[ip][_P_][sh],0.0);
				}
			}		
		}
	}
};

//	computes the gradient of the velocity for the export parameter
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
ex_velocity_grad(MathMatrix<dim, dim> vValue[],
        const MathVector<dim> vGlobIP[],
        number time, int si,
        const LocalVector& u,
        GridObject* elem,
        const MathVector<dim> vCornerCoords[],
        const MathVector<TFVGeom::dim> vLocIP[],
        const size_t nip,
        bool bDeriv,
        std::vector<std::vector<MathMatrix<dim, dim> > > vvvDeriv[])
{
// 	Get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

//  reference dimension
	static const int refDim = ref_elem_type::dim;

//  number of shape functions
	static const size_t numSH = ref_elem_type::numCorners;	
    
    bool m_scvf=false;
    bool m_scv=false;

    if (!(vGlobIP!=NULL && vLocIP!=NULL && vCornerCoords!=NULL) && nip!=0) UG_THROW("Error in ExportParameter: Navier-Stokes: VelocityGrad");

    Inter->integration_points(m_scvf, m_scv, elem,   vCornerCoords,   vGlobIP, nip);
        
        
    /*
    
    std::vector<number> vVolumeFraction(geo.num_scvf());
    //m_imDensitySCVF.user_data()->compute(&u, elem, vCornerCoords, bDeriv);
    
    (*m_imDensitySCVF.user_data())(&vVolumeFraction[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, NULL);
    //(*(m_imDensitySCVF.user_data()))(&vVolumeFraction[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, &u, NULL);
    

    //m_vDependentData[i]->compute(&u, elem, vCornerCoords, bDeriv);*/



//	FV1 SCVF ip
	if(m_scvf)
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		//  Loop dimensions for direction
			for(int d1 = 0; d1 < dim; ++d1)
			{
				//  Loop dimensions for derivative
				for(int d2 = 0; d2 <dim; ++d2)
				{
					vValue[ip](d1, d2) = 0.0;
					//	sum up contributions of each shape
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					{
						vValue[ip](d1, d2) += u(d1, sh)*scvf.global_grad(sh)[d2];		
                        //vValue[ip](d1, d2) +=  (m_imDensitySCV[sh] / m_imDensitySCVF[ip]) * u(d1, sh)*scvf.global_grad(sh)[d2];
                        if(bDeriv)
							vvvDeriv[ip][d1][sh](d1,d2) = scvf.global_grad(sh)[d2];
                            //vvvDeriv[ip][d1][sh](d1,d2) = (m_imDensitySCV[sh] / m_imDensitySCVF[ip]) * scvf.global_grad(sh)[d2];
					}
				}
			}
			if(bDeriv)
			{
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					MatSet(vvvDeriv[ip][_P_][sh],0.0);
				}
			}
		}
	}
//    FV1 SCVF ip
    else if(m_scv)
    {
    //    Loop Sub Control Volume Faces (SCVF)
        for (size_t ip = 0; ip < geo.num_scv(); ++ip)
        {
        //     Get current SCVF
            const typename TFVGeom::SCV& scv = geo.scv(ip);

        //  Loop dimensions for direction
            for(int d1 = 0; d1 < dim; ++d1)
            {
                //  Loop dimensions for derivative
                for(int d2 = 0; d2 <dim; ++d2)
                {
                    vValue[ip](d1, d2) = 0.0;
                    //    sum up contributions of each shape
                    for(size_t sh = 0; sh < scv.num_sh(); ++sh)
                    {
                        vValue[ip](d1, d2) += u(d1, sh)*scv.global_grad(sh)[d2];
                        //vValue[ip](d1, d2) += (m_imDensitySCV[sh] / m_imDensitySCV[ip]) * u(d1, sh)*scv.global_grad(sh)[d2];
                        if(std::isnan(vValue[ip](d1, d2)))
                        {
                            UG_LOG("Error SCV\n");
                        }
                        if(bDeriv)
                            vvvDeriv[ip][d1][sh](d1,d2) = scv.global_grad(sh)[d2];
                            //vvvDeriv[ip][d1][sh](d1,d2) = (m_imDensitySCV[sh] / m_imDensitySCV[ip]) * scv.global_grad(sh)[d2];
                    }
                }
            }
            if(bDeriv)
            {
                for(size_t sh = 0; sh < scv.num_sh(); ++sh)
                {
                    MatSet(vvvDeriv[ip][_P_][sh],0.0);
                }
            }
        }
    }
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		MathVector<refDim> vLocGrad[numSH];
        MathVector<refDim> locGrad;
        MathVector<refDim> vGlobGrad;

	//	Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		ReferenceMapping<ref_elem_type, dim> mapping(vCornerCoords);

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.grads(vLocGrad, vLocIP[ip]);
            
            //    compute global grad
            mapping.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
            
            MatSet(vValue[ip],0.0);
            for(size_t sh = 0; sh < numSH; ++sh)
            {
                MatVecMult(vGlobGrad, JTInv, vLocGrad[sh]);
                
            //  Loop dimensions for direction
                for(int d1 = 0; d1 < dim; ++d1)
                {
                    for(int d2 = 0; d2 <dim; ++d2)
                    {
                        vValue[ip](d1, d2) += u(d1, sh)*vGlobGrad[d2];
                        
                        if(bDeriv)
                            vvvDeriv[ip][d1][sh](d1,d2) = vGlobGrad[d2];
                        
                    }
                    
                    
                }
                if(bDeriv)
                {

                    MatSet(vvvDeriv[ip][_P_][sh],0.0);

                }
                
                
            }


		}
        
        /*if( true && nip >0)
        {
            

            
            printf("NO (SCVF and SCV)\n");
            printf("Nip = %zu\n",nip);
            
            printf("vVelocity[0] = %f     %f\n",u(0, 0),u(1, 0));
            printf("vVelocity[1] = %f     %f\n",u(0, 1),u(1, 1));
            printf("vVelocity[2] = %f     %f\n",u(0, 2),u(1, 2));
            printf("vVelocity[3] = %f     %f\n",u(0, 3),u(1, 3));
            
            printf("Coordinates[0] = %f     %f\n",vCornerCoords[0][0],vCornerCoords[0][1]);
            printf("Coordinates[1] = %f     %f\n",vCornerCoords[1][0],vCornerCoords[1][1]);
            printf("Coordinates[2] = %f     %f\n",vCornerCoords[2][0],vCornerCoords[2][1]);
            printf("Coordinates[3] = %f     %f\n",vCornerCoords[3][0],vCornerCoords[3][1]);
            
            printf("vGlobIP[0] = %f     %f\n",vGlobIP[0][0],vGlobIP[0][1]);
            printf("vGlobIP[1] = %f     %f\n",vGlobIP[1][0],vGlobIP[1][1]);
            
            printf("vLocIP[0] = %f     %f\n",vLocIP[0][0],vLocIP[0][1]);
            printf("vLocIP[1] = %f     %f\n",vLocIP[1][0],vLocIP[1][1]);
            
            printf("VelGrad[0]  = \n");
            printf("   %f     %f\n", vValue[0](0, 0),vValue[0](0, 1));
            printf("   %f     %f\n", vValue[0](1, 0),vValue[0](1, 1));

            printf("VelGrad[1]  = \n");
            printf("   %f     %f\n", vValue[1](0, 0),vValue[1](0, 1));
            printf("   %f     %f\n", vValue[1](1, 0),vValue[1](1, 1));

            
            

        }*/
	}
};

//    computes the velocities at scvf ips for the export parameter
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
ex_velocity_ip(MathVector<dim> vValue[],
        const MathVector<dim> vGlobIP[],
        number time, int si,
        const LocalVector& u,
        GridObject* elem,
        const MathVector<dim> vCornerCoords[],
        const MathVector<TFVGeom::dim> vLocIP[],
        const size_t nip,
        bool bDeriv,
        std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
//     Only first order implemented
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");
	
	if(bDeriv)
	{
		for(size_t ip = 0; ip < nip; ++ip)
			for(size_t c = 0; c < vvvDeriv[ip].size(); ++c)
				for(size_t sh = 0; sh < vvvDeriv[ip][c].size(); ++sh)
					for(size_t d = 0; d < vvvDeriv[ip][c][sh].size(); ++d)
						vvvDeriv[ip][c][sh][d] = 0.0;

	}

//     get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
//    reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
//  reference dimension
	static const int refDim = ref_elem_type::dim;
	static const size_t numSCVF = TFVGeom::numSCVF;
	static const size_t numSh = reference_element_traits<TElem>::reference_element_type::numCorners;
   
	if(vLocIP == geo.scvf_local_ips())
	{


	//    check for solutions to pass to stabilization in time-dependent case
		const LocalVector *pSol = &u, *pOldSol = NULL;
		number dt = 0.0;
		if(this->is_time_dependent())
		{
		//    get and check current and old solution
			const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
			if(vLocSol->size() != 2)
				UG_THROW("NavierStokes::add_def_A_elem: "
								" Stabilization needs exactly two time points.");

		//    remember local solutions
			pSol = &vLocSol->solution(0);
			pOldSol = &vLocSol->solution(1);
			dt = vLocSol->time(0) - vLocSol->time(1);
		}

	//    interpolate velocity at ip with standard lagrange interpolation
		
		MathVector<dim> StdVel[numSCVF];
		MathVector<dim> Vel_ip[numSCVF];
		number Rho_up[numSCVF];
		number Rho_d[numSCVF];
	   
		
		std_vel(  u,  geo,  StdVel, m_imDensitySCV,Rho_up,Rho_d);
		
		
	//    compute stabilized velocities and shapes for continuity equation
		// \todo: (optional) Here we can skip the computation of shapes, implement?
		bool Phase2[numSCVF];
		
		number mu_l, rho_l;
		number mu_g, rho_g;
		MathVector<dim> Source_l;
		MathVector<dim> Source_g;
		
		m_spStab->update(&geo, *pSol, StdVel, m_bStokes, m_imKinViscosity, m_imKinViscositySCV, m_imDensitySCVF, NULL, m_imDensitySCV, NULL, m_imSourceSCVF, m_imSourceSCV, pOldSol, dt);
		
		
		//    get a const (!!) reference to the stabilization
		const INavierStokesFV1Stabilization<dim>& stab = *m_spStab;
		


		
	//    Loop Sub Control Volume Faces (SCVF)
		for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		//     Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

			vValue[ip] = stab.stab_vel(ip);

			if(bDeriv)
			{
				//    Loop the shape functions
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					
					//    Add derivative of stabilized flux w.r.t velocity comp to local matrix
					
					for(int d1 = 0; d1 < dim; ++d1)
					{
						if(stab.vel_comp_connected())
						{
							for(int d2 = 0; d2 < dim; ++d2)
							{
								vvvDeriv[ip][d2][sh][d1] +=  stab.stab_shape_vel(ip, d1, d2, sh);
								
							}
							

						}
						else
						{
							vvvDeriv[ip][d1][sh][d1] +=  stab.stab_shape_vel(ip, d1, d1, sh);
						}
						
						//    Add derivative of stabilized flux w.r.t pressure to local matrix
						vvvDeriv[ip][_P_][sh][d1] +=  stab.stab_shape_p(ip, d1, sh);
						
						

					}

					
				}
			}
			
		}
	}
//     general case
	else
	{
		std::vector<number> vViscosity(nip);
		LocalVector* uu = const_cast<LocalVector*>(&u);
		(*m_imKinViscosity.user_data())(&vViscosity[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, uu, NULL);

		
		const number D = ElementDiameter<GridObject, TDomain>(*elem, *this->domain());

	//    get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//    storage for shape function at ip
		number vLocShape[numSh];

	//    Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		ReferenceMapping<ref_elem_type, dim> mapping(vCornerCoords);

	//    loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//    evaluate at shapes at ip
			rTrialSpace.shapes(vLocShape, vLocIP[ip]);
			
			number MaxShape = 0.0;
			size_t node_id = 0;
			MathVector<dim> stdVel(0.0);
			bool boolShape = false;
			
			if (!m_bStokes)
			{
				for(size_t sh = 0; sh < numSh; ++sh)
				{
					
					if(vLocShape[sh] > MaxShape)
					{
						node_id = sh;
						MaxShape = vLocShape[sh];
						boolShape = true;
					}
					for(int d = 0; d < dim; ++d)
						stdVel[d] += u(d, sh) * vLocShape[sh];
					

				}

			}
			const number ViscScale = vViscosity[ip]/pow(D,2);
			const number ConvScale = VecLength(stdVel)/D;
			number inv_scale = ViscScale;
			if (boolShape) inv_scale += ConvScale;

		//  Loop dimensions
			for(int d = 0; d < dim; ++d)
			{
			//    Loop the shape functions
				vValue[ip][d] = 0.0;
				for(size_t sh = 0; sh < numSh; ++sh)
				{
				//    Inerpolate the value
					vValue[ip][d] += (ViscScale / inv_scale) * u(d, sh) * vLocShape[sh];
		
					if(bDeriv)
						vvvDeriv[ip][d][sh] += (ViscScale / inv_scale) * vLocShape[sh];
				}
				
				if ( boolShape)
				{
				//    Inerpolate the value
					vValue[ip][d] += (ConvScale / inv_scale) * u(d, node_id) ;
		
					if(bDeriv)
						vvvDeriv[ip][d][node_id] += (ConvScale / inv_scale) ;
				}

			}
			
			if(bDeriv)
			{
				for(size_t sh = 0; sh < numSh; ++sh)
				{
					VecSet(vvvDeriv[ip][_P_][sh],0.0);
				}
			}
		}
	}
};
//    computes the nodal pressure for the export parameter
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
ex_nodal_pressure(number vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         const LocalVector& u,
         GridObject* elem,
         const MathVector<dim> vCornerCoords[],
         const MathVector<TFVGeom::dim> vLocIP[],
         const size_t nip,
         bool bDeriv,
         std::vector<std::vector<number> > vvvDeriv[])
{
//  get finite volume geometry
    static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//    reference element
    typedef typename reference_element_traits<TElem>::reference_element_type
            ref_elem_type;

//    number of shape functions
    static const size_t numSH =    ref_elem_type::numCorners;


//    FV1 SCVF ip
    if(vLocIP == geo.scvf_local_ips())
    {
    //    Loop Sub Control Volume Faces (SCVF)
        for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
        {
        //     Get current SCVF
            const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

        //    compute pressure at ip
            vValue[ip] = 0.0;
            for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                vValue[ip] += u(_P_, sh) * scvf.shape(sh);

        //    compute derivative w.r.t. to unknowns iff needed
            if(bDeriv)
            {
                for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                    vvvDeriv[ip][_P_][sh] = scvf.shape(sh);

                // do not forget that number of DoFs (== vvvDeriv[ip][_P_])
                // might be > scvf.num_sh() in case of hanging nodes!
                size_t ndof = vvvDeriv[ip][_P_].size();
                for (size_t sh = scvf.num_sh(); sh < ndof; ++sh)
                    vvvDeriv[ip][_P_][sh] = 0.0;
            }
        }
    }
//    FV1 SCV ip
    else if(vLocIP == geo.scv_local_ips())
    {
    //    Loop Sub Control Volumes (SCV)
        for(size_t ip = 0; ip < geo.num_scv(); ++ip)
        {
        //     Get current SCV
            const typename TFVGeom::SCV& scv = geo.scv(ip);

        //    get corner of SCV
            const size_t co = scv.node_id();

        //    solution at ip
            vValue[ip] = u(_P_, co);

        //    set derivatives if needed
            if(bDeriv)
            {
                size_t ndof = vvvDeriv[ip][_P_].size();
                for(size_t sh = 0; sh < ndof; ++sh)
                    vvvDeriv[ip][_P_][sh] = (sh==co) ? 1.0 : 0.0;
            }
        }
    }
//     general case
    else
    {
    //    get trial space
        LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

    //    storage for shape function at ip
        number vShape[numSH];

    //    loop ips
        for(size_t ip = 0; ip < nip; ++ip)
        {
        //    evaluate at shapes at ip
            rTrialSpace.shapes(vShape, vLocIP[ip]);

        //    compute concentration at ip
            vValue[ip] = 0.0;
            for(size_t sh = 0; sh < numSH; ++sh)
                vValue[ip] += u(_P_, sh) * vShape[sh];

        //    compute derivative w.r.t. to unknowns iff needed
        //    \todo: maybe store shapes directly in vvvDeriv
            if(bDeriv)
            {
                for(size_t sh = 0; sh < numSH; ++sh)
                    vvvDeriv[ip][_P_][sh] = vShape[sh];

                // beware of hanging nodes!
                size_t ndof = vvvDeriv[ip][_P_].size();
                for (size_t sh = numSH; sh < ndof; ++sh)
                    vvvDeriv[ip][_P_][sh] = 0.0;
            }
        }
    }
};

//    computes the gradient of the pressure for the export parameter
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
ex_pressure_grad(MathVector<dim> vValue[],
        const MathVector<dim> vGlobIP[],
        number time, int si,
        const LocalVector& u,
        GridObject* elem,
        const MathVector<dim> vCornerCoords[],
        const MathVector<TFVGeom::dim> vLocIP[],
        const size_t nip,
        bool bDeriv,
        std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
//     Get finite volume geometry
    static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//    reference element
    typedef typename reference_element_traits<TElem>::reference_element_type
            ref_elem_type;

//    reference dimension
    static const int refDim = ref_elem_type::dim;

//    number of shape functions
    static const size_t numSH =    ref_elem_type::numCorners;

//    FV1 SCVF ip
    if(vLocIP == geo.scvf_local_ips())
    {
    //    Loop Sub Control Volume Faces (SCVF)
        for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
        {
        //     Get current SCVF
            const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

            VecSet(vValue[ip], 0.0);

            for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                VecScaleAppend(vValue[ip], u(_P_, sh), scvf.global_grad(sh));

            if(bDeriv)
            {
                for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                    vvvDeriv[ip][_P_][sh] = scvf.global_grad(sh);

                // beware of hanging nodes!
                size_t ndof = vvvDeriv[ip][_P_].size();
                for (size_t sh = scvf.num_sh(); sh < ndof; ++sh)
                    vvvDeriv[ip][_P_][sh] = 0.0;
            }
        }
    }
//     general case
    else
    {
    //    get trial space
        LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

    //    storage for shape function at ip
        MathVector<refDim> vLocGrad[numSH];
        MathVector<refDim> locGrad;

    //    Reference Mapping
        MathMatrix<dim, refDim> JTInv;
        ReferenceMapping<ref_elem_type, dim> mapping(vCornerCoords);

    //    loop ips
        for(size_t ip = 0; ip < nip; ++ip)
        {
        //    evaluate at shapes at ip
            rTrialSpace.grads(vLocGrad, vLocIP[ip]);

        //    compute grad at ip
            VecSet(locGrad, 0.0);
            for(size_t sh = 0; sh < numSH; ++sh)
                VecScaleAppend(locGrad, u(_P_, sh), vLocGrad[sh]);

        //    compute global grad
            mapping.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
            MatVecMult(vValue[ip], JTInv, locGrad);

        //    compute derivative w.r.t. to unknowns iff needed
            if(bDeriv)
            {
                for(size_t sh = 0; sh < numSH; ++sh)
                    MatVecMult(vvvDeriv[ip][_P_][sh], JTInv, vLocGrad[sh]);

                // beware of hanging nodes!
                size_t ndof = vvvDeriv[ip][_P_].size();
                for (size_t sh = numSH; sh < ndof; ++sh)
                    vvvDeriv[ip][_P_][sh] = 0.0;
            }
        }
    }
};


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void NavierStokesFV1<Domain1d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_2
template<>
void NavierStokesFV1<Domain2d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_3
template<>
void NavierStokesFV1<Domain3d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, FV1Geometry<Prism, dim> >();
		register_func<Pyramid, FV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, FV1Geometry<Hexahedron, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesFV1: Hanging Nodes not implemented.")
	}
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void
NavierStokesFV1<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;
	static const int refDim = reference_element_traits<TElem>::dim;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 	id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	id, &T::template add_rhs_elem<TElem, TFVGeom>);
    //    set computation of linearized defect w.r.t
    
    m_imDensitySCV.     set_fct(id, this, &T::template lin_def_densitySCV<TElem, TFVGeom>);
    m_imDensitySCVF.    set_fct(id, this, &T::template lin_def_densitySCVF<TElem, TFVGeom>);
    m_imKinViscosity.   set_fct(id, this, &T::template lin_def_viscosity<TElem, TFVGeom>);
    m_imSourceSCV.      set_fct(id, this, &T::template lin_def_sourceSCV<TElem, TFVGeom>);
    
	m_exVelocity->    template set_fct<T,refDim>(id, this, &T::template ex_nodal_velocity<TElem, TFVGeom>);
	m_exVelocityGrad->template set_fct<T,refDim>(id, this, &T::template ex_velocity_grad<TElem, TFVGeom>);
    m_exVelocity_ip->template set_fct<T,refDim>(id, this, &T::template ex_velocity_ip<TElem, TFVGeom>);
    m_exPressure->    template set_fct<T,refDim>(id, this, &T::template ex_nodal_pressure<TElem, TFVGeom>);
    m_exPressureGrad->template set_fct<T,refDim>(id, this, &T::template ex_pressure_grad<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_2
template class NavierStokesFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesFV1<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
