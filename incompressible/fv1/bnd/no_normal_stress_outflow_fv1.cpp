/*
 * Copyright (c) 2012-2014:  G-CSC, Goethe University Frankfurt
 * Authors: Dmitry Logashenko, Andreas Vogel
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


#include "no_normal_stress_outflow_fv1.h"

#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesNoNormalStressOutflowFV1<TDomain>::
NavierStokesNoNormalStressOutflowFV1(SmartPtr< IncompressibleNavierStokesBase<TDomain> > spMaster)
: NavierStokesNoNormalStressOutflowBase<TDomain>(spMaster)
{
    m_imKinViscosity.set_comp_lin_defect(false);
    m_imDensity.set_comp_lin_defect(false);
    m_imSource.set_comp_lin_defect(false);
//	register imports
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensity);
    this->register_import(m_imSource);

//	initialize the imports from the master discretization
	set_kinematic_viscosity(spMaster->kinematic_viscosity ());
	set_density(spMaster->density ());
    set_source(spMaster->source ());
    
    m_imSource.set_rhs_part();

	//	update assemble functions
	register_all_funcs(false);
};


template<typename TDomain>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
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

////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

/**
 * Prepares the element loop for a given element type: computes the FV-geo, ...
 * Note that there are separate loops for every type of the grid elements.
 */
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	register subsetIndex at Geometry
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	Only first order implementation
	if(!(TFVGeom::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

//	check if kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("NavierStokesNoNormalStressOutflowFV1::prep_elem_loop:"
						" Kinematic Viscosity has not been set, but is required.\n");

//	check if Density has been set
	if(!m_imDensity.data_given())
		UG_THROW("NavierStokesNoNormalStressOutflowFV1::prep_elem_loop:"
						" Density has not been set, but is required.\n");

//	extract indices of boundary
	this->extract_scheduled_data();

//	request the subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
		geo.add_boundary_subset(*subsetIter);
}

/**
 * Finalizes the element loop for a given element type.
 */
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
fsh_elem_loop()
{
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	remove the bnd subsets
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
		geo.remove_boundary_subset(*subsetIter);
}


/**
 * General initializations of a given grid element for the assembling.
 */
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("NavierStokesNoNormalStressOutflowFV1::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	find and set the local and the global positions of the IPs for imports
	typedef typename TFVGeom::BF BF;
	typename std::vector<int>::const_iterator subsetIter;
	
	m_vLocIP.clear(); m_vGloIP.clear();
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
		const int bndSubset = *subsetIter;
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);
		for(size_t i = 0; i < vBF.size(); ++i)
		{
			m_vLocIP.push_back(vBF[i].local_ip());
			m_vGloIP.push_back(vBF[i].global_ip());
		}
	}
	// REMARK: The loop above determines the ordering of the integration points:
	// The "outer ordering" corresponds to the ordering of the subsets in
	// m_vBndSubSetIndex, and "inside" of this ordering, the ip's are ordered
	// according to the order of the boundary faces in the FV geometry structure.

	m_imKinViscosity.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imKinViscosity.set_global_ips(&m_vGloIP[0], m_vGloIP.size());
	
	m_imDensity.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imDensity.set_global_ips(&m_vGloIP[0], m_vGloIP.size());
    
    m_imSource.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
    m_imSource.set_global_ips(&m_vGloIP[0], m_vGloIP.size());
}

/// Assembling of the diffusive flux (due to the viscosity) in the Jacobian of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
diffusive_flux_Jac
(
	const size_t ip, // index of the integration point (for the viscosity)
	const BF& bf, // boundary face to assemble
	LocalMatrix& J, // local Jacobian to update
	const LocalVector& u // local solution
)
{
	MathMatrix<dim,dim> diffFlux, tang_diffFlux;
	MathVector<dim> normalStress;

	//UG_LOG("ip: " << ip << "\n");

	for(size_t sh = 0; sh < bf.num_sh(); ++sh) // loop shape functions
	{
		//UG_LOG("sh: " << sh << "\n");
	//	1. Compute the total flux
	//	- add \nabla u
		MatSet (diffFlux, 0);
		MatDiagSet (diffFlux, VecDot (bf.global_grad(sh), bf.normal()));
	
	//	- add (\nabla u)^T
		if(!m_spMaster->laplace())
			for (size_t d1 = 0; d1 < (size_t)dim; ++d1)
				for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
					diffFlux(d1,d2) += bf.global_grad(sh)[d1] * bf.normal()[d2];
	
	//	2. Subtract the normal part:
		tang_diffFlux = diffFlux;
		TransposedMatVecMult(normalStress, diffFlux, bf.normal ());
		for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
			for (size_t d1 = 0; d1 < (size_t)dim; ++d1)
				tang_diffFlux(d1,d2) -= bf.normal()[d1] * normalStress[d2];
	
	//	3. Scale by viscosity
		tang_diffFlux *= - m_imKinViscosity[ip] * m_imDensity [ip];
	
	//	4. Add flux to local Jacobian
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
				J(d1, bf.node_id(), d2, sh) += tang_diffFlux (d1, d2);
	}
}

/// Assembling of the diffusive flux (due to the viscosity) in the defect of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
diffusive_flux_defect
(
	const size_t ip, // index of the integration point (for the viscosity)
	const BF& bf, // boundary face to assemble
	LocalVector& d, // local defect to update
	const LocalVector& u // local solution
)
{
	MathMatrix<dim, dim> gradVel;
	MathVector<dim> diffFlux;

// 	1. Get the gradient of the velocity at ip
	for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
		{
		//	sum up contributions of each shape
			gradVel(d1, d2) = 0.0;
			for(size_t sh = 0; sh < bf.num_sh(); ++sh)
				gradVel(d1, d2) += bf.global_grad(sh)[d2] * u(d1, sh);
		}

//	2. Compute the total flux

//	- add (\nabla u) \cdot \vec{n}
	MatVecMult(diffFlux, gradVel, bf.normal());

//	- add (\nabla u)^T \cdot \vec{n}
	if(!m_spMaster->laplace())
		TransposedMatVecMultAdd(diffFlux, gradVel, bf.normal());

//	3. Subtract the normal part:
	VecScaleAppend (diffFlux, - VecDot (diffFlux, bf.normal()), bf.normal());

//	A4. Scale by viscosity
	diffFlux *= - m_imKinViscosity[ip] * m_imDensity [ip];
    
//	5. Add flux to local defect
	for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		d(d1, bf.node_id()) += diffFlux[d1];
}

/// Assembling of the convective flux (due to the quadratic inertial term) in the Jacobian of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
convective_flux_Jac
(
	const size_t ip, // index of the integration point (for the density)
	const BF& bf, // boundary face to assemble
	LocalMatrix& J, // local Jacobian to update
	const LocalVector& u, // local solution
    const MathVector<dim>& StdVel // velocity at ip
)
{
// The convection velocity according to the current approximation:
	number old_momentum_flux = VecDot (StdVel, bf.normal ()) * m_imDensity [ip];
	
// We assume that there should be no inflow through the outflow boundary:
	if (old_momentum_flux < 0)
		old_momentum_flux = 0;
	
//	Add flux to local Jacobian
	for(size_t sh = 0; sh < bf.num_sh(); ++sh)
	{
		number t = old_momentum_flux * bf.shape(sh);
		for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
			J(d1, bf.node_id(), d1, sh) += t;
	}
}

/// Assembling of the convective flux (due to the quadratic inertial term) in the defect of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
convective_flux_defect
(
	const size_t ip, // index of the integration point (for the density)
	const BF& bf, // boundary face to assemble
	LocalVector& d, // local defect to update
	const LocalVector& u, // local solution
	const MathVector<dim>& StdVel // velocity at ip
)
{
// The convection velocity according to the current approximation:
	number old_momentum_flux = VecDot (StdVel, bf.normal ()) * m_imDensity [ip];
	
// We assume that there should be no inflow through the outflow boundary:
	if (old_momentum_flux < 0)
		old_momentum_flux = 0;
	
// Add the flux to the defect:
	for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
		d(d1, bf.node_id()) += old_momentum_flux * StdVel[d1];
}
/// Assembling of the convective flux (due to the quadratic inertial term) in the Jacobian of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
pressure_flux_Jac
(
    const size_t ip, // index of the integration point (for the density)
    const BF& bf, // boundary face to assemble
    LocalMatrix& J, // local Jacobian to update
    const LocalVector& u // local solution
)
{
    MathVector<dim> n;
    number mag=sqrt(VecDot(bf.normal(),bf.normal()));
    VecScale(n,bf.normal(),mag);
    number flux;
    

    for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
        for(size_t sh = 0; sh < bf.num_sh(); ++sh)
        {
            flux=VecDot(bf.global_grad(sh),n);
            J(d1, bf.node_id(), _P_, sh) -= flux* n[d1] * bf.volume ();
            //J(d1, bf.node_id(), _P_, sh) -=  bf.shape(sh) * bf.normal ()[d1];
        }

}
/// Assembling of the convective flux (due to the quadratic inertial term) in the defect of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
pressure_flux_defect
(
    const size_t ip, // index of the integration point (for the density)
    const BF& bf, // boundary face to assemble
    LocalVector& d, // local defect to update
    const LocalVector& u, // local solution
    const number& pressure, // pressure at ip
    const MathVector<dim>& PressureGrad // velocity at ip
)
{
    //number static_pressure = pressure;//-0.5 * m_imDensity[ip] * VecProd(StdVel,StdVel);
// Add the flux to the defect:
    MathVector<dim> n;
    MathVector<dim> NormalGrad(0);
    number mag=sqrt(VecDot(bf.normal(),bf.normal()));
    VecScale(n,bf.normal(),mag);
    
    VecScaleAppend (NormalGrad, VecDot (PressureGrad, n), n);
    for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
        d(d1, bf.node_id()) -= NormalGrad[d1] * bf.volume();
        //d(d1, bf.node_id()) -= pressure* bf.normal ()[d1];
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	typedef typename TFVGeom::BF BF;

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	size_t ip = 0;
	for(subsetIter = m_vBndSubSetIndex.begin();
		subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;
		
	//	get the list of the ip's:
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);

	// 	loop the boundary faces
		typename std::vector<BF>::const_iterator bf;
		for(bf = vBF.begin(); bf != vBF.end(); ++bf, ++ip)
		{
            
        // A. Compute Velocity at ip
            MathVector<dim> stdVel(0.0);
            for(size_t sh = 0; sh < bf->num_sh(); ++sh)
                for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
                    stdVel[d1] += u(d1, sh) * bf->shape(sh);

            
		//	A. The momentum equation:
			diffusive_flux_Jac<BF> (ip, *bf, J, u);
			if (!m_spMaster->stokes ())
				convective_flux_Jac<BF> (ip, *bf, J, u, stdVel);
            //pressure_flux_Jac<BF> (ip, *bf, J, u);
			
		//	B. The continuity equation
            
            //const number D =   ElementDiameter<GridObject, TDomain>(*elem, *this->domain());
            //number scale = 1.0 / (m_imDensity[ip] * (VecLength(stdVel)/D + m_imKinViscosity[ip]/pow(D,2)));

			for(size_t sh = 0; sh < bf->num_sh(); ++sh) // loop shape functions
            {
                for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
                    J(_P_, bf->node_id (), d2, sh) += bf->shape(sh) * bf->normal()[d2]; //* m_imDensity [ip];
                
                //J(_P_, bf->node_id (), _P_, sh) += scale*VecProd(bf->global_grad(sh) , bf->normal());
            }
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implemented
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	typedef typename TFVGeom::BF BF;

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	size_t ip = 0;
	for(subsetIter = m_vBndSubSetIndex.begin();
		subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;
		
	//	get the list of the ip's:
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);

	// 	loop the boundary faces
		typename std::vector<BF>::const_iterator bf;
		for(bf = vBF.begin(); bf != vBF.end(); ++bf, ++ip)
		{
		// A. Compute Velocity at ip
			MathVector<dim> stdVel(0.0);
            MathVector<dim> PressureGrad(0.0);
            number pressure = 0;
            
            for(size_t sh = 0; sh < bf->num_sh(); ++sh)
            {
                for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
                {
                    stdVel[d1] += u(d1, sh) * bf->shape(sh);
                    PressureGrad[d1] += u(_P_, sh) * bf->global_grad(sh)[d1];
                }
                pressure += u(_P_, sh) * bf->shape(sh);
            }
                    

		// B. Momentum equation:
			diffusive_flux_defect<BF> (ip, *bf, d, u);
			if (!m_spMaster->stokes ())
				convective_flux_defect<BF> (ip, *bf, d, u, stdVel);
            //pressure_flux_defect<BF> (ip, *bf, d, u, pressure , PressureGrad);
		
		// c. Continuity equation:
            
            //const number D =   ElementDiameter<GridObject, TDomain>(*elem, *this->domain());
            //number scale = 1.0 / (m_imDensity[ip] * (VecLength(stdVel)/D + m_imKinViscosity[ip]/pow(D,2)));
            d(_P_, bf->node_id()) += VecDot (stdVel, bf->normal());// * m_imDensity[ip];
            //d(_P_, bf->node_id()) += scale*VecDot (PressureGrad, bf->normal());
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//     Only first order implemented
    UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");
    
//    if zero data given, return
    if(!m_imSource.data_given()) return;

//     get finite volume geometry
    static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    typedef typename TFVGeom::BF BF;

//     loop registered boundary segments
    typename std::vector<int>::const_iterator subsetIter;
    size_t ip = 0;
    for(subsetIter = m_vBndSubSetIndex.begin();
        subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
    {
    //    get subset index corresponding to boundary
        const int bndSubset = *subsetIter;
        
    //    get the list of the ip's:
        if(geo.num_bf(bndSubset) == 0) continue;
        const std::vector<BF>& vBF = geo.bf(bndSubset);

    //     loop the boundary faces
        typename std::vector<BF>::const_iterator bf;
        for(bf = vBF.begin(); bf != vBF.end(); ++bf, ++ip)
        {
            number pgh = VecDot(bf->global_ip(),m_imSource[ip]);

            for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
                d(d1, bf->node_id()) += - pgh * bf->normal ()[d1];
        }
    }
}

/// Assembling of the diffusive flux (due to the viscosity) in the defect of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
diffusive_flux_lin_defect
 (
     const size_t ip, // index of the integration point (for the density)
     const BF& bf, // boundary face to assemble
  std::vector<std::vector<number> > vvvLinDef[], // local defect to update
     const LocalVector& u // local solution
 )
{
    MathMatrix<dim, dim> gradVel;
    MathVector<dim> diffFlux;

//     1. Get the gradient of the velocity at ip
    for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
        for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
        {
        //    sum up contributions of each shape
            gradVel(d1, d2) = 0.0;
            for(size_t sh = 0; sh < bf.num_sh(); ++sh)
                gradVel(d1, d2) += bf.global_grad(sh)[d2] * u(d1, sh);
        }

//    2. Compute the total flux

//    - add (\nabla u) \cdot \vec{n}
    MatVecMult(diffFlux, gradVel, bf.normal());

//    - add (\nabla u)^T \cdot \vec{n}
    if(!m_spMaster->laplace())
        TransposedMatVecMultAdd(diffFlux, gradVel, bf.normal());

//    3. Subtract the normal part:
    VecScaleAppend (diffFlux, - VecDot (diffFlux, bf.normal()), bf.normal());

//    A4. Scale by viscosity
    //diffFlux *= - m_imKinViscosity[ip];

//    5. Add flux to local defect
    for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
        vvvLinDef[ip][d1][bf.node_id()]+= diffFlux[d1];
}
/// Assembling of the convective flux (due to the quadratic inertial term) in the defect of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
convective_flux_lin_defect
(
    const size_t ip, // index of the integration point (for the density)
    const BF& bf, // boundary face to assemble
 std::vector<std::vector<number> > vvvLinDef[], // local defect to update
    const LocalVector& u // local solution
)
{
    MathVector<dim> StdVel;
    number old_momentum_flux;
    VecSet(StdVel, 0);

// The convection velocity according to the current approximation:
    for(size_t sh = 0; sh < bf.num_sh(); ++sh)
        for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
            StdVel[d1] += u(d1, sh) * bf.shape(sh);
    old_momentum_flux = VecDot (StdVel, bf.normal ()) ;

// We assume that there should be no inflow through the outflow boundary:
    if (old_momentum_flux < 0)
        old_momentum_flux = 0;

// Add the flux to the defect:
    for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
        vvvLinDef[ip][d1][bf.node_id()]+= old_momentum_flux * StdVel[d1];
}
//    computes the linearized defect w.r.t to the Density
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
lin_def_density(const LocalVector& u, std::vector<std::vector<number> > vvvLinDef[], const size_t nip)
{
    
    UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

//     get finite volume geometry
    static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    typedef typename TFVGeom::BF BF;
    
    //    reset the values for the linearized defect
    for(size_t ip = 0; ip < nip; ++ip)
        for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
            for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
                vvvLinDef[ip][c][sh] = 0.0;
    //UG_LOG("Anfang add_def_A_elem");

//     loop registered boundary segments
    typename std::vector<int>::const_iterator subsetIter;
    size_t ip = 0;
    for(subsetIter = m_vBndSubSetIndex.begin();
        subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
    {
    //    get subset index corresponding to boundary
        const int bndSubset = *subsetIter;

    //    get the list of the ip's:
        if(geo.num_bf(bndSubset) == 0) continue;
        const std::vector<BF>& vBF = geo.bf(bndSubset);

    //     loop the boundary faces
        typename std::vector<BF>::const_iterator bf;
        for(bf = vBF.begin(); bf != vBF.end(); ++bf)
        {
            if (!m_spMaster->stokes ())
                convective_flux_lin_defect<BF> (ip, *bf, vvvLinDef, u);
                
        // Next IP:
            ip++;
        }
    }
}
//    computes the linearized defect w.r.t to the Density
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
lin_def_viscosity(const LocalVector& u, std::vector<std::vector<number> > vvvLinDef[], const size_t nip)
{
    
    UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

//     get finite volume geometry
    static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    typedef typename TFVGeom::BF BF;
    
    //    reset the values for the linearized defect
    for(size_t ip = 0; ip < nip; ++ip)
        for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
            for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
                vvvLinDef[ip][c][sh] = 0.0;
    //UG_LOG("Anfang add_def_A_elem");

//     loop registered boundary segments
    typename std::vector<int>::const_iterator subsetIter;
    size_t ip = 0;
    for(subsetIter = m_vBndSubSetIndex.begin();
        subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
    {
    //    get subset index corresponding to boundary
        const int bndSubset = *subsetIter;

    //    get the list of the ip's:
        if(geo.num_bf(bndSubset) == 0) continue;
        const std::vector<BF>& vBF = geo.bf(bndSubset);

    //     loop the boundary faces
        typename std::vector<BF>::const_iterator bf;
        for(bf = vBF.begin(); bf != vBF.end(); ++bf)
        {
            if (!m_spMaster->stokes ())
                diffusive_flux_lin_defect<BF> (ip, *bf, vvvLinDef, u);
                
        // Next IP:
            ip++;
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void NavierStokesNoNormalStressOutflowFV1<Domain1d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesNoNormalStressOutflowFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_2
template<>
void NavierStokesNoNormalStressOutflowFV1<Domain2d>::
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
		UG_THROW("NavierStokesNoNormalStressOutflowFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_3
template<>
void NavierStokesNoNormalStressOutflowFV1<Domain3d>::
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
		UG_THROW("NavierStokesNoNormalStressOutflowFV1: Hanging Nodes not implemented.")
	}
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void
NavierStokesNoNormalStressOutflowFV1<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 	id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	    id, &T::template add_rhs_elem<TElem, TFVGeom>);
}


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class NavierStokesNoNormalStressOutflowFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NavierStokesNoNormalStressOutflowFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesNoNormalStressOutflowFV1<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
