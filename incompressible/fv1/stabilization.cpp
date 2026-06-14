/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include <string>
#include <locale>

#include "stabilization.h"
#include "diffusion_length.h"

#include "common/math/math_vector_matrix/math_vector_functions.h"
#include "common/math/math_vector_matrix/math_matrix_functions.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

template <int dim>
SmartPtr<INavierStokesSRFV1Stabilization<dim> > CreateNavierStokesStabilization(const std::string& name)
{
    std::string n = TrimString(name);
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    
    if(n == "fields") return SmartPtr<NavierStokesFIELDSStabilization<dim> >(new NavierStokesFIELDSStabilization<dim>());
    if(n == "flow") return SmartPtr<NavierStokesFLOWStabilization<dim> >(new NavierStokesFLOWStabilization<dim>());
    if(n == "fields_2") return SmartPtr<NavierStokesFIELDS_2_Stabilization<dim> >(new NavierStokesFIELDS_2_Stabilization<dim>());
    if(n == "flow_2") return SmartPtr<NavierStokesFLOW_2_Stabilization<dim> >(new NavierStokesFLOW_2_Stabilization<dim>());
    if(n == "karimian") return SmartPtr<NavierStokesKARIMIANStabilization<dim> >(new NavierStokesKARIMIANStabilization<dim>());
    if(n == "no") return SmartPtr<NavierStokesNOStabilization<dim> >(new NavierStokesNOStabilization<dim>());
    
    UG_THROW("NavierStokes: stabilization type '"<<name<<"' not a valid name of a Schneider-Raw stabilization."
             " Options are: fields, flow");
}

/////////////////////////////////////////////////////////////////////////////
// Interface for Stabilization
/////////////////////////////////////////////////////////////////////////////

//	register a update function for a Geometry
template <int dim>

template <typename TFVGeom, typename TAssFunc>
void
INavierStokesFV1Stabilization<dim>::
register_update_func(TAssFunc func)
{
    //	get unique geometry id
    size_t id = GetUniqueFVGeomID<TFVGeom>();
    
    //	make sure that there is enough space
    if((size_t)id >= m_vUpdateFunc.size())
        m_vUpdateFunc.resize(id+1, NULL);
    
    //	set pointer
    m_vUpdateFunc[id] = (UpdateFunc)func;
}

/////////////////////////////////////////////////////////////////////////////
// Common functions for the Schneider-Raw-type Stabilizations
/////////////////////////////////////////////////////////////////////////////

template <int dim>
void
INavierStokesSRFV1Stabilization<dim>::
set_diffusion_length(std::string diffLength)
{
    std::string n = TrimString(diffLength);
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    
    if      (n == "raw")        m_diffLengthType = RAW;
    else if (n == "fivepoint")  m_diffLengthType = FIVEPOINT;
    else if (n == "cor")        m_diffLengthType = COR;
    else
        UG_THROW("Diffusion Length calculation method not found."
                 " Use one of [Raw, Fivepoint, Cor].");
}

template <int dim>
template <typename TFVGeom>
void
INavierStokesSRFV1Stabilization<dim>::
compute_diff_length(const TFVGeom& geo)
{
    // 	Compute Diffusion Length in corresponding IPs
    switch(m_diffLengthType)
    {
        case FIVEPOINT: NSDiffLengthFivePoint(m_vDiffLengthSqInv, geo); return;
        case RAW:       NSDiffLengthRaw(m_vDiffLengthSqInv, geo); return;
        case COR:       NSDiffLengthCor(m_vDiffLengthSqInv, geo); return;
        default: UG_THROW(" Diffusion Length type not found.");
    }
}

/////////////////////////////////////////////////////////////////////////////
// FIELDS
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesFIELDSStabilization<TDim>::
update(const FV1Geometry<TElem, dim>* geo,
       const LocalVector& vCornerValue,
       const MathVector<dim> vStdVel[],
       const bool bStokes,
       const DataImport<number, dim>& kinVisco,
       const DataImport<number, dim>& kinViscoSCV,
       const DataImport<number, dim>& density,
       const DataImport<number, dim>& density_old,
       const DataImport<number, dim>& densitySCV,
	   const DataImport<number, dim>& densitySCV_old,
	   const number ps[],
	   const MathVector<dim> vStdRelVel[],
	   const DataImport<MathVector<dim>, dim>& RelVelSCVF,
       const DataImport<MathVector<dim>, dim>& Source,
       const DataImport<MathVector<dim>, dim>& SourceSCV,
	   const DataImport<MathVector<dim>, dim>& PressGrad,
       const LocalVector* pvCornerValueOldTime, number dt,
	   const int jump_shape[],
	   const bool phase_2[],
	   const bool multiphase)
{
    
    //	abbreviation for pressure
    static const size_t _P_ = dim;
    
    //	Some constants
    static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
    static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
    
    //	compute upwind (no convective terms for the Stokes eq. => no upwind)
    if (! bStokes) this->compute_upwind(geo, vStdVel);
    
    //	compute diffusion length
    this->compute_diff_length(*geo);
    
    //	cache values
    number vViscoPerDiffLenSq[numIp];
    for(size_t ip = 0; ip < numIp; ++ip)
        vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);
    
    number vNormStdVelPerConvLen[numIp];
    if(!bStokes)
        for(size_t ip = 0; ip < numIp; ++ip)
            vNormStdVelPerConvLen[ip] = VecTwoNorm(vStdVel[ip]) / upwind_conv_length(ip);
    
    //	Find out if upwinded velocities depend on other ip velocities. In that case
    //	we have to solve a matrix system. Else the system is diagonal and we can
    //	compute the inverse directly
    
    //	diagonal case (i.e. upwind vel depend only on corner vel or no upwind)
    if(bStokes || !non_zero_shape_ip())
    {
        //	We can solve the systems ip by ip
        for(size_t ip = 0; ip < numIp; ++ip)
        {
            //	get SubControlVolumeFace
            const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
            
            //	First, we compute the contributions to the diagonal
            //	Note: - There is no contribution of the upwind vel to the diagonal
            //		    in this case, only for non-diag problems
            //		  - The diag does not depend on the dimension
            
            //	Diffusion part
            number diag = vViscoPerDiffLenSq[ip];
            
            //	Time part
            if(pvCornerValueOldTime != NULL)
                diag += 1./dt;
            
            //	Convective Term (no convective terms in the Stokes eq.)
            if (! bStokes)
                diag += vNormStdVelPerConvLen[ip];
            
            // 	Loop components of velocity
            for(size_t d = 0; d < (size_t)dim; d++)
            {
                //	Now, we can assemble the rhs. This rhs is assembled by all
                //	terms, that are non-dependent on the ip vel.
                //	Note, that we can compute the stab_shapes on the fly when setting
                //	up the system.
                
                //	Source
                number rhs = 0.0;
                if(Source.data_given())
                    rhs = (   density[ip] / density[ip]) * Source[ip][d];
                
                //	Time
                if(pvCornerValueOldTime != NULL)
                {
                    //	interpolate old time step
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                    
                    //	add to rhs
                    rhs += oldIPVel / dt;
                }
                
                //	loop shape functions
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                    //	Diffusion part
                    number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                    
                    //	Convective term
                    if (! bStokes) // no convective terms in the Stokes eq.
                        sumVel += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
                    
                    //	Add to rhs
                    rhs += sumVel * vCornerValue(d, k);
                    
                    //	set stab shape
                    stab_shape_vel(ip, d, d, k) = sumVel / diag;
                    
                    //	Pressure part
                    const number sumP = -1.0 * (scvf.global_grad(k))[d] / density[ip];
                    
                    //	Add to rhs
                    rhs += sumP * vCornerValue(_P_, k);
                    
                    //	set stab shape
                    stab_shape_p(ip, d, k) = sumP / diag;
                }
                
                //	Finally, the can invert this row
                stab_vel(ip)[d] = rhs / diag;
            }
        }
    }
    /// need to solve system
    else
    {
        // 	For the FIELDS stabilization, there is no connection between the
        //	velocity components. Thus we can solve a system of size=numIP for
        //	each component of the velocity separately. This results in smaller
        //	matrices, that we have to invert.
        
        //	First, we have to assemble the Matrix, that includes all connections
        //	between the ip velocity component. Note, that in this case the
        //	matrix is non-diagonal and we must invert it.
        //	The Matrix is the same for all dim-components, thus we assemble it only
        //	once
        
        //	size of the system
        static const size_t N = numIp;
        
        //	a fixed size matrix
        DenseMatrix< FixedArray2<number, N, N> > mat;
        
        //	reset all values of the matrix to zero
        mat = 0.0;
        
        //	Loop integration points
        for(size_t ip = 0; ip < numIp; ++ip)
        {
            //	Time part
            if(pvCornerValueOldTime != NULL)
                mat(ip, ip) += 1./dt;
            
            //	Diffusion part
            mat(ip, ip) += vViscoPerDiffLenSq[ip];
            
            //	cache this value
            const number scale = vNormStdVelPerConvLen[ip];
            
            //	Convective Term (standard)
            mat(ip, ip) += scale;
            
            //	Convective Term by upwind
            for(size_t ip2 = 0; ip2 < numIp; ++ip2)
                mat(ip, ip2) -= upwind_shape_ip(ip, ip2) * scale;
        }
        
        //	we now create a matrix, where we store the inverse matrix
        typename block_traits<DenseMatrix< FixedArray2<number, N, N> > >::inverse_type inv;
        
        //	get the inverse
        if(!GetInverse(inv, mat))
            UG_THROW("Could not compute inverse.");
        
        //	loop dimensions (i.e. components of the velocity)
        for(int d = 0; d < dim; ++d)
        {
            //	create two vectors
            DenseVector< FixedArray1<number, N> > contVel[numSh];
            DenseVector< FixedArray1<number, N> > contP[numSh];
            
            //	Now, we can create several vector that describes the contribution of the
            //	corner velocities and the corner pressure. For each of this contribution
            //	components, we will apply the inverted matrix to get the stab_shapes
            
            //	Loop integration points
            for(size_t ip = 0; ip < numIp; ++ip)
            {
                //	get SubControlVolumeFace
                const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
                
                //	loop shape functions
                for(size_t k = 0; k < numSh; ++k)
                {
                    //	Diffusion part
                    contVel[k][ip] = vViscoPerDiffLenSq[ip]	* scvf.shape(k);
                    
                    //	Convection part
                    contVel[k][ip] += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
                    
                    //	Pressure part
                    contP[k][ip] = -1.0 * (scvf.global_grad(k))[d] / density[ip];
                }
            }
            
            //	solution vector
            DenseVector< FixedArray1<number, N> > xVel, xP;
            
            //	compute all stab_shapes
            for(size_t k = 0; k < numSh; ++k)
            {
                //	apply for vel stab_shape
                MatMult(xVel, 1.0, inv, contVel[k]);
                
                //	apply for pressure stab_shape
                MatMult(xP, 1.0, inv, contP[k]);
                
                //	write values in data structure
                //\todo: can we optimize this, e.g. without copy?
                for(size_t ip = 0; ip < numIp; ++ip)
                {
                    //	write stab_shape for vel
                    stab_shape_vel(ip, d, d, k) = xVel[ip];
                    
                    //	write stab_shape for pressure
                    stab_shape_p(ip, d, k) = xP[ip];
                }
            }
            
            //	Finally, we can compute the values of the stabilized velocity for each
            //	integration point
            
            //	vector of all contributions
            DenseVector< FixedArray1<number, N> > f;
            
            //	Loop integration points
            for(size_t ip = 0; ip < numIp; ++ip)
            {
                //	get SubControlVolumeFace
                const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
                
                //	Source
                f[ip] = 0.0;
                if(Source.data_given())
                    f[ip] = (   density[ip] / density[ip]) * Source[ip][d];
                
                //	Time
                if(pvCornerValueOldTime != NULL)
                {
                    //	interpolate old time step
                    //	\todo: Is this ok? Or do we need the old stabilized vel ?
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                    
                    //	add to rhs
                    f[ip] += oldIPVel / dt;
                }
            }
            
            //	sum up all contributions of vel and p for rhs.
            for(size_t k = 0; k < numSh; ++k)
            {
                //	add velocity contribution
                VecScaleAdd(f, 1.0, f, vCornerValue(d, k), contVel[k]);
                
                //	add pressure contribution
                VecScaleAdd(f, 1.0, f, vCornerValue(_P_, k), contP[k]);
            }
            
            //	invert the system for all contributions
            DenseVector< FixedArray1<number, N> > x;
            MatMult(x, 1.0, inv, f);
            
            //	write values in data structure
            //\todo: can we optimize this, e.g. without copy?
            for(size_t ip = 0; ip < numIp; ++ip)
            {
                //	write stab_shape for vel
                stab_vel(ip)[d] = x[ip];
            }
        } // end dim loop
        
    } // end switch for non-diag
}

template <>
void NavierStokesFIELDSStabilization<1>::register_func()
{
    register_func<RegularEdge>();
}

template <>
void NavierStokesFIELDSStabilization<2>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
}

template <>
void NavierStokesFIELDSStabilization<3>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
    register_func<Tetrahedron>();
    register_func<Pyramid>();
    register_func<Prism>();
    register_func<Hexahedron>();
}

/////////////////////////////////////////////////////////////////////////////
// FLOW
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesFLOWStabilization<TDim>::
update(const FV1Geometry<TElem, dim>* geo,
       const LocalVector& vCornerValue,
       const MathVector<dim> vStdVel[],
       const bool bStokes,
       const DataImport<number, dim>& kinVisco,
       const DataImport<number, dim>& kinViscoSCV,
       const DataImport<number, dim>& density,
       const DataImport<number, dim>& density_old,
       const DataImport<number, dim>& densitySCV,
	   const DataImport<number, dim>& densitySCV_old,
	   const number ps[],
	   const MathVector<dim> vStdRelVel[],
	   const DataImport<MathVector<dim>, dim>& RelVelSCVF,
       const DataImport<MathVector<dim>, dim>& Source,
       const DataImport<MathVector<dim>, dim>& SourceSCV,
	   const DataImport<MathVector<dim>, dim>& PressGrad,
       const LocalVector* pvCornerValueOldTime, number dt,
	   const int jump_shape[],
	   const bool phase_2[],
	   const bool multiphase)
{
    //	abbreviation for pressure
    static const size_t _P_ = dim;
    
    //	Some constants
    static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
    static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
    
    //	compute upwind and downwind (no convective terms for the Stokes eq. => no upwind)
    if (! bStokes)
    {
        this->compute_upwind(geo, vStdVel);
        this->compute_downwind(geo, vStdVel);
    }
    
    //	compute diffusion length
    this->compute_diff_length(*geo);
    
    //	cache values
    number vViscoPerDiffLenSq[numIp];
    for(size_t ip = 0; ip < numIp; ++ip)
        vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);
    
    number vNormStdVelPerConvLen[numIp];
    number vNormStdVelPerDownLen[numIp];
    if(!bStokes)
        for(size_t ip = 0; ip < numIp; ++ip)
        {
            const number norm = VecTwoNorm(vStdVel[ip]);
            vNormStdVelPerConvLen[ip] = norm / upwind_conv_length(ip);
            vNormStdVelPerDownLen[ip] = norm / (downwind_conv_length(ip) + upwind_conv_length(ip));
        }
    
    //	Find out if upwinded velocities depend on other ip velocities. In that case
    //	we have to solve a matrix system. Else the system is diagonal and we can
    //	compute the inverse directly
    
    //	diagonal case (i.e. upwind vel depend only on corner vel or no upwind)
    if(bStokes || !non_zero_shape_ip())
    {
        //	Loop integration points
        for(size_t ip = 0; ip < numIp; ++ip)
        {
            //	get SubControlVolumeFace
            const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
            
            //	First, we compute the contributions to the diagonal
            //	Note: - There is no contribution of the upwind vel to the diagonal
            //		    in this case, only for non-diag problems
            //		  - The diag does not depend on the dimension
            
            //	the diagonal entry
            number diag = vViscoPerDiffLenSq[ip];
            
            //	Time part
            if(pvCornerValueOldTime != NULL)
                diag += 1./dt;
            
            //	Convective Term  (no convective terms in the Stokes eq.)
            if (! bStokes)
                diag += vNormStdVelPerConvLen[ip];
            
            // 	Loop components of velocity
            for(int d = 0; d < dim; d++)
            {
                //	Now, we can assemble the rhs. This rhs is assembled by all
                //	terms, that are non-dependent on the ip vel.
                //	Note, that we can compute the stab_shapes on the fly when setting
                //	up the system.
                
                //	Source
                number rhs = 0.0;
                if(Source.data_given())
                    rhs = (   density[ip] / density[ip]) * Source[ip][d];
                
                //	Time
                if(pvCornerValueOldTime != NULL)
                {
                    //	interpolate old time step
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                    
                    //	add to rhs
                    rhs += oldIPVel / dt;
                }
                
                //	loop shape functions
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                    //	Diffusion part
                    number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                    
                    //	Convective term (no convective terms in the Stokes eq.)
                    if (! bStokes)
                    {
                        sumVel += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
                        
                        sumVel += vNormStdVelPerDownLen[ip] *
                        (downwind_shape_sh(ip, k) - upwind_shape_sh(ip, k));
                    }
                    
                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        if(d2 == d) continue;
                        
                        sumVel -= vStdVel[ip][d2] * (scvf.global_grad(k))[d2];
                    }
                    
                    //	Add to rhs
                    rhs += sumVel * vCornerValue(d, k);
                    
                    //	set stab shape
                    stab_shape_vel(ip, d, d, k) = sumVel / diag;
                    
                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        if(d2 == d) continue;
                        
                        const number sumVel2 = vStdVel[ip][d] * (scvf.global_grad(k))[d2];
                        
                        rhs += sumVel2 * vCornerValue(d2, k);
                        
                        stab_shape_vel(ip, d, d2, k) = sumVel2 / diag;
                    }
                    
                    //	Pressure part
                    const number sumP = -1.0 * (scvf.global_grad(k))[d]  / density[ip];
                    
                    //	Add to rhs
                    rhs += sumP * vCornerValue(_P_, k);
                    
                    //	set stab shape
                    stab_shape_p(ip, d, k) = sumP / diag;
                }
                
                //	Finally, the can invert this row
                stab_vel(ip)[d] = rhs / diag;
            }
        }
    }
    /// need to solve system
    else
    {
        // 	For the FLOW stabilization, there is no connection between the
        //	velocity components. Thus we can solve a system of size=numIP for
        //	each component of the velocity separately. This results in smaller
        //	matrices, that we have to invert.
        
        //	First, we have to assemble the Matrix, that includes all connections
        //	between the ip velocity component. Note, that in this case the
        //	matrix is non-diagonal and we must invert it.
        //	The Matrix is the same for all dim-components. Thus we invert it only
        //	once.
        
        //	size of the system
        static const size_t N = numIp;
        
        //	a fixed size matrix
        DenseMatrix< FixedArray2<number, N, N> > mat;
        
        //	reset all values of the matrix to zero
        mat = 0.0;
        
        //	Loop integration points
        for(size_t ip = 0; ip < numIp; ++ip)
        {
            //	Time part
            if(pvCornerValueOldTime != NULL)
                mat(ip, ip) += 1./dt;
            
            //	Diffusion part
            mat(ip, ip) += vViscoPerDiffLenSq[ip];
            
            //	Convective Term (standard)
            mat(ip, ip) += vNormStdVelPerConvLen[ip];
            
            for(size_t ip2 = 0; ip2 < numIp; ++ip2)
            {
                //	Convective Term by upwind
                mat(ip, ip2) -= upwind_shape_ip(ip, ip2) * vNormStdVelPerConvLen[ip];
                
                //	correction of divergence error
                mat(ip,ip2) += vNormStdVelPerDownLen[ip] *
                (upwind_shape_ip(ip, ip2) - downwind_shape_ip(ip, ip2));
            }
        }
        
        //	we now create a matrix, where we store the inverse matrix
        typename block_traits<DenseMatrix< FixedArray2<number, N, N> > >::inverse_type inv;
        
        //	get the inverse
        if(!GetInverse(inv, mat))
            UG_THROW("Could not compute inverse.");
        
        
        //	create vectors
        DenseVector< FixedArray1<number, N> > contVel[dim][numSh];
        DenseVector< FixedArray1<number, N> > contP[numSh];
        DenseVector< FixedArray1<number, N> > xP;
        DenseVector< FixedArray1<number, N> > xVel;
        
        //	Now, we can create several vector that describes the contribution of the
        //	corner velocities. For each of this contribution
        //	components, we will apply the inverted matrix to get the stab_shapes
        
        //	Loop integration points
        for(int d = 0; d < dim; ++d)
        {
            for(size_t ip = 0; ip < numIp; ++ip)
            {
                //	get SubControlVolumeFace
                const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
                
                //	loop shape functions
                for(size_t k = 0; k < numSh; ++k)
                {
                    //	Pressure part
                    contP[k][ip] = -1.0 * (scvf.global_grad(k))[d] / density[ip];
                    
                    //	Diffusion part
                    contVel[d][k][ip] = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                    
                    //	Convection part
                    contVel[d][k][ip] += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
                    
                    //	terms for correction of divergence error
                    contVel[d][k][ip] += vNormStdVelPerDownLen[ip] *
                    (downwind_shape_sh(ip, k) - upwind_shape_sh(ip, k));
                    
                    
                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        if(d2 == d) continue;
                        
                        contVel[d][k][ip] -= vStdVel[ip][d2]
                        * (scvf.global_grad(k))[d2];
                        
                        contVel[d2][k][ip] = vStdVel[ip][d]
                        * (scvf.global_grad(k))[d2];
                    }
                }
            } // end ip
            
            //	compute all stab_shapes
            for(size_t k = 0; k < numSh; ++k)
            {
                //	apply for pressure stab_shape
                MatMult(xP, 1.0, inv, contP[k]);
                
                //	write stab_shape for pressure
                //\todo: can we optimize this, e.g. without copy?
                for(size_t ip = 0; ip < numIp; ++ip)
                    stab_shape_p(ip, d, k) = xP[ip];
                
                //	compute vel stab_shapes
                for(int d2 = 0; d2 < dim; ++d2)
                {
                    //	apply for vel stab_shape
                    MatMult(xVel, 1.0, inv, contVel[d2][k]);
                    
                    //	write stab_shape for vel
                    //\todo: can we optimize this, e.g. without copy?
                    for(size_t ip = 0; ip < numIp; ++ip)
                        stab_shape_vel(ip, d, d2, k) = xVel[ip];
                }
            }
            
            //	Finally, we can compute the values of the stabilized velocity for each
            //	integration point
            
            //	vector of all contributions
            DenseVector< FixedArray1<number, N> > f;
            
            //	sum up all contributions of vel and p for rhs.
            f = 0.0;
            for(size_t k = 0; k < numSh; ++k)
            {
                //	add velocity contribution
                for(int d2 = 0; d2 < dim; ++d2)
                    VecScaleAdd(f, 1.0, f, vCornerValue(d2, k), contVel[d2][k]);
                
                //	add pressure contribution
                VecScaleAdd(f, 1.0, f, vCornerValue(_P_, k), contP[k]);
            }
            
            //	Loop integration points
            for(size_t ip = 0; ip < numIp; ++ip)
            {
                //	get SubControlVolumeFace
                const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
                
                //	Source
                if(Source.data_given())
                    f[ip] += (   density[ip] / density[ip]) * Source[ip][d];
                
                //	Time
                if(pvCornerValueOldTime != NULL)
                {
                    //	interpolate old time step
                    //	\todo: Is this ok? Or do we need the old stabilized vel ?
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                    
                    //	add to rhs
                    f[ip] += oldIPVel / dt;
                }
            }
            
            //	invert the system for all contributions
            DenseVector< FixedArray1<number, N> > x;
            MatMult(x, 1.0, inv, f);
            
            //	write values in data structure
            //\todo: can we optimize this, e.g. without copy?
            for(size_t ip = 0; ip < numIp; ++ip)
            {
                //	write stab_shape for vel
                stab_vel(ip)[d] = x[ip];
            }
        } // end dim
        
    } // end switch for non-diag
}


template <>
void NavierStokesFLOWStabilization<1>::register_func()
{
    register_func<RegularEdge>();
}

template <>
void NavierStokesFLOWStabilization<2>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
}

template <>
void NavierStokesFLOWStabilization<3>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
    register_func<Tetrahedron>();
    register_func<Pyramid>();
    register_func<Prism>();
    register_func<Hexahedron>();
}
/////////////////////////////////////////////////////////////////////////////
// FLOW 2
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesFLOW_2_Stabilization<TDim>::
update(const FV1Geometry<TElem, dim>* geo,
       const LocalVector& vCornerValue,
       const MathVector<dim> vStdVel[],
       const bool bStokes,
       const DataImport<number, dim>& kinVisco,
       const DataImport<number, dim>& kinViscoSCV,
       const DataImport<number, dim>& density,
       const DataImport<number, dim>& density_old,
       const DataImport<number, dim>& densitySCV,
	   const DataImport<number, dim>& densitySCV_old,
	   const number ps[],
	   const MathVector<dim> vStdRelVel[],
	   const DataImport<MathVector<dim>, dim>& RelVelSCVF,
       const DataImport<MathVector<dim>, dim>& Source,
       const DataImport<MathVector<dim>, dim>& SourceSCV,
	   const DataImport<MathVector<dim>, dim>& PressGrad,
       const LocalVector* pvCornerValueOldTime, number dt,
	   const int jump_shape[],
	   const bool phase_2[],
	   const bool multiphase)
{
	//	abbreviation for pressure
	static const size_t _P_ = dim;
	//	Some constants
	static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
	static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
	
	//	compute upwind and downwind (no convective terms for the Stokes eq. => no upwind)
	if (! bStokes)
	{
		this->compute_upwind(geo, vStdVel);
		this->compute_downwind(geo, vStdVel);
	}
	
	
	//	compute diffusion length
	this->compute_diff_length(*geo);
	
	//	cache values
	number vViscoPerDiffLenSq[numIp];
	for(size_t ip = 0; ip < numIp; ++ip)
		vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);
	
	number vNormStdVelPerConvLen[numIp];
	number vNormStdVelPerDownLen[numIp];
	if(!bStokes)
		for(size_t ip = 0; ip < numIp; ++ip)
		{
			const number norm = VecTwoNorm(vStdVel[ip]);
			vNormStdVelPerConvLen[ip] = norm / upwind_conv_length(ip);
			vNormStdVelPerDownLen[ip] = norm / (downwind_conv_length(ip) + upwind_conv_length(ip));
		}
	
	//	Find out if upwinded velocities depend on other ip velocities. In that case
	//	we have to solve a matrix system. Else the system is diagonal and we can
	//	compute the inverse directly
	
	//	diagonal case (i.e. upwind vel depend only on corner vel or no upwind)
	if(bStokes || !non_zero_shape_ip())
	{
		//	Loop integration points
		for(size_t ip = 0; ip < numIp; ++ip)
		{
			//	get SubControlVolumeFace
			const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
			
			//	First, we compute the contributions to the diagonal
			//	Note: - There is no contribution of the upwind vel to the diagonal
			//		    in this case, only for non-diag problems
			//		  - The diag does not depend on the dimension
			
			//	the diagonal entry
			number diag = vViscoPerDiffLenSq[ip];
			
			//	Time part
			if(pvCornerValueOldTime != NULL)
				diag += 1./dt;
			
			//	Convective Term  (no convective terms in the Stokes eq.)
			if (! bStokes)
				diag += vNormStdVelPerConvLen[ip];
			
			// 	Loop components of velocity
			for(int d = 0; d < dim; d++)
			{
				//	Now, we can assemble the rhs. This rhs is assembled by all
				//	terms, that are non-dependent on the ip vel.
				//	Note, that we can compute the stab_shapes on the fly when setting
				//	up the system.
				
				//	Source
				number rhs = 0.0;
				if(Source.data_given())
					rhs = (   density[ip] / density[ip]) * Source[ip][d];
				
				//	Time
				if(pvCornerValueOldTime != NULL)
				{
					//	interpolate old time step
					number oldIPVel = 0.0;
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
					
					//	add to rhs
					rhs += oldIPVel / dt;
				}
				
				//	loop shape functions
				for(size_t k = 0; k < scvf.num_sh(); ++k)
				{
					//	Diffusion part
					number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
					
					//	Convective term (no convective terms in the Stokes eq.)
					if (! bStokes)
					{
						sumVel += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
						
						sumVel += vNormStdVelPerDownLen[ip] *
						(downwind_shape_sh(ip, k) - upwind_shape_sh(ip, k));
					}
					
					for(int d2 = 0; d2 < dim; ++d2)
					{
						if(d2 == d) continue;
						
						sumVel -= vStdVel[ip][d2] * (scvf.global_grad(k))[d2];
					}
					
					//	Add to rhs
					rhs += sumVel * vCornerValue(d, k);
					
					//	set stab shape
					stab_shape_vel(ip, d, d, k) = sumVel / diag;
					
					for(int d2 = 0; d2 < dim; ++d2)
					{
						if(d2 == d) continue;
						
						const number sumVel2 = vStdVel[ip][d] * (scvf.global_grad(k))[d2];
						
						rhs += sumVel2 * vCornerValue(d2, k);
						
						stab_shape_vel(ip, d, d2, k) = sumVel2 / diag;
					}
					
					//	Pressure part
					const number sumP = -1.0 * (scvf.global_grad(k))[d]  / density[ip];
					
					//	Add to rhs
					rhs += sumP * vCornerValue(_P_, k);
					
					//	set stab shape
					stab_shape_p(ip, d, k) = sumP / diag;
				}
				
				//	Finally, the can invert this row
				stab_vel(ip)[d] = rhs / diag;
			}
		}
	}
	/// need to solve system
	else
	{
		// 	For the FLOW stabilization, there is no connection between the
		//	velocity components. Thus we can solve a system of size=numIP for
		//	each component of the velocity separately. This results in smaller
		//	matrices, that we have to invert.
		
		//	First, we have to assemble the Matrix, that includes all connections
		//	between the ip velocity component. Note, that in this case the
		//	matrix is non-diagonal and we must invert it.
		//	The Matrix is the same for all dim-components. Thus we invert it only
		//	once.
		
		//	size of the system
		static const size_t N = numIp;
		
		//	a fixed size matrix
		DenseMatrix< FixedArray2<number, N, N> > mat;
		
		//	reset all values of the matrix to zero
		mat = 0.0;
		
		//	Loop integration points
		for(size_t ip = 0; ip < numIp; ++ip)
		{
			//	Time part
			if(pvCornerValueOldTime != NULL)
				mat(ip, ip) += 1./dt;
			
			//	Diffusion part
			mat(ip, ip) += vViscoPerDiffLenSq[ip];
			
			//	Convective Term (standard)
			mat(ip, ip) += vNormStdVelPerConvLen[ip];
			
			for(size_t ip2 = 0; ip2 < numIp; ++ip2)
			{
				//	Convective Term by upwind
				mat(ip, ip2) -= upwind_shape_ip(ip, ip2) * vNormStdVelPerConvLen[ip];
				
				//	correction of divergence error
				mat(ip,ip2) += vNormStdVelPerDownLen[ip] *
				(upwind_shape_ip(ip, ip2) - downwind_shape_ip(ip, ip2));
			}
		}
		
		//	we now create a matrix, where we store the inverse matrix
		typename block_traits<DenseMatrix< FixedArray2<number, N, N> > >::inverse_type inv;
		
		//	get the inverse
		if(!GetInverse(inv, mat))
			UG_THROW("Could not compute inverse.");
		
		
		//	create vectors
		DenseVector< FixedArray1<number, N> > contVel[dim][numSh];
		DenseVector< FixedArray1<number, N> > contP[numSh];
		DenseVector< FixedArray1<number, N> > xP;
		DenseVector< FixedArray1<number, N> > xVel;
		
		//	Now, we can create several vector that describes the contribution of the
		//	corner velocities. For each of this contribution
		//	components, we will apply the inverted matrix to get the stab_shapes
		
		//	Loop integration points
		for(int d = 0; d < dim; ++d)
		{
			for(size_t ip = 0; ip < numIp; ++ip)
			{
				//	get SubControlVolumeFace
				const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
				
				//	loop shape functions
				for(size_t k = 0; k < numSh; ++k)
				{
					//	Pressure part
					contP[k][ip] = -1.0 * (scvf.global_grad(k))[d] / density[ip];
					
					//	Diffusion part
					contVel[d][k][ip] = vViscoPerDiffLenSq[ip] * scvf.shape(k);
					
					//	Convection part
					contVel[d][k][ip] += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
					
					//	terms for correction of divergence error
					contVel[d][k][ip] += vNormStdVelPerDownLen[ip] *
					(downwind_shape_sh(ip, k) - upwind_shape_sh(ip, k));
					
					
					for(int d2 = 0; d2 < dim; ++d2)
					{
						if(d2 == d) continue;
						
						contVel[d][k][ip] -= vStdVel[ip][d2]
						* (scvf.global_grad(k))[d2];
						
						contVel[d2][k][ip] = vStdVel[ip][d]
						* (scvf.global_grad(k))[d2];
					}
				}
			} // end ip
			
			//	compute all stab_shapes
			for(size_t k = 0; k < numSh; ++k)
			{
				//	apply for pressure stab_shape
				MatMult(xP, 1.0, inv, contP[k]);
				
				//	write stab_shape for pressure
				//\todo: can we optimize this, e.g. without copy?
				for(size_t ip = 0; ip < numIp; ++ip)
					stab_shape_p(ip, d, k) = xP[ip];
				
				//	compute vel stab_shapes
				for(int d2 = 0; d2 < dim; ++d2)
				{
					//	apply for vel stab_shape
					MatMult(xVel, 1.0, inv, contVel[d2][k]);
					
					//	write stab_shape for vel
					//\todo: can we optimize this, e.g. without copy?
					for(size_t ip = 0; ip < numIp; ++ip)
						stab_shape_vel(ip, d, d2, k) = xVel[ip];
				}
			}
			
			//	Finally, we can compute the values of the stabilized velocity for each
			//	integration point
			
			//	vector of all contributions
			DenseVector< FixedArray1<number, N> > f;
			
			//	sum up all contributions of vel and p for rhs.
			f = 0.0;
			for(size_t k = 0; k < numSh; ++k)
			{
				//	add velocity contribution
				for(int d2 = 0; d2 < dim; ++d2)
					VecScaleAdd(f, 1.0, f, vCornerValue(d2, k), contVel[d2][k]);
				
				//	add pressure contribution
				VecScaleAdd(f, 1.0, f, vCornerValue(_P_, k), contP[k]);
			}
			
			//	Loop integration points
			for(size_t ip = 0; ip < numIp; ++ip)
			{
				//	get SubControlVolumeFace
				const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
				
				//	Source
				if(Source.data_given())
					f[ip] += (   density[ip] / density[ip]) * Source[ip][d];
				
				//	Time
				if(pvCornerValueOldTime != NULL)
				{
					//	interpolate old time step
					//	\todo: Is this ok? Or do we need the old stabilized vel ?
					number oldIPVel = 0.0;
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
					
					//	add to rhs
					f[ip] += oldIPVel / dt;
				}
			}
			
			//	invert the system for all contributions
			DenseVector< FixedArray1<number, N> > x;
			MatMult(x, 1.0, inv, f);
			
			//	write values in data structure
			//\todo: can we optimize this, e.g. without copy?
			for(size_t ip = 0; ip < numIp; ++ip)
			{
				//	write stab_shape for vel
				stab_vel(ip)[d] = x[ip];
			}
		} // end dim
		
	} // end switch for non-diag
}

template <>
void NavierStokesFLOW_2_Stabilization<1>::register_func()
{
    register_func<RegularEdge>();
}

template <>
void NavierStokesFLOW_2_Stabilization<2>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
}

template <>
void NavierStokesFLOW_2_Stabilization<3>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
    register_func<Tetrahedron>();
    register_func<Pyramid>();
    register_func<Prism>();
    register_func<Hexahedron>();
}
/////////////////////////////////////////////////////////////////////////////
// FIELDS 2
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesFIELDS_2_Stabilization<TDim>::
update(const FV1Geometry<TElem, dim>* geo,
	   const LocalVector& vCornerValue,
	   const MathVector<dim> vStdVel[],
	   const bool bStokes,
	   const DataImport<number, dim>& kinVisco,
	   const DataImport<number, dim>& kinViscoSCV,
	   const DataImport<number, dim>& density,
	   const DataImport<number, dim>& density_old,
	   const DataImport<number, dim>& densitySCV,
	   const DataImport<number, dim>& densitySCV_old,
	   const number ps[],
	   const MathVector<dim> vStdRelVel[],
	   const DataImport<MathVector<dim>, dim>& RelVelSCVF,
	   const DataImport<MathVector<dim>, dim>& Source,
	   const DataImport<MathVector<dim>, dim>& SourceSCV,
	   const DataImport<MathVector<dim>, dim>& PressGrad,
	   const LocalVector* pvCornerValueOldTime, number dt,
	   const int jump_shape[],
	   const bool phase_2[],
	   const bool multiphase)
{
	if( non_zero_shape_ip())
	{
		UG_THROW("Not implemented for ip velocities depending on other ip.");
	}
	
	//    abbreviation for pressure
	static const size_t _P_ = dim;
	//    abbreviation for VolumeFraction
	static const size_t _C_ = dim+1;
	
	//    Some constants
	static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
	static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
	
	for(size_t ip = 0; ip < numIp; ++ip)
	{
		for(size_t k = 0; k < numSh; ++k)
		{
			for(int d1 = 0; d1 < dim; d1++)
			{
				for(int d2 = 0; d2 < dim; d2++)
				{
					stab_shape_vel(ip, d2, d1, k) = 0.0;
				}
				stab_shape_p(ip, d1, k) = 0.0;
				stab_shape_c(ip, d1, k) = 0.0;
			}
		}
	}
	
	//    compute diffusion length
	this->compute_diff_length(*geo);

	MathVector<dim> vStdVel_ip_old[numIp];
	
	
	/*if(pvCornerValueOldTime != NULL )
	{
		for(size_t ip = 0; ip < numIp; ++ip)
		{
			const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
			
			VecSet(vStdVel_ip_old[ip],0.0);
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				for(int d = 0; d < dim; d++)
					vStdVel_ip_old[ip][d] +=  scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
		}
		
		if ( !bStokes )
		{
			
			this->compute_upwind(geo, vStdVel_ip_old);
			this->compute_downwind(geo, vStdVel_ip_old);
			

			for(size_t ip = 0; ip < numIp; ++ip)
			{
				const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
				
				number vViscoPerDiffLenSq_old = kinViscoSCV[ip] * diff_length_sq_inv(ip);
				number vNormStdVelPerConvLen_old = VecTwoNorm(vStdVel_ip_old[ip]) / upwind_conv_length(ip);
				
				number Rho_up = 0.0;
				number Rho_do = 0.0;


				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					Rho_up += upwind_shape_sh(ip, sh) * densitySCV_old[sh];
					Rho_do += downwind_shape_sh(ip, sh) * densitySCV_old[sh];
				}
				number Ratio_rho = factor*pow(fmin(Rho_up , Rho_do) / fmax (Rho_up , Rho_do), power);
				
				number R_rho_up = Rho_up * downwind_conv_length(ip) / ( Rho_up * downwind_conv_length(ip) + Rho_do * upwind_conv_length(ip));
				number Ratio_rho_do = 1.0 - R_rho_up;
				
				MathVector<dim> Vel_ip;
				VecSet(Vel_ip,0.0);
				for(int d = 0; d < dim; d++)
				{
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					{
						Vel_ip[d] += (Ratio_rho*upwind_shape_sh(ip, sh)  +  (1.0-Ratio_rho)*(R_rho_up * upwind_shape_sh(ip, sh) + Ratio_rho_do * downwind_shape_sh(ip, sh))) * (*pvCornerValueOldTime)(d, sh);
					}
					
				}
				number diag_old = vViscoPerDiffLenSq_old + vNormStdVelPerConvLen_old;
				VecScaleAdd(vStdVel_ip_old[ip], vViscoPerDiffLenSq_old/diag_old, vStdVel_ip_old[ip], vNormStdVelPerConvLen_old/diag_old,Vel_ip);
				
				
			}
			

			
		}

	}*/
	
	
	
	
	//    compute upwind and downwind (no convective terms for the Stokes eq. => no upwind)
	if (! bStokes)
	{
		this->compute_upwind(geo, vStdVel);
		this->compute_downwind(geo, vStdVel);
		
		/*if(RelVelSCVF.data_given())
		{
			this->compute_upwind_rel(geo, vStdRelVel);
			this->compute_downwind_rel(geo, vStdRelVel);
		}*/
			
	}
	
	MathVector<dim> RhoGrad[numIp];
	//MathVector<dim> ViscGrad[numIp];
	//MathVector<dim> Vel[numIp];
	//number DenMomentum[numIp];
	number RHO_up[numIp];
	number RHO_do[numIp];
	number MU_scvf[numIp];
	number Ratio_rho_up[numIp];
	number Ratio_rho_do[numIp];
	number Ratio[numIp];
	number power = 1.0;
	number factor = 1.0;
	if (! bStokes)
	{
		for(size_t ip = 0; ip < numIp; ++ip)
		{
			const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
			const size_t from = scvf.from();
			const size_t to = scvf.to();
			VecSet(RhoGrad[ip], 0.0);
			//VecSet(ViscGrad[ip], 0.0);
			//DenMomentum[ip] = 0.0;
			RHO_up[ip] = 0.0;
			RHO_do[ip] = 0.0;
			MU_scvf[ip] = 0.0;
			Ratio_rho_up[ip] = 0.0;
			number dRho = fabs(densitySCV[to]- densitySCV[from]);
			number theta = pow( dRho/(1.0+dRho), power);
			//const number Val = +VecTwoNorm(vStdVel[ip]) / (downwind_conv_length(ip) + upwind_conv_length(ip));
			//MU_scvf[ip] = 0.5*(densitySCV[from] * kinViscoSCV[from] + densitySCV[to] * kinViscoSCV[to]);
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				VecScaleAppend(RhoGrad[ip], densitySCV[sh], scvf.global_grad(sh));
				//VecScaleAppend(ViscGrad[ip], kinViscoSCV[sh]*densitySCV[sh], scvf.global_grad(sh));
				//DenMomentum[ip] += Val * (downwind_shape_sh(ip, sh) - upwind_shape_sh(ip, sh)) * densitySCV[sh];

				RHO_up[ip] += upwind_shape_sh(ip, sh) * densitySCV[sh];
				RHO_do[ip] += downwind_shape_sh(ip, sh) * densitySCV[sh];
				
				MU_scvf[ip] += upwind_shape_sh(ip, sh) *(densitySCV[sh] * kinViscoSCV[sh]);
				
				
			}
			MU_scvf[ip] = theta * MU_scvf[ip]   +(1.0- theta) * density[ip] * kinVisco[ip];
			number diff = fabs(RHO_up[ip] - RHO_do[ip]);
			Ratio[ip] = 0.0;//factor*pow(diff /(1.0 + diff), power);
			//if(Ratio[ip] < 0.8)printf("Ratio[%zu] = %f \n", ip,Ratio[ip]);
			
			Ratio_rho_up[ip] = RHO_up[ip] * downwind_conv_length(ip) / ( RHO_up[ip] * downwind_conv_length(ip) + RHO_do[ip] * upwind_conv_length(ip));
			Ratio_rho_do[ip] = 1.0 - Ratio_rho_up[ip];
			/*for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				for(int d = 0; d < dim; d++)
				{
					Vel[ip][d] += vCornerValue(d, sh)*(Ratio_rho_up[ip] * upwind_shape_sh(ip, sh) + Ratio_rho_do[ip] * downwind_shape_sh(ip, sh));
				}
			}*/
			//DenMomentum[ip]=VecProd(RhoGrad[ip],vStdVel[ip]);
		}
	}
	

	bool boolSource = (SourceSCV.data_given()) ? true : false;
	
	//    cache values
	number vViscoPerDiffLenSq[numIp];
	MathVector<dim> SOURCE[numIp];
	
	number vNormStdVelPerConvLen[numIp];
	number vNormStdVelPerDownLen[numIp];
	number vNormRelVelPerConvLen[numIp];
	
	number vPecletScale[numIp];
	
	MathVector<dim> vStdVel_stab[numIp];
	
	for(size_t ip = 0; ip < numIp; ++ip)
	{
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
		
		vViscoPerDiffLenSq[ip] = density[ip] * kinVisco[ip] * diff_length_sq_inv(ip);
		if(boolSource) SOURCE[ip] = Source[ip];
		
		VecSet(vStdVel_stab[ip],0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			for(int d = 0; d < dim; d++)
			{
				vStdVel_stab[ip][d] += scvf.shape(sh) * vCornerValue(d, sh);
			}
			  
		}
		
		if(!bStokes)
		{
			number Value1 = density[ip] * VecTwoNorm(vStdVel[ip]) / upwind_conv_length(ip);
			
			vNormStdVelPerConvLen[ip] = Value1;
			
			vNormStdVelPerDownLen[ip] = density[ip] * VecTwoNorm(vStdVel[ip]) / (downwind_conv_length(ip) + upwind_conv_length(ip));
			
			
			number Value2 = density[ip] * VecTwoNorm(vStdVel_stab[ip]) / upwind_conv_length(ip);
			vPecletScale[ip] = Value2 / (Value2 + MU_scvf[ip]*diff_length_sq_inv(ip));
			
			
			
			/*if(RelVelSCVF.data_given())
			{
				const number rhos = Inter->Density_max();
				const number alpha_max = Inter->Alpha_max();
				vNormRelVelPerConvLen[ip] = (rhos/alpha_max) * VecTwoNorm(vStdRelVel[ip]) / (upwind_conv_length_rel(ip) + downwind_conv_length_rel(ip));

			}*/
			
		}


	}
	
	MathVector<dim> vConsGravitySCVF[numIp];
	if(Inter->boolConsistentGravity())
	{
		Inter-> template ConsistentGravitySCVF<TElem>(vConsGravitySCVF, *geo, geo->corners(), numIp, densitySCV.values());
	}
	
	
	//    Find out if upwinded velocities depend on other ip velocities. In that case
	//    we have to solve a matrix system. Else the system is diagonal and we can
	//    compute the inverse directly
	
	//    diagonal case (i.e. upwind vel depend only on corner vel or no upwind)

	//    Loop integration points
	for(size_t ip = 0; ip < numIp; ++ip)
	{
		//    get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
		const size_t from=scvf.from();
		const size_t to=scvf.to();
		MathVector<dim> Direction = 0.0;
		VecSubtract(Direction,geo->scv_global_ips()[to],geo->scv_global_ips()[from]);
		VecScale(Direction,Direction,1.0/VecLengthSq(Direction));
		
		
		if (true)
		{

			
			//    First, we compute the contributions to the diagonal
			//    Note: - There is no contribution of the upwind vel to the diagonal
			//            in this case, only for non-diag problems
			//          - The diag does not depend on the dimension
			
			//    the diagonal entry
			number diag;
			if(true)//!multiphase)
			{
				diag = vViscoPerDiffLenSq[ip];
			}
			else
			{
				diag = diff_factor(ip) * diff_length_sq_inv(ip) /density[ip];
			}
			
			//    Time part
			//if(pvCornerValueOldTime != NULL)
				//diag += density[ip]/dt;
			
			//    Convective Term  (no convective terms in the Stokes eq.)
			if (! bStokes)
			{
				diag += vNormStdVelPerConvLen[ip];
				
				//diag += DenMomentum[ip];
			}
			
			//diag += DenMomentum[ip];
			//     Loop components of velocity
			for(int d = 0; d < dim; d++)
			{
				//    Now, we can assemble the rhs. This rhs is assembled by all
				//    terms, that are non-dependent on the ip vel.
				//    Note, that we can compute the stab_shapes on the fly when setting
				//    up the system.
				
				//    Source
				number rhs = 0.0;
				number rhs_mu = 0.0;
				if(Inter->boolConsistentGravity())
				{
					rhs =  vConsGravitySCVF[ip][d];
				}
				else
				{
					if(boolSource)
					{
						rhs =  SOURCE[ip][d];
					}
				}
				
				/*if(PressGrad.data_given())
				{
					rhs += PressGrad[ip][d];
					
				}*/

				
				//    Time
				/*if(pvCornerValueOldTime != NULL)
				{
					//	interpolate old time step
					number oldIPVel = 0.0;
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
					//    add to rhs
					rhs += density[ip] * oldIPVel/ dt; //(density_old[ip] / density[ip])
				}*/

				
				/*if (! bStokes)
				{
					rhs +=  vStdVel[ip][d] * vNormStdVelPerDownLen[ip] * (RHO_do[ip]-RHO_up[ip]);//
					
				}*/
				
				//    loop shape functions
				for(size_t k = 0; k < scvf.num_sh(); ++k)
				{
					//    Diffusion part
					number sumVel = 0.0;
					if(!multiphase)
					{
						
						sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k) ;
						//sumVel = vViscoPerDiffLenSq[ip] * (densitySCV[k] / density[ip]) * scvf.shape(k);
						
					}
					else
					{
						number SumVelRHS = 0.0;
						number SumVel2;
						for(int d1 = 0; d1 < dim; ++d1)
						{
							SumVel2 = vViscoPerDiffLenSq[ip] * slip_vel_shape_vel( ip,  d,  d1, k)  ;
							SumVelRHS += SumVel2 * vCornerValue(d1, k);
							stab_shape_vel(ip, d, d1, k) += SumVel2 / diag;
							
						}
						rhs += SumVelRHS;
						
						
					}
					//number sumVel = vViscoPerDiffLenSq[ip] * (densitySCV[k] / density[ip]) * scvf.shape(k);
					//number sumVel = vViscoPerDiffLenSq[ip] *(densitySCV[k])  * scvf.shape(k);
					
					//sumVel += -2.0 * kinVisco[ip] * (densitySCV[k]/density[ip]) * VecProd( RhoGrad[ip], scvf.global_grad(k));
					//sumVel +=  (densitySCV[k]/density[ip]) * VecProd( ViscGrad[ip], scvf.global_grad(k));
					//sumVel +=  (densitySCV[k]/pow(density[ip],2)) * scvf.shape(k) * kinVisco[ip]*VecLengthSq(RhoGrad[ip]);
					//sumVel +=  -(densitySCV[k]/pow(density[ip],2)) * scvf.shape(k) * VecProd(RhoGrad[ip],ViscGrad[ip]);
					
					//    Convective term (no convective terms in the Stokes eq.)
					if (! bStokes)
					{
						//sumVel += densitySCV[k] * vNormStdVelPerConvLen[ip] * (upwind_shape_sh(ip, k) - downwind_shape_sh(ip, k) );
						//sumVel += (densitySCV[k]/density[ip]) * vNormStdVelPerConvLen[ip] * (upwind_shape_sh(ip, k) );
						//sumVel += density[ip] *vNormStdVelPerConvLen[ip] * ( (1.0-Ratio[ip])*upwind_shape_sh(ip, k)  +  Ratio[ip]*(Ratio_rho_up[ip] * upwind_shape_sh(ip, k) + Ratio_rho_do[ip] * downwind_shape_sh(ip, k)));
						//sumVel +=  vNormStdVelPerConvLen[ip] * scvf.shape(k);
						sumVel +=  vNormStdVelPerConvLen[ip] * (vPecletScale[ip] * upwind_shape_sh(ip, k) + (1.0-vPecletScale[ip])*scvf.shape(k));
						
						
						//sumVel += -DenMomentum[ip]*scvf.shape(k) * ;
						
						
						sumVel += vPecletScale[ip] *vNormStdVelPerDownLen[ip] *(downwind_shape_sh(ip, k) - upwind_shape_sh(ip, k));
						
						
						
					}
					for(int d2 = 0; d2 < dim; ++d2)
					{
						if(d2 == d) continue;
						
						sumVel -= density[ip] * vStdVel[ip][d2] * (scvf.global_grad(k))[d2];
						
					}
					
					
					//    Add to rhs
					rhs += sumVel * vCornerValue(d, k);
					rhs_mu += scvf.shape(k) * vCornerValue(d, k);
					
					//    set stab shape
					stab_shape_vel(ip, d, d, k) += sumVel / diag;
					
					for(int d2 = 0; d2 < dim; ++d2)
					{
						if(d2 == d) continue;
						
						const number sumVel2 = density[ip]*vStdVel[ip][d] * (scvf.global_grad(k))[d2];
						
						rhs += sumVel2 * vCornerValue(d2, k);
						
						stab_shape_vel(ip, d, d2, k) = sumVel2 / diag;
					}

					
					if(true)//!multiphase)
					{
						//    Pressure part
						number sumP = -1.0 * scvf.global_grad(k)[d]  ;// - 1.0 * scvf.shape(k) *
						
						//    Add to rhs
						rhs += sumP * vCornerValue(_P_, k);
						
						stab_shape_p(ip, d, k) += sumP / diag;
						
						
						if(Inter->ParticleGradientForce())
							rhs += sumP * ps[k];
						
					}
					else
					{
						if((k == to) || (k == from))
						{
							const number sign = (k==to)? 1.0:-1.0;
							//    Pressure part
							number sumP = -1.0 * Direction[d] * sign  ;// - 1.0 * scvf.shape(k) *
							
							//    Add to rhs
							rhs += sumP * vCornerValue(_P_, k);
							
							stab_shape_p(ip, d, k) += sumP / diag;
							
						}

						
					}
					
					
					//    set stab shape
					
					if( false&&multiphase )
					{
						number sumPJump =0.0;


						/*if ((phase_2[ip] && jump_shape[k]<0) || (!phase_2[ip] && jump_shape[k]>0) )
							
						{
							
							if( true)
							{
								
								sumPJump =  jump_shape[k] * (scvf.global_grad(k)[d] ) ;
								//sumPJump +=  alpha3 * jump_shape[k] * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]  ;
								//sumPJump +=  alpha3 * jump_shape[k] * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) ;
								
							}
							else
							{
								
								sumPJump =  jump_shape[k] * (scvf.global_grad(k)[d] ) ;
								//sumPJump +=  alpha2 * jump_shape[k] * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) ;
							}

							//sumPJump +=  alpha2 * jump_shape[k] * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]  / density[ip];
							//sumPJump +=  alpha2 * jump_shape[k] * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) ;
						}

							
						
						//    Add to rhs
						rhs += sumPJump * pressure_jump_value(k);
						
						//stab_shape_p_jump(ip, d, k) = 0.0*sumPJump/diag;
						//    Add to contributions of the pressure jump at corner k with respect the corner pressure at node SH2
						for(size_t k2 = 0; k2 < scvf.num_sh(); ++k2)
						{
							stab_shape_p(ip, d, k2) += sumPJump * pressure_jump_shape_p(k,k2) / diag;
							for(int d1 = 0; d1 < dim; ++d1)
							{
								stab_shape_vel(ip, d, d1, k2) += sumPJump * pressure_jump_shape_vel( k,  d1, k2) / diag;

							}
						}*/
						
						number sumSlipVel = 0.0;
						
						if (  ((phase_2[ip] && jump_shape[k]<0) || (!phase_2[ip] && jump_shape[k]>0) )   )
						{
							sumSlipVel = -1.0 * vViscoPerDiffLenSq[ip] * scvf.shape(k) * jump_shape[k];
							//sumSlipVel = -1.0 * vViscoPerDiffLenSq[ip] * scvf.shape(k) * (densitySCV[k] / RHO[ip]) * jump_shape[k];
						}
						
						rhs += sumSlipVel * tang_vel(k,d);
						//stab_shape_slip_vel(ip, d, d, k) = 0.0*sumSlipVel/diag;
						
						for(size_t k2 = 0; k2 < scvf.num_sh(); ++k2)
						{
							for(int d1 = 0; d1 < dim; ++d1)
							{
								stab_shape_vel(ip, d, d1, k2) += sumSlipVel * tang_vel_shape_vel( k,  d,  d1,  k2) / diag;

							}
						}
			
					}
					 
										
					
					
					/*if (! bStokes && RelVelSCVF.data_given())
					{
						//    Pressure part
						number SumRelVel =  vNormRelVelPerConvLen[ip] * RelVelSCVF[k][d] * (downwind_shape_sh(ip, k) - upwind_shape_sh(ip, k) );
						//    Add to rhs
						rhs += SumRelVel * vCornerValue(_C_, k);
						
						//    set stab shape
						stab_shape_c(ip, d, k) += SumRelVel / diag;
					}*/
					
					
					
				}
				
				//    Finally, the can invert this row
				stab_vel(ip)[d] = rhs / diag;
			}
		}
	}
}
										

template <>
void NavierStokesFIELDS_2_Stabilization<1>::register_func()
{
	register_func<RegularEdge>();
}

template <>
void NavierStokesFIELDS_2_Stabilization<2>::register_func()
{
	register_func<RegularEdge>();
	register_func<Triangle>();
	register_func<Quadrilateral>();
}

template <>
void NavierStokesFIELDS_2_Stabilization<3>::register_func()
{
	register_func<RegularEdge>();
	register_func<Triangle>();
	register_func<Quadrilateral>();
	register_func<Tetrahedron>();
	register_func<Pyramid>();
	register_func<Prism>();
	register_func<Hexahedron>();
}
/////////////////////////////////////////////////////////////////////////////
// FLOW
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesKARIMIANStabilization<TDim>::
update(const FV1Geometry<TElem, dim>* geo,
       const LocalVector& vCornerValue,
       const MathVector<dim> vStdVel[],
       const bool bStokes,
       const DataImport<number, dim>& kinVisco,
       const DataImport<number, dim>& kinViscoSCV,
       const DataImport<number, dim>& density,
       const DataImport<number, dim>& density_old,
       const DataImport<number, dim>& densitySCV,
	   const DataImport<number, dim>& densitySCV_old,
	   const number ps[],
	   const MathVector<dim> vStdRelVel[],
	   const DataImport<MathVector<dim>, dim>& RelVelSCVF,
       const DataImport<MathVector<dim>, dim>& Source,
       const DataImport<MathVector<dim>, dim>& SourceSCV,
	   const DataImport<MathVector<dim>, dim>& PressGrad,
       const LocalVector* pvCornerValueOldTime, number dt,
	   const int jump_shape[],
	   const bool phase_2[],
	   const bool multiphase)
{
	if( non_zero_shape_ip())
	{
		UG_THROW("Not implemented for ip velocities depending on other ip.");
	}
	
	//    abbreviation for pressure
	static const size_t _P_ = dim;
	//    abbreviation for VolumeFraction
	static const size_t _C_ = dim+1;
	
	//    Some constants
	static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
	static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
	
	for(size_t ip = 0; ip < numIp; ++ip)
	{
		stab_vel(ip) = 0.0;
		conv_vel(ip) = 0.0;
		for(size_t k = 0; k < numSh; ++k)
		{
			for(int d1 = 0; d1 < dim; d1++)
			{
				for(int d2 = 0; d2 < dim; d2++)
				{
					stab_shape_vel(ip, d2, d1, k) = 0.0;
					conv_shape_vel(ip, d2, d1, k) = 0.0;
				}
				stab_shape_p(ip, d1, k) = 0.0;
				stab_shape_c(ip, d1, k) = 0.0;
				
				conv_shape_p(ip, d1, k) = 0.0;
				conv_shape_c(ip, d1, k) = 0.0;
			}
		}
		
	}
	
	//    compute diffusion length
	this->compute_diff_length(*geo);

	MathVector<dim> vStdVel_ip_old[numIp];
	
	
	
	
	//    compute upwind and downwind (no convective terms for the Stokes eq. => no upwind)
	if (! bStokes)
	{
		this->compute_upwind(geo, vStdVel);
		this->compute_downwind(geo, vStdVel);
		
	}
	
	MathVector<dim> RhoGrad[numIp];
	//MathVector<dim> ViscGrad[numIp];
	//MathVector<dim> Vel[numIp];
	//number DenMomentum[numIp];
	number RHO_up[numIp];
	number RHO_do[numIp];
	number Ratio_rho_up[numIp];
	number Ratio_rho_do[numIp];
	number Ratio[numIp];
	number power = 1.0;
	number factor = 1.0;
	if (! bStokes)
	{
		for(size_t ip = 0; ip < numIp; ++ip)
		{
			const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
			const size_t from = scvf.from();
			const size_t to = scvf.to();
			VecSet(RhoGrad[ip], 0.0);
			//VecSet(ViscGrad[ip], 0.0);
			//DenMomentum[ip] = 0.0;
			RHO_up[ip] = 0.0;
			RHO_do[ip] = 0.0;
			Ratio_rho_up[ip] = 0.0;
			//const number Val = +VecTwoNorm(vStdVel[ip]) / (downwind_conv_length(ip) + upwind_conv_length(ip));
			//MU_scvf[ip] = 0.5*(densitySCV[from] * kinViscoSCV[from] + densitySCV[to] * kinViscoSCV[to]);
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				VecScaleAppend(RhoGrad[ip], densitySCV[sh], scvf.global_grad(sh));
				//VecScaleAppend(ViscGrad[ip], kinViscoSCV[sh]*densitySCV[sh], scvf.global_grad(sh));
				//DenMomentum[ip] += Val * (downwind_shape_sh(ip, sh) - upwind_shape_sh(ip, sh)) * densitySCV[sh];

				RHO_up[ip] += upwind_shape_sh(ip, sh) * densitySCV[sh];
				RHO_do[ip] += downwind_shape_sh(ip, sh) * densitySCV[sh];
				
				
			}
			number diff = fabs(RHO_up[ip] - RHO_do[ip]);
			Ratio[ip] = 0.0;//factor*pow(diff /(1.0 + diff), power);
			//if(Ratio[ip] < 0.8)printf("Ratio[%zu] = %f \n", ip,Ratio[ip]);
			
			Ratio_rho_up[ip] = RHO_up[ip] * downwind_conv_length(ip) / ( RHO_up[ip] * downwind_conv_length(ip) + RHO_do[ip] * upwind_conv_length(ip));
			Ratio_rho_do[ip] = 1.0 - Ratio_rho_up[ip];
			/*for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				for(int d = 0; d < dim; d++)
				{
					Vel[ip][d] += vCornerValue(d, sh)*(Ratio_rho_up[ip] * upwind_shape_sh(ip, sh) + Ratio_rho_do[ip] * downwind_shape_sh(ip, sh));
				}
			}*/
			//DenMomentum[ip]=VecProd(RhoGrad[ip],vStdVel[ip]);
		}
	}
	

	bool boolSource = (SourceSCV.data_given()) ? true : false;
	
	//    cache values
	number vViscoPerDiffLenSq[numIp];
	MathVector<dim> SOURCE[numIp];
	
	number vNormStdVelPerConvLen[numIp];
	number vNormStdVelPerDownLen[numIp];
	number vNormRelVelPerConvLen[numIp];
	
	number theta[numIp];
	
	MathVector<dim> vStdVel_stab[numIp];
	
	for(size_t ip = 0; ip < numIp; ++ip)
	{
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
		
		vViscoPerDiffLenSq[ip] = density[ip] * kinVisco[ip] * diff_length_sq_inv(ip);
		if(boolSource) SOURCE[ip] = Source[ip];
		
		VecSet(vStdVel_stab[ip],0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			for(int d = 0; d < dim; d++)
			{
				vStdVel_stab[ip][d] += scvf.shape(sh) * vCornerValue(d, sh);
			}
			  
		}
		
		number dRho = fabs(densitySCV[scvf.to()]- densitySCV[scvf.from()]);
		theta[ip] = pow( dRho/(1.0+dRho), power);
		
		if(!bStokes)
		{
			number Value1 = density[ip] * VecTwoNorm(vStdVel[ip]) / upwind_conv_length(ip);
			
			vNormStdVelPerConvLen[ip] = Value1;
			
			vNormStdVelPerDownLen[ip] = density[ip] * VecTwoNorm(vStdVel[ip]) / (downwind_conv_length(ip) + upwind_conv_length(ip));
			
			
		}


	}
	
	MathVector<dim> vConsGravitySCVF[numIp];
	if(Inter->boolConsistentGravity())
	{
		Inter-> template ConsistentGravitySCVF<TElem>(vConsGravitySCVF, *geo, geo->corners(), numIp, densitySCV.values());
	}
	
	
	//    Find out if upwinded velocities depend on other ip velocities. In that case
	//    we have to solve a matrix system. Else the system is diagonal and we can
	//    compute the inverse directly
	
	//    diagonal case (i.e. upwind vel depend only on corner vel or no upwind)

	//    Loop integration points
	for(size_t ip = 0; ip < numIp; ++ip)
	{
		//    get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
		const size_t from=scvf.from();
		const size_t to=scvf.to();
		
		if (true)
		{

			
			//    First, we compute the contributions to the diagonal
			//    Note: - There is no contribution of the upwind vel to the diagonal
			//            in this case, only for non-diag problems
			//          - The diag does not depend on the dimension
			
			//    the diagonal entry
			number diag, diag2;

			diag = vViscoPerDiffLenSq[ip];
			diag2 = vViscoPerDiffLenSq[ip];
			
			//    Time part
			if(pvCornerValueOldTime != NULL)
				diag += density[ip]/dt;
			
			//    Convective Term  (no convective terms in the Stokes eq.)
			if (! bStokes)
			{
				diag += vNormStdVelPerConvLen[ip];
				
				diag2 += vNormStdVelPerConvLen[ip];
			}
			
			//diag += DenMomentum[ip];
			//     Loop components of velocity
			for(int d = 0; d < dim; d++)
			{
				//    Now, we can assemble the rhs. This rhs is assembled by all
				//    terms, that are non-dependent on the ip vel.
				//    Note, that we can compute the stab_shapes on the fly when setting
				//    up the system.
				
				//    Source
				number rhs = 0.0;
				number rhs_convected = 0.0;
				if(Inter->boolConsistentGravity())
				{
					rhs =  vConsGravitySCVF[ip][d];
				}
				else
				{
					if(boolSource)
					{
						rhs =  SOURCE[ip][d];
					}
				}
				
				/*if(PressGrad.data_given())
				{
					rhs += PressGrad[ip][d];
					
				}*/

				
				//    Time
				if(pvCornerValueOldTime != NULL)
				{
					//	interpolate old time step
					number oldIPVel = 0.0;
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
					//    add to rhs
					rhs += density[ip] * oldIPVel/ dt; //(density_old[ip] / density[ip])
				}
				
				/*if (! bStokes)
				{
					rhs +=  vStdVel[ip][d] * vNormStdVelPerDownLen[ip] * (RHO_do[ip]-RHO_up[ip]);//
					
				}*/
				
				//    loop shape functions
				for(size_t k = 0; k < scvf.num_sh(); ++k)
				{
					//    Diffusion part
					number sumVel = 0.0;
					number sumVel_c = 0.0;
						
					sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k) ;
					sumVel_c = vViscoPerDiffLenSq[ip] * scvf.shape(k) ;
					//sumVel = vViscoPerDiffLenSq[ip] * (densitySCV[k] / density[ip]) * scvf.shape(k);


					//number sumVel = vViscoPerDiffLenSq[ip] * (densitySCV[k] / density[ip]) * scvf.shape(k);
					//number sumVel = vViscoPerDiffLenSq[ip] *(densitySCV[k])  * scvf.shape(k);
					
					//sumVel += -2.0 * kinVisco[ip] * (densitySCV[k]/density[ip]) * VecProd( RhoGrad[ip], scvf.global_grad(k));
					//sumVel +=  (densitySCV[k]/density[ip]) * VecProd( ViscGrad[ip], scvf.global_grad(k));
					//sumVel +=  (densitySCV[k]/pow(density[ip],2)) * scvf.shape(k) * kinVisco[ip]*VecLengthSq(RhoGrad[ip]);
					//sumVel +=  -(densitySCV[k]/pow(density[ip],2)) * scvf.shape(k) * VecProd(RhoGrad[ip],ViscGrad[ip]);
					
					//    Convective term (no convective terms in the Stokes eq.)
					if (! bStokes)
					{
						//sumVel += densitySCV[k] * vNormStdVelPerConvLen[ip] * (upwind_shape_sh(ip, k) - downwind_shape_sh(ip, k) );
						//sumVel += (densitySCV[k]/density[ip]) * vNormStdVelPerConvLen[ip] * (upwind_shape_sh(ip, k) );
						//sumVel += density[ip] *vNormStdVelPerConvLen[ip] * ( (1.0-Ratio[ip])*upwind_shape_sh(ip, k)  +  Ratio[ip]*(Ratio_rho_up[ip] * upwind_shape_sh(ip, k) + Ratio_rho_do[ip] * downwind_shape_sh(ip, k)));
						//sumVel +=  vNormStdVelPerConvLen[ip] * scvf.shape(k);
						sumVel +=  vNormStdVelPerConvLen[ip] * ((1.0-theta[ip]) * upwind_shape_sh(ip, k) + theta[ip]*scvf.shape(k));
						sumVel += vNormStdVelPerDownLen[ip] *(1.0-theta[ip])*(downwind_shape_sh(ip, k) - upwind_shape_sh(ip, k));
						
						sumVel_c += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k) ;
						
						
						
					}
					for(int d2 = 0; d2 < dim; ++d2)
					{
						if(d2 == d) continue;
						
						sumVel -= density[ip] * vStdVel[ip][d2] * (scvf.global_grad(k))[d2];
						
					}
					
					
					//    Add to rhs
					rhs += sumVel * vCornerValue(d, k);
					rhs_convected += sumVel_c * vCornerValue(d, k);
					
					//    set stab shape
					stab_shape_vel(ip, d, d, k) += sumVel / diag;
					conv_shape_vel(ip, d, d, k) += sumVel_c / diag2;
					
					for(int d2 = 0; d2 < dim; ++d2)
					{
						if(d2 == d) continue;
						
						const number sumVel2 = density[ip]*vStdVel[ip][d] * (scvf.global_grad(k))[d2];
						
						rhs += sumVel2 * vCornerValue(d2, k);
						
						stab_shape_vel(ip, d, d2, k) = sumVel2 / diag;
					}

					
					//    Pressure part
					number sumP = -1.0 * scvf.global_grad(k)[d]  ;// - 1.0 * scvf.shape(k) *
					
					//    Add to rhs
					rhs += sumP * vCornerValue(_P_, k);
					rhs_convected += sumP * vCornerValue(_P_, k);
					
					stab_shape_p(ip, d, k) += sumP / diag;
					conv_shape_p(ip, d, k) += sumP / diag2;
					
					
					if(Inter->ParticleGradientForce())
						rhs += sumP * ps[k];
						
				}
				
				//    Finally, the can invert this row
				stab_vel(ip)[d] = rhs / diag;
				conv_vel(ip)[d] = rhs_convected / diag2;
			}
		}
	}
}


template <>
void NavierStokesKARIMIANStabilization<1>::register_func()
{
    register_func<RegularEdge>();
}

template <>
void NavierStokesKARIMIANStabilization<2>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
}

template <>
void NavierStokesKARIMIANStabilization<3>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
    register_func<Tetrahedron>();
    register_func<Pyramid>();
    register_func<Prism>();
    register_func<Hexahedron>();
}

/////////////////////////////////////////////////////////////////////////////
// NO
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesNOStabilization<TDim>::
update(const FV1Geometry<TElem, dim>* geo,
       const LocalVector& vCornerValue,
       const MathVector<dim> vStdVel[],
       const bool bStokes,
       const DataImport<number, dim>& kinVisco,
       const DataImport<number, dim>& kinViscoSCV,
       const DataImport<number, dim>& density,
       const DataImport<number, dim>& density_old,
       const DataImport<number, dim>& densitySCV,
	   const DataImport<number, dim>& densitySCV_old,
	   const number ps[],
	   const MathVector<dim> vStdRelVel[],
	   const DataImport<MathVector<dim>, dim>& RelVelSCVF,
       const DataImport<MathVector<dim>, dim>& Source,
       const DataImport<MathVector<dim>, dim>& SourceSCV,
	   const DataImport<MathVector<dim>, dim>& PressGrad,
       const LocalVector* pvCornerValueOldTime, number dt,
	   const int jump_shape[],
	   const bool phase_2[],
	   const bool multiphase)
{
//    abbreviation for pressure
    static const size_t _P_ = dim;

//    Some constants
    static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;

//    We can solve the systems ip by ip
    for(size_t ip = 0; ip < numIp; ++ip)
    {
    //    get SubControlVolumeFace
        const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);



    //     Loop components of velocity
        for(size_t d = 0; d < (size_t)dim; d++)
        {
            number rhs=0.0;
        //    loop shape functions
            for(size_t k = 0; k < scvf.num_sh(); ++k)
            {
            //    Diffusion part
                number sumVel = scvf.shape(k);

            //    Add to rhs
                rhs += sumVel * vCornerValue(d, k);

            //    set stab shape
                stab_shape_vel(ip, d, d, k) = sumVel;

            //    Pressure part
                const number sumP = -1.0 * 1 * (scvf.global_grad(k))[d] / density[ip];

            //    Add to rhs
                rhs += sumP * vCornerValue(_P_, k);

            //    set stab shape
                stab_shape_p(ip, d, k) = sumP;
            }

        //    Finally, the can invert this row
            stab_vel(ip)[d] = rhs;
        }
    }
}

template <>
void NavierStokesNOStabilization<1>::register_func()
{
    register_func<RegularEdge>();
}

template <>
void NavierStokesNOStabilization<2>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
}

template <>
void NavierStokesNOStabilization<3>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
    register_func<Tetrahedron>();
    register_func<Pyramid>();
    register_func<Prism>();
    register_func<Hexahedron>();
}


/////////////////////////////////////////////////////////////////////////////
// NO STABILIZATION (Note: The discretization is then unstable!)
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesFV1WithoutStabilization<TDim>::
update(const FV1Geometry<TElem, dim>* geo,
       const LocalVector& vCornerValue,
       const MathVector<dim> vStdVel[],
       const bool bStokes,
       const DataImport<number, dim>& kinVisco,
       const DataImport<number, dim>& kinViscoSCV,
       const DataImport<number, dim>& density,
       const DataImport<number, dim>& density_old,
       const DataImport<number, dim>& densitySCV,
	   const DataImport<number, dim>& densitySCV_old,
	   const number ps[],
	   const MathVector<dim> vStdRelVel[],
	   const DataImport<MathVector<dim>, dim>& RelVelSCVF,
       const DataImport<MathVector<dim>, dim>& Source,
       const DataImport<MathVector<dim>, dim>& SourceSCV,
	   const DataImport<MathVector<dim>, dim>& PressGrad,
       const LocalVector* pvCornerValueOldTime, number dt,
	   const int jump_shape[],
	   const bool phase_2[],
	   const bool multiphase)
{
//	Some constants
	static const size_t numIP = FV1Geometry<TElem, dim>::numSCVF;
	static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;

//	Compute upwind (no convective terms for the Stokes eq. => no upwind)
	if (! bStokes) this->compute_upwind(geo, vStdVel);

//	No dependence on the pressure:
	for (size_t ip = 0; ip < numIP; ip++)
		for (size_t i = 0; i < dim; i++)
			for (size_t sh = 0; sh < numSh; sh++)
				stab_shape_p (ip, i, sh) = 0;
	
//	The velocities are interpolated according to the FE shape functions:
	for (size_t ip = 0; ip < numIP; ip++)
	{
	//	get the shapes
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
		for (size_t sh = 0; sh < numSh; sh++)
		{
			number val = scvf.shape (sh);
			for (size_t i = 0; i < dim; i++)
			{
				for (size_t j = 0; j < dim; j++)
					stab_shape_vel (ip, i, j, sh) = 0;
				stab_shape_vel (ip, i, i, sh) = val;
			}
		}
		
	//	the interpolation is given as the argument
		stab_vel (ip) = vStdVel [ip];
	}
}

template <>
void NavierStokesFV1WithoutStabilization<1>::register_func()
{
	register_func<RegularEdge>();
}

template <>
void NavierStokesFV1WithoutStabilization<2>::register_func()
{
	register_func<RegularEdge>();
	register_func<Triangle>();
	register_func<Quadrilateral>();
}

template <>
void NavierStokesFV1WithoutStabilization<3>::register_func()
{
	register_func<RegularEdge>();
	register_func<Triangle>();
	register_func<Quadrilateral>();
	register_func<Tetrahedron>();
	register_func<Pyramid>();
	register_func<Prism>();
	register_func<Hexahedron>();
}

////////////////////////////////////////////////////////////////////////////////
//	explicit instantiations
////////////////////////////////////////////////////////////////////////////////

/*#ifdef UG_DIM_1
template class INavierStokesFV1Stabilization<1>;
template class NavierStokesFIELDSStabilization<1>;
template class NavierStokesFLOWStabilization<1>;

template SmartPtr<INavierStokesFV1Stabilization<1> >CreateNavierStokesStabilization<1>(const std::string& name);
#endif*/
#ifdef UG_DIM_2
template class INavierStokesFV1Stabilization<2>;
template class INavierStokesSRFV1Stabilization<2>;
template class NavierStokesFIELDSStabilization<2>;
template class NavierStokesFLOWStabilization<2>;
template class NavierStokesFIELDS_2_Stabilization<2>;
template class NavierStokesFLOW_2_Stabilization<2>;
template class NavierStokesKARIMIANStabilization<2>;

template SmartPtr<INavierStokesSRFV1Stabilization<2> >CreateNavierStokesStabilization<2>(const std::string& name);
#endif
#ifdef UG_DIM_3
template class INavierStokesFV1Stabilization<3>;
template class INavierStokesSRFV1Stabilization<3>;
template class NavierStokesFIELDSStabilization<3>;
template class NavierStokesFLOWStabilization<3>;
template class NavierStokesFIELDS_2_Stabilization<3>;
template class NavierStokesFLOW_2_Stabilization<3>;
template class NavierStokesKARIMIANStabilization<3>;

template SmartPtr<INavierStokesSRFV1Stabilization<3> >CreateNavierStokesStabilization<3>(const std::string& name);
#endif

} // namespace NavierStokes
} // end namespace ug
