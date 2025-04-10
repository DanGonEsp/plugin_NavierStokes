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
    if(n == "viscosity") return SmartPtr<NavierStokesVISCOSITY_Stabilization<dim> >(new NavierStokesVISCOSITY_Stabilization<dim>());
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
       const DataImport<number, dim>& densitySCV,
       const number pressure_jump[],
       const MathVector<dim> SlipVel[],
       const number jump_shape[],
       const MathVector<dim> normal[],
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt,
       const number density_ref,
       const bool multiphase,
       const bool phase_2[],
       const number theta,
       number** SCVFinterShape)
{
    if(multiphase)
        UG_THROW("Pressure Jump Not implemented for FIELDS stabilization.");
    
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
                if(pSource != NULL)
                    rhs = (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];
                
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
                if(pSource != NULL)
                    f[ip] = (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];
                
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
       const DataImport<number, dim>& densitySCV,
       const number pressure_jump[],
       const MathVector<dim> SlipVel[],
       const number jump_shape[],
       const MathVector<dim> normal[],
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt,
       const number density_ref,
       const bool multiphase,
       const bool phase_2[],
       const number theta,
       number** SCVFinterShape)
{
    if(multiphase)
        UG_THROW("Pressure Jump Not implemented for FLOW stabilization.");
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
                if(pSource != NULL)
                    rhs = (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];
                
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
                if(pSource != NULL)
                    f[ip] += (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];
                
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
       const DataImport<number, dim>& densitySCV,
       const number pressure_jump[],
       const MathVector<dim> SlipVel[],
       const number jump_shape[],
       const MathVector<dim> normal[],
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt,
       const number density_ref,
       const bool multiphase,
       const bool phase_2[],
       const number theta,
       number** SCVFinterShape)
{
    if( non_zero_shape_ip())
    {
        UG_THROW("Not implemented for ip velocities depending on other ip.");
    }
    
    
    
    //    abbreviation for pressure
    static const size_t _P_ = dim;
    
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
                    stab_shape_slip_vel(ip, d2, d1, k) = 0.0;
                }
                stab_shape_p(ip, d1, k) = 0.0;
                stab_shape_p_jump(ip, d1, k) =  0.0;
                
            }
        }
    }
    
    //    compute upwind and downwind (no convective terms for the Stokes eq. => no upwind)
    if (! bStokes)
    {
        this->compute_upwind(geo, vStdVel);
        //this->compute_downwind(geo, vStdVel);
    }
    
    //    compute diffusion length
    this->compute_diff_length(*geo);
    
    MathVector<dim> DenGrad; VecSet(DenGrad,0.0);
    number DenMomentum[numIp];
    for(size_t sh = 0; sh < numSh; ++sh)
    {
        //     get current SCV
        const typename FV1Geometry<TElem, dim>::SCV& scv = geo->scv(sh);
        VecScaleAppend(DenGrad, densitySCV[sh], scv.global_grad(sh));
    }
    for(size_t ip = 0; ip < numIp; ++ip)
    {
        DenMomentum[ip]=VecProd(DenGrad,vStdVel[ip])/density[ip];
        
    }

    number mu_2 = 0.0, mu_1 = 0.0;
    number rho_2 = 0.0, rho_1 = 0.0;
    
    number alpha1 = ((!multiphase && jump_shape[0]>0))? 1.0 : 1.0;
    number alpha2 = ((!multiphase && jump_shape[0]>0))? 1.0 : 1.0;
    number alpha3 = 1.0;
    
    //printf("Bug1 \n");
    
    //number fase1 = 0.0;
    
    if (multiphase)
    {

        
        for(size_t sh = 0; sh < numSh; ++sh)
        {
            if (jump_shape[sh]>0)
            {
                mu_2 = fmax( mu_2, densitySCV[sh] * kinViscoSCV[sh]);
                rho_2 = fmax( rho_2, densitySCV[sh] );
            }
            else
            {
                mu_1 = fmax( mu_1, densitySCV[sh] * kinViscoSCV[sh]);
                rho_1 = fmax( rho_1, densitySCV[sh] );
            }
            
            //fase1 +=jump_shape[sh];
        }
        
        
        if ((mu_2 < mu_1)||(mu_2<0.0)||(mu_1<0.0))
            UG_THROW("Viscosity in phase 1 is lower that phase 2");
        if ((rho_2 < rho_1)||(rho_2<0.0)||(rho_1<0.0))
            UG_THROW("Density in phase 1 is lower that phase 2");
    }
    

    //    cache values
    number vViscoPerDiffLenSq[numIp];
    number RHO[numIp];
    bool interface_change[numIp];
    
    number vNormStdVelPerConvLen[numIp];
    number w_pe[numIp];
    
    for(size_t ip = 0; ip < numIp; ++ip)
    {
        const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
        if (multiphase)
        {
            
            if(jump_shape[scvf.to()]*jump_shape[scvf.from()]<0.0)
            {
                interface_change[ip] = true;
            }
                
            else
            {
                interface_change[ip] = false;

            }
            
            if(phase_2[ip])
            {
                vViscoPerDiffLenSq[ip] = (mu_2/rho_2) * diff_length_sq_inv(ip);
                RHO[ip]=rho_2;
            }
            else
            {
                vViscoPerDiffLenSq[ip] = (mu_1/rho_1) * diff_length_sq_inv(ip);
                RHO[ip]=rho_1;
            }

            
        }
        else
        {
            vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);
            RHO[ip]=density[ip];
        }
        
        if(!bStokes)
        {
            const number norm = VecTwoNorm(vStdVel[ip]);
            vNormStdVelPerConvLen[ip] = norm / upwind_conv_length(ip);
            
            
        //    compute peclet number
            number Pe = VecProd(vStdVel[ip], scvf.normal())/VecTwoNormSq(scvf.normal())
             * VecDistance(geo->corners() [scvf.to()], geo->corners() [scvf.from()]) / kinVisco[ip];

        //    compute weight
            const number Pe2 = Pe * Pe;
            w_pe[ip] = Pe2 / (5.0 + Pe2);
            
        }


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
        
        if (true)
        {

            
            //    First, we compute the contributions to the diagonal
            //    Note: - There is no contribution of the upwind vel to the diagonal
            //            in this case, only for non-diag problems
            //          - The diag does not depend on the dimension
            
            //    the diagonal entry
            number diag = vViscoPerDiffLenSq[ip];
            
            //    Time part
            /*if(pvCornerValueOldTime != NULL)
                diag += 1./dt;*/
            
            //    Convective Term  (no convective terms in the Stokes eq.)
            if (! bStokes)
                diag += vNormStdVelPerConvLen[ip];
            
            
            //     Loop components of velocity
            for(int d = 0; d < dim; d++)
            {
                //    Now, we can assemble the rhs. This rhs is assembled by all
                //    terms, that are non-dependent on the ip vel.
                //    Note, that we can compute the stab_shapes on the fly when setting
                //    up the system.
                
                //    Source
                number rhs = 0.0;
                if(pSource != NULL)
                    rhs = (  ( RHO[ip]-density_ref) / RHO[ip]) * (*pSource)[ip][d];
                
                //    Time
                /*if(pvCornerValueOldTime != NULL)
                {
                    //    interpolate old time step
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                    
                    //    add to rhs
                    rhs += oldIPVel / dt;
                }*/
                
                //    loop shape functions
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                    //    Diffusion part
                    number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                    
                    //    Convective term (no convective terms in the Stokes eq.)
                    if (! bStokes)
                    {
                        sumVel += vNormStdVelPerConvLen[ip] * (w_pe[ip]*upwind_shape_sh(ip, k)+(1.0-w_pe[ip])*scvf.shape(k));
                        sumVel +=- DenMomentum[ip]*scvf.shape(k)
                    }
                    
                    
                    //    Add to rhs
                    rhs += sumVel * vCornerValue(d, k);
                    
                    //    set stab shape
                    stab_shape_vel(ip, d, d, k) = sumVel / diag;
                    
                    
                    //    Pressure part
                    number sumP = 0.0;
                    if(multiphase)
                    {
                        if( (interface_change[ip]))
                        {
                            sumP = -1.0 * alpha3 * (scvf.global_grad(k))[d]  / RHO[ip];
                            //sumP += -1.0 * alpha3 * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) / RHO[ip];
                            //sumP += -1.0 * alpha3 * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]  / RHO[ip];
                            
                        }
                        else
                        {
                            sumP = -1.0 * alpha2 * (scvf.global_grad(k))[d]  / RHO[ip];
                            //sumP += -1.0 * alpha2 * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) / RHO[ip];
                        }
                        
                        //sumP += -1.0 * alpha3 * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]  / rho_2;
                        //sumP += -1.0 * alpha2 * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) / rho_2;
                        
                    }
                    else
                    {
                        sumP = -1.0 * alpha1 * (scvf.global_grad(k))[d]  / density[ip];
                    }
                    
                    
                    
                    //    Add to rhs
                    rhs += sumP * vCornerValue(_P_, k);
                    
                    //    set stab shape
                     
                    /*if ( multiphase&&((jump_shape[scvf.from()]>0 && jump_shape[k]<0) || (jump_shape[scvf.from()]<0 && jump_shape[k]>0)))
                        sumP = 0.0;*/
                        
                    stab_shape_p(ip, d, k) = sumP / diag;
                    
                    number sumPJump =0.0;
                    number sumSlipVel = 0.0;
                    
                    if(multiphase)
                    {
                        
                        

                        if ((phase_2[ip] && jump_shape[k]<0) || (!phase_2[ip] && jump_shape[k]>0) )
                            
                        {
                            
                            if( interface_change[ip])
                            {
                                
                                sumPJump =  alpha3 * jump_shape[k] * (scvf.global_grad(k)[d] ) / RHO[ip];
                                //sumPJump +=  alpha3 * jump_shape[k] * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]  / RHO[ip];
                                //sumPJump +=  alpha3 * jump_shape[k] * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) / RHO[ip];
                                
                            }
                            else
                            {
                                
                                sumPJump =  alpha2 * jump_shape[k] * (scvf.global_grad(k)[d] ) / RHO[ip];
                                //sumPJump +=  alpha2 * jump_shape[k] * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) / RHO[ip];
                            }

                            //sumPJump +=  alpha2 * jump_shape[k] * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]  / density[ip];
                            //sumPJump +=  alpha2 * jump_shape[k] * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) / rho_1;
                        }

                            
                        
                        //    Add to rhs
                        rhs += sumPJump *(pressure_jump[k]);
                        
                
                        if (  ((phase_2[ip] && jump_shape[k]<0) || (!phase_2[ip] && jump_shape[k]>0) )   )
                            sumSlipVel = -1.0 * vViscoPerDiffLenSq[ip] * scvf.shape(k) * jump_shape[k];
                        
                        rhs += sumSlipVel * SlipVel[k][d];
            
                    }

                    
                    stab_shape_p_jump(ip, d, k) = sumPJump/diag;
                        
                    stab_shape_slip_vel(ip, d, d, k) = sumSlipVel/diag;
                        

                    
                    
                }
                
                //    Finally, the can invert this row
                stab_vel(ip)[d] = rhs / diag;
                //printf("Bug7 \n");
                
            }
            if(multiphase && (interface_change[ip]))
            {

                
                
                
                MathVector<dim> DX = normal[0];
                number Factor1 = 1.0;//Magnitud;//Vec2/Vec1;

                
                VecScale(stab_vel(ip)  ,  DX,  Factor1 * VecProd(DX,stab_vel(ip)));

                
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                    
                    
                    number ValueV = 0.0;
                    number ValueP = 0.0;
                    number ValueJumpP = 0.0;
                    number ValueSlipV = 0.0;
                    
                    
                    for(int d1 = 0; d1 < dim; d1++)
                    {
                        
                        ValueP += stab_shape_p(ip, d1, k) * DX[d1];
                        ValueJumpP += stab_shape_p_jump(ip, d1, k) * DX[d1];
                        
                        
                        
                        ValueV = 0.0;
                        ValueSlipV = 0.0;
                        
                        
                        for(int d2 = 0; d2 < dim; d2++)
                        {
                            ValueV       += stab_shape_vel(ip, d2, d1, k)      * DX[d2];
                            ValueSlipV   += stab_shape_slip_vel(ip, d2, d1, k) * DX[d2];
                            
                        }
                        for(int d2 = 0; d2 < dim; d2++)
                        {
                            
                            stab_shape_vel(ip, d2, d1, k) = Factor1 * ValueV * DX[d2];
                            stab_shape_slip_vel(ip, d2, d1, k) = Factor1 * ValueSlipV * DX[d2];
                            
                            
                            //stab_shape_vel(ip, d2, d1, k) = Factor1 * stab_shape_vel(ip, d2, d1, k);
                            //stab_shape_slip_vel(ip, d2, d1, k) = Factor1 * stab_shape_slip_vel(ip, d2, d1, k);
                            
                            

                            
                        }
                        
                    }
                    
                    for(int d1 = 0; d1 < dim; d1++)
                    {

                        stab_shape_p(ip, d1, k) = Factor1 * ValueP * DX[d1];
                        stab_shape_p_jump(ip, d1, k) =  Factor1 * ValueJumpP * DX[d1];
                        
                        
                        //stab_shape_p(ip, d1, k) = Factor1 * stab_shape_p(ip, d1, k);
                        //stab_shape_p_jump(ip, d1, k) =  Factor1 * stab_shape_p_jump(ip, d1, k);
                        
                        
                    }

                    
                }
                
                
            }
            
        }
    }
    /*
     number sumVel2;
     if ( d == d2)
         sumVel2 = scvf.shape(k);
     else
         sumVel2 = 0.0;
     //sumVel2 = scvf.shape(k) * normal[k][d2] *normal[k][d];
     */

}
/*{
    if(multiphase)
        UG_THROW("Slip Vel Jump Not implemented for fields_2 stabilization.");
    if( non_zero_shape_ip())
    {
        UG_THROW("Not implemented for ip velocities depending on other ip.");
    }
    
    //    abbreviation for pressure
    static const size_t _P_ = dim;
    
    //    Some constants
    static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
    static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
    
    //    compute upwind (no convective terms for the Stokes eq. => no upwind)
    if (! bStokes) this->compute_upwind(geo, vStdVel);
    
    //    compute diffusion length
    this->compute_diff_length(*geo);
    
    if(!multiphase)
    {
        
        //    cache values
        number vViscoPerDiffLenSq[numIp];
        for(size_t ip = 0; ip < numIp; ++ip)
        {
            vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);
        }
        
        number vNormStdVelPerConvLen[numIp];
        if(!bStokes)
            for(size_t ip = 0; ip < numIp; ++ip)
                vNormStdVelPerConvLen[ip] = VecTwoNorm(vStdVel[ip]) / upwind_conv_length(ip);
        
        //    Find out if upwinded velocities depend on other ip velocities. In that case
        //    we have to solve a matrix system. Else the system is diagonal and we can
        //    compute the inverse directly
        
        //    diagonal case (i.e. upwind vel depend only on corner vel or no upwind)
        
        //    We can solve the systems ip by ip
        for(size_t ip = 0; ip < numIp; ++ip)
        {
            //    get SubControlVolumeFace
            const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
            
            //    First, we compute the contributions to the diagonal
            //    Note: - There is no contribution of the upwind vel to the diagonal
            //            in this case, only for non-diag problems
            //          - The diag does not depend on the dimension
            
            //    Diffusion part
            number diag = vViscoPerDiffLenSq[ip];
            
            //    Time part
            if(pvCornerValueOldTime != NULL)
            {
                    diag += 1./dt;
            }
            
            
            //    Convective Term (no convective terms in the Stokes eq.)
            if (! bStokes)
                diag += vNormStdVelPerConvLen[ip];
            
            //     Loop components of velocity
            for(size_t d = 0; d < (size_t)dim; d++)
            {
                //    Now, we can assemble the rhs. This rhs is assembled by all
                //    terms, that are non-dependent on the ip vel.
                //    Note, that we can compute the stab_shapes on the fly when setting
                //    up the system.
                
                //    Source
                number rhs = 0.0;
                if(pSource != NULL)
                    rhs = (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];
                
                //    Time
                if(pvCornerValueOldTime != NULL)
                {
                    //    interpolate old time step
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                    
                    //    add to rhs
                    rhs += oldIPVel / dt;
                }
                
                //    loop shape functions
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                    //    Diffusion part
                    number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                    
                    //    Convective term
                    if (! bStokes) // no convective terms in the Stokes eq.
                        sumVel += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
                    
                    //    Add to rhs
                    rhs += sumVel * vCornerValue(d, k);
                    
                    //    set stab shape
                    stab_shape_vel(ip, d, d, k) = sumVel / diag;
                    for(size_t d2 = 0; d2 < dim; ++d2)
                    {
                        if(d2 == d) continue;
                        stab_shape_vel(ip, d, d2, k) = 0.0;
                    }
                    
                    
                    //    Pressure part
                    const number sumP = -1.0 * (scvf.global_grad(k))[d] / density[ip];
                    
                    //    Add to rhs
                    rhs += sumP * vCornerValue(_P_, k);
                    
                    //    set stab shape
                    stab_shape_p(ip, d, k) = sumP / diag;
                    
                    
                }
                
                //    Finally, the can invert this row
                stab_vel(ip)[d] = rhs / diag;
            }
        }
        // end switch for non-diag
    }
    else
    {
                
        number mu_2 = 0.0, mu_1 = 0.0;
        number rho_2 = 0.0, rho_1 = 0.0;
        
        
        number visc1 = -10000000, visc2 = -10000000;
        number rho1 = -10000000, rho2 = -10000000;
        
        for(size_t sh = 0; sh < numSh; ++sh)
        {
            if (jump_shape[sh]>0)
            {
                visc2 = fmax( visc2, densitySCV[sh] * kinViscoSCV[sh]);
                rho2 = fmax( rho2, densitySCV[sh] );
            }
            else
            {
                visc1 = fmax( visc1, -densitySCV[sh] * kinViscoSCV[sh]);
                rho1 = fmax( rho1, -densitySCV[sh] );
            }
            
        }
        
        mu_2 = visc2;
        mu_1 = -visc1;
        rho_2 = rho2;
        rho_1 = -rho1;
        
        if ((mu_2 < mu_1)||(mu_2<0.0)||(mu_1<0.0))
            UG_THROW("Viscosity in phase 1 is lower that phase 2");
        if ((rho_2 < rho_1)||(rho_2<0.0)||(rho_1<0.0))
            UG_THROW("Density in phase 1 is lower that phase 2");
        
        
        //    cache values
        number vViscoPerDiffLenSq[numIp];
        number vNormStdVelPerConvLen[numIp];
        bool interface_change[numIp];
        for(size_t ip = 0; ip < numIp; ++ip)
        {
            const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
            if(jump_shape[scvf.to()]*jump_shape[scvf.from()]<0.0)
            {
                interface_change[ip] = true;
                vViscoPerDiffLenSq[ip] = 1.0;
            }
            else
            {
                vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);
                interface_change[ip] = false;
            }
            if(!bStokes)
                vNormStdVelPerConvLen[ip] = VecTwoNorm(vStdVel[ip]) / upwind_conv_length(ip);
            
        }
        
        
        
        //    We can solve the systems ip by ip
        for(size_t ip = 0; ip < numIp; ++ip)
        {
            //    get SubControlVolumeFace
            const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
            
            //    First, we compute the contributions to the diagonal
            //    Note: - There is no contribution of the upwind vel to the diagonal
            //            in this case, only for non-diag problems
            //          - The diag does not depend on the dimension
            if(!interface_change[ip])
            {
                //    Diffusion part
                    number diag = vViscoPerDiffLenSq[ip];
                
                //    Time part
                if(pvCornerValueOldTime != NULL)
                {
                    diag += 1./dt;
                }
                
                
                //    Convective Term (no convective terms in the Stokes eq.)
                if (! bStokes)
                    diag += vNormStdVelPerConvLen[ip];
                
                //     Loop components of velocity
                for(size_t d = 0; d < (size_t)dim; d++)
                {
                    //    Now, we can assemble the rhs. This rhs is assembled by all
                    //    terms, that are non-dependent on the ip vel.
                    //    Note, that we can compute the stab_shapes on the fly when setting
                    //    up the system.
                    
                    //    Source
                    number rhs = 0.0;
                    if(pSource != NULL)
                    {

                        rhs = (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];
                    }
                    
                    //    Time
                    if(pvCornerValueOldTime != NULL)
                    {
                        //    interpolate old time step
                        number oldIPVel = 0.0;
                        for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                            oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                        
                        //    add to rhs
                        rhs += oldIPVel / dt;
                    }
            
                    //    loop shape functions
                    for(size_t k = 0; k < scvf.num_sh(); ++k)
                    {
                        //    Diffusion part
                        number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                        
                        //    Convective term
                        if (! bStokes ) // no convective terms in the Stokes eq.
                            sumVel += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
                        
                        //    Add to rhs
                        rhs += sumVel * vCornerValue(d, k);
                        
                        //    set stab shape
                        stab_shape_vel(ip, d, d, k) = sumVel / diag;
                        
                        for(size_t d2 = 0; d2 < dim; ++d2)
                        {
                            if(d2 == d) continue;
                            stab_shape_vel(ip, d, d2, k) = 0.0;
                        }
                        
                        
                        
                        
                        
                        number sumP;
                        
                        number alpha = 1.0;
                        
                        if (jump_shape[scvf.from()] * jump_shape[k]<0.0)
                            alpha = 1.0;
                        
                        sumP = -1.0 * alpha * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d])/ density[ip];
                            

                        
                        //    Add to rhs
                        rhs += sumP * vCornerValue(_P_, k);
                        
                        //    set stab shape
                        stab_shape_p(ip, d, k) = sumP / diag;
                        
                        
                        number sumPJump =0.0;
                         
                        if (jump_shape[scvf.from()]>0.0)
                        {
                            
                            if (jump_shape[k]<0.0)
                                sumPJump =  jump_shape[k] * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) / density[ip];
                            else
                                sumPJump = 0.0;
                            
                        }
                        else
                        {
                            if (jump_shape[k]>0.0)
                                sumPJump =  jump_shape[k] * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) / density[ip];
                            else
                                sumPJump = 0.0;
                        }
                        
                        sumPJump *= alpha;
                        //    Add to rhs
                        rhs += sumPJump * (pressure_jump[k]);
                        stab_shape_p_jump(ip, d, k) = sumPJump/diag;
                        
                        
                    }
                    
                    //    Finally, the can invert this row
                    stab_vel(ip)[d] = rhs / diag;
                }
            }
            else
            {

                            
                //     Loop components of velocity
                for(size_t d = 0; d < (size_t)dim; d++)
                {
                    //    Now, we can assemble the rhs. This rhs is assembled by all
                    //    terms, that are non-dependent on the ip vel.
                    //    Note, that we can compute the stab_shapes on the fly when setting
                    //    up the system.
                    
                    //    Source
                    number rhs = 0.0;

                    
                    //    loop shape functions
                    for(size_t k = 0; k < scvf.num_sh(); ++k)
                    {
                        
                        for(size_t d2 = 0; d2 < dim; ++d2)
                        {
                            number sumVel2;
                            if ( d == d2)
                                sumVel2 = scvf.shape(k);
                            else
                                sumVel2 = 0.0;
                            sumVel2 = scvf.shape(k) * normal[k][d2] *normal[k][d];

                            
                            rhs += sumVel2 * vCornerValue(d2, k);
                            
                            
                            stab_shape_vel(ip, d, d2, k) = sumVel2;
                        }
                        
                        
                        number sumP;
                        
                        //    Pressure part

                        sumP = 0.0 ;
                        
                        //    set stab shape
                        stab_shape_p(ip, d, k) = sumP ;
                        
                        
                        number sumPJump =0.0;

                        //    Add to rhs

                        stab_shape_p_jump(ip, d, k) = sumPJump;
                        
                    }
                    
                    //    Finally, the can invert this row
                    stab_vel(ip)[d] = rhs ;
                }
            }
        }
        
    }
}*/
/*
 
 {
 //    abbreviation for pressure
     static const size_t _P_ = dim;

 //    Some constants
     static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
     static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;

 //    compute upwind (no convective terms for the Stokes eq. => no upwind)
     if (! bStokes) this->compute_upwind(geo, vStdVel);

 //    compute diffusion length
     this->compute_diff_length(*geo);
     
     number mu_2 = 0.0, mu_1 = 0.0;
     number rho_2 = 0.0, rho_1 = 0.0;
     
     
     if(multiphase)
     {
         
         number visc1 = -10000000, visc2 = -10000000;
         number rho1 = -10000000, rho2 = -10000000;
         
         for(size_t sh = 0; sh < numSh; ++sh)
         {
             if (jump_shape[sh]>0)
             {
                 visc2 = fmax( visc2, densitySCV[sh] * kinViscoSCV[sh]);
                 rho2 = fmax( rho2, densitySCV[sh] );
             }
             else
             {
                 visc1 = fmax( visc1, -densitySCV[sh] * kinViscoSCV[sh]);
                 rho1 = fmax( rho1, -densitySCV[sh] );
             }
             
         }
         
         mu_2 = visc2;
         mu_1 = -visc1;
         rho_2 = rho2;
         rho_1 = -rho1;
         
         if ((mu_2 < mu_1)||(mu_2<0.0)||(mu_1<0.0))
             UG_THROW("Viscosity in phase 1 is lower that phase 2");
         if ((rho_2 < rho_1)||(rho_2<0.0)||(rho_1<0.0))
             UG_THROW("Density in phase 1 is lower that phase 2");
     }

 //    cache values
     number vViscoPerDiffLenSq[numIp];
     bool interface_change[numIp];
     for(size_t ip = 0; ip < numIp; ++ip)
     {
         if(!multiphase)
             vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);
         else
         {
             const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
             if(jump_shape[scvf.to()]*jump_shape[scvf.from()]<0.0)
             {
                 interface_change[ip] = true;
                 vViscoPerDiffLenSq[ip] = (mu_2/rho_2- mu_1/rho_1) * diff_length_sq_inv(ip);
             }
             else
             {
                 if(jump_shape[scvf.to()]>0.0)
                     vViscoPerDiffLenSq[ip] = (mu_2/rho_2) * diff_length_sq_inv(ip);
                 else
                     vViscoPerDiffLenSq[ip] = (mu_1/rho_1) * diff_length_sq_inv(ip);
                 
                 interface_change[ip] = false
             }
                 
             
         }
             
     }



      number vNormStdVelPerConvLen[numIp];
     if(!bStokes)
         for(size_t ip = 0; ip < numIp; ++ip)
             vNormStdVelPerConvLen[ip] = VecTwoNorm(vStdVel[ip]) / upwind_conv_length(ip);

 //    Find out if upwinded velocities depend on other ip velocities. In that case
 //    we have to solve a matrix system. Else the system is diagonal and we can
 //    compute the inverse directly

 //    diagonal case (i.e. upwind vel depend only on corner vel or no upwind)

         
     if(bStokes || !non_zero_shape_ip())
     {
         if(!multiphase)
         {
             //    We can solve the systems ip by ip
             for(size_t ip = 0; ip < numIp; ++ip)
             {
                 //    get SubControlVolumeFace
                 const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
                 
                 //    First, we compute the contributions to the diagonal
                 //    Note: - There is no contribution of the upwind vel to the diagonal
                 //            in this case, only for non-diag problems
                 //          - The diag does not depend on the dimension
                 
                 //    Diffusion part
                 number diag = vViscoPerDiffLenSq[ip];
                 
                 //    Time part
                 if(pvCornerValueOldTime != NULL)
                     diag += 1./dt;
                 
                 
                 //    Convective Term (no convective terms in the Stokes eq.)
                 if (! bStokes)
                     diag += vNormStdVelPerConvLen[ip];
                 
                 //     Loop components of velocity
                 for(size_t d = 0; d < (size_t)dim; d++)
                 {
                     //    Now, we can assemble the rhs. This rhs is assembled by all
                     //    terms, that are non-dependent on the ip vel.
                     //    Note, that we can compute the stab_shapes on the fly when setting
                     //    up the system.
                     
                     //    Source
                     number rhs = 0.0;
                     if(pSource != NULL)
                         rhs = (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];
                     
                     //    Time
                     if(pvCornerValueOldTime != NULL)
                     {
                         //    interpolate old time step
                         number oldIPVel = 0.0;
                         for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                             oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                         
                         //    add to rhs
                         rhs += oldIPVel / dt;
                     }
                     
                     //    loop shape functions
                     for(size_t k = 0; k < scvf.num_sh(); ++k)
                     {
                         //    Diffusion part
                         number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                         
                         //    Convective term
                         if (! bStokes) // no convective terms in the Stokes eq.
                             sumVel += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
                         
                         //    Add to rhs
                         rhs += sumVel * vCornerValue(d, k);
                         
                         //    set stab shape
                         stab_shape_vel(ip, d, d, k) = sumVel / diag;
                         
                         
                         //    Pressure part
                         const number sumP = -1.0 * (scvf.global_grad(k))[d] / density[ip];
                         
                         //    Add to rhs
                         rhs += sumP * vCornerValue(_P_, k);
                         
                         //    set stab shape
                         stab_shape_p(ip, d, k) = sumP / diag;
                         
                         
                         if ( multiphase)
                         {
                             size_t from = scvf.from();
                             size_t to = scvf.to();
                             number alpha = 1;
                             number sumPJump;
                             if (jump_shape[from]*jump_shape[to]<0 || (jump_shape[from]<0 && jump_shape[from]<0))
                             {
                                 
                                 if (jump_shape[k]>0.0)
                                     sumPJump =  alpha*jump_shape[k] * (scvf.global_grad(k))[d] / density[ip];
                                 else
                                     sumPJump = 0.0;
                                 
                             }
                             else
                             {
                                 if (jump_shape[k]<0.0)
                                     sumPJump =  alpha*jump_shape[k] * (scvf.global_grad(k))[d] / density[ip];
                                 else
                                     sumPJump = 0.0;
                             }
                             
                             
                             //    Add to rhs
                             rhs += sumPJump * (pressure_jump[k]);
                             stab_shape_p_jump(ip, d, k) = sumPJump/diag;
                             
                         }
                         
                     }
                     
                     //    Finally, the can invert this row
                     stab_vel(ip)[d] = rhs / diag;
                 }
             }
         }
         else
         {
             number mu_2 = 0.0, mu_1 = 0.0;
             number rho_2 = 0.0, rho_1 = 0.0;

             number visc1 = -10000000, visc2 = -10000000;
             number rho1 = -10000000, rho2 = -10000000;
             
             for(size_t sh = 0; sh < numSh; ++sh)
             {
                 if (jump_shape[sh]>0)
                 {
                     visc2 = fmax( visc2, densitySCV[sh] * kinViscoSCV[sh]);
                     rho2 = fmax( rho2, densitySCV[sh] );
                 }
                 else
                 {
                     visc1 = fmax( visc1, -densitySCV[sh] * kinViscoSCV[sh]);
                     rho1 = fmax( rho1, -densitySCV[sh] );
                 }
                 
             }
             
             mu_2 = visc2;
             mu_1 = -visc1;
             rho_2 = rho2;
             rho_1 = -rho1;
             
             if ((mu_2 < mu_1)||(mu_2<0.0)||(mu_1<0.0))
                 UG_THROW("Viscosity in phase 1 is lower that phase 2");
             if ((rho_2 < rho_1)||(rho_2<0.0)||(rho_1<0.0))
                 UG_THROW("Density in phase 1 is lower that phase 2");
             
             
             for(size_t ip = 0; ip < numIp; ++ip)
             {
                 const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
                 
                 if (jump_shape[from]*jump_shape[to]<0
                     vViscoPerDiffLenSq[ip] = (mu_2/rho_2- mu_1/rho_1) * diff_length_sq_inv(ip);
             }
             
             //    We can solve the systems ip by ip
             for(size_t ip = 0; ip < numIp; ++ip)
             {
                 //    get SubControlVolumeFace
                 const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
                 const size_t from=scvf.from();
                 const size_t to=scvf.to();
                 //    First, we compute the contributions to the diagonal
                 //    Note: - There is no contribution of the upwind vel to the diagonal
                 //            in this case, only for non-diag problems
                 //          - The diag does not depend on the dimension
                 
                 //    Diffusion part
                 number diag = vViscoPerDiffLenSq[ip];
                 
                 
                 //     Loop components of velocity
                 for(size_t d = 0; d < (size_t)dim; d++)
                 {
                     //    Now, we can assemble the rhs. This rhs is assembled by all
                     //    terms, that are non-dependent on the ip vel.
                     //    Note, that we can compute the stab_shapes on the fly when setting
                     //    up the system.
                     if (jump_shape[from]*jump_shape[to]<0)
                     {
                         //    Source
                         number rhs = 0.0;
                         if(pSource != NULL)
                             rhs =   - ( 1.0/rho_2-1.0/rho1 ) * density_ref * (*pSource)[ip][d];
                         
                         
                         //    loop shape functions
                         for(size_t k = 0; k < scvf.num_sh(); ++k)
                         {
                             //    Diffusion part
                             number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                             
                             //    Add to rhs
                             rhs += sumVel * vCornerValue(d, k);
                             
                             //    set stab shape
                             stab_shape_vel(ip, d, d, k) = sumVel / diag;
                             
                             
                             //    Pressure part
                             const number sumP = +1.0 * (scvf.global_grad(k))[d] * (1.0/rho_2-1.0/rho_1);
                             
                             //    Add to rhs
                             rhs += sumP * vCornerValue(_P_, k);
                             
                             //    set stab shape
                             stab_shape_p(ip, d, k) = sumP / diag;
                             
                             
                             if ( multiphase)
                             {
                                 size_t from = scvf.from();
                                 size_t to = scvf.to();
                                 number alpha = 1;
                                 number sumPJump = 1.0/rho_2;
                                 if (jump_shape[k]>0)
                                     sumPJump += -1.0 * jump_shape[k]* (1.0/rho_2 - 1.0/rho_1);
                                 sumPJump *=  scvf.global_grad(k)[d];
                                 
                                 //    Add to rhs
                                 rhs += sumPJump * (pressure_jump[k]);
                                 stab_shape_p_jump(ip, d, k) = sumPJump/diag;
                                 
                             }
                             
                         }
                     }
                     else
                     {
                         number rho;

                         rho= (jump_shape[from])<0.0)? rho_1:rho_2;
                             
                         
                         number rhs = 0.0;
                         if(pSource != NULL)
                             rhs = (  ( rho-density_ref) / rho) * (*pSource)[ip][d];
                         
                         //    Time
                         if(pvCornerValueOldTime != NULL)
                         {
                             //    interpolate old time step
                             number oldIPVel = 0.0;
                             for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                                 oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                             
                             //    add to rhs
                             rhs += oldIPVel / dt;
                         }
                         
                         //    loop shape functions
                         for(size_t k = 0; k < scvf.num_sh(); ++k)
                         {
                             //    Diffusion part
                             number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                             
                             //    Convective term
                             if (! bStokes) // no convective terms in the Stokes eq.
                                 sumVel += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
                             
                             //    Add to rhs
                             rhs += sumVel * vCornerValue(d, k);
                             
                             //    set stab shape
                             stab_shape_vel(ip, d, d, k) = sumVel / diag;
                             
                             
                             //    Pressure part
                             const number sumP = -1.0 * (scvf.global_grad(k))[d] / rho;
                             
                             //    Add to rhs
                             rhs += sumP * vCornerValue(_P_, k);
                             
                             //    set stab shape
                             stab_shape_p(ip, d, k) = sumP / diag;
                             
                             
                             if ( multiphase)
                             {
     
                                 number sumPJump;
                                 if (jump_shape[from]*jump_shape[to]<0 || (jump_shape[from]<0 && jump_shape[from]<0))
                                 {
                                     
                                     if (jump_shape[k]>0.0)
                                         sumPJump =  alpha*jump_shape[k] * (scvf.global_grad(k))[d] / density[ip];
                                     else
                                         sumPJump = 0.0;
                                     
                                 }
                                 else
                                 {
                                     if (jump_shape[k]<0.0)
                                         sumPJump =  alpha*jump_shape[k] * (scvf.global_grad(k))[d] / density[ip];
                                     else
                                         sumPJump = 0.0;
                                 }
                                 
                                 
                                 //    Add to rhs
                                 rhs += sumPJump * (pressure_jump[k]);
                                 stab_shape_p_jump(ip, d, k) = sumPJump/diag;
                                 
                             }
                             
                         }
                         
                     }
                     
                     //    Finally, the can invert this row
                     stab_vel(ip)[d] = rhs / diag;
                 }
             }
             
             
             
         
         }
     }
     
     /// need to solve system
     else
     {
         UG_THROW("Not implemented for ip velocities depending on other ip.");
     }
      // end switch for non-diag
 }
 
 */
                    
                    
                    
/* {
 //    abbreviation for pressure
     static const size_t _P_ = dim;

 //    Some constants
     static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
     static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;

 //    compute upwind (no convective terms for the Stokes eq. => no upwind)
     if (! bStokes) this->compute_upwind(geo, vStdVel);

 //    compute diffusion length
     this->compute_diff_length(*geo);

 //    cache values
     number vViscoPerDiffLenSq[numIp];
     for(size_t ip = 0; ip < numIp; ++ip)
     {
         vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);
     }



      number vNormStdVelPerConvLen[numIp];
     if(!bStokes)
         for(size_t ip = 0; ip < numIp; ++ip)
             vNormStdVelPerConvLen[ip] = VecTwoNorm(vStdVel[ip]) / upwind_conv_length(ip);

 //    Find out if upwinded velocities depend on other ip velocities. In that case
 //    we have to solve a matrix system. Else the system is diagonal and we can
 //    compute the inverse directly

 //    diagonal case (i.e. upwind vel depend only on corner vel or no upwind)

         
     if(bStokes || !non_zero_shape_ip())
     {
     //    We can solve the systems ip by ip
         for(size_t ip = 0; ip < numIp; ++ip)
         {
         //    get SubControlVolumeFace
             const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

         //    First, we compute the contributions to the diagonal
         //    Note: - There is no contribution of the upwind vel to the diagonal
         //            in this case, only for non-diag problems
         //          - The diag does not depend on the dimension

         //    Diffusion part
             number diag = vViscoPerDiffLenSq[ip];

         //    Time part
             if(pvCornerValueOldTime != NULL)
                 diag += 1./dt;
                 

         //    Convective Term (no convective terms in the Stokes eq.)
             if (! bStokes)
                 diag += vNormStdVelPerConvLen[ip];

         //     Loop components of velocity
             for(size_t d = 0; d < (size_t)dim; d++)
             {
             //    Now, we can assemble the rhs. This rhs is assembled by all
             //    terms, that are non-dependent on the ip vel.
             //    Note, that we can compute the stab_shapes on the fly when setting
             //    up the system.

             //    Source
                 number rhs = 0.0;
                 if(pSource != NULL)
                     rhs = (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];

             //    Time
                 if(pvCornerValueOldTime != NULL)
                 {
                 //    interpolate old time step
                     number oldIPVel = 0.0;
                     for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                         oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);

                 //    add to rhs
                     rhs += oldIPVel / dt;
                 }

             //    loop shape functions
                 for(size_t k = 0; k < scvf.num_sh(); ++k)
                 {
                 //    Diffusion part
                     number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);

                 //    Convective term
                     if (! bStokes) // no convective terms in the Stokes eq.
                         sumVel += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);

                 //    Add to rhs
                     rhs += sumVel * vCornerValue(d, k);

                 //    set stab shape
                     stab_shape_vel(ip, d, d, k) = sumVel / diag;
                     

                 //    Pressure part
                     const number sumP = -1.0 * (scvf.global_grad(k))[d] / density[ip];

                 //    Add to rhs
                     rhs += sumP * vCornerValue(_P_, k);

                 //    set stab shape
                     stab_shape_p(ip, d, k) = sumP / diag;
                     
                     
                     if ( multiphase)
                     {
                         size_t from = scvf.from();
                         size_t to = scvf.to();
                         number alpha = 1;
                         number sumPJump;
                         if (jump_shape[from]*jump_shape[to]<0 || (jump_shape[from]<0 && jump_shape[to]))
                         {
                             
                             if (jump_shape[k]>0.0)
                                 sumPJump =  alpha*jump_shape[k] * (scvf.global_grad(k))[d] / density[ip];
                             else
                                 sumPJump = 0.0;
                             
                         }
                         else
                         {
                             if (jump_shape[k]<0.0)
                                 sumPJump =  alpha*jump_shape[k] * (scvf.global_grad(k))[d] / density[ip];
                             else
                                 sumPJump = 0.0;
                         }

                             
                         //    Add to rhs
                         rhs += sumPJump * (pressure_jump[k]);
                         stab_shape_p_jump(ip, d, k) = sumPJump/diag;

                     }

                 }

             //    Finally, the can invert this row
                 stab_vel(ip)[d] = rhs / diag;
             }
         }
     }
     /// need to solve system
     else
     {
         UG_THROW("Not implemented for ip velocities depending on other ip.");
     }
      // end switch for non-diag
 }*/


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
// FIELDS 3
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesVISCOSITY_Stabilization<TDim>::
update(const FV1Geometry<TElem, dim>* geo,
       const LocalVector& vCornerValue,
       const MathVector<dim> vStdVel[],
       const bool bStokes,
       const DataImport<number, dim>& kinVisco,
       const DataImport<number, dim>& kinViscoSCV,
       const DataImport<number, dim>& density,
       const DataImport<number, dim>& densitySCV,
       const number pressure_jump[],
       const MathVector<dim> SlipVel[],
       const number jump_shape[],
       const MathVector<dim> normal[],
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt,
       const number density_ref,
       const bool multiphase,
       const bool phase_2[],
       const number theta,
       number** SCVFinterShape)
{
    if( non_zero_shape_ip())
    {
        UG_THROW("Not implemented for ip velocities depending on other ip.");
    }
    //    abbreviation for pressure
    static const size_t _P_ = dim;
    
    //    Some constants
    static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
    static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
    
    //    compute upwind and downwind (no convective terms for the Stokes eq. => no upwind)
    if (! bStokes)
    {
        this->compute_upwind(geo, vStdVel);
        this->compute_downwind(geo, vStdVel);
    }
    
    //    compute diffusion length
    this->compute_diff_length(*geo);
    
    number mu_2 = 0.0, mu_1 = 0.0;
    number rho_2 = 0.0, rho_1 = 0.0;
    number visc1 = -10000000, visc2 = -10000000;
    number rho1 = -10000000, rho2 = -10000000;
    
    number alpha1 = 1.0;
    number alpha2 = 1.0;
    number alpha3 = 1.0;
    

    
    if (multiphase)
    {

        
        for(size_t sh = 0; sh < numSh; ++sh)
        {
            if (jump_shape[sh]>0)
            {
                visc2 = fmax( visc2, densitySCV[sh] * kinViscoSCV[sh]);
                rho2 = fmax( rho2, densitySCV[sh] );
            }
            else
            {
                visc1 = fmax( visc1, -densitySCV[sh] * kinViscoSCV[sh]);
                rho1 = fmax( rho1, -densitySCV[sh] );
            }
            
        }
        
        mu_2 = visc2;
        mu_1 = -visc1;
        rho_2 = rho2;
        rho_1 = -rho1;
        
        if ((mu_2 < mu_1)||(mu_2<0.0)||(mu_1<0.0))
            UG_THROW("Viscosity in phase 1 is lower that phase 2");
        if ((rho_2 < rho_1)||(rho_2<0.0)||(rho_1<0.0))
            UG_THROW("Density in phase 1 is lower that phase 2");
    }
    

    //    cache values
    number vViscoPerDiffLenSq[numIp];
    bool interface_change[numIp];
    
    number vNormStdVelPerConvLen[numIp];
    number vNormStdVelPerDownLen[numIp];
    
    for(size_t ip = 0; ip < numIp; ++ip)
    {
        if (multiphase)
        {
            const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
            if(jump_shape[scvf.to()]*jump_shape[scvf.from()]<0.0)
            {
                interface_change[ip] = true;
            }
                
            else
            {
                interface_change[ip] = false;
            }
            if (phase_2[ip])
                vViscoPerDiffLenSq[ip] = (mu_2/rho_2) * diff_length_sq_inv(ip);
            else
                vViscoPerDiffLenSq[ip] = (mu_1/rho_1) * diff_length_sq_inv(ip);
            
        }
        else
            vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);
        
        if(!bStokes)
        {
            const number norm = VecTwoNorm(vStdVel[ip]);
            vNormStdVelPerConvLen[ip] = norm / upwind_conv_length(ip);
            vNormStdVelPerDownLen[ip] = norm / (downwind_conv_length(ip) + upwind_conv_length(ip));
            
        }


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
        
        if (!multiphase )
        {

            
            //    First, we compute the contributions to the diagonal
            //    Note: - There is no contribution of the upwind vel to the diagonal
            //            in this case, only for non-diag problems
            //          - The diag does not depend on the dimension
            
            //    the diagonal entry
            number diag = vViscoPerDiffLenSq[ip];
            
            //    Time part
            if(pvCornerValueOldTime != NULL)
                diag += 1./dt;
            
            //    Convective Term  (no convective terms in the Stokes eq.)
            if (! bStokes)
                diag += vNormStdVelPerConvLen[ip];
            
            //     Loop components of velocity
            for(int d = 0; d < dim; d++)
            {
                //    Now, we can assemble the rhs. This rhs is assembled by all
                //    terms, that are non-dependent on the ip vel.
                //    Note, that we can compute the stab_shapes on the fly when setting
                //    up the system.
                
                //    Source
                number rhs = 0.0;
                if(pSource != NULL)
                    rhs = (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];
                
                //    Time
                if(pvCornerValueOldTime != NULL)
                {
                    //    interpolate old time step
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                    
                    //    add to rhs
                    rhs += oldIPVel / dt;
                }
                
                //    loop shape functions
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                    //    Diffusion part
                    number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                    
                    //    Convective term (no convective terms in the Stokes eq.)
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
                    
                    //    Add to rhs
                    rhs += sumVel * vCornerValue(d, k);
                    
                    //    set stab shape
                    stab_shape_vel(ip, d, d, k) = sumVel / diag;
                    
                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        if(d2 == d) continue;
                        
                        const number sumVel2 = vStdVel[ip][d] * (scvf.global_grad(k))[d2];
                        
                        rhs += sumVel2 * vCornerValue(d2, k);
                        
                        stab_shape_vel(ip, d, d2, k) = sumVel2 / diag;
                    }
                    
                    //    Pressure part
                    number sumP;
                    
                    sumP = -1.0 * alpha1 * (scvf.global_grad(k))[d]  / density[ip];
                    
                    //    Add to rhs
                    rhs += sumP * vCornerValue(_P_, k);
                    
                    //    set stab shape
                    stab_shape_p(ip, d, k) = sumP / diag;
                    
                    number sumPJump =0.0;

                    stab_shape_p_jump(ip, d, k) = sumPJump/diag;
                }
                
                //    Finally, the can invert this row
                stab_vel(ip)[d] = rhs / diag;
                
            }
        }
        else
        {
            
            
            //    First, we compute the contributions to the diagonal
            //    Note: - There is no contribution of the upwind vel to the diagonal
            //            in this case, only for non-diag problems
            //          - The diag does not depend on the dimension
            
            //    the diagonal entry
            number diag = vViscoPerDiffLenSq[ip];
            
            //    Time part
            if(pvCornerValueOldTime != NULL)
                diag += 1./dt;
            
            //    Convective Term  (no convective terms in the Stokes eq.)
            if (! bStokes)
                diag += vNormStdVelPerConvLen[ip];
            

            
            //     Loop components of velocity
            for(int d = 0; d < dim; d++)
            {
                //    Now, we can assemble the rhs. This rhs is assembled by all
                //    terms, that are non-dependent on the ip vel.
                //    Note, that we can compute the stab_shapes on the fly when setting
                //    up the system.
                
                //    Source
                number rhs = 0.0;
                if(pSource != NULL)
                    rhs = (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];
                //rhs = (  ( density[ip]-density_ref) / density[ip]) * VecProd((*pSource)[ip], normal[ip]) * normal[ip][d];
                
                //    Time
                if(pvCornerValueOldTime != NULL)
                {
                    //    interpolate old time step
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        
                            oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                    //for(int d2 = 0; d2 < dim; d2++)
                    //oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d2, sh) * normal [ip][d2] * normal[ip][d];
                    
                    
                    //    add to rhs
                    rhs += oldIPVel / dt;
                }
                
                //    loop shape functions
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                    //    Diffusion part
                    number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                    
                    //    Convective term (no convective terms in the Stokes eq.)
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
                    
                    //    Add to rhs
                    rhs += sumVel * vCornerValue(d, k);
                    
                    //    set stab shape
                    stab_shape_vel(ip, d, d, k) = sumVel / diag;
                    
                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        if(d2 == d) continue;
                        
                        const number sumVel2 = vStdVel[ip][d] * (scvf.global_grad(k))[d2];
                        
                        rhs += sumVel2 * vCornerValue(d2, k);
                        
                        stab_shape_vel(ip, d, d2, k) = sumVel2 / diag;
                    }
                    
                    //    Pressure part
                    number sumP = 0.0;
                    if (interface_change[ip])
                    {
                        sumP = -1.0 * alpha3 * (scvf.global_grad(k)[d] ) / density[ip];
                        //sumP += -1.0 * alpha3 * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]  / density[ip];
                        //sumP += -1.0 * alpha3 * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d])/ density[ip];
                        
                    }
                    else
                    {
                        sumP =  -1.0 * alpha2 * (scvf.global_grad(k)[d] ) / density[ip];
                        //sumP += -1.0 * alpha2 * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]  / density[ip];
                        //sumP = -1.0 * alpha2 * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d])/ density[ip];
                    }
                    
                    //    Add to rhs
                    rhs += sumP * vCornerValue(_P_, k);
                    
                    //    set stab shape
                    stab_shape_p(ip, d, k) = sumP / diag;
                    
                    number sumPJump =0.0;
                    
                    
                    if (interface_change[ip] )
                    {
                        if((phase_2[ip] && jump_shape[k]<0.0) || (!phase_2[ip] && jump_shape[k]>0.0) )
                        {
                            //sumPJump =  alpha3 * jump_shape[k] * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d] / rho_1;
                            //sumPJump +=  alpha3 * jump_shape[k] * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d])/ density[ip];
                            sumPJump =  alpha3 * jump_shape[k] * (scvf.global_grad(k)[d] ) / density[ip];
                        }

                    }
                    else
                    {
                        if ((phase_2[ip] && jump_shape[k]<0.0) || (!phase_2[ip] && jump_shape[k]>0.0))
                        {
                            sumPJump =  alpha2 * jump_shape[k] * (scvf.global_grad(k)[d] ) / density[ip];
                            //sumPJump =  alpha2 * jump_shape[k] * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]  / density[ip];
                            //sumPJump =  alpha2 * jump_shape[k] * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) / rho_1;
                        }
                        
                    }
                        
                    //    Add to rhs
                    rhs += sumPJump *(pressure_jump[k]);
                    
                    
                    stab_shape_p_jump(ip, d, k) = sumPJump/diag;
                    
                    
                    
                    number sumSlipVel = 0.0;
                    if ((phase_2[ip] && jump_shape[k]<0.0) )
                        sumSlipVel = vViscoPerDiffLenSq[ip] * scvf.shape(k) * jump_shape[k];
                    
                    rhs += sumSlipVel * SlipVel[k][d];
                    stab_shape_slip_vel(ip, d, d, k) = sumSlipVel/diag;

                }
                
                //    Finally, the can invert this row
                stab_vel(ip)[d] = rhs / diag;
                
            }
            
            /*if (interface_change[ip])
            {
                VecScale(stab_vel(ip),normal[ip],VecProd(stab_vel(ip), normal[ip]));
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                    for(int d1 = 0; d1 < dim; d1++)
                    {
                        number ValueV = 0.0;
                        number ValueP = 0.0;
                        number ValueJumpP = 0.0;
                        number ValueSlipV = 0.0;
                        for(int d2 = 0; d2 < dim; d2++)
                        {
                            ValueV += stab_shape_vel(ip, d2, d1, k) * normal[ip][d2];
                            ValueP += stab_shape_p(ip, d2, k) * normal[ip][d2];
                            ValueJumpP += stab_shape_p_jump(ip, d2, k) * normal[ip][d2];
                            ValueSlipV += stab_shape_slip_vel(ip, d2, d1, k) * normal[ip][d2];
                            
                        }
                        for(int d2 = 0; d2 < dim; d2++)
                        {
                            stab_shape_vel(ip, d2, d1, k) = ValueV * normal[ip][d2];
                            stab_shape_slip_vel(ip, d2, d1, k) = ValueSlipV * normal[ip][d2];
                            
                        }
                        stab_shape_p(ip, d1, k) = ValueP * normal[ip][d1];
                        stab_shape_p_jump(ip, d1, k) = ValueJumpP * normal[ip][d1];
                    }
                
            
            }
            else
            {
                MathVector<dim> T = SlipVel[0];
                VecScale(T,T,sqrt(VecProd(T,T)));
                VecScale(stab_vel(ip),T,VecProd(stab_vel(ip), T));
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                    for(int d1 = 0; d1 < dim; d1++)
                    {
                        number ValueV = 0.0;
                        number ValueP = 0.0;
                        number ValueJumpP = 0.0;
                        number ValueSlipV = 0.0;
                        for(int d2 = 0; d2 < dim; d2++)
                        {
                            ValueV += stab_shape_vel(ip, d2, d1, k) * T[d2];
                            ValueP += stab_shape_p(ip, d2, k) * T[d2];
                            ValueJumpP += stab_shape_p_jump(ip, d2, k) * T[d2];
                            ValueSlipV += stab_shape_slip_vel(ip, d2, d1, k) * T[d2];
                            
                        }
                        for(int d2 = 0; d2 < dim; d2++)
                        {
                            stab_shape_vel(ip, d2, d1, k) = ValueV * T[d2];
                            
                        }
                        stab_shape_p(ip, d1, k) = ValueP * T[d1];
                        stab_shape_p_jump(ip, d1, k) = ValueJumpP * T[d1];
                    }
                
                
            }*/
            
            
                
        }
    }
    /*
     number sumVel2;
     if ( d == d2)
         sumVel2 = scvf.shape(k);
     else
         sumVel2 = 0.0;
     //sumVel2 = scvf.shape(k) * normal[k][d2] *normal[k][d];
     */

}
/*{
    if( non_zero_shape_ip())
    {
        UG_THROW("Not implemented for ip velocities depending on other ip.");
    }
    //    abbreviation for pressure
    static const size_t _P_ = dim;
    
    //    Some constants
    static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
    static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
    
    //    compute upwind and downwind (no convective terms for the Stokes eq. => no upwind)
    if (! bStokes)
    {
        this->compute_upwind(geo, vStdVel);
        this->compute_downwind(geo, vStdVel);
    }
    
    //    compute diffusion length
    this->compute_diff_length(*geo);
    
    number mu_2 = 0.0, mu_1 = 0.0;
    number rho_2 = 0.0, rho_1 = 0.0;
    number visc1 = -10000000, visc2 = -10000000;
    number rho1 = -10000000, rho2 = -10000000;
    
    number alpha1 = 1.0;
    number alpha2 = 1.0;
    number alpha3 = 0.0;
    

    
    if (multiphase)
    {

        
        for(size_t sh = 0; sh < numSh; ++sh)
        {
            if (jump_shape[sh]>0)
            {
                visc2 = fmax( visc2, densitySCV[sh] * kinViscoSCV[sh]);
                rho2 = fmax( rho2, densitySCV[sh] );
            }
            else
            {
                visc1 = fmax( visc1, -densitySCV[sh] * kinViscoSCV[sh]);
                rho1 = fmax( rho1, -densitySCV[sh] );
            }
            
        }
        
        mu_2 = visc2;
        mu_1 = -visc1;
        rho_2 = rho2;
        rho_1 = -rho1;
        
        if ((mu_2 < mu_1)||(mu_2<0.0)||(mu_1<0.0))
            UG_THROW("Viscosity in phase 1 is lower that phase 2");
        if ((rho_2 < rho_1)||(rho_2<0.0)||(rho_1<0.0))
            UG_THROW("Density in phase 1 is lower that phase 2");
    }
    

    //    cache values
    number vViscoPerDiffLenSq[numIp];
    bool interface_change[numIp];
    
    number vNormStdVelPerConvLen[numIp];
    number vNormStdVelPerDownLen[numIp];
    
    for(size_t ip = 0; ip < numIp; ++ip)
    {
        if (multiphase)
        {
            const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
            if(jump_shape[scvf.to()]*jump_shape[scvf.from()]<0.0)
            {
                interface_change[ip] = true;
            }
                
            else
            {
                interface_change[ip] = false;
            }
            if (phase_2[ip])
                vViscoPerDiffLenSq[ip] = (mu_2/rho_2) * diff_length_sq_inv(ip);
            else
                vViscoPerDiffLenSq[ip] = (mu_1/rho_1) * diff_length_sq_inv(ip);
            
        }
        else
            vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);
        
        if(!bStokes)
        {
            const number norm = VecTwoNorm(vStdVel[ip]);
            vNormStdVelPerConvLen[ip] = norm / upwind_conv_length(ip);
            vNormStdVelPerDownLen[ip] = norm / (downwind_conv_length(ip) + upwind_conv_length(ip));
            
        }


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
        
        if (!multiphase )
        {

            
            //    First, we compute the contributions to the diagonal
            //    Note: - There is no contribution of the upwind vel to the diagonal
            //            in this case, only for non-diag problems
            //          - The diag does not depend on the dimension
            
            //    the diagonal entry
            number diag = vViscoPerDiffLenSq[ip];
            
            //    Time part
            if(pvCornerValueOldTime != NULL)
                diag += 1./dt;
            
            //    Convective Term  (no convective terms in the Stokes eq.)
            if (! bStokes)
                diag += vNormStdVelPerConvLen[ip];
            
            //     Loop components of velocity
            for(int d = 0; d < dim; d++)
            {
                //    Now, we can assemble the rhs. This rhs is assembled by all
                //    terms, that are non-dependent on the ip vel.
                //    Note, that we can compute the stab_shapes on the fly when setting
                //    up the system.
                
                //    Source
                number rhs = 0.0;
                if(pSource != NULL)
                    rhs = (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];
                
                //    Time
                if(pvCornerValueOldTime != NULL)
                {
                    //    interpolate old time step
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                    
                    //    add to rhs
                    rhs += oldIPVel / dt;
                }
                
                //    loop shape functions
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                    //    Diffusion part
                    number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                    
                    //    Convective term (no convective terms in the Stokes eq.)
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
                    
                    //    Add to rhs
                    rhs += sumVel * vCornerValue(d, k);
                    
                    //    set stab shape
                    stab_shape_vel(ip, d, d, k) = sumVel / diag;
                    
                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        if(d2 == d) continue;
                        
                        const number sumVel2 = vStdVel[ip][d] * (scvf.global_grad(k))[d2];
                        
                        rhs += sumVel2 * vCornerValue(d2, k);
                        
                        stab_shape_vel(ip, d, d2, k) = sumVel2 / diag;
                    }
                    
                    //    Pressure part
                    number sumP;
                    
                    sumP = -1.0 * alpha1 * (scvf.global_grad(k))[d]  / density[ip];
                    
                    //    Add to rhs
                    rhs += sumP * vCornerValue(_P_, k);
                    
                    //    set stab shape
                    stab_shape_p(ip, d, k) = sumP / diag;
                    
                    number sumPJump =0.0;

                    stab_shape_p_jump(ip, d, k) = sumPJump/diag;
                }
                
                //    Finally, the can invert this row
                stab_vel(ip)[d] = rhs / diag;
                
            }
        }
        else
        {
            
            
            //    First, we compute the contributions to the diagonal
            //    Note: - There is no contribution of the upwind vel to the diagonal
            //            in this case, only for non-diag problems
            //          - The diag does not depend on the dimension
            
            //    the diagonal entry
            number diag = vViscoPerDiffLenSq[ip];
            
            //    Time part
            if(pvCornerValueOldTime != NULL)
                diag += 1./dt;
            
            //    Convective Term  (no convective terms in the Stokes eq.)
            if (! bStokes)
                diag += vNormStdVelPerConvLen[ip];
            
            //     Loop components of velocity
            for(int d = 0; d < dim; d++)
            {
                //    Now, we can assemble the rhs. This rhs is assembled by all
                //    terms, that are non-dependent on the ip vel.
                //    Note, that we can compute the stab_shapes on the fly when setting
                //    up the system.
                
                //    Source
                number rhs = 0.0;
                if(pSource != NULL)
                    rhs = (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];
                //rhs = (  ( density[ip]-density_ref) / density[ip]) * VecProd((*pSource)[ip], normal[ip]) * normal[ip][d];
                
                //    Time
                if(pvCornerValueOldTime != NULL)
                {
                    //    interpolate old time step
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        
                            oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                    //for(int d2 = 0; d2 < dim; d2++)
                    //oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d2, sh) * normal [ip][d2] * normal[ip][d];
                    
                    
                    //    add to rhs
                    rhs += oldIPVel / dt;
                }
                
                //    loop shape functions
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                    //    Diffusion part
                    number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                    
                    //    Convective term (no convective terms in the Stokes eq.)
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
                    
                    //    Add to rhs
                    rhs += sumVel * vCornerValue(d, k);
                    
                    //    set stab shape
                    stab_shape_vel(ip, d, d, k) = sumVel / diag;
                    
                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        if(d2 == d) continue;
                        
                        const number sumVel2 = vStdVel[ip][d] * (scvf.global_grad(k))[d2];
                        
                        rhs += sumVel2 * vCornerValue(d2, k);
                        
                        stab_shape_vel(ip, d, d2, k) = sumVel2 / diag;
                    }
                    
                    //    Pressure part
                    number sumP = 0.0;
                    if (interface_change[ip])
                    {
                        //sumP = -1.0 * alpha3 * (scvf.global_grad(k)[d] ) / density[ip];
                        //sumP = -1.0 * alpha3 * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]  / density[ip];
                        sumP = -1.0 * alpha3 * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d])/ density[ip];
                        
                    }
                    else
                    {
                        //sumP =  -1.0 * alpha2 * (scvf.global_grad(k)[d] ) / density[ip];
                        //sumP += -1.0 * alpha3 * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]  / density[ip];
                        sumP = -1.0 * alpha2 * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d])/ density[ip];
                    }
                    
                    //    Add to rhs
                    rhs += sumP * vCornerValue(_P_, k);
                    
                    //    set stab shape
                    stab_shape_p(ip, d, k) = sumP / diag;
                    
                    number sumPJump =0.0;
                    
                    
                    if (interface_change[ip] )
                    {
                        if((phase_2[ip] && jump_shape[k]<0.0) || (!phase_2[ip] && jump_shape[k]>0.0) )
                        {
                            //sumPJump =  alpha3 * jump_shape[k] * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d] / rho_1;
                            sumPJump =  alpha3 * jump_shape[k] * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d])/ density[ip];
                            //sumPJump =  alpha3 * jump_shape[k] * (scvf.global_grad(k)[d] ) / density[ip];
                        }

                    }
                    else
                    {
                        if ((phase_2[ip] && jump_shape[k]<0.0) || (!phase_2[ip] && jump_shape[k]>0.0))
                        {
                            //sumPJump =  alpha2 * jump_shape[k] * (scvf.global_grad(k)[d] ) / density[ip];
                            //sumPJump =  alpha2 * jump_shape[k] * VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]  / density[ip];
                            sumPJump =  alpha2 * jump_shape[k] * (scvf.global_grad(k)[d]  - VecProd(scvf.global_grad(k), normal[k] ) * normal[k][d]) / rho_1;
                        }
                        
                    }
                        
                    //    Add to rhs
                    rhs += sumPJump *(pressure_jump[k]);
                    
                    
                    stab_shape_p_jump(ip, d, k) = sumPJump/diag;
                    
                    
                    
                    number sumSlipVel = 0.0;
                    if ((phase_2[ip] && jump_shape[k]<0.0) || (!phase_2[ip] && jump_shape[k]>0.0))
                        sumSlipVel = -vViscoPerDiffLenSq[ip] * scvf.shape(k) * jump_shape[k];
                    
                    rhs += sumSlipVel * SlipVel[k][d];
                    stab_shape_slip_vel(ip, d, d, k) = sumSlipVel/diag;

                }
                
                //    Finally, the can invert this row
                stab_vel(ip)[d] = rhs / diag;
                
            }
            
            if (interface_change[ip])
            {
                VecScale(stab_vel(ip),normal[ip],VecProd(stab_vel(ip), normal[ip]));
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                    for(int d1 = 0; d1 < dim; d1++)
                    {
                        number ValueV = 0.0;
                        number ValueP = 0.0;
                        number ValueJumpP = 0.0;
                        number ValueSlipV = 0.0;
                        for(int d2 = 0; d2 < dim; d2++)
                        {
                            ValueV += stab_shape_vel(ip, d2, d1, k) * normal[ip][d2];
                            ValueP += stab_shape_p(ip, d2, k) * normal[ip][d2];
                            ValueJumpP += stab_shape_p_jump(ip, d2, k) * normal[ip][d2];
                            ValueSlipV += stab_shape_slip_vel(ip, d2, d1, k) * normal[ip][d2];
                            
                        }
                        for(int d2 = 0; d2 < dim; d2++)
                        {
                            stab_shape_vel(ip, d2, d1, k) = ValueV * normal[ip][d2];
                            stab_shape_slip_vel(ip, d2, d1, k) = ValueSlipV * normal[ip][d2];
                            
                        }
                        stab_shape_p(ip, d1, k) = ValueP * normal[ip][d1];
                        stab_shape_p_jump(ip, d1, k) = ValueJumpP * normal[ip][d1];
                    }
                
            
            }
            else
            {
                MathVector<dim> T = SlipVel[0];
                VecScale(T,T,sqrt(VecProd(T,T)));
                VecScale(stab_vel(ip),T,VecProd(stab_vel(ip), T));
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                    for(int d1 = 0; d1 < dim; d1++)
                    {
                        number ValueV = 0.0;
                        number ValueP = 0.0;
                        number ValueJumpP = 0.0;
                        number ValueSlipV = 0.0;
                        for(int d2 = 0; d2 < dim; d2++)
                        {
                            ValueV += stab_shape_vel(ip, d2, d1, k) * T[d2];
                            ValueP += stab_shape_p(ip, d2, k) * T[d2];
                            ValueJumpP += stab_shape_p_jump(ip, d2, k) * T[d2];
                            ValueSlipV += stab_shape_slip_vel(ip, d2, d1, k) * T[d2];
                            
                        }
                        for(int d2 = 0; d2 < dim; d2++)
                        {
                            stab_shape_vel(ip, d2, d1, k) = ValueV * T[d2];
                            
                        }
                        stab_shape_p(ip, d1, k) = ValueP * T[d1];
                        stab_shape_p_jump(ip, d1, k) = ValueJumpP * T[d1];
                    }
                
                
            }
            
            
                
        }
    }
    
     number sumVel2;
     if ( d == d2)
         sumVel2 = scvf.shape(k);
     else
         sumVel2 = 0.0;
     //sumVel2 = scvf.shape(k) * normal[k][d2] *normal[k][d];
     

}*/
/*
{
    if(multiphase)
        UG_THROW("Pressure Jump Not implemented for Viscosity stabilization.");
//    abbreviation for pressure
    static const size_t _P_ = dim;

//    Some constants
    static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
    static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;

//    compute upwind (no convective terms for the Stokes eq. => no upwind)
    if (! bStokes) this->compute_upwind(geo, vStdVel);

//    compute diffusion length
    this->compute_diff_length(*geo);

//    cache values
    number vViscoPerDiffLenSq[numIp];
    
    for(size_t ip = 0; ip < numIp; ++ip)
        vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);

     number vNormStdVelPerConvLen[numIp];
    if(!bStokes)
        for(size_t ip = 0; ip < numIp; ++ip)
            vNormStdVelPerConvLen[ip] = VecTwoNorm(vStdVel[ip]) / upwind_conv_length(ip);
    
    MathVector<dim> ViscosityGrad[numIp];
    for(size_t ip = 0; ip < numIp; ++ip)
    {
        const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
        VecSet(ViscosityGrad[ip],0.0);
        for(size_t sh = 0; sh < numSh; ++sh)
            for(size_t d = 0; d < (size_t)dim; d++)
                ViscosityGrad[ip][d] += density[ip] * (scvf.global_grad(sh))[d] * kinViscoSCV[sh];
    }
    
    for(size_t ip = 0; ip < numIp; ++ip)
        for(int d1 = 0; d1 < dim; ++d1)
            for(int d2 = 0; d2 < dim; ++d2)
                for(size_t k = 0; k < numSh; ++k)
                    stab_shape_vel(ip, d1, d2, k) = 0.0;

//    Find out if upwinded velocities depend on other ip velocities. In that case
//    we have to solve a matrix system. Else the system is diagonal and we can
//    compute the inverse directly

//    diagonal case (i.e. upwind vel depend only on corner vel or no upwind)
    if(bStokes || !non_zero_shape_ip())
    {
        

    //    We can solve the systems ip by ip
        for(size_t ip = 0; ip < numIp; ++ip)
        {
        //    get SubControlVolumeFace
            const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

        //    First, we compute the contributions to the diagonal
        //    Note: - There is no contribution of the upwind vel to the diagonal
        //            in this case, only for non-diag problems
        //          - The diag does not depend on the dimension

            


            
        //    Diffusion part
            number diag = vViscoPerDiffLenSq[ip];

        //    Time part
            if(pvCornerValueOldTime != NULL)
                //diag += 1.0 / dt;

        //    Convective Term (no convective terms in the Stokes eq.)
            if (! bStokes)
                diag +=  vNormStdVelPerConvLen[ip];

        //     Loop components of velocity
            for(size_t d = 0; d < (size_t)dim; d++)
            {
            //    Now, we can assemble the rhs. This rhs is assembled by all
            //    terms, that are non-dependent on the ip vel.
            //    Note, that we can compute the stab_shapes on the fly when setting
            //    up the system.

            //    Source
                number rhs = 0.0;
                if(pSource != NULL)
                    rhs =   ( density[ip]-density_ref)  * (*pSource)[ip][d] / density[ip];

            //    Time
                if(pvCornerValueOldTime != NULL)
                {
                //    interpolate old time step
                    number oldIPVel = 0.0;
                    
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                    {
                        oldIPVel += scvf.shape(sh) *  (*pvCornerValueOldTime)(d, sh);
                        //oldIPressure += -1.0 * (scvf.global_grad(sh))[d] *  (*pvCornerValueOldTime)(dim, sh) ;
                    }
                    
                        
                //    add to rhs
                    //rhs += oldIPVel / dt;
                    //rhs += oldIPressure;
                //   add density time derivative part from the continuity equation
                    //rhs -= vStdVel[ip][d] * (density[ip]-(*pDensity_OLD_SCVF)[ip]) /  dt ;
                }
                
                
 
            //    loop shape functions
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                //    Diffusion part
                    
                    
                    number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                    
                    //for(size_t d2 = 0; d2 < dim; ++d2) todo esto está comentado
                    {
                        sumVel += 1.0 * ViscosityGrad[ip][d2] * (scvf.global_grad(k))[d2] / density[ip];

                        const number sumVel2 = ViscosityGrad[ip][d2] * (scvf.global_grad(k))[d] / density[ip];
                        
                        if (d2!=d)
                        {
                            rhs += sumVel2 * vCornerValue(d2, k);
                            stab_shape_vel(ip, d, d2, k) += sumVel2 / diag;
                        }
                    }//
                                        
                //    Convective term
                    if (! bStokes) // no convective terms in the Stokes eq.
                    {
                        sumVel += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
                        //rhs -= vStdVel[ip][d]  * vNormStdVelPerConvLen[ip] * (density[ip]-upwind_shape_sh(ip, k) * densitySCV[k]);
                    }
                    

                //    Add to rhs
                    rhs += sumVel * vCornerValue(d, k);

                //    set stab shape
                    stab_shape_vel(ip, d, d, k) += sumVel / diag;

                //    Pressure part
                    const number sumP = - 1.0 * (scvf.global_grad(k))[d] / density[ip] ;

                //    Add to rhs
                    rhs += sumP * vCornerValue(_P_, k);

                //    set stab shape
                    stab_shape_p(ip, d, k) = sumP / diag;
                }

            //    Finally, the can invert this row
                stab_vel(ip)[d] = rhs / diag  ;
            }
        }
    }
    /// need to solve system
    else
    {
        UG_THROW("Not implemented for ip velocities depending on other ip.");
  

    } // end switch for non-diag
}*/
/*{
//    abbreviation for pressure
    static const size_t _P_ = dim;

//    Some constants
    static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
    static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;

//    compute upwind (no convective terms for the Stokes eq. => no upwind)
    if (! bStokes) this->compute_upwind(geo, vStdVel);

//    compute diffusion length
    this->compute_diff_length(*geo);

//    cache values
    number vViscoPerDiffLenSq[numIp];
    number upwind_density[numIp];
    for(size_t ip = 0; ip < numIp; ++ip)
    {
        vViscoPerDiffLenSq[ip] =  kinVisco[ip] * diff_length_sq_inv(ip);
        upwind_density[ip]=0.0;
        for(size_t k = 0; k < numSh; ++k)
            upwind_density[ip] += upwind_shape_sh(ip, k) * densitySCV[k];
    }
    
    

     number vNormStdVelPerConvLen[numIp];
    if(!bStokes)
        for(size_t ip = 0; ip < numIp; ++ip)
            vNormStdVelPerConvLen[ip] = VecTwoNorm(vStdVel[ip]) / upwind_conv_length(ip);

//    Find out if upwinded velocities depend on other ip velocities. In that case
//    we have to solve a matrix system. Else the system is diagonal and we can
//    compute the inverse directly

//    diagonal case (i.e. upwind vel depend only on corner vel or no upwind)
    if(bStokes || !non_zero_shape_ip())
    {
    //    We can solve the systems ip by ip
        for(size_t ip = 0; ip < numIp; ++ip)
        {
        //    get SubControlVolumeFace
            const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

        //    First, we compute the contributions to the diagonal
        //    Note: - There is no contribution of the upwind vel to the diagonal
        //            in this case, only for non-diag problems
        //          - The diag does not depend on the dimension

        //    Diffusion part
            number diag = vViscoPerDiffLenSq[ip];

        //    Time part
            if(pvCornerValueOldTime != NULL)
                diag += 1.0 / dt;

        //    Convective Term (no convective terms in the Stokes eq.)
            if (! bStokes)
                diag +=  vNormStdVelPerConvLen[ip];

        //     Loop components of velocity
            for(size_t d = 0; d < (size_t)dim; d++)
            {
            //    Now, we can assemble the rhs. This rhs is assembled by all
            //    terms, that are non-dependent on the ip vel.
            //    Note, that we can compute the stab_shapes on the fly when setting
            //    up the system.

            //    Source
                number rhs = 0.0;
                if(pSource != NULL)
                    rhs =   ( 1.0-density_ref/density[ip])  * (*pSource)[ip][d];

            //    Time
                if(pvCornerValueOldTime != NULL)
                {
                //    interpolate old time step
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPVel += scvf.shape(sh) *  (*pvCornerValueOldTime)(d, sh);

                //    add to rhs
                    rhs += oldIPVel / dt;
                //   add density time derivative part from the continuity equation
                    rhs -= ( vStdVel[ip][d] / density[ip]) * ( ( density[ip]-(*pDensity_OLD_SCVF)[ip] ) /  dt  + vNormStdVelPerConvLen[ip] * (density[ip]- upwind_density[ip]));
                }
                
                

            //    loop shape functions
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                //    Diffusion part
                    number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);

                //    Convective term
                    if (! bStokes) // no convective terms in the Stokes eq.
                    {
                        sumVel += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
                        //rhs -= vStdVel[ip][d]  * vNormStdVelPerConvLen[ip] * (density[ip]-upwind_shape_sh(ip, k) * densitySCV[k]);
                    }
                    

                //    Add to rhs
                    rhs += sumVel * vCornerValue(d, k);

                //    set stab shape
                    stab_shape_vel(ip, d, d, k) = sumVel / diag;

                //    Pressure part
                    const number sumP = -1.0 * (scvf.global_grad(k))[d] / density[ip];

                //    Add to rhs
                    rhs += sumP * vCornerValue(_P_, k);

                //    set stab shape
                    stab_shape_p(ip, d, k) = sumP / diag;
                }

            //    Finally, the can invert this row
                stab_vel(ip)[d] = rhs / diag;
            }
        }
    }
    /// need to solve system
    else
    {
        UG_THROW("Not implemented for ip velocities depending on other ip.");

    } // end switch for non-diag
}*/

template <>
void NavierStokesVISCOSITY_Stabilization<1>::register_func()
{
    register_func<RegularEdge>();
}

template <>
void NavierStokesVISCOSITY_Stabilization<2>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
}

template <>
void NavierStokesVISCOSITY_Stabilization<3>::register_func()
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
       const DataImport<number, dim>& densitySCV,
       const number pressure_jump[],
       const MathVector<dim> SlipVel[],
       const number jump_shape[],
       const MathVector<dim> normal[],
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt, 
       const number density_ref,
       const bool multiphase,
       const bool phase_2[],
       const number theta,
       number** SCVFinterShape)
{
    if(multiphase)
        UG_THROW("Pressure Jump Not implemented for KARIMIAN stabilization.");
    //    abbreviation for pressure
    static const size_t _P_ = dim;

    //    Some constants
    static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
    static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;

//    compute upwind and downwind (no convective terms for the Stokes eq. => no upwind)
    if (! bStokes)
    {
        this->compute_upwind(geo, vStdVel);
    }

//    compute diffusion length
    this->compute_diff_length(*geo);

//    cache values
    number vViscoPerDiffLenSq[numIp];
    number DenGrad[numIp];
    
    for(size_t ip = 0; ip < numIp; ++ip)
    {
        vViscoPerDiffLenSq[ip] =  kinVisco[ip] * diff_length_sq_inv(ip);
        
        const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
        DenGrad[ip]=0;
        for(int d = 0; d < dim; ++d)
            for(size_t sh = 0; sh < numSh; ++sh)
                DenGrad[ip] += densitySCV[sh]*scvf.global_grad(sh)[d] * vStdVel[ip][d];
    }
    
    number vNormStdVelPerConvLen[numIp];
    
    if(!bStokes)
        for(size_t ip = 0; ip < numIp; ++ip)
        {
            const number norm = VecTwoNorm(vStdVel[ip]);
            vNormStdVelPerConvLen[ip] = norm / upwind_conv_length(ip);
        }

    //    Find out if upwinded velocities depend on other ip velocities. In that case
    //    we have to solve a matrix system. Else the system is diagonal and we can
    //    compute the inverse directly

    //    diagonal case (i.e. upwind vel depend only on corner vel or no upwind)
    if(bStokes || !non_zero_shape_ip())
    {
    //    Loop integration points
        for(size_t ip = 0; ip < numIp; ++ip)
        {
        //    get SubControlVolumeFace
            const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

        //    First, we compute the contributions to the diagonal
        //    Note: - There is no contribution of the upwind vel to the diagonal
        //            in this case, only for non-diag problems
        //          - The diag does not depend on the dimension

        //    the diagonal entry
            number diag = vViscoPerDiffLenSq[ip];


        //    Convective Term  (no convective terms in the Stokes eq.)
            if (! bStokes)
                diag += vNormStdVelPerConvLen[ip];
            number stdVelCoeff=diag;
            
            if(pvCornerValueOldTime != NULL)
                diag += 1.0 / dt;


        //     Loop components of velocity
            for(int d = 0; d < dim; d++)
            {
            //    Now, we can assemble the rhs. This rhs is assembled by all
            //    terms, that are non-dependent on the ip vel.
            //    Note, that we can compute the stab_shapes on the fly when setting
            //    up the system.

            //    Source
                number rhs = 0.0;
                if(pSource != NULL)
                    rhs = (  ( density[ip]-density_ref) / density[ip]) * (*pSource)[ip][d];
                
                
                if(pvCornerValueOldTime != NULL)
                {
                //    interpolate old time step
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPVel += scvf.shape(sh) *  (*pvCornerValueOldTime)(d, sh) ;

                //    add to rhs
                    rhs += oldIPVel /  dt ;
                //   add density time derivative part from the continuity equation
                    //rhs -= vStdVel[ip][d] * (density[ip]-(*pDensity_OLD_SCVF)[ip]) /  dt ;
                }


            //    loop shape functions
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                //    Diffusion part
                    number sumVel = stdVelCoeff * scvf.shape(k);
                    
                   /* if(pvCornerValueOldTime != NULL)
                    {
                    //    interpolate old time step
 
                        sumVel -= scvf.shape(k)  / dt ;


                    }*/


                //    Convective term (no convective terms in the Stokes eq.)
                    if (! bStokes)
                    {
                        sumVel += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);

                    }

                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        if(d2 == d) continue;

                        sumVel -= vStdVel[ip][d2] * (scvf.global_grad(k))[d2];
                    }

                //    Add to rhs
                    rhs += sumVel * vCornerValue(d, k);

                //    set stab shape
                    stab_shape_vel(ip, d, d, k) = sumVel / diag;

                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        if(d2 == d) continue;

                        const number sumVel2 = vStdVel[ip][d] * (scvf.global_grad(k))[d2];

                        rhs += sumVel2 * vCornerValue(d2, k);

                        stab_shape_vel(ip, d, d2, k) = sumVel2 / diag;
                    }

                //    Pressure part
                    const number sumP = -1.0 * (scvf.global_grad(k))[d]  / density[ip];

                //    Add to rhs
                    rhs += sumP * vCornerValue(_P_, k);

                //    set stab shape
                    stab_shape_p(ip, d, k) = sumP / diag;
                }

            //    Finally, the can invert this row
                stab_vel(ip)[d] = rhs / diag;
            }
        }
    }
    /// need to solve system
    else
    {
        UG_THROW("Not implemented for ip velocities depending on other ip.");

    } // end switch for non-diag
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
       const DataImport<number, dim>& densitySCV,
       const number pressure_jump[],
       const MathVector<dim> SlipVel[],
       const number jump_shape[],
       const MathVector<dim> normal[],
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt,
       const number density_ref,
       const bool multiphase,
       const bool phase_2[],
       const number theta,
       number** SCVFinterShape)
{
    if(multiphase)
        UG_THROW("Pressure Jump Not implemented for NoStabiliation stabilization.");
//    abbreviation for pressure
    static const size_t _P_ = dim;

//    Some constants
    static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
    static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;

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
       const DataImport<number, dim>& densitySCV,
       const number pressure_jump[],
       const MathVector<dim> SlipVel[],
       const number jump_shape[],
       const MathVector<dim> normal[],
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt, 
       const number density_ref,
       const bool multiphase,
       const bool phase_2[],
       const number theta,
       number** SCVFinterShape)
{
    if(multiphase)
        UG_THROW("Pressure Jump Not implemented for  stabilization.");
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
template class NavierStokesVISCOSITY_Stabilization<2>;
template class NavierStokesKARIMIANStabilization<2>;

template SmartPtr<INavierStokesSRFV1Stabilization<2> >CreateNavierStokesStabilization<2>(const std::string& name);
#endif
#ifdef UG_DIM_3
template class INavierStokesFV1Stabilization<3>;
template class INavierStokesSRFV1Stabilization<3>;
template class NavierStokesFIELDSStabilization<3>;
template class NavierStokesFLOWStabilization<3>;
template class NavierStokesFIELDS_2_Stabilization<3>;
template class NavierStokesVISCOSITY_Stabilization<3>;
template class NavierStokesKARIMIANStabilization<3>;

template SmartPtr<INavierStokesSRFV1Stabilization<3> >CreateNavierStokesStabilization<3>(const std::string& name);
#endif

} // namespace NavierStokes
} // end namespace ug
