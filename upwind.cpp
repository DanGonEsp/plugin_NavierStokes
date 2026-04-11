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

// for minimum
#include <limits>
#include <algorithm>
#include <locale>

// function space, reference element
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "common/util/provider.h"

#include "upwind.h"

namespace ug {
namespace NavierStokes{

/////////////////////////////////////////////////////////////////////////////
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
template <typename TElem>
void
NavierStokesNoUpwind<dim>::
compute(const FV1Geometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
//	set shapes
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
		//	set upwind shape
			vUpShapeSh[ip][sh] = scvf.shape(sh);
		}

	//	compute convection length
	//  \todo: (optional) A convection length is not really defined for no upwind.
	//	       but in the computation of a stabilization the term cancels, so
	//   	   we only have to ensure that the conv_lengh is non-zero
        vConvLength[ip] = 1.0;
	}
}

template <int dim>
template <typename TElem>
void
NavierStokesNoUpwind<dim>::
compute(const CRFVGeometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
//	set shapes
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename CRFVGeometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
		//	set upwind shape
			vUpShapeSh[ip][sh] = scvf.shape(sh);
		}
	}
}

template <int dim>
template <typename TElem>
void
NavierStokesNoUpwind<dim>::
compute(const HCRFVGeometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
//	set shapes
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename HCRFVGeometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
		//	set upwind shape
			vUpShapeSh[ip][sh] = scvf.shape(sh);
		}
	}
}
/////////////////////////////////////////////////////////////////////////////
// Full Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
template <typename TElem>
void
NavierStokesFullUpwind<dim>::
compute(const FV1Geometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
//	two help vectors
	MathVector<dim> dist;

// 	get corners of elem
    const MathVector<dim>* corners = geo->corners();

// 	set shapes
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

    //  reset shapes to zero for all IPs
        for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
        	vUpShapeSh[ip][sh]=0.0;

    // 	switch upwind
        const number flux = VecDot(scvf.normal(), vIPVel[ip]);
        if(flux >= 0.0)
        {
        	vUpShapeSh[ip][scvf.from()] = 1.0;
        	vConvLength[ip] = VecDistance(scvf.global_ip(), corners[scvf.from()]);
        }
        else
        {
        	vUpShapeSh[ip][scvf.to()] = 1.0;
        	vConvLength[ip] = VecDistance(scvf.global_ip(), corners[scvf.to()]);
        }
    }
}

template <int dim>
template <typename TElem>
void
NavierStokesFullUpwind<dim>::
compute(const CRFVGeometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
//	two help vectors
	MathVector<dim> dist;

// 	get corners of elem
    const MathVector<dim>* elementfaces = geo->scv_global_ips();

// 	set shapes
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
    //	get SubControlVolumeFace
		const typename CRFVGeometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

    //  reset shapes to zero for all IPs
        for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
        	vUpShapeSh[ip][sh]=0.0;

    // 	switch upwind
        const number flux = VecDot(scvf.normal(), vIPVel[ip]);
        if(flux > 0.0)
        {
        	vUpShapeSh[ip][scvf.from()] = 1.0;
        	vConvLength[ip] = VecDistance(scvf.global_ip(), elementfaces[scvf.from()]);
        }
        else
        {
        	vUpShapeSh[ip][scvf.to()] = 1.0;
        	vConvLength[ip] = VecDistance(scvf.global_ip(), elementfaces[scvf.to()]);
        }
    }
}

// upwind shape in constraining dof is shape in constrained dof
template <int dim>
template <typename TElem>
void
NavierStokesFullUpwind<dim>::
compute(const HCRFVGeometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
	MathVector<dim> dist;

    if (geo->num_constrained_dofs()>0){
    	std::vector<size_t> constrainedShape(geo->num_scv()+geo->num_constrained_dofs());
    	for (size_t i=0;i<geo->num_sh();i++) constrainedShape[i] = i;
    	for (size_t i=0;i<geo->num_constrained_dofs();i++){
    		const typename HCRFVGeometry<TElem, dim>::CONSTRAINED_DOF& cd = geo->constrained_dof(i);
    		const size_t index = cd.index();
   			for (size_t j=0;j<cd.num_constraining_dofs();j++){
    			constrainedShape[cd.constraining_dofs_index(j)]=index;
    		}
    	}
    	// 	set shapes
    	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
    	{
    	    //	get SubControlVolumeFace
    		const typename HCRFVGeometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

    	    //  reset shapes to zero for all IPs
    	    for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
    	       	vUpShapeSh[ip][sh]=0.0;

    	    // 	switch upwind
    	    const number flux = VecDot(scvf.normal(), vIPVel[ip]);
    	    if(flux > 0.0)
    	    {
    	      	vUpShapeSh[ip][constrainedShape[scvf.from()]] = 1.0;
    	    }
    	    else
    	    {
    	      	vUpShapeSh[ip][constrainedShape[scvf.to()]] = 1.0;
    	    }
    	}
    }

// 	set shapes
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
    //	get SubControlVolumeFace
		const typename HCRFVGeometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

    //  reset shapes to zero for all IPs
        for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
        	vUpShapeSh[ip][sh]=0.0;

    // 	switch upwind
        const number flux = VecDot(scvf.normal(), vIPVel[ip]);
        if(flux > 0.0)
        {
        	vUpShapeSh[ip][scvf.from()] = 1.0;
        }
        else
        {
        	vUpShapeSh[ip][scvf.to()] = 1.0;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////
// Weighted Upwind on Crouzeix-Raviart type elements
// upwinding between full and no upwind
// shapes computed as (1-m_weight)*no_upwind_shape + m_weight*full_upwind_shape
//////////////////////////////////////////////////////////////////////////////////
/*
template <int dim>
template <typename TElem>
void
NavierStokesWeightedUpwind<dim>::
compute(const CRFVGeometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
	MathVector<dim> dist;

// 	set full upwind shapes
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
    //	get SubControlVolumeFace
		const typename CRFVGeometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

    //  reset shapes to zero for all IPs
        for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
        	vUpShapeSh[ip][sh]=0.0;

	// 	compute upwind shapes
        const number flux = VecDot(scvf.normal(), vIPVel[ip]);
        if(flux > 0.0)
        {
        	vUpShapeSh[ip][scvf.from()] = this->m_weight;
        }
        else
        {
        	vUpShapeSh[ip][scvf.to()] = this->m_weight;
        }

	//	add no upwind shapes
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			//	set upwind shape
			vUpShapeSh[ip][sh] += (1-this->m_weight) * scvf.shape(sh);
		}
    }
}
*/
/////////////////////////////////////////////////////////////////////////////
// Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

/// computes the closest node to a elem side ray intersection
template <typename TRefElem, int TWorldDim>
void GetNodeNextToCut(size_t& coOut,
                      const MathVector<TWorldDim>& IP,
                      const MathVector<TWorldDim>& IPVel,
                      const MathVector<TWorldDim>* vCornerCoords)
{
//	help variables
	size_t side = 0;
	MathVector<TWorldDim> globalIntersection;
	MathVector<TRefElem::dim> localIntersection;

//	compute intersection of ray in direction of ip velocity with elem side
//	we search the ray only in upwind direction
	if(!ElementSideRayIntersection<TRefElem, TWorldDim>
		(	side, globalIntersection, localIntersection,
			IP, IPVel, false /* i.e. search upwind */, vCornerCoords))
		UG_THROW("GetNodeNextToCut: Cannot find cut side.");

//	get reference element
	static const TRefElem& rRefElem = Provider<TRefElem>::get();
	const int dim = TRefElem::dim;

// 	reset minimum
	number min = std::numeric_limits<number>::max();

// 	loop corners of side
	for(size_t i = 0; i < rRefElem.num(dim-1, side, 0); ++i)
	{
	// 	get corner
		const size_t co = rRefElem.id(dim-1, side, 0, i);

	// 	Compute Distance to intersection
		number dist = VecDistanceSq(globalIntersection, vCornerCoords[co]);

	// 	if distance is smaller, choose this node
		if(dist < min)
		{
			min = dist;
			coOut = co;
		}
	}
}

template <int dim>
template <typename TElem>
void
NavierStokesSkewedUpwind<dim>::
compute(const FV1Geometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	loop all scvf
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

	// 	reset shapes to zero
 		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
 			vUpShapeSh[ip][sh] = 0.0;

 	//	if the velocity is zero, there will be no possibility to find the
 	//	cutted side. In this case we have no velocity and therefore there is
 	//	no convection. We set all upwind shapes to zero.
 		if(VecTwoNorm(vIPVel[ip]) < 1e-14) {
 			//  \todo: (optional) A convection length is not really defined.
 			//	       but in the computation of a stabilization the term cancels, so
 			//   	   we only have to ensure that the conv_lengh is non-zero
			vConvLength[ip] = 1.0;
			continue;
 		}

	// 	upwind corner
		size_t sh = 0;

	// 	find upwind node
		try{
			GetNodeNextToCut<typename FV1Geometry<TElem, dim>::ref_elem_type, dim>
					(sh, scvf.global_ip(), vIPVel[ip], vCornerCoords);
		}UG_CATCH_THROW("GetSkewedUpwindShapes: Cannot find upwind node.");

	// 	set upwind corner
		vUpShapeSh[ip][sh] = 1.0;

	//	compute convection length
	    vConvLength[ip] = VecDistance(scvf.global_ip(), vCornerCoords[sh]);
	}
}

template <int dim>
template <typename TElem>
void
NavierStokesSkewedUpwind<dim>::
compute(const CRFVGeometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	loop all scvf
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
    //	get SubControlVolumeFace
		const typename CRFVGeometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

		size_t num_sh =  scvf.num_sh();

	// 	reset shapes to zero
 		for(size_t sh = 0; sh < num_sh; ++sh)
 			vUpShapeSh[ip][sh] = 0.0;

 	//	if the velocity is zero, there will be no possibility to find the
 	//	cutted side. In this case we have no velocity and therefore there is
 	//	no convection. We set all upwind shapes to zero.
 		if(VecTwoNorm(vIPVel[ip]) < 1e-14) continue;

 	// 	side and intersection vectors
 		static const int refDim = CRFVGeometry<TElem, dim>::dim;
 		size_t side = 0;
 		MathVector<dim> globalIntersection;
 		MathVector<refDim> localIntersection;

 	// 	find local intersection and side
 		try{
 			ElementSideRayIntersection<typename CRFVGeometry<TElem, dim>::ref_elem_type, dim>
 			(	side, globalIntersection, localIntersection,
 				scvf.global_ip(), vIPVel[ip], false /* search upwind */, vCornerCoords);
 		}UG_CATCH_THROW("GetSkewedUpwindShapes: Cannot find cut side.");

 	// 	get linear trial space
 		static const ReferenceObjectID roid = reference_element_traits<TElem>::reference_element_type::REFERENCE_OBJECT_ID;
 		const LocalShapeFunctionSet<refDim>& TrialSpace =
 				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

 	// 	get Reference Element
 		typedef typename CRFVGeometry<TElem, dim>::ref_elem_type ref_elem_type;

 		number max = -1000;
 		size_t maxind=0;

 	// 	loop shape functions
 		for(size_t sh=0;sh < num_sh;sh++){
 			number shape = TrialSpace.shape(sh, localIntersection);
 			if (shape>max){
 				max=shape;
 				maxind = sh;
 			}
 		}
 		vUpShapeSh[ip][maxind] = 1;

 	//	compute conv length
		vConvLength[ip] = VecDistance(scvf.global_ip(), geo->scv(maxind).global_ip());
	}
}

/////////////////////////////////////////////////////////////////////////////
// Linear Profile Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
template <typename TElem>
void
NavierStokesLinearProfileSkewedUpwind<dim>::
compute(const FV1Geometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	loop all scvf
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

	// 	reset shapes to zero
 		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
 			vUpShapeSh[ip][sh] = 0.0;

 	//	if the velocity is zero, there will be no possibility to find the
 	//	cutted side. In this case we have no velocity and therefore there is
 	//	no convection. We set all upwind shapes to zero.
 		if(VecTwoNorm(vIPVel[ip]) < 1e-14) {
 			//  \todo: (optional) A convection length is not really defined.
 			//	       but in the computation of a stabilization the term cancels, so
 			//   	   we only have to ensure that the conv_lengh is non-zero
			vConvLength[ip] = 1.0;
			continue;
 		}

 	// 	side and intersection vectors
 		static const int refDim = FV1Geometry<TElem, dim>::dim;
 		size_t side = 0;
 		MathVector<dim> globalIntersection;
 		MathVector<refDim> localIntersection;

 	// 	find local intersection and side
 		try{
 			ElementSideRayIntersection<typename FV1Geometry<TElem, dim>::ref_elem_type, dim>
 			(	side, globalIntersection, localIntersection,
 				scvf.global_ip(), vIPVel[ip], false /* search upwind */, vCornerCoords);
 		}UG_CATCH_THROW("GetLinearProfileSkewedUpwindShapes: Cannot find cut side.");

 	// 	get linear trial space
 		static const ReferenceObjectID roid = reference_element_traits<TElem>::reference_element_type::REFERENCE_OBJECT_ID;
 		const LocalShapeFunctionSet<refDim>& TrialSpace =
 				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

 	// 	get Reference Element
 		typedef typename FV1Geometry<TElem, dim>::ref_elem_type ref_elem_type;
 		static const ref_elem_type& rRefElem
 			= Provider<ref_elem_type>::get();

 	// 	loop corners of side
 		for(size_t j = 0; j < rRefElem.num(dim-1, side, 0); ++j)
 		{
 		// 	get corner
 			const size_t co = rRefElem.id(dim-1, side, 0, j);

 		//	evaluate trial space
 			vUpShapeSh[ip][co] = TrialSpace.shape(co, localIntersection);
 		}

 	//	compute conv length
		vConvLength[ip] = VecDistance(scvf.global_ip(), globalIntersection);
	}
}

template <int dim>
template <typename TElem>
void
NavierStokesLinearProfileSkewedUpwind<dim>::
compute(const CRFVGeometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	loop all scvf
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
    //	get SubControlVolumeFace
		const typename CRFVGeometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

		size_t num_sh =  scvf.num_sh();

	// 	reset shapes to zero
 		for(size_t sh = 0; sh < num_sh; ++sh)
 			vUpShapeSh[ip][sh] = 0.0;

 	//	if the velocity is zero, there will be no possibility to find the
 	//	cutted side. In this case we have no velocity and therefore there is
 	//	no convection. We set all upwind shapes to zero.
 		if(VecTwoNorm(vIPVel[ip]) == 0.0) continue;

 	// 	side and intersection vectors
 		static const int refDim = CRFVGeometry<TElem, dim>::dim;
 		size_t side = 0;
 		MathVector<dim> globalIntersection;
 		MathVector<refDim> localIntersection;

 	// 	find local intersection and side
 		try{
 			ElementSideRayIntersection<typename CRFVGeometry<TElem, dim>::ref_elem_type, dim>
 			(	side, globalIntersection, localIntersection,
 				scvf.global_ip(), vIPVel[ip], false /* search upwind */, vCornerCoords);
 		}UG_CATCH_THROW("GetLinearProfileSkewedUpwindShapes: Cannot find cut side.");

 	// 	get linear trial space
 		static const ReferenceObjectID roid = reference_element_traits<TElem>::reference_element_type::REFERENCE_OBJECT_ID;
 		const LocalShapeFunctionSet<refDim>& TrialSpace =
 				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

 	// 	get Reference Element
 		typedef typename CRFVGeometry<TElem, dim>::ref_elem_type ref_elem_type;

  	// 	loop shape functions
 		for(size_t sh=0;sh < num_sh;sh++){
 			vUpShapeSh[ip][sh] = TrialSpace.shape(sh, localIntersection);
 		}

 	//	compute conv length
 		vConvLength[ip] = VecDistance(scvf.global_ip(), globalIntersection);
	}
}


/////////////////////////////////////////////////////////////////////////////
// Positive Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
template <typename TElem>
void
NavierStokesPositiveUpwind<dim>::
compute(const FV1Geometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{

//	1. Reset values and compute ip velocities and Compute mass fluxes at ip's

//	vector for flux values
	std::vector<number> vMassFlux(geo->num_scvf(), 0.0);
	std::vector<bool> vHasFlux(geo->num_scvf(), true);

//	loop all scvf
	const number eps = std::numeric_limits<number>::epsilon() * 10;
	size_t bNumNoFlux = 0;
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

	// 	reset shapes w.r.t corner value to zero
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			vUpShapeSh[ip][sh] = 0.0;

	// 	reset shapes w.r.t. ip value to zero and extract ip vel
		for(size_t j = 0; j < geo->num_scvf(); ++j)
			vUpShapeIp[ip][j] = 0.0;

		const number normSq = VecTwoNormSq(vIPVel[ip]);
		if(fabs(normSq) <= eps)
		{
			vUpShapeSh[ip][scvf.from()] = 0.5;
			vUpShapeSh[ip][scvf.to()] = 0.5;
			vHasFlux[ip] = false;
			bNumNoFlux++;
			continue;
		}

		vMassFlux[ip] = VecProd(vIPVel[ip], scvf.normal());
		const number vel = std::sqrt(normSq);
		const number len = VecTwoNorm(scvf.normal());
		if(fabs(vMassFlux[ip] / std::sqrt(vel*len))  <= eps)
		{
			vUpShapeSh[ip][scvf.from()] = 0.5;
			vUpShapeSh[ip][scvf.to()] = 0.5;
			vHasFlux[ip] = false;
			bNumNoFlux++;
			continue;
		}
	}


//	2. Handle each SCV separately
	if(bNumNoFlux != geo->num_scvf())
	{
		for(size_t sh = 0; sh < geo->num_sh(); ++sh)
		{
		//	reset inflow, outflow
			number m_in = 0, m_out = 0;
			std::vector<size_t> vSCVIP;
			std::vector<number> vFlux;

		//	loop subcontrol volume faces
			for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
			{
			//	if no flux skip
				if(!vHasFlux[ip]) continue;

			//	get SubControlVolumeFace
				const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

			//	if scvf is part of the scv, add fluxes
				if(scvf.from() == sh)
				{
				//	normal directed outwards
					vSCVIP.push_back(ip);
					vFlux.push_back( vMassFlux[ip] );
					m_in += -1.0 * std::min(vMassFlux[ip], 0.0);
					m_out += std::max(vMassFlux[ip], 0.0);
				}
				else if (scvf.to() == sh)
				{
				//	normal directed inwards
					vSCVIP.push_back(ip);
					vFlux.push_back( -1.0 * vMassFlux[ip] );
					m_in += -1.0 * std::min(-1.0 *  vMassFlux[ip], 0.0);
					m_out += std::max(-1.0 * vMassFlux[ip], 0.0);
				}
			}

		//	compute F
			number F = std::max(m_in, m_out);

		//	set shapes
			for(size_t i = 0; i < vSCVIP.size(); ++i)
			{
				if(vFlux[i] > 0)
				{
					number sum = 0.0;
					for(size_t j = 0; j < vSCVIP.size(); ++j)
					{
						if(vFlux[j] < 0)
						{
						//	set ip shapes
							sum += vUpShapeIp[vSCVIP[i]][vSCVIP[j]] = -1.0 * vFlux[j] / F;
						}
					}
				//	set nodal shapes
					vUpShapeSh[vSCVIP[i]][sh] = 1.0 - sum;
				}
			}
		}
	}

//	3. compute convection length

// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	compute upwind point
	MathVector<dim> upPos;
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

	//	reset upwind point
		VecSet(upPos, 0.0);

	//	sum up contributions
        for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
        	VecScaleAppend(upPos, vUpShapeSh[ip][sh], vCornerCoords[sh]);
        for (size_t j = 0; j < geo->num_scvf(); ++j)
        	VecScaleAppend(upPos, vUpShapeIp[ip][j], geo->scvf(j).global_ip());

    //	save convection length
		vConvLength[ip] = VecDistance(scvf.global_ip(), upPos);
	}
}

/////////////////////////////////////////////////////////////////////////////
// Regular Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
template <typename TElem>
void
NavierStokesRegularUpwind<dim>::
compute(const FV1Geometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
	const number eps = std::numeric_limits<number>::epsilon() * 10;
	const int nc = geo->num_scv();

	// ONLY 2D IMPLEMENTED
	if(dim != 2) UG_THROW("RegularUpwind only implemented for 2d.");

//	loop all scvf
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

	// 	reset shapes to zero
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			vUpShapeSh[ip][sh] = 0.0;
		for(size_t ip2 = 0; ip2 < scvf.num_sh(); ++ip2)
			vUpShapeIp[ip][ip2] = 0.0;

	//	if the velocity is zero, there will be no possibility to find the
	//	cutted side. In this case we have no velocity and therefore there is
	//	no convection. We set all upwind shapes to zero.
		const number normSq = VecTwoNormSq(vIPVel[ip]);
		if(fabs(normSq) <= eps) continue;

	//	flux over scvf
		number flux = VecProd(vIPVel[ip], scvf.normal());

	//	check if flux is very small (e.g. due to orthogonality of normal and
	//	flux direction)
	//	todo: Think about the small constants.
		if(fabs(flux) <= 100*eps)
		{
			// TODO: THIS IS 2D ONLY !!!!
			// the convection is parallel to the subcontrol volume surface
			flux = vIPVel[ip][0]*scvf.normal()[1] - vIPVel[ip][1]*scvf.normal()[0];
			if (flux>0)
			{
				// the velocity is pointing to the element boundary
				// take lin comb of pred and succ ip
				vUpShapeIp[ip][(ip+nc-1)%nc] = vUpShapeIp[ip][(ip+1)%nc] = 0.5;
			}
			else
			{
				// the velocity is pointing to the element midpoint
				// take lin comb of pred and succ node
				vUpShapeSh[ip][ip] = vUpShapeSh[ip][(ip+1)%nc] = 0.5;
			}
			continue;
	 	}

	//	get upwind scv
		int upwindSCV = -1;
		if(flux > 0.0) upwindSCV = scvf.from();
		else upwindSCV = scvf.to();

	//	get corresponding scv
		const typename FV1Geometry<TElem, dim>::SCV& scv = geo->scv(upwindSCV);

	// 	side and intersection vectors
		static const int refDim = FV1Geometry<TElem, dim>::dim;
		size_t side = 0;
		number lambda = 0.0;
		MathVector<dim> globalIntersection;
		MathVector<refDim> localIntersection;

		// TODO: THIS IS 2D ONLY !!!!
		//a) check for scvf intersection
		if(SCVFofSCVRayIntersection<refDim, dim>(side, lambda, globalIntersection, localIntersection,
		                            scvf.global_ip(), vIPVel[ip], false, scv.global_corners()))
		{
			int ipIntersect = (ip + 1) % nc;
			int ipIntersectOppose = (ip + nc - 1) % nc;
			if(flux > 0.0)
			{
				ipIntersect = (ip + nc - 1) % nc;
				ipIntersectOppose = (ip + 1) % nc;
			}

		// 	a.1) between two ip points
			if(lambda <= 0.5)
			{
				vUpShapeIp[ip][ipIntersectOppose] = lambda - 0.5;
				vUpShapeIp[ip][ipIntersect] = 1.0-(lambda-0.5);
			}

		// 	a.2) between ip point and elem side
			else
			{
			//	interpolation between corners and ip
				vUpShapeSh[ip][scvf.from()] = 0.5*2*(lambda-0.5);
				vUpShapeSh[ip][scvf.to()] = 0.5*2*(lambda-0.5);
				vUpShapeIp[ip][ipIntersect] = 1.0-2*(lambda-0.5);
			}

			continue;
		}
		//b) if not, on elem side intersection
		else
		{
			if(flux > 0.0)
			{
				switch (side)
				{
					case 0:
						// take linear profile between nodes on element side ip+1
						vUpShapeSh[ip][ip]		  = 1.0-0.5*lambda;
						vUpShapeSh[ip][(ip+1)%nc] = 0.5*lambda;
						break;
					case 3:
						// take linear profile between nodes on element side ip
						vUpShapeSh[ip][(ip+nc-1)%nc] = 1.0-0.5*(lambda+1.0);
						vUpShapeSh[ip][ip]			 = 0.5*(lambda+1.0);
						break;
					default:
						UG_THROW("This should not happen.");
				}
			}
			else
			{
				switch (side)
				{
					case 0:
						// take linear profile between nodes on element side ip+1
						vUpShapeSh[ip][(ip+1)%nc] = 1.0-0.5*lambda;
						vUpShapeSh[ip][(ip+2)%nc] = 0.5*lambda;
						break;
					case 3:
						// take linear profile between nodes on element side ip
						vUpShapeSh[ip][ip]		  = 1.0-0.5*(lambda+1.0);
						vUpShapeSh[ip][(ip+1)%nc] = 0.5*(lambda+1.0);
						break;
					default:
						UG_THROW("This should not happen.");
				}
			}
		}
	}

//	compute convection length

// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	compute upwind point
	MathVector<dim> upPos;
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

	//	reset upwind point
		VecSet(upPos, 0.0);

	//	sum up contributions
		for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(upPos, vUpShapeSh[ip][sh], vCornerCoords[sh]);
		for (size_t j = 0; j < geo->num_scvf(); ++j)
			VecScaleAppend(upPos, vUpShapeIp[ip][j], geo->scvf(j).global_ip());

	//	save convection length
		vConvLength[ip] = VecDistance(scvf.global_ip(), upPos);
	}
}

////////////////////////////////////////////////////////////////////////////////
//	explicit instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_2
template class NavierStokesNoUpwind<2>;
template class NavierStokesFullUpwind<2>;
//template class NavierStokesWeightedUpwind<2>;
template class NavierStokesSkewedUpwind<2>;
template class NavierStokesLinearProfileSkewedUpwind<2>;
template class NavierStokesPositiveUpwind<2>;
template class NavierStokesRegularUpwind<2>;
#endif

#ifdef UG_DIM_3
template class NavierStokesNoUpwind<3>;
template class NavierStokesFullUpwind<3>;
//template class NavierStokesWeightedUpwind<3>;
template class NavierStokesSkewedUpwind<3>;
template class NavierStokesLinearProfileSkewedUpwind<3>;
template class NavierStokesPositiveUpwind<3>;
template class NavierStokesRegularUpwind<3>;
#endif

} // namespace NavierStokes
} // end namespace ug
