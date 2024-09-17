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
       const DataImport<number, dim>* pDensity_OLD_SCVF,
       const DataImport<number, dim>* pDensity_OLD_SCV,
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt)
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
				if(pSource != NULL)
					rhs = (*pSource)[ip][d];

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
					f[ip] = (*pSource)[ip][d];

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
       const DataImport<number, dim>* pDensity_OLD_SCVF,
       const DataImport<number, dim>* pDensity_OLD_SCV,
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt)
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
				if(pSource != NULL)
					rhs = (*pSource)[ip][d];

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
					f[ip] += (*pSource)[ip][d];

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
       const DataImport<number, dim>* pDensity_OLD_SCVF,
       const DataImport<number, dim>* pDensity_OLD_SCV,
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt)
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

//    cache values
    number vViscoPerDiffLenSq[numIp];
    number densityOLD_SCVF[numIp];
    for(size_t ip = 0; ip < numIp; ++ip)
    {
        vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);
        densityOLD_SCVF[ip] = (*pDensity_OLD_SCVF)[ip];
    }
    number densityOLD_SCV[numSh];
    for(size_t sh = 0; sh < numSh; ++sh)
        densityOLD_SCV[sh] = (*pDensity_OLD_SCV)[sh];


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
                //diag += 1./dt + (density[ip] - densityOLD_SCVF[ip]) / ( dt * density[ip]);

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
                    rhs = (*pSource)[ip][d];

            //    Time
                if(pvCornerValueOldTime != NULL)
                {
                //    interpolate old time step
                    number oldIPMomentum = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPMomentum += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh) ;
                        //oldIPMomentum += densityOLD_SCV[sh] * scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);

                //    add to rhs
                    rhs += oldIPMomentum / dt;
                    //rhs += oldIPMomentum / ( dt * density[ip]);
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
    //     For the FIELDS stabilization, there is no connection between the
    //    velocity components. Thus we can solve a system of size=numIP for
    //    each component of the velocity separately. This results in smaller
    //    matrices, that we have to invert.

    //    First, we have to assemble the Matrix, that includes all connections
    //    between the ip velocity component. Note, that in this case the
    //    matrix is non-diagonal and we must invert it.
    //    The Matrix is the same for all dim-components, thus we assemble it only
    //    once

    //    size of the system
        static const size_t N = numIp;

    //    a fixed size matrix
        DenseMatrix< FixedArray2<number, N, N> > mat;

    //    reset all values of the matrix to zero
        mat = 0.0;

    //    Loop integration points
        for(size_t ip = 0; ip < numIp; ++ip)
        {
        //    Time part
            if(pvCornerValueOldTime != NULL)
                mat(ip, ip) += 1./dt;

        //    Diffusion part
            mat(ip, ip) += vViscoPerDiffLenSq[ip];

        //    cache this value
            const number scale = vNormStdVelPerConvLen[ip];

        //    Convective Term (standard)
            mat(ip, ip) += scale;

        //    Convective Term by upwind
            for(size_t ip2 = 0; ip2 < numIp; ++ip2)
                mat(ip, ip2) -= upwind_shape_ip(ip, ip2) * scale;
        }

    //    we now create a matrix, where we store the inverse matrix
        typename block_traits<DenseMatrix< FixedArray2<number, N, N> > >::inverse_type inv;

    //    get the inverse
        if(!GetInverse(inv, mat))
            UG_THROW("Could not compute inverse.");

    //    loop dimensions (i.e. components of the velocity)
        for(int d = 0; d < dim; ++d)
        {
        //    create two vectors
            DenseVector< FixedArray1<number, N> > contVel[numSh];
            DenseVector< FixedArray1<number, N> > contP[numSh];

        //    Now, we can create several vector that describes the contribution of the
        //    corner velocities and the corner pressure. For each of this contribution
        //    components, we will apply the inverted matrix to get the stab_shapes

        //    Loop integration points
            for(size_t ip = 0; ip < numIp; ++ip)
            {
            //    get SubControlVolumeFace
                const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

            //    loop shape functions
                for(size_t k = 0; k < numSh; ++k)
                {
                //    Diffusion part
                    contVel[k][ip] = vViscoPerDiffLenSq[ip]    * scvf.shape(k);

                //    Convection part
                    contVel[k][ip] += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);

                //    Pressure part
                    contP[k][ip] = -1.0 * (scvf.global_grad(k))[d] / density[ip];
                }
            }

        //    solution vector
            DenseVector< FixedArray1<number, N> > xVel, xP;

        //    compute all stab_shapes
            for(size_t k = 0; k < numSh; ++k)
            {
            //    apply for vel stab_shape
                MatMult(xVel, 1.0, inv, contVel[k]);

            //    apply for pressure stab_shape
                MatMult(xP, 1.0, inv, contP[k]);

            //    write values in data structure
                //\todo: can we optimize this, e.g. without copy?
                for(size_t ip = 0; ip < numIp; ++ip)
                {
                //    write stab_shape for vel
                    stab_shape_vel(ip, d, d, k) = xVel[ip];

                //    write stab_shape for pressure
                    stab_shape_p(ip, d, k) = xP[ip];
                }
            }

        //    Finally, we can compute the values of the stabilized velocity for each
        //    integration point

        //    vector of all contributions
            DenseVector< FixedArray1<number, N> > f;

        //    Loop integration points
            for(size_t ip = 0; ip < numIp; ++ip)
            {
            //    get SubControlVolumeFace
                const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

            //    Source
                f[ip] = 0.0;
                if(pSource != NULL)
                    f[ip] = (*pSource)[ip][d];

            //    Time
                if(pvCornerValueOldTime != NULL)
                {
                //    interpolate old time step
                //    \todo: Is this ok? Or do we need the old stabilized vel ?
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);

                //    add to rhs
                    f[ip] += oldIPVel / dt;
                }
            }

        //    sum up all contributions of vel and p for rhs.
            for(size_t k = 0; k < numSh; ++k)
            {
            //    add velocity contribution
                VecScaleAdd(f, 1.0, f, vCornerValue(d, k), contVel[k]);

            //    add pressure contribution
                VecScaleAdd(f, 1.0, f, vCornerValue(_P_, k), contP[k]);
            }

        //    invert the system for all contributions
            DenseVector< FixedArray1<number, N> > x;
            MatMult(x, 1.0, inv, f);

        //    write values in data structure
            //\todo: can we optimize this, e.g. without copy?
            for(size_t ip = 0; ip < numIp; ++ip)
            {
            //    write stab_shape for vel
                stab_vel(ip)[d] = x[ip];
            }
        } // end dim loop

    } // end switch for non-diag
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
NavierStokesFLOW_2_Stabilization<TDim>::
update(const FV1Geometry<TElem, dim>* geo,
       const LocalVector& vCornerValue,
       const MathVector<dim> vStdVel[],
       const bool bStokes,
       const DataImport<number, dim>& kinVisco,
       const DataImport<number, dim>& kinViscoSCV,
       const DataImport<number, dim>& density,
       const DataImport<number, dim>& densitySCV,
       const DataImport<number, dim>* pDensity_OLD_SCVF,
       const DataImport<number, dim>* pDensity_OLD_SCV,
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt)
{
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

//    cache values
    number vViscoPerDiffLenSq[numIp];
    number densityOLD_SCVF[numIp];
    MathVector<dim> viscGrad[numIp];
    for(size_t ip = 0; ip < numIp; ++ip)
    {
        vViscoPerDiffLenSq[ip] = kinVisco[ip] * diff_length_sq_inv(ip);
        densityOLD_SCVF[ip] = (*pDensity_OLD_SCVF)[ip];
        const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);
        VecSet(viscGrad[ip],0);
        for(size_t sh = 0; sh < numSh; ++sh)
            for(int d = 0; d < dim; ++d)
                viscGrad[ip][d] += kinViscoSCV[sh]*densitySCV[sh]*scvf.global_grad(sh)[d];
    }

    number densityOLD_SCV[numSh];
    for(size_t sh = 0; sh < numSh; ++sh)
        densityOLD_SCV[sh] = (*pDensity_OLD_SCV)[sh];
    
    number vNormStdVelPerConvLen[numIp];
    number vNormStdVelPerDownLen[numIp];
    
    if(!bStokes)
        for(size_t ip = 0; ip < numIp; ++ip)
        {
            const number norm = VecTwoNorm(vStdVel[ip]);
            vNormStdVelPerConvLen[ip] = norm / upwind_conv_length(ip);
            vNormStdVelPerDownLen[ip] = norm / (downwind_conv_length(ip) + upwind_conv_length(ip));
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

        //    Time part
        if(pvCornerValueOldTime != NULL)
        {
            diag += 1./dt;
            //diag += 1./dt + (density[ip] - densityOLD_SCVF[ip]) / ( dt * density[ip]);
        }

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
                    rhs = (*pSource)[ip][d];

            //    Time
                if(pvCornerValueOldTime != NULL)
                {
                //    interpolate old time step
                    number oldIPMomentum = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                    {
                        oldIPMomentum += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                        //oldIPMomentum += densityOLD_SCV[sh] * scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);
                    }

                //    add to rhs
                    rhs += oldIPMomentum / dt;
                    //rhs += oldIPMomentum / ( dt * density[ip]);
                }

            //    loop shape functions
                for(size_t k = 0; k < scvf.num_sh(); ++k)
                {
                //    Diffusion part
                    number sumVel = vViscoPerDiffLenSq[ip] * scvf.shape(k);
                    
                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        //sumVel += (viscGrad[ip][d2] / density[ip]) * (scvf.global_grad(k))[d2];
                    }

                //    Convective term (no convective terms in the Stokes eq.)
                    if (! bStokes)
                    {
                        //sumVel += (densitySCV[k]/density[ip]) * vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);
                        sumVel += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);

                        //sumVel += (densitySCV[k]/density[ip]) * vNormStdVelPerDownLen[ip] * (downwind_shape_sh(ip, k) - upwind_shape_sh(ip, k));
                        sumVel += vNormStdVelPerDownLen[ip] * (downwind_shape_sh(ip, k) - upwind_shape_sh(ip, k));
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
    //     For the FLOW stabilization, there is no connection between the
    //    velocity components. Thus we can solve a system of size=numIP for
    //    each component of the velocity separately. This results in smaller
    //    matrices, that we have to invert.

    //    First, we have to assemble the Matrix, that includes all connections
    //    between the ip velocity component. Note, that in this case the
    //    matrix is non-diagonal and we must invert it.
    //    The Matrix is the same for all dim-components. Thus we invert it only
    //    once.

    //    size of the system
        static const size_t N = numIp;

    //    a fixed size matrix
        DenseMatrix< FixedArray2<number, N, N> > mat;

    //    reset all values of the matrix to zero
        mat = 0.0;

    //    Loop integration points
        for(size_t ip = 0; ip < numIp; ++ip)
        {
        //    Time part
            if(pvCornerValueOldTime != NULL)
                mat(ip, ip) += 1./dt;

        //    Diffusion part
            mat(ip, ip) += vViscoPerDiffLenSq[ip];

        //    Convective Term (standard)
            mat(ip, ip) += vNormStdVelPerConvLen[ip];

            for(size_t ip2 = 0; ip2 < numIp; ++ip2)
            {
            //    Convective Term by upwind
                mat(ip, ip2) -= upwind_shape_ip(ip, ip2) * vNormStdVelPerConvLen[ip];

            //    correction of divergence error
                mat(ip,ip2) += vNormStdVelPerDownLen[ip] *
                                (upwind_shape_ip(ip, ip2) - downwind_shape_ip(ip, ip2));
            }
        }

    //    we now create a matrix, where we store the inverse matrix
        typename block_traits<DenseMatrix< FixedArray2<number, N, N> > >::inverse_type inv;

    //    get the inverse
        if(!GetInverse(inv, mat))
            UG_THROW("Could not compute inverse.");


    //    create vectors
        DenseVector< FixedArray1<number, N> > contVel[dim][numSh];
        DenseVector< FixedArray1<number, N> > contP[numSh];
        DenseVector< FixedArray1<number, N> > xP;
        DenseVector< FixedArray1<number, N> > xVel;

    //    Now, we can create several vector that describes the contribution of the
    //    corner velocities. For each of this contribution
    //    components, we will apply the inverted matrix to get the stab_shapes

    //    Loop integration points
        for(int d = 0; d < dim; ++d)
        {
            for(size_t ip = 0; ip < numIp; ++ip)
            {
            //    get SubControlVolumeFace
                const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

            //    loop shape functions
                for(size_t k = 0; k < numSh; ++k)
                {
                //    Pressure part
                    contP[k][ip] = -1.0 * (scvf.global_grad(k))[d] / density[ip];

                //    Diffusion part
                    contVel[d][k][ip] = vViscoPerDiffLenSq[ip] * scvf.shape(k);

                //    Convection part
                    contVel[d][k][ip] += vNormStdVelPerConvLen[ip] * upwind_shape_sh(ip, k);

                //    terms for correction of divergence error
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

        //    compute all stab_shapes
            for(size_t k = 0; k < numSh; ++k)
            {
            //    apply for pressure stab_shape
                MatMult(xP, 1.0, inv, contP[k]);

            //    write stab_shape for pressure
                //\todo: can we optimize this, e.g. without copy?
                for(size_t ip = 0; ip < numIp; ++ip)
                        stab_shape_p(ip, d, k) = xP[ip];

            //    compute vel stab_shapes
                for(int d2 = 0; d2 < dim; ++d2)
                {
                //    apply for vel stab_shape
                    MatMult(xVel, 1.0, inv, contVel[d2][k]);

                //    write stab_shape for vel
                    //\todo: can we optimize this, e.g. without copy?
                    for(size_t ip = 0; ip < numIp; ++ip)
                        stab_shape_vel(ip, d, d2, k) = xVel[ip];
                }
            }

        //    Finally, we can compute the values of the stabilized velocity for each
        //    integration point

        //    vector of all contributions
            DenseVector< FixedArray1<number, N> > f;

        //    sum up all contributions of vel and p for rhs.
            f = 0.0;
            for(size_t k = 0; k < numSh; ++k)
            {
            //    add velocity contribution
                for(int d2 = 0; d2 < dim; ++d2)
                    VecScaleAdd(f, 1.0, f, vCornerValue(d2, k), contVel[d2][k]);

            //    add pressure contribution
                VecScaleAdd(f, 1.0, f, vCornerValue(_P_, k), contP[k]);
            }

        //    Loop integration points
            for(size_t ip = 0; ip < numIp; ++ip)
            {
            //    get SubControlVolumeFace
                const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

            //    Source
                if(pSource != NULL)
                    f[ip] += (*pSource)[ip][d];

            //    Time
                if(pvCornerValueOldTime != NULL)
                {
                //    interpolate old time step
                //    \todo: Is this ok? Or do we need the old stabilized vel ?
                    number oldIPVel = 0.0;
                    for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                        oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);

                //    add to rhs
                    f[ip] += oldIPVel / dt;
                }
            }

        //    invert the system for all contributions
            DenseVector< FixedArray1<number, N> > x;
            MatMult(x, 1.0, inv, f);

        //    write values in data structure
            //\todo: can we optimize this, e.g. without copy?
            for(size_t ip = 0; ip < numIp; ++ip)
            {
            //    write stab_shape for vel
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
       const DataImport<number, dim>* pDensity_OLD_SCVF,
       const DataImport<number, dim>* pDensity_OLD_SCV,
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt)
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

template SmartPtr<INavierStokesSRFV1Stabilization<2> >CreateNavierStokesStabilization<2>(const std::string& name);
#endif
#ifdef UG_DIM_3
template class INavierStokesFV1Stabilization<3>;
template class INavierStokesSRFV1Stabilization<3>;
template class NavierStokesFIELDSStabilization<3>;
template class NavierStokesFLOWStabilization<3>;
template class NavierStokesFIELDS_2_Stabilization<3>;
template class NavierStokesFLOW_2_Stabilization<3>;

template SmartPtr<INavierStokesSRFV1Stabilization<3> >CreateNavierStokesStabilization<3>(const std::string& name);
#endif

} // namespace NavierStokes
} // end namespace ug
