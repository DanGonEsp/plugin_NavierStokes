/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__NAVIER_STOKES__UPWIND_INTERFACE__
#define __H__UG__NAVIER_STOKES__UPWIND_INTERFACE__

//#define UG_NSUPWIND_ASSERT(cond, exp)
// include define below to assert arrays used in stabilization
#define UG_NSUPWIND_ASSERT(cond, exp) UG_ASSERT((cond), exp)

#include <vector>

#include "lib_disc/common/local_algebra.h"
#include "lib_disc/spatial_disc/disc_util/fv_geom_base.h"
#include "lib_disc/spatial_disc/disc_util/fv_util.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

/////////////////////////////////////////////////////////////////////////////
// Interface for Upwinds on collocated grid (stabilization method)
/////////////////////////////////////////////////////////////////////////////


template <int dim>
class INavierStokesUpwind
{
	protected:
	///	used traits
		typedef fv1_dim_traits<dim, dim> traits;

	public:
	///	max number of SubControlVolumeFaces
		static const size_t maxNumSCVF = traits::maxNumSCVF + 10;

	/// max number of shape functions
		static const size_t maxNumSH = traits::maxNSH + 10;

	public:
	/// Abbreviation for own type
		typedef INavierStokesUpwind<dim> this_type;

	public:
	///	constructor
		INavierStokesUpwind()
			: m_numScvf(0), m_numSh(0),
			  m_bNonZeroShapeIp(true)
		{
			m_vComputeFunc.clear();
		}

	///	returns number of shapes
		size_t num_sh() const {return m_numSh;}

	///	returns number of sub control volume faces
		size_t num_scvf() const {return m_numScvf;}

	///	Convection Length
		number upwind_conv_length(size_t scvf) const
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			return m_vUpConvLength[scvf];
		}

	///	Convection Length
		number downwind_conv_length(size_t scvf) const
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			return m_vDownConvLength[scvf];
		}

	/// returns the upwind velocity
		MathVector<dim> upwind_vel(const size_t scvf, const LocalVector& CornerVel, const MathVector<dim> vStdVel[]) const;
    /// returns the upwind velocity
        MathVector<dim> upwind_momentum(const size_t scvf, const LocalVector& CornerVel, const number DensitySCV[], const MathVector<dim> vStdVel[]) const;
    /// returns the upwind velocity
        MathVector<dim> upwind_rel_vel( const size_t scvf,
                                        const MathVector<dim> vRelVelCorners[],
                                        const MathVector<dim> vRelVel[]) const;

	///	upwind shape for corner vel
		number upwind_shape_sh(size_t scvf, size_t sh) const
		{
			// UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			// UG_NSUPWIND_ASSERT(sh < m_numSh, "Invalid index");
			return m_vvUpShapeSh[scvf][sh];
		}

	///	upwind shape for corner vel
		number downwind_shape_sh(size_t scvf, size_t sh) const
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(sh < m_numSh, "Invalid index");
			return m_vvDownShapeSh[scvf][sh];
		}

	///	returns if upwind shape w.r.t. ip vel is non-zero
		bool non_zero_shape_ip() const {return m_bNonZeroShapeIp;}

	///	upwind shapes for ip vel
		number upwind_shape_ip(size_t scvf, size_t scvf2) const
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index: " << scvf);
			UG_NSUPWIND_ASSERT(scvf2 < m_numScvf, "Invalid index2: " << scvf2);
			return m_vvUpShapeIp[scvf][scvf2];
		}

	///	upwind shapes for ip vel
		number downwind_shape_ip(size_t scvf, size_t scvf2) const
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(scvf2 < m_numScvf, "Invalid index");
			return m_vvDownShapeIp[scvf][scvf2];
		}

	///	compute values for new geometry and corner velocities
		void update(const FVGeometryBase* geo,
		            const MathVector<dim> vStdVel[])
		{
			update_upwind(geo, vStdVel);
		}

	///	compute values for new geometry and corner velocities
		void update_upwind(const FVGeometryBase* geo,
						   const MathVector<dim> vStdVel[])
		{
			compute(geo, vStdVel, m_vvUpShapeSh, m_vvUpShapeIp, m_vUpConvLength);
		}

	///	compute values for new geometry and corner velocities
		void update_downwind(const FVGeometryBase* geo,
		                     const MathVector<dim> vStdVel[])
		{
			MathVector<dim> vDownIPVel[maxNumSCVF];
			for(size_t ip = 0; ip < m_numScvf; ++ip)
				VecScale(vDownIPVel[ip], vStdVel[ip], -1.0);

			compute(geo, vDownIPVel, m_vvDownShapeSh, m_vvDownShapeIp, m_vDownConvLength);
		}

	//////////////////////////
	// internal handling
	//////////////////////////

	protected:
	///	sets the shape ip flag
		void set_shape_ip_flag(bool flag) {m_bNonZeroShapeIp = flag;}

	///	non-const access to upwind shapes for corner vel
		number& upwind_shape_sh(size_t scvf, size_t sh)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(sh < m_numSh, "Invalid index");
			return m_vvUpShapeSh[scvf][sh];
		}

	///	non-const access to upwind shapes for corner vel
		number& downwind_shape_sh(size_t scvf, size_t sh)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(sh < m_numSh, "Invalid index");
			return m_vvDownShapeSh[scvf][sh];
		}

	///	non-const access to upwind shapes for ip vel
		number& upwind_shape_ip(size_t scvf, size_t scvf2)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(scvf2 < m_numScvf, "Invalid index");
			return m_vvUpShapeIp[scvf][scvf2];
		}

	///	non-const access to upwind shapes for ip vel
		number& downwind_shape_ip(size_t scvf, size_t scvf2)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(scvf2 < m_numScvf, "Invalid index");
			return m_vvDownShapeIp[scvf][scvf2];
		}

	///	non-const access to Convection Length
		number& upwind_conv_length(size_t scvf)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			return m_vUpConvLength[scvf];
		}

	///	non-const access to Convection Length
		number& down_upwind_conv_length(size_t scvf)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			return m_vDownConvLength[scvf];
		}

	///	number of current scvf
		size_t m_numScvf;

	///	number of current shape functions (usually in corners)
		size_t m_numSh;

	///	convection length
		number m_vUpConvLength[maxNumSCVF];
		number m_vDownConvLength[maxNumSCVF];

	///	upwind shapes for corners shape functions
		number m_vvUpShapeSh[maxNumSCVF][maxNumSH];
		number m_vvDownShapeSh[maxNumSCVF][maxNumSH];

	///	flag if ip shapes are non-zero
		bool m_bNonZeroShapeIp;

	///	upwind shapes for ip vels
		number m_vvUpShapeIp[maxNumSCVF][maxNumSCVF];
		number m_vvDownShapeIp[maxNumSCVF][maxNumSCVF];

	///	compute values for new geometry and corner velocities
		void compute(const FVGeometryBase* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF])
		{
			(*(m_vComputeFunc[m_id]))(geo, vIPVel, vUpShapeSh, vUpShapeIp, vConvLength);
		}

	//////////////////////////
	// registering process
	//////////////////////////
	public:
	///	register a update function for a Geometry
		template <typename TFVGeom, typename TAssFunc>
		void register_update_func(TAssFunc func);

	///	register a update function for a Geometry
		template <typename TFVGeom, typename TAssFunc>
		void register_update_ip_vel_func(TAssFunc func);

	///	set the Geometry type to use for next updates
		template <typename TFVGeom>
		void set_geometry_type();

	protected:

	///	type of update function
//		typedef void (this_type::*ComputeFunc)(
		typedef void (*ComputeFunc)(
								const FVGeometryBase* obj,
								const MathVector<dim> vIPVel[maxNumSCVF],
								number vUpShapeSh[maxNumSCVF][maxNumSH],
								number vUpShapeIp[maxNumSCVF][maxNumSCVF],
								number vConvLength[maxNumSCVF]);

	///	Vector holding all update functions
		std::vector<ComputeFunc> m_vComputeFunc;

	///	id of current geometry type
		int m_id;
};

//	register a update function for a Geometry
template <int dim>
template <typename TFVGeom, typename TAssFunc>
void
INavierStokesUpwind<dim>::
register_update_func(TAssFunc func)
{
//	get unique geometry id
	size_t id = GetUniqueFVGeomID<TFVGeom>();

//	make sure that there is enough space
	if((size_t)id >= m_vComputeFunc.size())
		m_vComputeFunc.resize(id+1, NULL);

//	set pointer
	m_vComputeFunc[id] = reinterpret_cast<ComputeFunc>(func);
}


//	set the Geometry type to use for next updates
template <int dim>
template <typename TFVGeom>
void
INavierStokesUpwind<dim>::
set_geometry_type()
{

//	get unique geometry id
	size_t id = GetUniqueFVGeomID<TFVGeom>();

//	check that function exists
	if(id >= m_vComputeFunc.size() || m_vComputeFunc[id] == NULL)
		UG_THROW("No update function registered for Geometry "<<id);

//	set current geometry
	m_id = id;

//	set sizes
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	m_numScvf = geo.num_scvf();
	m_numSh = geo.num_sh();
	UG_NSUPWIND_ASSERT(m_numScvf <= maxNumSCVF, "Invalid index");
	UG_NSUPWIND_ASSERT(m_numSh <= maxNumSH, "Invalid index");

}


///	upwind velocity
template <int dim>
MathVector<dim>
INavierStokesUpwind<dim>::
upwind_vel(const size_t scvf,
           const LocalVector& CornerVel,
           const MathVector<dim> vStdVel[]) const
{
//	reset result
	MathVector<dim> vel; VecSet(vel, 0.0);

//	add corner shapes
	for(size_t sh = 0; sh < num_sh(); ++sh)
		for(int d = 0; d < dim; ++d)
			vel[d] += upwind_shape_sh(scvf, sh) * CornerVel(d, sh);

//	done if only depending on shapes
	if(!non_zero_shape_ip()) return vel;

//	compute ip vel
	for(size_t scvf2 = 0; scvf2 < num_scvf(); ++scvf2)
		VecScaleAppend(vel, upwind_shape_ip(scvf, scvf2), vStdVel[scvf2]);

//	return value
	return vel;
}
///    upwind velocity
template <int dim>
MathVector<dim>
INavierStokesUpwind<dim>::
upwind_momentum(const size_t scvf,
           const LocalVector& CornerVel,
           const number DensitySCV[],
           const MathVector<dim> vStdVel[]) const
{
//    reset result
    MathVector<dim> momentum; VecSet(momentum, 0.0);

//    add corner shapes
    for(size_t sh = 0; sh < num_sh(); ++sh)
        for(int d = 0; d < dim; ++d)
            momentum[d] += upwind_shape_sh(scvf, sh) * CornerVel(d, sh) * DensitySCV[sh];

/*//    done if only depending on shapes
    if(!non_zero_shape_ip()) return vel;

//    compute ip vel
    for(size_t scvf2 = 0; scvf2 < num_scvf(); ++scvf2)
        VecScaleAppend(vel, upwind_shape_ip(scvf, scvf2), vStdVel[scvf2]);*/

//    return value
    return momentum;
}

///    upwind velocity
template <int dim>
MathVector<dim>
INavierStokesUpwind<dim>::
upwind_rel_vel(const size_t scvf,
               const MathVector<dim> vRelVelCorners[],
           const MathVector<dim> vRelVel[]) const
{
//    reset result
    MathVector<dim> vel; VecSet(vel, 0.0);

//    add corner shapes
    for(size_t sh = 0; sh < num_sh(); ++sh)
        for(int d = 0; d < dim; ++d)
            vel[d] += upwind_shape_sh(scvf, sh) * vRelVelCorners[sh][d];

//    done if only depending on shapes
    if(!non_zero_shape_ip()) return vel;

//    compute ip vel
    for(size_t scvf2 = 0; scvf2 < num_scvf(); ++scvf2)
        VecScaleAppend(vel, upwind_shape_ip(scvf, scvf2), vRelVel[scvf2]);

//    return value
    return vel;
}

template <template <class Elem, int WorldDim> class TFVGeom, int dim, typename TImpl>
class NavierStokesUpwindRegister
{
	public:
	///	Base class
		typedef INavierStokesUpwind<dim> base_type;

	protected:
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesUpwindRegister() {register_func(Int2Type<dim>());}

	private:
		void register_func(Int2Type<1>)
		{
			register_func<RegularEdge>();
		}

		void register_func(Int2Type<2>)
		{
			register_func<Triangle>();
			register_func<Quadrilateral>();
		}

		void register_func(Int2Type<3>)
		{
			register_func<Tetrahedron>();
			register_func<Pyramid>();
			register_func<Prism>();
			register_func<Hexahedron>();
		}

		template <typename TElem>
		void register_func()
		{
			typedef TFVGeom<TElem, dim> TGeom;
//			typedef void (TImpl::*TFunc)(
			typedef void (*TFunc)(
					const TGeom* obj,
					const MathVector<dim> vIPVel[maxNumSCVF],
					number vUpShapeSh[maxNumSCVF][maxNumSH],
					number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					number vConvLength[maxNumSCVF]);

			getImpl().template register_update_func<TGeom, TFunc>(&TImpl::template compute<TElem>);
		}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}
};

/// creates upwind based on a string identifier
template <int dim>
SmartPtr<INavierStokesUpwind<dim> > CreateNavierStokesUpwind(const std::string& name);

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__UPWIND_INTERFACE__ */
