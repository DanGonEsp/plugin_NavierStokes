/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__STABILIZATION__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__STABILIZATION__

//#define UG_NSSTAB_ASSERT(cond, exp)
// include define below to assert arrays used in stabilization
#define UG_NSSTAB_ASSERT(cond, exp) UG_ASSERT((cond), (exp))

#include "../../upwind_interface.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"


namespace ug{
namespace NavierStokes{

/**
 * Class for a general stabilization method for FV1-discretizations of the NS equations.
 *
 * A general stabilization is a method that expresses the additional DoFs of
 * the velocity (for the continuity equation) placed in the integration points
 * in terms of the corner values of the velocity and the pressure. Thus, the
 * the stabilization is a special interpolation of the velocity (involving the
 * pressure) from the corners of the primary grid elements into the integration
 * points of the secondary grid.
 *
 * \tparam	dim		dimension of the world
 */
template <int dim>
class INavierStokesFV1Stabilization
{
	protected:
	///	used traits
		typedef fv1_dim_traits<dim, dim> traits;

	public:
	///	number of SubControlVolumes
		static const size_t maxNumSCV = traits::maxNumSCV;

	///	max number of SubControlVolumeFaces
		static const size_t maxNumSCVF = traits::maxNumSCVF;

	/// max number of shape functions
		static const size_t maxNumSH = traits::maxNSH;

	private:
	/// abbreviation for own type
		typedef INavierStokesFV1Stabilization<dim> this_type;

	public:
	///	constructor
		INavierStokesFV1Stabilization()
		:	m_numScvf(0), m_numSh(0),
			m_spUpwind(NULL)
		{
			m_vUpdateFunc.clear();
		}

	///	sets the upwind method
		void set_upwind(SmartPtr<INavierStokesUpwind<dim> > spUpwind)
		{
			m_spUpwind = spUpwind;
			m_spConstUpwind = spUpwind;
		}

	///	returns the upwind
		const ConstSmartPtr<INavierStokesUpwind<dim> >& upwind() const {return m_spConstUpwind;}
		
	/// returns if the upwind pointer is valid
		bool upwind_valid() const {return m_spConstUpwind.valid();}

	///	set the FV1 Geometry type to use for next updates
		template <typename TFVGeom>
		void set_geometry_type()
		{
		//	get unique geometry id
			size_t id = GetUniqueFVGeomID<TFVGeom>();
		//	check that function exists
			if(id >= m_vUpdateFunc.size() || m_vUpdateFunc[id] == NULL)
				UG_THROW("No update function registered for Geometry "<<id);
		//	set current geometry
			m_id = id;

		//	set sizes
			TFVGeom& geo = GeomProvider<TFVGeom>::get();
			m_numScvf = geo.num_scvf();
			m_numSh = geo.num_sh();
		//	set sizes in upwind
			if(m_spUpwind.valid()) m_spUpwind->template set_geometry_type<TFVGeom>();
		}
	
	/////////////////////////////////////////
	// interface to the upwind
	/////////////////////////////////////////
	
	protected:
	
	///	computes the upwind shapes
		template <typename TFVGeom>
		void compute_upwind(const TFVGeom& geo,
		                    const MathVector<dim> vStdVel[])
		{
			UG_NSSTAB_ASSERT(m_spUpwind.valid(), "No upwind object");
			m_spUpwind->update_upwind(geo, vStdVel);
		};

	///	computes the downwind shapes
		template <typename TFVGeom>
		void compute_downwind(const TFVGeom& geo,
		                      const MathVector<dim> vStdVel[])
		{
			UG_NSSTAB_ASSERT(m_spUpwind.valid(), "No upwind object");
			m_spUpwind->update_downwind(geo, vStdVel);
		};
    
		
	/////////////////////////////////////////
	// the data interface (for the NS discretization)
	/////////////////////////////////////////

	public:
	
	///	number of integration points
		size_t num_ip () const {return m_numScvf;}
		
	///	number of shapes (corners)
		size_t num_sh () const {return m_numSh;}
	
	/// stabilized velocity
		const MathVector<dim>& stab_vel(size_t scvf) const
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			return m_vStabVel[scvf];
		}

	///	returns if stab velocity comp depends on other vel components
		bool vel_comp_connected() const {return m_bVelCompConnected;}

	/// computed stab shape for velocity. This is: The stab_vel derivative
	/// w.r.t velocity unknowns in the corner for each component
		number stab_shape_vel(size_t scvf, size_t compOut, size_t compIn, size_t sh) const
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
			UG_NSSTAB_ASSERT(compIn < dim, "Invalid index.");
			UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
			return m_vvvvStabShapeVel[scvf][compOut][compIn][sh];
		}

	///	computed stab shape for pressure.
		number stab_shape_p(size_t scvf, size_t compOut, size_t sh) const
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
			UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
			return m_vvvStabShapePressure[scvf][compOut][sh];
		}
    ///    computed stab shape for pressure.
        number stab_shape_p_jump(size_t scvf, size_t compOut, size_t sh) const
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vvvStabShapePressureJump[scvf][compOut][sh];
        }
        number stab_shape_slip_vel(size_t scvf, size_t compOut, size_t compIn, size_t sh) const
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(compIn < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vvvvStabShapeSlipVel[scvf][compOut][compIn][sh];
        }
    
    // Second Vel
    /// stabilized velocity
        const MathVector<dim>& stab_vel_2(size_t scvf) const
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            return m_vStabVel_2[scvf];
        }

    /// computed stab shape for velocity. This is: The stab_vel derivative
    /// w.r.t velocity unknowns in the corner for each component
        number stab_shape_vel_2(size_t scvf, size_t compOut, size_t compIn, size_t sh) const
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(compIn < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vvvvStabShapeVel_2[scvf][compOut][compIn][sh];
        }

    ///    computed stab shape for pressure.
        number stab_shape_p_2(size_t scvf, size_t compOut, size_t sh) const
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vvvStabShapePressure_2[scvf][compOut][sh];
        }
    ///    computed stab shape for pressure.
        number stab_shape_p_jump_2(size_t scvf, size_t compOut, size_t sh) const
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vvvStabShapePressureJump_2[scvf][compOut][sh];
        }
        number stab_shape_slip_vel_2(size_t scvf, size_t compOut, size_t compIn, size_t sh) const
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(compIn < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vvvvStabShapeSlipVel_2[scvf][compOut][compIn][sh];
        }
    


	///	compute values for new geometry and corner velocities
		void update(const FVGeometryBase* geo,
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
                    const number theta)
			{(this->*(m_vUpdateFunc[m_id]))(geo, vCornerValue, vStdVel,
											bStokes, kinVisco, kinViscoSCV, density, densitySCV, pressure_jump, SlipVel, jump_shape, normal, pSource,
											pvCornerValueOldTime, dt, density_ref, multiphase, phase_2, theta);}

	/////////////////////////////////////////
	// the data interface (for the implementation)
	/////////////////////////////////////////

	protected:
	
	/// stabilized velocity
		MathVector<dim>& stab_vel(size_t scvf)
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			return m_vStabVel[scvf];
		}

	///	sets the vel comp connected flag
		void set_vel_comp_connected(bool bVelCompConnected) {m_bVelCompConnected = bVelCompConnected;}

	/// computed stab shape for velocity. This is: The stab_vel derivative
	/// w.r.t velocity unknowns in the corner for each component
		number& stab_shape_vel(size_t scvf, size_t compOut, size_t compIn, size_t sh)
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
			UG_NSSTAB_ASSERT(compIn < dim, "Invalid index.");
			UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
			return m_vvvvStabShapeVel[scvf][compOut][compIn][sh];
		}

	///	computed stab shape for pressure.
		number& stab_shape_p(size_t scvf, size_t compOut, size_t sh)
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
			UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
			return m_vvvStabShapePressure[scvf][compOut][sh];
		}
    ///    computed stab shape for pressure jump.
        number& stab_shape_p_jump(size_t scvf, size_t compOut, size_t sh)
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vvvStabShapePressureJump[scvf][compOut][sh];
        }
        number& stab_shape_slip_vel(size_t scvf, size_t compOut, size_t compIn, size_t sh)
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(compIn < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vvvvStabShapeSlipVel[scvf][compOut][compIn][sh];
        }
    
    
    /// stabilized velocity For second phase
        MathVector<dim>& stab_vel_2(size_t scvf)
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            return m_vStabVel_2[scvf];
        }


    /// computed stab shape for velocity. This is: The stab_vel derivative
    /// w.r.t velocity unknowns in the corner for each component
        number& stab_shape_vel_2(size_t scvf, size_t compOut, size_t compIn, size_t sh)
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(compIn < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vvvvStabShapeVel_2[scvf][compOut][compIn][sh];
        }

    ///    computed stab shape for pressure.
        number& stab_shape_p_2(size_t scvf, size_t compOut, size_t sh)
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vvvStabShapePressure_2[scvf][compOut][sh];
        }
    ///    computed stab shape for pressure jump.
        number& stab_shape_p_jump_2(size_t scvf, size_t compOut, size_t sh)
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vvvStabShapePressureJump_2[scvf][compOut][sh];
        }
        number& stab_shape_slip_vel_2(size_t scvf, size_t compOut, size_t compIn, size_t sh)
        {
            UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(compIn < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vvvvStabShapeSlipVel_2[scvf][compOut][compIn][sh];
        }

	/////////////////////////////////////////
	// the data
	/////////////////////////////////////////
	
	private:

	///	number of current scvf
		size_t m_numScvf;

	///	number of current shape functions (usually in corners)
		size_t m_numSh;

	///	values of stabilized velocity at ip
		MathVector<dim> m_vStabVel[maxNumSCVF];

	///	flag if velocity components are interconnected
		bool m_bVelCompConnected;

	///	stab shapes w.r.t vel
		number m_vvvvStabShapeVel[maxNumSCVF][dim][dim][maxNumSH];

	///	stab shapes w.r.t pressure
		number m_vvvStabShapePressure[maxNumSCVF][dim][maxNumSH];
    
    ///    stab shapes w.r.t pressure jump
        number m_vvvStabShapePressureJump[maxNumSCVF][dim][maxNumSH];
    
    ///    stab shapes w.r.t Slip vel
        number m_vvvvStabShapeSlipVel[maxNumSCVF][dim][dim][maxNumSH];
    
    //-------------------------------
    
    ///    values of stabilized velocity at ip
        MathVector<dim> m_vStabVel_2[maxNumSCVF];

    ///    stab shapes w.r.t vel
        number m_vvvvStabShapeVel_2[maxNumSCVF][dim][dim][maxNumSH];

    ///    stab shapes w.r.t pressure
        number m_vvvStabShapePressure_2[maxNumSCVF][dim][maxNumSH];
    
    ///    stab shapes w.r.t pressure jump
        number m_vvvStabShapePressureJump_2[maxNumSCVF][dim][maxNumSH];
    
    ///    stab shapes w.r.t Slip vel
        number m_vvvvStabShapeSlipVel_2[maxNumSCVF][dim][dim][maxNumSH];

	///	id of current geometry type
		int m_id;

	///	Upwind object (if set)
		SmartPtr<INavierStokesUpwind<dim> > m_spUpwind;
	///	Upwind object (if set), the same as m_spUpwind, but as a const
		ConstSmartPtr<INavierStokesUpwind<dim> > m_spConstUpwind;
	
	//////////////////////////
	// registering mechanism
	//////////////////////////

	protected:
	///	type of update function
		typedef void (this_type::*UpdateFunc)(	const FVGeometryBase* geo,
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
                                                const number theta);

	public:
	///	register a update function for a Geometry
		template <typename TFVGeom, typename TAssFunc>
		void register_update_func(TAssFunc func);

	protected:
	///	Vector holding all update functions
		std::vector<UpdateFunc> m_vUpdateFunc;
};

/*-------- Schneider-Raw-type stabilizations --------*/

/**
 * The base class for the Schneider-Raw-type stabilizations for FV1-discretizations of the NS equations.
 *
 * \tparam	dim		dimension of the world
 */
template <int dim>
class INavierStokesSRFV1Stabilization
	:	public INavierStokesFV1Stabilization<dim>
{
	private:
	/// Abbreviation for own type
		typedef INavierStokesSRFV1Stabilization<dim> this_type;

	protected:
		enum DiffusionLength
		{
		    RAW = 0,
		    FIVEPOINT,
		    COR
		};

	public:
	///	constructor
		INavierStokesSRFV1Stabilization()
		{
		//	default setup
			set_diffusion_length("RAW");
		}

	///	sets the type of diff length used for evaluation
		void set_diffusion_length(std::string diffLength);

	protected:
	
	///	diff length
		number diff_length_sq_inv(size_t scvf) const
		{
			UG_NSSTAB_ASSERT(scvf < this_type::num_ip(), "Invalid index.");
			return m_vDiffLengthSqInv[scvf];
		}

	/////////////////////////////////////////
	// forward methods of Upwind Velocity
	/////////////////////////////////////////

	///	Convection Length
		number upwind_conv_length(size_t scvf) const
		{
			UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
			return this_type::upwind()->upwind_conv_length(scvf);
		}

	///	Convection Length
		number downwind_conv_length(size_t scvf) const
		{
			UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
			return this_type::upwind()->downwind_conv_length(scvf);
		}

	///	upwind shape for corner vel
		number upwind_shape_sh(size_t scvf, size_t sh) const
		{
			UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
			return this_type::upwind()->upwind_shape_sh(scvf, sh);
		}

	///	upwind shape for corner vel
		number downwind_shape_sh(size_t scvf, size_t sh) const
		{
			UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
			return this_type::upwind()->downwind_shape_sh(scvf, sh);
		}

	///	returns if upwind shape w.r.t. ip vel is non-zero
		bool non_zero_shape_ip() const
		{
			UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
			return this_type::upwind()->non_zero_shape_ip();
		}

	///	upwind shapes for ip vel
		number upwind_shape_ip(size_t scvf, size_t scvf2) const
		{
			UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
			return this_type::upwind()->upwind_shape_ip(scvf, scvf2);
		}

	///	upwind shapes for ip vel
		number downwind_shape_ip(size_t scvf, size_t scvf2) const
		{
			UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
			return this_type::upwind()->downwind_shape_ip(scvf, scvf2);
		}

	//////////////////////////
	// internal handling
	//////////////////////////

	protected:
	///	computes the diffusion length
		template <typename TFVGeom>
		void compute_diff_length(const TFVGeom& geo);

	//////////////////////////
	// data
	//////////////////////////

	protected:
	///	type of diffusion length computation
		DiffusionLength m_diffLengthType;

	///	vector holding diffusion Length squared and inverted
		number m_vDiffLengthSqInv[this_type::maxNumSCVF];

};

/// creates upwind based on a string identifier
template <int dim>
SmartPtr<INavierStokesSRFV1Stabilization<dim> > CreateNavierStokesStabilization(const std::string& name);

/////////////////////////////////////////////////////////////////////////////
// FIELDS
/////////////////////////////////////////////////////////////////////////////

/**
 * Implementation of the FIELDS stabilization
 */
template <int TDim>
class NavierStokesFIELDSStabilization
	: public INavierStokesSRFV1Stabilization<TDim>
{
	public:
	///	Base class
		typedef INavierStokesSRFV1Stabilization<TDim> base_type;

	///	This class
		typedef NavierStokesFIELDSStabilization<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::register_update_func;
		using base_type::diff_length_sq_inv;
		using base_type::stab_shape_vel;
		using base_type::stab_shape_p;
		using base_type::stab_vel;
		using base_type::set_vel_comp_connected;

	//	functions from upwind
		using base_type::upwind_conv_length;
		using base_type::downwind_conv_length;
		using base_type::upwind_shape_sh;
		using base_type::downwind_shape_sh;
		using base_type::non_zero_shape_ip;
		using base_type::upwind_shape_ip;
		using base_type::downwind_shape_ip;

	public:
	///	constructor
		NavierStokesFIELDSStabilization()
		{
		//	vel comp not coupled
			set_vel_comp_connected(false);

		//	register evaluation function
			register_func();
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		void update(const FV1Geometry<TElem, dim>* geo,
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
                    const number theta);

	private:
		void register_func();

		template <typename TElem>
		void register_func()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef void (this_type::*TFunc)(const TGeom* geo,
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
                                             const number theta);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};


/////////////////////////////////////////////////////////////////////////////
// FLOW
/////////////////////////////////////////////////////////////////////////////

/**
 * Implementation of the FLOW stabilization
 */
template <int TDim>
class NavierStokesFLOWStabilization
	: public INavierStokesSRFV1Stabilization<TDim>
{
	public:
	///	Base class
		typedef INavierStokesSRFV1Stabilization<TDim> base_type;

	///	This class
		typedef NavierStokesFLOWStabilization<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::register_update_func;
		using base_type::set_vel_comp_connected;
		using base_type::diff_length_sq_inv;
		using base_type::stab_shape_vel;
		using base_type::stab_shape_p;
		using base_type::stab_vel;

	//	functions from upwind
		using base_type::upwind_conv_length;
		using base_type::downwind_conv_length;
		using base_type::upwind_shape_sh;
		using base_type::downwind_shape_sh;
		using base_type::non_zero_shape_ip;
		using base_type::upwind_shape_ip;
		using base_type::downwind_shape_ip;

	public:
	///	constructor
		NavierStokesFLOWStabilization()
		{
		//	vel comp coupled
			set_vel_comp_connected(true);

		//	register evaluation function
			register_func();
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		void update(const FV1Geometry<TElem, dim>* geo,
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
                    const number theta);

	private:
		void register_func();

		template <typename TElem>
		void register_func()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef void (this_type::*TFunc)(const TGeom* geo,
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
                                             const number theta);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};

/////////////////////////////////////////////////////////////////////////////
// FIELDS
/////////////////////////////////////////////////////////////////////////////

/**
 * Implementation of the FIELDS stabilization
 */
template <int TDim>
class NavierStokesFIELDS_2_Stabilization
    : public INavierStokesSRFV1Stabilization<TDim>
{
    public:
    ///    Base class
        typedef INavierStokesSRFV1Stabilization<TDim> base_type;

    ///    This class
        typedef NavierStokesFIELDS_2_Stabilization<TDim> this_type;

    ///    Dimension
        static const int dim = TDim;

    protected:
    //    explicitly forward some function
        using base_type::register_update_func;
        using base_type::diff_length_sq_inv;
        using base_type::stab_shape_vel;
        using base_type::stab_shape_p;
        using base_type::stab_vel;
        using base_type::set_vel_comp_connected;
    
        using base_type::stab_shape_p_jump;
        using base_type::stab_shape_slip_vel;
    
    
        using base_type::stab_shape_vel_2;
        using base_type::stab_shape_p_2;
        using base_type::stab_vel_2;

        using base_type::stab_shape_p_jump_2;
        using base_type::stab_shape_slip_vel_2;

    //    functions from upwind
        using base_type::upwind_conv_length;
        using base_type::downwind_conv_length;
        using base_type::upwind_shape_sh;
        using base_type::downwind_shape_sh;
        using base_type::non_zero_shape_ip;
        using base_type::upwind_shape_ip;
        using base_type::downwind_shape_ip;

    public:
    ///    constructor
        NavierStokesFIELDS_2_Stabilization()
        {
        //    vel comp not coupled
            set_vel_comp_connected(true);

        //    register evaluation function
            register_func();
        }

    ///    update of values for FV1Geometry
        template <typename TElem>
        void update(const FV1Geometry<TElem, dim>* geo,
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
                    const number theta);
    

    private:
        void register_func();

        template <typename TElem>
        void register_func()
        {
            typedef FV1Geometry<TElem, dim> TGeom;
            typedef void (this_type::*TFunc)(const TGeom* geo,
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
                                             const number theta);

            this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
        }
};

/////////////////////////////////////////////////////////////////////////////
// FIELDS  3
/////////////////////////////////////////////////////////////////////////////

/**
 * Implementation of the FIELDS stabilization
 */
template <int TDim>
class NavierStokesVISCOSITY_Stabilization
    : public INavierStokesSRFV1Stabilization<TDim>
{
    public:
    ///    Base class
        typedef INavierStokesSRFV1Stabilization<TDim> base_type;

    ///    This class
        typedef NavierStokesVISCOSITY_Stabilization<TDim> this_type;

    ///    Dimension
        static const int dim = TDim;

    protected:
    //    explicitly forward some function
        using base_type::register_update_func;
        using base_type::diff_length_sq_inv;
        using base_type::stab_shape_vel;
        using base_type::stab_shape_p;
        using base_type::stab_vel;
        using base_type::set_vel_comp_connected;
    
        using base_type::stab_shape_p_jump;
        using base_type::stab_shape_slip_vel;

    //    functions from upwind
        using base_type::upwind_conv_length;
        using base_type::downwind_conv_length;
        using base_type::upwind_shape_sh;
        using base_type::downwind_shape_sh;
        using base_type::non_zero_shape_ip;
        using base_type::upwind_shape_ip;
        using base_type::downwind_shape_ip;

    public:
    ///    constructor
        NavierStokesVISCOSITY_Stabilization()
        {
        //    vel comp not coupled
            set_vel_comp_connected(true);

        //    register evaluation function
            register_func();
        }

    ///    update of values for FV1Geometry
        template <typename TElem>
        void update(const FV1Geometry<TElem, dim>* geo,
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
                    const number theta);

    private:
        void register_func();

        template <typename TElem>
        void register_func()
        {
            typedef FV1Geometry<TElem, dim> TGeom;
            typedef void (this_type::*TFunc)(const TGeom* geo,
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
                                             const number theta);

            this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
        }
};



/////////////////////////////////////////////////////////////////////////////
// FLOW
/////////////////////////////////////////////////////////////////////////////

/**
 * Implementation of the KARIMIAN stabilization
 */
template <int TDim>
class NavierStokesKARIMIANStabilization
    : public INavierStokesSRFV1Stabilization<TDim>
{
    public:
    ///    Base class
        typedef INavierStokesSRFV1Stabilization<TDim> base_type;

    ///    This class
        typedef NavierStokesKARIMIANStabilization<TDim> this_type;

    ///    Dimension
        static const int dim = TDim;

    protected:
    //    explicitly forward some function
        using base_type::register_update_func;
        using base_type::set_vel_comp_connected;
        using base_type::diff_length_sq_inv;
        using base_type::stab_shape_vel;
        using base_type::stab_shape_p;
        using base_type::stab_vel;

    //    functions from upwind
        using base_type::upwind_conv_length;
        using base_type::downwind_conv_length;
        using base_type::upwind_shape_sh;
        using base_type::downwind_shape_sh;
        using base_type::non_zero_shape_ip;
        using base_type::upwind_shape_ip;
        using base_type::downwind_shape_ip;

    public:
    ///    constructor
        NavierStokesKARIMIANStabilization()
        {
        //    vel comp coupled
            set_vel_comp_connected(true);

        //    register evaluation function
            register_func();
        }

    ///    update of values for FV1Geometry
        template <typename TElem>
        void update(const FV1Geometry<TElem, dim>* geo,
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
                    const number theta);

    private:
        void register_func();

        template <typename TElem>
        void register_func()
        {
            typedef FV1Geometry<TElem, dim> TGeom;
            typedef void (this_type::*TFunc)(const TGeom* geo,
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
                                             const number theta);

            this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
        }
};
/////////////////////////////////////////////////////////////////////////////
// NO
/////////////////////////////////////////////////////////////////////////////

/**
 * Implementation of the KARIMIAN stabilization
 */
template <int TDim>
class NavierStokesNOStabilization
    : public INavierStokesSRFV1Stabilization<TDim>
{
    public:
    ///    Base class
        typedef INavierStokesSRFV1Stabilization<TDim> base_type;

    ///    This class
        typedef NavierStokesNOStabilization<TDim> this_type;

    ///    Dimension
        static const int dim = TDim;

    protected:
    //    explicitly forward some function
        using base_type::register_update_func;
        using base_type::set_vel_comp_connected;
        using base_type::diff_length_sq_inv;
        using base_type::stab_shape_vel;
        using base_type::stab_shape_p;
        using base_type::stab_vel;

    //    functions from upwind
        using base_type::upwind_conv_length;
        using base_type::downwind_conv_length;
        using base_type::upwind_shape_sh;
        using base_type::downwind_shape_sh;
        using base_type::non_zero_shape_ip;
        using base_type::upwind_shape_ip;
        using base_type::downwind_shape_ip;

    public:
    ///    constructor
        NavierStokesNOStabilization()
        {
        //    vel comp coupled
            set_vel_comp_connected(false);

        //    register evaluation function
            register_func();
        }

    ///    update of values for FV1Geometry
        template <typename TElem>
        void update(const FV1Geometry<TElem, dim>* geo,
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
                    const number theta);

    private:
        void register_func();

        template <typename TElem>
        void register_func()
        {
            typedef FV1Geometry<TElem, dim> TGeom;
            typedef void (this_type::*TFunc)(const TGeom* geo,
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
                                             const number theta);

            this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
        }
};

/////////////////////////////////////////////////////////////////////////////
// NO STABILIZATION (Note: The discretization is then unstable!)
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesFV1WithoutStabilization
	: public INavierStokesFV1Stabilization<TDim>
{
	public:
	///	Base class
		typedef INavierStokesFV1Stabilization<TDim> base_type;

	///	This class
		typedef NavierStokesFV1WithoutStabilization<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::register_update_func;
		using base_type::set_vel_comp_connected;
		using base_type::stab_shape_vel;
		using base_type::stab_shape_p;
		using base_type::stab_vel;
		
	public:
	///	Constructor
		NavierStokesFV1WithoutStabilization ()
		{
		//	vel comp not interconnected
			set_vel_comp_connected(false);

		//	register evaluation function
			register_func();
		}
	
	///	update of values for FV1Geometry
		template <typename TElem>
		void update(const FV1Geometry<TElem, dim>* geo,
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
                    const number theta);

	private:
		void register_func();

		template <typename TElem>
		void register_func()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef void (this_type::*TFunc)(const TGeom* geo,
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
                                             const number theta);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__FV__STABILIZATION__ */
