/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
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

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__COMPRESSIBLE_NAVIER_STOKES_BASE__
#define __H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__COMPRESSIBLE_NAVIER_STOKES_BASE__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../navier_stokes_base.h"

namespace ug{
namespace NavierStokes{

/// \ingroup lib_disc_elem_disc
/// @{

/// Finite Volume Element Discretization for the compressible Navier-Stokes Equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the compressible Navier-Stokes equation for an thermical and caloric
 * ideal gas.
 *
 * The unknowns of the equation are:
 * <ul>
 * <li> \f$ \vec{u} \f$ velocity
 * <li> \f$ p \f$ pressure
 * <li> \f$ \rho \f$ density.
 * </ul>
 *
 * The equation takes the form
 * \f{align*}
 * 	\frac{\partial \rho \vec{u}}{\partial t}
 * 	- \nabla \left( \rho \nu (\nabla \vec{u} + (\nabla \vec{u})^T)
 * 		- \frac{2}{3} Id * \nabla \vec{u} \right)
 * 	+ \nabla \cdot \left( \rho \vec{u} \vec{u}^T \right)
 *  + \nabla p
 *  &= \vec{f}_m\\
 *
 *  \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) &= f_c \\
 *
 *  \frac{1}{\gamma - 1}\frac{\partial \rho}{\partial t}
 *  + \frac{1}{2} \frac{\partial \rho \vec{u}^2}{\partial t}
 *  + \frac{\gamma}{\gamma - 1} \nabla \cdot (\rho \vec{u}^T)
 *  + \frac{1}{2} \nabla \cdot (\rho \vec{u}^2 \vec{u}^T)
 *  - \nabla \left( \vec{u}^T \cdot \tau \right) &= f_e \\
 *
 * \f}
 *
 * with
 * <ul>
 * <li> \f$ \gamma \f$ is the adiabatic index
 * <li> \f$ \tau \f$ is the shear stress tensor \f$ \tau = \rho \nu (\nabla \vec{u} + (\nabla \vec{u})^T)
 * 		- \frac{2}{3} Id * \nabla \vec{u} \f$
 * <li>	\f$ \nu \f$ is the kinematic viscosity (temporarily constant, implement)
 * <li>	\f$ \vec{f}_m \equiv f_m(\vec{x},t) \f$
 * 		is a Source Term for the momentum equation (not implemented yet)
 * <li> \f$ f_c \f$ is a Source Term for the continuity equation (not implemented yet)
 * <li> \f$ f_e \f$ is a Source Term for the energy preserving equation (not implemented yet)
 * </ul>
 *
 * The first equation models the conservation of momentum and is therefore
 * referred to as the <b>Momentum equation</b>. The second equation models the
 * conservation of mass and is known as the <b>Mass continuity equation</b> or
 * simply <b>Continuity equation</b>. The third equation models the conservation
 * of energy and is therefore referred to as the <b>Energy equation</b>.
 *
 * These three equations are normalized by reference quantities (c.f. A. Gordner
 * 'Numerische Simulation nichtlinearer Aeroakustik bei kleinen Machzahlen'),
 * <ul>
 * <li> \f$ M_{ref}^2 \f$ the squared reference Mach number
 * 		\f$ M_{ref}^2 = \frac{\rho_{ref} \vec{u}_{ref}^2}{p_{ref}} \f$
  * <li> \f$ Re \f$ the Reynolds number
 * 		\f$ Re = \frac{\rho_{ref} L_{ref} \vec{u}_{ref}}{\nu} \f$
 * </ul>
 *
 * with
 * <ul>
 * <li> \f$ \vec{u}_{ref} = \max_{x,t} |\vec{u}(x,t)|\f$
 * <li> \f$ p_{ref} = \max_{x,t} |p(x,t)|\f$
 * <li> \f$ \rho_{ref} = \max_{x,t} |\rho(x,t)|\f$
 * <li> \f$ L_{ref} \f$ is a reference length.
 * </ul>
 *
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */
template<	typename TDomain>
class CompressibleNavierStokesBase
	: public NavierStokesBase<TDomain>
{
	protected:
	///	Base class type
		typedef NavierStokesBase<TDomain> base_type;

	///	own type
		typedef CompressibleNavierStokesBase<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
	/// \{
		CompressibleNavierStokesBase(const char* functions, const char* subsets);
		CompressibleNavierStokesBase(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
	/// \}

	///	sets the kinematic viscosity
	/**
	 * This method sets the kinematic viscosity value.
	 *
	 * \param[in]	data		kinematic Viscosity
	 */
	///	\{
		virtual void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user) = 0;

	///	returns kinematic viscosity
		virtual SmartPtr<CplUserData<number, dim> > kinematic_viscosity() = 0;

	///	sets the adiabatic index (also known as 'heat capacity ratio' or 'ratio of specific heats')
		virtual void set_adiabatic_index(SmartPtr<CplUserData<number, dim> > user) = 0;
		void set_adiabatic_index(number val);
#ifdef UG_FOR_LUA
		void set_adiabatic_index(const char* fctName);
#endif

	///	returns adiabatic index
		virtual SmartPtr<CplUserData<number, dim> > adiabatic_index() = 0;

	///	sets the source function
	/**
	 * This method sets the source value. A zero value is assumed as default.
	 * \param[in]	data		source data
	 */
	///	\{
		virtual void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user) = 0;
    
    ///    returns source
        virtual SmartPtr<CplUserData<MathVector<dim>, dim> > source() = 0;

	///	sets if Mach-number blending is used in momentum equation
		void set_mach_number_blend(bool machNrBlend) {m_bMachNrBlend = machNrBlend;}

  	public:
	///	returns if local time series is needed
		virtual bool requests_local_time_series() {return true;}

	///	returns string identifying disc type
		virtual std::string disc_type() const = 0;

	protected:
	///	flag if using Mach-number Blending
		bool m_bMachNrBlend;

	///	factor for exact jacobian, (1 for exact jacobian, 0 for fix point)
		using base_type::m_bFullNewtonFactor;
};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__COMPRESSIBLE_NAVIER_STOKES_BASE__*/


