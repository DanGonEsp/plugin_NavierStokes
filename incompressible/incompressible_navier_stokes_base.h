/*
 * Copyright (c) 2010-2017:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__INCOMPRESSIBLE_NAVIER_STOKES_BASE__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__INCOMPRESSIBLE_NAVIER_STOKES_BASE__

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

/// Finite Volume Element Discretization for the incompressible Navier-Stokes Equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the incompressible Navier-Stokes equation.
 *
 * The unknowns of the equation are:
 * <ul>
 * <li> \f$ \vec{u} \f$ velocity
 * <li> \f$ p \f$ pressure.
 * </ul>
 *
 * The equation takes the form
 * \f{align*}
 * 	\frac{\partial \rho \vec{u}}{\partial t}
 * 	- \nabla \left( \rho \nu (\nabla \vec{u} + (\nabla \vec{u})^T) \right)
 * 	+ \nabla \cdot \left( \rho \vec{u} \vec{u}^T \right)
 *  + \nabla p
 *  &= \vec{f}\\
 *
 *  \nabla \cdot (\rho \vec{u}) &= 0
 * \f}
 *
 * with
 * <ul>
 * <li>	\f$ \rho \f$ is the constant density
 * <li>	\f$ \nu \f$ is the kinematic viscosity (temporarily constant, implement)
 * <li>	\f$ \vec{f} \equiv f(\vec{x},t) \f$ is a Source Term (not implemented yet)
 * </ul>
 *
 * The first equation models the conservation of momentum and is therefore
 * referred to as the <b>Momentum equation</b>. The second equation models the
 * conservation of mass and is known as the <b>Mass continuity equation</b> or
 * simply <b>Continuity equation</b>.
 *
 * In component-wise notation, the equation reads
 *
 * \f{align*}
 * \frac{\partial \rho u_{d_1}}{\partial t}
 * - \sum_{d_2 = 1}^{\text{dim}} \frac{\partial}{\partial x_{d_2}}
 * 		\left( \rho \nu \left( \frac{\partial u_{d_1}}{\partial x_{d_2}}
 * 			+ \frac{\partial u_{d_2}}{\partial x_{d_1}} \right)\right)
 * + \sum_{d_2 = 1}^{\text{dim}}  \frac{\partial}{\partial x_{d_2}}
 * 		\left( \rho u_{d_2} u_{d_1} \right)
 * + \frac{\partial p}{\partial x_{d_1}}
 * &= f_{d_1} \qquad \forall \quad d_1 = 1, \dots, \text{dim}\\
 *
 * \sum_{d = 1}^{\text{dim}}  \frac{\partial \rho u_d}{\partial x_{d}} &= 0
 * \f}
 *
 * Since the density is assumed to be constant, it set to \f$\rho \equiv 1 \f$.
 *
 * The finite volume element discretization uses a dual mesh consisting of
 * Control Volumes \f$\mathcal{B}_h\f$, that cover the domain. For each
 * control volume \f$B \in \mathcal{B}_h\f$ the following equation is solved
 *
 * \f{align*}
 * 	\frac{\partial}{\partial t} \int_{B} \begin{pmatrix} \vec{u} \\ 0 \end{pmatrix} dV
 *  + \int_{\partial B}
 *  	\begin{pmatrix}
 *  		- \rho \nu \left(\nabla \vec{u} + (\nabla \vec{u})^T \right) \vec{n}
 *  		+ \vec{u} \vec{u}^T \vec{n} + p \vec{n} \\
 * 		 	\vec{u} \vec{n}
 * 		\end{pmatrix} dS
 * = \int_B \begin{pmatrix} \vec{f} \\ 0 \end{pmatrix} dV
 * \f}
 *
 * where \f$ \vec{n}\f$ is the outer normal on the boundary of the control volume
 * and \f$\int_B \cdot \; dV \f$ and \f$\int_{\partial B} \cdot \; dS \f$ are the
 * integration over the control volume and the integration over the boundary of
 * the control volume. By \$f B_{sh} \$f is denoted the control volume associated
 * to a shape function (i.e. vertices, since we use P1 Lagrange ansatz functions),
 * T usually denotes the Element.
 *
 * In order to number the local unknowns, we use the following notation.
 * <ul>
 * <li> \f$ sh = 0, \dots, n_{sh} - 1\f$ loop the \f$n_{sh}\f$ shape functions. For
 * 		this implementation these are the corners of the element
 * <li> \f$ d = 0, \dots, \text{dim} - 1 \f$ are the velocity component of each
 * 		corner
 * <li> _P_ := dim indicates the pressure component
 * </ul>
 * The access to local unknowns of local vectors and matrices is now given by
 * \f{align*}
 * 	&\mathbf{d}(d_1, sh), \mathbf{d}(\text{\_P\_}, sh),\\
 *  &\mathcal{J}(d_1, sh_1, d_2, sh_2), \mathcal{J}(d_1, sh_1, \text{\_P\_}, sh_2),
 * \f}
 *
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */
template<	typename TDomain>
class IncompressibleNavierStokesBase
	: public NavierStokesBase<TDomain>
{
	protected:
	///	Base class type
		typedef NavierStokesBase<TDomain> base_type;

	///	own type
		typedef IncompressibleNavierStokesBase<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
	/// \{
		IncompressibleNavierStokesBase(const char* functions, const char* subsets);
		IncompressibleNavierStokesBase(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
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

	///	sets the density
	/**
	 * This method sets the density value.
	 *
	 * \param[in]	data		density
	 */
	///	\{
		virtual void set_density(SmartPtr<CplUserData<number, dim> > user) = 0;
		void set_density(number val);
#ifdef UG_FOR_LUA
		void set_density(const char* fctName);
#endif
	///	\}

	///	returns density
		virtual SmartPtr<CplUserData<number, dim> > density() = 0;

	///	sets the source function
	/**
	 * This method sets the source value. A zero value is assumed as default.
	 * \param[in]	data		source data
	 */
	///	\{
		virtual void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user) = 0;
    ///    returns source
        virtual SmartPtr<CplUserData<MathVector<dim>, dim> > source() = 0;

	/// switches the convective terms off (to solve the Stokes equation)
	/**
	 * \param[in]	Stokes		true to solve Stokes (i.e. without the convective terms)
	 */
		void set_stokes(bool Stokes) {m_bStokes = Stokes;}
		bool stokes() {return m_bStokes;}
	
	
		void set_div_correction(bool bDiv)
		{
			m_div_correction = bDiv;
		}
	///    returns div_correction
		bool div_correction() {return m_div_correction;}

	///	sets assembling of diffusive term to laplace
	/**
	 * Flag to indicate, that in the diffusive term only the laplacian should
	 * be computed. This is valid only in cases, where the viscosity is
	 * constant, since then it holds that
	 *
	 *  div ( \nu (grad) v + (grad v)^T )
	 *  	=  \nu laplace v + \nu grad( div v )
	 *  	=  \nu laplace v
	 *
	 * for incompressible flow (i.e. div v = 0).
	 */
		void set_laplace(bool bLaplace) {m_bLaplace = bLaplace;}
		bool laplace() {return m_bLaplace;}

    ///	sets if peclet blending is used in momentum equation
        void set_peclet_blend(bool pecletBlend) {m_bPecletBlend = pecletBlend;}

        void set_grad_div(number factor){ m_gradDivFactor = factor; }

  	public:
	///	returns if local time series is needed
		virtual bool requests_local_time_series() {return true;}

	///	returns string identifying disc type
		virtual std::string disc_type() const = 0;

	///	returns the export of the velocity
		SmartPtr<CplUserData<MathVector<dim>, dim> > velocity() {return m_exVelocity;}

	///	returns the export of the velocity gradient
		SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > velocity_grad() {return m_exVelocityGrad;}
    
    ///    returns the export of the velocity at ip
        SmartPtr<CplUserData<MathVector<dim>, dim> > velocity_ip() {return m_exVelocity_ip;}
    
    ///    returns the export of the Pressure
        SmartPtr<CplUserData<number, dim> > pressure() {return m_exPressure;}
    
    ///    returns the export of the pressure gradient
        SmartPtr<CplUserData<MathVector<dim>, dim> > pressure_grad() {return m_exPressureGrad;}


	protected:
	///	flag if using Peclet Blending
		bool m_bPecletBlend;

	///	factor for exact jacobian, (1 for exact jacobian, 0 for fix point)
		using base_type::m_bFullNewtonFactor;
        using base_type:: m_boolFullNewton;

	/// factor for div grad stabilization
		number m_gradDivFactor;
	
	/// factor for div correction for trasport equation
		bool m_div_correction;

	/// flag if solving the Stokes equation
		bool m_bStokes;

	///	flag if using only laplace term
		bool m_bLaplace;

	///	Export for the velocity
		SmartPtr<DataExport<MathVector<dim>,dim> > m_exVelocity;

	///	Export for the velocity gradient
		SmartPtr<DataExport<MathMatrix<dim, dim>,dim> > m_exVelocityGrad;
    
    ///    Export for the velocity at ip
        SmartPtr<DataExport<MathVector<dim>,dim> > m_exVelocity_ip;
    
    ///    Export for the pressure
        SmartPtr<DataExport<number,dim> > m_exPressure;
    
    ///    Export for the pressure gradient
        SmartPtr<DataExport<MathVector<dim>,dim> > m_exPressureGrad;


};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__INCOMPRESSIBLE_NAVIER_STOKES_BASE__*/

