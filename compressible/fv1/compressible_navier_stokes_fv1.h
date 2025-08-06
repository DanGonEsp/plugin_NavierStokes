/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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

/*
 *     main parts are copied from the discretization of the incompressible Navier-Stokes Equations
 *     of Andreas Vogel and Christian Wehner
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__COMPRESSIBLE_NAVIER_STOKES_FV1_H__
#define __H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__COMPRESSIBLE_NAVIER_STOKES_FV1_H__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../compressible_navier_stokes_base.h"
#include "../../upwind_interface.h"

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
class CompressibleNavierStokesFV1
	: public CompressibleNavierStokesBase<TDomain>
{
	protected:
	///	Base class type
		typedef CompressibleNavierStokesBase<TDomain> base_type;

	///	own type
		typedef CompressibleNavierStokesFV1<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
	/// \{
		CompressibleNavierStokesFV1(const char* functions, const char* subsets);
		CompressibleNavierStokesFV1(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
	/// \}

	///	sets the kinematic viscosity
		void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user);

	///	returns kinematic viscosity
		SmartPtr<CplUserData<number, dim> > kinematic_viscosity() {return m_imKinViscosity.user_data ();}

	///	sets the adiabatic index
		void set_adiabatic_index(SmartPtr<CplUserData<number, dim> > user);

	///	returns adiabatic index
		SmartPtr<CplUserData<number, dim> > adiabatic_index() {return m_imAdiabaticIndex.user_data ();}

	///	sets the source function
		void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user);

	/// returns scvf source
		DataImport<MathVector<dim>, dim> sourceSCVF(){ return m_imSourceSCVF;}
    
    ///    returns source
        SmartPtr<CplUserData<MathVector<dim>, dim> > source() {return m_imSourceSCVF.user_data ();}
    

	///	sets an upwinding for the convective term of momentum equation
		void set_upwind(SmartPtr<INavierStokesUpwind<dim> > spUpwind)
			{m_spConvUpwind = spUpwind;}

	///	sets the upwind based on a string identifier
		void set_upwind(const std::string& name)
			{m_spConvUpwind = CreateNavierStokesUpwind<dim>(name);}

	public:
	///	type of trial space for each function used
		void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	///	returns string identifying disc type
		virtual std::string disc_type() const {return "fv1";};

	public:
	///	prepares the element loop
		/**
		 * This function is used to set all reference quantities (like e.g. the reynolds number)
		 * which are necessary the formulate the dimensionless Navier-Stokes-equations
		 */
		template<typename TElem, typename TFVGeom>
		void prep_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	prepares the element loop
	/**
	 * This function is used to prepare the element loop. It gets called, before
	 * the loop over all elements starts. Here, general checks are performed:
	 * - Is the correct finite volume geometry used
	 * - Has the stabilization been specified
	 *
	 * In addition, local coordinates of evaluation points are requested by
	 * the DataImports in case of element-fixed points.
	 */
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for evaluation
	/**
	 * This function prepare a specific element for assembling. This function
	 * is always called before any of the assemble_... functions is performed.
	 * Here, the Finite Volume Geometry is prepared and global positions of
	 * the evaluation points of the DataImports are requested (and local
	 * positions in case of hanging node assembling)
	 *
	 * \param[in]	elem		element to prepare the data for
	 * \param[in]	u			current local solution vector
	 * \param[in]	glob_ind	global indices of the local vector components
	 */
		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

	///	finishes the element loop
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();

	///	adds the stiffness part to the local jacobian
	/**
	 * This function adds the local contributions of the stiffness part to
	 * the local jacobian.
	 *
	 * For the definition of \f$ \vec{f}^{\text{diff}}|_{ip}\f$ and
	 * \f$ \vec{f}^{\text{conv}}|_{ip}\f$ see add_def_A_elem.
	 *
	 * The derivative of the diffusive flux is given by
	 * \f{align*}
	 * 	\frac{\partial \vec{f}^{\text{diff}}_{d_1}|_{ip}}{\partial \vec{u}_{d_2}^{sh}}
	 * 	&= - \nu \left( \delta_{d_1, d_2} \sum_{k=1}^{\text{dim}}
	 * 		\left. \frac{\partial \phi^{sh}}{\partial x_k}\right|_{ip} \vec{n_k}|_{ip}
	 * 		+ \left. \frac{\partial \phi^{sh}}{\partial x_{d_1}}\right|_{ip} \vec{n_{d_2}}|_{ip}
	 * 		\right)\\
	 * 	\frac{\partial \vec{f}^{\text{diff}}_{d_1}|_{ip}}{\partial p^{sh}}
	 * &= 0
	 * \f}
	 *
	 * For the derivative of the convective term, we can use a fixpoint linearization
	 * if we only differentiate the first factor of
	 * \f$
	 * 	\vec{f}^{\text{conv}}_{d_1}
	 * 	= (\vec{u}^{\text{blend}} (\vec{u}^{\text{blend}})^T)_{d_1 k} \vec{n}_k
	 *  = \vec{u}^{\text{blend}}_{d_1} \sum_{k=1}^{\text{dim}} \vec{u}^{\text{blend}}_k \vec{n}_k
	 *  = \vec{u}^{\text{blend}}_{d_1} \langle \vec{u}^{\text{blend}}, \vec{n} \rangle
	 * \f$
	 *
	 * This gives
	 * \f{align*}
	 * 	\frac{\partial \vec{f}^{\text{conv}}_{d_1}}{\partial \vec{u}_{d_2}^{sh}}
	 *  &= \langle \vec{u}^{\text{blend}}, \vec{n} \rangle
	 *  	\frac{\partial \vec{u}^{\text{blend}}_{d_1}}{\partial \vec{u}_{d_2}^{sh}}
	 *  &= \langle \vec{u}^{\text{blend}}, \vec{n} \rangle \left(
	 *  	(1-\omega) \delta_{d_1 d_2} \phi^{sh}
	 *  	+ \omega \frac{\partial \vec{u}^{\text{up}}_{d_1}}{\partial \vec{u}_{d_2}^{sh}}
	 *  \right)
	 * \f}
	 *
	 * The derivative of the pressure term is given by
	 * \f{align*}
	 * 	\frac{\partial p \vec{n}|_{ip}}{\partial \vec{u}_{d_2}^{sh}}
	 * 	&= 0\\
	 * 	\frac{\partial p \vec{n}|_{ip}}{\partial p^{sh}}
	 * &= \phi^{sh}|_{ip} \vec{n}
	 * \f}
	 *
	 * The derivative of the continuity term is given by
	 * \f{align*}
	 * 	\frac{\partial \vec{u}^{\text{stab}}|_{ip} \vec{n}|_{ip}}{\partial \vec{u}_{d_1}^{sh}}
	 * 	&= \sum_{d_2=1}^{dim} \left. \frac{\partial u_{d_2}^{\text{stab}}}
	 * 	{\partial \vec{u}_{d_1}^{sh}} \right|_{ip} \vec{n}_{d_2}\\
	 * 	\frac{\partial \vec{u}^{\text{stab}}|_{ip} \vec{n}|_{ip}}{\partial p^{sh}}
	 * &=\sum_{d_2=1}^{dim} \left. \frac{\partial u_{d_2}^{\text{stab}}}
	 * 	{\partial p^{sh}} \right|_{ip} \vec{n}_{d_2}\\
	 * \f}
	 *
	 * \param[in,out]	J		local jacobian with values added on output
	 * \param[in]		u		local vector of current solution
	 * \param[in]		time	current time step
	 *
	 * \tparam	TElem 	Element type
	 * \tparam	TFVGeom	Finite Volume Geometry
	 */
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	adds the stiffness part to the local defect
	/**
	 * This function adds the contribution of one element of the stiffness part to
	 * the local defect.
	 *
	 * The local solution \f$ \mathbf{u} \f$ passes the p1 conform lagrange
	 * dofs of the element, i.e. the unknowns in the corners of the element.
	 *
	 * Define the functional matrix of the velocity at some integration point by
	 * \f{align*}
	 * 	\nabla \vec{u}|_{ip} := \left(\left.\frac{\partial u_i}{\partial x_j}
	 * 							\right|_{ip}\right)_{i,j = 1,\dots,dim}
	 * \f}
	 *
	 * and the diffusive flux to
	 * \f{align*}
	 * 	\vec{f}^{\text{diff}}|_{ip} := - \nu \left( \nabla \vec{u}|_{ip} +
	 * 			(\nabla \vec{u}|_{ip})^T) \right) \cdot \vec{n}_{ip}
	 * \f}
	 *
	 * Suppose that a procedure producing some upwind velocity \f$ \vec{u}^{up}\f$
	 * is available for each ip. Then, we define the convective flux to
	 * \f{align*}
	 * 	\vec{f}^{\text{conv}}|_{ip} := \left( \vec{u}^{up}|_{ip}  \vec{u}^{up}|_{ip}^T
	 * 									\right) \cdot \vec{n}_{ip}
	 * \f}
	 *
	 * Define the pressure at the ip to be interpolated as
	 * \f{align*}
	 * 	p|_{ip} := \sum_{sh=0}^{n_{sh}-1} \mathbf{u}(\text{\_P\_}, sh) \phi_{sh}
	 * \f}
	 *
	 * The contribution added for the <b>momentum equation</b> is:
	 * \f{align*}
	 *  \forall d=0,\dots, \text{dim}-1, \quad
	 *  \forall sh=0,\dots, n_{sh}-1: \qquad
	 *	\mathbf{d}(d, sh) &\text{ +=}
	 *		\int_{\partial B_{sh} \cap T} \vec{f}^{\text{diff}}_d
	 *			+ \vec{f}^{\text{conv}}_d
	 *			+ p \vec{n}_d \;dS
	 *
	 *		\approx \sum_{i=1}^{n_{\partial B_{sh}}}
	 *			 	|\partial B^i_{sh} \cap T|
	 *			 	\left( \left.\vec{f}^{\text{diff}}_d\right|_{ip}
	 *			 	+ \left.\vec{f}^{\text{conv}}_d\right|_{ip}
	 *			 	+ p|_{ip} \vec{n}|_{ip} \right)\\
	 * \f}
	 *
	 * Assume, that some stabilized velocity \f$ \vec{u}^{stab}\f$ is available
	 * for each ip.
	 * Then, the contribution added for the <b>continuity equation</b> is:
	 * \f{align*}
	 *  \forall sh=0,\dots, n_{sh}-1: \qquad
	 *	\mathbf{d}(\text{\_P\_}, sh) &\text{ +=}
	 *		\int_{\partial B_{sh} \cap T} \vec{u}^{stab} \cdot \vec{n} \;dS
	 *
	 *		\approx \sum_{i=1}^{n_{\partial B_{sh}}}
	 *			 	|\partial B^i_{sh} \cap T| \; \vec{u}^{stab}|_{ip} \cdot \vec{n}|_{ip}\\
	 * \f}
	 *
	 * \param[in,out]	d		local defect with values added on output
	 * \param[in]		u		local vector of current solution
	 * \param[in]		time	current time step
	 *
	 * \tparam	TElem 	Element type
	 * \tparam	TFVGeom	Finite Volume Geometry
	 */
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	adds the mass part to the local jacobian
	/**
	 * This function adds the contribution of one element of the mass part to
	 * the local jacobian.
	 *
	 * The local solution \f$ \mathbf{u} \f$ passes the p1 conform lagrange
	 * dofs of the element, i.e. the unknowns in the corners of the element.
	 *
	 * The contribution added is:
	 * \f{align*}
	 *  \forall d=0,\dots, \text{dim}-1, \quad
	 *  \forall sh=0,\dots, n_{sh}-1: \qquad
	 *	\mathcal{J}(d, sh, d, sh) &\text{ +=} \int_{B_{sh} \cap T} dV
	 *							  = |B_{sh} \cap T|\\
	 * \f}
	 *
	 * \param[in,out]	J		local jacobian with values added on output
	 * \param[in]		u		local vector of current solution
	 * \param[in]		time	current time step
	 *
	 * \tparam	TElem 	Element type
	 * \tparam	TFVGeom Finite Volume Geometry
	 */
		template <typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	adds the mass part to the local defect
	/**
	 * This function adds the contribution of one element of the mass part to
	 * the local defect.
	 *
	 * The local solution \f$ \mathbf{u} \f$ passes the p1 conform lagrange
	 * dofs of the element, i.e. the unknowns in the corners of the element.
	 *
	 * The contribution added is:
	 * \f{align*}
	 *  \forall d=0,\dots, \text{dim}-1, \quad
	 *  \forall sh=0,\dots, n_{sh}-1: \qquad
	 *	\mathbf{d}(d, sh) &\text{ +=} \int_{B_{sh} \cap T} \vec{u}_d dV
	 *					  \approx |B_{sh} \cap T| \mathbf{u}(d, sh)\\
	 *
	 *	\mathbf{d}(\text{\_P\_, sh}) &\text{ +=} 0
	 * \f}
	 *
	 * \param[in,out]	d		local defect with values added on output
	 * \param[in]		u		local vector of current solution
	 * \param[in]		time	current time step
	 *
	 * \tparam	TElem 	Element type
	 * \tparam	TFVGeom	Finite Volume Geometry
	 */
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	adds the source part to the local defect
	/**
	 * This function adds the contribution of one element of the source part to
	 * the local defect, if a source function \f$ \vec{f} \f$ has been specified.
	 * Otherwise a zero source function is assumed.
	 *
	 * The contribution added is:
	 * \f{align*}
	 *  \forall d=0,\dots, \text{dim}-1, \quad
	 *  \forall sh=0,\dots, n_{sh}-1: \qquad
	 *	\mathbf{d}(d, sh) &\text{ +=} \int_{B_{sh} \cap T} \vec{f}_d \; dV
	 *					  \approx |B_{sh} \cap T| \; \vec{f}_d|_{ip}
	 * \f}
	 *
	 * \param[in,out]	d		local defect with values added on output
	 * \param[in]		time	current time step
	 *
	 * \tparam	TElem 	Element type
	 * \tparam	TFVGeom	Finite Volume Geometry
	 */
		template <typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	computes the Mach-number blended Upwind velocity
	/**
	 * This function compute the peclet blended velocity. First, the
	 * standard interpolated (no upwind) velocity is computed at the integration
	 * point of the scvf:
	 * \f{align*}
	 * 	\vec{u}^{\text{std}}_d|_{ip}
	 * 			:= \sum_{sh=0}^{n_{sh}-1} \mathbf{u}(d, sh) \; \phi_{sh}
	 * \f}
	 *
	 * Defining the Peclet Number
	 * \f{align*}
	 * 		Pe := \frac{\vec{u}^{\text{std}}|_{ip} \cdot \vec{n}|_{ip} \; L}{\nu}
	 * \f}
	 * with the Viscosity \f$ \nu \f$ and the diffusion length \f$ L \f$, a
	 * weighting factor is computed:
	 * \f{align*}
	 * 	\omega := \frac{Pe^2}{5 + Pe^2}
	 * \f}
	 *
	 * and the blended velocity is computed to
	 * \f{align*}
	 * 	\vec{u}^{\text{blend}}|_{ip} := \omega \vec{u}^{\text{up}}
	 * 									+ (1-\omega) \vec{u}^{\text{std}}
	 * \f}
	 *
	 * \param[in,out]	UpwindVel		upwind vel on entry, blended vel on exit
	 * \param[in]		scvf			SubControlVolumeFace of the Velocity
	 * \param[in]		StdVel			standard interpolation vel
	 * \param[in]		kinVisco		kinematic Viscosity at scvf
	 * \return			\f$\omega\f$ 	weighting factor
	 */
		template <typename TFVGeom>
		inline number mach_number_blending(MathVector<dim>& UpwindVel, const TFVGeom& geo, size_t ip,
		                           const MathVector<dim>& StdVel, number kinVisco);

	protected:
	///	Data import for source
		DataImport<MathVector<dim>, dim> m_imSourceSCV;
		DataImport<MathVector<dim>, dim> m_imSourceSCVF;

	///	Data import for kinematic viscosity
		DataImport<number, dim> m_imKinViscosity;

	///	Data import for adiabatic index
		DataImport<number, dim> m_imAdiabaticIndex;

	///	Upwinding for velocity in convective term of momentum equation
		SmartPtr<INavierStokesUpwind<dim> > m_spConvUpwind;

	/// abbreviation for pressure
		static const size_t _P_ = dim;
	/// abbreviation for density
		static const size_t _Rho_ = dim + 1;

		using base_type::m_bMachNrBlend;
		using base_type::m_bFullNewtonFactor;

		void init();

	///	reference quantities for getting dimensionless equations
		number m_maxVel;
		number m_maxPressure;
		number m_maxDensity;
		number m_refLength;

		number m_refMachNrSq;
		number m_numeratorOfRefReynoldsNr;

	private:
	///	register utils
	///	\{
		void register_all_funcs(bool bHang);
		template <typename TElem, typename TFVGeom> void register_func();
	/// \}
};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__COMPRESSIBLE_NAVIER_STOKES_FV1_H__*/
