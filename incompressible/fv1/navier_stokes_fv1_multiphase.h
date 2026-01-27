/*
 * Copyright (c) 2010-2014:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1M__NAVIER_STOKES_FV1M__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1M__NAVIER_STOKES_FV1M__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../incompressible_navier_stokes_base.h"
#include "../../upwind_interface.h"
#include "stabilization.h"
#include "../../properties_interface.h"

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
class NavierStokesFV1M
	: public IncompressibleNavierStokesBase<TDomain>
{
	protected:
	///	Base class type
		typedef IncompressibleNavierStokesBase<TDomain> base_type;

	///	own type
		typedef NavierStokesFV1M<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
	/// \{
		NavierStokesFV1M(const char* functions, const char* subsets);
		NavierStokesFV1M(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
	/// \}

	///	sets the kinematic viscosity
		void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user);

	///	returns kinematic viscosity
		SmartPtr<CplUserData<number, dim> > kinematic_viscosity() {return m_imKinViscosity.user_data ();}

	///	sets the density
		void set_density(SmartPtr<CplUserData<number, dim> > user);

	///	returns density
		SmartPtr<CplUserData<number, dim> > density() {return m_imDensitySCVF.user_data ();}

	///	sets the source function
		void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
    
    ///    returns source
        SmartPtr<CplUserData<MathVector<dim>, dim> > source() {return m_imSourceSCV.user_data ();}
    
        void set_relative_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
		void set_slip_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
        void set_diffusion(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > user);
    
        void set_average_gamma(SmartPtr<CplUserData<number, dim> > user);
        void set_interface_normal(SmartPtr<CplUserData<MathVector<dim>, dim> > user);


	///	sets the stabilization used to compute the stabilized velocities
		void set_stabilization(SmartPtr<INavierStokesFV1Stabilization<dim> > spStab)
			{m_spStab = spStab;}

	///	sets stabilization based on string identifier
		void set_stabilization(const std::string& name)
			{m_spStab = CreateNavierStokesStabilization<dim>(name);
			 if(m_spConvUpwind.valid()) m_spStab->set_upwind(m_spConvUpwind);
			 if(m_spConvUpwind_rel.valid()) m_spStab->set_upwind_rel(m_spConvUpwind_rel);}

	///	sets stabilization and diff length method based on string identifier (only for Schneider-Raw stabilizations)
		void set_stabilization(const std::string& name, const std::string& diffLength)
			{SmartPtr<INavierStokesSRFV1Stabilization<dim> > spStab = CreateNavierStokesStabilization<dim>(name);
			 spStab->set_diffusion_length(diffLength);
			 m_spStab = spStab;
			 if(m_spConvUpwind.valid()) m_spStab->set_upwind(m_spConvUpwind);
			 if(m_spConvUpwind_rel.valid()) m_spStab->set_upwind_rel(m_spConvUpwind_rel);}
			 
    /// returns stabilization	
		SmartPtr<INavierStokesFV1Stabilization<dim> > stabilization(){ return m_spStab;}

	///	sets a stabilization for upwinding (Physical Advection Correction)
        void set_upwind(SmartPtr<INavierStokesFV1Stabilization<dim> > spStab)
        	{m_spConvStab = spStab; m_spConvUpwind = SPNULL;}

	///	sets an upwinding for the convective term of momentum equation
		void set_upwind(SmartPtr<INavierStokesUpwind<dim> > spUpwind)
			{m_spConvStab = SPNULL; m_spConvUpwind = spUpwind;}

	///	sets the upwind based on a string identifier
		void set_upwind(const std::string& name)
			{m_spConvStab = SPNULL; m_spConvUpwind = CreateNavierStokesUpwind<dim>(name);
			if(m_spStab.valid() && m_spStab->upwind().invalid()) m_spStab->set_upwind(m_spConvUpwind);}
    ///    sets the upwind based on a string identifier
        void set_upwind_vol(const std::string& name)
            {m_spConvUpwind_vol = CreateNavierStokesUpwind<dim>(name);}
		void set_upwind_rel(const std::string& name)
			{m_spConvUpwind_rel = CreateNavierStokesUpwind<dim>(name);
			if(m_spStab.valid() && m_spStab->upwind_rel().invalid()) m_spStab->set_upwind_rel(m_spConvUpwind_rel);}

		void set_pac_upwind(bool bPac)
		{
			if (bPac==true){
				if (m_spConvUpwind.valid()) m_spStab->set_upwind(m_spConvUpwind);
				else UG_THROW("Upwind must be specified previously.\n");
				if (m_spStab.valid()) set_upwind(m_spStab);
				else UG_THROW("Stabilization must be specified previously.\n");
			}
		}
        void set_phase_parameters(Interface<dim>* user)
        {
            if(m_spStab.invalid())
                UG_THROW("Stabilization has not been set or must me set before.");
            if (!user->valid())
                UG_THROW("Interface parameters has not been initialized");
            m_spStab->set_phase_parameters(user);
            Inter = user;
        }
    
    ///    returns the export of the VolumeFractions
        SmartPtr<CplUserData<number, dim> > volume_fraction() {return m_exVolumeFraction;}
    ///    returns the export of the VolumeFraction gradient
        SmartPtr<CplUserData<MathVector<dim>, dim> > volume_fraction_grad() {return m_exVolumeFractionGrad;}
    ///    returns the export of the Einstein Viscosity
        SmartPtr<CplUserData<number, dim> > einstein_viscosity() {return m_exEinsteinViscosity;}
    ///    returns the export of the PariclePressure
        SmartPtr<CplUserData<number, dim> > particle_pressure() {return m_exPsPressure;}
    ///    returns the export of the ParticlePressure gradient
        SmartPtr<CplUserData<MathVector<dim>, dim> > particle_pressure_grad() {return m_exPsPressureGrad;}
    
			
	protected:

	using base_type::m_exVelocity;
	using base_type::m_exVelocityGrad;
	using base_type::m_exVelocity_ip;
    using base_type::m_exPressure;
    using base_type::m_exPressureGrad;
	
    ///    Export for the PariclePressure
    SmartPtr<DataExport<number,dim> > m_exPsPressure;
    ///    Export for the Particlepressure gradient
    SmartPtr<DataExport<MathVector<dim>,dim> > m_exPsPressureGrad;
    ///    Export for the VolumeFraction
    SmartPtr<DataExport<number,dim> > m_exVolumeFraction;
    ///    Export for the VolFraction gradient
    SmartPtr<DataExport<MathVector<dim>,dim> > m_exVolumeFractionGrad;
    ///    Export for the Einstein Viscosity
    SmartPtr<DataExport<number,dim> > m_exEinsteinViscosity;
    


	public:
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	///	returns string identifying disc type
		virtual std::string disc_type() const {return "fv1";};

	public:
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

	///	computes the pecled blended Upwind veloctity
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
		inline number peclet_blend(MathVector<dim>& UpwindVel, const TFVGeom& geo, size_t ip,
                                   const MathVector<dim>& StdVel, number kinVisco, number densitySCVF);
        template <typename TFVGeom>
        inline void std_vel(const LocalVector& u, const TFVGeom& geo, MathVector<dim>* StdVel, number* StdVol, MathVector<dim>* Vel, const DataImport<number, dim>& densitySCV, number* Rho_up, number* Rho_do, number* ConvRatio);
    
        template <typename TFVGeom>
        inline void std_rel_vel(const LocalVector& u, const TFVGeom& geo, MathVector<dim>* StdRelVel_div, MathVector<dim>* StdRelVel_total_t, MathVector<dim>* StdRelVel_t, const DataImport<number, dim>& densitySCV,  const DataImport<MathVector<dim>, dim>& RelVelSCV, const bool m_bStokes, const bool Jac);
    
        template <typename TFVGeom>
        inline void vel_grad(const LocalVector& u, const TFVGeom& geo, number* Gamma);

	///	export value of the velocity
		template <typename TElem, typename TFVGeom>
		void ex_nodal_velocity(MathVector<dim> vValue[],
							  const MathVector<dim> vGlobIP[],
							  number time, int si,
							  const LocalVector& u,
							  GridObject* elem,
							  const MathVector<dim> vCornerCoords[],
							  const MathVector<TFVGeom::dim> vLocIP[],
							  const size_t nip,
							  bool bDeriv,
							  std::vector<std::vector<MathVector<dim> > > vvvDeriv[]);
    ///    export value of the velocity
        template <typename TElem, typename TFVGeom>
        void ex_div_velocity(MathVector<dim> vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si,
                              const LocalVector& u,
                              GridObject* elem,
                              const MathVector<dim> vCornerCoords[],
                              const MathVector<TFVGeom::dim> vLocIP[],
                              const size_t nip,
                              bool bDeriv,
                              std::vector<std::vector<MathVector<dim> > > vvvDeriv[]);

	///	computes the value of the gradient of the velocity
		template <typename TElem, typename TFVGeom>
		void ex_velocity_grad(MathMatrix<dim,dim> vValue[],
							  const MathVector<dim> vGlobIP[],
							  number time, int si,
							  const LocalVector& u,
							  GridObject* elem,
							  const MathVector<dim> vCornerCoords[],
							  const MathVector<TFVGeom::dim> vLocIP[],
							  const size_t nip,
							  bool bDeriv,
							  std::vector<std::vector<MathMatrix<dim,dim> > > vvvDeriv[]);
    
    ///    export value of the pressure
        template <typename TElem, typename TFVGeom>
        void ex_nodal_pressure(number vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si,
                              const LocalVector& u,
                              GridObject* elem,
                              const MathVector<dim> vCornerCoords[],
                              const MathVector<TFVGeom::dim> vLocIP[],
                              const size_t nip,
                              bool bDeriv,
                              std::vector<std::vector<number > > vvvDeriv[]);
    
    ///    export value of the pressure gradient
        template <typename TElem, typename TFVGeom>
        void ex_pressure_grad(MathVector<dim> vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si,
                              const LocalVector& u,
                              GridObject* elem,
                              const MathVector<dim> vCornerCoords[],
                              const MathVector<TFVGeom::dim> vLocIP[],
                              const size_t nip,
                              bool bDeriv,
                              std::vector<std::vector<MathVector<dim> > > vvvDeriv[]);
    
    ///    export value of the volume fraction
        template <typename TElem, typename TFVGeom>
        void ex_nodal_volfraction(number vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si,
                              const LocalVector& u,
                              GridObject* elem,
                              const MathVector<dim> vCornerCoords[],
                              const MathVector<TFVGeom::dim> vLocIP[],
                              const size_t nip,
                              bool bDeriv,
                              std::vector<std::vector<number > > vvvDeriv[]);
    ///    export value of the VolumeFraction gradient
        template <typename TElem, typename TFVGeom>
        void ex_volume_fraction_grad(MathVector<dim> vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si,
                              const LocalVector& u,
                              GridObject* elem,
                              const MathVector<dim> vCornerCoords[],
                              const MathVector<TFVGeom::dim> vLocIP[],
                              const size_t nip,
                              bool bDeriv,
                              std::vector<std::vector<MathVector<dim> > > vvvDeriv[]);
    ///    export value of the einstein viscosity
        template <typename TElem, typename TFVGeom>
        void ex_nodal_einstein_viscosity(number vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si,
                              const LocalVector& u,
                              GridObject* elem,
                              const MathVector<dim> vCornerCoords[],
                              const MathVector<TFVGeom::dim> vLocIP[],
                              const size_t nip,
                              bool bDeriv,
                              std::vector<std::vector<number > > vvvDeriv[]);
    ///    export value of the partilce pressure
        template <typename TElem, typename TFVGeom>
        void ex_nodal_particle_pressure(number vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si,
                              const LocalVector& u,
                              GridObject* elem,
                              const MathVector<dim> vCornerCoords[],
                              const MathVector<TFVGeom::dim> vLocIP[],
                              const size_t nip,
                              bool bDeriv,
                              std::vector<std::vector<number > > vvvDeriv[]);
    ///    export value of the ParticlePressure gradient
        template <typename TElem, typename TFVGeom>
        void ex_particle_pressure_grad(MathVector<dim> vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si,
                              const LocalVector& u,
                              GridObject* elem,
                              const MathVector<dim> vCornerCoords[],
                              const MathVector<TFVGeom::dim> vLocIP[],
                              const size_t nip,
                              bool bDeriv,
                              std::vector<std::vector<MathVector<dim> > > vvvDeriv[]);
    
    
		
	protected:
	///	Data import for source
		DataImport<MathVector<dim>, dim> m_imSourceSCV;
		DataImport<MathVector<dim>, dim> m_imSourceSCVF;

	///	Data import for kinematic viscosity
		DataImport<number, dim> m_imKinViscosity;
        DataImport<number, dim> m_imKinViscosity_old;

	///	Data import for density
		DataImport<number, dim> m_imDensitySCVF;
        DataImport<number, dim> m_imDensitySCVF_old;
		DataImport<number, dim> m_imDensitySCV;
		DataImport<number, dim> m_imDensitySCV_old;
    
    ///    Data import for Diffusion
        DataImport<MathMatrix<dim,dim>, dim> m_imDiffusion;
    
    ///    Data import for Relative velocity
        DataImport<MathVector<dim>, dim> m_imRelativeVelocitySCV;
		DataImport<MathVector<dim>, dim> m_imRelativeVelocitySCVF;
	///    Data import for Slip velocity
		DataImport<MathVector<dim>, dim> m_imSlipVelocitySCVF;
	
	///    Data import for Nodal ShearRate
		DataImport<number, dim> m_imAverageGammaSCV;
    
    ///    Data import for multiphase flow
        DataImport<MathVector<dim>, dim> m_imSurfaceNormal;
        Interface<dim>* Inter = NULL;

	///	Stabilization for velocity in continuity equation
		SmartPtr<INavierStokesFV1Stabilization<dim> > m_spStab;
    
	///	Stabilization for velocity in convective term of momentum equation
	///	Here, the stabilization is used as an upwinding
		SmartPtr<INavierStokesFV1Stabilization<dim> > m_spConvStab;

	///	Upwinding for velocity in convective term of momentum equation
		SmartPtr<INavierStokesUpwind<dim> > m_spConvUpwind;
	///    Upwinding for Relative Velocity in convective term of momentum equation
		SmartPtr<INavierStokesUpwind<dim> > m_spConvUpwind_rel;
    ///    Upwinding for VolFraction in convective term of Transport equation
        SmartPtr<INavierStokesUpwind<dim> > m_spConvUpwind_vol;


	/// abbreviation for pressure and volume fraction
		static const size_t _P_ = dim;
        static const size_t _C_ = dim+1;

		using base_type::m_bPecletBlend;
		using base_type::m_bFullNewtonFactor;
        using base_type:: m_boolFullNewton;
		using base_type::m_bStokes;
		using base_type::m_bLaplace;
        using base_type::m_gradDivFactor;

		virtual void init();

    ///    computes the linearized defect w.r.t to the density SCV
        template <typename TElem, typename TFVGeom>
        void lin_def_densitySCV(const LocalVector& u,
                              std::vector<std::vector<number> > vvvLinDef[],
                              const size_t nip);
    ///    computes the linearized defect w.r.t to the density SCVF
        template <typename TElem, typename TFVGeom>
        void lin_def_densitySCVF(const LocalVector& u,
                              std::vector<std::vector<number> > vvvLinDef[],
                              const size_t nip);
    ///    computes the linearized defect w.r.t to the viscosity
        template <typename TElem, typename TFVGeom>
        void lin_def_viscosity(const LocalVector& u,
                              std::vector<std::vector<number> > vvvLinDef[],
                              const size_t nip);
    ///    computes the linearized defect w.r.t to the velocity
        template <typename TElem, typename TFVGeom>
        void lin_def_sourceSCV(const LocalVector& u,
                              std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                              const size_t nip);
    ///    computes the linearized defect w.r.t to the velocity
        template <typename TElem, typename TFVGeom>
        void lin_def_sourceSCVF(const LocalVector& u,
                              std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                              const size_t nip);
    
    
	protected:
	///	register utils
	///	\{
		virtual void register_all_funcs(bool bHang);
		template <typename TElem, typename TFVGeom> void register_func();
	/// \}
};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1M__NAVIER_STOKES_FV1M__*/
