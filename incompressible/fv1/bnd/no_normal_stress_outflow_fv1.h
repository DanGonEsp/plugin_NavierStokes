/*
 * Copyright (c) 2012-2014:  G-CSC, Goethe University Frankfurt
 * Authors: Dmitry Logashenko, Andreas Vogel
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

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__BND__NO_NORMAL_STRESS_OUTFLOW_FV1_
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__BND__NO_NORMAL_STRESS_OUTFLOW_FV1_

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../../bnd/no_normal_stress_outflow_base.h"

namespace ug{
namespace NavierStokes{


/// \ingroup lib_disc_elem_disc
/// @{

/// The zero-stress (neutral) outflow boundary condition for the incompressible NS equation
/**
 * This class implements the so-called neutral boundary condition for outflow
 * boundaries. Note that this class can be used only with the stabilized
 * vertex-centered discretization of the Navier-Stokes equations.
 *
 * This boundary condition imposes two equations on the unknown functions
 * at the outflow boundary \f$ F \f$:
 * <ul>
 * <li> \f$ \int_F p ds = 0 \f$, and
 * <li> \f$ \sigma \vec{n} \cdot \vec{n} = 0 \f$ on \f$ F \f$.
 * </ul>
 * where \f$ \vec{n} \f$ is the outer normal at \f$ F \f$ and
 * \f$ \sigma = \mu (\nabla \vec{u} + (\nabla \vec{u})^T) \f$ the stress tensor.
 */
template<	typename TDomain>
class NavierStokesNoNormalStressOutflowFV1
	: public NavierStokesNoNormalStressOutflowBase<TDomain>
{
	private:
	///	Base class type
		typedef NavierStokesNoNormalStressOutflowBase<TDomain> base_type;

	///	own type
		typedef NavierStokesNoNormalStressOutflowFV1<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
		NavierStokesNoNormalStressOutflowFV1(SmartPtr< IncompressibleNavierStokesBase<TDomain> > spMaster);

	protected:
	///	sets the kinematic viscosity
		virtual void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data)
			{m_imKinViscosity.set_data(data);}

	///	sets the density
		virtual void set_density(SmartPtr<CplUserData<number, dim> > data)
			{m_imDensity.set_data(data);}

	public:
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	public:
	///	prepares the element loop
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for evaluation
		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

	///	finishes the element loop
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();

	///	adds the stiffness part to the local jacobian
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	adds the stiffness part to the local defect
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	public:
	///	dummy implementations
	///	\{
		template <typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]){}
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]){}
		template <typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]){}
	/// \}

    
    ///    computes the linearized defect w.r.t to the density SCV
        template <typename TElem, typename TFVGeom>
        void lin_def_density(const LocalVector& u,
                              std::vector<std::vector<number> > vvvLinDef[],
                              const size_t nip);
    ///    computes the linearized defect w.r.t to the density SCV
        template <typename TElem, typename TFVGeom>
        void lin_def_viscosity(const LocalVector& u,
                              std::vector<std::vector<number> > vvvLinDef[],
                              const size_t nip);
	private:
	/// adds the diffusive part of the local Jacobian of the momentum equation
		template <typename BF>
		inline void diffusive_flux_Jac
		(
			const size_t ip,
			const BF& bf,
			LocalMatrix& J,
			const LocalVector& u
		);
	/// adds the diffusive part of the local defect of the momentum equation
		template <typename BF>
		inline void diffusive_flux_defect
		(
			const size_t ip,
			const BF& bf,
			LocalVector& d,
			const LocalVector& u
		);
	/// adds the convective part of the local Jacobian of the momentum equation
		template <typename BF>
		inline void convective_flux_Jac
		(
			const size_t ip,
			const BF& bf,
			LocalMatrix& J,
			const LocalVector& u,
            const MathVector<dim>& StdVel
		);
	/// adds the convective part of the local defect of the momentum equation
		template <typename BF>
		inline void convective_flux_defect
		(
			const size_t ip,
			const BF& bf,
			LocalVector& d,
			const LocalVector& u,
			const MathVector<dim>& StdVel
		);
    /// adds the pressure part of the local Jacobian of the momentum equation
        template <typename BF>
        inline void pressure_flux_Jac
        (
            const size_t ip,
            const BF& bf,
            LocalMatrix& J,
            const LocalVector& u
        );
    /// adds the pressure  part of the local defect of the momentum equation
        template <typename BF>
        inline void pressure_flux_defect
        (
            const size_t ip,
            const BF& bf,
            LocalVector& d,
            const LocalVector& u,
            const number& pressure,
            const MathVector<dim>& StdVel
        );
    /// adds the convective part of the local defect of the momentum equation
        template <typename BF>
        inline void convective_flux_lin_defect
        (
            const size_t ip,
            const BF& bf,
         std::vector<std::vector<number> > vvvLinDef[],
            const LocalVector& u
        );
    /// adds the convective part of the local defect of the momentum equation
        template <typename BF>
        inline void diffusive_flux_lin_defect
        (
            const size_t ip,
            const BF& bf,
         std::vector<std::vector<number> > vvvLinDef[],
            const LocalVector& u
        );
    
	protected:
	/// abbreviation for pressure
		static const size_t _P_ = dim;

		using base_type::m_spMaster;
		using base_type::m_vBndSubSetIndex;

	///	Data import for kinematic viscosity
		DataImport<number, dim> m_imKinViscosity;

	/// Data import for density
		DataImport<number, dim> m_imDensity;

	/// Boundary integration points of the viscosity and the density
		std::vector<MathVector<dim> > m_vLocIP;
		std::vector<MathVector<dim> > m_vGloIP;

	protected:
		void register_all_funcs(bool bHang);
		template<typename TElem, typename TFVGeom>
		void register_func();

};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__NO_NORMAL_STRESS_OUTFLOW_FV1_*/
