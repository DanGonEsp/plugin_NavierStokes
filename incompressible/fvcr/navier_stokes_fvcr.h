/*
 * Copyright (c) 2010-2017:  G-CSC, Goethe University Frankfurt
 * Author: Christian Wehner
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

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FVCR__NAVIER_STOKES_FVCR__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FVCR__NAVIER_STOKES_FVCR__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../incompressible_navier_stokes_base.h"
#include "../../upwind_interface.h"

namespace ug{
namespace NavierStokes{

/// \ingroup lib_disc_elem_disc
/// @{

template<	typename TDomain>
class NavierStokesFVCR
	: public IncompressibleNavierStokesBase<TDomain>
{
	private:
	///	Base class type
		typedef IncompressibleNavierStokesBase<TDomain> base_type;

	///	own type
		typedef NavierStokesFVCR<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
	/// \{
		NavierStokesFVCR(const char* functions, const char* subsets);
		NavierStokesFVCR(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
	/// \}

	///	sets the kinematic viscosity
		void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user);

	///	returns kinematic viscosity
		SmartPtr<CplUserData<number, dim> > kinematic_viscosity() {
            return m_imKinViscosity.user_data (); }

	///	sets the density
		void set_density(SmartPtr<CplUserData<number, dim> > user);

	///	returns density
		SmartPtr<CplUserData<number, dim> > density() {return m_imDensitySCVF.user_data ();}

	///	sets the source function
		void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
        void set_relative_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
        void set_divergence_rel_vel(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
        void set_mass_change(SmartPtr<CplUserData<number, dim> > user);
        void set_source_surface(SmartPtr<CplUserData<number, dim> > user);
		void set_defect_upwind(bool defectUpwind) { m_bDefectUpwind = defectUpwind;}
		bool get_defect_upwind() {return m_bDefectUpwind; }

	///	sets an upwinding for the convective term of momentum equation
		void set_upwind(SmartPtr<INavierStokesUpwind<dim> > spUpwind)
			{m_spConvUpwind = spUpwind;}

	///	sets the upwind based on a string identifier
		void set_upwind(const std::string& name)
			{m_spConvUpwind = CreateNavierStokesUpwind<dim>(name);}

	public:
	///	type of trial space for each function used
		void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);
		
	///	returns if hanging nodes are needed
		virtual bool use_hanging() const;

	///	returns string identifying disc type
		virtual std::string disc_type() const {return "fvcr";};

	protected:
		void init();

	public:
	/// mixed upwind for Crouzeix-Raviart elements
		template <typename TFVGeom>
		inline number peclet_blend(MathVector<dim>& UpwindVel, const TFVGeom& geo, size_t ip,
		                           const MathVector<dim>& StdVel, number kinVisco);

	protected:
	/// flag if using upwind in defect computation
		bool m_bDefectUpwind;

	///	Data import for source
		DataImport<MathVector<dim>, dim> m_imSource;
    
    ///    Data import for Relative velocity
        DataImport<MathVector<dim>, dim> m_imRelativeVelocity;
        DataImport<MathVector<dim>, dim> m_imDivergenceFlux;
    
    ///    Data import for mass term
        DataImport<number, dim> m_imMass;
        DataImport<number, dim> m_imSourceSurface;

	///	Data import for kinematic viscosity
		DataImport<number, dim> m_imKinViscosity;

	///	Data import for density
		DataImport<number, dim> m_imDensitySCVF;
		DataImport<number, dim> m_imDensitySCV;

	/// Data import for central gradient
		DataImport<MathMatrix<dim,dim>, dim> m_imCentralGradient;

	/// Data import for central gradient
		DataImport<MathVector<dim>, dim> m_imPressureGradient;

	///	Upwinding for velocity in convective term of momentum equation
		SmartPtr<INavierStokesUpwind<dim> > m_spConvUpwind;

	/// abbreviation for pressure
		static const size_t _P_ = dim;

		using base_type::m_bPecletBlend;
		using base_type::m_bFullNewtonFactor;
		using base_type::m_bStokes;
		using base_type::m_bLaplace;
		using base_type::m_gradDivFactor;
    
    protected:

    using base_type::m_exVelocity;
    using base_type::m_exVelocityGrad;
    using base_type::m_exPressure;
    using base_type::m_exPressureGrad;

	public:
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();
		
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);
		
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template <typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template <typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);
    
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
    ///    computes the linearized defect w.r.t to the viscosity
        template <typename TElem, typename TFVGeom>
        void lin_def_sourceSCV(const LocalVector& u,
                              std::vector<std::vector<MathVector<dim>> > vvvLinDef[],
                              const size_t nip);
    ///    computes the linearized defect w.r.t to the viscosity
        template <typename TElem, typename TFVGeom>
        void lin_def_sourceSCVF(const LocalVector& u,
                              std::vector<std::vector<number> > vvvLinDef[],
                              const size_t nip);
        template<typename TElem, typename TFVGeom>
        void lin_def_RelativeVelocity(const LocalVector& u,
                             std::vector<std::vector<MathVector<dim>> > vvvLinDef[],
                                const size_t nip);
    ///    computes the linearized defect w.r.t to the Relative Mass
        template <typename TElem, typename TFVGeom>
        void lin_def_relativeMass(const LocalVector& u,
                              std::vector<std::vector<number> > vvvLinDef[],
                              const size_t nip);
        template <typename TElem, typename TFVGeom>
        void lin_def_divergence_flux(const LocalVector& u,
                              std::vector<std::vector<MathVector<dim>> > vvvLinDef[],
                              const size_t nip);
    
    ///    export value of the velocity
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

    ///    computes the value of the gradient of the velocity
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

	private:
		template <typename TElem, typename TFVGeom>
		void register_func();
		void register_all_funcs(bool bHang);
};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FVCR__NAVIER_STOKES_FVCR__*/
