/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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


#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "navier_stokes_fvcr.h"

#include "pressure_separation.h"

#include "disc_constraint_fvcr.h"

#include "turbulent_viscosity_fvcr.h"

#include "bnd/inflow_fvcr.h"
#include "bnd/no_normal_stress_outflow_fvcr.h"

#include "bnd/symmetric_boundary_fvcr.h"
#include "cr_reorder.h"
#include "cr_ilut.h"
#include "pcr_ilut.h"

#include "lib_disc/function_spaces/grid_function.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace NavierStokes{

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct FunctionalityFVCR
{

/**
 * Function called for the registration of Domain and Algebra dependent parts
 * of the plugin. All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	//static const int dim = TDomain::dim;
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, TAlgebra> function_type;

	//	NavierStokesInflow FVCR
	{
		typedef NavierStokesInflowFVCR<TDomain, TAlgebra> T;
		typedef NavierStokesInflowBase<TDomain, TAlgebra> TBase;
		string name = string("NavierStokesInflowFVCR").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokesFVCR<TDomain> >)>("MasterElemDisc")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesInflowFVCR", tag);
	}
/*-- The following segment has been commented out (by the author?). Use with care! --*/
	//	Order CR-Cuthill-McKee
	{
		reg.add_function("OrderCRCuthillMcKee", static_cast<void (*)(approximation_space_type&,function_type&, bool)>(&OrderCRCuthillMcKee), grp);
	}
	
	//	Order CR-Cuthill-McKee
	{
		reg.add_function("CROrderCuthillMcKee", static_cast<void (*)(approximation_space_type&,function_type&, bool,bool,bool,bool)>(&CROrderCuthillMcKee), grp);
	}
	
	//	Order CR-Sloan
	{
		reg.add_function("CROrderSloan", static_cast<void (*)(approximation_space_type&,function_type&, bool,bool,bool)>(&CROrderSloan), grp);
	}
	
	//	Order CR-King
	{
		reg.add_function("CROrderKing", static_cast<void (*)(approximation_space_type&,function_type&, bool,bool,bool,bool)>(&CROrderKing), grp);
	}
	
	//	Order CR-Minimum Degree
	{
		reg.add_function("CROrderMinimumDegree", static_cast<void (*)(approximation_space_type&,function_type&, bool,bool,bool)>(&CROrderMinimumDegree), grp);
	}
/*-- End of the segment. --*/
	typedef ug::GridFunction<TDomain, TAlgebra> TFct;
	static const int dim = TDomain::dim;

	// Turbulent viscosity data
	// Smagorinsky model
	{
		string name = string("CRSmagorinskyTurbViscData").append(suffix);
		typedef CRSmagorinskyTurbViscData<TFct> T;
		typedef CplUserData<number, dim> TBase;
		typedef INewtonUpdate TBase2;
		reg.add_class_<T, TBase,TBase2>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,SmartPtr<TFct>,number)>("Approximation space, grid function, model parameter")
			.add_method("set_model_parameter", &T::set_model_parameter)
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
		#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
		#endif
			.add_method("set_turbulence_zero_bnd", &T::setTurbulenceZeroBoundaries)
			.add_method("update", &T::update)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRSmagorinskyTurbViscData", tag);
	}

	// Dynamic model
	{
		string name = string("CRDynamicTurbViscData").append(suffix);
		typedef CRDynamicTurbViscData<TFct> T;
		typedef CplUserData<number, dim> TBase;
		typedef INewtonUpdate TBase2;
		reg.add_class_<T, TBase,TBase2>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,SmartPtr<TFct>)>("Approximation space, grid function")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
		#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
		#endif
			.add_method("set_turbulence_zero_bnd", &T::setTurbulenceZeroBoundaries)
			.add_method("set_time_filter", &T::set_time_filter)
			.add_method("set_time_filter_eps", &T::set_time_filter_eps)
			.add_method("set_space_filter", &T::set_space_filter)
			.add_method("update", &T::update)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRDynamicTurbViscData", tag);
	}

	// SeparatedPressureSource
	{
		string name = string("SeparatedPressureSource").append(suffix);
		typedef SeparatedPressureSource<TFct> T;
		typedef CplUserData<MathVector<dim>, dim> TBase;
		typedef INewtonUpdate TBase2;
		reg.add_class_<T, TBase,TBase2>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,SmartPtr<TFct>)>("Approximation space, grid function")
				.add_method("set_source", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_source), "", "Source")
				.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "F_x")
				.add_method("set_source", static_cast<void (T::*)(number,number)>(&T::set_source), "", "F_x, F_y")
				.add_method("set_source", static_cast<void (T::*)(number,number,number)>(&T::set_source), "", "F_x, F_y, F_z")
			#ifdef UG_FOR_LUA
				.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source Vector")
			#endif
				.add_method("update", &T::update)
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SeparatedPressureSource", tag);
	}
	
	//	DiscConstraintFVCR
	{
		typedef DiscConstraintFVCR<TFct> T;
		typedef IDomainConstraint<typename TFct::domain_type,typename TFct::algebra_type> TBase;
		string name = string("DiscConstraintFVCR").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
		.template add_constructor<void (*)(SmartPtr<TFct>)>("Grid function")
		// bool bLinUpConvDefect,bool bLinUpConvJacobian,bool bLinPressureDefect,bool bLinPressureJacobian,bool bAdaptive
		.template add_constructor<void (*)(SmartPtr<TFct>,bool,bool,bool,bool,bool)>("Grid function,lin up def,lin up jac,lin p def,lin p jac,adaptivity")
		.template add_constructor<void (*)(SmartPtr<TFct>,bool,bool,bool,bool,bool,const char*)>("Grid function,lin up def,lin up jac,lin p def,lin p jac,adaptivity,bnd subsets")
		.template add_constructor<void (*)(SmartPtr<TFct>,bool,bool,bool,bool,bool,bool)>("Grid function,lin up def,lin up jac,lin p def,lin p jac,adaptivity,limiter")
		.template add_constructor<void (*)(SmartPtr<TFct>,bool,bool,bool,bool,bool,bool,const char*)>("Grid function,lin up def,lin up jac,lin p def,lin p jac,adaptivity,limiter,bnd subsets")
		.add_method("set_zero_grad_bnd", &T::set_zero_grad_bnd)
		.add_method("set_limiter", &T::set_limiter)
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DiscConstraintFVCR", tag);
	}
}

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();
	
	//	CR ILU Threshold
	{
		typedef CRILUTPreconditioner<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("CRILUT").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Incomplete LU Decomposition with threshold")
		.add_constructor()
		.template add_constructor<void (*)(number)>("threshold parameter")
		.template add_constructor<void (*)(number,number)>("threshold parameters")
		.template add_constructor<void (*)(number,number,number,number)>("threshold parameters")
		.template add_constructor<void (*)(number,bool)>("threshold parameter,storage info output")
		.template add_constructor<void (*)(number,number,bool)>("threshold vv,threshold vp/pv/pp,storage info output")
		.template add_constructor<void (*)(number,number,number,number,bool)>("threshold vv,vp,pv,pp,storage info output")
		.add_method("set_threshold",static_cast<void (T::*)(number,number,number,number)>(&T::set_threshold),
					"", "threshold", "sets threshold of incomplete LU factorisation")
		.add_method("set_threshold",static_cast<void (T::*)(number,number)>(&T::set_threshold),
					"", "threshold", "sets threshold of incomplete LU factorisation")
		.add_method("set_threshold",static_cast<void (T::*)(number)>(&T::set_threshold),
					"", "threshold", "sets threshold of incomplete LU factorisation")
		.add_method("set_info", &T::set_info,
					"", "info", "sets storage information output")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRILUT", tag);
	}
	
	//	parallel CR ILU Threshold
	{
		typedef PCRILUTPreconditioner<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("PCRILUT").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Incomplete LU Decomposition with threshold")
		.add_constructor()
		.template add_constructor<void (*)(number)>("threshold parameter")
		.template add_constructor<void (*)(number,number)>("threshold parameters")
		.template add_constructor<void (*)(number,number,number,number)>("threshold parameters")
		.template add_constructor<void (*)(number,bool)>("threshold parameter,storage info output")
		.template add_constructor<void (*)(number,number,bool)>("threshold vv,threshold vp/pv/pp,storage info output")
		.template add_constructor<void (*)(number,number,number,number,bool)>("threshold vv,vp,pv,pp,storage info output")
		.add_method("set_threshold",static_cast<void (T::*)(number,number,number,number)>(&T::set_threshold),
					"", "threshold", "sets threshold of incomplete LU factorisation")
		.add_method("set_threshold",static_cast<void (T::*)(number,number)>(&T::set_threshold),
					"", "threshold", "sets threshold of incomplete LU factorisation")
		.add_method("set_threshold",static_cast<void (T::*)(number)>(&T::set_threshold),
					"", "threshold", "sets threshold of incomplete LU factorisation")
		.add_method("set_info", &T::set_info,
					"", "info", "sets storage information output")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PCRILUT", tag);
	}
}

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	static const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

	//	Navier-Stokes FVCR
	{
		typedef NavierStokesFVCR<TDomain> T;
		typedef IncompressibleNavierStokesBase<TDomain> TBase;
		string name = string("NavierStokesFVCR").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Functions#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesUpwind<dim> >)>(&T::set_upwind))
			.add_method("set_upwind",  static_cast<void (T::*)(const std::string&)>(&T::set_upwind))
			.add_method("set_defect_upwind", &T::set_defect_upwind)
            .add_method("set_mass_change", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_mass_change))
            .add_method("set_source_surface", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_source_surface))
            .add_method("set_density_upwind", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_density_upwind))
            .add_method("set_relative_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_relative_velocity))
            .add_method("set_divergence_rel_vel", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_divergence_rel_vel))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFVCR", tag);
	}

	//	NavierStokesNoNormalStressOutflow FVCR
	{
		typedef NavierStokesNoNormalStressOutflowFVCR<TDomain> T;
		typedef NavierStokesNoNormalStressOutflowBase<TDomain> TBase;
		string name = string("NavierStokesNoNormalStressOutflowFVCR").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< IncompressibleNavierStokesBase<TDomain> >)>("MasterDisc")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesNoNormalStressOutflowFVCR", tag);
	}

	//	CRNavierStokesSymBC
	{
		typedef CRNavierStokesSymBC<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("CRNavierStokesSymBC").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< IncompressibleNavierStokesBase<TDomain> >)>("MasterDisc")
			.add_method("add", &T::add, "", "Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRNavierStokesSymBC", tag);
	}
}

/**
 * Function called for the registration of Dimension dependent parts
 * of the plugin. All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();
}

}; // end Functionality
} // end namespace NavierStokes


/**
 * This function is called when the plugin is loaded.
 */
void Init___NavierStokes___FVCR(Registry* reg, string grp)
{
	grp.append("SpatialDisc/NavierStokes/");
	typedef NavierStokes::FunctionalityFVCR Functionality;

	try{
//		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
		RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug

