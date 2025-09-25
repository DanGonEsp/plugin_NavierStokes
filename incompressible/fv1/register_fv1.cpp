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


#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "navier_stokes_fv1.h"
#include "bnd/inflow_fv1.h"
#include "bnd/no_normal_stress_outflow_fv1.h"
#include "bnd/symmetric_boundary_fv1.h"
#include "bnd/inflow_stress_fv1.h"
#include "bnd/wall_sliding_fv1.h"
#include "stabilization.h"
#include "pressure_jump.h"
#include "nodal_pressure_grad.h"

#include "turbulent_viscosity_fv1.h"

#include "lib_disc/function_spaces/grid_function.h"

#include "navier_stokes_fv1_multiphase.h"
#include "bnd/inflow_fv1_multiphase.h"
#include "bnd/no_normal_stress_outflow_fv1_multiphase.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace NavierStokes{

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct FunctionalityFV1
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

	//	NavierStokesInflow FV1
	{
		typedef NavierStokesInflowFV1<TDomain, TAlgebra> T;
		typedef NavierStokesInflowBase<TDomain, TAlgebra> TBase;
		string name = string("NavierStokesInflowFV1").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokesFV1<TDomain> >)>("MasterElemDisc")
 			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesInflowFV1", tag);
	}
    //    NavierStokesInflow FV1M
    {
        typedef NavierStokesInflowFV1M<TDomain, TAlgebra> T;
        typedef NavierStokesInflowBase<TDomain, TAlgebra> TBase;
        string name = string("NavierStokesInflowFV1M").append(suffix);
        reg.add_class_<T, TBase>(name, grp)
            .template add_constructor<void (*)(SmartPtr< NavierStokesFV1M<TDomain> >)>("MasterElemDisc")
             .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "NavierStokesInflowFV1M", tag);
    }
	
	typedef ug::GridFunction<TDomain, TAlgebra> TFct;
	static const int dim = TDomain::dim;

	// Turbulent viscosity data
	// Smagorinsky model
	{
		string name = string("FV1SmagorinskyTurbViscData").append(suffix);
		typedef FV1SmagorinskyTurbViscData<TFct> T;
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
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FV1SmagorinskyTurbViscData", tag);
	}

	// Dynamic model
	{
		string name = string("FV1DynamicTurbViscData").append(suffix);
		typedef FV1DynamicTurbViscData<TFct> T;
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
			.add_method("set_space_filter", &T::set_time_filter)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FV1DynamicTurbViscData", tag);
	}
    // NodalPressureGradient
    {
        string name = string("NodalPressureGradient").append(suffix);
        typedef NodalPressureGradient<TFct> T;
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
        reg.add_class_to_group(name, "NodalPressureGradient", tag);
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

	//	Navier-Stokes FV1
	{
		typedef NavierStokesFV1<TDomain> T;
		typedef IncompressibleNavierStokesBase<TDomain> TBase;
		string name = string("NavierStokesFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Functions#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_stabilization",  static_cast<void (T::*)(SmartPtr<INavierStokesFV1Stabilization<dim> >)>(&T::set_stabilization))
			.add_method("set_stabilization",  static_cast<void (T::*)(const std::string&)>(&T::set_stabilization))
			.add_method("set_stabilization",  static_cast<void (T::*)(const std::string&, const std::string&)>(&T::set_stabilization))
			.add_method("set_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesFV1Stabilization<dim> >)>(&T::set_upwind))
			.add_method("set_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesUpwind<dim> >)>(&T::set_upwind))
			.add_method("set_upwind",  static_cast<void (T::*)(const std::string&)>(&T::set_upwind))
            .add_method("set_relative_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_relative_velocity), "", "RelativeVel")
            .add_method("set_divergence_rel_vel", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_divergence_rel_vel), "", "RelVel Divergence")
            .add_method("set_mass_change", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_mass_change), "", "Mass")
            .add_method("set_particle_pressure", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_particle_pressure), "", "ParticlePressure")
			.add_method("set_pac_upwind", &T::set_pac_upwind, "", "Set pac upwind")
            .add_method("set_density_ref", static_cast<void (T::*)(number)>(&T::set_density_ref), "", "density_ref")
            .add_method("set_interface_value", static_cast<void (T::*)(number)>(&T::set_interface_value), "", "set_interface_value")
            .add_method("set_jump_shape", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>, const std::string&)>(&T::set_jump_shape), "", "JumpShape")
            .add_method("set_interface_normal", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_interface_normal), "", "SurfaceNormal")
            .add_method("set_vol_fraction", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> > )>(&T::set_vol_fraction), "", "VolFraction")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFV1", tag);
	}
    //    Navier-Stokes FV1M
    {
        typedef NavierStokesFV1M<TDomain> T;
        typedef IncompressibleNavierStokesBase<TDomain> TBase;
        string name = string("NavierStokesFV1M").append(suffix);
        reg.add_class_<T, TBase >(name, grp)
            .template add_constructor<void (*)(const char*,const char*)>("Functions#Subset(s)")
            .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
            .add_method("set_stabilization",  static_cast<void (T::*)(SmartPtr<INavierStokesFV1Stabilization<dim> >)>(&T::set_stabilization))
            .add_method("set_stabilization",  static_cast<void (T::*)(const std::string&)>(&T::set_stabilization))
            .add_method("set_stabilization",  static_cast<void (T::*)(const std::string&, const std::string&)>(&T::set_stabilization))
            .add_method("set_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesFV1Stabilization<dim> >)>(&T::set_upwind))
            .add_method("set_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesUpwind<dim> >)>(&T::set_upwind))
            .add_method("set_upwind",  static_cast<void (T::*)(const std::string&)>(&T::set_upwind))
            .add_method("set_upwind_vol",  static_cast<void (T::*)(const std::string&)>(&T::set_upwind_vol))
            .add_method("set_upwind",  static_cast<void (T::*)(SmartPtr<IConvectionShapes<dim> >)>(&T::set_upwind))
            .add_method("set_relative_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_relative_velocity), "", "RelativeVel")
            .add_method("set_divergence_rel_vel", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_divergence_rel_vel), "", "RelVel Divergence")
            .add_method("set_mass_change", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_mass_change), "", "Mass")
            .add_method("set_particle_pressure", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_particle_pressure), "", "ParticlePressure")
            .add_method("set_pac_upwind", &T::set_pac_upwind, "", "Set pac upwind")
            .add_method("set_density_ref", static_cast<void (T::*)(number)>(&T::set_density_ref), "", "density_ref")
            .add_method("set_interface_value", static_cast<void (T::*)(number)>(&T::set_interface_value), "", "set_interface_value")
            .add_method("set_jump_shape", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>, const std::string&)>(&T::set_jump_shape), "", "JumpShape")
            .add_method("set_interface_normal", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_interface_normal), "", "SurfaceNormal")
            .add_method("set_vol_fraction", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> > )>(&T::set_vol_fraction), "", "VolFraction")
            .add_method("volume_fraction", &T::volume_fraction)
            .add_method("volume_fraction_average", &T::volume_fraction_average)
            .add_method("volume_fraction_grad", &T::volume_fraction_grad)
            .add_method("particle_pressure", &T::particle_pressure)
            .add_method("set_phase_parameters", &T::set_phase_parameters)
            .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "NavierStokesFV1M", tag);
    }


	//	NavierStokesNoNormalStressOutflow FV1
	{
		typedef NavierStokesNoNormalStressOutflowFV1<TDomain> T;
		typedef NavierStokesNoNormalStressOutflowBase<TDomain> TBase;
		string name = string("NavierStokesNoNormalStressOutflowFV1").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< IncompressibleNavierStokesBase<TDomain> >)>("MasterDisc")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesNoNormalStressOutflowFV1", tag);
	}
    //    NavierStokesNoNormalStressOutflow FV1M
    {
        typedef NavierStokesNoNormalStressOutflowFV1M<TDomain> T;
        typedef NavierStokesNoNormalStressOutflowBase<TDomain> TBase;
        string name = string("NavierStokesNoNormalStressOutflowFV1M").append(suffix);
        reg.add_class_<T, TBase>(name, grp)
            .template add_constructor<void (*)(SmartPtr< IncompressibleNavierStokesBase<TDomain> >)>("MasterDisc")
            .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "NavierStokesNoNormalStressOutflowFV1M", tag);
    }

	//	NavierStokesSymBCFV1
	{
		typedef NavierStokesSymBCFV1<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("NavierStokesSymBCFV1").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< IncompressibleNavierStokesBase<TDomain> >)>("MasterDisc")
			.add_method("add", &T::add, "", "Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesSymBCFV1", tag);
	}
    //    NavierStokesSymBCFV1
    {
        typedef NavierStokesInflowStressFV1<TDomain> T;
        typedef IElemDisc<TDomain> TBase;
        string name = string("NavierStokesInflowStressFV1").append(suffix);
        reg.add_class_<T, TBase>(name, grp)
            .template add_constructor<void (*)(SmartPtr< IncompressibleNavierStokesBase<TDomain> >)>("MasterDisc")
            .add_method("add", &T::add, "", "Subset(s)")
            .add_method("set_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_velocity), "", "VelocityBC")
            .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "NavierStokesInflowStressFV1", tag);
    }
    
    

	//	NavierStokesWSBCFV1
	{
		typedef NavierStokesWSBCFV1<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("NavierStokesWSBCFV1").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< IncompressibleNavierStokesBase<TDomain> >)>("MasterDisc")
			.add_method("add", &T::add, "", "Subset(s)")
			.add_method("set_sliding_limit", static_cast<void (T::*)(number)>(&T::set_sliding_limit), "", "SlidingLimit")
			.add_method("set_sliding_factor", static_cast<void (T::*)(number)>(&T::set_sliding_factor), "", "SlidingFactor")
#ifdef UG_FOR_LUA
			.add_method("set_sliding_limit", static_cast<void (T::*)(const char*)>(&T::set_sliding_limit), "", "SlidingLimit")
			.add_method("set_sliding_factor", static_cast<void (T::*)(const char*)>(&T::set_sliding_factor), "", "SlidingFactor")
#endif
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesWSBCFV1", tag);
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

//	INavierStokesFV1Stabilization
	{
		typedef INavierStokesFV1Stabilization<dim> T;
		string name = string("INavierStokesFV1Stabilization").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_method("set_upwind", &T::set_upwind);
		reg.add_class_to_group(name, "INavierStokesFV1Stabilization", tag);
	}

//	INavierStokesSRFV1Stabilization
	{
		typedef INavierStokesSRFV1Stabilization<dim> T;
		typedef INavierStokesFV1Stabilization<dim> TBase;
		string name = string("INavierStokesSRFV1Stabilization").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_diffusion_length", &T::set_diffusion_length);
		reg.add_class_to_group(name, "INavierStokesSRFV1Stabilization", tag);
	}

//	NavierStokesFIELDSStabilization
	{
		typedef NavierStokesFIELDSStabilization<dim> T;
		typedef INavierStokesSRFV1Stabilization<dim> TBase;
		string name = string("NavierStokesFIELDSStabilization").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFIELDSStabilization", tag);
	}

//	NavierStokesFLOWStabilization
	{
		typedef NavierStokesFLOWStabilization<dim> T;
		typedef INavierStokesSRFV1Stabilization<dim> TBase;
		string name = string("NavierStokesFLOWStabilization").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFLOWStabilization", tag);
	}

//	NavierStokesFV1WithoutStabilization
	{
		typedef NavierStokesFV1WithoutStabilization<dim> T;
		typedef INavierStokesFV1Stabilization<dim> TBase;
		string name = string("NavierStokesFV1WithoutStabilization").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFV1WithoutStabilization", tag);
	}
//    INavierStokesPressureJump
    {
        typedef INavierStokesPressureJump<dim> T;
        string name = string("INavierStokesPressureJump").append(suffix);
        reg.add_class_<T>(name, grp);
        reg.add_class_to_group(name, "INavierStokesPressureJump", tag);
    }

//    NavierStokesFIELDSStabilization
    {
        typedef NavierStokesViscousPressureJump<dim> T;
        typedef INavierStokesPressureJump<dim> TBase;
        string name = string("NavierStokesViscousPressureJump").append(suffix);
        reg.add_class_<T, TBase>(name, grp)
            .add_constructor()
            .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "NavierStokesViscousPressureJump", tag);
    }
}

}; // end Functionality
} // end namespace NavierStokes


/**
 * This function is called when the plugin is loaded.
 */
void Init___NavierStokes___FV1(Registry* reg, string grp)
{
	grp.append("SpatialDisc/NavierStokes/");
	typedef NavierStokes::FunctionalityFV1 Functionality;

	try{
		RegisterDimension2d3dDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
//		RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
