/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "register_navier_stokes.h"
#include "navier_stokes_base.h"
#include "upwind.h"

#include "properties_interface.h"
#include "linker/granular_density_linker.h"
#include "linker/granular_viscosity_linker.h"
#include "linker/granular_source_linker.h"
#include "linker/relative_velocity_linker.h"
#include "linker/granular_diffusion_linker.h"
#include "linker/scale_linker.h"


using namespace std;
using namespace ug::bridge;

namespace ug{
namespace NavierStokes{

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
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
{}


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
{}

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

//    check whether those classes have already been registered
   {
	   string name = string("NavierStokesBase").append(suffix);
	   if(reg.get_class(name))
		   return;
   }

//	Navier-Stokes Base
	{
		typedef NavierStokesBase<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("NavierStokesBase").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
#endif

			.add_method("set_source", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(const std::vector<number>&)>(&T::set_source), "", "Source")
#ifdef UG_FOR_LUA
			.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source")
#endif

			.add_method("disc_type", &T::disc_type)
			.add_method("set_exact_jacobian", static_cast<void (T::*)(bool)>(&T::set_exact_jacobian), "", "ExactJacobian")
			.add_method("set_exact_jacobian", static_cast<void (T::*)(number)>(&T::set_exact_jacobian), "", "ExactJacobianFactor");
		reg.add_class_to_group(name, "NavierStokesBase", tag);
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
	
//    PropertiesInterface
	{
		typedef Interface<dim> T;
		string name = string("Interface").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_method("set_particle_density",  &T::set_particle_density)
			.add_method("set_air_density",  &T::set_air_density)
			.add_method("set_fluid_Visc",  &T::set_fluid_Visc)
			.add_method("set_particle_diameter",  &T::set_particle_diameter)
			.add_method("set_alpha_max",  &T::set_alpha_max)
			.add_method("set_alpha_min",  &T::set_alpha_min)
			.add_method("set_FR",  &T::set_FR)
			.add_method("set_B_phi",  &T::set_B_phi)
			.add_method("set_deltaGamma",  &T::set_deltaGamma)
			.add_method("set_limit",  &T::set_limit)
			.add_method("set_time_step_factor",  &T::set_time_step_factor)
			.add_method("set_interface_volume_fraction",  &T::set_interface_volume_fraction)
			.add_method("set_bool_particle_pressure_force",  &T::set_bool_particle_pressure_force)
			.add_method("set_bool_consistent_gravity",  &T::set_bool_consistent_gravity)
			.add_method("set_bool_initialized",  &T::set_bool_initialized)
			.add_method("set_drag_model",  &T::set_drag_model)
			.add_method("set_reference_pressure",  &T::set_reference_pressure)
			.add_method("RelVel_ext",  &T::RelVel_ext)
			.add_method("RE",  &T::RE)
			.add_method("CD",  &T::CD)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "Interface", tag);
	}
//    LineWriter
	{
		typedef LineWriter<dim> T;
		string name = string("LineWriter").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_method("write_line",  &T::write_line)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LineWriter", tag);
	}
//    Granular Density Linker
	{
		typedef GranularDensityLinker<dim> T;
		typedef DependentUserData<number, dim> TBase;
		string name = string("GranularDensityLinker").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_volume_fraction", static_cast<void (T::*)(number)>(&T::set_volume_fraction))
			.add_method("set_volume_fraction", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_volume_fraction))
			.add_method("set_particle_density",  &T::set_particle_density)
			.add_method("set_fluid_density",  &T::set_fluid_density)
			.add_method("set_interface_volume_fraction",  &T::set_interface_volume_fraction)
			.add_method("set_model",  &T::set_model)
		
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GranularDensityLinker", tag);
	}
//    GranularViscosityLinker
	{
		typedef GranularViscosityLinker<dim> T;
		typedef DependentUserData<number, dim> TBase;
		string name = string("GranularViscosityLinker").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_velocity_gradient", &T::set_velocity_gradient)
			.add_method("set_volume_fraction", static_cast<void (T::*)(number)>(&T::set_volume_fraction))
			.add_method("set_volume_fraction", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_volume_fraction))
			.add_method("set_eins_viscosity", static_cast<void (T::*)(number)>(&T::set_eins_viscosity))
			.add_method("set_eins_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_eins_viscosity))
			.add_method("set_granular_model",  &T::set_granular_model)
			.add_method("set_particle_density",  &T::set_particle_density)
			.add_method("set_particle_kinematicVisc",  &T::set_particle_kinematicVisc)
			.add_method("set_fluid_density",  &T::set_fluid_density)
			.add_method("set_fluid_Visc",  &T::set_fluid_Visc)
			.add_method("set_particle_diameter",  &T::set_particle_diameter)
			.add_method("set_alpha_max",  &T::set_alpha_max)
			.add_method("set_alpha_min",  &T::set_alpha_min)
			//.add_method("set_packing_factor",  &T::set_packing_factor)
			.add_method("set_mix_density", static_cast<void (T::*)(number)>(&T::set_mix_density))
			.add_method("set_mix_density", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_mix_density))
			.add_method("set_particle_pressure", static_cast<void (T::*)(number)>(&T::set_particle_pressure))
			.add_method("set_particle_pressure", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_particle_pressure))
			.add_method("set_interface_volume_fraction",  &T::set_interface_volume_fraction)
			.add_method("set_deltaPs",  &T::set_deltaPs)
			.add_method("set_deltaI",  &T::set_deltaI)
			.add_method("set_FricMu_1",  &T::set_FricMu_1)
			.add_method("set_FricMu_2",  &T::set_FricMu_2)
			.add_method("set_I_0",  &T::set_I_0)
			.add_method("set_deltaGamma",  &T::set_deltaGamma)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GranularViscosityLinker", tag);
	}
//    GranularSourceLinker
	{
		typedef GranularSourceLinker<dim> T;
		typedef DependentUserData<MathVector<dim>, dim> TBase;
		string name = string("GranularSourceLinker").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_mix_density", static_cast<void (T::*)(number)>(&T::set_mix_density))
			.add_method("set_mix_density", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_mix_density))
			.add_method("set_particle_density",  &T::set_particle_density)
			.add_method("set_fluid_density",  &T::set_fluid_density)
			.add_method("set_gravity",  &T::set_gravity)
			.add_method("set_cons_gravity",  &T::set_cons_gravity)
			.add_method("set_phase_parameters", &T::set_phase_parameters)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GranularSourceLinker", tag);
	}
//    RelativeVelocityLinker
	{
		typedef RelativeVelocityLinker<dim> T;
		typedef DependentUserData<MathVector<dim>, dim> TBase;
		string name = string("RelativeVelocityLinker").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_volume_fraction", static_cast<void (T::*)(number)>(&T::set_volume_fraction))
			.add_method("set_volume_fraction", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_volume_fraction))
			.add_method("set_volume_grad", &T::set_volume_grad)
			.add_method("set_mix_density", static_cast<void (T::*)(number)>(&T::set_mix_density))
			.add_method("set_mix_density", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_mix_density))
			.add_method("set_mix_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_mix_kinematic_viscosity))
			.add_method("set_mix_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_mix_kinematic_viscosity))
			.add_method("set_einstein_visc", static_cast<void (T::*)(number)>(&T::set_einstein_visc))
			.add_method("set_einstein_visc", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_einstein_visc))
			.add_method("set_ps_grad", &T::set_ps_grad)
			.add_method("set_pressure_grad", &T::set_pressure_grad)
			.add_method("set_derivative_mask",  &T::set_derivative_mask)
			.add_method("set_particle_density",  &T::set_particle_density)
			.add_method("set_fluid_density",  &T::set_fluid_density)
			.add_method("set_particle_diameter",  &T::set_particle_diameter)
			.add_method("set_alpha_max",  &T::set_alpha_max)
			.add_method("set_fluid_viscosity",  &T::set_fluid_viscosity)
			.add_method("set_gravity",  &T::set_gravity)
			.add_method("set_vol_limit",  &T::set_vol_limit)
			.add_method("set_rel_vel",  &T::set_rel_vel)
			.add_method("set_dragCoeff",  &T::set_dragCoeff)
			.add_method("activate_relative_vel",  &T::activate_relative_vel)
			.add_method("set_phase_parameters", &T::set_phase_parameters)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "RelativeVelocityLinker", tag);
	}
//    GranularDiffusionLinker
	{
		typedef GranularDiffusionLinker<dim> T;
		typedef DependentUserData<MathMatrix<dim,dim>, dim> TBase;
		string name = string("GranularDiffusionLinker").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_gamma", static_cast<void (T::*)(number)>(&T::set_gamma))
			.add_method("set_gamma", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_gamma))
			.add_method("set_velocity_gradient", &T::set_velocity_gradient)
			.add_method("set_mix_viscosity", static_cast<void (T::*)(number)>(&T::set_mix_viscosity))
			.add_method("set_mix_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_mix_viscosity))
			.add_method("set_diff_factor",  &T::set_diff_factor)
			.add_method("const_kinematic_visc",  &T::const_kinematic_visc)
			.add_method("set_phase_parameters", &T::set_phase_parameters)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GranularDiffusionLinker", tag);
	}
//    ScalePressureLinker
	{
		typedef ScaleLinker<dim> T;
		typedef DependentUserData<number, dim> TBase;
		string name = string("ScaleLinker").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_import_1", static_cast<void (T::*)(number)>(&T::set_import_1))
			.add_method("set_import_1", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_import_1))
			.add_method("set_import_2", static_cast<void (T::*)(number)>(&T::set_import_2))
			.add_method("set_import_2", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_import_2))
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ScaleLinker", tag);
	}

//    check whether those classes have already been registered
   {
	   string name = string("INavierStokesUpwind").append(suffix);
	   if(reg.get_class(name))
		   return;
   }

/////////////////////////////////////////////////////////////////////////////
// Upwind
/////////////////////////////////////////////////////////////////////////////

//	INavierStokesUpwind
	{
		typedef INavierStokesUpwind<dim> T;
		string name = string("INavierStokesUpwind").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "INavierStokesUpwind", tag);
	}

//	NavierStokesNoUpwind
	{
		typedef NavierStokesNoUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesNoUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesNoUpwind", tag);
	}

//	NavierStokesFullUpwind
	{
		typedef NavierStokesFullUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesFullUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFullUpwind", tag);
	}

//	NavierStokesSkewedUpwind
	{
		typedef NavierStokesSkewedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesSkewedUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesSkewedUpwind", tag);
	}

//	NavierStokesLinearProfileSkewedUpwind
	{
		typedef NavierStokesLinearProfileSkewedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesLinearProfileSkewedUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesLinearProfileSkewedUpwind", tag);
	}

//	NavierStokesPositiveUpwind
	{
		typedef NavierStokesPositiveUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesPositiveUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesPositiveUpwind", tag);
	}

//	NavierStokesRegularUpwind
	{
		typedef NavierStokesRegularUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesRegularUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesRegularUpwind", tag);
	}

//	NavierStokesWeightedUpwind
/*	{
		typedef NavierStokesWeightedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesWeightedUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
		.template add_constructor<void (*)(number)>("Weight factor")
		.add_method("set_weight", &T::set_weight, "", "")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesWeightedUpwind", tag);
	}*/
}

}; // end Functionality
} // end namespace NavierStokes

/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_NavierStokes(Registry* reg, string grp)
{
	grp.append("SpatialDisc/NavierStokes/");
	typedef NavierStokes::Functionality Functionality;

	try{
		RegisterDimension2d3dDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
		Init___CompressibleNavierStokes(reg, grp);
		Init___IncompressibleNavierStokes(reg, grp);
//		RegisterAlgebraDependent<Functionality>(*reg,grp);
		//RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}


}// namespace ug
