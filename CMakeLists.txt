# Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

################################################################################
# SAMPLE
#
# This script shows how easy it is to add a build script for your plugins.
# Simply specify the plugins name and a list of your source files.
################################################################################

# Name of your plugin and sources.
set(pluginName	NavierStokes)
set(SOURCES		sample_plugin.cpp)

set(SOURCES	navier_stokes_base.cpp
			upwind_interface.cpp
			upwind.cpp
			register_navier_stokes.cpp
			
			compressible/compressible_navier_stokes_base.cpp
			compressible/fv1/compressible_navier_stokes_fv1.cpp
			compressible/fv1/register_fv1.cpp
			compressible/compressible_navier_stokes_plugin.cpp


			incompressible/incompressible_navier_stokes_base.cpp	
			
			incompressible/fv1/navier_stokes_fv1.cpp
			incompressible/fe/navier_stokes_fe.cpp
			incompressible/fvcr/navier_stokes_fvcr.cpp
			incompressible/fv/navier_stokes_fv.cpp

			incompressible/bnd/no_normal_stress_outflow_base.cpp
			incompressible/fv1/bnd/no_normal_stress_outflow_fv1.cpp
			incompressible/fv/bnd/no_normal_stress_outflow_fv.cpp
			incompressible/fvcr/bnd/no_normal_stress_outflow_fvcr.cpp

			incompressible/fv1/bnd/wall_sliding_fv1.cpp
			incompressible/fvcr/bnd/symmetric_boundary_fvcr.cpp
			incompressible/fv1/bnd/symmetric_boundary_fv1.cpp

			incompressible/fvcr/cr_reorder.cpp

			incompressible/fv1/stabilization.cpp
			
			incompressible/fv1/register_fv1.cpp
			incompressible/fv/register_fv.cpp
			incompressible/fvcr/register_fvcr.cpp
			incompressible/fe/register_fe.cpp	
			
			incompressible/incompressible_navier_stokes_plugin.cpp)


################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 2.8.12...3.20)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4)
endif(buildEmbeddedPlugins)
