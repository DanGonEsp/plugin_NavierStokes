/*
 * Copyright (c) 2012-2013:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__NAVIER_STOKES__BND__INFLOW_BASE_IMPL_MULTIPHASE__
#define __H__UG__NAVIER_STOKES__BND__INFLOW_BASE_IMPL_MULTIPHASE__

#include "inflow_base_multiphase.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace NavierStokes{

template <typename TDomain, typename TAlgebra>
void NavierStokesInflowBaseMultiphase<TDomain,TAlgebra>::
add(const std::vector<number>& vVel, const std::vector<number>& vMassFlux, const char* subsetsBND)
{
	if((int)vVel.size() != dim || (int)vMassFlux.size() != dim)
		UG_THROW("NavierStokesInflow: Setting velocity vector of dimension "<<
		         vVel.size()<<" to a Discretization for world dim " << dim);

	add(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)), SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vMassFlux)), subsetsBND);
}


#ifdef UG_FOR_LUA
template <typename TDomain, typename TAlgebra>
void NavierStokesInflowBaseMultiphase<TDomain,TAlgebra>::
add(const char* fctName, const char* fctName2, const char* subsetsBND)
{
	add(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName), LuaUserDataFactory<MathVector<dim>, dim>::create(fctName2), subsetsBND);
}
#endif

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__BND__INFLOW_BASE_IMPL_MULTIPHASE__ */
