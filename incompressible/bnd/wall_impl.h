/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__NAVIER_STOKES__BND__WALL_IMPL__
#define __H__UG__NAVIER_STOKES__BND__WALL_IMPL__

#include "wall.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace NavierStokes{

template <typename TDomain, typename TAlgebra>
NavierStokesWall<TDomain,TAlgebra>::
NavierStokesWall(SmartPtr< IncompressibleNavierStokesBase<TDomain> > spMaster)
	: m_spDirichletConstraint(new DirichletBoundary<TDomain,TAlgebra>)
{
	m_vFctName = spMaster->symb_fcts();

	/*if(m_vFctName.size() != TDomain::dim + 1)
		UG_THROW("NavierStokesWall::set_functions: This Boundary "
				"Condition works on exactly dim+1 (velocity+pressure) "
				"components, but "<<m_vFctName.size()<<"components given.");*/
}

template <typename TDomain, typename TAlgebra>
void NavierStokesWall<TDomain,TAlgebra>::
add(const char* subsetsBND)
{
	if(m_vFctName.empty())
		UG_THROW("NavierStokesWall::add: Symbolic names for"
				" velocity and pressure not set. Please set them first.");

	for(int i = 0; i < TDomain::dim; ++i)
	{
		m_spDirichletConstraint->add(0.0, m_vFctName[i].c_str(), subsetsBND);
	}
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__BND__WALL_IMPL__ */
