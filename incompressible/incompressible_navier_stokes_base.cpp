/*
 * Copyright (c) 2012-2017:  G-CSC, Goethe University Frankfurt
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

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

#include "incompressible_navier_stokes_base.h"
#include <vector>
#include <string>

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
IncompressibleNavierStokesBase<TDomain>::IncompressibleNavierStokesBase(const char* functions,
                                    const char* subsets)
: NavierStokesBase<TDomain>(functions, subsets),
  m_bPecletBlend(false),
  m_bStokes(false),
  m_bLaplace(false)
{
      m_exVelocity = make_sp(new DataExport<MathVector<dim>, dim>(functions));
      m_exVelocityGrad = make_sp(new DataExport<MathMatrix<dim, dim>, dim>(functions));
      m_exVelocity_ip = make_sp(new DataExport<MathVector<dim>, dim>(functions));
      m_exPressure = make_sp(new DataExport<number, dim>(functions));
      m_exPressureGrad = make_sp(new DataExport<MathVector<dim>, dim>(functions));

};

template<typename TDomain>
IncompressibleNavierStokesBase<TDomain>::IncompressibleNavierStokesBase(const std::vector<std::string>& vFct,
                                    const std::vector<std::string>& vSubset)
: NavierStokesBase<TDomain>(vFct, vSubset),
  m_bPecletBlend(false),
  m_bStokes(false),
  m_bLaplace(false)
{
  std::string functions;
  for(size_t i = 0; i < vFct.size(); ++i){
    if(i > 0) functions.append(",");
    functions.append(vFct[i]);
  }
      m_exVelocity = make_sp(new DataExport<MathVector<dim>, dim>(functions.c_str()));
      m_exVelocityGrad = make_sp(new DataExport<MathMatrix<dim, dim>, dim>(functions.c_str()));
      m_exVelocity_ip = make_sp(new DataExport<MathVector<dim>, dim>(functions.c_str()));
      m_exPressure = make_sp(new DataExport<number, dim>(functions.c_str()));
      m_exPressureGrad = make_sp(new DataExport<MathVector<dim>, dim>(functions.c_str()));

};


/////////// density

template<typename TDomain>
void IncompressibleNavierStokesBase<TDomain>::
set_density(number val)
{
	set_density(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void IncompressibleNavierStokesBase<TDomain>::
set_density(const char* fctName)
{
	set_density(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class IncompressibleNavierStokesBase<Domain1d>;
#endif
#ifdef UG_DIM_2
template class IncompressibleNavierStokesBase<Domain2d>;
#endif
#ifdef UG_DIM_3
template class IncompressibleNavierStokesBase<Domain3d>;
#endif

} // namespace NavierStokes
} // end namespace ug



