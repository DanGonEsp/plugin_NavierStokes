/*
 * Copyright (c) 2013-2017:  G-CSC, Goethe University Frankfurt
 * Author: Jonas Simon
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


#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

#include "wall_sliding_fv1.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace NavierStokes{

/**
 * The add method for the boundary subsets:
 */
template<typename TDomain>
void
NavierStokesWSBCFV1<TDomain>::
add
(
	const char* subsets
)
{
	m_vScheduledBndSubSets.push_back(subsets);
}

template<typename TDomain>
void NavierStokesWSBCFV1<TDomain>::
set_sliding_factor(number val)
{
	m_imSlidingFactor.set_data(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void NavierStokesWSBCFV1<TDomain>::
set_sliding_factor(const char* fctName)
{
	m_imSlidingFactor.set_data(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif

template<typename TDomain>
void NavierStokesWSBCFV1<TDomain>::
set_sliding_limit(number val)
{
	m_imSlidingLimit.set_data(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void NavierStokesWSBCFV1<TDomain>::
set_sliding_limit(const char* fctName)
{
	m_imSlidingLimit.set_data(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif

template<typename TDomain>
void NavierStokesWSBCFV1<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	if(bNonRegularGrid)
		UG_THROW("NavierStokesWSBCFV1: only regular grid implemented.");

//	check number
	if(vLfeID.size() != dim+1)
		UG_THROW("NavierStokesWSBCFV1: Needs exactly "<<dim+1<<" functions.");

//	check that Lagrange 1st order
	for(size_t i = 0; i < vLfeID.size(); ++i)
		if(vLfeID[i] != LFEID(LFEID::LAGRANGE, dim, 1))
			UG_THROW("NavierStokesWSBCFV1: only first order Lagrange supported.");
}

/**
 * Prepares the element loop for a given element type: computes the FV-geo, ...
 * Note that there are separate loops for every type of the grid elements.
 */
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesWSBCFV1<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	register subsetIndex at Geometry
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	Only first order implementation
	if(!(TFVGeom::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

//	check if kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("NavierStokesWSBCFV1::prep_elem_loop:"
						" Kinematic Viscosity has not been set, but is required.\n");

//	check if Density has been set
	if(!m_imDensity.data_given())
		UG_THROW("NavierStokesWSBCFV1::prep_elem_loop:"
						" Density has not been set, but is required.\n");

//  check if sliding factor has been set
	if(!m_imSlidingFactor.data_given())
		UG_THROW("NavierStokesWSBCFV1::prep_elem_loop:"
						" Sliding Factor has not been set, but is required.\n");

//  check if sliding limit has been set
	if(!m_imSlidingLimit.data_given())
		UG_THROW("NavierStokesWSBCFV1::prep_elem_loop:"
						" Sliding Limit has not been set. It can be set to 0 for deactivation.\n");

//	extract indices of boundary
	extract_scheduled_data();

//	request the subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
		geo.add_boundary_subset(*subsetIter);
}

/**
 * Finalizes the element loop for a given element type.
 */
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesWSBCFV1<TDomain>::
fsh_elem_loop()
{
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	remove the bnd subsets
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
		geo.remove_boundary_subset(*subsetIter);
}


/**
 * General initializations of a given grid element for the assembling.
 */
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesWSBCFV1<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("NavierStokesWSBCFV1::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	find and set the local and the global positions of the IPs for imports
	typedef typename TFVGeom::BF BF;
	typename std::vector<int>::const_iterator subsetIter;

	m_vLocIP.clear(); m_vGloIP.clear();
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
		const int bndSubset = *subsetIter;
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);
		for(size_t i = 0; i < vBF.size(); ++i)
		{
			m_vLocIP.push_back(vBF[i].local_ip());
			m_vGloIP.push_back(vBF[i].global_ip());
		}
	}
	// REMARK: The loop above determines the ordering of the integration points:
	// The "outer ordering" corresponds to the ordering of the subsets in
	// m_vBndSubSetIndex, and "inside" of this ordering, the ip's are ordered
	// according to the order of the boundary faces in the FV geometry structure.

//	Set ips for imports
	m_imKinViscosity.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imKinViscosity.set_global_ips(&m_vGloIP[0], m_vGloIP.size());

	m_imDensity.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imDensity.set_global_ips(&m_vGloIP[0], m_vGloIP.size());

	m_imSlidingLimit.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imSlidingLimit.set_global_ips(&m_vGloIP[0], m_vGloIP.size());

	m_imSlidingFactor.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imSlidingFactor.set_global_ips(&m_vGloIP[0], m_vGloIP.size());
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesWSBCFV1<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	typedef typename TFVGeom::BF BF;

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	size_t ip = 0;
	for(subsetIter = m_vBndSubSetIndex.begin();
		subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;
		
	//	get the list of the ip's:
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);

	// 	loop the boundary faces
		typename std::vector<BF>::const_iterator bf;
		for(bf = vBF.begin(); bf != vBF.end(); ++bf, ++ip)
		{

			////////////////////////////////////////////////////
			// Add convective flux
			////////////////////////////////////////////////////

			MathVector<dim> normal;

            VecNormalize(normal, bf->normal());

			if (!m_spMaster->stokes ())
			{
            //  Interpolate velocity with standard lagrange interpolation
				MathVector<dim> StdVel(0.0);
				for(size_t sh = 0; sh < bf->num_sh(); ++sh)
					for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
						StdVel[d1] += u(d1, sh) * bf->shape(sh);
				
            //  Calculate tangential part of the velocity
				VecScaleAppend(StdVel, -VecDot(StdVel, normal), normal);

				number old_momentum_flux = VecDot (StdVel, bf->normal ()) * m_imDensity [ip];
			
			//	Add flux to local Jacobian
				for(size_t sh = 0; sh < bf->num_sh(); ++sh)
				{
					number t = old_momentum_flux * bf->shape(sh);
					for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
						J(d1, bf->node_id(), d1, sh) += t;
				}
			}
			
			////////////////////////////////////////////////////
			// Add sliding coefficient part
			////////////////////////////////////////////////////
			
			MathMatrix<dim,dim> parallelVel;
			MatSet(parallelVel, 0.0);
			
		//  Evaluate velocity parallel to boundary
			for(size_t sh = 0; sh < bf->num_sh(); ++sh) {
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				{
					parallelVel(d1,d1) += bf->shape(sh);
				}
			}

		//  Calculate u_{||} = (u - (u \cdot n) \cdot n)
			for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
					for(size_t sh = 0; sh < (size_t)dim; ++sh)
						parallelVel(d1,d2) -= bf->shape(sh)*normal[d1]*normal[d2];
			
		//  Add k \cdot u_{||} to local Jacobian
			for(int d1 = 0; d1 < dim; ++d1)
				for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
					for(size_t sh = 0; sh < bf->num_sh(); ++sh)
					{
						J(d1, bf->node_id(), d2, sh) += (-1.0) * parallelVel(d1,d2) * m_imSlidingFactor[ip];
					}

			////////////////////////////////////////////////////
			// Add sliding limit part
			////////////////////////////////////////////////////
		
			MathVector<dim> normalizedParallelVel;
			VecSet(normalizedParallelVel, 0.0);

			for(size_t sh = 0; sh < bf->num_sh(); ++sh) {
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				{
					normalizedParallelVel[d1] += bf->shape(sh);
				}
			}

			if(VecLength(normalizedParallelVel) != 0)
			{
				VecNormalize(normalizedParallelVel, normalizedParallelVel);
				for(size_t sh = 0; sh < bf->num_sh(); ++sh)
					for(int d1 = 0; d1 < dim; ++d1)
						J(d1, bf->node_id(), d1, sh) += (-1.0) * normalizedParallelVel[d1] * m_imSlidingLimit[ip];
			}
            
            //    pressure term
            for(size_t sh = 0; sh < bf->num_sh(); ++sh)
                for(int d1 = 0; d1 < dim; ++d1)
                {
                    J(d1, bf->node_id(), _P_, sh) += bf->shape(sh) * bf->normal()[d1];
                }
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesWSBCFV1<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

	// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	typedef typename TFVGeom::BF BF;
	
	// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	size_t ip = 0;
	for(subsetIter = m_vBndSubSetIndex.begin();
		subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;
		
	//	get the list of the ip's:
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);

	// 	loop the boundary faces
		typename std::vector<BF>::const_iterator bf;
		for(bf = vBF.begin(); bf != vBF.end(); ++bf, ++ip)
		{

			////////////////////////////////////////////////////
			// Add convective flux
			////////////////////////////////////////////////////

            MathVector<dim> normal;
                        
            VecNormalize(normal, bf->normal());

			if (!m_spMaster->stokes ())
			{
            //  Interpolate velocity with standard lagrange interpolation
				MathVector<dim> StdVel(0.0);
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
					for(size_t sh = 0; sh < bf->num_sh(); ++sh)
						StdVel[d1] += u(d1, sh) * bf->shape(sh);

            //  Calculate parallel part of velocity
				VecScaleAppend(StdVel, -VecDot(StdVel, normal), normal);

				number old_momentum_flux = VecDot (StdVel, bf->normal ()) * m_imDensity [ip];
			
			// Add the flux to the defect
				for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
					d(d1, bf->node_id()) += old_momentum_flux * StdVel[d1];
			}

			////////////////////////////////////////////////////
			// Add sliding koefficient part
			////////////////////////////////////////////////////
			
			MathVector<dim> parallelVel;
			for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			{
				parallelVel[d1] = 0.0;
				for(size_t sh = 0; sh < (size_t)dim; ++sh)
				{
					parallelVel[d1] += u(d1, sh) * bf->shape(sh);
				}
			}

			VecScaleAppend(parallelVel, -VecDot(parallelVel, normal), normal);

			for(int d1 = 0; d1 < dim; ++d1)
				d(d1, bf->node_id()) += (-1.0) * parallelVel[d1] * m_imSlidingFactor[ip];


			////////////////////////////////////////////////////
			// Add sliding limit part
			////////////////////////////////////////////////////
			
			if(VecLength(parallelVel) != 0)
			{
				VecNormalize(parallelVel, parallelVel);
				for(int d1 = 0; d1 < dim; ++d1)
					d(d1, bf->node_id()) += (1.0) * parallelVel[d1] * m_imSlidingLimit[ip];
			}
            
            
            // pressure term
            number pressure = 0.0;
            for(size_t sh = 0; sh < bf->num_sh(); ++sh)
                pressure += bf->shape(sh) * u(_P_, sh);
                
            for(int d1 = 0; d1 < dim; ++d1)
                d(d1, bf->node_id()) += pressure * bf->normal()[d1];
		}
	}
}

/**
 * converts the subset names where the BC is imposed to the corresponding subset
 * indices (i.e. m_vScheduledBndSubSets -> m_vBndSubSetIndex):
 */
template<typename TDomain>
void
NavierStokesWSBCFV1<TDomain>::
extract_scheduled_data()
{
//	clear all extracted data
	m_vBndSubSetIndex.clear();

//	loop all scheduled subsets
	for(size_t i = 0; i < m_vScheduledBndSubSets.size(); ++i)
	{
	//	create Subset Group
		SubsetGroup subsetGroup;

	//	convert strings
		try{
			subsetGroup = this->approx_space()->subset_grp_by_name(m_vScheduledBndSubSets[i].c_str());
		}UG_CATCH_THROW("'NavierStokesWSBCFV1:extract_scheduled_data':"
						" Subsets '" <<m_vScheduledBndSubSets[i].c_str() <<"' not"
						" all contained in ApproximationSpace.");

	//	get subsethandler
		const ISubsetHandler& rSH = *this->function_pattern()->subset_handler();

	// 	loop subsets
		for(size_t si = 0; si < subsetGroup.size(); ++si)
		{
		//	get subset index
			const int subsetIndex = subsetGroup[si];

		//	check that subsetIndex is valid
			if(subsetIndex < 0 || subsetIndex >= rSH.num_subsets())
			{
				UG_LOG("ERROR in 'NavierStokesWSBCFV1:extract_scheduled_data':"
						" Invalid subset Index " << subsetIndex <<
						". (Valid is 0, .. , " << rSH.num_subsets() <<").\n");
				return;
			}

		// save the index
			m_vBndSubSetIndex.push_back(subsetIndex);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesWSBCFV1<TDomain>::
NavierStokesWSBCFV1(SmartPtr< IncompressibleNavierStokesBase<TDomain> > spMaster)
: IElemDisc<TDomain>(spMaster->symb_fcts(), spMaster->symb_subsets()), m_spMaster (spMaster)
{
//	check number of functions
	if(this->num_fct() != dim+1)
		UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");

//	yet no boundary subsets
	m_vBndSubSetIndex.clear ();
    
    m_imKinViscosity.set_comp_lin_defect(false);
    m_imDensity.set_comp_lin_defect(false);
    m_imSlidingFactor.set_comp_lin_defect(false);
    m_imSlidingLimit.set_comp_lin_defect(false);

//	register imports
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensity);
	this->register_import(m_imSlidingFactor);
	this->register_import(m_imSlidingLimit);

//	initialize the imports from the master discretization
	m_imKinViscosity.set_data(spMaster->kinematic_viscosity ());
	m_imDensity.set_data(spMaster->density ());

//	register assemble functions
	this->register_all_funcs(false);
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void NavierStokesWSBCFV1<Domain1d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesNoNormalStressOutflowFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_2
template<>
void NavierStokesWSBCFV1<Domain2d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesNoNormalStressOutflowFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_3
template<>
void NavierStokesWSBCFV1<Domain3d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, FV1Geometry<Prism, dim> >();
		register_func<Pyramid, FV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, FV1Geometry<Hexahedron, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesNoNormalStressOutflowFV1: Hanging Nodes not implemented.")
	}
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void
NavierStokesWSBCFV1<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	   	id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct(	id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	id, &T::template add_rhs_elem<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class NavierStokesWSBCFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NavierStokesWSBCFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesWSBCFV1<Domain3d>;
#endif

} // namespace NavierStokes
} // end namespace ug
