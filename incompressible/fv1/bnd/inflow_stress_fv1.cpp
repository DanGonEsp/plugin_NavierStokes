/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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

#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

#include "inflow_stress_fv1.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace NavierStokes{

template<typename TDomain>
void NavierStokesInflowStressFV1<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	if(bNonRegularGrid)
		UG_THROW("NavierStokesSymBCFV1: only regular grid implemented.");

//	check number
	if(vLfeID.size() != dim+1)
		UG_THROW("NavierStokesSymBCFV1: Needs exactly "<<dim+1<<" functions.");

//	check that Lagrange 1st order
	for(size_t i = 0; i < vLfeID.size(); ++i)
		if(vLfeID[i] != LFEID(LFEID::LAGRANGE, dim, 1))
			UG_THROW("NavierStokesSymBCFV1: only first order Lagrange supported.");
}
/**
 * converts the subset names where the BC is imposed to the corresponding subset
 * indices (i.e. m_vScheduledBndSubSets -> m_vBndSubSetIndex):
 */
template<typename TDomain>
void
NavierStokesInflowStressFV1<TDomain>::
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
		}UG_CATCH_THROW("'NavierStokesInflowStressFV1:extract_scheduled_data':"
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
				UG_LOG("ERROR in 'NavierStokesInflowStressFV1:extract_scheduled_data':"
						" Invalid subset Index " << subsetIndex <<
						". (Valid is 0, .. , " << rSH.num_subsets() <<").\n");
				return;
			}

		// save the index
			m_vBndSubSetIndex.push_back(subsetIndex);
		}
	}
}

/**
 * The add method for the boundary subsets:
 */
template<typename TDomain>
void
NavierStokesInflowStressFV1<TDomain>::
add
(
	const char* subsets // string with the ','-separated names of the subsets
)
{
	m_vScheduledBndSubSets.push_back(subsets);
}

/**
 * Prepares the element loop for a given element type: computes the FV-geo, ...
 * Note that there are separate loops for every type of the grid elements.
 */
template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokesInflowStressFV1<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	register subsetIndex at Geometry
	static TFVGeom<TElem, dim>& geo = GeomProvider<TFVGeom<TElem,dim> >::get();

// 	Only first order implementation
	if(!(TFVGeom<TElem, dim>::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

//	check if kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("NavierStokesInflowStressFV1::prep_elem_loop:"
						" Kinematic Viscosity has not been set, but is required.\n");

//	check if Density has been set
	if(!m_imDensity.data_given())
		UG_THROW("NavierStokesInflowStressFV1::prep_elem_loop:"
						" Density has not been set, but is required.\n");

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
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokesInflowStressFV1<TDomain>::
fsh_elem_loop()
{
	static TFVGeom<TElem, dim>& geo = GeomProvider<TFVGeom<TElem,dim> >::get();

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
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokesInflowStressFV1<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	static TFVGeom<TElem, dim>& geo = GeomProvider<TFVGeom<TElem,dim> >::get();
	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("NavierStokesInflowStressFV1::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	find and set the local and the global positions of the IPs for imports
	typedef typename TFVGeom<TElem, dim>::BF BF;
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

	m_imKinViscosity.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imKinViscosity.set_global_ips(&m_vGloIP[0], m_vGloIP.size());

	m_imDensity.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imDensity.set_global_ips(&m_vGloIP[0], m_vGloIP.size());
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokesInflowStressFV1<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// 	Only first order implementation
		UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

	// 	get finite volume geometry
		static const TFVGeom<TElem, dim>& geo = GeomProvider<TFVGeom<TElem,dim> >::get();

	//	check for source term to pass to the stabilization
//		const DataImport<MathVector<dim>, dim>* pSource = NULL;
//		if(m_imSourceSCVF.data_given())	pSource = &m_imSourceSCVF;

	//	check for solutions to pass to stabilization in time-dependent case
/*		const LocalVector *pSol = &u, *pOldSol = NULL;
		number dt = 0.0;
		if(this->is_time_dependent())
		{
		//	get and check current and old solution
			const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
			if(vLocSol->size() != 2)
				UG_THROW("NavierStokes::add_jac_A_elem: "
								" Stabilization needs exactly two time points.");

		//	remember local solutions
			pSol = &vLocSol->solution(0);
			pOldSol = &vLocSol->solution(1);
			dt = vLocSol->time(0) - vLocSol->time(1);
		}
*/
	//	interpolate velocity at ip with standard lagrange interpolation
		static const size_t numSCVF = TFVGeom<TElem, dim>::numSCVF;
		MathVector<dim> StdVel[numSCVF];
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
			const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);
			VecSet(StdVel[ip], 0.0);

			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				for(int d1 = 0; d1 < dim; ++d1)
					StdVel[ip][d1] += u(d1, sh) * scvf.shape(sh);
		}

	typedef typename TFVGeom<TElem, dim>::BF BF;

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
		for(bf = vBF.begin(); bf != vBF.end(); ++bf)
		{
			// 	loop shape functions
			for(size_t sh = 0; sh < bf->num_sh(); ++sh)
			{
				////////////////////////////////////////////////////
				// Add normal flux
				////////////////////////////////////////////////////

			// 	Compute flux derivative at IP
				/*const number flux_sh =  -1.0 * m_imKinViscosity[ip]  * m_imDensity[ip]
										* VecDot(bf->global_grad(sh), bf->normal());

			// 	Add flux derivative  to local matrix
				for(int d1 = 0; d1 < dim; ++d1)
				{
					J(d1, bf->node_id(), d1, sh) += flux_sh;
				}
				if(!m_spMaster->laplace())
				{
					for(int d1 = 0; d1 < dim; ++d1)
						for(int d2 = 0; d2 < dim; ++d2)
						{
							const number flux2_sh = -1.0 * m_imKinViscosity[ip]  * m_imDensity[ip]
													* bf->global_grad(sh)[d1] * bf->normal()[d2];
							J(d1, bf->node_id(), d2, sh) += flux2_sh;
						}
				}
				//	pressure term
				for(int d1 = 0; d1 < dim; ++d1)
				{
					J(d1, bf->node_id(), _P_, sh) += bf->shape(sh) * bf->normal()[d1];
				}*/
				// do nothing for continuity equation because un = 0
                const number D =   ElementDiameter<GridObject, TDomain>(*elem, *this->domain());
                number scale = 1.0 / (m_imDensity[ip] * (VecLength(StdVel[ip])/D + m_imKinViscosity[ip]/pow(D,2)));

                /*for(size_t sh = 0; sh < bf->num_sh(); ++sh) // loop shape functions
                {
                    for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
                        J(_P_, bf->node_id (), d2, sh) += bf->shape(sh) * bf->normal()[d2]; //* m_imDensity [ip];
                    
                    //J(_P_, bf->node_id (), _P_, sh) += scale*VecProd(bf->global_grad(sh) , bf->normal());
                }*/
			}// end of loop shape functions
			
		}
	}
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokesInflowStressFV1<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implemented
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = GeomProvider<TFVGeom<TElem,dim> >::get();
	typedef typename TFVGeom<TElem, dim>::BF BF;

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;

	size_t ip = 0;

	for(subsetIter = m_vBndSubSetIndex.begin();subsetIter != m_vBndSubSetIndex.end(); ++subsetIter){
		//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;

	//	get the list of the ip's:
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);

	// 	loop the boundary faces
		typename std::vector<BF>::const_iterator bf;
		for(bf = vBF.begin(); bf != vBF.end(); ++bf)
		{
		// 	1. Interpolate Functional Matrix of velocity at ip
			MathMatrix<dim, dim> gradVel;
			for(int d1 = 0; d1 < dim; ++d1)
				for(int d2 = 0; d2 <dim; ++d2)
				{
				//	sum up contributions of each shape
					gradVel(d1, d2) = 0.0;
					for(size_t sh = 0; sh < bf->num_sh(); ++sh)
						gradVel(d1, d2) += bf->global_grad(sh)[d2]
			                 * u(d1, sh);
				}
/*
		//	2. Compute flux
			MathVector<dim> diffFlux;
		//	Add (\nabla u) \cdot \vec{n}
			MatVecMult(diffFlux, gradVel, bf->normal());

		//	Add (\nabla u)^T \cdot \vec{n}
			if(!m_spMaster->laplace())
				TransposedMatVecMultAdd(diffFlux, gradVel, bf->normal());
    

		//	scale by viscosity
			VecScale(diffFlux, diffFlux, (-1.0) * m_imKinViscosity[ip] * m_imDensity[ip]);

			for(int d1 = 0; d1 < dim; ++d1)
			{
				d(d1, bf->node_id()) += diffFlux[d1];
			}
			// pressure term
			number pressure = 0.0;
			for(size_t sh = 0; sh < bf->num_sh(); ++sh)
				pressure += bf->shape(sh) * u(_P_, sh);
				
			for(int d1 = 0; d1 < dim; ++d1)
				d(d1, bf->node_id()) += pressure * bf->normal()[d1];*/
            
            
        // A. Compute Velocity at ip
            MathVector<dim> stdVel(0.0);
            MathVector<dim> PressureGrad(0.0);
            
            for(size_t sh = 0; sh < bf->num_sh(); ++sh)
            {
                for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
                {
                    stdVel[d1] += u(d1, sh) * bf->shape(sh);
                    PressureGrad[d1] += u(_P_, sh) * bf->global_grad(sh)[d1];
                }
            }
			
			// do nothing for continuity equation because un = 0
            const number D =   ElementDiameter<GridObject, TDomain>(*elem, *this->domain());
            number scale = 1.0 / (m_imDensity[ip] * (VecLength(stdVel)/D + m_imKinViscosity[ip]/pow(D,2)));
            d(_P_, bf->node_id()) -= VecDot (stdVel, bf->normal());// * m_imDensity[ip];
            //printf("Defect error = %f\n",d(_P_, bf->node_id()));
            //d(_P_, bf->node_id()) += scale*VecDot (PressureGrad, bf->normal());
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesInflowStressFV1<TDomain>::
NavierStokesInflowStressFV1(SmartPtr< IncompressibleNavierStokesBase<TDomain> > spMaster)
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

//	register imports
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensity);

//	initialize the imports from the master discretization
	m_imKinViscosity.set_data(spMaster->kinematic_viscosity ());
	m_imDensity.set_data(spMaster->density ());

//	register assemble functions
	this->register_all_funcs(false);
}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void
NavierStokesInflowStressFV1<TDomain>::
register_all_funcs(bool bHang)
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::DimElemList ElemList;
	// REMARK: Note that we register this boundary condition only
	// for the full-dimensional elements (DimElemList instead of AllElemList).

//	switch assemble functions
	if(!bHang) boost::mpl::for_each<ElemList>( RegisterFV1<FV1Geometry>(this) );
	else UG_THROW("Not implemented.");
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void
NavierStokesInflowStressFV1<TDomain>::
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
template class NavierStokesInflowStressFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NavierStokesInflowStressFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesInflowStressFV1<Domain3d>;
#endif

} // namespace NavierStokes
} // end namespace ug
