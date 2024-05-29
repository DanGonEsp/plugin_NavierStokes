/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__NAVIER_STOKES_TOOLS__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__NAVIER_STOKES_TOOLS__

#include <vector>
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "common/profiler/profiler.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
#include "lib_disc/quadrature/quadrature_provider.h"
#include "lib_disc/function_spaces/grid_function_global_user_data.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

namespace ug{

// Crouzeix-Raviart function is interpolated to Lagrange 1 function
template <typename TGridFunction>
void interpolateCRToLagrange(TGridFunction& uLagrange,TGridFunction& uCR){
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	///	world dimension
	static const int dim = domain_type::dim;

	// get grid
	grid_type& grid = *uLagrange.domain()->grid();

	// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// side type
	typedef typename elem_type::side side_type;

	//  volume attachment
	typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;
	aVertexNumber m_acVolume;
	ANumber m_aVolume;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type& posAcc = uLagrange.domain()->position_accessor();
	
	// scalar function or vector function
	int spaceDim = uLagrange.num_fct();
	if ((int)uLagrange.num_fct()==dim) spaceDim=dim;

	for (int i=0;i<spaceDim;i++){
		if (uLagrange.local_finite_element_id(i) != LFEID(LFEID::LAGRANGE, TGridFunction::dim, 1)){
			UG_THROW("First parameter must be of Lagrange 1 type.");
		}
		if (uCR.local_finite_element_id(i) != LFEID(LFEID::CROUZEIX_RAVIART, TGridFunction::dim, 1)){
			UG_THROW("Second parameter must be of CR type.");
		}
	}
/*	if (uLagrange.num_fct()==dim+1){
		if (uLagrange.local_finite_element_id(dim) != LFEID(LFEID::PIECEWISE_CONSTANT, TGridFunction::dim, 0)){
					UG_THROW("Parameter dim must be of piecewise type.");
		}
	}
	if (uCR.num_fct()==dim+1){
		if (uCR.local_finite_element_id(dim) != LFEID(LFEID::PIECEWISE_CONSTANT, TGridFunction::dim, 0)){
					UG_THROW("Parameter dim must be of piecewise constant 1 type.");
		}
	}*/

	DimFV1Geometry<dim> geo;

	grid.template attach_to<Vertex>(m_aVolume);
	m_acVolume.access(grid,m_aVolume);

	SetAttachmentValues(m_acVolume,grid.template begin<Vertex>(), grid.template end<Vertex>(), 0);

	PeriodicBoundaryManager* pbm = grid.periodic_boundary_manager();

	// coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	// create Multiindex
	std::vector<DoFIndex> multInd;

	static const size_t MaxNumSidesOfElem = 10;

	typedef MathVector<dim> MVD;
	std::vector<MVD> uValue(MaxNumSidesOfElem);
	
	// set lagrange function to zero
	for(int si = 0; si < uLagrange.num_subsets(); ++si){
		//	get iterators
		VertexIterator iter = uLagrange.template begin<Vertex>(si);
		VertexIterator iterEnd = uLagrange.template end<Vertex>(si);
		for(  ;iter !=iterEnd; ++iter){
			Vertex* vrt = *iter;
			if (pbm && pbm->is_slave(vrt)) continue;
			for (int d=0;d<spaceDim;d++){
				uLagrange.inner_dof_indices(vrt, d, multInd);
				DoFRef(uLagrange,multInd[0])=0;
			}
		}
	}

	for(int si = 0; si < uLagrange.num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = uLagrange.template begin<elem_type>(si);
		ElemIterator iterEnd = uLagrange.template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			// get Crouzeix-Raviart trial space
			const LocalShapeFunctionSet<dim>& crTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
			}

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			grid.associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			size_t nofsides = sides.size();
			for (size_t s=0;s < nofsides;s++)
			{
				for (int d=0;d<spaceDim;d++){
					//	get indices of function fct on vertex
					uCR.inner_dof_indices(sides[s], d, multInd);
					//	read value of index from vector
					uValue[s][d]=DoFRef(uCR,multInd[0]);
				}
			}

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), uLagrange.domain()->subset_handler().get());
			for(size_t i = 0; i < numVertices; ++i){
				const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(i);
				number scvVol = scv.volume();
				m_acVolume[vVrt[i]]+=scvVol;
				std::vector<number> vShape;
				crTrialSpace.shapes(vShape, scv.local_ip());
				MathVector<dim> localValue = 0;
				for (size_t s=0;s<nofsides;s++)
					for (int d=0;d<spaceDim;d++)
						localValue[d]+=vShape[s]*uValue[s][d];
				for (int d=0;d<spaceDim;d++){
					uLagrange.inner_dof_indices(vVrt[i], d, multInd);
					DoFRef(uLagrange,multInd[0])+=scvVol*localValue[d];
				}
			}
			if (uLagrange.num_fct()==dim+1){
				uLagrange.inner_dof_indices(elem,dim,multInd);
				DoFRef(uLagrange,multInd[0])=DoFRef(uCR,multInd[0]);
			}
		}
	}
	// finish computation by averaging
	for(int si = 0; si < uLagrange.num_subsets(); ++si){
		//	get iterators
		VertexIterator iter = uLagrange.template begin<Vertex>(si);
		VertexIterator iterEnd = uLagrange.template end<Vertex>(si);
		for(  ;iter !=iterEnd; ++iter){
			Vertex* vrt = *iter;
			if (pbm && pbm->is_slave(vrt)) continue;
			for (int d=0;d<spaceDim;d++){
				uLagrange.inner_dof_indices(vrt, d, multInd);
				DoFRef(uLagrange,multInd[0])/=m_acVolume[vrt];
			}
		}
	}
}

// compute vorticity
// for velocity field (u,v,w) the vorticity is \partial_x v - \partial_y u
template <typename TGridFunction>
void vorticityFVCR(TGridFunction& vort,TGridFunction& u)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	///	world dimension
	static const int dim = domain_type::dim;

	// get grid
	grid_type& grid = *u.domain()->grid();

	// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

    /// side type
	typedef typename elem_type::side side_type;

	//  volume attachment
	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	aSideNumber m_acVolume;
	ANumber m_aVolume;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type& posAcc = u.domain()->position_accessor();

	DimCRFVGeometry<dim> geo;

	grid.template attach_to<side_type>(m_aVolume);
	m_acVolume.access(grid,m_aVolume);

	SetAttachmentValues(m_acVolume,grid.template begin<side_type>(), grid.template end<side_type>(), 0);

	if (vort.local_finite_element_id(0) != LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
				UG_THROW("Component " << 0 << " in approximation space of parameter 1 must be of Crouzeix-Raviart type.");
	}

	for (int d=0;d<dim;d++){
		if (u.local_finite_element_id(d) != LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
			UG_THROW("Component " << d << " in approximation space of parameter 2 must be of Crouzeix-Raviart type.");
		}
	}

	PeriodicBoundaryManager* pbm = grid.periodic_boundary_manager();

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	create Multiindex
	std::vector<DoFIndex> multInd;

	for(int si = 0; si < u.num_subsets(); ++si)
	{
		if (si>0) continue;
		//	get iterators
		ElemIterator iter = u.template begin<elem_type>(si);
		ElemIterator iterEnd = u.template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
			}
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), u.domain()->subset_handler().get());

			size_t nofsides = geo.num_scv();

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			grid.associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			static const size_t MaxNumSidesOfElem = 18;

			typedef MathVector<dim> MVD;
			std::vector<MVD> uValue(MaxNumSidesOfElem);

			for (size_t s=0;s < nofsides;s++)
			{
				for (int d=0;d<dim;d++){
					//	get indices of function fct on vertex
					u.dof_indices(sides[s], d, multInd);
					//	read value of index from vector
					uValue[s][d]=DoFRef(u,multInd[0]);
				}
			}

			for (size_t s=0;s < nofsides;s++)
			{
				number localvort = 0;
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);

				for(size_t sh = 0 ; sh < nofsides; ++sh){
					localvort += uValue[sh][1] * scv.global_grad(sh)[0] - uValue[sh][0] * scv.global_grad(sh)[1];
				};
				
				number vol = scv.volume();

				localvort*=vol;

				vort.dof_indices(sides[s], 0, multInd);
				DoFRef(vort,multInd[0])+=localvort;
				m_acVolume[sides[s]] += vol;
			}
		}
	}
	// average vorticity
	number maxvort = 0;
	for(int si = 0; si < u.num_subsets(); ++si)
	{
		SideIterator sideIter = vort.template begin<side_type>(si);
		SideIterator sideIterEnd = vort.template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			//	get Elem
			side_type* elem = *sideIter;
			// if periodic slave continue
			if (pbm && pbm->is_slave(elem)) continue;
			vort.dof_indices(elem, 0, multInd);
			DoFRef(vort,multInd[0])/=m_acVolume[elem];
			//UG_LOG("[" << 0.5*(posAcc[elem->vertex(0)][0] + posAcc[elem->vertex(0)][0]) << "," << 
			//	   0.5*(posAcc[elem->vertex(0)][1] + posAcc[elem->vertex(0)][1]) << "]" << " " << DoFRef(vort,multInd[0]) << "\n");
			if (DoFRef(vort,multInd[0])*DoFRef(vort,multInd[0])>maxvort*maxvort){
				maxvort = DoFRef(vort,multInd[0]);
			}
		}
	}
	grid.template detach_from<side_type>(m_aVolume);
}

template <typename TGridFunction>
void vorticityFV1(TGridFunction& vort,TGridFunction& u)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	///	world dimension
	static const int dim = domain_type::dim;

	// get grid
	grid_type& grid = *u.domain()->grid();

	// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	//  volume attachment
	typedef PeriodicAttachmentAccessor<Vertex,ANumber > aSideNumber;
	aSideNumber m_acVolume;
	ANumber m_aVolume;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// side iterator
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexConstIterator;

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type& posAcc = u.domain()->position_accessor();

	DimFV1Geometry<dim> geo;

	grid.template attach_to<Vertex>(m_aVolume);
	m_acVolume.access(grid,m_aVolume);

	SetAttachmentValues(m_acVolume,grid.template begin<Vertex>(), grid.template end<Vertex>(), 0);

	if (vort.local_finite_element_id(0) != LFEID(LFEID::LAGRANGE, dim, 1)){
				UG_THROW("Component " << 0 << " in approximation space of parameter 1 must be of Crouzeix-Raviart type.");
	}

	for (int d=0;d<dim;d++){
		if (u.local_finite_element_id(d) != LFEID(LFEID::LAGRANGE, dim, 1)){
			UG_THROW("Component " << d << " in approximation space of parameter 2 must be of Crouzeix-Raviart type.");
		}
	}

	PeriodicBoundaryManager* pbm = grid.periodic_boundary_manager();

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	create Multiindex
	std::vector<DoFIndex> multInd;

	for(int si = 0; si < u.num_subsets(); ++si)
	{
		if (si>0) continue;
		//	get iterators
		ElemIterator iter = u.template begin<elem_type>(si);
		ElemIterator iterEnd = u.template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				// UG_LOG("co_coord(" << i<< "+1,:)=" << coCoord[i] << "\n");
			}
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), u.domain()->subset_handler().get());

			static const size_t MaxNumSidesOfElem = 18;

			typedef MathVector<dim> MVD;
			std::vector<MVD> uValue(MaxNumSidesOfElem);

			for (size_t co=0;co < numVertices;co++)
			{
				for (int d=0;d<dim;d++){
					//	get indices of function fct on vertex
					u.dof_indices(elem->vertex(co), d, multInd);
					//	read value of index from vector
					uValue[co][d]=DoFRef(u,multInd[0]);
				}
			}

			for (size_t co=0;co < numVertices;co++)
			{
				number localvort = 0;
				const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(co);

				for(size_t sh = 0 ; sh < scv.num_sh(); ++sh){
					localvort += uValue[sh][1] * scv.global_grad(sh)[0] - uValue[sh][0] * scv.global_grad(sh)[1];
				};

				number vol = scv.volume();

				localvort*=vol;

				vort.dof_indices(elem->vertex(co), 0, multInd);
				DoFRef(vort,multInd[0])+=localvort;
				m_acVolume[elem->vertex(co)] += vol;
			}
		}
	}
	// average vorticity
	number maxvort = 0;
	for(int si = 0; si < u.num_subsets(); ++si)
	{
		// average vorticity
		VertexConstIterator vertexIter = vort.template begin<Vertex>(si);
		VertexConstIterator vertexIterEnd = vort.template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			//	get Elem
			Vertex* vrt = *vertexIter;
			// if periodic slave continue
			if (pbm && pbm->is_slave(vrt)) continue;
			vort.dof_indices(vrt, 0, multInd);
			DoFRef(vort,multInd[0])/=m_acVolume[vrt];
			if (DoFRef(vort,multInd[0])*DoFRef(vort,multInd[0])>maxvort*maxvort){
				maxvort = DoFRef(vort,multInd[0]);
			}
		}
	}
	grid.template detach_from<Vertex>(m_aVolume);
}

template <typename TGridFunction>
void vorticity(TGridFunction& vort,TGridFunction& u){
	if (vort.local_finite_element_id(0) == LFEID(LFEID::LAGRANGE, TGridFunction::dim, 1)){
		vorticityFV1<TGridFunction>(vort,u);
	} else
		if (vort.local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, TGridFunction::dim, 1)){
			vorticityFVCR<TGridFunction>(vort,u);
		} else {
			UG_LOG("Function type " << vort.local_finite_element_id(0) << " not supported in vorticity computation.");
		}
}

template <typename TGridFunction>
void DrivenCavityEvalAtPoints(const std::vector<MathVector<2> >& vPos,
                              GlobalGridFunctionNumberData<TGridFunction>& GFEval,
                              const number vReferenceValue[])
{
#ifndef UG_DIM_3
	number maxdiff = 0, diffsum = 0;
	ug::Table<std::stringstream> table;
	table(0,0)<<"#"; table(0,1)<<"Position";table(0,2)<<"Measure";table(0,3)<<"Reference";table(0,4)<<"Difference";
	for(size_t i = 0; i < vPos.size(); ++i)
	{
		number val;
		const number ref = vReferenceValue[i];
		GFEval.evaluate_global(val, vPos[i]);
		const number localdiff = std::abs(ref-val);
		maxdiff = std::max(localdiff, maxdiff);
		diffsum += localdiff;

		table(i+1, 0) << std::setw(2) << i+1;
		table(i+1, 1) << std::fixed << std::setprecision(4) << vPos[i];
		table(i+1, 2) << std::fixed << std::setprecision(8) << val;
		table(i+1, 3) << std::fixed << std::setprecision(7) << ref;
		table(i+1, 4) << std::scientific << std::setprecision( 3 ) << localdiff;
	}
	std::cout << table;
	UG_LOG("\t     Max Diff: " << maxdiff << "\n")
	UG_LOG("\t Average Diff: " << diffsum/vPos.size() << "\n\n");
#else
	UG_LOG("DrivenCavityEvalAtPoints: Failed. Compile without dimension 3 to use this command.\n");
#endif
}

template <typename TGridFunction>
void DrivenCavityLinesEval(SmartPtr<TGridFunction> u, std::vector<std::string> vVelCmp, size_t Re)
{
	// data from Ghia paper, see also cast3m implementation
	const size_t numGhiaPoints = 17;

	const number vertGhiaPosX = 0.5;
	const number vertGhiaPosY[17]={0.0000,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5000,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1.0000};
	const number vertGhia_100[17]  ={0.,-0.03717,-0.04192,-0.04775,-0.06434,-0.10150,-0.15662,-0.2109,-0.20581,-0.13641,0.00332,0.23151,0.68717,0.73722,0.78871,0.84123,1.};
	const number vertGhia_400[17]  ={0.,-0.08186,-0.09266,-0.10338,-0.14612,-0.24299,-0.32726,-0.17119,-0.11477,0.02135,0.16256,0.29093,0.55892,0.61756,0.68439,0.75837,1.};
	const number vertGhia_1000[17] ={0.,-0.18109,-0.20196,-0.22220,-0.29730,-0.38289,-0.27805,-0.10648,-0.06080,0.05702,0.18719,0.33304,0.46604,0.51117,0.57492,0.65928,1.};
	const number vertGhia_3200[17] ={0.00000,-0.32407,-0.35344,-0.37827,-0.41933,-0.34323,-0.24427,-0.086636,-0.04272,0.07156,0.19791,0.34682,0.46101,0.46547,0.48296,0.53236,1.00000};
	const number vertGhia_5000[17] ={0.00000,-0.41165,-0.42901,-0.43643,-0.40435,-0.33050,-0.22855,-0.07404,-0.03039,0.08183,0.20087,0.33556,0.46036,0.45992,0.46120,0.48223,1.00000};
	const number vertGhia_7500[17] ={0.00000,-0.43154,-0.43590,-0.43025,-0.38324,-0.32393,-0.23176,-0.07503,-0.03800,0.08342,0.20591,0.34228,0.47167,0.47323,0.47048,0.47244,1.00000};
	const number vertGhia_10000[17]={0.00000,-0.42735,-0.42537,-0.41657,-0.38000,-0.32709,-0.23186,-0.07540,-0.03111,0.08344,0.20673,0.34635,0.47804,0.48070,0.47783,0.47221,1.00000};

	const number horizGhiaPosX[17]=   {0.0000,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5000,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1.0000};
	const number horizGhiaPosY = 0.5;
	const number horizGhia100[17]  ={0.00000,0.09233,0.10091,0.10890,0.12317,0.16077,0.17507,0.17527,0.05454,-0.24533,-0.22445,-0.16914,-0.10313,-0.08864,-0.07391,-0.05906,0.00000};
	const number horizGhia400[17]  ={0.00000,0.18360,0.19713,0.20920,0.22965,0.28124,0.30203,0.30174,0.05186,-0.38598,-0.44993,-0.3827,-0.22847,-0.19254,-0.15663,-0.12146,0.00000};
	const number horizGhia1000[17] ={0.00000,0.27485,0.29012,0.30353,0.32627,0.37095,0.33075,0.32235,0.02526,-0.31966,-0.42665,-0.51550,-0.39188,-0.33714,-0.27669,-0.21388,0.00000};
	const number horizGhia3200[17] ={0.00000,0.39560,0.40917,0.41906,0.42768,0.37119,0.29030,0.28188,0.00999,-0.31184,-0.37401,-0.44307,-0.54053,-0.52357,-0.47425,-0.39017,0.00000};
	const number horizGhia5000[17] ={0.00000,0.42447,0.43329,0.43648,0.42951,0.35368,0.28066,0.27280,0.00945,-0.30018,-0.36214,-0.41442,-0.52876,-0.55408,-0.55069,-0.49774,0.00000};
	const number horizGhia7500[17] ={0.00000,0.43979,0.44030,0.43564,0.41824,0.35060,0.28117,0.27348,0.00824,-0.30448,-0.36213,-0.41050,-0.48590,-0.52347,-0.55216,-0.53858,0.00000};
	const number horizGhia10000[17]={0.00000,0.43983,0.43733,0.43124,0.41487,0.35070,0.28003,0.27224,0.00831,-0.30719,-0.36737,-0.41496,-0.45863,-0.49099,-0.52987,-0.54302,0.00000};

	const number* vertGhia = NULL;
	const number* horizGhia = NULL;
	switch (Re)
	{
		case 100:	vertGhia = vertGhia_100;   horizGhia = horizGhia100; break;
		case 400:	vertGhia = vertGhia_400;   horizGhia = horizGhia400; break;
		case 1000:	vertGhia = vertGhia_1000;  horizGhia = horizGhia1000; break;
		case 3200:	vertGhia = vertGhia_3200;  horizGhia = horizGhia3200; break;
		case 5000:	vertGhia = vertGhia_5000;  horizGhia = horizGhia5000; break;
		case 7500:	vertGhia = vertGhia_7500;  horizGhia = horizGhia7500; break;
		case 10000:	vertGhia = vertGhia_10000; horizGhia = horizGhia10000; break;
		default: break;
	}

	// Botella reference data for Re=1000
	const number vertBotella_1000[17]={0.0000000 , -0.1812881 , -0.2023300 , -0.2228955 , -0.3004561 , -0.3885691 , -0.2803696 , -0.1081999 , -0.0620561 ,
	                                     0.0570178 , 0.1886747 , 0.3372212 , 0.4723329 , 0.5169277 , 0.5808359 , 0.6644227 , 1.0000000};
	const number horizBotella_1000[17]={0.0000000 , 0.2807056 , 0.2962703 , 0.3099097 , 0.3330442 , 0.3769189 , 0.3339924 , 0.3253592 , 0.0257995 ,
	                                     -0.3202137 , -0.4264545 , -0.5264392 , -0.4103754 , -0.3553213 , -0.2936869 , -0.2279225 , 0.0000000};

	// ug4 fvcr / linear upwind/linear pressure reference data for Re=3200 (computed on level 6)
	const number vertUG4_3200[17]={0, -3.565562e-01, -3.854442e-01,-4.083199e-01,-4.329174e-01, -3.455453e-01,-2.425736e-01,-8.123961e-02,-3.689807e-02,7.721821e-02,2.018428e-01,
	                                 3.482477e-01,4.618284e-01,4.654109e-01,4.811885e-01,5.280763e-01,1};
	const number horizUG4_3200[17]={0,3.961922e-01,4.113146e-01,4.224736e-01,4.327940e-01,3.771819e-01,2.964934e-01,2.881224e-01,1.425516e-02,-3.136940e-01,
	                                 -3.787908e-01,-4.436374e-01,-5.672544e-01,-5.610925e-01,-5.200556e-01,-4.390652e-01,0};

	// data from Erturk, Corke, Gkcl paper
	const size_t numErturkPoints = 23;
	const number vertErturkPosX = 0.5;
	const number vertErturkPosY[23]={1.00000,0.99000,0.98000,0.97000,0.96000,0.95000,0.94000,0.93000,0.92000,0.91000,0.90000,0.50000,0.20000,0.18000,0.16000,0.14000,0.12000,0.10000,0.08000,0.06000,0.04000,0.02000,0.00000};
	const number horizErturkPosX[23]={1.00000,0.98500,0.97000,0.95500,0.94000,0.92500,0.91000,0.89500,0.88000,0.86500,0.85000,0.50000,0.15000,0.13500,0.12000,0.10500,0.09000,0.07500,0.06000,0.04500,0.03000,0.01500,0.00000};
	const number horizErturkPosY = 0.5;
	const number horizErturk_1000[23] ={0.00000,-0.09730,-0.21730,-0.34000,-0.44170,-0.50520,-0.52630,-0.51320,-0.48030,-0.44070,-0.40280,0.02580,0.37560,0.37050,0.36050,0.34600,0.32730,0.30410,0.27460,0.23490,0.17920,0.10190,0.00000};
    const number horizErturk_2500[23] ={0.0000,-0.1675,-0.3725,-0.5192,-0.5603,-0.5268,-0.4741,-0.4321,-0.4042,-0.3843,-0.3671,0.0160,0.3918,0.4078,0.4187,0.4217,0.4142,0.3950,0.3649,0.3238,0.2633,0.1607,0.0000};
	const number horizErturk_5000[23] ={0.0000,-0.2441,-0.5019,-0.5700,-0.5139,-0.4595,-0.4318,-0.4147,-0.3982,-0.3806,-0.3624,0.0117,0.3699,0.3878,0.4070,0.4260,0.4403,0.4426,0.4258,0.3868,0.3263,0.2160,0.0000};
	const number horizErturk_7500[23] ={0.0000,-0.2991,-0.5550,-0.5434,-0.4748,-0.4443,-0.4283,-0.4118,-0.3938,-0.3755,-0.3574,0.0099,0.3616,0.3779,0.3950,0.4137,0.4337,0.4495,0.4494,0.4210,0.3608,0.2509,0.0000};
	const number horizErturk_10000[23]={0.0000,-0.3419,-0.5712,-0.5124,-0.4592,-0.4411,-0.4256,-0.4078,-0.3895,-0.3715,-0.3538,0.0088,0.3562,0.3722,0.3885,0.4056,0.4247,0.4449,0.4566,0.4409,0.3844,0.2756,0.0000};
	const number horizErturk_12500[23]={0.0000,-0.3762,-0.5694,-0.4899,-0.4534,-0.4388,-0.4221,-0.4040,-0.3859,-0.3682,-0.3508,0.0080,0.3519,0.3678,0.3840,0.4004,0.4180,0.4383,0.4563,0.4522,0.4018,0.2940,0.0000};
	const number horizErturk_15000[23]={0.0000,-0.4041,-0.5593,-0.4754,-0.4505,-0.4361,-0.4186,-0.4005,-0.3828,-0.3654,-0.3481,0.0074,0.3483,0.3641,0.3801,0.3964,0.4132,0.4323,0.4529,0.4580,0.4152,0.3083,0.0000};
	const number horizErturk_17500[23]={0.0000,-0.4269,-0.5460,-0.4664,-0.4482,-0.4331,-0.4153,-0.3975,-0.3800,-0.3627,-0.3457,0.0069,0.3452,0.3608,0.3767,0.3929,0.4093,0.4273,0.4484,0.4602,0.4254,0.3197,0.0000};
	const number horizErturk_20000[23]={0.0000,-0.4457,-0.5321,-0.4605,-0.4459,-0.4300,-0.4122,-0.3946,-0.3774,-0.3603,-0.3434,0.0065,0.3423,0.3579,0.3736,0.3897,0.4060,0.4232,0.4438,0.4601,0.4332,0.3290,0.0000};
	const number horizErturk_21000[23]={0.0000,-0.4522,-0.5266,-0.4588,-0.4449,-0.4287,-0.4110,-0.3936,-0.3764,-0.3593,-0.3425,0.0063,0.3413,0.3567,0.3725,0.3885,0.4048,0.4218,0.4420,0.4596,0.4357,0.3323,0.0000};
	const number vertErturk_1000[23] ={1.00000,0.84860,0.70650,0.59170,0.51020,0.45820,0.42760,0.41010,0.39930,0.39130,0.38380,-0.06200,-0.37560,-0.38690,-0.38540,-0.36900,-0.33810,-0.29600,-0.24720,-0.19510,-0.13920,-0.07570,0.00000};
    const number vertErturk_2500[23] ={1.0000,0.7704,0.5924,0.4971,0.4607,0.4506,0.4470,0.4424,0.4353,0.4256,0.4141,-0.0403,-0.3228,-0.3439,-0.3688,-0.3965,-0.4200,-0.4250,-0.3979,-0.3372,-0.2547,-0.1517,0.0000};
	const number vertErturk_5000[23] ={1.0000,0.6866,0.5159,0.4749,0.4739,0.4738,0.4683,0.4582,0.4452,0.4307,0.4155,-0.0319,-0.3100,-0.3285,-0.3467,-0.3652,-0.3876,-0.4168,-0.4419,-0.4272,-0.3480,-0.2223,0.0000};
	const number vertErturk_7500[23] ={1.0000,0.6300,0.4907,0.4817,0.4860,0.4824,0.4723,0.4585,0.4431,0.4275,0.4123,-0.0287,-0.3038,-0.3222,-0.3406,-0.3587,-0.3766,-0.3978,-0.4284,-0.4491,-0.3980,-0.2633,0.0000};
	const number vertErturk_10000[23]={1.0000,0.5891,0.4837,0.4891,0.4917,0.4843,0.4711,0.4556,0.4398,0.4243,0.4095,-0.0268,-0.2998,-0.3179,-0.3361,-0.3543,-0.3721,-0.3899,-0.4142,-0.4469,-0.4259,-0.2907,0.0000};
	const number vertErturk_12500[23]={1.0000,0.5587,0.4833,0.4941,0.4937,0.4833,0.4684,0.4523,0.4366,0.4216,0.4070,-0.0256,-0.2967,-0.3146,-0.3326,-0.3506,-0.3685,-0.3859,-0.4054,-0.4380,-0.4407,-0.3113,0.0000};
	const number vertErturk_15000[23]={1.0000,0.5358,0.4850,0.4969,0.4937,0.4811,0.4653,0.4492,0.4338,0.4190,0.4047,-0.0247,-0.2942,-0.3119,-0.3297,-0.3474,-0.3652,-0.3827,-0.4001,-0.4286,-0.4474,-0.3278,0.0000};
	const number vertErturk_17500[23]={1.0000,0.5183,0.4871,0.4982,0.4925,0.4784,0.4622,0.4463,0.4312,0.4166,0.4024,-0.0240,-0.2920,-0.3096,-0.3271,-0.3446,-0.3622,-0.3797,-0.3965,-0.4206,-0.4490,-0.3412,0.0000};
	const number vertErturk_20000[23]={1.0000,0.5048,0.4889,0.4985,0.4906,0.4754,0.4592,0.4436,0.4287,0.4142,0.4001,-0.0234,-0.2899,-0.3074,-0.3248,-0.3422,-0.3595,-0.3769,-0.3936,-0.4143,-0.4475,-0.3523,0.0000};
	const number vertErturk_21000[23]={1.0000,0.5003,0.4895,0.4983,0.4897,0.4742,0.4580,0.4425,0.4277,0.4132,0.3992,-0.0232,-0.2892,-0.3066,-0.3239,-0.3412,-0.3585,-0.3758,-0.3925,-0.4121,-0.4463,-0.3562,0.0000};


	const number* horizErturk = NULL;
	const number* vertErturk = NULL;
	switch (Re)
	{
		case 1000:	horizErturk = horizErturk_1000;  vertErturk = vertErturk_1000; break;
		case 2500:	horizErturk = horizErturk_2500;  vertErturk = vertErturk_2500; break;
		case 5000:	horizErturk = horizErturk_5000;  vertErturk = vertErturk_5000; break;
		case 7500:	horizErturk = horizErturk_7500;  vertErturk = vertErturk_7500; break;
		case 10000:	horizErturk = horizErturk_10000; vertErturk = vertErturk_10000; break;
		case 15000:	horizErturk = horizErturk_15000; vertErturk = vertErturk_15000; break;
		case 12500:	horizErturk = horizErturk_12500; vertErturk = vertErturk_12500; break;
		case 17500:	horizErturk = horizErturk_17500; vertErturk = vertErturk_17500; break;
		case 20000:	horizErturk = horizErturk_20000; vertErturk = vertErturk_20000; break;
		case 21000:	horizErturk = horizErturk_21000; vertErturk = vertErturk_21000; break;
		default: break;
	}

	// create Evaluation UserData
	GlobalGridFunctionNumberData<TGridFunction> uGFEval(u, vVelCmp[0].c_str());
	GlobalGridFunctionNumberData<TGridFunction> vGFEval(u, vVelCmp[1].c_str());
	std::vector<MathVector<2> > vPos(numGhiaPoints);

	if(vertGhia != NULL && horizGhia != NULL)
	{
		UG_LOG("  ------ Ghia, Re = " << Re << ": u values on a vertical line through x = 0.5  ------\n");
		for(size_t i = 0; i < numGhiaPoints; ++i) vPos[i] = MathVector<2>(vertGhiaPosX, vertGhiaPosY[i]);
		DrivenCavityEvalAtPoints<TGridFunction>(vPos, uGFEval, vertGhia);

		UG_LOG("  ------ Ghia, Re = " << Re << ": v values on a horizontal line through y = 0.5  ------\n");
		for(size_t i = 0; i < numGhiaPoints; ++i) vPos[i] = MathVector<2>(horizGhiaPosX[i], horizGhiaPosY);
		DrivenCavityEvalAtPoints<TGridFunction>(vPos, vGFEval, horizGhia);
	}

	if (Re==1000)
	{
		UG_LOG("  ------ Botella/Peyret, Re = " << Re << ": u values on a vertical line through x = 0.5  ------\n");
		for(size_t i = 0; i < numGhiaPoints; ++i) vPos[i] = MathVector<2>(vertGhiaPosX, vertGhiaPosY[i]);
		DrivenCavityEvalAtPoints<TGridFunction>(vPos, uGFEval, vertBotella_1000);

		UG_LOG("  ------ Botella/Peyret, Re = " << Re << ": v values on a horizontal line through y = 0.5  ------\n");
		for(size_t i = 0; i < numGhiaPoints; ++i) vPos[i] = MathVector<2>(horizGhiaPosX[i], horizGhiaPosY);
		DrivenCavityEvalAtPoints<TGridFunction>(vPos, vGFEval, horizBotella_1000);
	}

	if (Re==3200)
	{
		UG_LOG("  ------ ug4, Re = " << Re << ": u values on a vertical line through x = 0.5  ------\n");
		for(size_t i = 0; i < numGhiaPoints; ++i) vPos[i] = MathVector<2>(vertGhiaPosX, vertGhiaPosY[i]);
		DrivenCavityEvalAtPoints<TGridFunction>(vPos, uGFEval, vertUG4_3200);

		UG_LOG("  ------ ug4, Re = " << Re << ": v values on a horizontal line through y = 0.5  ------\n");
		for(size_t i = 0; i < numGhiaPoints; ++i) vPos[i] = MathVector<2>(horizGhiaPosX[i], horizGhiaPosY);
		DrivenCavityEvalAtPoints<TGridFunction>(vPos, vGFEval, horizUG4_3200);
	}

	if(horizErturk != NULL && vertErturk != NULL)
	{
		vPos.resize(numErturkPoints);
		UG_LOG("  ------ Erturk, Re = " << Re << ": u values on a vertical line through x = 0.5  ------\n");
		for(size_t i = 0; i < numErturkPoints; ++i) vPos[i] = MathVector<2>(vertErturkPosX, vertErturkPosY[i]);
		DrivenCavityEvalAtPoints<TGridFunction>(vPos, uGFEval, vertErturk);

		UG_LOG("  ------ Erturk, Re = " << Re << ": v values on a horizontal line through y = 0.5  ------\n");
		for(size_t i = 0; i < numErturkPoints; ++i) vPos[i] = MathVector<2>(horizErturkPosX[i], horizErturkPosY);
		DrivenCavityEvalAtPoints<TGridFunction>(vPos, vGFEval, horizErturk);
	}
}

/*
// Computation of maximum CFL-number over elements
// General CFL number is given as |u|*deltaT/h where u is velocity, deltaT time step length
// and h mesh width.
// Here following approximation for element local CFL-number is used:
// In element compute interpolated velocity v_bary in barycenter
// go over all "edges" between Crouzeix-Raviart nodes c_i and c_j (middle of edge/face)
// approximate velocity between nodes is 1/|c_i-c_j|*|(c_i-c_j)*v_bary|
// step length is |c_i-c_j| Therefore overall local cfl number is
// 1/|c_i-c_j|^2*|(c_i-c_j)*v_bary|
// Function computes maximum over all elements.
 */
template<typename TGridFunction>
number cflNumber(TGridFunction& u,number deltaT){
	/// domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	world dimension
	static const int dim = domain_type::dim;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// side type
	typedef typename elem_type::side side_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	// cfl number
	number maxCfl=0;
	//	get domain
	domain_type& domain = *u.domain().get();

	DimCRFVGeometry<dim> geo;

	//	create Multiindex
	std::vector<DoFIndex > multInd;

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = u.template begin<elem_type>(si);
		ElemIterator iterEnd = u.template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			MathVector<dim> bary,localBary;
			bary = 0;

			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				bary += coCoord[i];
			};
			bary /= numVertices;

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			//	memory for shapes
			std::vector<number> vShape;

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			rTrialSpace.shapes(vShape, localBary);

			size_t nofsides = geo.num_scv();

			MathVector<dim> baryV;
			baryV = 0;
			for (size_t s=0;s<nofsides;s++){
				MathVector<dim> localbaryV;
				for (int d=0;d<dim;d++){
					u.dof_indices(sides[s], d, multInd);
					localbaryV[d]=DoFRef(u,multInd[0]);
				}
				localbaryV *= vShape[s];
				baryV += localbaryV;
			}
			for (size_t i=0;i<nofsides;i++){
				const typename DimCRFVGeometry<dim>::SCV& scvi = geo.scv(i);
				MathVector<dim> iCoord = scvi.global_ip();
				for (size_t j=i+1;j<nofsides;j++){
					const typename DimCRFVGeometry<dim>::SCV& scvj = geo.scv(j);
					MathVector<dim> jCoord = scvj.global_ip();
					MathVector<dim> subVec;
					VecSubtract(subVec,iCoord,jCoord);
					number localCfl=deltaT*(number)1.0/VecTwoNormSq(subVec)*std::abs(subVec*baryV);
					if (localCfl>maxCfl) maxCfl=localCfl;
				}
			}
		}
	}
	UG_LOG("Max CFL number is " << maxCfl << "\n");
    return maxCfl;
}

// Compute 1/|\Omega| \int_{\Omega} \hat{u_i} \hat{u_i} dx , where \Omega is the computational domain.
// Elementwise approximation of integral by interpolation to barycenter.
template<typename TGridFunction>
number kineticEnergy(TGridFunction& u){
	// total kinetic energy
	number totalE=0;
	// total volume
	number totalVol=0;

	/// domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	world dimension
	static const int dim = domain_type::dim;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// side type
	typedef typename elem_type::side side_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	//	get domain of grid function
	domain_type& domain = *u.domain().get();
	DimCRFVGeometry<dim> geo;

	 std::vector<DoFIndex> multInd;

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	// assemble deformation tensor fluxes
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = u.template begin<elem_type>(si);
		ElemIterator iterEnd = u.template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			MathVector<dim> bary,localBary;
			bary = 0;

			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				bary += coCoord[i];
			};
			bary /= numVertices;

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			//	memory for shapes
			std::vector<number> vShape;

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			rTrialSpace.shapes(vShape, localBary);

			size_t nofsides = geo.num_scv();

			MathVector<dim> value;
			value = 0;
			number elementVolume = 0;
			for (size_t s=0;s<nofsides;s++){
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
				MathVector<dim> localValue;
				for (int d=0;d<dim;d++){
					u.dof_indices(sides[s], d, multInd);
					localValue[d]=DoFRef(u,multInd[0]);
				}
				localValue *= vShape[s];
				value += localValue;
				elementVolume += scv.volume();
			}
			for (int d=0;d<dim;d++){
				totalE += elementVolume*value[d]*value[d];
			}

			totalVol+=elementVolume;
		}
	}
	// average
	totalE/=(number)totalVol;
	UG_LOG("Total kinetic energy in domain is " << totalE << "\n");
	return totalE;
}

// clear file content
/* void clearFile(std::string filename){
	std::fstream file(filename.c_str(), std::fstream::out | std::fstream::trunc);
}
*/
// write numbers into file
/* void writeNumbers(std::string filename,const size_t step,const number t,const number data){
	std::fstream file(filename.c_str(), std::fstream::out | std::fstream::app);
	file << "t(" << step << ")=" << t << ";d(" << step << ")=" << data << ";" << std::endl;
}
*/


template <typename TGridFunction>
std::vector<number> DragLift(SmartPtr<TGridFunction> spGridFct,
                              const char* vCmp,
                              const char* BndSubsets, const char* InnerSubsets,
                              number kinVisco, number density,
                              int quadOrder)
{
	static const int dim = TGridFunction::dim;
	static const int WorldDim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;


//	read subsets
	SubsetGroup innerSSGrp(spGridFct->domain()->subset_handler());
	if(InnerSubsets != NULL){
		innerSSGrp.add(TokenizeString(InnerSubsets));
		if(!SameDimensionsInAllSubsets(innerSSGrp))
			UG_THROW("DragDrift: Subsets '"<<InnerSubsets<<"' do not have same dimension."
					 "Can not integrate on subsets of different dimensions.");
	}
	else{
		innerSSGrp.add_all();
		RemoveLowerDimSubsets(innerSSGrp);
	}

//	read subsets
	SubsetGroup bndSSGrp(spGridFct->domain()->subset_handler());
	if(BndSubsets != NULL)
		bndSSGrp.add(TokenizeString(BndSubsets));
	else
		UG_THROW("DragDrift: No boundary subsets passed.");

//	get function group
	const FunctionGroup vFctID = spGridFct->fct_grp_by_name(vCmp);
	std::vector<LFEID> vLFEID;
	for(size_t fct = 0; fct < vFctID.size(); ++fct){
		vLFEID.push_back(spGridFct->lfeid(vFctID[fct]));
	}

//	reset the result
	number int_lift = 0;
	number int_drag = 0;

//	loop subsets
	for(size_t i = 0; i < innerSSGrp.size(); ++i)
	{
	//	get subset index
		const int si = innerSSGrp[i];

	//	skip empty subset
		if(innerSSGrp.dim(i) == DIM_SUBSET_EMPTY_GRID) continue;

	//	check dimension
		if(innerSSGrp.dim(i) != dim)
			UG_THROW("DragDrift: Dimension of inner subset is "<<
					 innerSSGrp.dim(i)<<", but only World Dimension "<<dim<<
					 " subsets can be used for inner subsets.");

	//	note: this iterator is for the base elements, e.g. Face and not
	//			for the special type, e.g. Triangle, Quadrilateral
		const_iterator iterBegin = spGridFct->template begin<grid_base_object>(si);
		const_iterator iterEnd = spGridFct->template end<grid_base_object>(si);
		const_iterator iter = iterBegin;

		typename domain_traits<TGridFunction::dim>::position_accessor_type& aaPos
			= spGridFct->domain()->position_accessor();
		const ISubsetHandler* ish = spGridFct->domain()->subset_handler().get();
		Grid& grid = *spGridFct->domain()->grid();

	//	this is the base element type (e.g. Face). This is the type when the
	//	iterators above are dereferenciated.
		typedef typename domain_traits<dim>::element_type Element;
		typedef typename domain_traits<dim>::side_type Side;

	//	vector of corner coordinates of element corners (to be filled for each elem)
		std::vector<MathVector<WorldDim> > vCorner;
		std::vector<int> vSubsetIndex;

	// 	iterate over all elements
		for(; iter != iterEnd; ++iter)
		{
		//	get element
			Element* pElem = *iter;

		//	get all corner coordinates
			CollectCornerCoordinates(vCorner, *pElem, aaPos, true);

		//	get reference object id
			const ReferenceObjectID elemRoid = pElem->reference_object_id();

		//	get sides
			typename Grid::traits<Side>::secure_container vSide;
			grid.associated_elements_sorted(vSide, pElem);
			vSubsetIndex.resize(vSide.size());
			for(size_t i = 0; i < vSide.size(); ++i)
				vSubsetIndex[i] = ish->get_subset_index(vSide[i]);

			DimReferenceMapping<dim, WorldDim>& rMapping
				= ReferenceMappingProvider::get<dim, WorldDim>(elemRoid, vCorner);

			const DimReferenceElement<dim>& rRefElem
				= ReferenceElementProvider::get<dim>(elemRoid);

		//	get element values
			std::vector<DoFIndex> vInd;
			std::vector<std::vector<number> > vvValue(vFctID.size());
			for(size_t fct = 0; fct < vvValue.size(); ++fct){
				spGridFct->dof_indices(pElem, vFctID[fct], vInd);
				vvValue[fct].resize(vInd.size());
				for(size_t sh = 0; sh < vInd.size(); ++sh)
					vvValue[fct][sh] = DoFRef(*spGridFct, vInd[sh]);
			}
			const static int _P_ = dim;

		//	loop sub elements
			for(size_t side = 0; side < vSide.size(); ++side)
			{
			//	check if side used
				if(!bndSSGrp.contains(vSubsetIndex[side])) continue;

			//	get side
				Side* pSide = vSide[side];

				std::vector<MathVector<WorldDim> > vSideCorner(rRefElem.num(dim-1, side, 0));
				std::vector<MathVector<dim> > vLocalSideCorner(rRefElem.num(dim-1, side, 0));
				for(size_t co = 0; co < vSideCorner.size(); ++co){
					vSideCorner[co] = vCorner[rRefElem.id(dim-1, side, 0, co)];
					vLocalSideCorner[co] = rRefElem.corner(rRefElem.id(dim-1, side, 0, co));
				}

			//	side quad rule
				const ReferenceObjectID sideRoid = pSide->reference_object_id();
				const QuadratureRule<dim-1>& rSideQuadRule
						= QuadratureRuleProvider<dim-1>::get(sideRoid, quadOrder);

			// 	normal
				MathVector<WorldDim> Normal;
				ElementNormal<WorldDim>(sideRoid, Normal, &vSideCorner[0]);
				VecNormalize(Normal, Normal);
				VecScale(Normal, Normal, -1); // inner normal

			//	a tangental
				MathVector<WorldDim> Tangental(0.0);
				Tangental[0] = Normal[dim-1];
				Tangental[dim-1] = -Normal[0];

			//	quadrature points
				const number* vWeight = rSideQuadRule.weights();
				const size_t nip = rSideQuadRule.size();
				std::vector<MathVector<dim> > vLocalIP(nip);
				std::vector<MathVector<dim> > vGlobalIP(nip);

				DimReferenceMapping<dim-1, dim>& map
					= ReferenceMappingProvider::get<dim-1, dim>(sideRoid, vLocalSideCorner);

				for(size_t ip = 0; ip < nip; ++ip)
					map.local_to_global(vLocalIP[ip], rSideQuadRule.point(ip));

				for(size_t ip = 0; ip < nip; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]);

			//	compute transformation matrices
				DimReferenceMapping<dim-1, dim>& map2
					= ReferenceMappingProvider::get<dim-1, dim>(sideRoid, vSideCorner);
				std::vector<MathMatrix<dim-1, WorldDim> > vJT(nip);
				map2.jacobian_transposed(&(vJT[0]), rSideQuadRule.points(), nip);

				std::vector<MathMatrix<dim, WorldDim> > vElemJT(nip);
				rMapping.jacobian_transposed(&(vElemJT[0]), &vLocalIP[0], nip);

			//	loop integration points
				for(size_t ip = 0; ip < nip; ++ip)
				{
				// 	1. Interpolate Functional Matrix of velocity at ip
					std::vector<MathVector<dim> > vvLocGradV[dim];
					std::vector<MathVector<dim> > vvGradV[dim];
					MathMatrix<dim, dim> JTInv;
					Inverse(JTInv, vElemJT[ip]);
					for(int d1 = 0; d1 < dim; ++d1){
						const LocalShapeFunctionSet<dim>& rTrialSpaceP =
								LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[d1]);
						rTrialSpaceP.grads(vvLocGradV[d1], vLocalIP[ip]);

						vvGradV[d1].resize(vvLocGradV[d1].size());
						for(size_t sh = 0; sh < vvGradV[d1].size(); ++sh)
							MatVecMult(vvGradV[d1][sh], JTInv, vvLocGradV[d1][sh]);
					}

					MathMatrix<dim, dim> gradVel;
					for(int d1 = 0; d1 < dim; ++d1){
						for(int d2 = 0; d2 <dim; ++d2){
							gradVel(d1, d2) = 0.0;
							for(size_t sh = 0; sh < vvValue[d1].size(); ++sh)
								gradVel(d1, d2) += vvValue[d1][sh] * vvGradV[d1][sh][d2];
						}
					}

				//	1. Interpolate pressure at ip
					const LocalShapeFunctionSet<dim>& rTrialSpaceP =
							LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[_P_]);
					std::vector<number> vShapeP;
					rTrialSpaceP.shapes(vShapeP, vLocalIP[ip]);

					number pressure = 0.0;
					for(size_t sh = 0; sh < vvValue[_P_].size(); ++sh)
						pressure += vShapeP[sh] * vvValue[_P_][sh];

				//	2. Compute flux
					MathVector<dim> diffFlux;
					MatVecMult(diffFlux, gradVel, Normal);

				//	get quadrature weight
					const number weightIP = vWeight[ip];

				//	get determinate of mapping
					const number det = SqrtGramDeterminant(vJT[ip]);

				//	add contribution of integration point
					int_drag +=  weightIP * det *
								(kinVisco*density *VecDot(diffFlux, Tangental) * Normal[dim-1]
								                                    - pressure * Normal[0]);
					int_lift -=  weightIP * det *
								(kinVisco*density *VecDot(diffFlux, Tangental) * Normal[0]
								                                    + pressure * Normal[dim-1]);
				}
			} // end bf
		} // end elem
	} // end subsets

#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = int_drag;
		com.allreduce(&local, &int_drag, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
		local = int_lift;
		com.allreduce(&local, &int_lift, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif

//	return the summed integral contributions of all elements
	std::vector<number> vals(2);
	vals[0] = int_drag;
	vals[1] = int_lift;
	return vals;
}


} // end namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__NAVIER_STOKES_TOOLS__ */
