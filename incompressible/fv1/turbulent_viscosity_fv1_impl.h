/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_FV1_DATA_IMPL_H_
#define __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_FV1_DATA_IMPL_H_

#include "turbulent_viscosity_fv1.h"

namespace ug{
namespace NavierStokes{

// transfer by injection
template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>::transferToLowerLevels(aVertexNumber& aaData,ApproximationSpace<domain_type>& approximationSpace){
	for(int si = 0; si < approximationSpace.domain()->subset_handler()->num_subsets(); ++si){
		// transfer to lower levels, averaging over child edges (2d) / child faces (3d)
		for (size_t lev=approximationSpace.num_levels()-2; ;lev--){
			const DoFDistribution& lDD = *approximationSpace.dof_distribution(GridLevel(lev, GridLevel::LEVEL));
			const MultiGrid& grid = *lDD.multi_grid();
			typedef typename DoFDistribution::traits<Vertex>::const_iterator coarseLevelVertexIter;
			coarseLevelVertexIter clvIter, clvIterEnd;
			clvIter = lDD.template begin<Vertex>(si);
			clvIterEnd = lDD.template end<Vertex>(si);
			for (;clvIter != clvIterEnd;clvIter++){
				Vertex* vertex = *clvIter;
				aaData[vertex] += aaData[grid.get_child<Vertex>(vertex, 0)];
			}
			if (lev==0) break;
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>::fillAttachment(aVertexDimVector& aaU,SmartPtr<TGridFunction> u){
	//	get domain
	domain_type& domain = *u->domain().get();
	//	create Multiindex
	std::vector<DoFIndex> multInd;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		ElemIterator iter = u->template begin<Vertex>(si);
		ElemIterator iterEnd = u->template end<Vertex>(si);
		for(  ;iter !=iterEnd; ++iter)
		{
			Vertex* vertex = iter;
			for (int d=0;d<dim;d++){
				u->dof_indices(vertex, d, multInd);
				aaU[vertex][d]=DoFRef(*u,multInd[0]);
			}
		}
	}
}

// go over all elements, interpolate data to barycenter, average by multiplying with corresponding element volume and deviding by complete adjacent element volume
template <typename TData, int dim, typename TImpl,typename TGridFunction>
template <typename VType>
void StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>::elementFilter(PeriodicAttachmentAccessor<Vertex,Attachment<VType> >& aaUHat,aVertexNumber& aaVol,const PeriodicAttachmentAccessor<Vertex,Attachment<VType> >& aaU){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	DimFV1Geometry<dim> geo;

	// set attachment values to zero
	SetAttachmentValues(aaUHat , m_uInfo->template begin<Vertex>(), m_uInfo->template end<Vertex>(), 0);
	SetAttachmentValues(aaVol , m_uInfo->template begin<Vertex>(), m_uInfo->template end<Vertex>(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

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
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			//	memory for shapes
			std::vector<number> vShape;

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			rTrialSpace.shapes(vShape, localBary);

			size_t noc = elem->num_vertices();

			VType value;
			value = 0;
			number elementVolume = 0;


			for (size_t co=0;co<noc;co++){
				const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(co);
				VType localValue = aaU[elem->vertex(co)];
				//for debug UG_LOG(localValue << "\n");
				localValue *= vShape[co];
				value += localValue;
				elementVolume += scv.volume();
			}
			value *= elementVolume;
			for (size_t co=0;co<noc;co++){
				aaUHat[elem->vertex(co)] += value;
				aaVol[elem->vertex(co)] += elementVolume;
			}
		}
	}
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		VertexIterator vertexIter = m_uInfo->template begin<Vertex>(si);
		VertexIterator vertexIterEnd = m_uInfo->template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			Vertex* vertex = *vertexIter;
			if (pbm && pbm->is_slave(vertex)) continue;
			aaUHat[vertex]/=(number)aaVol[vertex];
		}
	}
}

// go over all elements, interpolate data to barycenter, average by multiplying with corresponding element volume and deviding by complete adjacent element volume
template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>::elementFilter(aVertexDimVector& aaUHat,aVertexNumber& aaVol,SmartPtr<TGridFunction> u){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	DimFV1Geometry<dim> geo;

	 std::vector<DoFIndex> multInd;

	// set attachment values to zero
	SetAttachmentValues(aaUHat , m_uInfo->template begin<Vertex>(), m_uInfo->template end<Vertex>(), 0);
	SetAttachmentValues(aaVol , m_uInfo->template begin<Vertex>(), m_uInfo->template end<Vertex>(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

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
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			//	memory for shapes
			std::vector<number> vShape;

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			rTrialSpace.shapes(vShape, localBary);

			size_t noc = elem->num_vertices();

			MathVector<dim> value;
			value = 0;
			number elementVolume = 0;
			for (size_t co=0;co<noc;co++){
				const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(co);
				MathVector<dim> localValue;
				for (int d=0;d<dim;d++){
					u->dof_indices(elem->vertex(co), d, multInd);
					localValue[d]=DoFRef(*u,multInd[0]);
				}
				//for debug UG_LOG("localValue=" << localValue << "\n");
				//for debug UG_LOG("vShape=" << vShape[s] << "\n");
				localValue *= vShape[co];
				value += localValue;
				elementVolume += scv.volume();
			}
			//for debug UG_LOG("value=" << value << " vol=" << elementVolume << "\n");
			value *= elementVolume;
			for (size_t co=0;co<noc;co++){
				aaVol[elem->vertex(co)]+=elementVolume;
				aaUHat[elem->vertex(co)] += value;
			}
		}
	}
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		VertexIterator vertexIter = m_uInfo->template begin<Vertex>(si);
		VertexIterator vertexIterEnd = m_uInfo->template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			Vertex* vertex = *vertexIter;
			if (pbm && pbm->is_slave(vertex)) continue;
			aaUHat[vertex]/=(number)aaVol[vertex];
		}
	}
}

// go over all elements, interpolate data to scv barycenter, average by multiplying with corresponding scv volume and deviding by volume of complete control volume
template <typename TData, int dim, typename TImpl,typename TGridFunction>
template <typename VType>
void StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>::scvFilter(PeriodicAttachmentAccessor<Vertex,Attachment<VType> >& aaUHat,aVertexNumber& aaVol,const PeriodicAttachmentAccessor<Vertex,Attachment<VType> >& aaU){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	DimFV1Geometry<dim> geo;

	// set attachment values to zero
	SetAttachmentValues(aaUHat , m_uInfo->template begin<Vertex>(), m_uInfo->template end<Vertex>(), 0);
	SetAttachmentValues(aaVol , m_uInfo->template begin<Vertex>(), m_uInfo->template end<Vertex>(), 0);

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
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

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
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			size_t noc = elem->num_vertices();

			MathVector<dim> scvLocalBary;
			for (size_t co=0;co<noc;co++){
				const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(co);

				scvLocalBary = 0;
				// compute barycenter of scv
				for (size_t i=0;i<scv->num_corners();i++){
					scvLocalBary += scv->loal_corner(i);
				}
				scvLocalBary/=(number)(scv->num_corners());
				//	memory for shapes
				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, scvLocalBary);
				VType localValue = 0;

				for (size_t j=0;j<noc;j++){
					localValue += vShape[j]*aaU[elem->vertex(j)];
				}

				localValue *= scv.volume();
				aaVol[elem->vertex(co)]  += scv.volume();
				aaUHat[elem->vertex(co)] += localValue;
			}
		}
	}
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		VertexIterator vertexIter = m_uInfo->template begin<Vertex>(si);
		VertexIterator vertexIterEnd = m_uInfo->template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			Vertex* vertex = *vertexIter;
			if (pbm && pbm->is_slave(vertex)) continue;
			aaUHat[vertex]/=(number)aaVol[vertex];
		}
	}
}

// go over all elements, interpolate data to scv barycenter, average by multiplying with corresponding scv volume and deviding by volume of complete control volume
template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>::scvFilter(aVertexDimVector& aaUHat,aVertexNumber& aaVol,SmartPtr<TGridFunction> u){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	DimFV1Geometry<dim> geo;

	// set attachment values to zero
	SetAttachmentValues(aaUHat , m_uInfo->template begin<Vertex>(), m_uInfo->template end<Vertex>(), 0);
	SetAttachmentValues(aaVol , m_uInfo->template begin<Vertex>(), m_uInfo->template end<Vertex>(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	//	create Multiindex
	std::vector<DoFIndex> multInd;

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

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
					LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			size_t noc = elem->num_vertices();

			static const size_t MaxNumVerticesOfElem = 10;

			typedef MathVector<dim> MVD;
			std::vector<MVD> uValue(MaxNumVerticesOfElem);

			for (size_t co=0;co<noc;co++){
				for (int d=0;d<dim;d++){
					u->dof_indices(elem->vertex(co), d, multInd);
					uValue[co][d]=DoFRef(*u,multInd[0]);
				}
			};

			MathVector<dim> scvLocalBary;
			for (size_t co=0;co<noc;co++){
				const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(co);
				scvLocalBary = 0;
				// compute barycenter of scv
				for (size_t i=0;i<scv.num_corners();i++){
					scvLocalBary += scv.local_corner(i);
				}
				scvLocalBary/=(number)(scv.num_corners());
				//	memory for shapes
				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, scvLocalBary);
				MathVector<dim> localValue = 0;

				for (size_t j=0;j<noc;j++)
					for (int d=0;d<dim;d++)
						localValue += vShape[j]*uValue[j];

				localValue *= scv.volume();
				aaVol[elem->vertex(co)]  += scv.volume();
				aaUHat[elem->vertex(co)] += localValue;
			}
		}
	}
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		VertexIterator vertexIter = m_uInfo->template begin<Vertex>(si);
		VertexIterator vertexIterEnd = m_uInfo->template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			Vertex* vertex = *vertexIter;
			if (pbm && pbm->is_slave(vertex)) continue;
			aaUHat[vertex]/=(number)aaVol[vertex];
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>::assembleDeformationTensor(aVertexTensor& aaDefTensor,aVertexNumber& aaVol,SmartPtr<TGridFunction> u){
	//	get domain
	domain_type& domain = *u->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();

	//	create Multiindex
	std::vector<DoFIndex> multInd;

	DimFV1Geometry<dim> geo;

	// add boundary subsets to enforce boundary subset computations in geo.update()
	for(size_t i = 0; i < this->m_turbZeroSg.size(); ++i){
		geo.add_boundary_subset(this->m_turbZeroSg[i]);
	}

	// set attachment values to zero
	SetAttachmentValues(aaDefTensor , u->template begin<Vertex>(), u->template end<Vertex>(), 0);
	SetAttachmentValues(aaVol , u->template begin<Vertex>(), u->template end<Vertex>(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			//	get iterators
			ElemIterator iter = u->template begin<elem_type>(si);
			ElemIterator iterEnd = u->template end<elem_type>(si);

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
				};

				//	evaluate finite volume geometry
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

				static const size_t MaxNumVerticesOfElem = 10;

				typedef MathVector<dim> MVD;
				std::vector<MVD> uValue(MaxNumVerticesOfElem);
				MVD ipVelocity;

				size_t noc = elem->num_vertices();

				size_t nip = geo.num_scvf();

				for (size_t co=0;co < noc;co++)
				{
					const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(co);
					for (int d=0;d<dim;d++){
						u->dof_indices(elem->vertex(co), d, multInd);
						uValue[co][d]=DoFRef(*u,multInd[0]);
					}
					aaVol[elem->vertex(co)] += scv.volume();
				}

				for (size_t ip=0;ip<nip;ip++){
					// 	get current SCVF
					ipVelocity = 0;
					const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
					for (size_t co=0;co < noc;co++){
						for (int d=0;d<dim;d++){
						    ipVelocity[d] += scvf.shape(co)*uValue[co][d];
						};
					};
					dimMat ipDefTensorFlux;
					ipDefTensorFlux = 0;
					for (int i=0;i<dim;i++){
						for (int j=0;j<dim;j++){
							ipDefTensorFlux[i][j]= 0.5 * (ipVelocity[i] * scvf.normal()[j] + ipVelocity[j] * scvf.normal()[i]);
						}
					}
					aaDefTensor[elem->vertex(scvf.from())]+=ipDefTensorFlux;
					aaDefTensor[elem->vertex(scvf.to())]-=ipDefTensorFlux;
				}
				for(size_t sgi = 0; sgi < this->m_turbZeroSg.size(); ++sgi){
					const size_t sgsi=this->m_turbZeroSg[sgi];
					if (geo.num_bf(sgsi) == 0) continue;
					for(size_t bfi = 0; bfi < geo.num_bf(sgsi); ++bfi){
						const typename DimFV1Geometry<dim>::BF& bf = geo.bf(sgsi, bfi);
						const size_t nodeID = bf.node_id();
						for (int i=0;i<dim;i++)
							for (int j=0;j<dim;j++){
								aaDefTensor[elem->vertex(nodeID)][i][j] += 0.5 * (uValue[nodeID][i] * bf.normal()[j] + uValue[nodeID][j] * bf.normal()[i]);
							}
					}
				}
			}
	}
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		VertexIterator vertexIter = u->template begin<Vertex>(si);
		VertexIterator vertexIterEnd = u->template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			Vertex* vertex = *vertexIter;
			if (pbm && pbm->is_slave(vertex)) continue;
			aaDefTensor[vertex]/=(number)aaVol[vertex];
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>::scaleTensorByNorm(aVertexTensor& aaTensor){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		VertexIterator vertexIter = m_uInfo->template begin<Vertex>(si);
		VertexIterator vertexIterEnd = m_uInfo->template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			Vertex* vertex = *vertexIter;
			if (pbm && pbm->is_slave(vertex)) continue;
			aaTensor[vertex]*=(number)FNorm(aaTensor[vertex]);
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>::assembleDeformationTensor(aVertexTensor& aaDefTensor,aVertexNumber& aaVol,aVertexDimVector aaU){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	//	create Multiindex
	std::vector<DoFIndex> multInd;

	DimFV1Geometry<dim> geo;

	// add boundary subsets to enforce boundary subset computations in geo.update()
	for(size_t i = 0; i < this->m_turbZeroSg.size(); ++i){
		geo.add_boundary_subset(this->m_turbZeroSg[i]);
	}

	SetAttachmentValues(aaDefTensor , m_uInfo->template begin<Vertex>(), m_uInfo->template end<Vertex>(), 0);
	SetAttachmentValues(aaVol , m_uInfo->template begin<Vertex>(), m_uInfo->template end<Vertex>(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			//	get iterators
			ElemIterator iter = m_uInfo->template begin<elem_type>(si);
			ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

			//for debug UG_LOG("|||||||||||||||||||| si = " << si << "\n");

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
					// //for debug UG_LOG("co_coord(" << i<< "+1,:)=" << coCoord[i] << "\n");
				};

				//	evaluate finite volume geometry
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

				static const size_t MaxNumVerticesOfElem = 10;

				typedef MathVector<dim> MVD;
				std::vector<MVD> uValue(MaxNumVerticesOfElem);
				MVD ipVelocity;

				UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

				size_t noc = elem->num_vertices();

				size_t nip = geo.num_scvf();

				for (size_t co=0;co < noc;co++)
				{
					const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(co);
					uValue[co]=aaU[elem->vertex(co)];
					//for debug UG_LOG("uvalue(" << s << ")=" << uValue[s] << "\n");
					aaVol[elem->vertex(co)] += scv.volume();
				}

				for (size_t ip=0;ip<nip;ip++){
					// 	get current SCVF
					ipVelocity = 0;
					const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
					for (size_t co=0;co < noc;co++){
						for (int d=0;d<dim;d++){
						    ipVelocity[d] += scvf.shape(co)*uValue[co][d];
						};
					};
					dimMat ipDefTensorFlux;
					ipDefTensorFlux = 0;
					for (int d=0;d<dim;d++){
						for (int j=0;j<dim;j++){
							ipDefTensorFlux[d][j]= 0.5 * (ipVelocity[d] * scvf.normal()[j] + ipVelocity[j] * scvf.normal()[d]);
						}
					}
					aaDefTensor[elem->vertex(scvf.from())]+=ipDefTensorFlux;
					aaDefTensor[elem->vertex(scvf.to())]-=ipDefTensorFlux;
				}
				for(size_t sgi = 0; sgi < this->m_turbZeroSg.size(); ++sgi){
					const size_t sgsi=this->m_turbZeroSg[sgi];
					if (geo.num_bf(sgsi) == 0) continue;
					for(size_t bfi = 0; bfi < geo.num_bf(sgsi); ++bfi){
						const typename DimFV1Geometry<dim>::BF& bf = geo.bf(sgsi, bfi);
						const size_t nodeID = bf.node_id();
						for (int i=0;i<dim;i++)
							for (int j=0;j<dim;j++){
								aaDefTensor[elem->vertex(nodeID)][i][j] += 0.5 * (uValue[nodeID][i] * bf.normal()[j] + uValue[nodeID][j] * bf.normal()[i]);
							}
					}
				}
			}
	}
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		VertexIterator vertexIter = m_uInfo->template begin<Vertex>(si);
		VertexIterator vertexIterEnd = m_uInfo->template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			Vertex* vertex = *vertexIter;
			if (pbm && pbm->is_slave(vertex)) continue;
			aaDefTensor[vertex]/=(number)aaVol[vertex];
		}
	}
}

// Frobenius norm of dim x dim matrix
template <typename TData, int dim, typename TImpl,typename TGridFunction>
number StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>::FNorm(MathMatrix<dim,dim> M){
	number norm=0;
	for (int d1=0;d1<dim;d1++)
		for (int d2=0;d2<dim;d2++){
			norm +=  M[d1][d2] * M[d1][d2];
		}
	return sqrt(2.0*norm);
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>::addUiUjTerm(aVertexTensor& aaResult,const number factor,aVertexDimVector aaU){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		VertexIterator vertexIter = m_uInfo->template begin<Vertex>(si);
		VertexIterator vertexIterEnd = m_uInfo->template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			Vertex* vertex = *vertexIter;
			if (pbm && pbm->is_slave(vertex)) continue;
			dimMat Tij;
			for (int d1=0;d1 < dim;d1++)
				for (int d2=0;d2 < dim;d2++)
					Tij[d1][d2] = aaU[vertex][d1]*aaU[vertex][d2];
			Tij*=factor;
			aaResult[vertex]+=Tij;
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>::addUiUjTerm(aVertexTensor& aaResult,const number factor,SmartPtr<TGridFunction> u){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	//	create Multiindex
	std::vector<DoFIndex> multInd;

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		VertexIterator vertexIter = m_uInfo->template begin<Vertex>(si);
		VertexIterator vertexIterEnd = m_uInfo->template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			Vertex* vertex = *vertexIter;
			if (pbm && pbm->is_slave(vertex)) continue;
			dimMat Tij;
			MathVector<dim> uValue;
			for (int d=0;d<dim;d++){
				u->dof_indices(vertex, d, multInd);
				uValue[d]=DoFRef(*u,multInd[0]);
			}
			for (int d1=0;d1 < dim;d1++)
				for (int d2=0;d2 < dim;d2++)
					Tij[d1][d2] = uValue[d1]*uValue[d2];
			Tij*=factor;
			aaResult[vertex]+=Tij;
		}
	}
}

template<typename TGridFunction>
void FV1SmagorinskyTurbViscData<TGridFunction>::update(){
	//	get domain of grid function
	domain_type& domain = *m_u->domain().get();
	SetAttachmentValues(m_acTurbulentViscosity, m_grid->template begin<Vertex>(), m_grid->template end<Vertex>(), 0);
	//	coord and vertex array
//	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
//	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	// assemble deformation tensor fluxes
	this->assembleDeformationTensor(m_acDeformation,m_acVolume,m_u);
	// compute turbulent viscosity , loop over vertices
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		if ((this->m_turbZeroSg.size()!=0) && (this->m_turbZeroSg.contains(si)==true)) continue;
		VertexIterator vertexIter = m_u->template begin<Vertex>(si);
		VertexIterator vertexIterEnd = m_u->template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			//	get Elem
			Vertex* vertex = *vertexIter;
			if (m_pbm && m_pbm->is_slave(vertex)){
				continue;
			}
			number delta = m_acVolume[vertex];
			// for possible other choices of delta see Fr�hlich p 160
			delta = pow(delta,(number)1.0/(number)dim);
			number tensorNorm = this->FNorm(m_acDeformation[vertex]);
			m_acTurbulentViscosity[vertex] = m_c * delta*delta * tensorNorm;
		}
	}
	// transfer attachment data to lower levels
	this->transferToLowerLevels(m_acTurbulentViscosity,*m_spApproxSpace);
}

template<typename TGridFunction>
void FV1DynamicTurbViscData<TGridFunction>::update(){
	//	get domain of grid function
	domain_type& domain = *m_u->domain().get();
	//	get position accessor
	// for debug typedef typename domain_type::position_accessor_type position_accessor_type;
	// for debug const position_accessor_type& posAcc = domain.position_accessor();

	// initialize attachment values with 0
//	SetAttachmentValues(m_acDeformation , m_grid->template begin<Vertex>(), m_grid->template end<Vertex>(), 0);
	SetAttachmentValues(m_acTurbulentViscosity, m_grid->template begin<Vertex>(), m_grid->template end<Vertex>(), 0);
//	SetAttachmentValues(m_acVolume,m_grid->template begin<Vertex>(), m_grid->template end<Vertex>(), 0);
	SetAttachmentValues(m_acTurbulentC,m_grid->template begin<Vertex>(), m_grid->template end<Vertex>(), 0);
//	SetAttachmentValues(m_acVolumeHat,m_grid->template begin<Vertex>(), m_grid->template end<Vertex>(), 0);
//	SetAttachmentValues(m_acUHat,m_grid->template begin<Vertex>(), m_grid->template end<Vertex>(), 0);
//	SetAttachmentValues(m_acDeformationHat,m_grid->template begin<Vertex>(), m_grid->template end<Vertex>(), 0);
//	SetAttachmentValues(m_acLij,m_grid->template begin<Vertex>(), m_grid->template end<Vertex>(), 0);
	SetAttachmentValues(m_acMij,m_grid->template begin<Vertex>(), m_grid->template end<Vertex>(), 0);

	// compute Lij term \hat{u_i u_j} - \hat{u_i} \hat{u_j}
	// \hat{u}
	this->elementFilter(m_acUHat,m_acVolumeHat,m_u);
	// use Mij attachment to store first Lij part
	// u_i u_j
	this->addUiUjTerm(m_acMij,1.0,m_u);
	// \hat{u_i u_j}
	this->elementFilter(m_acLij,m_acVolumeHat,m_acMij);
	// \hat{u_i u_j} - \hat{u_i} \hat{u_j}
	this->addUiUjTerm(m_acLij,-1.0,m_acUHat);

	// Mij term
	// first term |\hat{S}| \hat{S}
	// assemble \hat{S} using \hat{u}
	this->assembleDeformationTensor(m_acDeformationHat,m_acVolume,m_acUHat);
	// normalize \hat{S}
	this->scaleTensorByNorm(m_acDeformationHat);
	// Mij second term \hat{|S|S}
	// compute S
	this->assembleDeformationTensor(m_acDeformation,m_acVolumeHat,m_u);
	// compute |S| S
	this->scaleTensorByNorm(m_acDeformation);
	// filter |S| S
	//for debug UG_LOG("------------------------------------------------------\n");
	this->elementFilter(m_acMij,m_acVolumeHat,m_acDeformation);

	bool use_filter = false;

	//	create Multiindex
	std::vector<DoFIndex> multInd;

	// complete Mij term computation by scaling and adding the two terms,
	// solve the local least squares problem and compute local c and local turbulent viscosity
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		if (use_filter==false)
			if ((this->m_turbZeroSg.size()!=0) && (this->m_turbZeroSg.contains(si)==true)) continue;
		VertexIterator vertexIter = m_u->template begin<Vertex>(si);
		VertexIterator vertexIterEnd = m_u->template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++){
			Vertex* vertex = *vertexIter;
			if (m_pbm && m_pbm->is_slave(vertex)) continue;
			// use c to compute turbulent viscosity
			number delta = m_acVolume[vertex];
			delta = pow(delta,(number)1.0/(number)dim);
			number deltaHat = m_acVolumeHat[vertex];
			deltaHat = pow(deltaHat,(number)1.0/(number)dim);
			m_u->dof_indices(vertex, 0, multInd);
			m_u->dof_indices(vertex, 1, multInd);
			m_acDeformationHat[vertex] *= -2*deltaHat*deltaHat;
			m_acMij[vertex] *= 2*delta*delta;
			m_acMij[vertex] += m_acDeformationHat[vertex];
			//for debug UG_LOG("Mij " << FNorm(m_acMij[vertex]) << "\n");
			for (int d1=0;d1<dim;d1++){
				for (int d2=0;d2<dim;d2++){
					//for debug UG_LOG(m_acMij[vertex][d1][d2] << " ");
				}
				//for debug UG_LOG("\n");
			}
			//for debug UG_LOG("Lij " << FNorm(m_acLij[vertex]) << "\n");
			for (int d1=0;d1<dim;d1++){
				for (int d2=0;d2<dim;d2++){
					//for debug UG_LOG(m_acLij[vertex][d1][d2] << " ");
				}
				//for debug UG_LOG("\n");
			}
			// compute local c
			// solve least squares problem
			number c = 0;
			for (int d1=0;d1<dim;d1++)
				for (int d2=0;d2<dim;d2++)
					c += m_acLij[vertex][d1][d2]*m_acMij[vertex][d1][d2];
			number denom=0;
			//for debug UG_LOG("c=" << c << "\n");
			for (int d1=0;d1<dim;d1++)
				for (int d2=0;d2<dim;d2++)
					denom += m_acMij[vertex][d1][d2]*m_acMij[vertex][d1][d2];
			if (denom>1e-15)
				c/=(number)denom;
			else c=0;

			if (m_spaceFilter==false){
				if (m_timeFilter==false){
					m_acTurbulentViscosity[vertex] = c * delta*delta * this->FNorm(m_acDeformation[vertex]);
				} else {
					m_acTurbulentC[vertex]= (m_timeFilterEps * c + (1-m_timeFilterEps)*m_acTurbulentC[vertex]);
					m_acTurbulentViscosity[vertex] = m_acTurbulentC[vertex] * delta*delta * this->FNorm(m_acDeformation[vertex]);
				}
				if (m_acTurbulentViscosity[vertex]+m_viscosityNumber<m_small) m_acTurbulentViscosity[vertex] = m_viscosityNumber + m_small;			}
			else{
				// store c in viscosity array
				m_acTurbulentViscosity[vertex] = c;
			}
		}
	}
	if (m_spaceFilter==true){
		// filter c
		if (m_timeFilter==false)
			this->elementFilter(m_acTurbulentC,m_acVolumeHat,m_acTurbulentViscosity);
		else
			// store c in volumeHat array
			this->elementFilter(m_acVolumeHat,m_acVolumeHat,m_acTurbulentViscosity);
		// compute turbulent viscosity
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
		{
			//for debug UG_LOG("si = " << si << "\n");
			if ((this->m_turbZeroSg.size()!=0) && (this->m_turbZeroSg.contains(si)==true)) continue;
			VertexIterator vertexIter = m_u->template begin<Vertex>(si);
			VertexIterator vertexIterEnd = m_u->template end<Vertex>(si);
			for(  ;vertexIter !=vertexIterEnd; vertexIter++){
				Vertex* vertex = *vertexIter;
				if (m_pbm && m_pbm->is_slave(vertex)) continue;
				number delta = m_acVolume[vertex];
				delta = pow(delta,(number)1.0/(number)dim);
				if (m_timeFilter==true)
					// time averaging, note that c has been stored in m_acVolumeHat
					m_acTurbulentC[vertex]= (m_timeFilterEps * m_acVolumeHat[vertex] + (1-m_timeFilterEps)*m_acTurbulentC[vertex]);
				m_acTurbulentViscosity[vertex] = m_acTurbulentC[vertex] * delta*delta * this->FNorm(m_acDeformation[vertex]);
				if (m_acTurbulentViscosity[vertex]+m_viscosityNumber<m_small) m_acTurbulentViscosity[vertex] = m_viscosityNumber+m_small;
				//for debug UG_LOG("nu_t = " << m_acTurbulentViscosity[vertex]  << " c = " << m_acTurbulentC[vertex] << " delta = " << delta << " co=[" << 0.5*(posAcc[vertex->vertex(0)][0] + posAcc[vertex->vertex(1)][0]) << "," << 0.5*(posAcc[vertex->vertex(0)][1] + posAcc[vertex->vertex(1)][1]) << "]\n");
			}
		}
	}
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_FV1_DATA_IMPL_H_ */
