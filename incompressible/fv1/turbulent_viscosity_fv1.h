/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__NAVIER_STOKES__INCOMPRESSIBLE__FV1__TURBULENT_VISCOSITY_DATA_FV1__
#define __H__UG__NAVIER_STOKES__INCOMPRESSIBLE__FV1__TURBULENT_VISCOSITY_DATA_FV1__

#include "common/common.h"

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton_update_interface.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_grid/tools/subset_group.h"
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "lib_grid/algorithms/attachment_util.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace NavierStokes{

/**
concept derived from grid_function_user_data.h
*/
template <typename TData, int dim, typename TImpl,typename TGridFunction>
class StdTurbulentViscosityDataFV1
	: 	public StdUserData<StdTurbulentViscosityDataFV1<TData,dim,TImpl,TGridFunction>, TData,dim>
{
		///	domain type
		typedef typename TGridFunction::domain_type domain_type;

		///	algebra type
		typedef typename TGridFunction::algebra_type algebra_type;

		/// position accessor type
		typedef typename domain_type::position_accessor_type position_accessor_type;

		///	grid type
		typedef typename domain_type::grid_type grid_type;

		/// element type
		typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

		/// element iterator
		typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

		/// vertex iterator
		typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

		/// attachment accessor types
		typedef MathMatrix<dim,dim> dimMat;
		typedef Attachment<dimMat> ATensor;

		typedef MathVector<dim> vecDim;
		typedef Attachment<vecDim> AMathVectorDim;

		typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;
		typedef PeriodicAttachmentAccessor<Vertex,ATensor > aVertexTensor;
		typedef PeriodicAttachmentAccessor<Vertex,AMathVectorDim > aVertexDimVector;
	public:
		////////////////
		// one value
		////////////////
		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si) const
		{
			UG_THROW("StdTurbulentViscosityDataFV1: Need element.");
		}

		virtual void operator() (TData vValue[],
		                         const MathVector<dim> vGlobIP[],
		                         number time, int si, const size_t nip) const
		{
			UG_THROW("StdTurbulentViscosityDataFV1: Need element.");
		}

		template <int refDim>
		void evaluate(TData vValue[],
						const MathVector<dim> vGlobIP[],
						number time, int si,
						GridObject* elem,
						const MathVector<dim> vCornerCoords[],
						const MathVector<refDim> vLocIP[],
						const size_t nip,
						LocalVector* u,
						const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<refDim>(vValue,vGlobIP,time,si,elem,
			                                    vCornerCoords,vLocIP,nip,u,vJT);
		}

		virtual void compute(LocalVector* u, GridObject* elem,
		                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
		{
			const int si = this->subset();
			for(size_t s = 0; s < this->num_series(); ++s)
				getImpl().template evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
			                  elem, vCornerCoords, this->template local_ips<dim>(s),
			                  this->num_ip(s), u);
		}

		virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
		                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
		{
			const int si = this->subset();
			for(size_t s = 0; s < this->num_series(); ++s)
				getImpl().template evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
			                  elem, vCornerCoords, this->template local_ips<dim>(s),
			                  this->num_ip(s), &(u->solution(this->time_point(s))));
		}

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool continuous() const {return false;}

	///	returns if grid function is needed for evaluation
		virtual bool requires_grid_fct() const {return true;}

		void assembleDeformationTensor(aVertexTensor& aaDefTensor,ATensor& aa_aDefTensor,aVertexNumber& aaVol,ANumber& aa_aVol,SmartPtr<TGridFunction> u);
		void assembleDeformationTensor(aVertexTensor& aaDefTensor,ATensor& aa_aDefTensor,aVertexNumber& aaVol,ANumber& aa_aVol,aVertexDimVector aaU);

		number FNorm(MathMatrix<dim,dim> M);

		void addUiUjTerm(aVertexTensor& aaDefTensor,const number factor,aVertexDimVector aaU);
		void addUiUjTerm(aVertexTensor& aaDefTensor,const number factor,SmartPtr<TGridFunction> u);

		template <typename VType>
		void elementFilter(PeriodicAttachmentAccessor<Vertex,Attachment<VType> >& aaUHat, Attachment<VType>& aa_aUHat,aVertexNumber& aaVol, ANumber& aa_aVol,const PeriodicAttachmentAccessor<Vertex,Attachment<VType> >& aaU);

		void elementFilter(aVertexDimVector& aaUHat, AMathVectorDim& aa_aUHat,aVertexNumber& aaVol, ANumber& aa_aVol,SmartPtr<TGridFunction> u);

		template <typename VType>
		void scvFilter(PeriodicAttachmentAccessor<Vertex,Attachment<VType> >& aaUHat,aVertexNumber& aaVol,const PeriodicAttachmentAccessor<Vertex,Attachment<VType> >& aaU);

		void scvFilter(aVertexDimVector& aaUHat,aVertexNumber& aaVol,SmartPtr<TGridFunction> u);

		void fillAttachment(aVertexDimVector& aaU,SmartPtr<TGridFunction> u);

		void transferToLowerLevels(aVertexNumber& aaData,ApproximationSpace<domain_type>& approximationSpace);

		void scaleTensorByNorm(aVertexTensor& aaTensor);

		// set non-periodic boundaries so that viscosity can be set to zero there
		void setTurbulenceZeroBoundaries(const char* subsets){
			try{
				m_turbZeroSg = m_uInfo->subset_grp_by_name(subsets);
			}UG_CATCH_THROW("ERROR while parsing Subsets.");
		}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}

		// grid function
		SmartPtr<TGridFunction> m_uInfo;

		// subset group
		SubsetGroup m_turbZeroSg;
};


template <typename TGridFunction>
class FV1SmagorinskyTurbViscData
	: public StdTurbulentViscosityDataFV1<number, TGridFunction::dim,
	  	  	  FV1SmagorinskyTurbViscData<TGridFunction>,TGridFunction >, virtual public INewtonUpdate
{
	///	domain type
		typedef typename TGridFunction::domain_type domain_type;

	///	algebra type
		typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
		typedef typename domain_type::position_accessor_type position_accessor_type;

	///	world dimension
		static const int dim = domain_type::dim;

	///	grid type
		typedef typename domain_type::grid_type grid_type;

	/// element type
		typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// element iterator
		typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
		typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

	/// attachment accessor types
		typedef MathMatrix<dim,dim> dimMat;
		typedef Attachment<dimMat> ATensor;

		typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;
		typedef PeriodicAttachmentAccessor<Vertex,ATensor > aVertexTensor;

	public:
		/**
		 * This method sets the kinematic viscosity value. Kinematic viscosity is added to turbulent viscosity in evaluation routine.
		 *
		 * \param[in]	data		kinematic Viscosity
		 */
		///	\{
		void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user){
			m_imKinViscosity = user;
		}
		void set_kinematic_viscosity(number val){
			set_kinematic_viscosity(make_sp(new ConstUserNumber<dim>(val)));
		}
	#ifdef UG_FOR_LUA
		void set_kinematic_viscosity(const char* fctName){
			set_kinematic_viscosity(LuaUserDataFactory<number, dim>::create(fctName));
		}
	#endif
		///	\}

	private:
		///	Data import for kinematic viscosity
		SmartPtr<CplUserData<number,dim> > m_imKinViscosity;

	private:
	// grid function
		SmartPtr<TGridFunction> m_u;
		
		grid_type* m_grid;

	//  turbulent viscosity attachment
		aVertexNumber m_acTurbulentViscosity;
		ANumber m_aTurbulentViscosity;
		
	//  volume attachment
		aVertexNumber m_acVolume;
		ANumber m_aVolume;

	//	deformation tensor attachment
		aVertexTensor m_acDeformation;
		ATensor m_aDeformation;
		
	//  Smagorinsky model parameter, typical values [0.01 0.1]
		number m_c;

		//	approximation space for level and surface grid
		SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

		//  periodic boundary manager
		PeriodicBoundaryManager* m_pbm;

	public:
	/// constructor
		FV1SmagorinskyTurbViscData(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct,number c = 0.05){
			m_c = c;
			m_u = spGridFct;
			this->m_uInfo = m_u;
			m_spApproxSpace = approxSpace;
			domain_type& domain = *m_u->domain().get();
			grid_type& grid = *domain.grid();
			m_grid = &grid;
			m_pbm = m_grid->periodic_boundary_manager();
			// attachments
			grid.template attach_to<Vertex>(m_aTurbulentViscosity);
			grid.template attach_to<Vertex>(m_aVolume);
			grid.template attach_to<Vertex>(m_aDeformation);
			// accessors
			m_acTurbulentViscosity.access(grid,m_aTurbulentViscosity);
			m_acVolume.access(grid,m_aVolume);
			m_acDeformation.access(grid,m_aDeformation);

		}
		
		virtual ~FV1SmagorinskyTurbViscData() {
			domain_type& domain = *m_u->domain().get();
			grid_type& grid = *domain.grid();
			grid.template detach_from<Vertex>(m_aTurbulentViscosity);
			grid.template detach_from<Vertex>(m_aVolume);
			grid.template detach_from<Vertex>(m_aDeformation);
		};

		void set_model_parameter(number c){
			m_c = c;
		}
		
		template <int refDim>
		inline void evaluate(number vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     GridObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     LocalVector* u,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
		//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();
			elem_type* element = static_cast<elem_type*>(elem);

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");

			size_t noc=element->num_vertices();

		//	get trial space
			try{
			const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));

		//	memory for shapes
			std::vector<number> vShape;

		//	loop ips
			for(size_t ip = 0; ip < nip; ++ip)
			{
			//	evaluate shapes at ip
				rTrialSpace.shapes(vShape, vLocIP[ip]);

			// 	compute solution at integration point
				vValue[ip] = 0.0;
				for(size_t sh = 0; sh < noc; ++sh)
				{
					 vValue[ip] += vShape[sh] * m_acTurbulentViscosity[element->vertex(sh)];
				}
			}

			}
			UG_CATCH_THROW("TurbulentViscosityData: trial space missing, Reference Object: "
				         <<roid<<", Trial Space: LAGRANGE 1, refDim="<<refDim);

			number kinViscValues[max_number_of_ips];
			(*m_imKinViscosity)(kinViscValues,
                    vGlobIP,
                    time, si,
                    elem,
                    vCornerCoords,
                    vLocIP,
                    nip,
                    u,
                    vJT);
			for (size_t ip=0;ip < nip;ip++){
				// UG_LOG("turbVis(" << ip << ")=" << vValue[ip] << "+" << kinViscValues[ip] << "\n");
				vValue[ip] += kinViscValues[ip];
			}
		}
		
		static const size_t max_number_of_ips = 20;
	
		int m_counter = 0;
		const int m_update = 2;
	

		void update();

};

template <typename TGridFunction>
class FV1DynamicTurbViscData
	: public StdTurbulentViscosityDataFV1<number, TGridFunction::dim,
	  	  	  FV1DynamicTurbViscData<TGridFunction>,TGridFunction >, virtual public INewtonUpdate
{
	///	domain type
		typedef typename TGridFunction::domain_type domain_type;

	///	algebra type
		typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
		typedef typename domain_type::position_accessor_type position_accessor_type;

	///	world dimension
		static const int dim = domain_type::dim;

	///	grid type
		typedef typename domain_type::grid_type grid_type;

	/// element type
		typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// element iterator
		typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
		typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

	/// attachment accessor types
		typedef MathMatrix<dim,dim> dimMat;
		typedef Attachment<dimMat> ATensor;

		typedef MathVector<dim> vecDim;
		typedef Attachment<vecDim> AMathVectorDim;

		typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;
		typedef PeriodicAttachmentAccessor<Vertex,ATensor > aVertexTensor;
		typedef PeriodicAttachmentAccessor<Vertex,AMathVectorDim > aVertexDimVector;

	public:
		/**
		 * This method sets the kinematic viscosity value. Kinematic viscosity is added to turbulent viscosity in evaluation routine.
		 *
		 * \param[in]	data		kinematic Viscosity
		 */
		///	\{
		void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user){
			m_imKinViscosity = user;
		}
		void set_kinematic_viscosity(number val){
			m_viscosityNumber = val;
			set_kinematic_viscosity(make_sp(new ConstUserNumber<dim>(val)));
		}
	#ifdef UG_FOR_LUA
		void set_kinematic_viscosity(const char* fctName){
			set_kinematic_viscosity(LuaUserDataFactory<number, dim>::create(fctName));
		}
	#endif
		///	\}

	private:
		///	Data import for kinematic viscosity
		SmartPtr<CplUserData<number,dim> > m_imKinViscosity;

		number m_viscosityNumber;

		static const number m_small;

	private:
	// grid function
		SmartPtr<TGridFunction> m_u;
		
		grid_type* m_grid;

	//  turbulent viscosity attachment
		aVertexNumber m_acTurbulentViscosity;
		ANumber m_aTurbulentViscosity;
		
	//  turbulent model parameter attachment
		aVertexNumber m_acTurbulentC;
		ANumber m_aTurbulentC;

	//  volume attachment
		aVertexNumber m_acVolume;
		ANumber m_aVolume;

	//  coarser grid volume attachment
		aVertexNumber m_acVolumeHat;
		ANumber m_aVolumeHat;

	//  filtered u attachment
		aVertexDimVector m_acUHat;
		AMathVectorDim m_aUHat;

	//	deformation tensor attachment
		aVertexTensor m_acDeformation;
		ATensor m_aDeformation;
		
	//	coarser grid deformation tensor attachment
		aVertexTensor m_acDeformationHat;
		ATensor m_aDeformationHat;

	//	Leonard tensor attachment
		aVertexTensor m_acLij;
		ATensor m_aLij;

	//	Mij tensor attachment
		aVertexTensor m_acMij;
		ATensor m_aMij;

	//	approximation space for level and surface grid
		SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  periodic boundary manager
		PeriodicBoundaryManager* m_pbm;

	public:
	/// constructor
		FV1DynamicTurbViscData(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
			m_u = spGridFct;
			this->m_uInfo = m_u;
			m_spApproxSpace = approxSpace;
			domain_type& domain = *m_u->domain().get();
			grid_type& grid = *domain.grid();
			m_grid = &grid;
			m_pbm = m_grid->periodic_boundary_manager();
			// attachments
			grid.template attach_to<Vertex>(m_aTurbulentViscosity);
			grid.template attach_to<Vertex>(m_aTurbulentC);
			grid.template attach_to<Vertex>(m_aVolume);
			grid.template attach_to<Vertex>(m_aVolumeHat);
			grid.template attach_to<Vertex>(m_aUHat);
			grid.template attach_to<Vertex>(m_aDeformation);
			grid.template attach_to<Vertex>(m_aDeformationHat);
			grid.template attach_to<Vertex>(m_aLij);
			grid.template attach_to<Vertex>(m_aMij);
			// accessors
			m_acTurbulentViscosity.access(grid,m_aTurbulentViscosity);
			m_acTurbulentC.access(grid,m_aTurbulentC);
			m_acVolume.access(grid,m_aVolume);
			m_acVolumeHat.access(grid,m_aVolumeHat);
			m_acUHat.access(grid,m_aUHat);
			m_acDeformation.access(grid,m_aDeformation);
			m_acDeformationHat.access(grid,m_aDeformationHat);
			m_acLij.access(grid,m_aLij);
			m_acMij.access(grid,m_aMij);
			// default setting for filtering of model constant c
			m_spaceFilter=true;
			m_timeFilter=false;
			m_timeFilterEps=1;
		}
		
		virtual ~FV1DynamicTurbViscData() {
			domain_type& domain = *m_u->domain().get();
			grid_type& grid = *domain.grid();
			grid.template detach_from<Vertex>(m_aTurbulentViscosity);
			grid.template detach_from<Vertex>(m_aTurbulentC);
			grid.template detach_from<Vertex>(m_aVolume);
			grid.template detach_from<Vertex>(m_aVolumeHat);
			grid.template detach_from<Vertex>(m_aUHat);
			grid.template detach_from<Vertex>(m_aDeformation);
			grid.template detach_from<Vertex>(m_aDeformationHat);
			grid.template detach_from<Vertex>(m_aLij);
			grid.template detach_from<Vertex>(m_aMij);
		};
		
		template <int refDim>
		inline void evaluate(number vValue[],
									 const MathVector<dim> vGlobIP[],
				                     number time, int si,
				                     GridObject* elem,
				                     const MathVector<dim> vCornerCoords[],
				                     const MathVector<refDim> vLocIP[],
				                     const size_t nip,
				                     LocalVector* u,
				                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
			elem_type* element = static_cast<elem_type*>(elem);

			size_t noc=element->num_vertices();

			//	get trial space
			try{
				const LocalShapeFunctionSet<refDim>& rTrialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));

				//	memory for shapes
				std::vector<number> vShape;

				//	loop ips
				for(size_t ip = 0; ip < nip; ++ip)
					{
					//	evaluate shapes at ip
						rTrialSpace.shapes(vShape, vLocIP[ip]);

					// 	compute solution at integration point
						vValue[ip] = 0.0;
						for(size_t sh = 0; sh < noc; ++sh)
						{
							 vValue[ip] += vShape[sh] * m_acTurbulentViscosity[element->vertex(sh)];
						}
					}

			}
			UG_CATCH_THROW("TurbulentViscosityData: trial space missing, Reference Object: "
			         <<roid<<", Trial Space: LAGRANGE 1, refDim="<<refDim);

			number kinViscValues[max_number_of_ips];
			(*m_imKinViscosity)(kinViscValues,
							vGlobIP,
		                    time, si,
		                    elem,
		                    vCornerCoords,
		                    vLocIP,
		                    nip,
		                    u,
		                    vJT);
			for (size_t ip=0;ip < nip;ip++){
				// UG_LOG("turbVis(" << ip << ")=" << vValue[ip] << "+" << kinViscValues[ip] << "\n");
				vValue[ip] += kinViscValues[ip];
			}
		}
		
		static const size_t max_number_of_ips = 20;

		bool m_spaceFilter;
		bool m_timeFilter;
		number m_timeFilterEps;

		void update();

		void set_space_filter(bool b){
			m_spaceFilter=b;
		}
		void set_time_filter(bool b){
			m_timeFilter=b;
			m_timeFilterEps=0.001;
		}
		void set_time_filter_eps(number eps){
			if (eps!=1)
				m_timeFilter=true;
			else
				m_timeFilter=false;
			m_timeFilterEps=eps;
		}
};


template <typename TGridFunction>
const number FV1DynamicTurbViscData<TGridFunction>::m_small = 1e-8;

} // namespace NavierStokes
} // end namespace ug

// include implementation
#include "turbulent_viscosity_fv1_impl.h"

#endif /* __H__UG__NAVIER_STOKES__INCOMPRESSIBLE__FV1__TURBULENT_VISCOSITY_DATA_FV1__ */
