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

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__SHEAR_STRESS__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__SHEAR_STRESS__

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
#include "../../properties_interface.h"
#include "common/math/ugmath.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/util/attachment_operations.hpp"
#endif

namespace ug{
namespace NavierStokes{



template<int dim>
void computeElemBarycenter(MathVector<dim>& bary,int NumSh , const MathVector<dim> Vec[]){
	bary = 0;
	for (size_t i=0;i<NumSh;i++){
		bary+=Vec[i];
	}
	bary/=NumSh;
}

/**
concept derived from grid_function_user_data.h
 */
template <typename TGridFunction>
class ShearStressFV1
:     public StdUserData<ShearStressFV1<TGridFunction>, number, TGridFunction::dim>,
      virtual public INewtonUpdate
{
    ///    domain type
    typedef typename TGridFunction::domain_type domain_type;

    ///    algebra type
    typedef typename TGridFunction::algebra_type algebra_type;

    /// position accessor type
    typedef typename domain_type::position_accessor_type position_accessor_type;

    ///    world dimension
    static const int dim = domain_type::dim;

    ///    grid type
    typedef typename domain_type::grid_type grid_type;

    /// element type
    typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

    /// MathVector<dim> attachment
    //        typedef MathVector<dim> vecDim;
    //        typedef Attachment<vecDim> AMathVectorDim;

    /// attachment accessor
    typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;

    /// element iterator
    typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

    /// vertex iterator
    typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

          private:

    //    ShearRate attachment accessor (interpolated ShearRate in vertices)
    ANumber m_aSR;
    aVertexNumber m_shear_rate;

    //  volume attachment accessor
    ANumber m_aVol;
    aVertexNumber m_vol;

    // level set grid function
    SmartPtr<TGridFunction> m_u;

    //    approximation space for level and surface grid
    SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

    //  grid
    grid_type* m_grid;

          private:

    ///    Data import for source
    SmartPtr<CplUserData<MathVector<dim>,dim> > m_imSource;

          public:
    /////////// Source

    void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
    {
        m_imSource = data;
    }

    void set_source(number f_x)
    {
        SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
        for (int i=0;i<dim;i++){
            f->set_entry(i, f_x);
        }
        set_source(f);
    }

    void set_source(number f_x, number f_y)
    {
        if (dim!=2){
            UG_THROW("NavierStokes: Setting source vector of dimension 2"
                    " to a Discretization for world dim " << dim);
        } else {
            SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
            f->set_entry(0, f_x);
            f->set_entry(1, f_y);
            set_source(f);
        }
    }

    void set_source(number f_x, number f_y, number f_z)
    {
        if (dim<3){
            UG_THROW("NavierStokes: Setting source vector of dimension 3"
                    " to a Discretization for world dim " << dim);
        }
        else
        {
            SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
            f->set_entry(0, f_x);
            f->set_entry(1, f_y);
            f->set_entry(2, f_z);
            set_source(f);
        }
    }

#ifdef UG_FOR_LUA
    void set_source(const char* fctName)
    {
        set_source(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
    }
#endif

          public:
    /// constructor
    ShearStressFV1(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
        
		if (spGridFct->num_fct() != dim+2)
			UG_THROW("NavierStokesMultiphase: Need exactly "<<dim+2<<" functions");
		for (int d=0;d<dim+2;d++){
			if (spGridFct->local_finite_element_id(d) != LFEID(LFEID::LAGRANGE, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Lagrange P1 type.");
			}
		}

		m_u = spGridFct;
        domain_type& domain = *m_u->domain().get();
        grid_type& grid = *domain.grid();
        m_grid = &grid;
        m_spApproxSpace = approxSpace;
        set_source(0.0);
        grid.template attach_to<Vertex>(m_aSR);
        grid.template attach_to<Vertex>(m_aVol);
        m_shear_rate.access(grid,m_aSR);
        m_vol.access(grid,m_aVol);
        // set all values to zero
        SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
        SetAttachmentValues(m_shear_rate, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
        //this->update();
    }

    virtual ~ShearStressFV1(){};

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
        UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
        elem_type* element = static_cast<elem_type*>(elem);

        //    reference object id
        ReferenceObjectID roid = elem->reference_object_id();

        const size_t numVertices = element->num_vertices();
        //    get domain of grid function
        const domain_type& domain = *m_u->domain().get();

        //    get position accessor
        typedef typename domain_type::position_accessor_type position_accessor_type;
        const position_accessor_type& posAcc = domain.position_accessor();

//        position_accessor_type aaPos = m_u->domain()->position_accessor();

        // coord and vertex array
        MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
        Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
        DimFV1Geometry<dim> geo;

        for(size_t i = 0; i < numVertices; ++i){
            vVrt[i] = element->vertex(i);
            coCoord[i] = posAcc[vVrt[i]];
        };

        // evaluate finite volume geometry
        geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

		// Lagrange 1 trial space
		const LocalShapeFunctionSet<refDim>& rTrialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));

        std::vector<number> shapes;
        for (size_t ip=0;ip<nip;ip++)
        {
            number value = 0.0;
			rTrialSpace.shapes(shapes,vLocIP[ip]);
            for (size_t sh=0;sh<numVertices;sh++)
                value += m_shear_rate[vVrt[sh]]*shapes[sh];
            
            vValue[ip] = value;
            
        }
        
        
            
    }; // evaluate

    void update(){
        //    get domain
		UG_LOG("Updating Shear Rate... \n");
        domain_type& domain = *m_u->domain().get();
        //    create Multiindex
        std::vector<DoFIndex> multInd;
        DimFV1Geometry<dim> geo;
        //    coord and vertex array
        MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
        MathVector<dim> coGrad[domain_traits<dim>::MaxNumVerticesOfElem];
        Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

        //    get position accessor
        typedef typename domain_type::position_accessor_type position_accessor_type;
        const position_accessor_type& posAcc = domain.position_accessor();

        // set volume and p values to zero
        SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
        SetAttachmentValues(m_shear_rate, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
        // compute pressure in vertices by averaging
        for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
            ElemIterator iter = m_u->template begin<elem_type>(si);
            ElemIterator iterEnd = m_u->template end<elem_type>(si);
            for(  ;iter !=iterEnd; ++iter)
            {
                elem_type* elem = *iter;
                const size_t numVertices = elem->num_vertices();
                for(size_t i = 0; i < numVertices; ++i){
                    vVrt[i] = elem->vertex(i);
                    coCoord[i] = posAcc[vVrt[i]];
                };
                geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
                for(size_t i = 0; i < numVertices; ++i)
				{
                    number scvVol = geo.scv(i).volume();
                    m_vol[vVrt[i]]+=scvVol;
                    
                    MathMatrix<dim,dim> VelGrad; MatSet(VelGrad,0.0);
                    
                    //    sum up contributions of each shape
                    for(size_t sh = 0; sh < numVertices; ++sh)
                    {
                        //  Loop dimensions for derivative
                        for(int d1 = 0; d1 <dim; ++d1)
                        {
                            m_u->dof_indices(elem->vertex(sh), d1, multInd);
                            //    read value of index from vector
                            number uVal = DoFRef(*m_u,multInd[0]);
                        //  Loop dimensions for direction
                            for(int d2 = 0; d2 < dim; ++d2)
                            {
                                VelGrad(d1, d2) += uVal*geo.scv(i).global_grad(sh)[d2];
                            }
                        }
                    }
                    number gamma=0.0;
                    // compute inner sum
                    for(int d1 = 0; d1 < dim; ++d1)
                    {
                        for(int d2 = 0; d2 < dim; ++d2)
                        {
                            gamma += pow((VelGrad(d1,d2) + VelGrad(d2,d1)),2);
                        }
                    }
                    
                    gamma =sqrt((0.5*gamma));
                    
                    
                    m_shear_rate[vVrt[i]] += gamma * scvVol;
                    
                }
            }
        }
		
		#ifdef UG_PARALLEL
			AttachmentAllReduce<Vertex> (*domain.grid(), m_aVol, PCL_RO_SUM);
			AttachmentAllReduce<Vertex> (*domain.grid(), m_aSR, PCL_RO_SUM);
		#endif
		
        PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
        // go over all vertices and average
        for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
            VertexIterator iter = m_u->template begin<Vertex>(si);
            VertexIterator iterEnd = m_u->template end<Vertex>(si);
            for(  ;iter !=iterEnd; ++iter)
            {
                Vertex* vrt = *iter;
                if (pbm && pbm->is_slave(vrt)) continue;
                    m_shear_rate[vrt] /= m_vol[vrt];
            }
        }
    }

          private:
    static const size_t max_number_of_ips = 20;

          public:
    virtual void operator() (number& value,
                             const MathVector<dim>& globIP,
                             number time, int si) const
    {
        UG_THROW("ShearStressUserData: Need element.");
    }

    virtual void operator() (number vValue[],
                             const MathVector<dim> vGlobIP[],
                             number time, int si, const size_t nip) const
    {
        UG_THROW("ShearStress: Need element.");
    }

    virtual void compute(LocalVector* u, GridObject* elem,
                         const MathVector<dim> vCornerCoords[], bool bDeriv = false)
    {
        const int si = this->subset();
        for(size_t s = 0; s < this->num_series(); ++s)
            evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
                          elem, vCornerCoords, this->template local_ips<dim>(s),
                          this->num_ip(s), u);
    }

    virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
                         const MathVector<dim> vCornerCoords[], bool bDeriv = false)
    {
        const int si = this->subset();
        for(size_t s = 0; s < this->num_series(); ++s)
            evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
                          elem, vCornerCoords, this->template local_ips<dim>(s),
                          this->num_ip(s), &(u->solution(this->time_point(s))));
    }

    ///    returns if provided data is continuous over geometric object boundaries
    virtual bool continuous() const {return false;}

    ///    returns if grid function is needed for evaluation
    virtual bool requires_grid_fct() const {return true;}
};

/**
concept derived from grid_function_user_data.h
 */
template <typename TGridFunction>
class ParticlePressureFV1
:     public StdUserData<ParticlePressureFV1<TGridFunction>, number, TGridFunction::dim>,
	  virtual public INewtonUpdate
{
	///    domain type
	typedef typename TGridFunction::domain_type domain_type;

	///    algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	///    world dimension
	static const int dim = domain_type::dim;
	///    Pressure
	static const int _P_ = domain_type::dim;
	///    Pressure
	static const int _C_ = domain_type::dim+1;

	///    grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// MathVector<dim> attachment
	//        typedef MathVector<dim> vecDim;
	//        typedef Attachment<vecDim> AMathVectorDim;

	/// attachment accessor
	typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

		  private:

	//    ShearRate attachment accessor (interpolated ShearRate in vertices)
	ANumber m_aSR;
	aVertexNumber m_shear_rate;

	//  volume attachment accessor
	ANumber m_aVol;
	aVertexNumber m_vol;

	// level set grid function
	SmartPtr<TGridFunction> m_u;

	//    approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;

		  private:

	///    Data import for source
	SmartPtr<CplUserData<MathVector<dim>,dim> > m_imSource;
	Interface<dim>* Inter;

		  public:
	/////////// Source

	void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
	{
		m_imSource = data;
	}

	void set_source(number f_x)
	{
		SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
		for (int i=0;i<dim;i++){
			f->set_entry(i, f_x);
		}
		set_source(f);
	}

	void set_source(number f_x, number f_y)
	{
		if (dim!=2){
			UG_THROW("NavierStokes: Setting source vector of dimension 2"
					" to a Discretization for world dim " << dim);
		} else {
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			set_source(f);
		}
	}

	void set_source(number f_x, number f_y, number f_z)
	{
		if (dim<3){
			UG_THROW("NavierStokes: Setting source vector of dimension 3"
					" to a Discretization for world dim " << dim);
		}
		else
		{
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			f->set_entry(2, f_z);
			set_source(f);
		}
	}

#ifdef UG_FOR_LUA
	void set_source(const char* fctName)
	{
		set_source(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
	}
#endif
	  void set_phase_parameters(Interface<dim>* user)
	  {
		  if (!user) UG_THROW("Interface pointer is null!");
		  if (!user->valid())
			  UG_THROW("Interface parameters has not been initialized");
		  Inter = user;
	  }

		  public:
	/// constructor
	ParticlePressureFV1(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
		
		if (spGridFct->num_fct() != dim+2)
			UG_THROW("NavierStokesMultiphase: Need exactly "<<dim+2<<" functions");
		for (int d=0;d<dim+2;d++)
		{
			if (spGridFct->local_finite_element_id(d) != LFEID(LFEID::LAGRANGE, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Lagrange P1 type.");
			}
		}
		m_u = spGridFct;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		set_source(0.0);
		grid.template attach_to<Vertex>(m_aSR);
		grid.template attach_to<Vertex>(m_aVol);
		m_shear_rate.access(grid,m_aSR);
		m_vol.access(grid,m_aVol);
		// set all values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_shear_rate, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//this->update();
	}

	virtual ~ParticlePressureFV1(){};

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
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);

		//    reference object id
		ReferenceObjectID roid = elem->reference_object_id();

		const size_t numVertices = element->num_vertices();
		const size_t MaxVertices = domain_traits<dim>::MaxNumVerticesOfElem;
		//    get domain of grid function
		const domain_type& domain = *m_u->domain().get();

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

//        position_accessor_type aaPos = m_u->domain()->position_accessor();

		// coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
		DimFV1Geometry<dim> geo;

		for(size_t i = 0; i < numVertices; ++i){
			vVrt[i] = element->vertex(i);
			coCoord[i] = posAcc[vVrt[i]];
		};

		// evaluate finite volume geometry
		geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

		// Lagrange 1 trial space
		const LocalShapeFunctionSet<refDim>& rTrialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));

		std::vector<number> shapes;
		number Ps[numVertices];
		number Gamma[numVertices];
		for (size_t sh=0;sh<numVertices;sh++)
			Gamma[sh] = m_shear_rate[vVrt[sh]];
		
		//Inter->Ps( Ps, NULL, Gamma, *u, _C_, numVertices, false);
		
		for (size_t ip=0;ip<nip;ip++)
		{
			number value = 0.0;

			rTrialSpace.shapes(shapes,vLocIP[ip]);
			
			
			for (size_t sh=0;sh<numVertices;sh++)
				value += Gamma[sh]*shapes[sh];
			
			vValue[ip] = value;
			
		}
		
		
			
	}; // evaluate

	void update(){
		//    get domain
		UG_LOG("Updating Particle Pressure Ps... \n");
		
		domain_type& domain = *m_u->domain().get();
		//    create Multiindex
		std::vector<DoFIndex> multInd;
		DimFV1Geometry<dim> geo;
		//    coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		MathVector<dim> coGrad[domain_traits<dim>::MaxNumVerticesOfElem];
		Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

		// set volume and p values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_shear_rate, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		// compute pressure in vertices by averaging
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			ElemIterator iter = m_u->template begin<elem_type>(si);
			ElemIterator iterEnd = m_u->template end<elem_type>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				elem_type* elem = *iter;
				const size_t numVertices = elem->num_vertices();
				for(size_t i = 0; i < numVertices; ++i){
					vVrt[i] = elem->vertex(i);
					coCoord[i] = posAcc[vVrt[i]];
				};
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				for(size_t i = 0; i < numVertices; ++i){
					number scvVol = geo.scv(i).volume();
					m_vol[vVrt[i]]+=scvVol;
					
					MathMatrix<dim,dim> VelGrad; MatSet(VelGrad,0.0);
					
					//    sum up contributions of each shape
					for(size_t sh = 0; sh < numVertices; ++sh)
					{
						//  Loop dimensions for derivative
						for(int d1 = 0; d1 <dim; ++d1)
						{
							m_u->dof_indices(elem->vertex(sh), d1, multInd);
							//    read value of index from vector
							number uVal = DoFRef(*m_u,multInd[0]);
						//  Loop dimensions for direction
							for(int d2 = 0; d2 < dim; ++d2)
							{
								VelGrad(d1, d2) += uVal*geo.scv(i).global_grad(sh)[d2];
							}
						}
					}
					number gamma=0.0;
					// compute inner sum
					for(int d1 = 0; d1 < dim; ++d1)
					{
						for(int d2 = 0; d2 < dim; ++d2)
						{
							gamma += pow((VelGrad(d1,d2) + VelGrad(d2,d1)),2);
						}
					}
					
					gamma =sqrt((0.5*gamma));
					
					
					m_shear_rate[vVrt[i]] += gamma * scvVol;
					
				}
			}
		}
		PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
		// go over all vertices and average
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			VertexIterator iter = m_u->template begin<Vertex>(si);
			VertexIterator iterEnd = m_u->template end<Vertex>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				Vertex* vrt = *iter;
				if (pbm && pbm->is_slave(vrt)) continue;
					m_shear_rate[vrt] /= m_vol[vrt];
			}
		}
	}

		  private:
	static const size_t max_number_of_ips = 20;

		  public:
	virtual void operator() (number& value,
							 const MathVector<dim>& globIP,
							 number time, int si) const
	{
		UG_THROW("LevelSetUserData: Need element.");
	}

	virtual void operator() (number vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si, const size_t nip) const
	{
		UG_THROW("LevelSetUserData: Need element.");
	}

	virtual void compute(LocalVector* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), &(u->solution(this->time_point(s))));
	}

	///    returns if provided data is continuous over geometric object boundaries
	virtual bool continuous() const {return false;}

	///    returns if grid function is needed for evaluation
	virtual bool requires_grid_fct() const {return true;}
};


/**
concept derived from grid_function_user_data.h
 */
template <typename TGridFunction>
class PressureGradientMean
:     public StdUserData<PressureGradientMean<TGridFunction>, MathVector<TGridFunction::dim>, TGridFunction::dim>,
	  virtual public INewtonUpdate
{
	///    domain type
	typedef typename TGridFunction::domain_type domain_type;

	///    algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	///    world dimension
	static const int dim = domain_type::dim;
	///    Pressure
	static const int _P_ = domain_type::dim;
	///    Pressure
	static const int _C_ = domain_type::dim+1;

	///    grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// MathVector<dim> attachment
	typedef MathVector<dim> vecDim;
	typedef Attachment<vecDim> AMathVectorDim;

	/// attachment accessor
	typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;
	typedef PeriodicAttachmentAccessor<Vertex,AMathVectorDim > aVertexDimVector;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

private:

	//    Normal attachment accessor (average normal in vertices)
	AMathVectorDim m_aGradient;
	aVertexDimVector m_gradient;
	
	//    Normal attachment accessor (average normal in vertices)
	//AMathVectorDim m_aTang;
	//aVertexDimVector m_tang;

	//  volume attachment accessor
	ANumber m_aVol;
	aVertexNumber m_vol;

	// level set grid function
	SmartPtr<TGridFunction> m_u;

	//    approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;
	
	number m_limit = 1e-03;
	number m_limit_grad = 1e-03;
	number m_theta_cr = 34.0*3.1416/180.0;

private:

	///    Data import for source
	SmartPtr<CplUserData<MathVector<dim>,dim> > m_imSource;
	Interface<dim>* Inter;

		  public:
	void set_theta(number data)
	{
		m_theta_cr = data*3.1416/180.0;
	}
	void set_gradient_limit(number data)
	{
		m_limit = data;
		m_limit_grad = 0.001 * data;
	}
	
	void set_phase_parameters(Interface<dim>* user)
	{
		if (!user) UG_THROW("Interface pointer is null!");
		if (!user->valid())
			UG_THROW("Interface parameters has not been initialized");
		Inter = user;
	}

public:
	/// constructor
	PressureGradientMean(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
		
		if (spGridFct->num_fct() != dim+2)
			UG_THROW("NavierStokesMultiphase: Need exactly "<<dim+2<<" functions");
		for (int d=0;d<dim+2;d++)
		{
			if (spGridFct->local_finite_element_id(d) != LFEID(LFEID::LAGRANGE, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Lagrange P1 type.");
			}
		}
		m_u = spGridFct;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		grid.template attach_to<Vertex>(m_aGradient);
		//grid.template attach_to<Vertex>(m_aTang);
		grid.template attach_to<Vertex>(m_aVol);
		m_gradient.access(grid,m_aGradient);
		//m_tang.access(grid,m_aTang);
		m_vol.access(grid,m_aVol);
		// set all values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_gradient, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//SetAttachmentValues(m_tang, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//this->update();
	}

	virtual ~PressureGradientMean(){};

	template <int refDim>
	inline void evaluate(MathVector<dim> vValue[],
						 const MathVector<dim> vGlobIP[],
						 number time, int si,
						 GridObject* elem,
						 const MathVector<dim> vCornerCoords[],
						 const MathVector<refDim> vLocIP[],
						 const size_t nip,
						 LocalVector* u,
						 const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);

		//    reference object id
		ReferenceObjectID roid = elem->reference_object_id();

		const size_t numVertices = element->num_vertices();
		const size_t MaxVertices = domain_traits<dim>::MaxNumVerticesOfElem;
		//    get domain of grid function
		const domain_type& domain = *m_u->domain().get();

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

//        position_accessor_type aaPos = m_u->domain()->position_accessor();

		// coord and vertex array
		MathVector<dim> coCoord[numVertices];
		Vertex* vVrt[numVertices];
		DimFV1Geometry<dim> geo;

		for(size_t i = 0; i < numVertices; ++i){
			vVrt[i] = element->vertex(i);
			coCoord[i] = posAcc[vVrt[i]];
		};

		// evaluate finite volume geometry
		geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

		// Lagrange 1 trial space
		const LocalShapeFunctionSet<refDim>& rTrialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));

		std::vector<number> shapes;
		

	//    storage for shape function at ip
		MathVector<refDim> vLocGrad[numVertices];
		MathVector<refDim> locGrad;

	//    Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		
		DimReferenceMapping<refDim, dim>& mapping = ReferenceMappingProvider::get<refDim, dim>(roid, coCoord);
		
		
		for (size_t ip=0;ip<nip;ip++)
		{
			MathVector<dim> gradientP = 0.0;
			
			//MathVector<dim> GradC = 0.0;
			rTrialSpace.shapes(shapes,vLocIP[ip]);
			
			
			for (size_t sh=0;sh<numVertices;sh++)
				for(int d = 0; d < refDim; ++d)
				{
					gradientP[d] += m_gradient[vVrt[sh]][d]*shapes[sh];
				}
			
			vValue[ip] = gradientP;
		}
		
		
			
	}; // evaluate

	void update(){
		//    get domain
		UG_LOG("Updating PressureGradientMean... \n");
		domain_type& domain = *m_u->domain().get();
		//    create Multiindex
		std::vector<DoFIndex> multInd;
		DimFV1Geometry<dim> geo;
		//    coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		MathVector<dim> coGrad[domain_traits<dim>::MaxNumVerticesOfElem];
		Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

		// set volume, tang and normal values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_gradient, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//SetAttachmentValues(m_tang, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		// compute pressure in vertices by averaging
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			ElemIterator iter = m_u->template begin<elem_type>(si);
			ElemIterator iterEnd = m_u->template end<elem_type>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				elem_type* elem = *iter;
				const size_t numVertices = elem->num_vertices();
				for(size_t i = 0; i < numVertices; ++i){
					vVrt[i] = elem->vertex(i);
					coCoord[i] = posAcc[vVrt[i]];
				};
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				for(size_t i = 0; i < numVertices; ++i){
					number scvVol = geo.scv(i).volume();
					
					MathVector<dim> GradP; VecSet(GradP,0.0);
					
					//    sum up contributions of each shape
					for(size_t sh = 0; sh < numVertices; ++sh)
					{
						m_u->dof_indices(elem->vertex(sh), _P_, multInd);
						//    read value of index from vector
						number uVal = DoFRef(*m_u,multInd[0]);
						//uVal = fmax(uVal, 0.0);
						//  Loop dimensions for derivative
						for(int d1 = 0; d1 <dim; ++d1)
						{
							GradP[d1] += uVal*geo.scv(i).global_grad(sh)[d1];
						}
					}

					m_vol[vVrt[i]]+=scvVol;
					for(int d1 = 0; d1 <dim; ++d1)
						m_gradient[vVrt[i]][d1] += GradP[d1] * scvVol;
					
					
				}
			}
		}
		
		#ifdef UG_PARALLEL
			AttachmentAllReduce<Vertex> (*domain.grid(), m_aVol, PCL_RO_SUM);
			AttachmentAllReduce<Vertex> (*domain.grid(), m_aGradient, PCL_RO_SUM);
		#endif
		
		PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
		// go over all vertices and average
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			VertexIterator iter = m_u->template begin<Vertex>(si);
			VertexIterator iterEnd = m_u->template end<Vertex>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				Vertex* vrt = *iter;
				if (pbm && pbm->is_slave(vrt)) continue;

				for(int d1 = 0; d1 <dim; ++d1)
					m_gradient[vrt][d1] /= m_vol[vrt];

				/*MathVector<dim> tang; VecSet(tang,0.0);
				if(dim == 2)
				{
					number ss;
					if(fabs(m_normal[vrt][0])<m_limit)
						ss = 1.0;
					else
						ss = (m_normal[vrt][0]*m_normal[vrt][1]>0.0)? 1.0 : -1.0;
					
					tang[0] =   ss * m_normal[vrt][1];
					tang[1] = - fabs(m_normal[vrt][0]);
				}
				else if (dim == 3)
				{
					tang[0] = -m_normal[vrt][0]*m_normal[vrt][1];
					tang[1] = -m_normal[vrt][1]*m_normal[vrt][2];
					tang[2] =  m_normal[vrt][0]*m_normal[vrt][0] + m_normal[vrt][1]*m_normal[vrt][1];
					UG_THROW("MeanPressureGradient: Works only for Dim = 2");
				}
				else UG_THROW("MeanPressureGradient: Works only for Dim = 2 , 3");

				//VecScale(tang, tang,1.0/VecTwoNorm(tang));
				m_tang[vrt] = tang;*/
				
			}
		}
	}

private:
	static const size_t max_number_of_ips = 20;

public:
	virtual void operator() (MathVector<dim>& value,
							 const MathVector<dim>& globIP,
							 number time, int si) const
	{
		UG_THROW("PressureGradientMean: Need element.");
	}

	virtual void operator() (MathVector<dim> vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si, const size_t nip) const
	{
		UG_THROW("PressureGradientMean: Need element.");
	}

	virtual void compute(LocalVector* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), &(u->solution(this->time_point(s))));
	}

	///    returns if provided data is continuous over geometric object boundaries
	virtual bool continuous() const {return false;}

	///    returns if grid function is needed for evaluation
	virtual bool requires_grid_fct() const {return true;}
};

/**
concept derived from grid_function_user_data.h
 */
template <typename TGridFunction>
class DuneNormal
:     public StdUserData<DuneNormal<TGridFunction>, MathVector<TGridFunction::dim>, TGridFunction::dim>,
	  virtual public INewtonUpdate
{
	///    domain type
	typedef typename TGridFunction::domain_type domain_type;

	///    algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	///    world dimension
	static const int dim = domain_type::dim;
	///    Pressure
	static const int _P_ = domain_type::dim;
	///    Pressure
	static const int _C_ = domain_type::dim+1;

	///    grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// MathVector<dim> attachment
	typedef MathVector<dim> vecDim;
	typedef Attachment<vecDim> AMathVectorDim;

	/// attachment accessor
	typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;
	typedef PeriodicAttachmentAccessor<Vertex,AMathVectorDim > aVertexDimVector;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

private:

	//    Normal attachment accessor (average normal in vertices)
	AMathVectorDim m_aNormal;
	aVertexDimVector m_normal;
	
	//    Normal attachment accessor (average normal in vertices)
	//AMathVectorDim m_aTang;
	//aVertexDimVector m_tang;

	//  volume attachment accessor
	ANumber m_aVol;
	aVertexNumber m_vol;

	// level set grid function
	SmartPtr<TGridFunction> m_u;

	//    approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;
	
	number m_limit = 1e-03;
	number m_limit_grad = 1e-03;
	number m_theta_cr = 34.0*3.1416/180.0;

private:

	///    Data import for source
	SmartPtr<CplUserData<MathVector<dim>,dim> > m_imSource;
	Interface<dim>* Inter;

		  public:
	void set_theta(number data)
	{
		m_theta_cr = data*3.1416/180.0;
	}
	void set_gradient_limit(number data)
	{
		m_limit = data;
		m_limit_grad = 0.001 * data;
	}
	
	void set_phase_parameters(Interface<dim>* user)
	{
		if (!user) UG_THROW("Interface pointer is null!");
		if (!user->valid())
			UG_THROW("Interface parameters has not been initialized");
		Inter = user;
	}

public:
	/// constructor
	DuneNormal(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
		
		if (spGridFct->num_fct() != dim+2)
			UG_THROW("NavierStokesMultiphase: Need exactly "<<dim+2<<" functions");
		for (int d=0;d<dim+2;d++)
		{
			if (spGridFct->local_finite_element_id(d) != LFEID(LFEID::LAGRANGE, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Lagrange P1 type.");
			}
		}
		m_u = spGridFct;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		grid.template attach_to<Vertex>(m_aNormal);
		//grid.template attach_to<Vertex>(m_aTang);
		grid.template attach_to<Vertex>(m_aVol);
		m_normal.access(grid,m_aNormal);
		//m_tang.access(grid,m_aTang);
		m_vol.access(grid,m_aVol);
		// set all values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_normal, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//SetAttachmentValues(m_tang, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//this->update();
	}

	virtual ~DuneNormal(){};

	template <int refDim>
	inline void evaluate(MathVector<dim> vValue[],
						 const MathVector<dim> vGlobIP[],
						 number time, int si,
						 GridObject* elem,
						 const MathVector<dim> vCornerCoords[],
						 const MathVector<refDim> vLocIP[],
						 const size_t nip,
						 LocalVector* u,
						 const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);

		//    reference object id
		ReferenceObjectID roid = elem->reference_object_id();

		const size_t numVertices = element->num_vertices();
		const size_t MaxVertices = domain_traits<dim>::MaxNumVerticesOfElem;
		//    get domain of grid function
		const domain_type& domain = *m_u->domain().get();

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

//        position_accessor_type aaPos = m_u->domain()->position_accessor();

		// coord and vertex array
		MathVector<dim> coCoord[numVertices];
		Vertex* vVrt[numVertices];
		DimFV1Geometry<dim> geo;

		for(size_t i = 0; i < numVertices; ++i){
			vVrt[i] = element->vertex(i);
			coCoord[i] = posAcc[vVrt[i]];
		};

		// evaluate finite volume geometry
		geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

		// Lagrange 1 trial space
		const LocalShapeFunctionSet<refDim>& rTrialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));

		std::vector<number> shapes;
		

	//    storage for shape function at ip
		MathVector<refDim> vLocGrad[numVertices];
		MathVector<refDim> locGrad;

	//    Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		
		DimReferenceMapping<refDim, dim>& mapping = ReferenceMappingProvider::get<refDim, dim>(roid, coCoord);
		
		
		for (size_t ip=0;ip<nip;ip++)
		{
			MathVector<dim> normal = 0.0;
			
			//MathVector<dim> GradC = 0.0;
			rTrialSpace.shapes(shapes,vLocIP[ip]);
			
			
			for (size_t sh=0;sh<numVertices;sh++)
				for(int d = 0; d < refDim; ++d)
				{
					normal[d] += m_normal[vVrt[sh]][d]*shapes[sh];
				}
			number normal_mag =  VecTwoNorm(normal);

			VecScale(normal,normal,1.0/normal_mag);
			
			vValue[ip] = normal;
		}
		
		
			
	}; // evaluate

	void update(){
		//    get domain
		UG_LOG("Updating normal... \n");
		domain_type& domain = *m_u->domain().get();
		//    create Multiindex
		std::vector<DoFIndex> multInd;
		DimFV1Geometry<dim> geo;
		//    coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		MathVector<dim> coGrad[domain_traits<dim>::MaxNumVerticesOfElem];
		Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

		// set volume, tang and normal values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_normal, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//SetAttachmentValues(m_tang, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		// compute pressure in vertices by averaging
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			ElemIterator iter = m_u->template begin<elem_type>(si);
			ElemIterator iterEnd = m_u->template end<elem_type>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				elem_type* elem = *iter;
				const size_t numVertices = elem->num_vertices();
				for(size_t i = 0; i < numVertices; ++i){
					vVrt[i] = elem->vertex(i);
					coCoord[i] = posAcc[vVrt[i]];
				};
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				for(size_t i = 0; i < numVertices; ++i){
					number scvVol = geo.scv(i).volume();
					
					MathVector<dim> GradC; VecSet(GradC,0.0);
					MathVector<dim> Normal; VecSet(Normal,0.0);
					MathVector<dim> nn; VecSet(Normal,0.0); nn[dim-1] = -1.0;
					
					//    sum up contributions of each shape
					for(size_t sh = 0; sh < numVertices; ++sh)
					{
						m_u->dof_indices(elem->vertex(sh), _C_, multInd);
						//    read value of index from vector
						number uVal = DoFRef(*m_u,multInd[0]);
						//uVal = fmax(uVal, 0.0);
						//  Loop dimensions for derivative
						for(int d1 = 0; d1 <dim; ++d1)
						{
							GradC[d1] += uVal*geo.scv(i).global_grad(sh)[d1];
						}
					}
					const number grad_c_mag = VecTwoNorm(GradC);
					const number eps = 5e-01;
					const number p = 2.0;
					const number alpha = fmin(pow(grad_c_mag,p)/(pow(grad_c_mag,p) + eps),1.0 );
					VecScaleAdd(GradC,alpha, GradC, 1.0-alpha,nn);
					const number GRAD_c_mag = VecTwoNorm(GradC);
					/*if(grad_c_mag<=m_limit_grad)
					{
						VecSet(Normal,0.0);
					}
					else
					{
						VecScale(Normal,GradC,-1.0/grad_c_mag);
						m_vol[vVrt[i]]+=scvVol;
					}*/
					VecScale(Normal,GradC,-1.0/GRAD_c_mag);
					m_vol[vVrt[i]]+=scvVol;
					
					for(int d1 = 0; d1 <dim; ++d1)
						m_normal[vVrt[i]][d1] += Normal[d1] * scvVol;
					
					
				}
			}
		}
		
		#ifdef UG_PARALLEL
			AttachmentAllReduce<Vertex> (*domain.grid(), m_aVol, PCL_RO_SUM);
			AttachmentAllReduce<Vertex> (*domain.grid(), m_aNormal, PCL_RO_SUM);
		#endif
		
		PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
		// go over all vertices and average
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			VertexIterator iter = m_u->template begin<Vertex>(si);
			VertexIterator iterEnd = m_u->template end<Vertex>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				Vertex* vrt = *iter;
				if (pbm && pbm->is_slave(vrt)) continue;

				for(int d1 = 0; d1 <dim; ++d1)
					m_normal[vrt][d1] /= m_vol[vrt];

				/*MathVector<dim> tang; VecSet(tang,0.0);
				if(dim == 2)
				{
					number ss;
					if(fabs(m_normal[vrt][0])<m_limit)
						ss = 1.0;
					else
						ss = (m_normal[vrt][0]*m_normal[vrt][1]>0.0)? 1.0 : -1.0;
					
					tang[0] =   ss * m_normal[vrt][1];
					tang[1] = - fabs(m_normal[vrt][0]);
				}
				else if (dim == 3)
				{
					tang[0] = -m_normal[vrt][0]*m_normal[vrt][1];
					tang[1] = -m_normal[vrt][1]*m_normal[vrt][2];
					tang[2] =  m_normal[vrt][0]*m_normal[vrt][0] + m_normal[vrt][1]*m_normal[vrt][1];
					UG_THROW("DuneNormal: Works only for Dim = 2");
				}
				else UG_THROW("DuneNormal: Works only for Dim = 2 , 3");

				//VecScale(tang, tang,1.0/VecTwoNorm(tang));
				m_tang[vrt] = tang;*/
				
			}
		}
	}

private:
	static const size_t max_number_of_ips = 20;

public:
	virtual void operator() (MathVector<dim>& value,
							 const MathVector<dim>& globIP,
							 number time, int si) const
	{
		UG_THROW("Normal: Need element.");
	}

	virtual void operator() (MathVector<dim> vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si, const size_t nip) const
	{
		UG_THROW("Normal: Need element.");
	}

	virtual void compute(LocalVector* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), &(u->solution(this->time_point(s))));
	}

	///    returns if provided data is continuous over geometric object boundaries
	virtual bool continuous() const {return false;}

	///    returns if grid function is needed for evaluation
	virtual bool requires_grid_fct() const {return true;}
};


/**
concept derived from grid_function_user_data.h
 */
template <typename TGridFunction>
class SlipVelocity
:     public StdUserData<SlipVelocity<TGridFunction>, MathVector<TGridFunction::dim>, TGridFunction::dim>,
	  virtual public INewtonUpdate
{
	///    domain type
	typedef typename TGridFunction::domain_type domain_type;

	///    algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	///    world dimension
	static const int dim = domain_type::dim;
	///    Pressure
	static const int _P_ = domain_type::dim;
	///    Pressure
	static const int _C_ = domain_type::dim+1;

	///    grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// MathVector<dim> attachment
	typedef MathVector<dim> vecDim;
	typedef Attachment<vecDim> AMathVectorDim;

	/// attachment accessor
	typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;
	typedef PeriodicAttachmentAccessor<Vertex,AMathVectorDim > aVertexDimVector;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

private:

	//    Normal attachment accessor (average normal in vertices)
	AMathVectorDim m_aNormal;
	aVertexDimVector m_normal;
	
	//    Normal attachment accessor (average normal in vertices)
	//AMathVectorDim m_aTang;
	//aVertexDimVector m_tang;

	//  volume attachment accessor
	ANumber m_aVol;
	aVertexNumber m_vol;

	// level set grid function
	SmartPtr<TGridFunction> m_u;

	//    approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;
	
	number m_limit = 1e-03;
	number m_limit_grad = 1e-03;
	number m_theta_cr = 34.0*3.1416/180.0;
	number m_vel = 0.15;

private:

	///    Data import for source
	SmartPtr<CplUserData<MathVector<dim>,dim> > m_imSource;
	Interface<dim>* Inter;

		  public:
	void set_theta(number data)
	{
		m_theta_cr = data*3.1416/180.0;
	}
	void set_vel(number data)
	{
		m_vel = data;
	}
	void set_gradient_limit(number data)
	{
		m_limit = data;
		m_limit_grad = 0.001 * data;
	}
	
	void set_phase_parameters(Interface<dim>* user)
	{
		if (!user) UG_THROW("Interface pointer is null!");
		if (!user->valid())
			UG_THROW("Interface parameters has not been initialized");
		Inter = user;
	}

public:
	/// constructor
	SlipVelocity(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
		
		if (spGridFct->num_fct() != dim+2)
			UG_THROW("NavierStokesMultiphase: Need exactly "<<dim+2<<" functions");
		for (int d=0;d<dim+2;d++)
		{
			if (spGridFct->local_finite_element_id(d) != LFEID(LFEID::LAGRANGE, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Lagrange P1 type.");
			}
		}
		m_u = spGridFct;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		grid.template attach_to<Vertex>(m_aNormal);
		//grid.template attach_to<Vertex>(m_aTang);
		grid.template attach_to<Vertex>(m_aVol);
		m_normal.access(grid,m_aNormal);
		//m_tang.access(grid,m_aTang);
		m_vol.access(grid,m_aVol);
		// set all values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_normal, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//SetAttachmentValues(m_tang, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//this->update();
	}

	virtual ~SlipVelocity(){};

	template <int refDim>
	inline void evaluate(MathVector<dim> vValue[],
						 const MathVector<dim> vGlobIP[],
						 number time, int si,
						 GridObject* elem,
						 const MathVector<dim> vCornerCoords[],
						 const MathVector<refDim> vLocIP[],
						 const size_t nip,
						 LocalVector* u,
						 const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);

		//    reference object id
		ReferenceObjectID roid = elem->reference_object_id();

		const size_t numVertices = element->num_vertices();
		const size_t MaxVertices = domain_traits<dim>::MaxNumVerticesOfElem;
		//    get domain of grid function
		const domain_type& domain = *m_u->domain().get();

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

//        position_accessor_type aaPos = m_u->domain()->position_accessor();

		// coord and vertex array
		MathVector<dim> coCoord[numVertices];
		Vertex* vVrt[numVertices];
		DimFV1Geometry<dim> geo;

		for(size_t i = 0; i < numVertices; ++i){
			vVrt[i] = element->vertex(i);
			coCoord[i] = posAcc[vVrt[i]];
		};

		// evaluate finite volume geometry
		geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

		// Lagrange 1 trial space
		const LocalShapeFunctionSet<refDim>& rTrialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));

		std::vector<number> shapes;
		

	//    storage for shape function at ip
		MathVector<refDim> vLocGrad[numVertices];
		MathVector<refDim> locGrad;

	//    Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		
		DimReferenceMapping<refDim, dim>& mapping = ReferenceMappingProvider::get<refDim, dim>(roid, coCoord);
		
		bool cut_elem=false;
		bool boolInside = false;
		Inter->cut_element(cut_elem,boolInside,  u,_C_);
		bool ComputeSlipVel =(cut_elem || boolInside || !boolInside)? true : false;
		
		for (size_t ip=0;ip<nip;ip++)
		{
			MathVector<dim> tang = 0.0;
			MathVector<dim> normal; VecSet(normal, 0.0);
			const number eps_n = 1e-10;
			const number eps_t = 1e-10;
			const number eps_slope = 1e-05;
			
			if(ComputeSlipVel)
			{
				MathVector<dim> GradC = 0.0;
				rTrialSpace.shapes(shapes,vLocIP[ip]);
				
				
				//    evaluate at shapes at ip
				rTrialSpace.grads(vLocGrad, vLocIP[ip]);
				//    compute grad at ip
				VecSet(locGrad, 0.0);
				for(size_t sh = 0; sh < numVertices; ++sh)
					VecScaleAppend(locGrad, (*u)(_C_, sh), vLocGrad[sh]);
				
				//    compute global grad
				mapping.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
				MatVecMult(GradC, JTInv, locGrad);
				
				number gradc_mag =  VecTwoNorm(GradC);
				
				for (size_t sh=0;sh<numVertices;sh++)
				{
					for(int d = 0; d < refDim; ++d)
					{
						normal[d] += m_normal[vVrt[sh]][d]*shapes[sh];
					}
				}
				number normal_mag =  sqrt(pow(VecTwoNorm(normal),2.0) + eps_n*eps_n);

				VecScale(normal,normal,1.0/normal_mag);
				
				number nxy2 = 0.0;
				for(int d = 0; d < dim-1; ++d)
					nxy2 += normal[d]*normal[d];

				number nxy = sqrt(nxy2 + eps_t*eps_t);
				number nz = sqrt(normal[dim-1]*normal[dim-1] + eps_n*eps_n);
				


				
				MathVector<dim> h = normal;
				h[dim-1] = 0.0;

				number hmag = sqrt(VecProd(h,h) + eps_t*eps_t);
				VecScale(h, h, 1.0 / hmag);
				
				MathVector<dim> z = 0.0;
				z[dim-1] = 1.0;
				VecScaleAdd(tang,cos(m_theta_cr),  h, - sin(m_theta_cr),z);
				
				
				
				
				number slope = nxy / nz;
				number Value = slope - tan(m_theta_cr);
				Value = (Value + sqrt(pow(Value,2.0) + eps_slope)) / 2.0;
				//Value = (Value + fabs(Value))/2.0;
				Value /= sqrt(1.0 + pow(slope,2.0));
				
				if(std::isnan(Value))
				{
					UG_LOG("  normal = " <<normal[0]<<"  "<<normal[1]<<".\n");
					
					UG_LOG("  normal0 = " <<m_normal[vVrt[0]][0]<<"  "<<m_normal[vVrt[0]][1]<<".\n");
					UG_LOG("  normal1 = " <<m_normal[vVrt[1]][0]<<"  "<<m_normal[vVrt[1]][1]<<".\n");
					UG_LOG("  normal2 = " <<m_normal[vVrt[2]][0]<<"  "<<m_normal[vVrt[2]][1]<<".\n");
					
					
					UG_THROW("Non valid number SlipVel  Value = " <<Value<<"  slope = "<< slope <<"  tang_mag = .\n");
				}
				VecScale(tang,tang,m_vel*Value);
			}
			//vValue[ip] = normal;
			//MathVector<dim> W = 0.0; W[dim-1] = -6.85985 * sin(30*3.1416/180); W[dim-2] = 6.85985 * cos(30*3.1416/180);
			//tang[dim-1]=0.0;
			vValue[ip] = tang;
			
			
			
		}
		
		
			
	}; // evaluate

	void update(){
		//    get domain
		UG_LOG("Updating Slip Velocity... \n");
		domain_type& domain = *m_u->domain().get();
		//    create Multiindex
		std::vector<DoFIndex> multInd;
		DimFV1Geometry<dim> geo;
		//    coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		MathVector<dim> coGrad[domain_traits<dim>::MaxNumVerticesOfElem];
		Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

		// set volume, tang and normal values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_normal, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//SetAttachmentValues(m_tang, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		// compute pressure in vertices by averaging
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			ElemIterator iter = m_u->template begin<elem_type>(si);
			ElemIterator iterEnd = m_u->template end<elem_type>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				elem_type* elem = *iter;
				const size_t numVertices = elem->num_vertices();
				for(size_t i = 0; i < numVertices; ++i){
					vVrt[i] = elem->vertex(i);
					coCoord[i] = posAcc[vVrt[i]];
				};
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				for(size_t i = 0; i < numVertices; ++i){
					number scvVol = geo.scv(i).volume();
					
					MathVector<dim> GradC; VecSet(GradC,0.0);
					MathVector<dim> Normal; VecSet(Normal,0.0);
					MathVector<dim> nn; VecSet(Normal,0.0); nn[dim-1] = -1.0;
					
					//    sum up contributions of each shape
					for(size_t sh = 0; sh < numVertices; ++sh)
					{
						m_u->dof_indices(elem->vertex(sh), _C_, multInd);
						//    read value of index from vector
						number uVal = DoFRef(*m_u,multInd[0]);
						//uVal = fmax(uVal, 0.0);
						//  Loop dimensions for derivative
						for(int d1 = 0; d1 <dim; ++d1)
						{
							GradC[d1] += uVal*geo.scv(i).global_grad(sh)[d1];
						}
					}
					const number grad_c_mag = VecTwoNorm(GradC);
					const number eps = 5e-01;
					const number p = 2.0;
					const number alpha = pow(grad_c_mag,p)/(pow(grad_c_mag,p) + eps) ;
					VecScaleAdd(GradC,alpha, GradC, 1.0-alpha,nn);
					const number GRAD_c_mag = VecTwoNorm(GradC);
					
					/*if(grad_c_mag<=m_limit_grad)
					{
						VecSet(Normal,0.0);
					}
					else
					{
						VecScale(Normal,GradC,-1.0/grad_c_mag);
						m_vol[vVrt[i]]+=scvVol;
					}*/
					
					VecScale(Normal,GradC,-1.0/GRAD_c_mag);
					m_vol[vVrt[i]]+=scvVol;
					
					for(int d1 = 0; d1 <dim; ++d1)
						m_normal[vVrt[i]][d1] += Normal[d1] * scvVol;
					
					
				}
			}
		}
		
		#ifdef UG_PARALLEL
			AttachmentAllReduce<Vertex> (*domain.grid(), m_aVol,    PCL_RO_SUM);
			AttachmentAllReduce<Vertex> (*domain.grid(), m_aNormal, PCL_RO_SUM);
		#endif
		
		PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
		// go over all vertices and average
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			VertexIterator iter = m_u->template begin<Vertex>(si);
			VertexIterator iterEnd = m_u->template end<Vertex>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				Vertex* vrt = *iter;
				if (pbm && pbm->is_slave(vrt)) continue;

				for(int d1 = 0; d1 <dim; ++d1)
					m_normal[vrt][d1] /= m_vol[vrt];

				/*MathVector<dim> tang; VecSet(tang,0.0);
				if(dim == 2)
				{
					number ss;
					if(fabs(m_normal[vrt][0])<m_limit)
						ss = 1.0;
					else
						ss = (m_normal[vrt][0]*m_normal[vrt][1]>0.0)? 1.0 : -1.0;
					
					tang[0] =   ss * m_normal[vrt][1];
					tang[1] = - fabs(m_normal[vrt][0]);
				}
				else if (dim == 3)
				{
					tang[0] = -m_normal[vrt][0]*m_normal[vrt][1];
					tang[1] = -m_normal[vrt][1]*m_normal[vrt][2];
					tang[2] =  m_normal[vrt][0]*m_normal[vrt][0] + m_normal[vrt][1]*m_normal[vrt][1];
					UG_THROW("SlipVelocity: Works only for Dim = 2");
				}
				else UG_THROW("SlipVelocity: Works only for Dim = 2 , 3");

				//VecScale(tang, tang,1.0/VecTwoNorm(tang));
				m_tang[vrt] = tang;*/
				
			}
		}
	}

private:
	static const size_t max_number_of_ips = 20;

public:
	virtual void operator() (MathVector<dim>& value,
							 const MathVector<dim>& globIP,
							 number time, int si) const
	{
		UG_THROW("SlipVel: Need element.");
	}

	virtual void operator() (MathVector<dim> vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si, const size_t nip) const
	{
		UG_THROW("SlipVel: Need element.");
	}

	virtual void compute(LocalVector* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), &(u->solution(this->time_point(s))));
	}

	///    returns if provided data is continuous over geometric object boundaries
	virtual bool continuous() const {return false;}

	///    returns if grid function is needed for evaluation
	virtual bool requires_grid_fct() const {return true;}
};

/**
concept derived from grid_function_user_data.h
 */
template <typename TGridFunction>
class SlipDiffusion
:     public StdUserData<SlipDiffusion<TGridFunction>, MathMatrix<TGridFunction::dim,TGridFunction::dim>, TGridFunction::dim>,
	  virtual public INewtonUpdate
{
	///    domain type
	typedef typename TGridFunction::domain_type domain_type;

	///    algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	///    world dimension
	static const int dim = domain_type::dim;
	///    Pressure
	static const int _P_ = domain_type::dim;
	///    Pressure
	static const int _C_ = domain_type::dim+1;

	///    grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// MathVector<dim> attachment
	typedef MathVector<dim> vecDim;
	typedef Attachment<vecDim> AMathVectorDim;

	/// attachment accessor
	typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;
	typedef PeriodicAttachmentAccessor<Vertex,AMathVectorDim > aVertexDimVector;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

private:

	//    Normal attachment accessor (average normal in vertices)
	AMathVectorDim m_aNormal;
	aVertexDimVector m_normal;
	
	//    Normal attachment accessor (average normal in vertices)
	//AMathVectorDim m_aTang;
	//aVertexDimVector m_tang;

	//  volume attachment accessor
	ANumber m_aVol;
	aVertexNumber m_vol;

	// level set grid function
	SmartPtr<TGridFunction> m_u;

	//    approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;
	
	//  grid
	grid_type* m_grid;
	
	number m_limit = 1e-03;
	number m_limit_grad = 1e-03;
	number m_theta_cr = 34.0*3.1416/180.0;
	number m_diff = 0.15;

private:

	///    Data import for Diffusion source
	SmartPtr<CplUserData<MathMatrix<dim,dim>,dim> > m_imDiffusion;
	///    Data import for normal
	SmartPtr<CplUserData<MathVector<dim>,dim> > m_imNormal;
	Interface<dim>* Inter;
	
public:
	/////////// Source

	void set_diffusion(SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > data)
	{
		m_imDiffusion = data;
	}
	void set_diffusion(number f_x)
	{
		SmartPtr<ConstUserMatrix<dim,dim> > f(new ConstUserMatrix<dim,dim>());
		f->set_diag_tensor(f_x);
		set_diffusion(f);
	}
	void set_normal(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
	{
		m_imNormal = data;
	}
	void set_normal(number f_x)
	{
		SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
		for (int i=0;i<dim;i++){
			f->set_entry(i, f_x);
		}
		set_normal(f);
	}
	void set_normal(number f_x, number f_y)
	{
		if (dim!=2){
			UG_THROW("NavierStokes: Setting source vector of dimension 2"
					" to a Discretization for world dim " << dim);
		} else {
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			set_normal(f);
		}
	}

	void set_normal(number f_x, number f_y, number f_z)
	{
		if (dim<3){
			UG_THROW("NavierStokes: Setting source vector of dimension 3"
					" to a Discretization for world dim " << dim);
		}
		else
		{
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			f->set_entry(2, f_z);
			set_normal(f);
		}
	}
#ifdef UG_FOR_LUA
	void set_normal(const char* fctName)
	{
		set_normal(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
	}
#endif

	void set_theta(number data)
	{
		m_theta_cr = data*3.1416/180.0;
	}
	void set_diff(number data)
	{
		m_diff = data;
	}
	void set_gradient_limit(number data)
	{
		m_limit = data;
		m_limit_grad = 0.001 * data;
	}
	
	void set_phase_parameters(Interface<dim>* user)
	{
		if (!user) UG_THROW("Interface pointer is null!");
		if (!user->valid())
			UG_THROW("Interface parameters has not been initialized");
		Inter = user;
	}

public:
	/// constructor
	SlipDiffusion(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
		
		if (spGridFct->num_fct() != dim+2)
			UG_THROW("NavierStokesMultiphase: Need exactly "<<dim+2<<" functions");
		for (int d=0;d<dim+2;d++)
		{
			if (spGridFct->local_finite_element_id(d) != LFEID(LFEID::LAGRANGE, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Lagrange P1 type.");
			}
		}
		m_u = spGridFct;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		set_diffusion(0.0);
		set_normal(0.0);

		grid.template attach_to<Vertex>(m_aNormal);
		//grid.template attach_to<Vertex>(m_aTang);
		grid.template attach_to<Vertex>(m_aVol);
		m_normal.access(grid,m_aNormal);
		//m_tang.access(grid,m_aTang);
		m_vol.access(grid,m_aVol);
		// set all values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_normal, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//SetAttachmentValues(m_tang, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		this->update();
	}

	virtual ~SlipDiffusion(){};

	template <int refDim>
	inline void evaluate(MathMatrix<dim,dim> vValue[],
						 const MathVector<dim> vGlobIP[],
						 number time, int si,
						 GridObject* elem,
						 const MathVector<dim> vCornerCoords[],
						 const MathVector<refDim> vLocIP[],
						 const size_t nip,
						 LocalVector* u,
						 const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);

		//    reference object id
		ReferenceObjectID roid = elem->reference_object_id();

		const size_t numVertices = element->num_vertices();
		const size_t MaxVertices = domain_traits<dim>::MaxNumVerticesOfElem;
		//    get domain of grid function
		const domain_type& domain = *m_u->domain().get();

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

//        position_accessor_type aaPos = m_u->domain()->position_accessor();

		// coord and vertex array
		MathVector<dim> coCoord[numVertices];
		Vertex* vVrt[numVertices];
		DimFV1Geometry<dim> geo;

		for(size_t i = 0; i < numVertices; ++i){
			vVrt[i] = element->vertex(i);
			coCoord[i] = posAcc[vVrt[i]];
		};

		// evaluate finite volume geometry
		geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

		// Lagrange 1 trial space
		const LocalShapeFunctionSet<refDim>& rTrialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));
		
		

		std::vector<number> shapes;
		

	//    storage for shape function at ip
		//MathVector<refDim> vLocGrad[numVertices];
		MathVector<refDim> locGrad;

	//    Reference Mapping
		//MathMatrix<dim, refDim> JTInv;
		
		DimReferenceMapping<refDim, dim>& mapping = ReferenceMappingProvider::get<refDim, dim>(roid, coCoord);
		
		(*m_imDiffusion)(vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
		
		std::vector<MathVector<dim> > vNormal(nip);
		(*m_imNormal)(&vNormal[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
		
		bool cut_elem = false;
		bool boolInside = false;
		Inter->cut_element(cut_elem,boolInside, u, _C_);
		bool ComputeSlipDiff =(cut_elem || boolInside)? true : false;
		
		for (size_t ip=0;ip<nip;ip++)
		{
			MathVector<dim> normal = 0.0;
			MathVector<dim> tang = 0.0;
			number Value = 0.0;
			
			if(ComputeSlipDiff)
			{
				//MathVector<dim> GradC = 0.0;
				rTrialSpace.shapes(shapes,vLocIP[ip]);
				
				//    evaluate at shapes at ip
				//rTrialSpace.grads(vLocGrad, vLocIP[ip]);
				//    compute grad at ip
				//VecSet(locGrad, 0.0);
				//for(size_t sh = 0; sh < numVertices; ++sh)
					//VecScaleAppend(locGrad, (*u)(_C_, sh), vLocGrad[sh]);
				
				//    compute global grad
				//mapping.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
				//MatVecMult(GradC, JTInv, locGrad);
				
				
				for (size_t sh=0;sh<numVertices;sh++)
					for(int d = 0; d < dim; ++d)
					{
						normal[d] += m_normal[vVrt[sh]][d]*shapes[sh];
						//tang[d] += m_tang[vVrt[sh]][d]*shapes[sh];
					}
				number normal_mag =  VecTwoNorm(normal);
				//number gradc_mag =  VecTwoNorm(GradC);
				if(normal_mag<m_limit)
				{
					VecSet(normal,0.0);
					normal[dim-1] = 1.0;
				}

				VecScale(normal,normal,1.0/normal_mag);
				
				number tang_mag =  0.0;
				for(int d = 0; d < dim-1; ++d)
					tang_mag += pow(normal[d],2.0);
				tang_mag = sqrt(tang_mag);
				
				VecSet(tang,0.0);
				if(tang_mag>m_limit)// && gradc_mag > m_limit)
				{
					for(int d = 0; d < dim-1; ++d)
						tang[d] = normal[d]*cos(m_theta_cr) / tang_mag;
					tang[dim-1] = -sin(m_theta_cr);
				}
				//if(gradc_mag > 1e-01)
				//if((cut_elem || boolInside) && VecProd(tang,-GradC) > 0)
				if(true)// VecProd(tang,-GradC) > 0)
				{
					number eps_slope = 1e-03;
					number theta = acos(normal[dim-1]);
					//UG_LOG("theta = "<< theta<<"\n");
					number slope = tan(theta);
					Value = (slope - tan(m_theta_cr))/sqrt(1.0 + pow(slope,2.0));
					//*sin(theta-m_theta_cr)
					Value = (Value + sqrt(pow(Value,2.0) + eps_slope)) / 2.0;
					//Value = (Value + fabs(Value))/2.0;
					Value *= m_diff;
				}
			}
			
			//for(int d1 = 0; d1 < dim; ++d1)
				//for(int d2 = 0; d2 < dim; ++d2)
					//vValue[ip](d1,d2) += Value * tang[d1]*tang[d2];
			for(int d1 = 0; d1 < dim; ++d1)
				vValue[ip](d1,d1) += Value;
		
		}
		
		
			
	}; // evaluate

	void update(){
		//    get domain
		UG_LOG("Updating Slip Diffusion... \n");
		domain_type& domain = *m_u->domain().get();
		//    create Multiindex
		std::vector<DoFIndex> multInd;
		DimFV1Geometry<dim> geo;
		//    coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		MathVector<dim> coGrad[domain_traits<dim>::MaxNumVerticesOfElem];
		Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

		// set volume, tang and normal values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_normal, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//SetAttachmentValues(m_tang, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		// compute pressure in vertices by averaging
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			ElemIterator iter = m_u->template begin<elem_type>(si);
			ElemIterator iterEnd = m_u->template end<elem_type>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				elem_type* elem = *iter;
				const size_t numVertices = elem->num_vertices();
				for(size_t i = 0; i < numVertices; ++i){
					vVrt[i] = elem->vertex(i);
					coCoord[i] = posAcc[vVrt[i]];
				};
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				for(size_t i = 0; i < numVertices; ++i){
					number scvVol = geo.scv(i).volume();
					
					MathVector<dim> GradC; VecSet(GradC,0.0);
					MathVector<dim> Normal; VecSet(Normal,0.0);
					
					//    sum up contributions of each shape
					for(size_t sh = 0; sh < numVertices; ++sh)
					{
						m_u->dof_indices(elem->vertex(sh), _C_, multInd);
						//    read value of index from vector
						number uVal = DoFRef(*m_u,multInd[0]);
						//uVal = fmax(uVal, 0.0);
						//  Loop dimensions for derivative
						for(int d1 = 0; d1 <dim; ++d1)
						{
							GradC[d1] += uVal*geo.scv(i).global_grad(sh)[d1];
						}
					}
					number grad_c_mag = VecTwoNorm(GradC);
					if(grad_c_mag<m_limit_grad)
					{
						VecSet(Normal,0.0);
					}
					else
					{
						VecScale(Normal,GradC,-1.0/grad_c_mag);
						m_vol[vVrt[i]]+=scvVol;
					}
					
					for(int d1 = 0; d1 <dim; ++d1)
						m_normal[vVrt[i]][d1] += Normal[d1] * scvVol;
					
					
				}
			}
		}
		
		#ifdef UG_PARALLEL
			AttachmentAllReduce<Vertex> (*domain.grid(), m_aVol, PCL_RO_SUM);
			AttachmentAllReduce<Vertex> (*domain.grid(), m_aNormal, PCL_RO_SUM);
		#endif
		
		PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
		// go over all vertices and average
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			VertexIterator iter = m_u->template begin<Vertex>(si);
			VertexIterator iterEnd = m_u->template end<Vertex>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				Vertex* vrt = *iter;
				if (pbm && pbm->is_slave(vrt)) continue;
				if(m_vol[vrt] > 1e-10)
				{
					for(int d1 = 0; d1 <dim; ++d1)
						m_normal[vrt][d1] /= m_vol[vrt];
				}
				else
				{
					VecSet(m_normal[vrt],0.0);
					m_normal[vrt][dim-1] = 1.0;
				}
				/*MathVector<dim> tang; VecSet(tang,0.0);
				if(dim == 2)
				{
					number ss;
					if(fabs(m_normal[vrt][0])<m_limit)
						ss = 1.0;
					else
						ss = (m_normal[vrt][0]*m_normal[vrt][1]>0.0)? 1.0 : -1.0;
					
					tang[0] =   ss * m_normal[vrt][1];
					tang[1] = - fabs(m_normal[vrt][0]);
				}
				else if (dim == 3)
				{
					tang[0] = -m_normal[vrt][0]*m_normal[vrt][1];
					tang[1] = -m_normal[vrt][1]*m_normal[vrt][2];
					tang[2] =  m_normal[vrt][0]*m_normal[vrt][0] + m_normal[vrt][1]*m_normal[vrt][1];
					UG_THROW("SlipDiffusion: Works only for Dim = 2");
				}
				else UG_THROW("SlipDiffusion: Works only for Dim = 2 , 3");

				//VecScale(tang, tang,1.0/VecTwoNorm(tang));
				m_tang[vrt] = tang;*/
				
			}
		}
	}

private:
	static const size_t max_number_of_ips = 20;

public:
	virtual void operator() (MathMatrix<dim,dim>& value,
							 const MathVector<dim>& globIP,
							 number time, int si) const
	{
		UG_THROW("SlipDiff: Need element.");
	}

	virtual void operator() (MathMatrix<dim,dim> vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si, const size_t nip) const
	{
		UG_THROW("SlipDiff: Need element.");
	}

	virtual void compute(LocalVector* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), &(u->solution(this->time_point(s))));
	}

	///    returns if provided data is continuous over geometric object boundaries
	virtual bool continuous() const {return false;}

	///    returns if grid function is needed for evaluation
	virtual bool requires_grid_fct() const {return true;}
};


/**
concept derived from grid_function_user_data.h
 */
template <typename TGridFunction>
class RelativeVelocity
:     public StdUserData<RelativeVelocity<TGridFunction>, MathVector<TGridFunction::dim>, TGridFunction::dim>,
	  virtual public INewtonUpdate
{
	///    domain type
	typedef typename TGridFunction::domain_type domain_type;

	///    algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	///    world dimension
	static const int dim = domain_type::dim;
	///    Pressure
	static const int _P_ = domain_type::dim;
	///    Pressure
	static const int _C_ = domain_type::dim+1;

	///    grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// MathVector<dim> attachment
	//typedef MathVector<dim> vecDim;
	//typedef Attachment<vecDim> AMathVectorDim;

	/// attachment accessor
	typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;
	//typedef PeriodicAttachmentAccessor<Vertex,AMathVectorDim > aVertexDimVector;
	//typedef Grid::AttachmentAccessor<elem_type,ANumber > aElementNumber;
	
	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

private:

	//    Normal attachment accessor (average normal in vertices)
	
	//  volume attachment accessor
	ANumber m_aRelVel;
	aVertexNumber m_rel_vel;

	// level set grid function
	SmartPtr<TGridFunction> m_u;

	//    approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;
	
	// subset group
	SubsetGroup m_relVelZeroSg;
	

	private:

	///    Data import for source
	//SmartPtr<CplUserData<MathVector<dim>,dim> > m_imPsGrad;
	
	///    Data import for Viscosity
	SmartPtr<CplUserData<number,dim> > m_imMixViscosity;
	//bool m_bvisc = false;
	
	Interface<dim> iface;
	Interface<dim>* Inter = &iface;

	public:
	/////////// Particle Pressure Grad import

	/*void set_ps_grad(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
	{
		m_imPsGrad = data;
	}

	void set_ps_grad(number f_x)
	{
		SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
		for (int i=0;i<dim;i++){
			f->set_entry(i, f_x);
		}
		set_ps_grad(f);
	}*/

	/*void set_viscosity(SmartPtr<CplUserData<number, dim> > data)
	{
		m_imMixViscosity = data;
		m_bvisc = true;
	}*/
	// set non-periodic boundaries so that viscosity can be set to zero there
	void setRelVelZeroBoundaries(const char* subsets){
		try{
			m_relVelZeroSg = m_u->subset_grp_by_name(subsets);
		}UG_CATCH_THROW("ERROR while parsing Subsets.");
	}
	
	void set_phase_parameters(Interface<dim>* user)
	{
		if (!user) UG_THROW("Interface pointer is null!");
		if (!user->valid())
			UG_THROW("Interface parameters has not been initialized");
		Inter = user;
	}

public:
	/// constructor
	RelativeVelocity(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
		
		if (spGridFct->num_fct() != dim+2)
			UG_THROW("NavierStokesMultiphase: Need exactly "<<dim+2<<" functions");
		for (int d=0;d<dim+2;d++)
		{
			if (spGridFct->local_finite_element_id(d) != LFEID(LFEID::LAGRANGE, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Lagrange P1 type.");
			}
		}
		if (!Inter) UG_THROW("Interface pointer is null!");
		
		m_u = spGridFct;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		
		//set_ps_grad(0.0);
		
		grid.template attach_to<Vertex>(m_aRelVel);
		
		m_rel_vel.access(grid,m_aRelVel);

		SetAttachmentValues(m_rel_vel, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 6.9);
		
		//this->update();
	}

	virtual ~RelativeVelocity(){};

	template <int refDim>
	inline void evaluate(MathVector<dim> vValue[],
						 const MathVector<dim> vGlobIP[],
						 number time, int si,
						 GridObject* elem,
						 const MathVector<dim> vCornerCoords[],
						 const MathVector<refDim> vLocIP[],
						 const size_t nip,
						 LocalVector* u,
						 const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		/*if(!m_bvisc)
			UG_THROW("RelativeVelocity StdUserData: Dynamic viscosity required");*/
			
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);

		//    reference object id
		ReferenceObjectID roid = elem->reference_object_id();

		const size_t numVertices = element->num_vertices();
		const size_t MaxVertices = domain_traits<dim>::MaxNumVerticesOfElem;
		//    get domain of grid function
		const domain_type& domain = *m_u->domain().get();

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

		// coord and vertex array
		MathVector<dim> coCoord[MaxVertices];
		Vertex* vVrt[MaxVertices];
		DimFV1Geometry<dim> geo;

		for(size_t i = 0; i < numVertices; ++i){
			vVrt[i] = element->vertex(i);
			coCoord[i] = posAcc[vVrt[i]];
		};

		// evaluate finite volume geometry
		try{
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
		}
		UG_CATCH_THROW("RelativeVel Export Parameter::evaluate:"
					   " Cannot update Finite Volume Geometry.");
		
		// Lagrange 1 trial space
		const LocalShapeFunctionSet<refDim>& rTrialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));

		std::vector<number> shapes;

		
		

		
		
		const number gy = Inter->gravity();
		const number rho_s = Inter->Density_s();
		const number rho_a =Inter->Density_a();
		const number dp =Inter->diameter();
		

		/*number  vViscosity[numVertices];
		(*m_imMixViscosity)(vViscosity, geo.scv_global_ips(), time, si, elem, coCoord, geo.scv_local_ips(), numVertices, u , NULL);*/
		
		number RelVelNodes[numVertices];
		
		for (size_t sh=0;sh<numVertices;sh++)
		{
			number W = m_rel_vel[vVrt[sh]];
			//size_t iter = 0;
			//Inter->RelVel(W,  iter,    (*u)(_C_,sh),   rho_a, rho_a, dp, rho_s,  fabs(gy));
			RelVelNodes[sh] = W;
		}
	
		
		for (size_t ip=0;ip<nip;ip++)
		{
			rTrialSpace.shapes(shapes,vLocIP[ip]);
			
			number W = 0.0;
			for (size_t sh=0;sh<numVertices;sh++)
				W += RelVelNodes[sh]*shapes[sh];

				
			MathVector<dim> RelVel; VecSet(RelVel,0.0);
			
			
			RelVel[dim-1] = -W;
			vValue[ip] = RelVel;
			
		}
		
	}; // evaluate

	void update(){
		//    get domain
		
		UG_LOG("Updating Relative Velocity... \n");
		domain_type& domain = *m_u->domain().get();
		//    create Multiindex
		std::vector<DoFIndex> multInd;
		DimFV1Geometry<dim> geo;
		//    coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		MathVector<dim> coGrad[domain_traits<dim>::MaxNumVerticesOfElem];
		Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

		
		const number gy = Inter->gravity();
		const number rho_s = Inter->Density_s();
		const number rho_a =Inter->Density_a();
		const number dp =Inter->diameter();


		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			ElemIterator iter = m_u->template begin<elem_type>(si);
			ElemIterator iterEnd = m_u->template end<elem_type>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				elem_type* elem = *iter;
				static const int refDim = elem_type::dim;
				const size_t numVertices = elem->num_vertices();
				for(size_t i = 0; i < numVertices; ++i){
					vVrt[i] = elem->vertex(i);
					coCoord[i] = posAcc[vVrt[i]];
				};
				
				//	create storage
				LocalIndices ind;
				LocalVector localU;

				// 	get global indices
				m_u->indices(elem, ind);
				// 	adapt local algebra
				localU.resize(ind);

				// 	read local values of u
				GetLocalVector(localU, *m_u);
				
				
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				
				/*number vViscosity[numVertices];
				(*m_imMixViscosity)(vViscosity, geo.scv_global_ips(), 0, si, elem, coCoord, geo.scv_local_ips(), geo.num_scv(), &localU, NULL);*/
				
				
				
				for(size_t i = 0; i < numVertices; ++i)
				{

					m_u->dof_indices(elem->vertex(i), _C_, multInd);
					//    read value of index from vector
					number uVal = DoFRef(*m_u,multInd[0]);
						

					
					number Ws = m_rel_vel[vVrt[i]];
					
					size_t iter = 0;

					Inter->RelVel(Ws,  iter,  uVal,   rho_a, rho_a, dp, rho_s,  fabs(gy));
					
					m_rel_vel[vVrt[i]] = Ws;
					
					
				}
			}
		}
		
	}

private:
	static const size_t max_number_of_ips = 20;

public:
	virtual void operator() (MathVector<dim>& value,
							 const MathVector<dim>& globIP,
							 number time, int si) const
	{
		UG_THROW("RelativeVelocity: Need element.");
	}

	virtual void operator() (MathVector<dim> vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si, const size_t nip) const
	{
		UG_THROW("RelativeVelocity: Need element.");
	}

	virtual void compute(LocalVector* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), &(u->solution(this->time_point(s))));
	}

	///    returns if provided data is continuous over geometric object boundaries
	virtual bool continuous() const {return false;}

	///    returns if grid function is needed for evaluation
	virtual bool requires_grid_fct() const {return true;}
};

template <typename TGridFunction>
class RelaxedParticleViscosity:
	public StdUserData<RelaxedParticleViscosity<TGridFunction>, number, TGridFunction::dim>,
	virtual public INewtonUpdate
{
	///    domain type
	typedef typename TGridFunction::domain_type domain_type;

	///    algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	///    world dimension
	static const int dim = domain_type::dim;
	///    Pressure
	static const int _P_ = domain_type::dim;
	///    Vol of fraction
	static const int _C_ = domain_type::dim+1;

	///    grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// MathVector<dim> attachment
	//        typedef MathVector<dim> vecDim;
	//        typedef Attachment<vecDim> AMathVectorDim;

	/// attachment accessor
	typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

private:

	
	//    Kinetic viscosity attachment accessor (New)
	ANumber m_aShearRate;
	aVertexNumber m_shear_rate;
	
	//    Kinetic viscosity attachment accessor (New)
	ANumber m_aViscNew;
	aVertexNumber m_visc_new;
	
	//    Kinetic viscosity attachment accessor (Old)
	ANumber m_aViscOld;
	aVertexNumber m_visc_old;

	//  volume attachment accessor
	ANumber m_aVol;
	aVertexNumber m_vol;

	// level set grid function
	SmartPtr<TGridFunction> m_u;

	//    approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;

private:

	
	Interface<dim> iface;
	Interface<dim>* Inter = &iface;

	int m_counter, m_skip_update;
	float m_relaxation_factor;
	bool m_bRelaxation;
	

public:
	
	
	void set_phase_parameters(Interface<dim>* user)
	{
		if (!user) UG_THROW("Interface pointer is null!");
		if (!user->valid())
			UG_THROW("Interface parameters has not been initialized");
		Inter = user;
	}
	void set_relaxation_parameters(bool bRelaxation, float relaxation_factor, int skip_update)
	{
		m_skip_update = skip_update;
		m_relaxation_factor = relaxation_factor;
		m_bRelaxation = bRelaxation;
		if (m_relaxation_factor < 0.0) UG_THROW("Invalid Relaxation factor!");
		if (!m_bRelaxation) m_relaxation_factor = 1.0;
		
	}

public:
	/// constructor
	RelaxedParticleViscosity(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct)
	{
		
		if (spGridFct->num_fct() != dim+2)
			UG_THROW("NavierStokesMultiphase: Need exactly "<<dim+2<<" functions");
		for (int d=0;d<dim+2;d++)
		{
			if (spGridFct->local_finite_element_id(d) != LFEID(LFEID::LAGRANGE, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Lagrange P1 type.");
			}
		}
		if (!Inter) UG_THROW("Interface pointer is null!");

		m_u = spGridFct;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		
		m_counter = -1;
		
		grid.template attach_to<Vertex>(m_aShearRate);
		grid.template attach_to<Vertex>(m_aViscNew);
		grid.template attach_to<Vertex>(m_aViscOld);
		grid.template attach_to<Vertex>(m_aVol);
		
		m_shear_rate.access(grid,m_aShearRate);
		m_visc_new.access(grid,m_aViscNew);
		m_visc_old.access(grid,m_aViscOld);
		m_vol.access(grid,m_aVol);
		
		// set all values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_shear_rate, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_visc_new, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_visc_old, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);

		this->update();
	}

	virtual ~RelaxedParticleViscosity(){};

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
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);

		//    reference object id
		ReferenceObjectID roid = elem->reference_object_id();

		const size_t numVertices = element->num_vertices();
		//    get domain of grid function
		const domain_type& domain = *m_u->domain().get();

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();


		// coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
		DimFV1Geometry<dim> geo;

		for(size_t i = 0; i < numVertices; ++i){
			vVrt[i] = element->vertex(i);
			coCoord[i] = posAcc[vVrt[i]];
		};

		// evaluate finite volume geometry
		try{
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
		}
		UG_CATCH_THROW("RelativeVel Export Parameter::evaluate:"
					   " Cannot update Finite Volume Geometry.");


		// Lagrange 1 trial space
		const LocalShapeFunctionSet<refDim>& rTrialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));

		std::vector<number> shapes;

		
		bool bAverageVisc = true;
		bool bAverageHarmonic = false;
		number AveValue = 0.0;
		
		if(bAverageVisc)
		{
			number value_new = 0.0;
			number value_old = 0.0;
			if(bAverageHarmonic)
			{
				UG_THROW("Harmonic Averaged Viscosity  non implemented.");
			}
			else
			{	number vol = 0.0;
				for (size_t sh=0;sh<numVertices;sh++)
				{
					value_new += geo.scv(sh).volume() * m_visc_new[vVrt[sh]];
					value_old += geo.scv(sh).volume() * m_visc_old[vVrt[sh]];
					vol += geo.scv(sh).volume();
				}
				value_new =  value_new / vol;
				value_old = value_old / vol;

			}
			
			AveValue = m_relaxation_factor * value_new +  (1.0 - m_relaxation_factor) * value_old;
		}

		for (size_t ip=0;ip<nip;ip++)
		{
			if(bAverageVisc)
				vValue[ip] = AveValue;
			else
			{
				
				number value_new = 0.0;
				number value_old = 0.0;
				rTrialSpace.shapes(shapes,vLocIP[ip]);
				for (size_t sh=0;sh<numVertices;sh++)
				{
					value_new += m_visc_new[vVrt[sh]]*shapes[sh];
					value_old += m_visc_old[vVrt[sh]]*shapes[sh];
				}

				vValue[ip] = m_relaxation_factor * value_new + (1.0 - m_relaxation_factor) * value_old;
			}
			
		}
		
		
		
			
	}; // evaluate

	void update(){
		//    get domain
		if(m_counter>=0)
		{
			++m_counter;
			if( m_counter % (m_skip_update+1) != 0)
			{
				//UG_LOG("Skipping...\n");
				return ;
			}
			else m_counter = 0;
		}
		
		UG_LOG("Updating Relaxed Viscosity...\n");
		
		domain_type& domain = *m_u->domain().get();
		//    create Multiindex
		std::vector<DoFIndex> multInd;
		DimFV1Geometry<dim> geo;
		//    coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		MathVector<dim> coGrad[domain_traits<dim>::MaxNumVerticesOfElem];
		Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();
		
		// set volume and p values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_shear_rate, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//SetAttachmentValues(m_visc_new, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		// compute pressure in vertices by averaging
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			ElemIterator iter = m_u->template begin<elem_type>(si);
			ElemIterator iterEnd = m_u->template end<elem_type>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				elem_type* elem = *iter;
				const size_t numVertices = elem->num_vertices();
				for(size_t i = 0; i < numVertices; ++i){
					vVrt[i] = elem->vertex(i);
					coCoord[i] = posAcc[vVrt[i]];
				};
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				for(size_t i = 0; i < numVertices; ++i)
				{
					number scvVol = geo.scv(i).volume();
					m_vol[vVrt[i]]+=scvVol;
					
					MathMatrix<dim,dim> VelGrad; MatSet(VelGrad,0.0);
					
					//    sum up contributions of each shape
					for(size_t sh = 0; sh < numVertices; ++sh)
					{
						//  Loop dimensions for derivative
						for(int d1 = 0; d1 <dim; ++d1)
						{
							m_u->dof_indices(elem->vertex(sh), d1, multInd);
							//    read value of index from vector
							number uVal = DoFRef(*m_u,multInd[0]);
						//  Loop dimensions for direction
							for(int d2 = 0; d2 < dim; ++d2)
							{
								VelGrad(d1, d2) += uVal*geo.scv(i).global_grad(sh)[d2];
							}
						}
					}
					number gamma=0.0;
					// compute inner sum
					for(int d1 = 0; d1 < dim; ++d1)
					{
						for(int d2 = 0; d2 < dim; ++d2)
						{
							gamma += pow((VelGrad(d1,d2) + VelGrad(d2,d1)),2);
						}
					}
					
					gamma =sqrt((0.5*gamma));
					
					
					m_shear_rate[vVrt[i]] += gamma * scvVol;
					
				}
			}
		}
		
		#ifdef UG_PARALLEL
			AttachmentAllReduce<Vertex> (*domain.grid(), m_aVol, PCL_RO_SUM);
			AttachmentAllReduce<Vertex> (*domain.grid(), m_aShearRate, PCL_RO_SUM);
		#endif
	
		PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
		// go over all vertices and average
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			VertexIterator iter = m_u->template begin<Vertex>(si);
			VertexIterator iterEnd = m_u->template end<Vertex>(si);
			
			for(  ;iter !=iterEnd; ++iter)
			{
				Vertex* vrt = *iter;
				
				if (pbm && pbm->is_slave(vrt)) continue;
				
				m_shear_rate[vrt] /= m_vol[vrt];
	
				number C_value, mu_s, Dmu_s, gamma_average;
				number mu_eins = 0.0;
				number Dmu_eins = 0.0;

				m_u->dof_indices(vrt, _C_, multInd);
				C_value=DoFRef(*m_u,multInd[0]);
				gamma_average = m_shear_rate[vrt];


				Inter->MU_I_Viscosity( mu_s, Dmu_s,gamma_average, C_value, false);
				//Inter->Einstein_viscosity( mu_eins, Dmu_eins, C_value, false);
				
				if (m_counter<0)
				{
					m_visc_old[vrt] = mu_s;
					m_counter = 0;
				}
				else
					m_visc_old[vrt] = m_visc_new[vrt];
				
				m_visc_new[vrt] = mu_s;

				
				
			}
			
		}
		
	}

		  private:
	static const size_t max_number_of_ips = 20;

		  public:
	virtual void operator() (number& value,
							 const MathVector<dim>& globIP,
							 number time, int si) const
	{
		UG_THROW("ShearStressUserData: Need element.");
	}

	virtual void operator() (number vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si, const size_t nip) const
	{
		UG_THROW("ShearStress: Need element.");
	}

	virtual void compute(LocalVector* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, vCornerCoords, this->template local_ips<dim>(s),
						  this->num_ip(s), &(u->solution(this->time_point(s))));
	}

	///    returns if provided data is continuous over geometric object boundaries
	virtual bool continuous() const {return false;}

	///    returns if grid function is needed for evaluation
	virtual bool requires_grid_fct() const {return true;}
};




} // namespace NavierStokes
} // end namespace ug


#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__SHEAR_STRESS__ */
