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

namespace ug{
namespace NavierStokes{

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
        this->update();
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
        //printf("Updating Velocity Grad... \n");
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
                          elem, NULL, this->template local_ips<dim>(s),
                          this->num_ip(s), u);
    }

    virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
                         const MathVector<dim> vCornerCoords[], bool bDeriv = false)
    {
        const int si = this->subset();
        for(size_t s = 0; s < this->num_series(); ++s)
            evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
                          elem, NULL, this->template local_ips<dim>(s),
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
		this->update();
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
		//printf("Updating Velocity Grad... \n");
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
						  elem, NULL, this->template local_ips<dim>(s),
						  this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, NULL, this->template local_ips<dim>(s),
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
	
	void set_phase_parameters(Interface<dim>* user)
	{
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
		this->update();
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
		for (size_t ip=0;ip<nip;ip++)
		{
			MathVector<dim> normal = 0.0;
			MathVector<dim> GradC = 0.0;
			MathVector<dim> tang = 0.0;

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
			
			
			for (size_t sh=0;sh<numVertices;sh++)
				for(int d = 0; d < refDim; ++d)
				{
					normal[d] += m_normal[vVrt[sh]][d]*shapes[sh];
					//tang[d] += m_tang[vVrt[sh]][d]*shapes[sh];
				}
			number normal_mag =  VecTwoNorm(normal);
			number gradc_mag =  VecTwoNorm(GradC);
			if(normal_mag<m_limit)
			{
				VecSet(normal,0.0);
				normal[dim-1] = 1.0;
			}
			else
				VecScale(normal,normal,1.0/normal_mag);
			
			number tang_mag =  0.0;
			for(int d = 0; d < refDim-1; ++d)
				tang_mag += pow(normal[d],2.0);
			tang_mag = sqrt(tang_mag);
	
			VecSet(tang,0.0);
			if(tang_mag>m_limit && gradc_mag > m_limit)
			{
				for(int d = 0; d < refDim-1; ++d)
					tang[d] = normal[d]*cos(m_theta_cr) / tang_mag;
				tang[dim-1] = -sin(m_theta_cr);
			}
			
			number theta = acos(normal[dim-1]);
			number slope = fabs(tan(theta));
			number Value = slope - tan(m_theta_cr);
			Value = (Value + fabs(Value))/2.0;
			Value /= sqrt(1.0 + pow(slope,2.0));
		

			VecScale(tang,tang,m_vel*Value);
			vValue[ip] = tang;
			
			
			
		}
		
		
			
	}; // evaluate

	void update(){
		//    get domain
		//printf("Updating Slip Velocity... \n");
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
						
						//  Loop dimensions for derivative
						for(int d1 = 0; d1 <dim; ++d1)
						{
							GradC[d1] += uVal*geo.scv(i).global_grad(sh)[d1];
						}
					}
					number grad_c_mag = VecTwoNorm(GradC);
					if(grad_c_mag<m_limit)
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
						  elem, NULL, this->template local_ips<dim>(s),
						  this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, NULL, this->template local_ips<dim>(s),
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
	
	number m_limit = 1e-02;
	number m_theta_cr = 34.0*3.1416/180.0;
	number m_vel = 0.15;

private:

	///    Data import for source
	SmartPtr<CplUserData<MathMatrix<dim,dim>,dim> > m_imDiffusion;
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

	void set_theta(number data)
	{
		m_theta_cr = data*3.1416/180.0;
	}
	void set_vel(number data)
	{
		m_vel = data;
	}
	
	void set_phase_parameters(Interface<dim>* user)
	{
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
		MathVector<refDim> vLocGrad[numVertices];
		MathVector<refDim> locGrad;

	//    Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		
		DimReferenceMapping<refDim, dim>& mapping = ReferenceMappingProvider::get<refDim, dim>(roid, coCoord);
		
		(*m_imDiffusion)(vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
		
		bool cut_elem = false;
		Inter->cut_element(cut_elem, u, _C_);
		
		for (size_t ip=0;ip<nip;ip++)
		{
			MathVector<dim> normal = 0.0;
			MathVector<dim> GradC = 0.0;
			MathVector<dim> tang = 0.0;

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
			
			
			for (size_t sh=0;sh<numVertices;sh++)
				for(int d = 0; d < dim; ++d)
				{
					normal[d] += m_normal[vVrt[sh]][d]*shapes[sh];
					//tang[d] += m_tang[vVrt[sh]][d]*shapes[sh];
				}
			number normal_mag =  VecTwoNorm(normal);
			number gradc_mag =  VecTwoNorm(GradC);
			if(normal_mag<m_limit)
			{
				VecSet(normal,0.0);
				normal[dim-1] = 1.0;
			}
			else
				VecScale(normal,normal,1.0/normal_mag);
			
			number tang_mag =  0.0;
			for(int d = 0; d < dim-1; ++d)
				tang_mag += pow(normal[d],2.0);
			tang_mag = sqrt(tang_mag);
	
			VecSet(tang,0.0);
			if(tang_mag>m_limit)
			{
				for(int d = 0; d < dim-1; ++d)
					tang[d] = normal[d]*cos(m_theta_cr) / tang_mag;
				tang[dim-1] = -sin(m_theta_cr);
			}
			number Value = 0.0;
			//if(gradc_mag > 1e-01)
			if(cut_elem && VecProd(tang,-GradC) > 0)
			{
				number theta = acos(normal[dim-1]);
				number slope = fabs(tan(theta));
				Value = slope - tan(m_theta_cr);
				Value = (Value + fabs(Value))/2.0;
				Value *= m_vel*sin(theta-m_theta_cr)/sqrt(1.0 + pow(slope,2.0));
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
		//printf("Updating Slip Diffusion... \n");
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
					m_vol[vVrt[i]]+=scvVol;
					
					MathVector<dim> GradC; VecSet(GradC,0.0);
					MathVector<dim> Normal; VecSet(Normal,0.0);
					
					//    sum up contributions of each shape
					for(size_t sh = 0; sh < numVertices; ++sh)
					{
						m_u->dof_indices(elem->vertex(sh), _C_, multInd);
						//    read value of index from vector
						number uVal = DoFRef(*m_u,multInd[0]);
						
						//  Loop dimensions for derivative
						for(int d1 = 0; d1 <dim; ++d1)
						{
							GradC[d1] += uVal*geo.scv(i).global_grad(sh)[d1];
						}
					}
					number grad_c_mag = VecTwoNorm(GradC);
					if(grad_c_mag<m_limit)
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
						  elem, NULL, this->template local_ips<dim>(s),
						  this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, NULL, this->template local_ips<dim>(s),
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
	//typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;
	//typedef PeriodicAttachmentAccessor<Vertex,AMathVectorDim > aVertexDimVector;
	typedef Grid::AttachmentAccessor<elem_type,ANumber > aElementNumber;
	
	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
	//typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

private:

	//    Normal attachment accessor (average normal in vertices)
	//AMathVectorDim m_aNormal;
	//aVertexDimVector m_normal;
	
	//  volume attachment accessor
	//ANumber m_aVol;
	//aVertexNumber m_vol;
	
	//  volume attachment accessor
	ANumber m_aRelVel;
	aElementNumber m_rel_vel;

	// level set grid function
	SmartPtr<TGridFunction> m_u;

	//    approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;
	
	number m_limit = 1e-03;
	number m_theta_cr = 34.0*3.1416/180.0;
	number m_vel = 0.15;

private:

	///    Data import for source
	SmartPtr<CplUserData<MathVector<dim>,dim> > m_imRelVel;
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
	
	void set_phase_parameters(Interface<dim>* user)
	{
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
		m_u = spGridFct;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		//grid.template attach_to<Vertex>(m_aNormal);
		//grid.template attach_to<Vertex>(m_aVol);
		grid.template attach_to<elem_type>(m_aRelVel);
		
		//m_normal.access(grid,m_aNormal);
		//m_vol.access(grid,m_aVol);
		m_rel_vel.access(grid,m_aRelVel);
		// set all values to zero
		//SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//SetAttachmentValues(m_normal, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_rel_vel, m_u->template begin<elem_type>(), m_u->template end<elem_type>(), 6.9);
		this->update();
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
			MathVector<dim> GradC = 0.0;
			MathVector<dim> tang = 0.0;

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
			
			
			/*for (size_t sh=0;sh<numVertices;sh++)
				for(int d = 0; d < refDim; ++d)
				{
					normal[d] += m_normal[vVrt[sh]][d]*shapes[sh];
					//tang[d] += m_tang[vVrt[sh]][d]*shapes[sh];
				}*/
			number normal_mag =  VecTwoNorm(normal);
			number gradc_mag =  VecTwoNorm(GradC);
			if(normal_mag<m_limit)
			{
				VecSet(normal,0.0);
				normal[dim-1] = 1.0;
			}
			else
				VecScale(normal,normal,1.0/normal_mag);
			
			number tang_mag =  0.0;
			for(int d = 0; d < refDim-1; ++d)
				tang_mag += pow(normal[d],2.0);
			tang_mag = sqrt(tang_mag);
	
			VecSet(tang,0.0);
			if(tang_mag>m_limit && gradc_mag > m_limit)
			{
				for(int d = 0; d < refDim-1; ++d)
					tang[d] = normal[d]*cos(m_theta_cr) / tang_mag;
				tang[dim-1] = -sin(m_theta_cr);
			}
			
			number theta = acos(normal[dim-1]);
			number slope = fabs(tan(theta));
			number Value = slope - tan(m_theta_cr);
			Value = (Value + fabs(Value))/2.0;
			Value /= sqrt(1.0 + pow(slope,2.0));
		

			VecScale(tang,tang,m_vel*Value);
			vValue[ip] = tang;
			
		}
		
	}; // evaluate

	void update(){
		//    get domain
		//printf("Updating Relative Velocity... \n");
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
		//SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		//SetAttachmentValues(m_normal, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
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
				number W1 = m_rel_vel[elem];
				//m_rel_vel[elem] = Inter->RelVel(W1,  iter,    MU2,   rho_a, rho_a, dp, rho_s,  fabs(gy), 5.0);
				//m_pOld[elem]+=DoFRef(*m_u,multInd[0]);
				
			}
		}
		/*PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
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
			}
		}*/
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
						  elem, NULL, this->template local_ips<dim>(s),
						  this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
						 const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
						  elem, NULL, this->template local_ips<dim>(s),
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
