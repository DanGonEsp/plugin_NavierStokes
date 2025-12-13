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
    typedef Grid::AttachmentAccessor<elem_type,ANumber > aElementNumber;

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
        const LocalShapeFunctionSet<dim>& lagrange1 =
                LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

        std::vector<number> shapes;
        for (size_t ip=0;ip<nip;ip++)
        {
            number value = 0.0;
            MathVector<refDim> LocalCoord_aux;
            for(int d = 0; d < refDim; ++d)
                LocalCoord_aux[d]=vLocIP[ip][d];
            lagrange1.shapes(shapes,LocalCoord_aux[ip]);
            for (size_t sh=0;sh<numVertices;sh++)
                value += m_shear_rate[vVrt[sh]]*shapes[sh];
            
            vValue[ip] = value;
            
        }
        
        
            
    }; // evaluate

    void update(){
        //    get domain
        printf("Updating Velocity Grad... \n");
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
	typedef Grid::AttachmentAccessor<elem_type,ANumber > aElementNumber;

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
		const LocalShapeFunctionSet<dim>& lagrange1 =
				LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

		std::vector<number> shapes;
		number Ps[numVertices];
		number Gamma[numVertices];
		for (size_t sh=0;sh<numVertices;sh++)
			Gamma[sh] = m_shear_rate[vVrt[sh]];
		
		//Inter->Ps( Ps, NULL, Gamma, *u, _C_, numVertices, false);
		
		for (size_t ip=0;ip<nip;ip++)
		{
			number value = 0.0;
			MathVector<refDim> LocalCoord_aux;
			for(int d = 0; d < refDim; ++d)
				LocalCoord_aux[d]=vLocIP[ip][d];
			lagrange1.shapes(shapes,LocalCoord_aux[ip]);
			
			
			for (size_t sh=0;sh<numVertices;sh++)
				value += Gamma[sh]*shapes[sh];
			
			vValue[ip] = value;
			
		}
		
		
			
	}; // evaluate

	void update(){
		//    get domain
		printf("Updating Velocity Grad... \n");
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



} // namespace NavierStokes
} // end namespace ug


#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__SHEAR_STRESS__ */
