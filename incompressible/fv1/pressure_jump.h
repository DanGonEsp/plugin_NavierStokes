/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__PRESSURE__JUMP__MODEL__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__PRESSURE__JUMP__MODEL__

//#define UG_NSSTAB_ASSERT(cond, exp)
// include define below to assert arrays used in stabilization
#define UG_NSSTAB_ASSERT(cond, exp) UG_ASSERT((cond), (exp))

#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "../../upwind_interface.h"

namespace ug{
namespace NavierStokes{

template <int dim>
class INavierStokesPressureJump
{
    protected:
    ///    used traits
        typedef fv1_dim_traits<dim, dim> traits;

        enum DiffusionLength
        {
            RAW = 0,
            FIVEPOINT,
            COR
        };

    public:
    ///    number of SubControlVolumes
        static const size_t maxNumSCV = traits::maxNumSCV;

    ///    max number of SubControlVolumeFaces
        static const size_t maxNumSCVF = traits::maxNumSCVF;

    /// max number of shape functions
        static const size_t maxNumSH = traits::maxNSH;

    private:
    /// abbreviation for own type
        typedef INavierStokesPressureJump<dim> this_type;

    public:
    ///    constructor
        INavierStokesPressureJump()
        :    m_numScv(0), m_numSh(0), m_spUpwind(NULL)
        {
            m_vUpdateFunc.clear();
            //    default setup
            set_diffusion_length("RAW");
        }
    
    ///    sets the type of diff length used for evaluation
        void set_diffusion_length(std::string diffLength);
    
    ///    sets the upwind method
        void set_upwind(SmartPtr<INavierStokesUpwind<dim> > spUpwind)
        {
            m_spUpwind = spUpwind;
            m_spConstUpwind = spUpwind;
        }

    ///    returns the upwind
        const ConstSmartPtr<INavierStokesUpwind<dim> >& upwind() const {return m_spConstUpwind;}
        
    /// returns if the upwind pointer is valid
        bool upwind_valid() const {return m_spConstUpwind.valid();}

    ///    set the FV1 Geometry type to use for next updates
        template <typename TFVGeom>
        void set_geometry_type()
        {
        //    get unique geometry id
            size_t id = GetUniqueFVGeomID<TFVGeom>();
        //    check that function exists
            if(id >= m_vUpdateFunc.size() || m_vUpdateFunc[id] == NULL)
                UG_THROW("No update function registered for Geometry "<<id);
        //    set current geometry
            m_id = id;

        //    set sizes
            TFVGeom& geo = GeomProvider<TFVGeom>::get();
            m_numScv = geo.num_scv();
            m_numSh = geo.num_sh();
            if(m_spUpwind.valid()) m_spUpwind->template set_geometry_type<TFVGeom>();

        }
    /////////////////////////////////////////
    // interface to the upwind
    /////////////////////////////////////////
    
    protected:
    
    ///    computes the upwind shapes
        template <typename TFVGeom>
        void compute_upwind(const TFVGeom& geo,
                            const MathVector<dim> vStdVel[])
        {
            UG_NSSTAB_ASSERT(m_spUpwind.valid(), "No upwind object");
            m_spUpwind->update_upwind(geo, vStdVel);
        };

    ///    computes the downwind shapes
        template <typename TFVGeom>
        void compute_downwind(const TFVGeom& geo,
                              const MathVector<dim> vStdVel[])
        {
            UG_NSSTAB_ASSERT(m_spUpwind.valid(), "No upwind object");
            m_spUpwind->update_downwind(geo, vStdVel);
        };
    
        
    /////////////////////////////////////////
    // the data interface (for the NS discretization)
    /////////////////////////////////////////

    public:
    
    ///    number of integration points
        size_t num_ip () const {return m_numScv;}
        
    ///    number of shapes (corners)
        size_t num_sh () const {return m_numSh;}
    
    ///    returns if stab velocity comp depends on other vel components
        bool press_jump_comp_connected() const {return m_bPressJumpCompConnected;}
    
    /// pressure jump
        const number& pressure_jump(size_t sh) const
        {
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vPressureJump[sh];
        }

    /// computed stab shape for velocity. This is: The pressure_jump derivative
    /// w.r.t velocity unknowns in the corner for each component
        number pressure_shape_vel(size_t sh1, size_t d, size_t sh2) const
        {
            UG_NSSTAB_ASSERT(sh1 < m_numSh, "Invalid index.");
            UG_NSSTAB_ASSERT(d < dim,       "Invalid index.");
            UG_NSSTAB_ASSERT(sh2 < m_numSh, "Invalid index.");
            return m_vvvvPressureShapeVel[sh1][d][sh2];
        }

    ///    computed stab shape for pressure.
        number pressure_shape_p(size_t sh1, size_t sh2) const
        {
            UG_NSSTAB_ASSERT(sh1 < m_numSh, "Invalid index.");
            UG_NSSTAB_ASSERT(sh2 < m_numSh, "Invalid index.");
            return m_vvvPressureShapeP[sh1][sh2];
        }
    
    /// Tangent vel
        const MathVector<dim>& tang_vel(size_t sh) const
        {
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vTangVel[sh];
        }
        number tang_vel_shape_vel(size_t sh1, size_t compOut, size_t compIn, size_t sh) const
        {
            UG_NSSTAB_ASSERT(sh1 < m_numSh, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(compIn < dim,  "Invalid index.");
            UG_NSSTAB_ASSERT(sh < m_numSh,  "Invalid index.");
            return m_vvvvTangVelShapeVel[sh1][compOut][compIn][sh];
        }



    ///    compute values for new geometry and corner velocities
        void update(const FVGeometryBase* geo,
                    const LocalVector& vCornerValue,
                    const MathVector<dim> vStdVel[],
                    const bool bStokes,
                    const DataImport<MathVector<dim>, dim>& n,
                    const DataImport<number, dim>& kinViscoSCV,
                    const DataImport<number, dim>& density,
                    const DataImport<number, dim>& densitySCV,
                    const number JumpShape[],
                    const LocalVector* pvCornerValueOldTime, number dt,
                    const number density_ref,
                    const number mu_l,
                    const number rho_l,
                    const MathVector<dim> Source_l,
                    const number mu_g,
                    const number rho_g,
                    const MathVector<dim> Source_g,
                    const number interface_value)
    {(this->*(m_vUpdateFunc[m_id]))(geo, vCornerValue, vStdVel, bStokes,n, kinViscoSCV, density, densitySCV, JumpShape, pvCornerValueOldTime, dt, density_ref, mu_l, rho_l, Source_l, mu_g, rho_g,Source_g, interface_value);}


    /////////////////////////////////////////
    // the data interface (for the implementation)
    /////////////////////////////////////////

    protected:
    
    ///    sets the pressure comp connected flag
        void set_pressure_jump_comp_connected(bool bPressJumpCompConnected) {m_bPressJumpCompConnected = bPressJumpCompConnected;}
    
    /// Pressure jump
        number& pressure_jump(size_t sh)
        {
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vPressureJump[sh];
        }

    /// computed stab shape for velocity. This is: The pressure_jump derivative
    /// w.r.t velocity unknowns in the corner for each component
        number& pressure_shape_vel(size_t sh1, size_t d, size_t sh2)
        {
            UG_NSSTAB_ASSERT(sh1 < m_numSh, "Invalid index.");
            UG_NSSTAB_ASSERT(d < dim,       "Invalid index.");
            UG_NSSTAB_ASSERT(sh2 < m_numSh, "Invalid index.");
            return m_vvvvPressureShapeVel[sh1][d][sh2];
        }

    ///    computed stab shape for pressure.
        number& pressure_shape_p(size_t sh1, size_t sh2)
        {
            UG_NSSTAB_ASSERT(sh1 < m_numSh, "Invalid index.");
            UG_NSSTAB_ASSERT(sh2 < m_numSh, "Invalid index.");
            return m_vvvPressureShapeP[sh1][sh2];
        }
    
    /// tangent vector
        MathVector<dim>& tang_vel(size_t sh)
        {
            UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
            return m_vTangVel[sh];
        }
    
        number& tang_vel_shape_vel(size_t sh1, size_t compOut, size_t compIn, size_t sh2)
        {
            UG_NSSTAB_ASSERT(sh1 < m_numSh, "Invalid index.");
            UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
            UG_NSSTAB_ASSERT(compIn < dim,  "Invalid index.");
            UG_NSSTAB_ASSERT(sh2 < m_numSh, "Invalid index.");
            return m_vvvvTangVelShapeVel[sh1][compOut][compIn][sh2];
        }
    ///    diff length
        number diff_length_sq_inv(size_t scvf) const
        {
            UG_NSSTAB_ASSERT(scvf < this_type::num_ip(), "Invalid index.");
            return m_vDiffLengthSqInv[scvf];
        }
    
    /////////////////////////////////////////
    // forward methods of Upwind Velocity
    /////////////////////////////////////////

    ///    Convection Length
        number upwind_conv_length(size_t scvf) const
        {
            UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
            return this_type::upwind()->upwind_conv_length(scvf);
        }

    ///    Convection Length
        number downwind_conv_length(size_t scvf) const
        {
            UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
            return this_type::upwind()->downwind_conv_length(scvf);
        }

    ///    upwind shape for corner vel
        number upwind_shape_sh(size_t scvf, size_t sh) const
        {
            UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
            return this_type::upwind()->upwind_shape_sh(scvf, sh);
        }

    ///    upwind shape for corner vel
        number downwind_shape_sh(size_t scvf, size_t sh) const
        {
            UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
            return this_type::upwind()->downwind_shape_sh(scvf, sh);
        }

    ///    returns if upwind shape w.r.t. ip vel is non-zero
        bool non_zero_shape_ip() const
        {
            UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
            return this_type::upwind()->non_zero_shape_ip();
        }

    ///    upwind shapes for ip vel
        number upwind_shape_ip(size_t scvf, size_t scvf2) const
        {
            UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
            return this_type::upwind()->upwind_shape_ip(scvf, scvf2);
        }

    ///    upwind shapes for ip vel
        number downwind_shape_ip(size_t scvf, size_t scvf2) const
        {
            UG_NSSTAB_ASSERT(this_type::upwind_valid(), "No upwind object");
            return this_type::upwind()->downwind_shape_ip(scvf, scvf2);
        }
    
    ///    computes the diffusion length
        template <typename TFVGeom>
        void compute_diff_length(const TFVGeom& geo);

    /////////////////////////////////////////
    // the data
    /////////////////////////////////////////
    
    private:

    ///    number of current scv
        size_t m_numScv;

    ///    number of current shape functions (usually in corners)
        size_t m_numSh;

    ///    values of pressure Jump at ip
        number m_vPressureJump[maxNumSH];
    
    ///    pressure jump shapes w.r.t vel
        number m_vvvvPressureShapeVel[maxNumSH][dim][maxNumSH];

    ///    pressure shapes w.r.t pressure
        number m_vvvPressureShapeP[maxNumSH][maxNumSH];
    
        MathVector<dim> m_vTangVel[maxNumSH];

    ///    velocity jump shapes w.r.t vel
        number m_vvvvTangVelShapeVel[maxNumSH][dim][dim][maxNumSH];
    
    ///    type of diffusion length computation
        DiffusionLength m_diffLengthType;

    ///    vector holding diffusion Length squared and inverted
        number m_vDiffLengthSqInv[this_type::maxNumSCVF];


    
    ///    flag if pressure jump components are interconnected
        bool m_bPressJumpCompConnected;
    


    ///    id of current geometry type
        int m_id;
    
    ///    Upwind object (if set)
        SmartPtr<INavierStokesUpwind<dim> > m_spUpwind;
    ///    Upwind object (if set), the same as m_spUpwind, but as a const
        ConstSmartPtr<INavierStokesUpwind<dim> > m_spConstUpwind;
    
    //////////////////////////
    // registering mechanism
    //////////////////////////

    protected:
    ///    type of update function
        typedef void (this_type::*UpdateFunc)(const FVGeometryBase* geo,
                                              const LocalVector& vCornerValue,
                                              const MathVector<dim> vStdVel[],
                                              const bool bStokes,
                                              const DataImport<MathVector<dim>, dim>& n,
                                              const DataImport<number, dim>& kinViscoSCV,
                                              const DataImport<number, dim>& density,
                                              const DataImport<number, dim>& densitySCV,
                                              const number JumpShape[],
                                              const LocalVector* pvCornerValueOldTime, number dt,
                                              const number density_ref,
                                              const number mu_l,
                                              const number rho_l,
                                              const MathVector<dim> Source_l,
                                              const number mu_g,
                                              const number rho_g,
                                              const MathVector<dim> Source_g,
                                              const number interface_value);

    public:
    ///    register a update function for a Geometry
        template <typename TFVGeom, typename TAssFunc>
        void register_update_func(TAssFunc func);

    protected:
    ///    Vector holding all update functions
        std::vector<UpdateFunc> m_vUpdateFunc;
};


/// creates upwind based on a string identifier
template <int dim>
SmartPtr<INavierStokesPressureJump<dim> > CreateNavierStokesPressureJump(const std::string& name);

/////////////////////////////////////////////////////////////////////////////
// PressureJump
/////////////////////////////////////////////////////////////////////////////


template <int TDim>
class NavierStokesViscousPressureJump
    : public INavierStokesPressureJump<TDim>
{
    public:
    ///    Base class
        typedef INavierStokesPressureJump<TDim> base_type;

    ///    This class
        typedef NavierStokesViscousPressureJump<TDim> this_type;

    ///    Dimension
        static const int dim = TDim;
	//    abbreviation for pressure
		static const size_t _P_ = dim;
	//    abbreviation for VolFraction
		static const size_t _C_ = dim+1;

    protected:
    //    explicitly forward some function
        using base_type::register_update_func;
        using base_type::set_pressure_jump_comp_connected;

        using base_type::pressure_shape_vel;
        using base_type::pressure_shape_p;
        using base_type::pressure_jump;
    
        using base_type::tang_vel_shape_vel;
        using base_type::tang_vel;
    
    
        using base_type::diff_length_sq_inv;
    
    //    functions from upwind
        using base_type::upwind_conv_length;
        using base_type::downwind_conv_length;
        using base_type::upwind_shape_sh;
        using base_type::downwind_shape_sh;
        using base_type::non_zero_shape_ip;
        using base_type::upwind_shape_ip;
        using base_type::downwind_shape_ip;
    

        


    public:
    ///    constructor
        NavierStokesViscousPressureJump()
        {
        //    vel comp not coupled
            set_pressure_jump_comp_connected(false);

        //    register evaluation function
            register_func();
        }
    
        //void Viscosity_jump( number& mu_2, number& mu_1, const DataImport<number, dim>& N_value,const DataImport<number, dim>& viscositySCV, const DataImport<number, dim>& densitySCV, const size_t nip);
    

    ///    update of values for FV1Geometry
        template <typename TElem>
        void update(const FV1Geometry<TElem, dim>* geo,
                    const LocalVector& vCornerValue,
                    const MathVector<dim> vStdVel[],
                    const bool bStokes,
                    const DataImport<MathVector<dim>, dim>& n,
                    const DataImport<number, dim>& kinViscoSCV,
                    const DataImport<number, dim>& density,
                    const DataImport<number, dim>& densitySCV,
                    const number JumpShape[],
                    const LocalVector* pvCornerValueOldTime, number dt,
                    const number density_ref,
                    const number mu_l,
                    const number rho_l,
                    const MathVector<dim> Source_l,
                    const number mu_g,
                    const number rho_g,
                    const MathVector<dim> Source_g,
                    const number interface_value);
    


    private:
        void register_func();

        template <typename TElem>
        void register_func()
        {
            typedef FV1Geometry<TElem, dim> TGeom;
            typedef void (this_type::*TFunc)(const TGeom* geo,
                                             const LocalVector& vCornerValue,
                                             const MathVector<dim> vStdVel[],
                                             const bool bStokes,
                                             const DataImport<MathVector<dim>, dim>& n,
                                             const DataImport<number, dim>& kinViscoSCV,
                                             const DataImport<number, dim>& density,
                                             const DataImport<number, dim>& densitySCV,
                                             const number JumpShape[],
                                             const LocalVector* pvCornerValueOldTime, number dt,
                                             const number density_ref,
                                             const number mu_l,
                                             const number rho_l,
                                             const MathVector<dim> Source_l,
                                             const number mu_g,
                                             const number rho_g,
                                             const MathVector<dim> Source_g,
                                             const number interface_value);

            this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
        }
};




} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__PRESSURE__JUMP__MODEL__ */
