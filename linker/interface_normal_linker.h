/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__INTERFACE_NORMAL_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__INTERFACE_NORMAL_LINKER__

#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif


namespace ug{


////////////////////////////////////////////////////////////////////////////////
// NORMAL INTERFACE LINKER
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class InterfaceNormalLinker
    : public StdDataLinker< InterfaceNormalLinker<dim>, MathVector<dim>, dim>
{
    ///    Base class type
        typedef StdDataLinker< InterfaceNormalLinker<dim>, MathVector<dim>, dim> base_type;

    public:
    InterfaceNormalLinker() :
            m_spVolumeFraction(NULL), m_spDVolumeFraction(NULL),
            m_spVolumeGrad(NULL), m_spDVolumeGrad(NULL),
            interface_volume_fraction(0.5)

        {
        //    this linker needs exactly five input
            this->set_num_input(2);
        }


        inline void evaluate (MathVector<dim>& value,
                              const MathVector<dim>& globIP,
                              number time, int si) const
        {
            UG_LOG("InterfaceNormalLinker::evaluate single called");
            number volume_fraction;
            MathVector<dim> volume_grad;
            
            (*m_spVolumeFraction)(volume_fraction, globIP, time, si);
            (*m_spVolumeGrad)(volume_grad, globIP, time, si);

                                                    
                                                                                                                    
            number vol_grad;
            const number eps = 1e-04;

                
            vol_grad=sqrt(VecProd(volume_grad, volume_grad));
            
            MathVector<dim> n;
            
            if (vol_grad> eps)
                VecScale(n, volume_grad,1.0/vol_grad);
            else
            {
                VecSet(n,0.0);
                n[dim-1] = -1.0;
            }
            value=n;

        }

        template <int refDim>
        inline void evaluate(MathVector<dim> vNormal[],
                             const MathVector<dim> vGlobIP[],
                             number time, int si,
                             GridObject* elem,
                             const MathVector<dim> vCornerCoords[],
                             const MathVector<refDim> vLocIP[],
                             const size_t nip,
                             LocalVector* u,
                             const MathMatrix<refDim, dim>* vJT = NULL) const
        {
            UG_LOG("InterfaceNormalLinker::evaluate single called");
            std::vector<number> vVolume(nip);
            std::vector<MathVector<dim>> vVolumeGrad(nip);


            (*m_spVolumeFraction)(&vVolume[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spVolumeGrad)(&vVolumeGrad[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);


            
            number vol_grad;
            const number eps = 1e-04;
            for(size_t ip = 0; ip < nip; ++ip)
            {
                
                vol_grad=sqrt(VecProd(vVolumeGrad[ip], vVolumeGrad[ip]));
                
                MathVector<dim> n;
                
                if (vol_grad> eps)
                    VecScale(n, vVolumeGrad[ip],1.0/vol_grad);
                else
                {
                    VecSet(n,0.0);
                    n[dim-1] = -1.0;
                }
                vNormal[ip]=n;
            }
        }

        template <int refDim>
        void eval_and_deriv(MathVector<dim> vNormal[],
                            const MathVector<dim> vGlobIP[],
                            number time, int si,
                            GridObject* elem,
                            const MathVector<dim> vCornerCoords[],
                            const MathVector<refDim> vLocIP[],
                            const size_t nip,
                            LocalVector* u,
                            bool bDeriv,
                            int s,
                            std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
                            const MathMatrix<refDim, dim>* vJT = NULL) const
    {        
        
        int s_VOL_ = base_type::series_id(_VOL_, s);
        int s_DVOL_ = base_type::series_id(_DVOL_, s);
        
        const number* vVolumeFraction   = m_spVolumeFraction->values(s_VOL_);
        const MathVector<dim>* vVolumeGrad = m_spVolumeGrad->values(s_DVOL_);
        
        std::vector<number> vVolumeFraction2(nip);
        (*m_spVolumeFraction)(&vVolumeFraction2[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
        
        for(size_t ip = 0; ip < nip; ++ip)
        {
            if(fabs(vVolumeFraction[ip]-vVolumeFraction2[ip]) > 1e-05 )
                UG_THROW("Volume fraction Values are not consistent in Granular Viscosity Linker");
                             
        }


        
        number vol_grad;
        const number eps = 1e-04;
        for(size_t ip = 0; ip < nip; ++ip)
        {

            vol_grad=sqrt(VecProd(vVolumeGrad[ip], vVolumeGrad[ip]));

            MathVector<dim> n;
            
            if (vol_grad> eps)
                VecScale(n, vVolumeGrad[ip],1.0/vol_grad);
            else
            {
                VecSet(n,0.0);
                n[dim-1] = -1.0;
            }
            vNormal[ip]=n;
            
            
        }
        
        //    Compute the derivatives at all ips     //
        /////////////////////////////////////////////
        
        //    check if something to do
        if(!bDeriv || this->zero_derivative()) return;
        
        //    clear all derivative values
        this->set_zero(vvvDeriv, nip);
        
    }
    

    public:
    
    ///    set density import
        void set_volume_fraction(SmartPtr<CplUserData<number, dim> > data)
        {
            m_spVolumeFraction = data;
            m_spDVolumeFraction = data.template cast_dynamic<DependentUserData<number, dim> >();
            base_type::set_input(_VOL_, data, data);
        }

        void set_volume_fraction(number val)
        {
            set_volume_fraction(make_sp(new ConstUserNumber<dim>(val)));
        }
    ///    set gravity import
        void set_volume_grad(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
        {
            m_spVolumeGrad = data;
            m_spDVolumeGrad = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim> >();
            base_type::set_input(_DVOL_, data, data);
        }

    protected:
    
    ///    import for density
        static const size_t _VOL_ = 0;
        SmartPtr<CplUserData<number, dim> > m_spVolumeFraction;
        SmartPtr<DependentUserData<number, dim> > m_spDVolumeFraction;
    ///    import for density
        static const size_t _DVOL_ = 1;
        SmartPtr<CplUserData<MathVector<dim>, dim> > m_spVolumeGrad;
        SmartPtr<DependentUserData<MathVector<dim>, dim> > m_spDVolumeGrad;
    
    


    public:

        void set_interface_volume_fraction(float R) {
            interface_volume_fraction = R;
        }

    protected:

        float interface_volume_fraction;


};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__INTERFACE_NORMAL_LINKER__ */

