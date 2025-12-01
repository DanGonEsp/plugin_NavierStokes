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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_SOURCE_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_SOURCE_LINKER__

#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "../properties_interface.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// Granular Source linker
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class GranularSourceLinker
    : public StdDataLinker< GranularSourceLinker<dim>, MathVector<dim>, dim>
{
    ///    Base class type
        typedef StdDataLinker< GranularSourceLinker<dim>, MathVector<dim>, dim> base_type;

    public:
    GranularSourceLinker() :
            m_spMixDensity(NULL), m_spDMixDensity(NULL),
            rho_a(1.2),
            rho_s(1000.0),
            //packing_factor(1.0),
            m_gravity(-9.81),
            m_bConstGravity(false)
        {
        //    this linker needs exactly 1 input
            this->set_num_input(1);
        }


        inline void evaluate (MathVector<dim>& value,
                              const MathVector<dim>& globIP,
                              number time, int si) const
        {
            UG_LOG("GranularSourceLinker::evaluate single called");
            number particle_diameter;
            
            (*m_spMixDensity)(particle_diameter, globIP, time, si);
            
                                                    
                                                                 
                                                          
                                                                   
                                                    
        //    Variables
            MathVector<dim> Vel;VecSet(Vel,0.0);
            Vel[dim-1] = m_gravity;



            VecScale(value, Vel,0);
        }

        template <int refDim>
        inline void evaluate(MathVector<dim> vGranularSource[],
                             const MathVector<dim> vGlobIP[],
                             number time, int si,
                             GridObject* elem,
                             const MathVector<dim> vCornerCoords[],
                             const MathVector<refDim> vLocIP[],
                             const size_t nip,
                             LocalVector* u,
                             const MathMatrix<refDim, dim>* vJT = NULL) const
        {
            bool m_scvf=false;
            bool m_scv=false;

            if (!(vGlobIP!=NULL && vLocIP!=NULL && vCornerCoords!=NULL) && nip>0) UG_THROW("Error in Granular Source");
            Inter->integration_points(m_scvf, m_scv, elem,   vCornerCoords,   vGlobIP, nip);
            
            DimFV1Geometry<dim> geo;
            geo.update(elem, vCornerCoords, NULL);
            
            if( m_bConstGravity  )
            {
                const size_t numSh = geo.num_scv();
                std::vector<number> vMixDensitySCV(geo.num_scv());
                (*m_spMixDensity)(&vMixDensitySCV[0], geo.scv_global_ips(), time, si,
                                  elem, vCornerCoords, geo.scv_local_ips(), numSh, u, NULL);
                number DensitySCV[numSh];
                for(size_t ip = 0; ip < geo.num_scv(); ++ip)
                {
                    DensitySCV[ip] = vMixDensitySCV[ip];
                }
                
                StdLinConsistentGravityX<refDim> RhoG;
                MathVector<refDim> vConsGravity[numSh];
                MathVector<dim> Gravity; VecSet(Gravity,0.0); Gravity[dim-1] = m_gravity;
                RhoG.template prepare<dim>(vConsGravity, numSh, vCornerCoords, DensitySCV, Gravity);
                if(m_scvf)
                {
                    for(size_t ip = 0; ip < nip; ++ip)
                    {
                        //     get current SCVF
                        const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
                        MathVector<dim> vRhoGravity; VecSet(vRhoGravity,0.0);
                        
                        MathVector<refDim> LocalCoord_aux;
                        for(int d = 0; d < refDim; ++d)
                            LocalCoord_aux[d]=scvf.local_ip()[d];
                        
                        MathMatrix<dim, refDim> JTInv_aux;
                        for(int d = 0; d <dim ; ++d)
                            for(int d1 = 0; d1 < refDim; ++d1)
                                JTInv_aux(d,d1)=scvf.JTInv()(d,d1);
                        
                        MathVector<refDim> vLocalGrad[numSh];
                        for(size_t sh = 0; sh < numSh; sh++)
                            for(size_t d = 0; d < refDim; d++)
                                vLocalGrad[sh][d] = scvf.local_grad_vector()[sh][d];
                        
                        RhoG.template compute<dim>(vRhoGravity,LocalCoord_aux, JTInv_aux ,vLocalGrad, vConsGravity);
                        
                        vGranularSource[ip] = vRhoGravity;
                    }
                    
                    
                }
                else if(m_scv)
                {
                    for(size_t ip = 0; ip < nip; ++ip)
                    {
                        //     get current SCV
                        const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(ip);
                        MathVector<dim> vRhoGravity; VecSet(vRhoGravity,0.0);
                        
                        
                        MathVector<refDim> LocalCoord_aux;
                        for(int d = 0; d < refDim; ++d)
                            LocalCoord_aux[d]=scv.local_ip()[d];
                        
                        MathMatrix<dim, refDim> JTInv_aux;
                        for(int d = 0; d <dim ; ++d)
                            for(int d1 = 0; d1 < refDim; ++d1)
                                JTInv_aux(d,d1)=scv.JTInv()(d,d1);
                        
                        MathVector<refDim> vLocalGrad[numSh];
                        for(size_t sh = 0; sh < numSh; sh++)
                            for(size_t d = 0; d < refDim; d++)
                                vLocalGrad[sh][d] = scv.local_grad_vector()[sh][d];
                        
                        
                        RhoG.template compute<dim>(vRhoGravity,LocalCoord_aux, JTInv_aux,vLocalGrad, vConsGravity);
                        
                        vGranularSource[ip] = vRhoGravity;
                    }
                }
                else 
                {
                    
                    ReferenceObjectID roid = elem->reference_object_id();
                    const LocalShapeFunctionSet<refDim>& rTrialSpace =
                            LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));
                    
                    MathVector<refDim> vLocGrad[numSh];
                    //MathVector<refDim> locGrad;
                    //MathVector<refDim> vGlobGrad;
                    
                    MathMatrix<dim, refDim> JTInv;
                    //    get Reference Mapping
                    DimReferenceMapping<refDim, dim>& map = ReferenceMappingProvider::get<refDim, dim>(roid, vCornerCoords);
                    
                    
                    for(size_t ip = 0; ip < nip; ++ip)
                    {
                        
                        rTrialSpace.grads(vLocGrad, vLocIP[ip]);
                        
                        //    compute global grad
                        map.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
                        
                        MathVector<dim> vRhoGravity; VecSet(vRhoGravity,0.0);
                        
                        
                        RhoG.template compute<dim>(vRhoGravity,vLocIP[ip], JTInv ,vLocGrad, vConsGravity);
                        
                        vGranularSource[ip] = vRhoGravity;
                        
                        /*if(vRhoGravity[1]<-12 || fabs(vRhoGravity[1]-vRhoGravity2[1])>1e-03)
                        {
                            printf("vRhoGravity   = [%f  %f]\n",vRhoGravity[0],vRhoGravity[1]);
                            printf("vRhoGravity2  = [%f  %f]\n",vRhoGravity2[0],vRhoGravity2[1]);
                            
                            
                            printf("JTInv_aux  = \n");
                            printf("   %f     %f    \n", JTInv_aux(0, 0),JTInv_aux(0, 1));
                            printf("   %f     %f    \n", JTInv_aux(1, 0),JTInv_aux(1, 1));
                            
                            
                            printf("JTInv  = \n");
                            printf("   %f     %f    \n", JTInv(0, 0),JTInv(0, 1));
                            printf("   %f     %f    \n", JTInv(1, 0),JTInv(1, 1));
                        }*/
                        
                        
                    }
                    

                    
                    
                    
                }
                
                
                
            }
            else
            {
                
                std::vector<number> vMixDensity(nip);
                (*m_spMixDensity)(&vMixDensity[0], vGlobIP, time, si,
                                elem, vCornerCoords, vLocIP, nip, u, vJT);
                                
                MathVector<dim> W;
                for(size_t ip = 0; ip < nip; ++ip)
                {
                    VecSet(W,0);
                    W[dim-1] = m_gravity;
                    VecScale(vGranularSource[ip],W,vMixDensity[ip]);
                }


            }
            
        }

        template <int refDim>
        void eval_and_deriv(MathVector<dim> vGranularSource[],
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
        
        
        bool m_scvf=false;
        bool m_scv=false;
        
        if (!(vGlobIP!=NULL && vLocIP!=NULL && vCornerCoords!=NULL) && nip>0) UG_THROW("Error in Granular Source");
        Inter->integration_points(m_scvf, m_scv, elem,   vCornerCoords,   vGlobIP, nip);
        
        DimFV1Geometry<dim> geo;
        geo.update(elem, vCornerCoords, NULL);
        
        
        
        if( m_bConstGravity  )
        {
            const size_t numSh = geo.num_scv();
            std::vector<number> vMixDensitySCV(geo.num_scv());
            (*m_spMixDensity)(&vMixDensitySCV[0], geo.scv_global_ips(), time, si,
                              elem, vCornerCoords, geo.scv_local_ips(), numSh, u, NULL);
            number DensitySCV[numSh];
            for(size_t ip = 0; ip < geo.num_scv(); ++ip)
            {
                DensitySCV[ip] = vMixDensitySCV[ip];
            }
            
            StdLinConsistentGravityX<refDim> RhoG;
            MathVector<refDim> vConsGravity[numSh];
            MathVector<dim> Gravity; VecSet(Gravity,0.0); Gravity[dim-1] = m_gravity;
            RhoG.template prepare<dim>(vConsGravity, numSh, vCornerCoords, DensitySCV, Gravity);
            if(m_scvf)
            {
                for(size_t ip = 0; ip < nip; ++ip)
                {
                    //     get current SCVF
                    const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
                    MathVector<dim> vRhoGravity; VecSet(vRhoGravity,0.0);
                    
                    MathVector<refDim> LocalCoord_aux;
                    for(int d = 0; d < refDim; ++d)
                        LocalCoord_aux[d]=scvf.local_ip()[d];
                    
                    MathMatrix<dim, refDim> JTInv_aux;
                    for(int d = 0; d <dim ; ++d)
                        for(int d1 = 0; d1 < refDim; ++d1)
                            JTInv_aux(d,d1)=scvf.JTInv()(d,d1);
                    
                    MathVector<refDim> vLocalGrad[numSh];
                    for(size_t sh = 0; sh < numSh; sh++)
                        for(size_t d = 0; d < refDim; d++)
                            vLocalGrad[sh][d] = scvf.local_grad_vector()[sh][d];
                    
                    RhoG.template compute<dim>(vRhoGravity,LocalCoord_aux, JTInv_aux ,vLocalGrad, vConsGravity);
                    
                    vGranularSource[ip] = vRhoGravity;
                }
                
                
            }
            else if(m_scv)
            {
                for(size_t ip = 0; ip < nip; ++ip)
                {
                    //     get current SCV
                    const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(ip);
                    MathVector<dim> vRhoGravity; VecSet(vRhoGravity,0.0);
                    
                    
                    MathVector<refDim> LocalCoord_aux;
                    for(int d = 0; d < refDim; ++d)
                        LocalCoord_aux[d]=scv.local_ip()[d];
                    
                    MathMatrix<dim, refDim> JTInv_aux;
                    for(int d = 0; d <dim ; ++d)
                        for(int d1 = 0; d1 < refDim; ++d1)
                            JTInv_aux(d,d1)=scv.JTInv()(d,d1);
                    
                    MathVector<refDim> vLocalGrad[numSh];
                    for(size_t sh = 0; sh < numSh; sh++)
                        for(size_t d = 0; d < refDim; d++)
                            vLocalGrad[sh][d] = scv.local_grad_vector()[sh][d];
                    
                    
                    RhoG.template compute<dim>(vRhoGravity,LocalCoord_aux, JTInv_aux,vLocalGrad, vConsGravity);
                    
                    vGranularSource[ip] = vRhoGravity;
                }
            }
            else
            {
                
                ReferenceObjectID roid = elem->reference_object_id();
                
                const LocalShapeFunctionSet<refDim>& rTrialSpace =
                        LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));
                
                MathVector<refDim> vLocGrad[numSh];
                //MathVector<refDim> locGrad;
                //MathVector<refDim> vGlobGrad;
                
                MathMatrix<dim, refDim> JTInv;
                //    get Reference Mapping
                DimReferenceMapping<refDim, dim>& map = ReferenceMappingProvider::get<refDim, dim>(roid, vCornerCoords);
                
                
                for(size_t ip = 0; ip < nip; ++ip)
                {
                    
                    rTrialSpace.grads(vLocGrad, vLocIP[ip]);
                    
                    //    compute global grad
                    map.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
                    
                    MathVector<dim> vRhoGravity; VecSet(vRhoGravity,0.0);
                    
                    
                    RhoG.template compute<dim>(vRhoGravity,vLocIP[ip], JTInv ,vLocGrad, vConsGravity);
                    
                    vGranularSource[ip] = vRhoGravity;
                    
                    /*if(vRhoGravity[1]<-12 || fabs(vRhoGravity[1]-vRhoGravity2[1])>1e-03)
                    {
                        printf("vRhoGravity   = [%f  %f]\n",vRhoGravity[0],vRhoGravity[1]);
                        printf("vRhoGravity2  = [%f  %f]\n",vRhoGravity2[0],vRhoGravity2[1]);
                        
                        
                        printf("JTInv_aux  = \n");
                        printf("   %f     %f    \n", JTInv_aux(0, 0),JTInv_aux(0, 1));
                        printf("   %f     %f    \n", JTInv_aux(1, 0),JTInv_aux(1, 1));
                        
                        
                        printf("JTInv  = \n");
                        printf("   %f     %f    \n", JTInv(0, 0),JTInv(0, 1));
                        printf("   %f     %f    \n", JTInv(1, 0),JTInv(1, 1));
                    }*/
                    
                    
                }
                

                
                
                
            }
            
            
            
        }
        else
        {   
            
            
            int s_RHO_ = base_type::series_id(_RHO_, s);
            //int s_P_   = base_type::series_id(_P_  , s);
            
            //    get the data of the ip series
            const number* vMixDensity = m_spMixDensity->values(s_RHO_);
            
            MathVector<dim> W;
            for(size_t ip = 0; ip < nip; ++ip)
            {
                VecSet(W,0);
                W[dim-1] = m_gravity;
                VecScale(vGranularSource[ip],W,vMixDensity[ip]);
            }


        }
        
    

        
        
        
        
        /*
        else
        {
            for(size_t ip = 0; ip < nip; ++ip)
                VecSet(vGranularSource[ip],0);
        }*/


        
        //    Compute the derivatives at all ips     //
        /////////////////////////////////////////////
        
        //    check if something to do
        if(!bDeriv || this->zero_derivative()) return;
        
        //    clear all derivative values
        this->set_zero(vvvDeriv, nip);
        /*
         //    Derivatives of Viscosity
         if(m_partialDerivMask == 0 && m_spDVolumeFraction.valid() && !m_spDVolumeFraction->zero_derivative())
         for(size_t ip = 0; ip < nip; ++ip)
         for(size_t fct = 0; fct < m_spDVolumeFraction->num_fct(); ++fct)
         {
         //    get derivative of viscosity w.r.t. to all functions
         const number* vDVolumeFraction = m_spDVolumeFraction->deriv(s, ip, fct);
         
         //    get common fct id for this function
         const size_t commonFct = this->input_common_fct(_VOL_, fct);
         
         //    loop all shapes and set the derivative
         for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
         {
         //  DarcyVel_fct[sh] -= mu_fct_sh / mu * q
         VecScaleAppend(vvvDeriv[ip][commonFct][sh], -vDViscosity[sh] / vViscosity[ip], vDarcyVel[ip]);
         }
         }
         
         //    Derivatives of Density
         if( m_partialDerivMask == 0 && m_spDDensity.valid() && !m_spDDensity->zero_derivative())
         for(size_t ip = 0; ip < nip; ++ip)
         for(size_t fct = 0; fct < m_spDDensity->num_fct(); ++fct)
         {
         //    get derivative of viscosity w.r.t. to all functions
         const number* vDDensity = m_spDDensity->deriv(s, ip, fct);
         
         //    get common fct id for this function
         const size_t commonFct = this->input_common_fct(_RHO_, fct);
         
         //    Precompute K/mu * g
         MathVector<dim> Kmug;
         
         //    a) compute K * g
         MatVecMult(Kmug, vPermeability[ip], vGravity[ip]);
         
         //    b) compute K* g / mu
         VecScale(Kmug, Kmug, 1./vViscosity[ip]);
         
         //    loop all shapes and set the derivative
         for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
         {
         UG_ASSERT(commonFct < vvvDeriv[ip].size(), commonFct<<", "<<vvvDeriv[ip].size());
         UG_ASSERT(sh < vvvDeriv[ip][commonFct].size(), sh<<", "<<vvvDeriv[ip][commonFct].size());
         //  DarcyVel_fct[sh] += K/mu * (rho_fct_sh * g)
         VecScaleAppend(vvvDeriv[ip][commonFct][sh],
         vDDensity[sh], Kmug);
         }
         }
         
        
        
        
        //    Derivatives of VolumeFraction
        if(m_partialDerivMask == 0 && m_spDVolumeFraction.valid() && !m_spDVolumeFraction->zero_derivative())
            for(size_t ip = 0; ip < nip; ++ip){
                
                vol=fmin(fmax(0.0,vVolume[ip]),1.0);
                phi=vPackingFactor[ip]*vol;
                Ratio=vPackingFactor[ip]*(-2*phi*phi-pow(phi,2/3)*(3*phi-8)+6*phi+1)*exp(5*phi/(3*(phi-1)))/(3*pow(1+pow(phi,1/3),2)*(phi-1)*pow(phi,2/3))-0.017811336463534444;
                VecSet(W,0);
                W[Dim-1]+=-vRelVel[ip]*Ratio;
                
                for(size_t fct = 0; fct < m_spDVolumeFraction->num_fct(); ++fct)
                {
                    //    get derivative of viscosity w.r.t. to all functions
                    const number* vDVolumeFraction = m_spDVolumeFraction->deriv(s, ip, fct);
                    
                    //    get common fct id for this function
                    const size_t commonFct = this->input_common_fct(_VOL_, fct);
                    
                    //    loop all shapes and set the derivative
                    for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                    {
                        if(vVolume[ip]<=m_VolLimit &&  vol_grad<=limit_grad)
                            VecScaleAppend(vvvDeriv[ip][commonFct][sh], vDVolumeFraction[sh], W);
                        else
                            VecScaleAppend(vvvDeriv[ip][commonFct][sh], 0, W);
    
                    }
                }
            }
            
         
        //    Derivatives of Pressure
            if(m_partialDerivMask ==0 && m_spDPressureGrad.valid() && !m_spDPressureGrad->zero_derivative()  )
            for(size_t ip = 0; ip < nip; ++ip)
                for(size_t fct = 0; fct < m_spDPressureGrad->num_fct(); ++fct)
                {
                //    get derivative of viscosity w.r.t. to all functions
                    const MathVector<dim>* vDPressureGrad = m_spDPressureGrad->deriv(s, ip, fct);

                //    get common fct id for this function
                    const size_t commonFct = this->input_common_fct(_DP_, fct);

                //    Precompute -K/mu
                    MathMatrix<dim,dim> Kmu;

                //    a) compute -K/mu
                    MatScale(Kmu, -1.0/vViscosity[ip],vPermeability[ip]);

                //    loop all shapes and set the derivative
                    for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                    {
                        MathVector<dim> tmp;
                        MatVecMult(tmp, Kmu, vDPressureGrad[sh]);

                        vvvDeriv[ip][commonFct][sh] += tmp;
                    }
                }

        //    Derivatives of Permeability
            if(m_spDPermeability.valid() && !m_spDPermeability->zero_derivative())
            for(size_t ip = 0; ip < nip; ++ip)
                for(size_t fct = 0; fct < m_spDPermeability->num_fct(); ++fct)
                {
                //    get derivative of viscosity w.r.t. to all functions
                    const MathMatrix<dim,dim>* vDPermeability = m_spDPermeability->deriv(s, ip, fct);

                //    get common fct id for this function
                    const size_t commonFct = this->input_common_fct(_K_, fct);

                //    Variables
                    MathVector<dim> Vel;

                //    compute rho*g
                    VecScale(Vel, vGravity[ip], vDensity[ip]);

                //     compute rho*g - \nabla p
                    VecSubtract(Vel, Vel, vPressureGrad[ip]);

                //    compute Darcy velocity q := K / mu * (rho*g - \nabla p)
                    VecScale(Vel, Vel, 1./vViscosity[ip]);

                //    loop all shapes and set the derivative
                    for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                    {
                        MathVector<dim> tmp;
                        MatVecMult(tmp, vDPermeability[sh], Vel);

                        vvvDeriv[ip][commonFct][sh] += tmp;
                    }
                }*/
        }

    public:
    
    ///    set density import
        void set_mix_density(SmartPtr<CplUserData<number, dim> > data)
        {
            m_spMixDensity = data;
            m_spDMixDensity = data.template cast_dynamic<DependentUserData<number, dim> >();
            base_type::set_input(_RHO_, data, data);
        }

        void set_mix_density(number val)
        {
            set_mix_density(make_sp(new ConstUserNumber<dim>(val)));
        }
    ///    set Particle Pressure Grad import
        /*void set_ps_grad(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
        {
            m_spPsGrad = data;
            m_spDPsGrad = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim> >();
            base_type::set_input(_DPS_, data, data);
        }*/

    protected:
    
    
    ///    import for density
        static const size_t _RHO_ = 0;
        SmartPtr<CplUserData<number, dim> > m_spMixDensity;
        SmartPtr<DependentUserData<number, dim> > m_spDMixDensity;
    
    ///    import for Particle Pressure Gradient
        /*static const size_t _DPS_ = 5;
        SmartPtr<CplUserData<MathVector<dim>, dim> > m_spPsGrad;
        SmartPtr<DependentUserData<MathVector<dim>, dim> > m_spDPsGrad;*/
    
        Interface<dim>* Inter;


    public:

        void set_particle_density(float R) {
            rho_s = R;
        }
        void set_fluid_density(float R) {
            rho_a = R;
        }
        void set_gravity(float R) {
            m_gravity = R;
        }
        void set_cons_gravity(bool R) {
            m_bConstGravity = R;
        }
        void set_phase_parameters(Interface<dim>* user)
        {
            if (!user->valid())
                UG_THROW("Interface parameters has not been initialized");
            Inter = user;
        }
    protected:
        // disable certain derivatives
        int m_partialDerivMask;
        float rho_a, rho_s;//, packing_factor;
        float m_gravity;
        bool m_bConstGravity;

};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_SOURCE_LINKER__ */

