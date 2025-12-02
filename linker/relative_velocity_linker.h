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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__RELATIVE_VELOCITY_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__RELATIVE_VELOCITY_LINKER__

#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "../properties_interface.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// Relative Velocity linker
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class RelativeVelocityLinker
    : public StdDataLinker< RelativeVelocityLinker<dim>, MathVector<dim>, dim>
{
    ///    Base class type
        typedef StdDataLinker< RelativeVelocityLinker<dim>, MathVector<dim>, dim> base_type;

    public:
    RelativeVelocityLinker() :
            m_spVolumeFraction(NULL), m_spDVolumeFraction(NULL),
            m_spVolumeGrad(NULL), m_spDVolumeGrad(NULL),
            m_spMixDensity(NULL), m_spDMixDensity(NULL),
            m_spMixKinViscosity(NULL), m_spDMixKinViscosity(NULL),
            m_spPsGrad(NULL), m_spDPsGrad(NULL),
            m_spPressureGrad(NULL), m_spDPressureGrad(NULL),
            m_spEinsVisc(NULL), m_spDEinsVisc(NULL),
            m_partialDerivMask(0),
            m_VolLimit(1.0),
            m_Wr(0.0),
            m_Cd(1.0),
            rho_s(2500),
            rho_a(1.2),
            dp(1e-03),
            mu_a(1.48e-05),
            alpha_max(0.635),
            m_gravitation(-9.81),
            m_BoolRelativeVel(true)
        {
        //    this linker needs exactly five input
            this->set_num_input(7);
        }


        inline void evaluate (MathVector<dim>& value,
                              const MathVector<dim>& globIP,
                              number time, int si) const
        {
            UG_LOG("RelativeVelocityLinker::evaluate single called");
            number volume_fraction;
            MathVector<dim> volume_grad;
            number particle_diameter;
            number EinsVisc;
            
            (*m_spVolumeFraction)(volume_fraction, globIP, time, si);
            (*m_spVolumeGrad)(volume_grad, globIP, time, si);
            (*m_spMixKinViscosity)(particle_diameter, globIP, time, si);
            (*m_spEinsVisc)(EinsVisc, globIP, time, si);
                                                    
                                                                 
                                                          
                                                                   
                                                    
        //    Variables
            MathVector<dim> Vel;VecSet(Vel,0.0);



            VecScale(value, Vel,0);
        }

        template <int refDim>
        inline void evaluate(MathVector<dim> vRelativeVelocity[],
                             const MathVector<dim> vGlobIP[],
                             number time, int si,
                             GridObject* elem,
                             const MathVector<dim> vCornerCoords[],
                             const MathVector<refDim> vLocIP[],
                             const size_t nip,
                             LocalVector* u,
                             const MathMatrix<refDim, dim>* vJT = NULL) const
        {
            std::vector<number> vVolume(nip);
            std::vector<number> vMixDensity(nip);
            std::vector<number> vMixKinVisc(nip);
            std::vector<MathVector<dim>> vPsGrad(nip);
            std::vector<number> vEinsVisc(nip);
            std::vector<MathVector<dim>> vGravityForce(nip);

            (*m_spVolumeFraction)(&vVolume[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spMixKinViscosity)(&vMixKinVisc[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spMixDensity)(&vMixDensity[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spPsGrad)(&vPsGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spEinsVisc)(&vEinsVisc[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);

            
            

            
            /*bool m_cut_element_scvf=false;
            size_t numSH=0;
            number interface=0.5;
            if (vGlobIP!=NULL && vLocIP!=NULL && vCornerCoords!=NULL)
                cut_element(m_cut_element_scvf,numSH ,u, elem, vCornerCoords, vGlobIP, interface);*/
            
            
            MathVector<dim> W;
            
            const number Dim=vRelativeVelocity[0].size();

            for(size_t ip = 0; ip < nip; ++ip)
            {
                VecSet(W, 0.0);
                number RelVel = 0.0;
                if(m_BoolRelativeVel )
                {
                    RelativeVel(RelVel, vPsGrad[ip][Dim-1], m_gravitation, vVolume[ip], vMixDensity[ip], vMixKinVisc[ip], vEinsVisc[ip], m_Wr, m_Cd, rho_a, rho_s, dp, mu_a, alpha_max, Inter,false);
					W[Dim-1] = RelVel;
                }
                
                vRelativeVelocity[ip] = W;
                
            }
        }

        template <int refDim>
        void eval_and_deriv(MathVector<dim> vRelativeVelocity[],
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
        //    get the data of the ip series

        int s_VOL_ = base_type::series_id(_VOL_, s);
        //int s_V_ = base_type::series_id(_DVOL_, s);
        int s_RHO_ = base_type::series_id(_RHO_, s);
        int s_MU_ = base_type::series_id(_MU_, s);
        int s_DPS_  = base_type::series_id(_DPS_ , s);
        int s_EMU_  = base_type::series_id(_EMU_ , s);
        
        const number* vVolume = m_spVolumeFraction->values(s_VOL_);
        //const MathVector<dim>* vVolumeGrad = m_spVolumeGrad->values(s);
        const number* vMixDensity = m_spMixDensity->values(s_RHO_);
        const number* vMixKinVisc = m_spDMixKinViscosity->values(s_MU_);
        const MathVector<dim>* vPsGrad = m_spPsGrad->values(s_DPS_);
        const number* vEinsVisc = m_spEinsVisc->values(s_EMU_);
        //const MathVector<dim>* vPressureGrad = m_spPressureGrad->values(s_DP_);
        
        

        
        /*bool m_cut_element_scvf=false;
        size_t numSH=0;
        number interface=0.5;
        if (vGlobIP!=NULL && vLocIP!=NULL && vCornerCoords!=NULL)
            cut_element(m_cut_element_scvf,numSH ,u, elem, vCornerCoords, vGlobIP, interface);*/
        
        
        MathVector<dim> W;
        
        const size_t Dim=vRelativeVelocity[0].size();

        for(size_t ip = 0; ip < nip; ++ip)
        {
            VecSet(W, 0.0);
            number RelVel = 0.0;
            if(m_BoolRelativeVel )
            {
                RelativeVel(RelVel, vPsGrad[ip][Dim-1], m_gravitation, vVolume[ip], vMixDensity[ip], vMixKinVisc[ip], vEinsVisc[ip], m_Wr, m_Cd, rho_a, rho_s, dp , mu_a, alpha_max, Inter, true);
                W[dim-1] = RelVel;
            }
            
            vRelativeVelocity[ip] = W;

            
        }

        /*bool m_scvf=false;
        bool m_scv=false;
        number VolFrac_surface;
        //if (!(vGlobIP!=NULL && vLocIP!=NULL && vCornerCoords!=NULL)) UG_THROW("Error in Granular Source");
        
        Inter->integration_points(m_scvf, m_scv, elem,   vCornerCoords,   vGlobIP, nip);
        bool f = m_scv;
        if(m_scv)
            for(size_t ip = 0; ip < nip; ++ip)
            {

                if (!((vGlobIP[ip][0] > 23.495) && (vGlobIP[ip][0] < 24.005) && (vGlobIP[ip][1] > 1.324) && (vGlobIP[ip][1] < 1.778) ))
                {
                    f = f & false;
                    
                }
            }
        if(f)
        {
            printf("Relative velocity Nip = %zu\n",nip);
            
            printf("Coordinates[0] = %f     %f\n",vCornerCoords[0][0],vCornerCoords[0][1]);
            printf("Coordinates[1] = %f     %f\n",vCornerCoords[1][0],vCornerCoords[1][1]);
            printf("Coordinates[2] = %f     %f\n",vCornerCoords[2][0],vCornerCoords[2][1]);
            //printf("Coordinates[3] = %f     %f\n",vCornerCoords[3][0],vCornerCoords[3][1]);
            
            printf("vGlobIP[0] = %f     %f\n",vGlobIP[0][0],vGlobIP[0][1]);
            printf("vGlobIP[1] = %f     %f\n",vGlobIP[1][0],vGlobIP[1][1]);
            printf("vGlobIP[2] = %f     %f\n",vGlobIP[2][0],vGlobIP[2][1]);
            //printf("vGlobIP[3] = %f     %f\n",vGlobIP[3][0],vGlobIP[3][1]);
            
            printf("vVolumeFraction[0] =    %f     \n",vVolume[0]);
            printf("vVolumeFraction[1] =    %f     \n",vVolume[1]);
            printf("vVolumeFraction[2] =    %f     \n",vVolume[2]);
            //printf("vVolumeFraction[3] =    %f     \n",vVolume[3]);
             
            printf("vRelativeVelocity[0] = %f     %f\n",vRelativeVelocity[0][0],vRelativeVelocity[0][1]);
            printf("vRelativeVelocity[1] = %f     %f\n",vRelativeVelocity[1][0],vRelativeVelocity[1][1]);
            printf("vRelativeVelocity[2] = %f     %f\n",vRelativeVelocity[2][0],vRelativeVelocity[2][1]);
            //printf("vRelativeVelocity[3] = %f     %f\n",vRelativeVelocity[3][0],vRelativeVelocity[3][1]);
             
            printf("vPsGrad[0] = %f     %f\n",vPsGrad[0][0],vPsGrad[0][1]);
            printf("vPsGrad[1] = %f     %f\n",vPsGrad[1][0],vPsGrad[1][1]);
            printf("vPsGrad[2] = %f     %f\n",vPsGrad[2][0],vPsGrad[2][1]);
            //printf("vPsGrad[3] = %f     %f\n",vPsGrad[3][0],vPsGrad[3][1]);
             
            printf("vGravityForce[0] = %f     %f\n",0.0,m_gravitation*vMixDensity[0]);
            printf("vGravityForce[1] = %f     %f\n",0.0,m_gravitation*vMixDensity[1]);
            printf("vGravityForce[2] = %f     %f\n",0.0,m_gravitation*vMixDensity[2]);
            //printf("vGravityForce[3] = %f     %f\n",vGravityForce[3][0],vGravityForce[3][1]);
            
            
            //printf("F[0] =    %f     \n",Fy[0]);
            //printf("F[1] =    %f     \n",Fy[1]);
            //printf("F[2] =    %f     \n",Fy[2]);
            
            if (m_scv) printf("SCV_____\n");
            if (m_scvf) printf("SCVF_____\n");
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
                        if(m_BoolRelativeVel && vVolume[ip]<=m_VolLimit &&  vol_grad<=limit_grad)
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
        static void RelativeVel(number& RelVel, const number vPsGrad, const number m_gravitation, const number vVolume, const number vMixDensity, const number vMixKinVisc, const number vEinstVisc, const number m_Wr, const number m_Cd, const number rho_a, const number rho_s, const number dp, const number mu_a, const number alpha_max, Interface<dim>* Inter, bool EvalAndDeriv)
        {
            const number Fy = -vPsGrad+m_gravitation * vMixDensity;
            number W = 0.0;
            number W1 = 0.0;
            //number W2 = 0.0;
            //number W3 = 0.0;
            //number W4 = 0.0;
            //number W5 = 0.0;
            //number W6 = 0.0;
            //number W7 = 0.0;
            //number Re,Cd;
            //const number phi=fmin(alpha_max-1e-04, fmax(vVolume,0.0));
            //Ratio = (1-phi)/((1+pow(phi,1/3))*exp(5*phi/(3*(1-phi))))-vol*0.017811336463534444;
            //const number Ratio=(1.0-phi)/((1.0+pow(phi,1.0/3.0))*exp(5.0*phi/(3.0*(1.0-phi))));
            const number MU2 = vMixDensity*vMixKinVisc;
            //Cd_phi=Cd;//*pow(Ratio,2)*(1-phi)/(1+pow(phi,1/3));
            //vol_grad=sqrt(VecProd(vVolumeGrad[ip], vVolumeGrad[ip]));
            bool ff= false;
            if( Fy < 0.0 )
            {
                
                W1 = m_Wr*Ratio;
                size_t iter = 0;
                Inter->RelVel(W1,  iter,      MU2,   vMixDensity, dp, rho_s,  fabs(Fy/vMixDensity), 1.0);
                //Inter->RelVel(W2,  iter,      MU2,   vMixDensity, dp, rho_s,  fabs(Fy/vMixDensity), 1e-03);
                //Re = Inter->RE( MU2, vMixDensity, dp,  W2);
                //Cd = Inter->CD( Re,  Inter->DragModel());
                //W2 = (rho_s-rho_a) * pow(dp,2.0) * fabs(Fy/vMixDensity) / (18.0 * MU2);
                //W3 = (rho_s-vMixDensity) * pow(dp,2.0) * fabs(Fy/vMixDensity) / (18.0 * MU2);
                
                //W4 = (rho_s-rho_a) * pow(dp,2.0) * fabs(Fy/vMixDensity) / (18.0 * vEinstVisc);
                //W5 = (rho_s-vMixDensity) * pow(dp,2.0) * fabs(Fy/vMixDensity) / (18.0 * vEinstVisc);
                
                
                //Inter->RelVel(W7,  iter, vEinstVisc, vMixDensity, dp, rho_s,  fabs(Fy/vMixDensity), 1e-03);
                //number W8 = 0.0;
                //Inter->RelVel(W8,  iter, mu_a, rho_a, dp, rho_s,  9.81, 1e-03);
                //printf("Vel = %f\n",W8);
                //W6=sqrt((4.0/3.0)*dp*(rho_s/vMixDensity-1.0)*fabs(Fy)/m_Cd);
                //W7=sqrt((4.0/3.0)*dp*(rho_s/vMixDensity-1.0)*fabs(Fy)/m_Cd);
                
                W = -W1;//fmin(W1,W2);
                //if(phi>0.1) ff= true;
                //VecScale(W, vol,W);
                
            }
            /*else
            {
                W1 = m_Wr*Ratio;
                W2 = (rho_s-rho_a) * pow(dp,2.0) * Fy / (18.0 * MU2);
                W3 = Ratio*(rho_s-vMixDensity) * pow(dp,2.0) * Fy / (18.0 * vEinstVisc);
                W4=sqrt((4.0/3.0)*dp*(rho_s/rho_a-1.0)*Fy/m_Cd);
                W = fmin(W1,W3);
                ff= true;
                
            }*/
            if(ff && EvalAndDeriv)
            {
                printf("W1 = %f,  W2 = %f, Log2 = %f,  W3 = %f, Log3 = %f, phi = %f,  Ratio = %f\n", W1, W2, log10(fabs(W2)), W3, log10(fabs(W3)), phi, Ratio);
            }
            
            RelVel=W;
            
            
            if(std::isnan(RelVel)) UG_THROW("Error in  RelVelLinker: Value = NaN" <<"  Ws = "<<RelVel<<".");
            
            /*if (vol_grad>=limit_grad &&  phi>=limit_vol)
             {
             number pi=3.14159265359;
             number theta_c=34*pi/180;
             MathVector<dim> n; VecScale(n, vVolumeGrad[ip],-1/vol_grad);
             number r_s=0;
             for(size_t n_i = 0; n_i < Dim-1; ++n_i)
             r_s+=pow(n[n_i],2);
             r_s=sqrt(r_s);
             number theta=pi/2-acos(r_s);
             
             if (theta>theta_c)
             {
             VecScale(n, n, cos(theta_c)/r_s);
             n[Dim-1]=-sin(theta_c);
             VecScale(vRelativeVelocity[ip], n, w*Ratio);
             }
             }*/
        
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

    ///    set density import
        void set_mix_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data)
        {
            m_spMixKinViscosity = data;
            m_spDMixKinViscosity = data.template cast_dynamic<DependentUserData<number, dim> >();
            base_type::set_input(_MU_, data, data);
        }

        void set_mix_kinematic_viscosity(number val)
        {
            set_mix_kinematic_viscosity(make_sp(new ConstUserNumber<dim>(val)));
        }

    ///    set Particle Pressure Grad import
        void set_ps_grad(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
        {
            m_spPsGrad = data;
            m_spDPsGrad = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim> >();
            base_type::set_input(_DPS_, data, data);
        }
    ///    set Pressure Grad import
        void set_pressure_grad(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
        {
            m_spPressureGrad = data;
            m_spDPressureGrad = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim> >();
            base_type::set_input(_DP_, data, data);
        }
    ///    set density import
        void set_einstein_visc(SmartPtr<CplUserData<number, dim> > data)
        {
            m_spEinsVisc = data;
            m_spDEinsVisc = data.template cast_dynamic<DependentUserData<number, dim> >();
            base_type::set_input(_EMU_, data, data);
        }

        void set_einstein_visc(number val)
        {
            set_einstein_visc(make_sp(new ConstUserNumber<dim>(val)));
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
    
    ///    import for density
        static const size_t _RHO_ = 2;
        SmartPtr<CplUserData<number, dim> > m_spMixDensity;
        SmartPtr<DependentUserData<number, dim> > m_spDMixDensity;

    ///    import for density
        static const size_t _MU_ = 3;
        SmartPtr<CplUserData<number, dim> > m_spMixKinViscosity;
        SmartPtr<DependentUserData<number, dim> > m_spDMixKinViscosity;
    
    ///    import for Particle Pressure Gradient
        static const size_t _DPS_ = 4;
        SmartPtr<CplUserData<MathVector<dim>, dim> > m_spPsGrad;
        SmartPtr<DependentUserData<MathVector<dim>, dim> > m_spDPsGrad;
    
    ///    import for Pressure Gradient
        static const size_t _DP_ = 5;
        SmartPtr<CplUserData<MathVector<dim>, dim> > m_spPressureGrad;
        SmartPtr<DependentUserData<MathVector<dim>, dim> > m_spDPressureGrad;
    
        static const size_t _EMU_ = 6;
        SmartPtr<CplUserData<number, dim> > m_spEinsVisc;
        SmartPtr<DependentUserData<number, dim> > m_spDEinsVisc;
    
    
        Interface<dim>* Inter;


    public:
        void set_derivative_mask(int mask) {
            m_partialDerivMask = mask;
            std::cerr << "Setting some derivatives: "<< m_partialDerivMask << "(" << this <<")" << std::endl;
        }
        void set_particle_density(float R) {
            rho_s = R;
        }
        void set_fluid_density(float R) {
            rho_a = R;
        }
        void set_particle_diameter(float R) {
            dp = R;
        }
        void set_fluid_viscosity(float mu) {
            mu_a = mu;
        }
        void set_alpha_max(float alpha) {
            alpha_max = alpha;
        }
        void set_gravity(float R) {
            m_gravitation = R;
        }
        void set_vol_limit(float VolLimit) {
            m_VolLimit = VolLimit;
        }
        void set_rel_vel(float w) {
            m_Wr = w;
        }
        void set_dragCoeff(float Cd) {
            m_Cd = Cd;
        }
        void activate_relative_vel(bool R) {
            m_BoolRelativeVel = R;
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
        float m_GradLimit, m_VolLimit, m_Wr, m_Cd, rho_s, rho_a, dp, mu_a, alpha_max, m_gravitation;
        bool m_BoolRelativeVel;
    


};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__RELATIVE_VELOCITY_LINKER__ */

