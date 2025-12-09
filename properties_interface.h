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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC_NV_CUT_ELEMENT_
#define __H__UG__LIB_DISC__SPATIAL_DISC_NV_CUT_ELEMENT_

// extern libraries
#include <cmath>
#include <map>
#include <vector>

// other ug4 modules
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"
#include "consistent_gravity_multiphase.h"

// library intern includes



namespace ug{

template <int TDim, int TWorldDim = TDim>
class Interface
{


	public:
	///	dimension of reference element
		static const int dim = TDim;

	///	dimension of world
		static const int worldDim = TWorldDim;

	public:
	/// construct object and initialize local values and sizes
        Interface():
            m_P0(0.0),
            rho_s(2400.0),
            rho_a(1.2),
            mu_a(1.48e-5),
            dp(0.001),
            Fr(0.05),
            B_phi(1.0),
            alpha_max(0.635),
            alpha_min(0.57),
            deltaGamma(1e-4),
            m_limit(1e6),
            m_dt(1.0),
            drag_model(1),
            m_bParticleGradientForce(false),
            m_bConsistentGravity(false),
            m_init(false)
        {
            
        }

        void Ps(number* ParticlePressure, number* DCParticlePressure, const MathMatrix<dim,dim> VelocityGrad[], const LocalVector& u, const size_t _C_, const size_t numSH, const bool deriv)
        {
            
            for(size_t sh = 0; sh < numSH; ++sh)
            {
                number gamma=0.0;
                const number phi=fmin(1.0, fmax(u(_C_,sh),0));
                
                
                // compute inner sum
                for(int d1 = 0; d1 < dim; ++d1)
                {
                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        gamma += pow(VelocityGrad[sh](d1,d2) + VelocityGrad[sh](d2,d1),2);
                    }
                }
                gamma =sqrt(pow(deltaGamma,2)+(0.5*gamma));
                
                
                //Stokes number
                const number St=gamma*rho_s*pow(dp,2)/mu_a;
                
                //Permanent contact pressure
                const number pff=(phi >= alpha_min) ? Fr *pow(  phi-alpha_min,3) /pow(alpha_max-phi,5) : 0.0;
                //Dynamic pressure
                const number pa = mu_a*(1.0+St)*pow(B_phi*phi/(alpha_max-phi),2)*gamma;
                
                ParticlePressure[sh] = fmax(pff+pa,0.0);
                if(deriv)
                {
                    number dpff =(phi > alpha_min) ? pff *(2*phi+3*alpha_max-5*alpha_min)/((  phi-alpha_min)*(alpha_max-phi)) : 0.0;
                    
                    const number dpa = 2*mu_a*(1.0+St)*pow(B_phi,2.0)*phi*alpha_max*gamma/(pow(alpha_max-phi,3));
                    
                    
                    DCParticlePressure[sh] = dpff + dpa;
                }
            }

        }

        void Einstein_viscosity(number& ss, number& Dss,const number phi, const bool deriv)
        {

            ss= (phi < alpha_max)? mu_a*pow(1.0-phi/alpha_max,-2.5*alpha_max) : m_limit;
            ss=fmin(ss,m_limit);
            if(deriv)
            {
                Dss = (phi < alpha_max)? 2.5*ss/(1.0-phi/alpha_max) : 0.0;
            }
            
            if(std::isnan(ss) || ss<0.0) UG_THROW("Error in Einstein ViscosityLinker: Value = NaN" <<"  Volume Fraction = "<<phi<<".");
            //if(phi > alpha_max) UG_LOG("Phi > phi_max in Einstein Viscosity\n");
        }
    
        number RelVel_ext(const number mu_a1, const number rho_a1, const number dp1, const number rho_s1, const number g1, const number E1)
        {
            number Vel; size_t iter;
            RelVel(Vel, iter,  mu_a1,  rho_a1, dp1,  rho_s1,  g1,  E1);
            number Re = RE(mu_a1,rho_a1,dp1,Vel);
            number cd = CD(Re,  drag_model);
            printf("Ws = %f\n",Vel);
            printf("Cd = %f\n",cd);
            printf("Iter = %zu\n",iter);
            return Vel;
            
        }
        void RelVel(number& Rel, size_t& iter, const number mu_a1, const number rho_a1, const number dp1, const number rho_s1, const number g1, const number E1)
        {
            size_t mod = drag_model;
            size_t i=0;
            number e=10.0;
            number w2=6.8598479912268;
            number w1=6.8598478663758;
            number re, c;
            while(e>E1)
            {
                
                w1=w2;
                re=RE(mu_a1,rho_a1,dp1,w1);
                c = CD(re,mod);
                w2 = sqrt((4.0/3.0)*dp1*(rho_s1/rho_a1-1.0)*g1/c);
                e = 100.0 * fabs(w2-w1)/w2;
                i=i+1;
                if(i > 40) UG_THROW("Error in RelativeVelocityLinker: Reached " << i <<" iterations in RelVel. Vel = "  << w2 <<"   Error = "<<e<<".");
            }
            iter = i;
            Rel = w2;

            if(std::isnan(Rel) || Rel<0.0) UG_THROW("Error in RelativeVelocityLinker: Relative Vel = "<< Rel <<".");
        }
        number RE(const number mu_a1, const number rho_a1, const number dp1, const number w1)
        {
            return rho_a1*dp1*w1/mu_a1;
        }
        number CD(const number RE, const size_t mod)
        {
            number c;
            if(mod == 0)
                c = 24.0 / RE;
            else if (mod == 1)
                c=pow(0.63+4.8/sqrt(RE),2.0);
            else if (mod == 2)
            {//Schiller-Naumann
                if(RE>1000.0)
                    c = 0.44;
                else
                    c = (24.0/RE)*(1.0+0.15*pow(RE,0.687));
                
            }
            else if (mod == 3)
            {//Turton and Levenspiel
                c=(24.0/RE)*(1.0+0.173*pow(RE,0.657))+0.413/(1.0+16300.0*pow(RE,-1.09));
            }
            else UG_THROW("Error in RelativeVelocityLinker: Drag Coefficient.");
            return c;
        }
    
        /*
         if(Inter->boolConsistenGravity())
         {
             StdLinConsistentGravity<refDim> RhoG;
             MathVector<refDim> vConsGravity[numSh];
             MathVector<dim> Gravity; VecSet(Gravity,0.0); Gravity[dim-1]=9.81;
             RhoG.template prepare<dim>(vConsGravity, numSh, vCornerCoords, m_imDensitySCV.values(), Gravity);
             
             for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
             {
                 //     get current SCVF
                 const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
                 
                 RhoG.template compute<dim>(vConsGravitySCVF[ip],scvf.local_ip(), scvf.JTInv(),scvf.local_grad_vector(), vConsGravity);
             }
             
             for(size_t ip = 0; ip < geo.num_scv(); ++ip)
             {
                 //     get current SCV
                 const typename TFVGeom::SCV& scv = geo.scv(ip);
                 
                 RhoG.template compute<dim>(vConsGravitySCV[ip],scv.local_ip(), scv.JTInv(),scv.local_grad_vector(), vConsGravity);
             }
         }
         */
        template <typename TElem, typename TFVGeom>
        inline
        void ConsistentGravitySCV( MathVector<dim>* vConsGravitySCV, const TFVGeom& geo, const MathVector<dim> vCornerCoords[], size_t nip, const number* DensitySCV)
        {
            static const int refDim = TElem::dim;
            static const size_t numSh = reference_element_traits<TElem>::reference_element_type::numCorners;
            StdLinConsistentGravityX<refDim> RhoG;
            MathVector<refDim> vConsGravity[numSh];
            MathVector<dim> Gravity; VecSet(Gravity,0.0); Gravity[dim-1]=-9.81;
            RhoG.template prepare<dim>(vConsGravity, numSh, vCornerCoords, DensitySCV, Gravity);
            
            for(size_t ip = 0; ip < nip; ++ip)
            {
                //     get current SCV
                const typename TFVGeom::SCV& scv = geo.scv(ip);
                
                RhoG.template compute<dim>(vConsGravitySCV[ip],scv.local_ip(), scv.JTInv(),scv.local_grad_vector(), vConsGravity);
            }

        }
    
        template <typename TElem, typename TFVGeom>
        inline
        void ConsistentGravitySCVF( MathVector<dim>* vConsGravitySCVF, const TFVGeom& geo, const MathVector<dim> vCornerCoords[], size_t nip, const number* DensitySCV)
        {
            static const int refDim = TElem::dim;
            static const size_t numSh = reference_element_traits<TElem>::reference_element_type::numCorners;
            StdLinConsistentGravityX<refDim> RhoG;
            MathVector<refDim> vConsGravity[numSh];
            MathVector<dim> Gravity; VecSet(Gravity,0.0); Gravity[dim-1]=-9.81;
            RhoG.template prepare<dim>(vConsGravity, numSh, vCornerCoords, DensitySCV, Gravity);
            
            for(size_t ip = 0; ip < nip; ++ip)
            {
                //     get current SCVF
                const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
                
                RhoG.template compute<dim>(vConsGravitySCVF[ip],scvf.local_ip(), scvf.JTInv(),scvf.local_grad_vector(), vConsGravity);
            }
            

        }
    
        
        void integration_points( bool &b_scvf, bool &b_scv, GridObject* elem, const MathVector<dim> vCornerCoords[], const MathVector<dim> vGlobIP[], const size_t nip)
        {
            
            if(nip == 0) return;
            b_scvf=true;
            b_scv=true;
            
            DimFV1Geometry<dim> geo;
            geo.update(elem, vCornerCoords, NULL);
            
            const size_t numSCVFip = geo.num_scvf_ips();
            const size_t numSCVip = geo.num_scv_ips();
            
            
            if(geo.num_scvf_ips()!=nip)
            {
                b_scvf=false;
            }
            if(geo.num_scv_ips()!=nip || geo.num_sh() != nip)
            {
                b_scv=false;
            }
            
            for(size_t ip = 0; ip < nip; ++ip)
            {
                if (ip<numSCVFip)
                {
                    b_scvf =(b_scvf&&(vGlobIP[ip] == geo.scvf_global_ips()[ip]));
                }
                if (ip<numSCVip)
                {
                    b_scv  =(b_scv&&(vGlobIP[ip]  == geo.scv_global_ips()[ip]));
                }
            }

            /*if( !b_scv && !b_scvf )
            {
                printf("Other ips ----------------------------------------------\n");
                for(size_t ip = 0; ip < nip; ++ip)
                {
                    printf("x1 vGlobIP  %f   vCornerCoords   %f      geo    %f\n" ,vGlobIP[ip][0], vCornerCoords[ip][0], geo.scvf_global_ips()[ip][0]   );
                    printf("y1 vGlobIP  %f   vCornerCoords   %f      geo    %f\n" ,vGlobIP[ip][1], vCornerCoords[ip][1], geo.scvf_global_ips()[ip][1]   );
                   
                }
            }*/
   
        }
    
        void cut_element(number &value ,bool &cut_elem, bool &phase2, LocalVector* u, const number interface)
        {
            
            bool cut = false;
            size_t _C_=(*u).num_all_fct() -1 ;
            (*u).access_all();

            const size_t numSH=(*u).num_all_dof(_C_);
            size_t inside=0;
            size_t outside=0;
            number c;
            value = 1000;
            for(size_t sh = 0; sh < numSH; ++sh)
            {
                c = fmin(1.0, fmax((*u)(_C_,sh),0));
                value = fmin(value, c);
                if (c>interface)
                    inside += 1;
                else
                    outside +=1;
            }
            
            if (inside==numSH || outside == numSH)
            {
                cut = false;
                phase2 = (inside==numSH)? true : false;
            }
            else
            {
                cut = true;
                phase2 = false;
            }
            
            cut_elem=cut;
        }
        /*void cut_elementSCV( bool &cut_scv , number* N_value, const number interface, size_t nip, LocalVector* u)
        {
            
            bool cut = false;

            size_t inside=0;
            number c;
            number N;
                        
            size_t _C_=(*u).num_all_fct() -1 ;
            (*u).access_all();
            
            for(size_t ip = 0; ip < nip; ++ip)
            {
                c = fmin(1.0, fmax((*u)(_C_,ip),0));
                
                if (c>=interface)
                {
                    inside += 1;
                    
                     N = 1.0;
                }
                else
                {
                    N = -1.0;
                }
                if (N_value != NULL)
                    N_value[ip]=N;
                
            }
            
            if (inside != nip && inside != 0)
                cut = true;

            cut_scv=cut;
        }*/
        
    
    
        void MatAddTraspose( MathMatrix<dim,dim>& mOut, const MathMatrix<dim,dim> m1)
        {
            for(size_t i = 0; i < dim; ++i)
                for(size_t j = 0; j < dim; ++j)
                {
                    mOut[i][j] = m1[i][j] + m1[j][i];
                }
        }
    
        number MatDotTraspose( MathMatrix<dim,dim> m1, const MathVector<dim> m2)
        {
            number s = 0.0;
            for(size_t i = 0; i < dim; ++i)
                for(size_t j = 0; j < dim; ++j)
                    s += m2[i] * m1[i][j] * m2[j];
                
            
            return s;
        }
    
        MathVector<2> MaxMin( number* N_value, const number* density, const size_t nip)
        {
            MathVector<2> DensityJump;
            DensityJump[1] = -100000;
            DensityJump[0] = 100000000;
            for(size_t ip = 0; ip < nip; ++ip)
            {
                if (N_value[ip]>0)  DensityJump[1] = fmax( DensityJump[1], density[ip] );
                else DensityJump[0] = fmin( DensityJump[0], density[ip] );
            }
            return DensityJump;
        }
        MathVector<2> MaxMin( number* N_value, std::vector<number> density, const size_t nip)
        {
            MathVector<2> DensityJump;
            DensityJump[1] = -100000;
            DensityJump[0] = 100000000;
            for(size_t ip = 0; ip < nip; ++ip)
            {
                if (N_value[ip]>0)  DensityJump[1] = fmax( DensityJump[1], density[ip] );
                else DensityJump[0] = fmin( DensityJump[0], density[ip] );
            }
            return DensityJump;
        }
        MathVector<2> MaxMin( number* N_value,const number* viscosity, const number* density, const size_t nip)
        {
            MathVector<2> ViscJump;
            ViscJump[1] = -100000;
            ViscJump[0] = 100000000;
            for(size_t ip = 0; ip < nip; ++ip)
            {
                if (N_value[ip]>0)  ViscJump[1] = fmax( ViscJump[1], density[ip] * viscosity[ip]);
                else ViscJump[0] = fmin( ViscJump[0], density[ip] * viscosity[ip]);
            }
            return ViscJump;
        }
        MathVector<2> MaxMin( number* N_value, std::vector<number> viscosity, std::vector<number> density, const size_t nip)
        {
            MathVector<2> ViscJump;
            ViscJump[1] = -100000;
            ViscJump[0] = 100000000;
            for(size_t ip = 0; ip < nip; ++ip)
            {
                if (N_value[ip]>0)  ViscJump[1] = fmax( ViscJump[1], density[ip] * viscosity[ip]);
                else ViscJump[0] = fmin( ViscJump[0], density[ip] * viscosity[ip]);
            }
            return ViscJump;
        }
        void interface_position( MathVector<dim>* x_interface, const MathVector<dim> vCornerCoords[], const number* N_value, const number rho_l, const number rho_g, const number interface ,DimFV1Geometry<dim> geo,const size_t nip, LocalVector*  u)
        {
                        
            size_t _C_=(*u).num_all_fct() -1 ;
            (*u).access_all();
            
            MathVector<dim> x, DX;
            size_t NumSCVF = geo.num_scvf();
            number n[nip];
            number theta_to, theta_from, c_to, c_from, DC;
            number rho;
            for(size_t ip = 0; ip < nip; ++ip)
            {
                n[ip] = 0.0;
                VecSet(x_interface[ip],0.0);
            }
            for(size_t ip = 0; ip < NumSCVF; ++ip)
            {
                const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
                size_t from=scvf.from();
                size_t to=scvf.to();
                
                if (N_value[from]*N_value[to] < 0.0)
                {

                    c_from = fmin(1.0, fmax((*u)(_C_,from),0));
                    c_to = fmin(1.0, fmax((*u)(_C_,to),0));

                    
                    DC=c_to-c_from;
                    VecSubtract(DX,vCornerCoords[to],vCornerCoords[from]);
                    
                    theta_to=  (c_to   - interface)/DC;
                    theta_from=(c_from - interface)/DC;
                    
                    VecScaleAppend(x_interface[to], theta_to,   DX);
                    VecScaleAppend(x_interface[from],   theta_from, DX);
                    
                    n[from] += 1.0;
                    n[to] += 1.0;
                }
            }
            for(size_t ip = 0; ip < nip; ++ip)
            {

                rho = (N_value[ip] > 0.0)? rho_g:rho_l;
                VecScale(x_interface[ip], x_interface[ip], (rho_l-rho_g) / (n[ip] * rho) );
                
            }
            
        }
        number interface_density(number VolFrac_surf, LocalVector* u, GridObject* elem, const MathVector<dim> vCornerCoords[], const size_t ip, const number rho_s, const number rho_a, const bool m_scvf, const bool m_scv)
        {
            
            size_t _C_=(*u).num_all_fct() -1 ;
            size_t numSH=(*u).num_all_dof(_C_);
            if (m_scvf)
            {
                number c1, c2, rho1, rho2;
                
                DimFV1Geometry<dim> geo;
                geo.update(elem, vCornerCoords, NULL);
                
                
                const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
                
                c1 = fmin(1.0, fmax((*u)(_C_,scvf.from()),0));
                c2 = fmin(1.0, fmax((*u)(_C_,scvf.to()),0));
                rho1=(rho_s-rho_a)*c1+rho_a;
                rho2=(rho_s-rho_a)*c2+rho_a;
                
                return fmin(rho1,rho2);
            }
            else
            {
                number rho=1e10, c;
                number rho_aux;
                for(size_t sh = 0; sh < numSH; ++sh)
                {
                    c = fmin(1.0, fmax((*u)(_C_,sh),0));
                    rho_aux=(rho_s-rho_a)*c+rho_a;
                    rho = fmin (rho_aux,rho);
                    
                }
                return rho;
                
                
            }


        }
    
        template <typename TElem, typename TFVGeom>
        inline
        void PropertiesJump(const TFVGeom& geo, const DataImport<number, dim>& VolFraction, const DataImport<number, dim>& JumpShape, const DataImport<number, dim>& DensitySCV, const DataImport<number, dim>& KinViscSCV, DataImport<MathVector<dim>, dim>& SourceSCV,size_t numSh, bool& interface, bool* Phase2, number& mu_l, number& mu_g, number& rho_l, number& rho_g, MathVector<dim>& Source_l, MathVector<dim>& Source_g, const number m_interface_vol_fraction)
        {
            
            UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");
            
            interface = cut_interface(JumpShape, numSh);
            
            
            
            if(!interface)
            {
                return;
            }
            
            /*if (true)
            {
                
                for(size_t sh = 0; sh < geo.num_scvf(); ++sh)
                    printf("Vol[%zu] = %f\n", sh, VolFraction[sh]);
                
                printf("Pos[0] = %f     %f\n", geo.scv_global_ips()[0][0],geo.scv_global_ips()[0][1]);
                printf("Pos[1] = %f     %f\n", geo.scv_global_ips()[1][0],geo.scv_global_ips()[1][1]);
                printf("Pos[2] = %f     %f\n", geo.scv_global_ips()[2][0],geo.scv_global_ips()[2][1]);
                

            }*/
            
            RhoMuSource(mu_l,  mu_g,  rho_l, rho_g, Source_l, Source_g, DensitySCV, KinViscSCV, SourceSCV, JumpShape, numSh);

            
            /*number Cval=0.0;
            number vol = 0.0;
            //     loop Sub Control Volumes (SCV)
            for(size_t ip = 0; ip < geo.num_scv(); ++ip)
            {
                //     get current SCV
                const typename TFVGeom::SCV& scv = geo.scv(ip);
                
                //     get associated node
                const int sh = scv.node_id();
                
                Cval += VolFraction[sh]*scv.volume();
                vol += scv.volume();

            }
            Cval *= 1.0/vol;*/
            this->template phase<TFVGeom>(geo,JumpShape.values(), Phase2);
            
        }
        bool cut_interface(const DataImport<number, dim>& JumpShape, const size_t numSh)
        {
            
            number max=0;
            number min=0;
            
            bool interface=false;
            for(size_t sh = 0; sh < numSh; ++sh)
            {
                max=fmax(JumpShape[sh],max);
                min=fmin(JumpShape[sh],min);
                
            }
            if (-0.5<(max+min) && (max+min)<0.5)
            {
                interface=true;
                
            }
            return interface;


        }
    
        template <typename TElem, typename TFVGeom>
        inline
        void phase(const TFVGeom& geo, const number JumpShape[], bool* Phase2)
        {
            
            for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
            {
                //     get current SCV
                const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
                if(JumpShape[scvf.to()]*JumpShape[scvf.from()]<0.0)
                {

                    Phase2[ip]=true;
                    
                }
                else
                {

                    if(JumpShape[scvf.to()]>0.0)
                    {
                        Phase2[ip]=true;
                    }
                    else
                    {
                        Phase2[ip]=false;
                    }
                }
                
            }


        }
        template <typename TElem, typename TFVGeom>
        inline
        bool cut_interface_concentration( const TFVGeom& geo, const LocalVector& u, const size_t _C_, const size_t numSh, number interface_value, size_t& num_inside, size_t& num_outside, number& value_in, number& value_out, number* JumpShape)
        {
            bool interface = false;

            size_t inside=0;
            size_t outside=0;
            value_in = 0.0;
            value_out = 0.0;
            
            for(size_t sh = 0; sh < numSh; ++sh)
            {

                if (u(_C_, sh)>interface_value)
                {
                    inside += 1;
                    JumpShape[sh] = 1;
                    value_in += u(_C_, sh);
                }
                else
                {
                    outside +=1;
                    JumpShape[sh] = -1;
                    value_out += u(_C_, sh);
                }
            }
            
            if (inside==numSh || outside == numSh)
            {
                interface = false;
                value_in = (value_in + value_out) / (inside+outside);
                value_out = value_in;
                inside = inside+outside;
                outside = inside;
            }
            else
            {
                interface = true;
                value_in = value_in / inside;
                value_out = value_out / outside;
            }
            
            num_inside = inside;
            num_outside = outside;
            
            
            return interface;


        }
    
        void RhoMuSource( number& mu_l, number& mu_g, number& rho_l, number& rho_g, MathVector<dim>& vSource_l, MathVector<dim>& vSource_g, const DataImport<number, dim>& DensitySCV, const DataImport<number, dim>& KinViscSCV, const DataImport<MathVector<dim>, dim>& SourceSCV, const DataImport<number, dim>& JumpShape, const size_t numSh)
        {
    
            bool boolSource = (SourceSCV.data_given()) ? true : false;
            
            number mu_2 = 0.0, mu_1 = 0.0;
            number rho_2 = 0.0, rho_1 = 0.0;
            

            MathVector<dim> Source_1;
            MathVector<dim> Source_2;
            
            VecSet(Source_1,0.0);
            VecSet(Source_2,0.0);
            
            
            int Count_1 = 0;
            int Count_2 = 0;
            
            for(size_t sh = 0; sh < numSh; ++sh)
            {
                if (JumpShape[sh]>0)
                {
                    mu_2  +=  (DensitySCV[sh] * KinViscSCV[sh]);
                    rho_2 +=  (DensitySCV[sh]) ;

                    if(boolSource) {VecScaleAppend(Source_2 , 1.0  ,SourceSCV[sh] );}
                    
                    Count_2 +=1;
                    
                }
                else
                {
                    mu_1  +=  (DensitySCV[sh] * KinViscSCV[sh]);
                    rho_1 +=  DensitySCV[sh] ;
                    if(boolSource) {VecScaleAppend(Source_1   , 1.0  ,SourceSCV[sh] );}
                    Count_1 +=1;
                    
                }
                
            }
            mu_1 =  mu_1 / Count_1 ;
            mu_2 =  mu_2 / Count_2 ;
            rho_1 = rho_1 / Count_1;
            rho_2 = rho_2 / Count_2;
            
            if(boolSource)
            {
                VecScale(Source_1,Source_1,1.0/Count_1);
                VecScale(Source_2,Source_2,1.0/Count_2);
            }
            
            if ((mu_2 < mu_1)||(mu_2<=0.0)||(mu_1<=0.0))
            {
                printf("Mu2 = %f    Mu1 = %f\n",mu_2,mu_1);
                printf("Rho2 = %f    Rho1 = %f\n",rho_2,rho_1);
                
                for(size_t sh = 0; sh < numSh; ++sh)
                    printf("mu[%zu] = %f\n", sh, DensitySCV[sh] * KinViscSCV[sh]);
                for(size_t sh = 0; sh < numSh; ++sh)
                    printf("rho[%zu] = %f\n", sh, DensitySCV[sh] );
                for(size_t sh = 0; sh < numSh; ++sh)
                    printf("JumpShape[%zu] = %f\n", sh, JumpShape[sh] );
                
                
                
                UG_THROW("Viscosity in phase 1 is lower that phase 2");
            }
            if ((rho_2 < rho_1)||(rho_2<=0.0)||(rho_1<=0.0))
                UG_THROW("Density in phase 1 is lower that phase 2");
            
            
            
            
            mu_l = mu_2;
            mu_g = mu_1;
            rho_l = rho_2;
            rho_g = rho_1;

            vSource_l = Source_2;
            vSource_g = Source_1;

        }
    
        template <typename TElem, typename TFVGeom>
        inline
        void InterfaceSCVFShapes( number** SCVFinterShape, const TFVGeom& geo, const DataImport<number, dim>& VolFraction, const DataImport<number, dim>& JumpShape, size_t numSh, const number m_interface_vol_fraction)
        {
            /////////////////////////////////////////////////////////////////////////////
            // Calculation X_interface
            /////////////////////////////////////////////////////////////////////////////
            
            /*number** SCVFinterShape = new number*[numSCVF];
            for(size_t count=0; count<numSCVF; count++)
                SCVFinterShape[count] = new number[numSh];*/
            
            UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");
            
            size_t numSCVF = geo.num_scvf();
            typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
            static const int refDim = ref_elem_type::dim;
            
            MathVector<refDim> vLocIP_inter;
            VecSet(vLocIP_inter,0.0);
            
            MathVector<refDim> vLocIP_SCVF[numSCVF];
            
            number theta_to, c_to, c_from, DC;
            

            number count_interface=0;

            for(size_t ip = 0; ip < numSCVF; ++ip)
            {
                const typename FV1Geometry<TElem, dim>::SCVF scvf = geo.scvf(ip);

                const size_t from=scvf.from();
                const size_t to=scvf.to();
                
                VecSet(vLocIP_SCVF[ip],0.0);
                
                if (JumpShape[from]*JumpShape[to] < 0.0)
                {
                    c_from = VolFraction[from];
                    c_to = VolFraction[to];

                    DC=c_to-c_from;
                    
                    theta_to=  (c_to   - m_interface_vol_fraction)/DC;
                    
                    count_interface += 1.0;
                    VecScaleAppend(vLocIP_SCVF[ip], 1.0  ,geo.scv_local_ips()[to],-1.0 * theta_to, geo.scv_local_ips()[to],theta_to,geo.scv_local_ips()[from]);
                    VecScaleAppend(vLocIP_inter   , 1.0  ,vLocIP_SCVF[ip]);
                    
                }
                else
                {
                    VecScaleAppend(vLocIP_SCVF[ip], 1.0  ,geo.scv_local_ips()[to],-1.0 * 0.5, geo.scv_local_ips()[to],0.5,geo.scv_local_ips()[from]);
                }

            }
            

            VecScale(vLocIP_inter,vLocIP_inter,1.0/count_interface);
            
            LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();


        //    Reference Mapping

            ReferenceMapping<ref_elem_type, dim> mapping(geo.scv_global_ips());
            
            
            

            
            for(size_t ip = 0; ip < numSCVF; ++ip)
            {
                
                VecScale(vLocIP_SCVF[ip],vLocIP_SCVF[ip],0.5);
                VecScaleAppend(vLocIP_SCVF[ip], 0.5  ,vLocIP_inter);
                               
                rTrialSpace.shapes(SCVFinterShape[ip], vLocIP_SCVF[ip]);
                
                /*const typename FV1Geometry<TElem, dim>::SCVF scvf = geo.scvf(ip);
                const size_t from=scvf.from();
                const size_t to=scvf.to();
                printf("vGlobIP_SCVF[from][0] =  %f \n", geo.scv_global_ips()[from][0]);
                printf("vGlobIP_SCVF[from][1] =  %f \n", geo.scv_global_ips()[from][1]);
                printf("vGlobIP_SCVF[to][0] =  %f \n",geo.scv_global_ips()[to][0]);
                printf("vGlobIP_SCVF[to][1] =  %f \n",geo.scv_global_ips()[to][1]);
                printf(" Shapes  ------------------------------------");
                for(size_t sh = 0; sh < numSh; ++sh)
                    printf("SCVF shape [%zu][%zu] = %f\n",ip,sh,interShape[ip][sh]);*/

                
            }
        }
    
        template <typename TElem, typename TFVGeom>
        inline
        void InterfaceShape(number* InterfaceShape, const TFVGeom& geo, const DataImport<number, dim>& VolFraction, const DataImport<number, dim>& JumpShape, size_t numSh, const number m_interface_vol_fraction)
        {
            /////////////////////////////////////////////////////////////////////////////
            // Calculation X_interface
            /////////////////////////////////////////////////////////////////////////////
            
            /*number** SCVFinterShape = new number*[numSCVF];
            for(size_t count=0; count<numSCVF; count++)
                SCVFinterShape[count] = new number[numSh];*/
            
            UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");
            
            size_t numSCVF = geo.num_scvf();
            typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
            static const int refDim = ref_elem_type::dim;
            
            MathVector<refDim> vLocIP_inter;
            VecSet(vLocIP_inter,0.0);
    
            
            number theta_to, c_to, c_from, DC;
            

            number count_interface=0;

            for(size_t ip = 0; ip < numSCVF; ++ip)
            {
                const typename FV1Geometry<TElem, dim>::SCVF scvf = geo.scvf(ip);

                const size_t from=scvf.from();
                const size_t to=scvf.to();
                
                
                if (JumpShape[from]*JumpShape[to] < 0.0)
                {
                    c_from = VolFraction[from];
                    c_to = VolFraction[to];

                    DC=c_to-c_from;
                    
                    theta_to=  (c_to   - m_interface_vol_fraction)/DC;
                    
                    count_interface += 1.0;
                    VecScaleAppend(vLocIP_inter, 1.0  ,geo.scv_local_ips()[to],-1.0 * theta_to, geo.scv_local_ips()[to],theta_to,geo.scv_local_ips()[from]);
                    
                }

            }
            

            VecScale(vLocIP_inter,vLocIP_inter,1.0/count_interface);
            
            LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();


        //    Reference Mapping

            ReferenceMapping<ref_elem_type, dim> mapping(geo.scv_global_ips());
            
            rTrialSpace.shapes(InterfaceShape, vLocIP_inter);
            

        }
    public:
        void set_particle_density(float R) {
            rho_s = R;
        }
        void set_air_density(float R) {
            rho_a = R;
        }
        void set_fluid_Visc(float R) {
            mu_a = R;
        }
        void set_particle_diameter(float R) {
            dp = R;
        }
        void set_alpha_max(float R) {
            alpha_max = R;
        }
        void set_alpha_min(float R) {
            alpha_min = R;
        }
        /*void set_packing_factor(float R) {
            m_packing_factor = R;
        }*/
        void set_FR(float R) {
            Fr = R;
        }
        void set_B_phi(float R) {
            B_phi = R;
        }
        void set_deltaGamma(float R) {
            deltaGamma = R;
        }
        void set_limit(float R) {
            m_limit = R;
        }
        void set_time_step_factor(float R) {
            m_dt = R;
        }
        void set_reference_pressure(float R) {
            m_P0 = R;
        }
        void set_bool_particle_pressure_force(bool R) {
            m_bParticleGradientForce = R;
        }
        void set_bool_consistent_gravity(bool R) {
            m_bConsistentGravity = R;
        }
        void set_bool_initialized(bool R) {
            m_init = R;
        }
        void set_drag_model(size_t R) {
            drag_model = R;
        }

        number Density_s(){ return rho_s;}
        number Density_a(){ return rho_a;}
        number Density_max(){ return (rho_s-rho_a) * alpha_max + rho_a;}
		number Viscosity_a(){ return mu_a;}
        number Alpha_max(){ return alpha_max;}
        number DT(){ return m_dt;}
        number ReferencePressure(){ return m_P0;}
        bool ParticleGradientForce(){ return m_bParticleGradientForce;}
        bool boolConsistentGravity(){ return m_bConsistentGravity;}
        
        size_t DragModel(){ return drag_model;}
        bool valid(){ return m_init;}
    
    
    protected:
        float m_P0, rho_s, rho_a, mu_a, dp, Fr, B_phi, alpha_max, alpha_min, deltaGamma, m_limit, m_dt;
        size_t drag_model;
        bool m_bParticleGradientForce, m_bConsistentGravity, m_init;
        

        


};



}

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC_NV_CUT_ELEMENT_ */


