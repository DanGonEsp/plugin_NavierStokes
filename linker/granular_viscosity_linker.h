
/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Jonas Simon
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_VISCOSITY_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_VISCOSITY_LINKER__

#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "../properties_interface.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// Granular Viscosity linker
////////////////////////////////////////////////////////////////////////////////

/// Linker for the Granular viscosity
/**
 * This linker computes the Bingham viscosity \f$ \mathbf{q} = - \frac{\mathbf{K}}{\mu} ( \nabla p - \rho \mathbf{g} ) \f$,
 * where
 * <ul>
 * <li> \f$ \mathbf{K} \f$        permeability
 * <li> \f$ \mu \f$                viscosity
 * <li> \f$ \rho \f$            density
 * <li> \f$ \mathbf{g} \f$        gravity
 * <li> \f$ \nabla p \f$        pressure gradient
 * </ul>
 * are input parameters.
 */

template <int dim>
class GranularViscosityLinker
    : public StdDataLinker< GranularViscosityLinker<dim>, number, dim>
{
    //    Base class type
    typedef StdDataLinker< GranularViscosityLinker<dim>, number, dim> base_type;

    //  Constructor
    public:
        GranularViscosityLinker() :
            m_spVolumeFraction(NULL), m_spDVolumeFraction(NULL),
            m_spVelocityGrad(NULL), m_spDVelocityGrad(NULL),
            m_spMixDensity(NULL), m_spDMixDensity(NULL),
            m_spParticlePressure(NULL), m_spDParticlePressure(NULL),
            m_spMixViscosity(NULL), m_spDMixViscosity(NULL),
            m_model(0),
            //m_packing_factor(1),
            rho_s(1600.0),
            nu_s(1.0e-06),
            rho_a(1.2),
            
            mu_a(1.48e-5),
            I_0(0.279),
            FricMu_1(0.38),
            FricMu_2(0.64),
            dp(0.001),

            alpha_max(0.635),
            alpha_min(0.57),
            deltaI(1e-3),
            deltaPs(1.48e-04),
            deltaGamma(1e-4),
            interface_volume_fraction(0.5)
        {
        //    this linker needs exactly four input
            this->set_num_input(5);
        }

        // function for evaluation at single ip?
        inline void evaluate (number& value,
                            const MathVector<dim>& globIP,
                            number time, int si) const
        {
            UG_LOG("GranularViscosityLinker::evaluate single called");
            number volume_fraction;
            MathMatrix<dim,dim> velocityGrad;

            (*m_spVolumeFraction)(volume_fraction, globIP, time, si);
            (*m_spVelocityGrad)(velocityGrad, globIP, time, si);

            // Compute second invariant
            number gamma = 0.0;
            number I = 0.0;
            number mu_friction = 0.0;
            number viscosity = 0.0;
            

            if (true){
                printf("Flag 1");
            }
            // compute inner sum
            for(int d1 = 0; d1 < dim; ++d1)
            {
                for(int d2 = 0; d2 < dim; ++d2)
                {
                    //gamma += pow(velocityGrad(d1,d2) + velocityGrad(d2,d1),2);
                }
            }
            gamma = 0.0001+sqrt((0.5*gamma));
            I=gamma*dp/sqrt((0.0001+abs(1))/rho_s);
            mu_friction=FricMu_1+(FricMu_2-FricMu_1)/(1.0+I_0/I);
            viscosity = mu_friction*1/gamma;
            
            
            //value = viscosity + yieldStress/;
            // compute mu = eta + sigma / \sqrt( delta + 1 / I)
            value = volume_fraction*viscosity + (1.0-volume_fraction)*mu_a;
        }

        // function for evaluation at multiple ips??
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
            

            std::vector<number> vVolumeFraction(nip);
            std::vector<MathMatrix<dim,dim> > vVelocityGrad(nip);
            std::vector<number> vMixDensity(nip);
            std::vector<number> vParticlePressure(nip);
            std::vector<number> vMixViscosity(nip);

            (*m_spVolumeFraction)(&vVolumeFraction[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spVelocityGrad)(&vVelocityGrad[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spMixDensity)(&vMixDensity[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spParticlePressure)(&vParticlePressure[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spMixViscosity)(&vMixViscosity[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);

            
            number VolFraction;
            number viscosity_granular[nip];
            number viscosity_granular_aux;
            
            number gamma[nip];
            number gamma_aux;
            
            //number mu_eins[nip];
            //number mu_eins_aux;
            
            number mu_sand[nip];
            number mu_sand_aux;
            
            for(size_t ip = 0; ip < nip; ++ip)
            {
                mu_sand_aux=0;
                //mu_eins_aux=mu_a;
                
                VolFraction=fmin(1.0, fmax(vVolumeFraction[ip],0.0));

                switch (m_model) {
                    case 0:
                        viscosity_granular_aux = Constant_viscosity(VolFraction, mu_a, nu_s*rho_s,interface_volume_fraction);
                        
                        break;
                    case 1:
                        viscosity_granular_aux = Proportional_viscosity(VolFraction, mu_a, nu_s*rho_s);
                        
                        break;
                    case 2:
                        viscosity_granular_aux= vMixViscosity[ip];
                        
                        break;
                    case 3:
                        Granular_viscosity_1(mu_sand_aux,vVelocityGrad[ip], vParticlePressure[ip], rho_s, dp, FricMu_1, FricMu_2, I_0, deltaGamma, deltaPs, deltaI,gamma_aux);

                        viscosity_granular_aux = vMixViscosity[ip] + mu_sand_aux ;

                        break;
                    case 4:
                        
                        if (VolFraction>interface_volume_fraction)
                            Granular_viscosity_1(mu_sand_aux,vVelocityGrad[ip], vParticlePressure[ip], rho_s, dp, FricMu_1, FricMu_2, I_0, deltaGamma, deltaPs, deltaI,gamma_aux);
                        
                        viscosity_granular_aux = vMixViscosity[ip] + mu_sand_aux;

                            
                        gamma[ip]=gamma_aux;
                        mu_sand[ip]=mu_sand_aux;
                        //mu_eins[ip]=mu_eins_aux;
                        

                        break;
                    default:
                        UG_THROW("Wrong model selected for granular viscosity: GranularViscosityLinker has options"
                                       " model= 0, 1 , 2, 3")
                        break;
                }
				if(std::isnan(viscosity_granular_aux) || viscosity_granular_aux<0.0) UG_THROW("Error ViscosityLinker: Value = NaN" <<".");
                viscosity_granular[ip]=viscosity_granular_aux;
                vValue[ip]=viscosity_granular[ip]/vMixDensity[ip];



            }

        }

        template <int refDim>
        void eval_and_deriv(number vValue[],
                             const MathVector<dim> vGlobIP[],
                             number time, int si,
                             GridObject* elem,
                             const MathVector<dim> vCornerCoords[],
                             const MathVector<refDim> vLocIP[],
                             const size_t nip,
                             LocalVector* u,
                             bool bDeriv,
                             int s,
                             std::vector<std::vector<number> > vvvDeriv[],
                             const MathMatrix<refDim, dim>* vJT = NULL) const
        {
            
            int s_VOL_ = base_type::series_id(_VOL_, s);
            int s_DV_  = base_type::series_id(_DV_ , s);
            int s_RHO_ = base_type::series_id(_RHO_, s);
            int s_P_   = base_type::series_id(_P_  , s);
            int s_MU_ = base_type::series_id(_MU_, s);
            
            const number*               vVolumeFraction   = m_spVolumeFraction->values(s_VOL_);
            const MathMatrix<dim,dim>*  vVelocityGrad     = m_spVelocityGrad->values(s_DV_);
            const number*               vMixDensity       = m_spMixDensity->values(s_RHO_);
            const number*               vParticlePressure = m_spParticlePressure->values(s_P_);
            const number*               vMixViscosity     = m_spMixViscosity->values(s_MU_);
            
            number VolFraction;
            //number MixViscosity;
            
            number viscosity_granular[nip];
            number viscosity_granular_aux;
            
            number gamma[nip];
            number gamma_aux;
            
            //number mu_eins[nip];
            //number mu_eins_aux;
            
            number mu_sand[nip];
            number mu_sand_aux;
            
            for(size_t ip = 0; ip < nip; ++ip)
            {
                mu_sand_aux=0;
                //mu_eins_aux=mu_a;

                VolFraction=fmin(1.0, fmax(vVolumeFraction[ip],0.0));


                switch (m_model) {
                    case 0:

                        viscosity_granular_aux = Constant_viscosity(VolFraction, mu_a, nu_s*rho_s,interface_volume_fraction);
                        
                        break;
                    case 1:

                        viscosity_granular_aux = Proportional_viscosity(VolFraction, mu_a, nu_s*rho_s);

                        
                        break;
                    case 2:
                        
                        viscosity_granular_aux = vMixViscosity[ip];
                        break;
                    case 3:

                        Granular_viscosity_1(mu_sand_aux,vVelocityGrad[ip], vParticlePressure[ip], rho_s, dp, FricMu_1, FricMu_2, I_0, deltaGamma, deltaPs, deltaI,gamma_aux);

                        viscosity_granular_aux = vMixViscosity[ip]+mu_sand_aux ;
                        
                        gamma[ip]=gamma_aux;
                        mu_sand[ip]=mu_sand_aux;
                        //mu_eins[ip]=mu_eins_aux;

                        break;
                    case 4:
                        
                        if (VolFraction>interface_volume_fraction)
                            Granular_viscosity_1(mu_sand_aux,vVelocityGrad[ip], vParticlePressure[ip], rho_s, dp, FricMu_1, FricMu_2, I_0, deltaGamma, deltaPs, deltaI,gamma_aux);
                        
                        viscosity_granular_aux = vMixViscosity[ip] + mu_sand_aux;

                            
                        gamma[ip]=gamma_aux;
                        mu_sand[ip]=mu_sand_aux;
                        //mu_eins[ip]=mu_eins_aux;

                        break;
                    default:
                        UG_THROW("Wrong model selected for granular viscosity: GranularViscosityLinker has options"
                                       " model= 0, 1 , 2, 3, 4")
                        break;
                }
				if(std::isnan(viscosity_granular_aux) || viscosity_granular_aux<0.0) UG_THROW("Error ViscosityLinker: Value = NaN" <<".");
                viscosity_granular[ip]=viscosity_granular_aux;
                vValue[ip]=viscosity_granular[ip]/vMixDensity[ip];
            

            }
            
            /*bool m_scvf=false;
            bool m_scv=false;
            bool cut_element=false;
            number VolFrac_surface;
            //if (!(vGlobIP!=NULL && vLocIP!=NULL && vCornerCoords!=NULL)) UG_THROW("Error in Granular Source");
            
            Inter->integration_points(m_scvf, m_scv, elem,   vCornerCoords,   vGlobIP, nip);
            Inter->cut_element(VolFrac_surface, cut_element, u, 0.9047);
            bool f = cut_element && m_scv;
            if(cut_element && m_scv)
                for(size_t ip = 0; ip < nip; ++ip)
                {
                
                    if (!((vGlobIP[ip][0] > 24.124) && (vGlobIP[ip][0] < 24.376)))
                    {
                        f = f && false;
                        
                    }
                }
            if(f)
            {
                for(size_t ip2 = 0; ip2 < nip; ++ip2)
                {
                    printf("Visc_Solid[%zu] = %f        %f      %f\n",ip2, vValue[ip2]*vMixDensity[ip2], vVolumeFraction[ip2], vParticlePressure[ip2]);
                }
                for(size_t ip2 = 0; ip2 < nip; ++ip2)
                {
                    printf("MUMU_Solid[%zu] = %f        %f\n",ip2, mu_sand[ip2]*vMixDensity[ip2], mu_eins[ip2]);
                }
                for(size_t ip2 = 0; ip2 < nip; ++ip2)
                {
                    printf("U_solid[%zu] = %lf      %lf\n",ip2, vCornerCoords[ip2][0],vCornerCoords[ip2][1]);
                }
                for(size_t ip2 = 0; ip2 < nip; ++ip2)
                {
                    printf("vGlobIP[%zu] = %lf      %lf\n",ip2, vGlobIP[ip2][0],vGlobIP[ip2][1]);
                }
                if (m_scv) printf("SCV_____\n");
                if (m_scvf) printf("SCVF_____\n");
            }*/

            //    Compute the derivatives at all ips     //
            /////////////////////////////////////////////

        //    check if something is left to do
            if(!bDeriv || this->zero_derivative()) return;

        //    clear all derivative values
            this->set_zero(vvvDeriv, nip);
            
        //  Derivatives of Volume Fraction
            /*if(m_spDVolumeFraction.valid() && !m_spDVolumeFraction->zero_derivative() && (m_model == 0 || m_model == 1 || m_model == 4))
            {

                for(size_t ip = 0; ip < nip; ++ip)
                {
                    VolFraction=fmin(1.0, fmax(vVolumeFraction[ip],0.0));
                    number deriv1=0.0;
                    number dmu1=0.0;
                    number deriv2=0.0;
                    number dmu2=0.0;
                    number deriv=0.0;
                    number mu_dev=0.0;
                    switch (m_model) {
                        case 1:
                            mu_dev= (nu_s*rho_s-mu_a);
                            break;
                        case 2:
                            mu_dev= (2.5*mu_eins[ip]/(1.0-VolFraction/alpha_max));
                            break;
                        case 3:
   
                            dmu2=2.5*(mu_eins[ip])/(1.0-VolFraction/alpha_max);
                            deriv2 = 1.0;
                            //deriv2 = pow(  viscosity_granular[ip] , 2) * ( VolFraction / pow(mu_eins[ip]+mu_sand[ip],2) + (1-VolFraction)/pow(mu_eins[ip],2)    );
                            //deriv = -pow(  viscosity_granular[ip] , 2) * ( 1 / (mu_eins[ip]+mu_sand[ip]) - 1 / mu_eins[ip]    );
                            
                            
                            mu_dev= deriv1 * dmu1  +  deriv2 * dmu2  +  deriv;

                            break;
                        case 4:
                            
                            dmu2=2.5*(mu_eins[ip])/(1-vVolumeFraction[ip]/alpha_max);
                            deriv2 = 1.0;
                            
                            deriv = 0.0;
                            mu_dev= deriv1 * dmu1  +  deriv2 * dmu2  + deriv;

                            break;
                        default:
                            UG_THROW("Error in GranularViscosityLinker model= 0, 1,2")
                            break;
                    }
                    
                    mu_dev = mu_dev/vMixDensity[ip];
                
                    
                    for(size_t fct = 0; fct < m_spDVolumeFraction->num_fct(); ++fct)
                    {
                    //    get derivative of volume fraction w.r.t. to all functions
                        const number* vDVolumeFraction = m_spDVolumeFraction->deriv(s_VOL_, ip, fct);

                    //    get common fct id for this function
                        const size_t commonFct = this->input_common_fct(_VOL_, fct);

                    //    loop all shapes and set the derivative
                        for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                        {
                            vvvDeriv[ip][commonFct][sh] += vDVolumeFraction[sh] * mu_dev;

                        }
                    }
                }
            }*/
            
        //  Derivatives of Volume Fraction
            /*if(m_spDMixViscosity.valid() && !m_spDMixViscosity->zero_derivative() && (m_model == 2 || m_model == 3 || m_model == 4))
            {

                for(size_t ip = 0; ip < nip; ++ip)
                {
                    VolFraction=fmin(1.0, fmax(vVolumeFraction[ip],0.0));
                    //number deriv1=0.0;
                    //number dmu1=0.0;
                    //number deriv2=0.0;
                    //number dmu2=0.0;
                    //number deriv=0.0;
                    number mu_dev=0.0;
                    switch (m_model) {
                        case 2:
                            mu_dev= 1.0;
                            break;
                        case 3:
                                                        
                            mu_dev= 1.0;

                            break;
                        case 4:
                            
                            //dmu2=1.0;
                            //deriv2 = 1.0;
                            
                            //dmu1=1.0;
                            //deriv1 = 1.0;
                            
                            //deriv = 0.0;
                            mu_dev= 0.0;
                            UG_THROW("Error in GranularViscosityLinker model= 0, 1,2")

                            break;
                        default:
                            UG_THROW("Error in GranularViscosityLinker model= 0, 1,2")
                            break;
                    }
                    
                    mu_dev = mu_dev/vMixDensity[ip];
                
                    
                    for(size_t fct = 0; fct < m_spDMixViscosity->num_fct(); ++fct)
                    {
                    //    get derivative of volume fraction w.r.t. to all functions
                        const number* vDMixViscosity = m_spDMixViscosity->deriv(s_MU_, ip, fct);

                    //    get common fct id for this function
                        const size_t commonFct = this->input_common_fct(_MU_, fct);

                    //    loop all shapes and set the derivative
                        for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                        {
                            vvvDeriv[ip][commonFct][sh] += vDMixViscosity[sh] * mu_dev;

                        }
                    }
                }
            }*/
            
        //  Derivatives of Density
            /*if(m_spDMixDensity.valid() && !m_spDMixDensity->zero_derivative() )
            {
                for(size_t ip = 0; ip < nip; ++ip)
                {
                    for(size_t fct = 0; fct < m_spDMixDensity->num_fct(); ++fct)
                    {
                    //    get derivative of volume fraction w.r.t. to all functions
                        const number* vDMixDensity = m_spDMixDensity->deriv(s_RHO_, ip, fct);

                    //    get common fct id for this function
                        const size_t commonFct = this->input_common_fct(_RHO_, fct);

                    //    loop all shapes and set the derivative
                        for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                        {
                            vvvDeriv[ip][commonFct][sh] += - vDMixDensity[sh]*vValue[ip]/vMixDensity[ip];

                        }
                    }
                }
            }*/
		//  Derivatives of MixViscosity
			
			/*if(m_spDMixViscosity.valid() && !m_spDMixViscosity->zero_derivative() && (m_model== 2  || m_model== 3 || m_model== 4) )
			{
				for(size_t ip = 0; ip < nip; ++ip)
				{
					
					for(size_t fct = 0; fct < m_spDMixViscosity->num_fct(); ++fct)
					{
					//    get derivative of volume fraction w.r.t. to all functions
						const number* vDMixViscosity = m_spDMixViscosity->deriv(s_MU_, ip, fct);

					//    get common fct id for this function
						const size_t commonFct = this->input_common_fct(_MU_, fct);

					//    loop all shapes and set the derivative
						for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
						{
							vvvDeriv[ip][commonFct][sh] += vDMixViscosity[sh] / vMixDensity[ip];
						}
					}
				}
			}
            
        //  Derivatives of ParticlePressure
            
            if(m_spDParticlePressure.valid() && !m_spDParticlePressure->zero_derivative() && (m_model== 3  || m_model== 4) )
            {
                for(size_t ip = 0; ip < nip; ++ip)
                {
                    VolFraction=fmin(1.0, fmax(vVolumeFraction[ip],0));
                    
                    number deriv1=0;
                    number dmu1=0;
                    number deriv2=0;
                    number dmu2=0;
                    number mu_dev=0;

                    switch (m_model) {
                        case 3:
                            DevPressure_1(dmu1, gamma[ip], vParticlePressure[ip], mu_a, rho_s, dp, FricMu_1, FricMu_2, I_0, deltaGamma,deltaPs, deltaI);
                            deriv1 = 1.0;//VolFraction * pow(  mu_eins[ip] , 2) / pow(  mu_eins[ip] + (1-VolFraction)*mu_sand[ip]  ,2);
                            mu_dev=   (deriv1 * dmu1  +  deriv2 * dmu2 );
                            
                            mu_dev = mu_dev / vMixDensity[ip];
                            
                            break;
                        case 4:
                            if ( VolFraction>interface_volume_fraction)
                                DevPressure_1(dmu1, gamma[ip], vParticlePressure[ip], mu_a, rho_s, dp, FricMu_1, FricMu_2, I_0, deltaGamma,deltaPs, deltaI);
                            deriv1 = 1.0;
                            mu_dev=   (deriv1 * dmu1  +  deriv2 * dmu2 );
                            
                            mu_dev = mu_dev / vMixDensity[ip];
                            break;
                        default:
                            UG_THROW("Error in GranularViscosityLinker model= 2, 3")
                            break;
                    }
                    
                    
                    
                    for(size_t fct = 0; fct < m_spDParticlePressure->num_fct(); ++fct)
                    {
                    //    get derivative of volume fraction w.r.t. to all functions
                        const number* vDParticlePressure = m_spDParticlePressure->deriv(s_P_, ip, fct);

                    //    get common fct id for this function
                        const size_t commonFct = this->input_common_fct(_P_, fct);

                    //    loop all shapes and set the derivative
                        for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                        {
                            vvvDeriv[ip][commonFct][sh] += vDParticlePressure[sh] * mu_dev  ;
                        }
                    }
                }
            }*/
        
        //  Derivatives of velocity gradient
            
            /*if(m_spDVelocityGrad.valid() && !m_spDVelocityGrad->zero_derivative() && (m_model==3 || m_model==4))
            {
                for(size_t ip = 0; ip < nip; ++ip)
                {
                    VolFraction=fmin(1.0, fmax(vVolumeFraction[ip],0));
                    
                    MathMatrix<dim,dim> Deriv;
                    MathMatrix<dim,dim> Dmu1;
                    
                    MathMatrix<dim,dim> Mu_deriv;
                    
                    MatSet(Dmu1,0);
                    MatSet(Dmu1,0);
                    number deriv1=0;
                    
                    switch (m_model) {
                        case 3:
                            DevGamma( Dmu1,  gamma[ip], mu_sand[ip],  vParticlePressure[ip], vVelocityGrad[ip], mu_a, rho_s, dp, FricMu_1,  FricMu_2,  I_0,  deltaGamma,  deltaPs,  deltaI, Inter);
                            //deriv1 =  VolFraction * pow(  mu_eins[ip] , 2) / pow(  mu_eins[ip] + (1-VolFraction)*mu_sand[ip]  ,2);
                            deriv1=1.0;
                            MatScale(Mu_deriv, deriv1 / vMixDensity[ip] ,Dmu1);
                            
                            break;
                        case 4:
                            if ( VolFraction>interface_volume_fraction)
                                DevGamma( Dmu1,  gamma[ip], mu_sand[ip],  vParticlePressure[ip], vVelocityGrad[ip], mu_a, rho_s, dp, FricMu_1,  FricMu_2,  I_0,  deltaGamma,  deltaPs,  deltaI, Inter);
                            deriv1 =  1.0;
                            MatScale(Mu_deriv, deriv1 / vMixDensity[ip] ,Dmu1);
                            break;
                        default:
                            UG_THROW("Error in GranularViscosityLinker model= 2, 3")
                            break;
                    }
                    
                    
                    for(size_t fct = 0; fct < m_spDVelocityGrad->num_fct(); ++fct)
                    {
                    //    get derivative of velocity gradient w.r.t. to all functions
                        const MathMatrix<dim,dim>* vDVelocityGrad = m_spDVelocityGrad->deriv(s_DV_, ip, fct);

                    //    get common fct id for this function
                        const size_t commonFct = this->input_common_fct(_DV_, fct);
                        //printf("num_fct  %d  \n",m_spDVelocityGrad->num_fct());

                    //    loop all shapes and set the derivative
                        for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                        {
                            

                            Inter->MatAddTraspose(Deriv,vDVelocityGrad[sh]);
                            vvvDeriv[ip][commonFct][sh] += MatMultiplyElment(Mu_deriv,Deriv) ;

                            
                        }
                    }
                }
            }*/
             
            
            
        }

    protected:
    
        static number MatMultiplyElment( const MathMatrix<dim,dim> m1, const MathMatrix<dim,dim> m2)
        {
            number Sum=0.0;
            for(size_t i = 0; i < dim; ++i)
                for(size_t j = 0; j < dim; ++j)
                {
                    Sum += m1[i][j] * m2[i][j];
                }
            return Sum;
        }


        static number Constant_viscosity(const number phi, const number mu_1, const number mu_2, const number interface_volume_fraction)
        {
            if (phi>=interface_volume_fraction) return (mu_2);
            else return (mu_1);
            //return 1/(phi/mu_1+(1-phi)/mu_2);
            //return (phi*mu_1+(1-phi)*mu_2);
        }
        static number Proportional_viscosity(const number phi, const number mu_1, const number mu_2)
        {
            const number MU = (phi*mu_2+(1-phi)*mu_1);
            if(std::isnan(MU)) UG_THROW("Error in Linear ViscosityLinker: NaN");
            return MU;
        }
        
        static void Granular_viscosity_1(number& mu_s, const MathMatrix<dim,dim> VelocityGrad, const number Ps, const number rho_s, const number dp, const number FricMu_1, const number FricMu_2, const number I_0, const number deltaGamma, const number deltaPs, const number deltaI, number& gamma)
        {
             
            //number grad_vel_mag=0.0;
            number I;
            number mu_friction;
            //number St;

            gamma=0.0;
            //number grad_vel=0;
            
            // compute inner sum
            for(int d1 = 0; d1 < dim; ++d1)
            {
                for(int d2 = 0; d2 < dim; ++d2)
                {
                    gamma += pow(VelocityGrad(d1,d2) + VelocityGrad(d2,d1),2);
                    //grad_vel += pow(VelocityGrad(d1,d2) + VelocityGrad(d2,d1)-div_factor*(Div1),2);
                }
            }
            
            

            gamma=   sqrt( pow(deltaGamma,2) + 0.5 * gamma    );
            //grad_vel=sqrt( pow(deltaGamma,2) + 0.5 * grad_vel );
            
            /*if (false)
            {
                printf("Gamma   =  %lf\n",gamma);
                printf("Grad vel=  %lf\n",grad_vel);
                printf("Grad rel=  %lf\n",gamma-grad_vel);
            }*/
            //grad_vel_mag =sqrt(0.5*grad_vel_mag);
            //grad_vel_mag=gamma;
            
            
            //Stokes number
            //St=gamma*rho_s*pow(dp,2)/mu_a;
            //Permanent contact pressure



            I=gamma*dp/pow((Ps+deltaPs)/rho_s,0.5)+deltaI;
            
            //I=mu_a*(1+St)*gamma/(Ps+deltaPs)+deltaI;//grad_vel_mag*Viscosity_fluid/((Pressure_s+0.001)*ParticleDensity);
            
            mu_friction=FricMu_1+(FricMu_2-FricMu_1)/(1.0+I_0/I);


            mu_s=mu_friction*Ps/gamma;
            
            if(std::isnan(mu_s)) 
				UG_THROW("Error in MU(I) ViscosityLinker: NaN"
                                     <<"\nPressure = "<< Ps
                                     << "\nGamma(I) = " << gamma
                                     << "\nMU(I) = " <<mu_friction <<".");
            
            


        }
        static void Granular_viscosity_2(number& ss, number& mu_s, number& mu_einstein, number& fac,const MathMatrix<dim,dim> VelocityGrad, const number Ps,const number mu_a, const number rho_s, const number dp, const number FricMu_1, const number FricMu_2, const number alpha_min, const number alpha_max, const number I_0, const number deltaGamma, const number deltaPs, const number deltaI, number vol, number& gamma, const number packing_factor, const number power)
        {
             
            //number grad_vel_mag=0.0;
            number I;
            number mu_friction;
            //number St;

            gamma=0.0;

            // compute inner sum
            for(int d1 = 0; d1 < dim; ++d1)
            {
                for(int d2 = 0; d2 < dim; ++d2)
                {
                    gamma += pow(VelocityGrad(d1,d2) + VelocityGrad(d2,d1),2);
                    //grad_vel_mag += pow((VelocityGrad(d1,d2) + VelocityGrad(d2,d1)),2);
                }
            }

            gamma=sqrt(pow(deltaGamma,2)+0.5*gamma);
         
            //grad_vel_mag =sqrt(0.5*grad_vel_mag);
            //grad_vel_mag=gamma;
            
            
            //Stokes number
            //St=gamma*rho_s*pow(dp,2)/mu_a;
            //Permanent contact pressure



            I=gamma*dp/pow((Ps+deltaPs)/rho_s,0.5)+deltaI;
            
            //I=mu_a*(1+St)*gamma/(Ps+deltaPs)+deltaI;//grad_vel_mag*Viscosity_fluid/((Pressure_s+0.001)*ParticleDensity);
            
            mu_friction=FricMu_1+(FricMu_2-FricMu_1)/(1.0+I_0/I);


            mu_s=mu_friction*Ps/gamma+mu_a;
            //Einstein_viscosity(mu_einstein, vol, mu_a, alpha_max, 1e5);
            
            fac=pow(vol,power);
            
            ss=1/((1-fac)/mu_einstein+fac/mu_s);
            
            if(std::isnan(ss)) UG_THROW("Error in MU(I) Mixed ViscosityLinker: NaN");

        }

    
        
        static void DevGamma(MathMatrix<dim,dim>& ss,  const number gamma, const number Mu_s, const number Ps,const MathMatrix<dim,dim> VelGrad, const number mu_a, const number rho_s, const number dp, const number FricMu_1, const number FricMu_2, const number I_0, const number deltaGamma, const number deltaPs, const number deltaI, Interface<dim>* Inter)
        {
            //number DI;
            number I;
            //number Dmu_friction;
            number mu_friction;
            I=gamma*dp/pow((Ps+deltaPs)/rho_s,0.5)+deltaI;
            //DI=dp/pow((Ps+deltaPs)/rho_s,0.5)+deltaI;
            mu_friction=FricMu_1+(FricMu_2-FricMu_1)/(1.0+I_0/I);
            //Dmu_friction=((FricMu_2-FricMu_1)*I_0/pow(I+I_0,2))*DI;
            
            Inter->MatAddTraspose(ss,VelGrad);
            
            //MatScale(ss, (-(Mu_s)+Ps*Dmu_friction)/(2*pow(gamma,2)) ,ss);
            //MatScale(ss, (Ps*Dmu_friction)/(2*pow(gamma,2)) ,ss);
            MatScale(ss, (-Ps*mu_friction)/(pow(gamma,3)) ,ss);

        }
    
        static void DevPressure_1(number& ss,  const number gamma, const number Ps, const number mu_a, const number rho_s, const number dp, const number FricMu_1, const number FricMu_2, const number I_0, const number deltaGamma, const number deltaPs, const number deltaI)
        {

            number I;
            //number DI;
            number mu_friction;
            //number Dmu_friction;
            
            /*
            I=mu_a*(1+St)*gamma/(Ps+deltaPs)+deltaI;//grad_vel_mag*Viscosity_fluid/((Pressure_s+0.001)*ParticleDensity);
            DI=-I/(Ps+deltaPs);*/
            
            I=gamma*dp/pow((Ps+deltaPs)/rho_s,0.5)+deltaI;
            //DI=-0.5*I/(Ps+deltaPs);
            
            mu_friction=FricMu_1+(FricMu_2-FricMu_1)/(1.0+I_0/I);
            //Dmu_friction=(FricMu_2-FricMu_1)*I_0/pow(I+I_0,2);

            //ss=(Ps*Dmu_friction*DI+mu_friction)/gamma;
            ss=(mu_friction)/gamma;

        }
     
        /*void Granular_viscosity_2(number& ss, const number& VolumeFraction, const MathMatrix<dim,dim>& VelocityGrad, const number delta, const number Viscosity_fluid, const number ParticleDensity, const number FrictionCoefficient_1, const number FrictionCoefficient_2)
        {
            number gamma=0.0;
            number grad_vel_mag=0.0;
            number I;
            number mu_friction;
            number Pressure_s;
            number pff;
            number pa;
            number fac;
            // compute inner sum
            for(int d1 = 0; d1 < dim; ++d1)
            {
                for(int d2 = 0; d2 < dim; ++d2)
                {   fac=0.5;
                    if ( d1==d2){fac=1/3;}
                    gamma += pow(fac*(VelocityGrad(d1,d2) + VelocityGrad(d2,d1)),2);
                    grad_vel_mag += pow(VelocityGrad(d1,d2) + VelocityGrad(d2,d1),2);
                }
            }
            gamma =sqrt(delta+(0.5*gamma));
            grad_vel_mag =sqrt(0.5*grad_vel_mag);
            
            pff=(VolumeFraction >= alpha_min) ? Fr*
            pow(        fmin(VolumeFraction,alpha_max-0.025)-alpha_min              ,3) /
            pow(        alpha_max-fmin(VolumeFraction,alpha_max-0.025)              ,5) : 0.0;
            
            pa = pow(B_phi*fmin(alpha_max,fmax(0,VolumeFraction)) /fmax(0.1*delta,alpha_max-fmin(alpha_max,fmax(0,VolumeFraction))),2)*Viscosity_fluid*gamma;

            Pressure_s = pff+pa+delta;
            I=grad_vel_mag*Viscosity_fluid/(Pressure_s*ParticleDensity);
            
            
            
            mu_friction=FrictionCoefficient_1+(FrictionCoefficient_2-FrictionCoefficient_1)/(1.0+I_0/I);
            //viscosity_mix=Viscosity_fluid*pow(1.0-alpha_limit/0.64,-2.5*0.64);


            ss=mu_friction*Pressure_s/(gamma)+Viscosity_fluid;
            

            if (ss<0){
                printf("-------------\n\n");
                printf("Viscosity  %f\n",ss);
                printf("Gamma      %f\n",gamma);
                printf("Mu         %f\n",mu_friction);
                printf("pff         %f\n",pff);
                printf("pa         %f\n",pa);
            }
            if (std::isnan(gamma)){printf("Flag 5");}
            if (std::isnan(mu_friction)){printf("Flag 4");}
            
            if (isinf(Pressure_s)){printf("Flag 3");}
            if (isinf(mu_friction)){printf("Flag 4");}
            if (isinf(gamma)){printf("Flag 5");}
        }*/



    public:
    // Setter functions for imports
    ///    set volume fraction import
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


    /// set velocity gradient import
        void set_velocity_gradient(SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > data)
        {
            m_spVelocityGrad = data;
            m_spDVelocityGrad = data.template cast_dynamic<DependentUserData<MathMatrix<dim,dim>, dim> >();
            base_type::set_input(_DV_, data, data);
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
        void set_particle_pressure(SmartPtr<CplUserData<number, dim> > data)
        {
            m_spParticlePressure = data;
            m_spDParticlePressure = data.template cast_dynamic<DependentUserData<number, dim> >();
            base_type::set_input(_P_, data, data);
        }
        void set_particle_pressure(number val)
        {
            set_particle_pressure(make_sp(new ConstUserNumber<dim>(val)));
        }
        void set_eins_viscosity(SmartPtr<CplUserData<number, dim> > data)
        {
            m_spMixViscosity = data;
            m_spDMixViscosity = data.template cast_dynamic<DependentUserData<number, dim> >();
            base_type::set_input(_MU_, data, data);
        }

        void set_eins_viscosity(number val)
        {
			set_eins_viscosity(make_sp(new ConstUserNumber<dim>(val)));
        }

    protected:
         //  variables for storing imports
        ///    import for volume fraction
            static const size_t _VOL_ = 0;
            SmartPtr<CplUserData<number, dim> > m_spVolumeFraction;
            SmartPtr<DependentUserData<number, dim> > m_spDVolumeFraction;

        ///    import for velocity gradient
            static const size_t _DV_ = 1;
            SmartPtr<CplUserData<MathMatrix<dim,dim>, dim>> m_spVelocityGrad;
            SmartPtr<DependentUserData<MathMatrix<dim,dim>,dim> > m_spDVelocityGrad;
    
        ///    import for density
            static const size_t _RHO_ = 2;
            SmartPtr<CplUserData<number, dim> > m_spMixDensity;
            SmartPtr<DependentUserData<number, dim> > m_spDMixDensity;
    
        ///    import for density
            static const size_t _P_ = 3;
            SmartPtr<CplUserData<number, dim> > m_spParticlePressure;
            SmartPtr<DependentUserData<number, dim> > m_spDParticlePressure;
    
        ///    import for volume fraction
            static const size_t _MU_ = 4;
            SmartPtr<CplUserData<number, dim> > m_spMixViscosity;
            SmartPtr<DependentUserData<number, dim> > m_spDMixViscosity;
    
            Interface<dim>* Inter;
    
    
    
        public:
            void set_granular_model(int model) {
                m_model = model;
            }
            void set_particle_density(float R) {
                rho_s = R;
            }
            void set_particle_kinematicVisc(float R) {
                nu_s = R;
            }
            void set_fluid_density(float R) {
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
            void set_interface_volume_fraction(float R) {
                interface_volume_fraction = R;
            }
            void set_deltaPs(float R) {
                deltaPs = R;
            }
            void set_deltaI(float R) {
                deltaI = R;
            }
            void set_FricMu_1(float R) {
                FricMu_1 = R;
            }
            void set_FricMu_2(float R) {
                FricMu_2 = R;
            }
            void set_I_0(float R) {
                I_0 = R;
            }
            void set_deltaGamma(float R) {
                deltaGamma = R;
            }
        protected:

            int m_model;
            float rho_s, nu_s, rho_a, mu_a, I_0, FricMu_1, FricMu_2, dp, alpha_max, alpha_min, deltaI, deltaPs, deltaGamma;
            float interface_volume_fraction; //m_packing_factor;
    
};

} // end of namespace ug

#endif // __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_VISCOSITY_LINKER__
