
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_PRESSURE_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_PRESSURE_LINKER__

#include "lib_disc/spatial_disc/user_data/linker/linker.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// Scale Pressure linker
////////////////////////////////////////////////////////////////////////////////

/// Linker for the Scale pressure
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
class ScalePressureLinker
    : public StdDataLinker< ScalePressureLinker<dim>, number, dim>
{
    //    Base class type
    typedef StdDataLinker< ScalePressureLinker<dim>, number, dim> base_type;

    //  Constructor
    public:
        ScalePressureLinker() :
            m_spVolumeFraction(NULL), m_spDVolumeFraction(NULL),
            m_spPressure(NULL), m_spDPressure(NULL)
        {
        //    this linker needs exactly four input
            this->set_num_input(2);
        }

        // function for evaluation at single ip?
        inline void evaluate (number& value,
                            const MathVector<dim>& globIP,
                            number time, int si) const
        {
            UG_LOG("ScalePressureLinker::evaluate single called");
            number volume_fraction;

            (*m_spVolumeFraction)(volume_fraction, globIP, time, si);

            // Compute second invariant
            number gamma = 0.0;
 
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

            
            
            //value = viscosity + yieldStress/;
            // compute mu = eta + sigma / \sqrt( delta + 1 / I)
            value = gamma*volume_fraction*viscosity + (1.0-volume_fraction) ;
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
            UG_LOG("ScalePressureLinker::evaluate single called");
            std::vector<number> vVolumeFraction(nip);


            (*m_spVolumeFraction)(&vVolumeFraction[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);

            
            if (true){
                printf("Flag 1");
            }

            for(size_t ip = 0; ip < nip; ++ip)
            {
                // Compute second invariant
                number gamma=0.0;
                number viscosity=0.0;
                

                // compute inner sum
                for(int d1 = 0; d1 < dim; ++d1)
                {
                    for(int d2 = 0; d2 < dim; ++d2)
                    {
                        //gamma += pow(vVelocityGrad[ip](d1,d2) + vVelocityGrad[ip](d2,d1),2);
                    }
                }
                gamma = 0.0001+sqrt((0.5*gamma));
                
                viscosity = 1/gamma;
                
                // compute mu = eta + sigma
                vValue[ip] = vVolumeFraction[ip]*viscosity+(1.0-vVolumeFraction[ip]);
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
        
        //UG_LOG("ParticlePressureLinker::eval_and_deriv called");
        //    get the data of the ip series
            const number* vVolumeFraction = m_spVolumeFraction->values(s);
            const number* vPressure = m_spPressure->values(s);
           


  
            
            number VolFraction;
            for(size_t ip = 0; ip < nip; ++ip)
            {

                VolFraction=0.6*fmin(1.0, fmax(vVolumeFraction[ip],0));

                vValue[ip]=-vPressure[ip]*VolFraction;
                

            }
            

            //    Compute the derivatives at all ips     //
            /////////////////////////////////////////////

        //    check if something is left to do
            if(!bDeriv || this->zero_derivative()) return;

        //    clear all derivative values
            this->set_zero(vvvDeriv, nip);
            
        //  Derivatives of Volume Fraction
            if(m_spDVolumeFraction.valid() && !m_spDVolumeFraction->zero_derivative())
            {
                for(size_t ip = 0; ip < nip; ++ip)
                {
                                        
                    for(size_t fct = 0; fct < m_spDVolumeFraction->num_fct(); ++fct)
                    {
                    //    get derivative of volume fraction w.r.t. to all functions
                        const number* vDVolumeFraction = m_spDVolumeFraction->deriv(s, ip, fct);

                    //    get common fct id for this function
                        const size_t commonFct = this->input_common_fct(_VOL_, fct);

                    //    loop all shapes and set the derivative
                        for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                        {

                            vvvDeriv[ip][commonFct][sh] += vDVolumeFraction[sh]*(vPressure[ip]);

                        }
                    }
                }
            }
            
        //  Derivatives of Volume Fraction
            if(m_spDPressure.valid() && !m_spDPressure->zero_derivative())
            {
                for(size_t ip = 0; ip < nip; ++ip)
                {
                    VolFraction=fmin(1.0, fmax(vVolumeFraction[ip],0));
                    
                    for(size_t fct = 0; fct < m_spDPressure->num_fct(); ++fct)
                    {
                    //    get derivative of volume fraction w.r.t. to all functions
                        const number* vDPressure = m_spDPressure->deriv(s, ip, fct);

                    //    get common fct id for this function
                        const size_t commonFct = this->input_common_fct(_P_, fct);

                    //    loop all shapes and set the derivative
                        for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                        {

                            vvvDeriv[ip][commonFct][sh] += vDPressure[sh]*(VolFraction);

                        }
                    }
                }
            }
             
            

 /*
            static void DevPressure_1(number& ss, const MathMatrix<dim,dim> VelocityGrad, const number mu_a, const number rho_s, const number dp, const number FricMu_1, const number FricMu_2, const number alpha_min, const number alpha_max, const number Fr, const number B_phi, const number I_0, const number DeltaGamma, const number deltaPs, const number deltaI)
            {
                number gamma=0.0;
                number grad_vel_mag=0.0;
                number I;
                number DI;
                number mu_friction;
                number Dmu_friction;
                number Ps;
                number DPs;
                number pff;
                number Dpff;
                number pa;
                number Dpa;
                number St;
                number Div=0;
                number fac;
                for(int d1 = 0; d1 < dim; ++d1)
                    Div += VelocityGrad(d1,d1);
                // compute inner sum
                for(int d1 = 0; d1 < dim; ++d1)
                {
                    for(int d2 = 0; d2 < dim; ++d2)
                    {   fac=0;
                        if ( d1==d2){ fac=2/3;}
                            
                        gamma += pow(VelocityGrad(d1,d2) + VelocityGrad(d2,d1)-fac*(Div),2);
                        grad_vel_mag += pow((VelocityGrad(d1,d2) + VelocityGrad(d2,d1)),2);
                    }
                }

                gamma =sqrt(DeltaGamma+(0.5*gamma));
                grad_vel_mag =sqrt(0.5*grad_vel_mag);
                
                
                
                //Stokes number
                St=gamma*rho_s*pow(dp,2)/mu_a;
                //Permanent contact pressure
                pff=(phi >= alpha_min) ? Fr *pow(  phi-alpha_min,3) /pow(alpha_max-phi,5) : 0.0;
                Dpff=(phi >= alpha_min) ? pff*(2*phi+3*alpha_max-5*alpha_min)/((  phi-alpha_min)*(alpha_max-phi)) : 0;
                //Dynamic pressure
                pa = mu_a*(1.0+St)*pow(B_phi*phi/(alpha_max-phi),2)*gamma;
                Dpa = 2*mu_a*(1.0+St)*pow(B_phi,2)*phi*alpha_max*gamma/pow(alpha_max-phi,3);//pow(B_phi*vol /(alpha_max-vol),2)*Viscosity_fluid*gamma;
                Ps=pff +pa;
                DPs = Dpff+Dpa;
                
                
                I=mu_a*(1+St)*gamma/(Ps+deltaPs)+deltaI;//grad_vel_mag*Viscosity_fluid/((Pressure_s+0.001)*ParticleDensity);
                DI=-I/(Ps+deltaPs);
                mu_friction=FricMu_1+(FricMu_2-FricMu_1)/(1.0+I_0/I);
                Dmu_friction=(FricMu_2-FricMu_1)*I_0/pow(I+I_0,2);

                ss=(DPs/gamma)*(Ps*Dmu_friction*DI+mu_friction);
                
            }*/
        }

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

    
    ///    set volume fraction import
        void set_pressure(SmartPtr<CplUserData<number, dim> > data)
        {
            m_spPressure = data;
            m_spDPressure = data.template cast_dynamic<DependentUserData<number, dim> >();
            base_type::set_input(_P_, data, data);
        }

        void set_pressure(number val)
        {
            set_pressure(make_sp(new ConstUserNumber<dim>(val)));
        }




    protected:
         //  variables for storing imports
        ///    import for volume fraction
            static const size_t _VOL_ = 0;
            SmartPtr<CplUserData<number, dim> > m_spVolumeFraction;
            SmartPtr<DependentUserData<number, dim> > m_spDVolumeFraction;
    
        ///    import for  pressure
            static const size_t _P_ = 1;
            SmartPtr<CplUserData<number, dim> > m_spPressure;
            SmartPtr<DependentUserData<number, dim> > m_spDPressure;


};

} // end of namespace ug

#endif // __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_VISCOSITY_LINKER__
