
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_DIFFUSION_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_DIFFUSION_LINKER__

#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "../properties_interface.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// Granular Diffusion linker
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
class GranularDiffusionLinker
    : public StdDataLinker< GranularDiffusionLinker<dim>, MathMatrix<dim,dim>, dim>
{
    //    Base class type
    typedef StdDataLinker< GranularDiffusionLinker<dim>, MathMatrix<dim,dim>, dim> base_type;

    //  Constructor
    public:
    GranularDiffusionLinker() :
            m_spGamma(NULL), m_spDGamma(NULL),
            m_Diff_factor(0.0),
			m_BoolConsKinVisc(true)
        {
        //    this linker needs exactly four input
            this->set_num_input(3);
        }

        // function for evaluation at single ip?
        inline void evaluate (MathMatrix<dim,dim>& value,
                            const MathVector<dim>& globIP,
                            number time, int si) const
        {
            UG_LOG("GranularDiffusionLinker::evaluate single called");
            // Compute second invariant
            //number gamma = 0.0;
            //number I = 0.0;
            //number mu_friction = 0.0;
            //number viscosity = 0.0;
            

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

            MatSet(value,1.0);
        }

        // function for evaluation at multiple ips??
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
            std::vector<number> vGamma(nip);
			std::vector<MathMatrix<dim,dim> > vVelocityGrad(nip);
			std::vector<number> vMixViscosity(nip);

            (*m_spGamma)(&vGamma[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
			(*m_spVelocityGrad)(&vVelocityGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
			(*m_spMixViscosity)(&vMixViscosity[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
			const number nu_a = (Inter-> Viscosity_a()) / Inter->Density_a();

            for(size_t ip = 0; ip < nip; ++ip)
            {
                
                MathMatrix<dim,dim> Diff;
                MatSet(Diff,0.0);
				if (m_BoolConsKinVisc)
				{
					MatDiagSet(Diff,m_Diff_factor * nu_a * vGamma[ip]);
					vValue[ip]=Diff;
				}
				else
				{
					number gamma=0.0;
					// compute inner sum
					for(int d1 = 0; d1 < dim; ++d1)
					{
						for(int d2 = 0; d2 < dim; ++d2)
						{
							gamma += pow((vVelocityGrad[ip](d1,d2) + vVelocityGrad[ip](d2,d1)),2);
						}
					}
					
					MatDiagSet(Diff,m_Diff_factor * vMixViscosity[ip] * gamma);
					vValue[ip]=Diff;
				}

            }
        }

        template <int refDim>
        void eval_and_deriv(MathMatrix<dim,dim> vValue[],
                             const MathVector<dim> vGlobIP[],
                             number time, int si,
                             GridObject* elem,
                             const MathVector<dim> vCornerCoords[],
                             const MathVector<refDim> vLocIP[],
                             const size_t nip,
                             LocalVector* u,
                             bool bDeriv,
                             int s,
                             std::vector<std::vector<MathMatrix<dim,dim> > > vvvDeriv[],
                             const MathMatrix<refDim, dim>* vJT = NULL) const
        {
        
            
            int s_GAMMA_ = base_type::series_id(_GAMMA_, s);
			int s_DV_  = base_type::series_id(_DV_ , s);
			int s_MU_ = base_type::series_id(_MU_, s);
           
        //    get the data of the ip series
            const number* 				vGamma 			= m_spGamma->values(s_GAMMA_);
			const MathMatrix<dim,dim>*  vVelocityGrad	= m_spVelocityGrad->values(s_DV_);
			const number*               vMixViscosity	= m_spMixViscosity->values(s_MU_);

            //number viscosity_granular;
            //number VolFraction;
            
            
            /*for(size_t ip = 0; ip < nip; ++ip)
            {
                VolFraction+=vMixKinViscosity[ip];
                
            }
            //VolFraction=;
            VolFraction=fmin(1.0, fmax(VolFraction/nip,0));*/
			
			const number nu_a = (Inter-> Viscosity_a()) / Inter->Density_a();
            for(size_t ip = 0; ip < nip; ++ip)
            {
                                
                MathMatrix<dim,dim> Diff;
                MatSet(Diff,0.0);
				if(m_BoolConsKinVisc)
				{
					MatDiagSet(Diff,m_Diff_factor * nu_a * vGamma[ip]);
					vValue[ip]=Diff;
				}
				else
				{
					
					number gamma=0.0;
					// compute inner sum
					for(int d1 = 0; d1 < dim; ++d1)
					{
						for(int d2 = 0; d2 < dim; ++d2)
						{
							gamma += pow((vVelocityGrad[ip](d1,d2) + vVelocityGrad[ip](d2,d1)),2);
						}
					}
					
					gamma =sqrt((0.5*gamma));
					
					MatDiagSet(Diff,m_Diff_factor * vMixViscosity[ip] * gamma);
					
					vValue[ip]=Diff;
				}
            }
            
            
            /*size_t _C_=(*u).num_all_fct() -1 ;
            size_t numSH=(*u).num_all_dof(_C_);
            bool f= true;
            for(size_t ip = 0; ip < numSH; ++ip)
            {
            
                if (!((vCornerCoords[ip][0] > 17.374) && (vCornerCoords[ip][0] < 17.626)) && !((vCornerCoords[ip][1] > 1.99484) && (vCornerCoords[ip][1] < 2.2166)))
                {
                    f = f && false;
                }
            }
            for(size_t ip = 0; ip < nip; ++ip)
            {
                if (f)
                {
                    if(ip==0) printf("------------------------------------------------- Diff at linker\n");
                    printf("------------------------------------------------- Tau[%zu] = %f\n",ip, Tau[ip]);
                }
            }*/

                
            

            //    Compute the derivatives at all ips     //
            /////////////////////////////////////////////

        //    check if something is left to do
            if(!bDeriv || this->zero_derivative()) return;

        //    clear all derivative values
            this->set_zero(vvvDeriv, nip);
            /*
        //  Derivatives of Volume Fraction
            if(m_spDMixKinViscosity.valid() && !m_spDMixKinViscosity->zero_derivative())
            {
                for(size_t ip = 0; ip < nip; ++ip)
                {
                    VolFraction=fmin(1.0, fmax(vMixKinViscosity[ip],0));
                    
                    for(size_t fct = 0; fct < m_spDMixKinViscosity->num_fct(); ++fct)
                    {
                    //    get derivative of volume fraction w.r.t. to all functions
                        const number* vDMixKinViscosity = m_spDMixKinViscosity->deriv(s, ip, fct);

                    //    get common fct id for this function
                        const size_t commonFct = this->input_common_fct(_MU_, fct);
                        number Dmu_einstein;
                        number Dmu_rheology;

                    //    loop all shapes and set the derivative
                        for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                        {
                            switch (m_model) {
                                case 0:
                                    vvvDeriv[ip][commonFct][sh]=vDMixKinViscosity[sh]*(nu_s*rho_s-nu_a*rho_s);
                                    break;
                                case 1:
                                    vvvDeriv[ip][commonFct][sh]=vDMixKinViscosity[sh]*(2.5*m_packing_factor*vValue[ip]/(1-vMixKinViscosity[ip]*m_packing_factor/alpha_max));
                                    break;
                                case 2:
                                    DevGranular_viscosity_1(Dmu_rheology,VolFraction*m_packing_factor, vVelocityGrad[ip], nu_a*rho_a, rho_s, dp, FricMu_1, FricMu_2, alpha_min, alpha_max, Fr, B_phi, I_0, nu_a, deltaPs, deltaI);
                                    vvvDeriv[ip][commonFct][sh]=vDMixKinViscosity[sh]*Dmu_rheology*m_packing_factor;
                                    break;
                                case 3:
                                    Einstein_viscosity(Dmu_einstein,VolFraction*m_packing_factor, nu_a*rho_a, alpha_max);
                                    Granular_viscosity_1(Dmu_rheology,VolFraction*m_packing_factor, vVelocityGrad[ip], nu_a*rho_a, rho_s, dp, FricMu_1, FricMu_2, alpha_min, alpha_max, Fr, B_phi, I_0, nu_a, deltaPs, deltaI);
                                    vvvDeriv[ip][commonFct][sh]=vDMixKinViscosity[sh]*(Dmu_rheology-Dmu_einstein+(1-VolFraction)*2.5*m_packing_factor*Dmu_einstein/(1-VolFraction*m_packing_factor/alpha_max));
                                    break;
                                default:
                                    UG_THROW("Wrong model selected for granular viscosity: GranularViscosityLinker has options"
                                                   " model= 0, 1 , 2, 3")
                                    break;
                            }
                            

                            

                        }
                    }
                }
            }
            */
        
            /*
            //  Derivatives of velocity gradient
            //UG_LOG("Derivatives of velocity gradient missing");
            if(m_spDVelocityGrad.valid() && !m_spDVelocityGrad->zero_derivative())
            {
                for(size_t ip = 0; ip < nip; ++ip)
                {
                    for(size_t fct = 0; fct < m_spDYieldStress->num_fct(); ++fct)
                    {
                    //    get derivative of velocity gradient w.r.t. to all functions
                        const number* vDYieldStress = m_spDYieldStress->deriv(s, ip, fct);

                    //    get common fct id for this function
                        const size_t commonFct = this->input_common_fct(_DV_, fct);

                    //    loop all shapes and set the derivative
                        for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                        {
                            VecScaleAppend(vvvDeriv[ip][commonFct][sh], vDYieldStress[sh], rInvariant);
                        }
                    }
                }
            }
             
            
        
            */
            
        }


    public:
    // Setter functions for imports
    
    ///    set gamma import
        void set_gamma(SmartPtr<CplUserData<number, dim> > data)
        {
            m_spGamma = data;
            m_spDGamma = data.template cast_dynamic<DependentUserData<number, dim> >();
            base_type::set_input(_GAMMA_, data, data);
        }
        void set_gamma(number val)
        {
            set_gamma(make_sp(new ConstUserNumber<dim>(val)));
        }
	
	/// set velocity gradient import
		void set_velocity_gradient(SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > data)
		{
			m_spVelocityGrad = data;
			m_spDVelocityGrad = data.template cast_dynamic<DependentUserData<MathMatrix<dim,dim>, dim> >();
			base_type::set_input(_DV_, data, data);
		}
	
		void set_mix_viscosity(SmartPtr<CplUserData<number, dim> > data)
		{
			m_spMixViscosity = data;
			m_spDMixViscosity = data.template cast_dynamic<DependentUserData<number, dim> >();
			base_type::set_input(_MU_, data, data);
		}

		void set_mix_viscosity(number val)
		{
			set_mix_viscosity(make_sp(new ConstUserNumber<dim>(val)));
		}
	
        void set_diff_factor(float R)
        {
            if(R < 0.0)
                UG_THROW("GranularDiffusionLinker:  Diffusion Factor must be greater or equal zero.");
            m_Diff_factor = R;
            
        }
		void const_kinematic_visc(bool R) {
			m_BoolConsKinVisc = R;
		}

		void set_phase_parameters(Interface<dim>* user)
		{
			if (!user->valid())
				UG_THROW("Interface parameters has not been initialized");
			Inter = user;
		}

    protected:
         //  variables for storing imports
    
        ///    import for gamma
            static const size_t _GAMMA_ = 0;
            SmartPtr<CplUserData<number, dim> > m_spGamma;
            SmartPtr<DependentUserData<number, dim> > m_spDGamma;
	
		///    import for velocity gradient
			static const size_t _DV_ = 1;
			SmartPtr<CplUserData<MathMatrix<dim,dim>, dim>> m_spVelocityGrad;
			SmartPtr<DependentUserData<MathMatrix<dim,dim>,dim> > m_spDVelocityGrad;
	
		///    import for volume fraction
			static const size_t _MU_ = 2;
			SmartPtr<CplUserData<number, dim> > m_spMixViscosity;
			SmartPtr<DependentUserData<number, dim> > m_spDMixViscosity;
	
			Interface<dim>* Inter;
	
    protected:
            float m_Diff_factor;
			bool m_BoolConsKinVisc;
    
};

} // end of namespace ug

#endif // __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_VISCOSITY_LINKER__
