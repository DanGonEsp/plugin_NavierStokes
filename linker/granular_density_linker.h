
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_DENSITY_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_DENSITY_LINKER__

#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// Granular Density linker
////////////////////////////////////////////////////////////////////////////////

/// Linker for the Granular density
/**
 * This linker computes the Bingham density \f$ \mathbf{q} = - \frac{\mathbf{K}}{\mu} ( \nabla p - \rho \mathbf{g} ) \f$,
 * where
 * <ul>
 * <li> \f$ \mathbf{K} \f$        permeability
 * <li> \f$ \mu \f$                density
 * <li> \f$ \rho \f$            density
 * <li> \f$ \mathbf{g} \f$        gravity
 * <li> \f$ \nabla p \f$        pressure gradient
 * </ul>
 * are input parameters.
 */

template <int dim>
class GranularDensityLinker
    : public StdDataLinker< GranularDensityLinker<dim>, number, dim>
{
    //    Base class type
    typedef StdDataLinker< GranularDensityLinker<dim>, number, dim> base_type;

    //  Constructor
    public:
        GranularDensityLinker() :
            m_spVolumeFraction(NULL), m_spDVolumeFraction(NULL),
            m_model(0)
        {
        //    this linker needs exactly four input
            this->set_num_input(1);
        }

        // function for evaluation at single ip?
        inline void evaluate (number& value,
                            const MathVector<dim>& globIP,
                            number time, int si) const
        {
            UG_LOG("GranularDensityLinker::evaluate single called");
            number volume_fraction;

            (*m_spVolumeFraction)(volume_fraction, globIP, time, si);
            printf("GranularDensityLinker::evaluate called");
            // Compute second invariant
            //number s=fmin(1.0, fmax(volume_fraction,0));
			
			number rho_s = Inter->Density_s();
			number rho_a = Inter->Density_a();
			number packing_factor = Inter->packing_factor();
			number s= packing_factor * volume_fraction;
            value = s*rho_s + (1.0-s)*rho_a ;
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
            (*m_spVolumeFraction)(&vVolumeFraction[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);
            
			const number rho_s = Inter->Density_s();
			const number rho_a = Inter->Density_a();
			const number interface_value = Inter->interface_value();
			const number packing_factor = Inter->packing_factor();
			number c[nip];
            number value;
            for(size_t ip = 0; ip < nip; ++ip)
            {
				
				//c[ip] = packing_factor * vVolumeFraction[ip];
				c[ip]=packing_factor*(vVolumeFraction[ip] + sqrt(vVolumeFraction[ip]*vVolumeFraction[ip] + m_eps*m_eps))/2.0;
				
                switch (m_model) {
                    case 0:
                        value =Constant_density( c[ip] ,   rho_a,   rho_s, packing_factor * interface_value);

                        break;

                    case 1:
                        value =Linear_density( c[ip] ,   rho_a,   rho_s);


                        break;
                    default:
                        UG_THROW("No model for Granular Density Linker: model= 0, 1")
                        break;
                }
                
                vValue[ip]=value;
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
            const number* vVolumeFraction   = m_spVolumeFraction->values(s_VOL_);
            
			const number rho_s = Inter->Density_s();
			const number rho_a = Inter->Density_a();
			const number interface_value = Inter->interface_value();
			const number packing_factor = Inter->packing_factor();
			number c[nip];
            number value;
            for(size_t ip = 0; ip < nip; ++ip)
            {
				//c[ip] = packing_factor * vVolumeFraction[ip];
				c[ip]=packing_factor*(vVolumeFraction[ip] + sqrt(vVolumeFraction[ip]*vVolumeFraction[ip] + m_eps*m_eps))/2.0;
                

                switch (m_model) {
                    case 0:
                        value =Constant_density( c[ip],   rho_a,   rho_s, packing_factor * interface_value);

                        break;

                    case 1:
                        value =Linear_density( c[ip],   rho_a,   rho_s );
                        
                        break;
                    default:
                        UG_THROW("No model for Granular Density Linker: model= 0, 1")
                        break;
                }
                
                vValue[ip]=value;
				if(std::isnan(value) || value<0.0) UG_THROW("Error in DensityLinker: NaN" << value<< "   phi = "<< c);
            }

            //    Compute the derivatives at all ips     //
            /////////////////////////////////////////////

        //    check if something is left to do
            if(!bDeriv || this->zero_derivative()) return;

        //    clear all derivative values
            this->set_zero(vvvDeriv, nip);
            
        //  Derivatives of Volume Fraction
            if(m_spDVolumeFraction.valid() && !m_spDVolumeFraction->zero_derivative() && m_model == 1)
            {
                for(size_t ip = 0; ip < nip; ++ip)
                {
					//const number vol2 = vVolumeFraction[ip]*vVolumeFraction[ip];
					const number deriv_value = packing_factor *(rho_s-rho_a);// (c[ip]/sqrt(vol2 + m_eps*m_eps));
					
                    for(size_t fct = 0; fct < m_spDVolumeFraction->num_fct(); ++fct)
                    {
                    //    get derivative of volume fraction w.r.t. to all functions
                        const number* vDVolumeFraction = m_spDVolumeFraction->deriv(s_VOL_, ip, fct);

                    //    get common fct id for this function
                        const size_t commonFct = this->input_common_fct(_VOL_, fct);

                    //    loop all shapes and set the derivative
                        for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                        {
                            for(size_t ip2 = 0; ip2 < nip; ++ip2)
                            {
                                vvvDeriv[ip2][commonFct][sh]=  deriv_value*vDVolumeFraction[sh];
                            }

                        }
                    }
                }
            }

        }
    protected:
        static number Constant_density(const number phi, const number rho_1, const number rho_2, const number interface)
        {
            if (phi>=interface)
                return (rho_2);
            else
                return (rho_1);
                
        }
        static number Linear_density(const number phi, const number rho_1, const number rho_2)
        {
            return (phi*rho_2+(1-phi)*rho_1);
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


    protected:
         //  variables for storing imports
        ///    import for volume fraction
            static const size_t _VOL_ = 0;
            SmartPtr<CplUserData<number, dim> > m_spVolumeFraction;
            SmartPtr<DependentUserData<number, dim> > m_spDVolumeFraction;
	
			Interface<dim>* Inter;
    
    
        public:

            void set_model(std::string density_model)
            {
                std::string n = TrimString(density_model);
                std::transform(n.begin(), n.end(), n.begin(), ::tolower);
                
                if (n == "constant")
                    m_model = 0;
                else if (n == "linear")
                    m_model = 1;
                else
                    UG_THROW("Density calculation method not found."
                             " Use one of [constant, linear].");
            }
			void set_phase_parameters(Interface<dim>* user)
			{
				if (!user) UG_THROW("Interface pointer is null!");
				if (!user->valid())
					UG_THROW("Interface parameters has not been initialized");
				Inter = user;
			}
        protected:
            int m_model;
			float m_eps = 1e-09;
    
};

} // end of namespace ug

#endif // __H__UG__LIB_DISC__SPATIAL_DISC__GRANULAR_DENSITY_LINKER__






