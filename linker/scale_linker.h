
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_LINKER__

#include "lib_disc/spatial_disc/user_data/linker/linker.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// Scale linker
////////////////////////////////////////////////////////////////////////////////

/// Linker for the Scale


template <int dim>
class ScaleLinker
    : public StdDataLinker< ScaleLinker<dim>, number, dim>
{
    //    Base class type
    typedef StdDataLinker< ScaleLinker<dim>, number, dim> base_type;

    //  Constructor
    public:
        ScaleLinker() :
            m_spImport1(NULL), m_spDImport1(NULL),
            m_spImport2(NULL), m_spDImport2(NULL)
        {
        //    this linker needs exactly four input
            this->set_num_input(2);
        }

        // function for evaluation at single ip?
        inline void evaluate (number& value,
                            const MathVector<dim>& globIP,
                            number time, int si) const
        {
            UG_LOG("ScaleLinker::evaluate single called");
            number volume_fraction;

            (*m_spImport1)(volume_fraction, globIP, time, si);

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

            std::vector<number> vImport1(nip);
            std::vector<number> vImport2(nip);


            (*m_spImport1)(&vImport1[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spImport2)(&vImport2[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);

            

            for(size_t ip = 0; ip < nip; ++ip)
            {
                vValue[ip]=vImport1[ip]*vImport2[ip];
            
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
        
            int s_IMP1_ = base_type::series_id(_IMP1_, s);
            int s_IMP2_ = base_type::series_id(_IMP2_, s);
            

            const number* vImport1 = m_spImport1->values(s_IMP1_);
            const number* vImport2 = m_spImport2->values(s_IMP2_);
                       
  
            
            for(size_t ip = 0; ip < nip; ++ip)
            {
                vValue[ip]=vImport1[ip]*vImport2[ip];
            
            }
            

            //    Compute the derivatives at all ips     //
            /////////////////////////////////////////////

        //    check if something is left to do
            if(!bDeriv || this->zero_derivative()) return;

        //    clear all derivative values
            this->set_zero(vvvDeriv, nip);
            
        //  Derivatives of Volume Fraction
            if(m_spDImport1.valid() && !m_spDImport1->zero_derivative())
            {
                for(size_t ip = 0; ip < nip; ++ip)
                {
                                        
                    for(size_t fct = 0; fct < m_spDImport1->num_fct(); ++fct)
                    {
                    //    get derivative of volume fraction w.r.t. to all functions
                        const number* vDImport1 = m_spDImport1->deriv(s_IMP1_, ip, fct);

                    //    get common fct id for this function
                        const size_t commonFct = this->input_common_fct(_IMP1_, fct);

                    //    loop all shapes and set the derivative
                        for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                        {

                            vvvDeriv[ip][commonFct][sh] += vDImport1[sh]*(vImport2[ip]);

                        }
                    }
                }
            }
            
        //  Derivatives of Volume Fraction
            if(m_spDImport2.valid() && !m_spDImport2->zero_derivative())
            {
                for(size_t ip = 0; ip < nip; ++ip)
                {
                    
                    for(size_t fct = 0; fct < m_spDImport2->num_fct(); ++fct)
                    {
                    //    get derivative of volume fraction w.r.t. to all functions
                        const number* vDImport2 = m_spDImport2->deriv(s_IMP2_, ip, fct);

                    //    get common fct id for this function
                        const size_t commonFct = this->input_common_fct(_IMP2_, fct);

                    //    loop all shapes and set the derivative
                        for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
                        {

                            vvvDeriv[ip][commonFct][sh] += vDImport2[sh]*(vImport1[ip]);

                        }
                    }
                }
            }
             
        }

    public:
    // Setter functions for imports
    ///    set  import 1
        void set_import_1(SmartPtr<CplUserData<number, dim> > data)
        {
            m_spImport1 = data;
            m_spDImport1 = data.template cast_dynamic<DependentUserData<number, dim> >();
            base_type::set_input(_IMP1_, data, data);
        }

        void set_import_1(number val)
        {
            set_import_1(make_sp(new ConstUserNumber<dim>(val)));
        }

    
    ///    set import 2
        void set_import_2(SmartPtr<CplUserData<number, dim> > data)
        {
            m_spImport2 = data;
            m_spDImport2 = data.template cast_dynamic<DependentUserData<number, dim> >();
            base_type::set_input(_IMP2_, data, data);
        }

        void set_import_2(number val)
        {
            set_import_2(make_sp(new ConstUserNumber<dim>(val)));
        }




    protected:
         //  variables for storing imports
        ///    import1
            static const size_t _IMP1_ = 0;
            SmartPtr<CplUserData<number, dim> > m_spImport1;
            SmartPtr<DependentUserData<number, dim> > m_spDImport1;
    
        ///    import2
            static const size_t _IMP2_ = 1;
            SmartPtr<CplUserData<number, dim> > m_spImport2;
            SmartPtr<DependentUserData<number, dim> > m_spDImport2;


};

} // end of namespace ug

#endif // __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_LINKER__
