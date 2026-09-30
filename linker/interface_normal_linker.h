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
            m_spVolumeGrad(NULL), m_spDVolumeGrad(NULL)

        {
        //    this linker needs exactly five input
            this->set_num_input(2);
        }


		inline void evaluate(MathVector<dim>& value,
							 const MathVector<dim>& globIP,
							 number time, int si) const
		{
			MathVector<dim> volumeGrad;
			(*m_spVolumeGrad)(volumeGrad, globIP, time, si);

			compute_normal(value, volumeGrad);
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
			std::vector<MathVector<dim> > vVolumeGrad(nip);

			(*m_spVolumeGrad)(&vVolumeGrad[0], vGlobIP, time, si,
							  elem, vCornerCoords, vLocIP, nip, u, vJT);

			for(size_t ip = 0; ip < nip; ++ip)
				compute_normal(vNormal[ip], vVolumeGrad[ip]);
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
			const int s_DVOL_ = base_type::series_id(_DVOL_, s);

			const MathVector<dim>* vVolumeGrad = m_spVolumeGrad->values(s_DVOL_);

			std::vector<MathMatrix<dim, dim> > vJ(nip);

			for(size_t ip = 0; ip < nip; ++ip)
				compute_normal_and_jacobian(vNormal[ip], vJ[ip], vVolumeGrad[ip]);

			if(!bDeriv || this->zero_derivative())
				return;

			this->set_zero(vvvDeriv, nip);

			if(m_spVolumeGrad->zero_derivative())
				return;

			for(size_t fct = 0; fct < base_type::input_num_fct(_DVOL_); ++fct)
			{
				const size_t commonFct = base_type::input_common_fct(_DVOL_, fct);

				for(size_t ip = 0; ip < nip; ++ip)
				{
					for(size_t sh = 0; sh < m_spDVolumeGrad->num_sh(fct); ++sh)
					{
						const MathVector<dim>& dGrad =
							m_spDVolumeGrad->deriv(s_DVOL_, ip, fct, sh);

						MathVector<dim> dNormal;
						VecSet(dNormal, 0.0);

						for(int i = 0; i < dim; ++i)
							for(int j = 0; j < dim; ++j)
								dNormal[i] += vJ[ip](i,j)*dGrad[j];

						vvvDeriv[ip][commonFct][sh] += dNormal;
					}
				}
			}
		}
	
		inline void compute_normal_and_jacobian(MathVector<dim>& n,
												MathMatrix<dim, dim>& J,
												const MathVector<dim>& grad) const
		{
			const number epsGrad = 1e-4;
			const number epsNorm = 1e-12;

			const number r2 = VecProd(grad, grad);
			const number R = r2 + epsGrad*epsGrad;
			const number s = sqrt(R);
			const number alpha = r2/R;

			MathVector<dim> vertical;
			VecSet(vertical, 0.0);
			vertical[dim-1] = 1.0;

			MathVector<dim> gradReg;
			VecScale(gradReg, grad, -1.0/s);

			MathVector<dim> q;
			VecScaleAdd(q, alpha, gradReg, 1.0-alpha, vertical);

			const number q2 = VecProd(q, q);
			const number Q = sqrt(q2 + epsNorm*epsNorm);

			VecScale(n, q, 1.0/Q);

			MathMatrix<dim, dim> DqDg;

			for(int i = 0; i < dim; ++i)
			{
				for(int j = 0; j < dim; ++j)
				{
					const number dalpha =
						2.0*epsGrad*epsGrad*grad[j]/(R*R);

					const number delta =
						(i == j) ? 1.0 : 0.0;

					// derivative of -grad/s
					const number dGradReg =
						-delta/s + grad[i]*grad[j]/(s*s*s);

					DqDg(i,j) =
						dalpha*(gradReg[i] - vertical[i])
						+ alpha*dGradReg;
				}
			}

			for(int i = 0; i < dim; ++i)
			{
				for(int j = 0; j < dim; ++j)
				{
					J(i,j) = 0.0;

					for(int k = 0; k < dim; ++k)
					{
						const number delta =
							(i == k) ? 1.0 : 0.0;

						const number DnDq =
							delta/Q - q[i]*q[k]/(Q*Q*Q);

						J(i,j) += DnDq*DqDg(k,j);
					}
				}
			}
		}
	
		inline void compute_normal(MathVector<dim>& n, const MathVector<dim>& grad) const
		{
			const number epsGrad = 1e-4;
			const number epsNorm = 1e-12;

			const number r2 = VecProd(grad, grad);
			const number R = r2 + epsGrad*epsGrad;
			const number s = sqrt(R);
			const number alpha = r2 / R;

			// Local interface normal: -grad(c)
			MathVector<dim> gradReg;
			VecScale(gradReg, grad, -1.0/s);

			// Upward fallback for vanishing gradient
			MathVector<dim> vertical;
			VecSet(vertical, 0.0);
			vertical[dim-1] = 1.0;

			MathVector<dim> q;
			VecScaleAdd(q, alpha, gradReg, 1.0-alpha, vertical);

			const number q2 = VecProd(q, q);
			const number qNorm = sqrt(q2 + epsNorm*epsNorm);

			VecScale(n, q, 1.0/qNorm);
		}
    

    public:
    
    ///    set gravity import
        void set_volume_grad(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
        {
            m_spVolumeGrad = data;
            m_spDVolumeGrad = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim> >();
            base_type::set_input(_DVOL_, data, data);
        }

    protected:
    
	///    import for volume fraction grad
        static const size_t _DVOL_ = 0;
        SmartPtr<CplUserData<MathVector<dim>, dim> > m_spVolumeGrad;
        SmartPtr<DependentUserData<MathVector<dim>, dim> > m_spDVolumeGrad;
    


};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__INTERFACE_NORMAL_LINKER__ */

