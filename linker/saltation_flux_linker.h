/*
 * Copyright (c) 2013-2015: G-CSC, Goethe University Frankfurt
 *
 * Saltation flux linker for the Navier-Stokes plugin.
 *
 * The linker computes a diffuse-interface volumetric saltation mass flux
 *
 *      J_salt = delta_Gamma q_s
 *
 * using a saturated saltation closure following the formulation used by
 * Ortiz and Smolarkiewicz.
 *
 * The surface saltation flux is
 *
 *      q_s = (C/g) (u_* - u_*t)_+ tau_Gamma,
 *
 * where
 *
 *      tau_Gamma = 2 mu_a P D n,
 *
 *      D = 1/2 (grad(u) + grad(u)^T),
 *
 *      P = I - n tensor n,
 *
 *      n = -grad(c)/sqrt(|grad(c)|^2 + eps_n^2),
 *
 *      u_* = sqrt(|tau_Gamma|/rho_a).
 *
 * The surface flux is mapped to the diffuse interface using
 *
 *      delta_Gamma = sqrt(|grad(c)|^2 + eps_delta^2) - eps_delta.
 *
 * Therefore,
 *
 *      J_salt = delta_Gamma (C/g) (u_* - u_*t)_+ tau_Gamma.
 *
 * The output has units kg/(m^2 s).
 *
 * Analytical derivatives are provided with respect to:
 *
 *      grad(c)
 *      grad(u)
 *
 * The threshold is treated piecewise:
 *
 *      u_* <= u_*t : J_salt = 0
 *      u_* >  u_*t : active saltation
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__SALTATION_FLUX_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__SALTATION_FLUX_LINKER__

#include <cmath>
#include <vector>

#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "../properties_interface.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Saltation Flux Linker
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class SaltationFluxLinker : public StdDataLinker<SaltationFluxLinker<dim>, MathVector<dim>, dim>
{
	typedef StdDataLinker<SaltationFluxLinker<dim>, MathVector<dim>, dim> base_type;

public:

	SaltationFluxLinker() : m_spVolumeGrad(NULL), m_spDVolumeGrad(NULL), m_spVelocityGrad(NULL), m_spDVelocityGrad(NULL), Inter(NULL), m_C(5.5), m_uStarThreshold(0.22), m_epsNormal(1e-8), m_epsDelta(1e-8)
	{
		this->set_num_input(2);
	}


private:

	////////////////////////////////////////////////////////////////////////////
	// Flux evaluation
	////////////////////////////////////////////////////////////////////////////

	void compute_flux(MathVector<dim>& flux, const MathVector<dim>& gradC, const MathMatrix<dim,dim>& gradU) const
	{
		VecSet(flux, 0.0);

		if(!Inter) UG_THROW("SaltationFluxLinker: Interface pointer is null.");

		const number rho_a = Inter->Density_a();
		const number mu_a = Inter->Viscosity_a();
		const number grav = std::fabs(Inter->gravity());

		if(rho_a <= 0.0) UG_THROW("SaltationFluxLinker: air density must be positive.");
		if(mu_a < 0.0) UG_THROW("SaltationFluxLinker: air viscosity must be non-negative.");
		if(grav <= 0.0) UG_THROW("SaltationFluxLinker: gravity magnitude must be positive.");

		////////////////////////////////////////////////////////////////////////////
		// Diffuse-interface geometry
		////////////////////////////////////////////////////////////////////////////

		number gradC2 = 0.0;
		for(size_t d = 0; d < dim; ++d) gradC2 += gradC[d]*gradC[d];

		const number sN = std::sqrt(gradC2 + m_epsNormal*m_epsNormal);
		const number sDelta = std::sqrt(gradC2 + m_epsDelta*m_epsDelta);
		const number deltaGamma = sDelta - m_epsDelta;

		if(deltaGamma <= 0.0) return;

		MathVector<dim> n;
		for(size_t d = 0; d < dim; ++d) n[d] = -gradC[d]/sN;

		////////////////////////////////////////////////////////////////////////////
		// Air strain-rate tensor
		//
		// D = 1/2 (grad(u) + grad(u)^T)
		////////////////////////////////////////////////////////////////////////////

		MathMatrix<dim,dim> D;

		for(size_t i = 0; i < dim; ++i)
		{
			for(size_t j = 0; j < dim; ++j) D(i,j) = 0.5*(gradU(i,j) + gradU(j,i));
		}

		////////////////////////////////////////////////////////////////////////////
		// D n
		////////////////////////////////////////////////////////////////////////////

		MathVector<dim> Dn;
		VecSet(Dn, 0.0);

		for(size_t i = 0; i < dim; ++i)
		{
			for(size_t j = 0; j < dim; ++j) Dn[i] += D(i,j)*n[j];
		}

		////////////////////////////////////////////////////////////////////////////
		// Tangential projection
		//
		// P D n = D n - n (n . D n)
		////////////////////////////////////////////////////////////////////////////

		number nDn = 0.0;
		for(size_t d = 0; d < dim; ++d) nDn += n[d]*Dn[d];

		////////////////////////////////////////////////////////////////////////////
		// Aerodynamic tangential traction
		//
		// tau_Gamma = 2 mu_a P D n
		////////////////////////////////////////////////////////////////////////////

		MathVector<dim> tau;

		for(size_t d = 0; d < dim; ++d) tau[d] = 2.0*mu_a*(Dn[d] - n[d]*nDn);

		////////////////////////////////////////////////////////////////////////////
		// Magnitude of tangential traction
		////////////////////////////////////////////////////////////////////////////

		number tau2 = 0.0;
		for(size_t d = 0; d < dim; ++d) tau2 += tau[d]*tau[d];

		const number tauMag = std::sqrt(tau2);

		if(tauMag <= 0.0) return;

		////////////////////////////////////////////////////////////////////////////
		// Friction velocity
		//
		// u_* = sqrt(|tau_Gamma|/rho_a)
		////////////////////////////////////////////////////////////////////////////

		const number uStar = std::sqrt(tauMag/rho_a);

		////////////////////////////////////////////////////////////////////////////
		// Threshold condition
		////////////////////////////////////////////////////////////////////////////

		if(uStar <= m_uStarThreshold) return;

		const number excess = uStar - m_uStarThreshold;

		////////////////////////////////////////////////////////////////////////////
		// Saturated saltation surface flux
		//
		// q_s = (C/g) (u_* - u_*t) tau_Gamma
		//
		// Diffuse volumetric flux
		//
		// J_salt = delta_Gamma q_s
		////////////////////////////////////////////////////////////////////////////

		const number factor = deltaGamma*(m_C/grav)*excess;

		for(size_t d = 0; d < dim; ++d) flux[d] = factor*tau[d];
	}


	////////////////////////////////////////////////////////////////////////////
	// Analytical directional derivative
	////////////////////////////////////////////////////////////////////////////

	void compute_flux_derivative(MathVector<dim>& dFlux, const MathVector<dim>& gradC, const MathMatrix<dim,dim>& gradU, const MathVector<dim>& dGradC, const MathMatrix<dim,dim>& dGradU) const
	{
		VecSet(dFlux, 0.0);

		if(!Inter) UG_THROW("SaltationFluxLinker: Interface pointer is null.");

		const number rho_a = Inter->Density_a();
		const number mu_a = Inter->Viscosity_a();
		const number grav = std::fabs(Inter->gravity());

		if(rho_a <= 0.0) UG_THROW("SaltationFluxLinker: air density must be positive.");
		if(mu_a < 0.0) UG_THROW("SaltationFluxLinker: air viscosity must be non-negative.");
		if(grav <= 0.0) UG_THROW("SaltationFluxLinker: gravity magnitude must be positive.");

		////////////////////////////////////////////////////////////////////////////
		// grad(c)
		////////////////////////////////////////////////////////////////////////////

		number gradC2 = 0.0;
		for(size_t d = 0; d < dim; ++d) gradC2 += gradC[d]*gradC[d];

		const number sN = std::sqrt(gradC2 + m_epsNormal*m_epsNormal);
		const number sDelta = std::sqrt(gradC2 + m_epsDelta*m_epsDelta);
		const number deltaGamma = sDelta - m_epsDelta;

		if(deltaGamma <= 0.0) return;

		////////////////////////////////////////////////////////////////////////////
		// Interface normal
		////////////////////////////////////////////////////////////////////////////

		MathVector<dim> n;
		for(size_t d = 0; d < dim; ++d) n[d] = -gradC[d]/sN;

		////////////////////////////////////////////////////////////////////////////
		// Derivative of interface normal
		//
		// dn = -[I/sN - grad(c) tensor grad(c)/sN^3] dgrad(c)
		////////////////////////////////////////////////////////////////////////////

		number gradCDotDGradC = 0.0;
		for(size_t d = 0; d < dim; ++d) gradCDotDGradC += gradC[d]*dGradC[d];

		const number sN3 = sN*sN*sN;

		MathVector<dim> dn;
		for(size_t d = 0; d < dim; ++d) dn[d] = -dGradC[d]/sN + gradC[d]*gradCDotDGradC/sN3;

		////////////////////////////////////////////////////////////////////////////
		// Derivative of diffuse-interface delta
		////////////////////////////////////////////////////////////////////////////

		const number dDeltaGamma = gradCDotDGradC/sDelta;

		////////////////////////////////////////////////////////////////////////////
		// D and dD
		////////////////////////////////////////////////////////////////////////////

		MathMatrix<dim,dim> D;
		MathMatrix<dim,dim> dD;

		for(size_t i = 0; i < dim; ++i)
		{
			for(size_t j = 0; j < dim; ++j)
			{
				D(i,j) = 0.5*(gradU(i,j) + gradU(j,i));
				dD(i,j) = 0.5*(dGradU(i,j) + dGradU(j,i));
			}
		}

		////////////////////////////////////////////////////////////////////////////
		// Dn and derivative
		//
		// d(Dn) = dD n + D dn
		////////////////////////////////////////////////////////////////////////////

		MathVector<dim> Dn;
		MathVector<dim> dDn;

		VecSet(Dn, 0.0);
		VecSet(dDn, 0.0);

		for(size_t i = 0; i < dim; ++i)
		{
			for(size_t j = 0; j < dim; ++j)
			{
				Dn[i] += D(i,j)*n[j];
				dDn[i] += dD(i,j)*n[j] + D(i,j)*dn[j];
			}
		}

		////////////////////////////////////////////////////////////////////////////
		// n . Dn and derivative
		////////////////////////////////////////////////////////////////////////////

		number nDn = 0.0;
		number dNDn = 0.0;

		for(size_t d = 0; d < dim; ++d)
		{
			nDn += n[d]*Dn[d];
			dNDn += dn[d]*Dn[d] + n[d]*dDn[d];
		}

		////////////////////////////////////////////////////////////////////////////
		// Tangential traction and derivative
		//
		// tau  = 2 mu_a (Dn - n(n.Dn))
		//
		// dtau = 2 mu_a [dDn - dn(n.Dn) - n d(n.Dn)]
		////////////////////////////////////////////////////////////////////////////

		MathVector<dim> tau;
		MathVector<dim> dTau;

		for(size_t d = 0; d < dim; ++d)
		{
			tau[d] = 2.0*mu_a*(Dn[d] - n[d]*nDn);
			dTau[d] = 2.0*mu_a*(dDn[d] - dn[d]*nDn - n[d]*dNDn);
		}

		////////////////////////////////////////////////////////////////////////////
		// Magnitude of traction
		////////////////////////////////////////////////////////////////////////////

		number tau2 = 0.0;
		for(size_t d = 0; d < dim; ++d) tau2 += tau[d]*tau[d];

		const number tauMag = std::sqrt(tau2);

		if(tauMag <= 0.0) return;

		////////////////////////////////////////////////////////////////////////////
		// Friction velocity
		////////////////////////////////////////////////////////////////////////////

		const number uStar = std::sqrt(tauMag/rho_a);

		if(uStar <= m_uStarThreshold) return;

		////////////////////////////////////////////////////////////////////////////
		// Derivative of friction velocity
		//
		// d|tau| = tau . dtau / |tau|
		//
		// du_* = (tau . dtau)/(2 rho_a u_* |tau|)
		////////////////////////////////////////////////////////////////////////////

		number tauDotDTau = 0.0;
		for(size_t d = 0; d < dim; ++d) tauDotDTau += tau[d]*dTau[d];

		const number dUStar = tauDotDTau/(2.0*rho_a*uStar*tauMag);

		////////////////////////////////////////////////////////////////////////////
		// Derivative of saltation flux
		//
		// J = delta (C/g) (u_* - u_*t) tau
		//
		// dJ = (C/g) [
		//          ddelta (u_* - u_*t) tau
		//        + delta du_* tau
		//        + delta (u_* - u_*t) dtau
		//      ]
		////////////////////////////////////////////////////////////////////////////

		const number excess = uStar - m_uStarThreshold;
		const number Cg = m_C/grav;

		for(size_t d = 0; d < dim; ++d) dFlux[d] = Cg*(dDeltaGamma*excess*tau[d] + deltaGamma*dUStar*tau[d] + deltaGamma*excess*dTau[d]);
	}


public:

	////////////////////////////////////////////////////////////////////////////
	// Single integration-point evaluation
	////////////////////////////////////////////////////////////////////////////

	inline void evaluate(MathVector<dim>& value, const MathVector<dim>& globIP, number time, int si) const
	{
		MathVector<dim> gradC;
		MathMatrix<dim,dim> gradU;

		(*m_spVolumeGrad)(gradC, globIP, time, si);
		(*m_spVelocityGrad)(gradU, globIP, time, si);

		compute_flux(value, gradC, gradU);
	}


	////////////////////////////////////////////////////////////////////////////
	// Multiple integration-point evaluation
	////////////////////////////////////////////////////////////////////////////

	template <int refDim>
	inline void evaluate(MathVector<dim> vValue[], const MathVector<dim> vGlobIP[], number time, int si, GridObject* elem, const MathVector<dim> vCornerCoords[], const MathVector<refDim> vLocIP[], const size_t nip, LocalVector* u, const MathMatrix<refDim,dim>* vJT = NULL) const
	{
		std::vector<MathVector<dim> > vVolumeGrad(nip);
		std::vector<MathMatrix<dim,dim> > vVelocityGrad(nip);

		(*m_spVolumeGrad)(&vVolumeGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
		(*m_spVelocityGrad)(&vVelocityGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);

		for(size_t ip = 0; ip < nip; ++ip) compute_flux(vValue[ip], vVolumeGrad[ip], vVelocityGrad[ip]);
	}


	////////////////////////////////////////////////////////////////////////////
	// Evaluation and analytical derivatives
	////////////////////////////////////////////////////////////////////////////

	template <int refDim>
	void eval_and_deriv(MathVector<dim> vValue[], const MathVector<dim> vGlobIP[], number time, int si, GridObject* elem, const MathVector<dim> vCornerCoords[], const MathVector<refDim> vLocIP[], const size_t nip, LocalVector* u, bool bDeriv, int s, std::vector<std::vector<MathVector<dim> > > vvvDeriv[], const MathMatrix<refDim,dim>* vJT = NULL) const
	{
		const int s_DC_ = base_type::series_id(_DC_, s);
		const int s_DU_ = base_type::series_id(_DU_, s);

		const MathVector<dim>* vVolumeGrad = m_spVolumeGrad->values(s_DC_);
		const MathMatrix<dim,dim>* vVelocityGrad = m_spVelocityGrad->values(s_DU_);

		for(size_t ip = 0; ip < nip; ++ip) compute_flux(vValue[ip], vVolumeGrad[ip], vVelocityGrad[ip]);

		if(!bDeriv || this->zero_derivative()) return;

		this->set_zero(vvvDeriv, nip);

		////////////////////////////////////////////////////////////////////////////
		// Derivatives with respect to grad(c)
		////////////////////////////////////////////////////////////////////////////

		if(m_spDVolumeGrad.valid() && !m_spDVolumeGrad->zero_derivative())
		{
			for(size_t ip = 0; ip < nip; ++ip)
			{
				for(size_t fct = 0; fct < m_spDVolumeGrad->num_fct(); ++fct)
				{
					const MathVector<dim>* vDVolumeGrad = m_spDVolumeGrad->deriv(s_DC_, ip, fct);
					const size_t commonFct = this->input_common_fct(_DC_, fct);

					for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
					{
						MathMatrix<dim,dim> zeroGradU;
						MathVector<dim> dFlux;

						MatSet(zeroGradU, 0.0);

						compute_flux_derivative(dFlux, vVolumeGrad[ip], vVelocityGrad[ip], vDVolumeGrad[sh], zeroGradU);

						vvvDeriv[ip][commonFct][sh] += dFlux;
					}
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////
		// Derivatives with respect to grad(u)
		////////////////////////////////////////////////////////////////////////////

		if(m_spDVelocityGrad.valid() && !m_spDVelocityGrad->zero_derivative())
		{
			for(size_t ip = 0; ip < nip; ++ip)
			{
				for(size_t fct = 0; fct < m_spDVelocityGrad->num_fct(); ++fct)
				{
					const MathMatrix<dim,dim>* vDVelocityGrad = m_spDVelocityGrad->deriv(s_DU_, ip, fct);
					const size_t commonFct = this->input_common_fct(_DU_, fct);

					for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
					{
						MathVector<dim> zeroGradC;
						MathVector<dim> dFlux;

						VecSet(zeroGradC, 0.0);

						compute_flux_derivative(dFlux, vVolumeGrad[ip], vVelocityGrad[ip], zeroGradC, vDVelocityGrad[sh]);

						vvvDeriv[ip][commonFct][sh] += dFlux;
					}
				}
			}
		}
	}


public:

	////////////////////////////////////////////////////////////////////////////
	// Input setters
	////////////////////////////////////////////////////////////////////////////

	void set_volume_grad(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
	{
		m_spVolumeGrad = data;
		m_spDVolumeGrad = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim> >();
		base_type::set_input(_DC_, data, data);
	}


	void set_velocity_gradient(SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > data)
	{
		m_spVelocityGrad = data;
		m_spDVelocityGrad = data.template cast_dynamic<DependentUserData<MathMatrix<dim,dim>, dim> >();
		base_type::set_input(_DU_, data, data);
	}


	////////////////////////////////////////////////////////////////////////////
	// Physical parameters
	////////////////////////////////////////////////////////////////////////////

	void set_phase_parameters(Interface<dim>* user)
	{
		if(!user) UG_THROW("SaltationFluxLinker: Interface pointer is null.");
		if(!user->valid()) UG_THROW("SaltationFluxLinker: Interface parameters have not been initialized.");

		Inter = user;
	}


	////////////////////////////////////////////////////////////////////////////
	// Saltation parameters
	////////////////////////////////////////////////////////////////////////////

	void set_saltation_coefficient(number C) {m_C = C;}

	void set_threshold_friction_velocity(number uStarThreshold) {m_uStarThreshold = uStarThreshold;}

	void set_normal_epsilon(number eps) {m_epsNormal = eps;}

	void set_delta_epsilon(number eps) {m_epsDelta = eps;}


protected:

	////////////////////////////////////////////////////////////////////////////
	// Input indices
	////////////////////////////////////////////////////////////////////////////

	static const size_t _DC_ = 0;
	static const size_t _DU_ = 1;


	////////////////////////////////////////////////////////////////////////////
	// grad(c)
	////////////////////////////////////////////////////////////////////////////

	SmartPtr<CplUserData<MathVector<dim>, dim> > m_spVolumeGrad;
	SmartPtr<DependentUserData<MathVector<dim>, dim> > m_spDVolumeGrad;


	////////////////////////////////////////////////////////////////////////////
	// grad(u)
	////////////////////////////////////////////////////////////////////////////

	SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > m_spVelocityGrad;
	SmartPtr<DependentUserData<MathMatrix<dim,dim>, dim> > m_spDVelocityGrad;


	////////////////////////////////////////////////////////////////////////////
	// Shared phase parameters
	////////////////////////////////////////////////////////////////////////////

	Interface<dim>* Inter;


	////////////////////////////////////////////////////////////////////////////
	// Saltation parameters
	////////////////////////////////////////////////////////////////////////////

	number m_C;
	number m_uStarThreshold;
	number m_epsNormal;
	number m_epsDelta;
};


} // namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__SALTATION_FLUX_LINKER__ */
