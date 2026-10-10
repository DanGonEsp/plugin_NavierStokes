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

	SaltationFluxLinker() : m_spVolumeGrad(NULL), m_spDVolumeGrad(NULL), m_spVelocityGrad(NULL), m_spDVelocityGrad(NULL), m_spTurbulentKinViscosity(NULL), m_spDTurbulentKinViscosity(NULL), Inter(NULL), m_C(5.5), m_uStarThreshold(0.22), m_epsNormal(1e-8), m_epsDelta(1e-8)
	{
		this->set_num_input(3);
	}
	struct SaltationData
	{
		number sN;
		number sDelta;
		number deltaGamma;
		number normal2;
		number nDn;
		number tauMag;
		number uStar;
		MathVector<dim> normal;
		MathVector<dim> Dn;
		MathVector<dim> tau;
		MathMatrix<dim,dim> D;
	};


private:

	////////////////////////////////////////////////////////////////////////////
	// Flux evaluation
	////////////////////////////////////////////////////////////////////////////

	void compute_flux(MathVector<dim>& flux, const MathVector<dim>& gradC, const MathMatrix<dim,dim>& gradU, number nuT) const
	{
		VecSet(flux, 0.0);

		if(!Inter) UG_THROW("SaltationFluxLinker: Interface pointer is null.");

		const number grav = std::fabs(Inter->gravity());

		if(grav <= 0.0) UG_THROW("SaltationFluxLinker: gravity magnitude must be positive.");

		SaltationData data;
		compute_saltation_data(data, gradC, gradU, nuT);

		if(data.deltaGamma <= 0.0) return;
		if(data.uStar <= m_uStarThreshold) return;

		const number excess = data.uStar - m_uStarThreshold;
		const number factor = data.deltaGamma*(m_C/grav)*excess;

		for(size_t d = 0; d < dim; ++d)
			flux[d] = factor*data.tau[d];
	}


	////////////////////////////////////////////////////////////////////////////
	// Analytical directional derivative
	////////////////////////////////////////////////////////////////////////////

	void compute_flux_derivative(MathVector<dim>& dFlux, const MathVector<dim>& gradC, const MathMatrix<dim,dim>& gradU, number nuT, const MathVector<dim>& dGradC, const MathMatrix<dim,dim>& dGradU, number dNuT = 0.0) const
	{
		VecSet(dFlux, 0.0);

		if(!Inter)
			UG_THROW("SaltationFluxLinker: Interface pointer is null.");

		const number rho_a = Inter->Density_a();
		const number mu_a = Inter->Viscosity_a();
		const number mu_t = rho_a*nuT;
		const number mu_eff = mu_a + mu_t;
		const number grav = std::fabs(Inter->gravity());

		if(rho_a <= 0.0)
			UG_THROW("SaltationFluxLinker: air density must be positive.");

		if(mu_a < 0.0)
			UG_THROW("SaltationFluxLinker: air viscosity must be non-negative.");

		if(grav <= 0.0)
			UG_THROW("SaltationFluxLinker: gravity magnitude must be positive.");

		SaltationData data;
		compute_saltation_data(data, gradC, gradU, nuT);

		if(data.deltaGamma <= 0.0)
			return;

		if(data.tauMag <= 0.0)
			return;

		if(data.uStar <= m_uStarThreshold)
			return;

		number gradCDotDGradC = 0.0;

		for(size_t d = 0; d < dim; ++d)
			gradCDotDGradC += gradC[d]*dGradC[d];

		const number sN3 = data.sN*data.sN*data.sN;

		MathVector<dim> dn;

		for(size_t d = 0; d < dim; ++d)
			dn[d] = -dGradC[d]/data.sN + gradC[d]*gradCDotDGradC/sN3;
		
		number dNormal2 = 0.0;

		for(size_t d = 0; d < dim; ++d)
			dNormal2 += 2.0*data.normal[d]*dn[d];

		const number dDeltaGamma = gradCDotDGradC/data.sDelta;

		MathMatrix<dim,dim> dD;

		for(size_t i = 0; i < dim; ++i)
		{
			for(size_t j = 0; j < dim; ++j)
				dD(i,j) = 0.5*(dGradU(i,j) + dGradU(j,i));
		}

		MathVector<dim> dDn;
		VecSet(dDn, 0.0);

		for(size_t i = 0; i < dim; ++i)
		{
			for(size_t j = 0; j < dim; ++j)
				dDn[i] += dD(i,j)*data.normal[j] + data.D(i,j)*dn[j];
		}

		number dNDn = 0.0;

		for(size_t d = 0; d < dim; ++d)
			dNDn += dn[d]*data.Dn[d] + data.normal[d]*dDn[d];

		MathVector<dim> dTau;
		const number dMuEff = rho_a*dNuT;

		for(size_t d = 0; d < dim; ++d)
		{
			const number tangent = data.normal2*data.Dn[d] - data.normal[d]*data.nDn;
			const number dTangent = dNormal2*data.Dn[d] + data.normal2*dDn[d] - dn[d]*data.nDn - data.normal[d]*dNDn;
			dTau[d] = 2.0*dMuEff*tangent + 2.0*mu_eff*dTangent;
		}

		number tauDotDTau = 0.0;

		for(size_t d = 0; d < dim; ++d)
			tauDotDTau += data.tau[d]*dTau[d];

		const number dUStar = tauDotDTau/(2.0*rho_a*data.uStar*data.tauMag);
		const number excess = data.uStar - m_uStarThreshold;
		const number Cg = m_C/grav;

		for(size_t d = 0; d < dim; ++d)
			dFlux[d] = Cg*(dDeltaGamma*excess*data.tau[d] + data.deltaGamma*dUStar*data.tau[d] + data.deltaGamma*excess*dTau[d]);
	}
	
	void compute_saltation_data(SaltationData& data, const MathVector<dim>& gradC, const MathMatrix<dim,dim>& gradU, number nuT) const
	{
		data.sN = 0.0;
		data.sDelta = 0.0;
		data.deltaGamma = 0.0;
		data.normal2 = 0.0;
		data.nDn = 0.0;
		data.tauMag = 0.0;
		data.uStar = 0.0;

		VecSet(data.normal, 0.0);
		VecSet(data.Dn, 0.0);
		VecSet(data.tau, 0.0);

		if(!Inter)
			UG_THROW("SaltationFluxLinker: Interface pointer is null.");

		const number rho_a = Inter->Density_a();
		const number mu_a = Inter->Viscosity_a();
		const number mu_t = rho_a*nuT;
		const number mu_eff = mu_a + mu_t;

		if(rho_a <= 0.0)
			UG_THROW("SaltationFluxLinker: air density must be positive.");

		if(mu_a < 0.0)
			UG_THROW("SaltationFluxLinker: air viscosity must be non-negative.");
		if(nuT < 0.0)
			UG_THROW("SaltationFluxLinker: turbulent kinematic viscosity must be non-negative.");

		number gradC2 = 0.0;

		for(size_t d = 0; d < dim; ++d)
			gradC2 += gradC[d]*gradC[d];

		data.sN = std::sqrt(gradC2 + m_epsNormal*m_epsNormal);
		data.sDelta = std::sqrt(gradC2 + m_epsDelta*m_epsDelta);
		data.deltaGamma = data.sDelta - m_epsDelta;

		for(size_t d = 0; d < dim; ++d)
			data.normal[d] = -gradC[d]/data.sN;

		for(size_t d = 0; d < dim; ++d)
			data.normal2 += data.normal[d]*data.normal[d];

		for(size_t i = 0; i < dim; ++i)
		{
			for(size_t j = 0; j < dim; ++j)
				data.D(i,j) = 0.5*(gradU(i,j) + gradU(j,i));
		}

		for(size_t i = 0; i < dim; ++i)
		{
			for(size_t j = 0; j < dim; ++j)
				data.Dn[i] += data.D(i,j)*data.normal[j];
		}

		for(size_t d = 0; d < dim; ++d)
			data.nDn += data.normal[d]*data.Dn[d];

		for(size_t d = 0; d < dim; ++d)
			data.tau[d] = 2.0*mu_eff*(data.normal2*data.Dn[d] - data.normal[d]*data.nDn);

		number tau2 = 0.0;

		for(size_t d = 0; d < dim; ++d)
			tau2 += data.tau[d]*data.tau[d];

		data.tauMag = std::sqrt(tau2);

		if(data.tauMag > 0.0)
			data.uStar = std::sqrt(data.tauMag/rho_a);
	}


public:

	////////////////////////////////////////////////////////////////////////////
	// Single integration-point evaluation
	////////////////////////////////////////////////////////////////////////////

	inline void evaluate(MathVector<dim>& value, const MathVector<dim>& globIP, number time, int si) const
	{
		MathVector<dim> gradC;
		MathMatrix<dim,dim> gradU;
		number nuT;
		
		(*m_spTurbulentKinViscosity)(nuT, globIP, time, si);
		(*m_spVolumeGrad)(gradC, globIP, time, si);
		(*m_spVelocityGrad)(gradU, globIP, time, si);

		compute_flux(value, gradC, gradU, nuT);
	}


	////////////////////////////////////////////////////////////////////////////
	// Multiple integration-point evaluation
	////////////////////////////////////////////////////////////////////////////

	template <int refDim>
	inline void evaluate(MathVector<dim> vValue[], const MathVector<dim> vGlobIP[], number time, int si, GridObject* elem, const MathVector<dim> vCornerCoords[], const MathVector<refDim> vLocIP[], const size_t nip, LocalVector* u, const MathMatrix<refDim,dim>* vJT = NULL) const
	{
		std::vector<MathVector<dim> > vVolumeGrad(nip);
		std::vector<MathMatrix<dim,dim> > vVelocityGrad(nip);
		std::vector<number> vNuT(nip);
		
		(*m_spTurbulentKinViscosity)(&vNuT[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
		(*m_spVolumeGrad)(&vVolumeGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
		(*m_spVelocityGrad)(&vVelocityGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);

		for(size_t ip = 0; ip < nip; ++ip)
			compute_flux(vValue[ip], vVolumeGrad[ip], vVelocityGrad[ip], vNuT[ip]);
	}


	////////////////////////////////////////////////////////////////////////////
	// Evaluation and analytical derivatives
	////////////////////////////////////////////////////////////////////////////

	template <int refDim>
	void eval_and_deriv(MathVector<dim> vValue[], const MathVector<dim> vGlobIP[], number time, int si, GridObject* elem, const MathVector<dim> vCornerCoords[], const MathVector<refDim> vLocIP[], const size_t nip, LocalVector* u, bool bDeriv, int s, std::vector<std::vector<MathVector<dim> > > vvvDeriv[], const MathMatrix<refDim,dim>* vJT = NULL) const
	{
		const int s_DC_ = base_type::series_id(_DC_, s);
		const int s_DU_ = base_type::series_id(_DU_, s);
		const int s_DNUT_ = base_type::series_id(_DNUT_, s);

		const MathVector<dim>* vVolumeGrad = m_spVolumeGrad->values(s_DC_);
		const MathMatrix<dim,dim>* vVelocityGrad = m_spVelocityGrad->values(s_DU_);
		const number* vNuT = m_spTurbulentKinViscosity->values(s_DNUT_);

		for(size_t ip = 0; ip < nip; ++ip)
			compute_flux(vValue[ip], vVolumeGrad[ip], vVelocityGrad[ip], vNuT[ip]);

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

						compute_flux_derivative(dFlux, vVolumeGrad[ip], vVelocityGrad[ip], vNuT[ip], vDVolumeGrad[sh], zeroGradU);

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

						compute_flux_derivative(dFlux, vVolumeGrad[ip], vVelocityGrad[ip], vNuT[ip], zeroGradC, vDVelocityGrad[sh]);

						vvvDeriv[ip][commonFct][sh] += dFlux;
					}
				}
			}
		}
		// Derivatives through the turbulent kinematic viscosity input.
		if(m_spDTurbulentKinViscosity.valid() && !m_spDTurbulentKinViscosity->zero_derivative())
		{
			MathVector<dim> zeroGradC;
			MathMatrix<dim,dim> zeroGradU;
			VecSet(zeroGradC, 0.0);
			MatSet(zeroGradU, 0.0);

			for(size_t ip = 0; ip < nip; ++ip)
			{
				for(size_t fct = 0; fct < m_spDTurbulentKinViscosity->num_fct(); ++fct)
				{
					const number* vDNuT = m_spDTurbulentKinViscosity->deriv(s_DNUT_, ip, fct);
					const size_t commonFct = this->input_common_fct(_DNUT_, fct);

					for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
					{
						MathVector<dim> dFlux;
						compute_flux_derivative(dFlux, vVolumeGrad[ip], vVelocityGrad[ip], vNuT[ip], zeroGradC, zeroGradU, vDNuT[sh]);
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
	
	void set_turbulent_kinematic_viscosity(SmartPtr<CplUserData<number,dim> > data)
	{
		m_spTurbulentKinViscosity = data;
		m_spDTurbulentKinViscosity = data.template cast_dynamic<DependentUserData<number,dim> >();
		base_type::set_input(_DNUT_, data, data);
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
	static const size_t _DNUT_ = 2;


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
	// Mu_t(u)
	////////////////////////////////////////////////////////////////////////////
	
	SmartPtr<CplUserData<number,dim> > m_spTurbulentKinViscosity;
	SmartPtr<DependentUserData<number,dim> > m_spDTurbulentKinViscosity;


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
