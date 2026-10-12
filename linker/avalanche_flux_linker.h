/*
 * Copyright (c) 2013-2015: G-CSC, Goethe University Frankfurt
 *
 * Avalanche flux linker for the Navier-Stokes plugin.
 *
 * The linker embeds the Degenerate Parabolic Sliding Model (DPSM)
 * into the diffuse-interface volume-fraction formulation.
 *
 * For a bed surface z = h(x), the DPSM reads
 *
 *   dh/dt = nu_av div_h[
 *              ((|grad_h h| - tan(theta_cr))_+) / (1 + |grad_h h|^2)
 *              * grad_h h / |grad_h h|
 *           ].
 *
 * In conservative form the downhill surface flux is
 *
 *   q_av = -nu_av
 *          ((|grad_h h| - tan(theta_cr))_+) / (1 + |grad_h h|^2)
 *          * grad_h h / |grad_h h|.
 *
 * With the upward diffuse-interface normal
 *
 *   n = -grad(c) / sqrt(|grad(c)|^2 + eps_n^2),
 *
 * and n_h = n - (n.e_z)e_z,
 *
 *   |grad_h h| = |n_h| / |n_z|,
 *   -grad_h h / |grad_h h| = n_h / |n_h|.
 *
 * Hence
 *
 *   q_av = nu_av
 *          ((s - tan(theta_cr))_+) / (1 + s^2)
 *          * n_h / |n_h|,
 *
 * where s = |n_h| / |n_z|.
 *
 * As in SaltationFluxLinker, the surface flux is localized through
 *
 *   delta_Gamma = sqrt(|grad(c)|^2 + eps_delta^2) - eps_delta,
 *
 * so that
 *
 *   J_av = delta_Gamma q_av.
 *
 * The output is the diffuse avalanche flux entering the dimensionless
 * concentration equation:
 *
 *   J_av = delta_Gamma q_av.
 *
 * With nu_av in [m^2/s] and delta_Gamma in [1/m], the output has units [m/s].
 *
 * Analytical derivatives are provided with respect to grad(c).
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__AVALANCHE_FLUX_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__AVALANCHE_FLUX_LINKER__

#include <cmath>
#include <vector>

#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "../properties_interface.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Avalanche Flux Linker
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class AvalancheFluxLinker : public StdDataLinker<AvalancheFluxLinker<dim>, MathVector<dim>, dim>
{
	typedef StdDataLinker<AvalancheFluxLinker<dim>, MathVector<dim>, dim> base_type;

public:

	AvalancheFluxLinker() : m_spVolumeGrad(NULL), m_spDVolumeGrad(NULL), Inter(NULL), m_nuAv(0.0), m_thetaCr(34.0*std::acos(-1.0)/180.0), m_criticalSlope(std::tan(34.0*std::acos(-1.0)/180.0)), m_epsNormal(1e-8), m_epsDelta(1e-8), m_epsVertical(1e-8)
	{
		this->set_num_input(1);
	}

	struct AvalancheData
	{
		number sN;
		number sDelta;
		number deltaGamma;
		number normal2;
		number horizontalNormalMag;
		number verticalNormalMag;
		number slope;
		number excessSlope;
		number mobility;
		MathVector<dim> normal;
		MathVector<dim> horizontalNormal;
		MathVector<dim> downhillDirection;
	};


private:

	////////////////////////////////////////////////////////////////////////////
	// Flux evaluation
	////////////////////////////////////////////////////////////////////////////

	void compute_flux(MathVector<dim>& flux, const MathVector<dim>& gradC) const
	{
		VecSet(flux, 0.0);

		if(!Inter) UG_THROW("AvalancheFluxLinker: Interface pointer is null.");

		AvalancheData data;
		compute_avalanche_data(data, gradC);

		if(data.deltaGamma <= 0.0) return;
		if(data.excessSlope <= 0.0) return;
		if(data.horizontalNormalMag <= 0.0) return;

		const number factor = m_nuAv*data.deltaGamma*data.mobility;

		for(size_t d = 0; d < dim; ++d)
			flux[d] = factor*data.downhillDirection[d];
	}


	////////////////////////////////////////////////////////////////////////////
	// Analytical directional derivative
	////////////////////////////////////////////////////////////////////////////

	void compute_flux_derivative(MathVector<dim>& dFlux, const MathVector<dim>& gradC, const MathVector<dim>& dGradC) const
	{
		VecSet(dFlux, 0.0);

		if(!Inter) UG_THROW("AvalancheFluxLinker: Interface pointer is null.");

		AvalancheData data;
		compute_avalanche_data(data, gradC);

		if(data.deltaGamma <= 0.0) return;
		if(data.excessSlope <= 0.0) return;
		if(data.horizontalNormalMag <= 0.0) return;

		number gradCDotDGradC = 0.0;

		for(size_t d = 0; d < dim; ++d)
			gradCDotDGradC += gradC[d]*dGradC[d];

		const number sN3 = data.sN*data.sN*data.sN;

		MathVector<dim> dNormal;

		for(size_t d = 0; d < dim; ++d)
			dNormal[d] = -dGradC[d]/data.sN + gradC[d]*gradCDotDGradC/sN3;

		const number dDeltaGamma = gradCDotDGradC/data.sDelta;

		MathVector<dim> dHorizontalNormal;
		VecSet(dHorizontalNormal, 0.0);

		for(size_t d = 0; d < dim-1; ++d)
			dHorizontalNormal[d] = dNormal[d];

		number horizontalDotDerivative = 0.0;

		for(size_t d = 0; d < dim-1; ++d)
			horizontalDotDerivative += data.horizontalNormal[d]*dHorizontalNormal[d];

		const number dHorizontalNormalMag = horizontalDotDerivative/data.horizontalNormalMag;

		const number nz = data.normal[dim-1];
		const number dNz = dNormal[dim-1];
		const number dVerticalNormalMag = nz*dNz/data.verticalNormalMag;

		const number dSlope = dHorizontalNormalMag/data.verticalNormalMag - data.horizontalNormalMag*dVerticalNormalMag/(data.verticalNormalMag*data.verticalNormalMag);

		const number slope2 = data.slope*data.slope;
		const number denom = 1.0 + slope2;
		const number dMobility_dSlope = (denom - 2.0*data.slope*data.excessSlope)/(denom*denom);
		const number dMobility = dMobility_dSlope*dSlope;

		MathVector<dim> dDownhillDirection;
		VecSet(dDownhillDirection, 0.0);

		const number horizontalMag2 = data.horizontalNormalMag*data.horizontalNormalMag;

		for(size_t d = 0; d < dim-1; ++d)
			dDownhillDirection[d] = dHorizontalNormal[d]/data.horizontalNormalMag - data.horizontalNormal[d]*dHorizontalNormalMag/horizontalMag2;

		for(size_t d = 0; d < dim; ++d)
		{
			dFlux[d] = m_nuAv*(dDeltaGamma*data.mobility*data.downhillDirection[d] + data.deltaGamma*dMobility*data.downhillDirection[d] + data.deltaGamma*data.mobility*dDownhillDirection[d]);
		}
	}


	////////////////////////////////////////////////////////////////////////////
	// DPSM data
	////////////////////////////////////////////////////////////////////////////

	void compute_avalanche_data(AvalancheData& data, const MathVector<dim>& gradC) const
	{
		data.sN = 0.0;
		data.sDelta = 0.0;
		data.deltaGamma = 0.0;
		data.normal2 = 0.0;
		data.horizontalNormalMag = 0.0;
		data.verticalNormalMag = 0.0;
		data.slope = 0.0;
		data.excessSlope = 0.0;
		data.mobility = 0.0;

		VecSet(data.normal, 0.0);
		VecSet(data.horizontalNormal, 0.0);
		VecSet(data.downhillDirection, 0.0);

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

		number horizontal2 = 0.0;

		for(size_t d = 0; d < dim-1; ++d)
		{
			data.horizontalNormal[d] = data.normal[d];
			horizontal2 += data.horizontalNormal[d]*data.horizontalNormal[d];
		}

		data.horizontalNormal[dim-1] = 0.0;
		data.horizontalNormalMag = std::sqrt(horizontal2);

		const number nz = data.normal[dim-1];
		data.verticalNormalMag = std::sqrt(nz*nz + m_epsVertical*m_epsVertical);

		if(data.horizontalNormalMag <= 0.0)
			return;

		data.slope = data.horizontalNormalMag/data.verticalNormalMag;
		data.excessSlope = data.slope - m_criticalSlope;

		if(data.excessSlope <= 0.0)
		{
			data.excessSlope = 0.0;
			return;
		}

		data.mobility = data.excessSlope/(1.0 + data.slope*data.slope);

		for(size_t d = 0; d < dim-1; ++d)
			data.downhillDirection[d] = data.horizontalNormal[d]/data.horizontalNormalMag;

		data.downhillDirection[dim-1] = 0.0;
	}


public:

	////////////////////////////////////////////////////////////////////////////
	// Single integration-point evaluation
	////////////////////////////////////////////////////////////////////////////

	inline void evaluate(MathVector<dim>& value, const MathVector<dim>& globIP, number time, int si) const
	{
		MathVector<dim> gradC;
		(*m_spVolumeGrad)(gradC, globIP, time, si);
		compute_flux(value, gradC);
	}


	////////////////////////////////////////////////////////////////////////////
	// Multiple integration-point evaluation
	////////////////////////////////////////////////////////////////////////////

	template <int refDim>
	inline void evaluate(MathVector<dim> vValue[], const MathVector<dim> vGlobIP[], number time, int si, GridObject* elem, const MathVector<dim> vCornerCoords[], const MathVector<refDim> vLocIP[], const size_t nip, LocalVector* u, const MathMatrix<refDim,dim>* vJT = NULL) const
	{
		std::vector<MathVector<dim> > vVolumeGrad(nip);

		(*m_spVolumeGrad)(&vVolumeGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);

		for(size_t ip = 0; ip < nip; ++ip)
			compute_flux(vValue[ip], vVolumeGrad[ip]);
	}


	////////////////////////////////////////////////////////////////////////////
	// Evaluation and analytical derivatives
	////////////////////////////////////////////////////////////////////////////

	template <int refDim>
	void eval_and_deriv(MathVector<dim> vValue[], const MathVector<dim> vGlobIP[], number time, int si, GridObject* elem, const MathVector<dim> vCornerCoords[], const MathVector<refDim> vLocIP[], const size_t nip, LocalVector* u, bool bDeriv, int s, std::vector<std::vector<MathVector<dim> > > vvvDeriv[], const MathMatrix<refDim,dim>* vJT = NULL) const
	{
		const int s_DC_ = base_type::series_id(_DC_, s);
		const MathVector<dim>* vVolumeGrad = m_spVolumeGrad->values(s_DC_);

		for(size_t ip = 0; ip < nip; ++ip)
			compute_flux(vValue[ip], vVolumeGrad[ip]);

		if(!bDeriv || this->zero_derivative()) return;

		this->set_zero(vvvDeriv, nip);

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
						MathVector<dim> dFlux;
						compute_flux_derivative(dFlux, vVolumeGrad[ip], vDVolumeGrad[sh]);
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


	////////////////////////////////////////////////////////////////////////////
	// Shared phase parameters
	////////////////////////////////////////////////////////////////////////////

	void set_phase_parameters(Interface<dim>* user)
	{
		if(!user) UG_THROW("AvalancheFluxLinker: Interface pointer is null.");
		if(!user->valid()) UG_THROW("AvalancheFluxLinker: Interface parameters have not been initialized.");

		Inter = user;
	}


	////////////////////////////////////////////////////////////////////////////
	// DPSM parameters
	////////////////////////////////////////////////////////////////////////////

	void set_avalanche_coefficient(number nuAv)
	{
		if(nuAv < 0.0)
			UG_THROW("AvalancheFluxLinker: avalanche coefficient must be non-negative.");

		m_nuAv = nuAv;
	}

	void set_angle_of_repose(number thetaCr)
	{
		if(thetaCr <= 0.0 || thetaCr >= 0.5*std::acos(-1.0))
			UG_THROW("AvalancheFluxLinker: angle of repose must be in radians and satisfy 0 < theta < pi/2.");

		m_thetaCr = thetaCr;
		m_criticalSlope = std::tan(thetaCr);
	}

	void set_angle_of_repose_degrees(number thetaCrDeg)
	{
		if(thetaCrDeg <= 0.0 || thetaCrDeg >= 90.0)
			UG_THROW("AvalancheFluxLinker: angle of repose must satisfy 0 < theta < 90 degrees.");

		set_angle_of_repose(thetaCrDeg*std::acos(-1.0)/180.0);
	}

	void set_normal_epsilon(number eps)
	{
		if(eps <= 0.0)
			UG_THROW("AvalancheFluxLinker: normal epsilon must be positive.");

		m_epsNormal = eps;
	}

	void set_delta_epsilon(number eps)
	{
		if(eps <= 0.0)
			UG_THROW("AvalancheFluxLinker: delta epsilon must be positive.");

		m_epsDelta = eps;
	}

	void set_vertical_epsilon(number eps)
	{
		if(eps <= 0.0)
			UG_THROW("AvalancheFluxLinker: vertical epsilon must be positive.");

		m_epsVertical = eps;
	}


protected:

	////////////////////////////////////////////////////////////////////////////
	// Input indices
	////////////////////////////////////////////////////////////////////////////

	static const size_t _DC_ = 0;


	////////////////////////////////////////////////////////////////////////////
	// grad(c)
	////////////////////////////////////////////////////////////////////////////

	SmartPtr<CplUserData<MathVector<dim>, dim> > m_spVolumeGrad;
	SmartPtr<DependentUserData<MathVector<dim>, dim> > m_spDVolumeGrad;


	////////////////////////////////////////////////////////////////////////////
	// Shared phase parameters
	////////////////////////////////////////////////////////////////////////////

	Interface<dim>* Inter;


	////////////////////////////////////////////////////////////////////////////
	// DPSM parameters
	////////////////////////////////////////////////////////////////////////////

	number m_nuAv;
	number m_thetaCr;
	number m_criticalSlope;
	number m_epsNormal;
	number m_epsDelta;
	number m_epsVertical;
};


} // namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__AVALANCHE_FLUX_LINKER__ */
