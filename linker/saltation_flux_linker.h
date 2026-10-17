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
 *
 * The surface saltation flux magnitude is
 *
 *      q_s = (C/g) (u_* - u_*t)_+ tau_salt,
 *
 * where
 *
 *      tau_salt = mu_eff sqrt(2 D:D),
 *
 *      mu_eff = mu_a + rho_a nu_t,
 *
 *      D = 1/2 (grad(u) + grad(u)^T),
 *
 *      u_* = sqrt(tau_salt/rho_a).
 *
 * The saltation direction is given by the normalized local velocity,
 *
 *      d_salt = u/sqrt(|u|^2 + eps_u^2).
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

	SaltationFluxLinker() : m_spVolumeGrad(NULL), m_spDVolumeGrad(NULL), m_spVelocityGrad(NULL), m_spDVelocityGrad(NULL), m_spTurbulentKinViscosity(NULL), m_spDTurbulentKinViscosity(NULL), m_spVelocity(NULL), m_spDVelocity(NULL), m_spNormal(NULL), m_spDNormal(NULL), Inter(NULL), m_C(5.5), m_uStarThreshold(0.22), m_epsNormal(1e-8), m_epsDelta(1e-8), m_epsVelocity(1e-4), m_directionMode(VELOCITY_DIRECTION)
	{
		this->set_num_input(5);
	}
	struct SaltationData
	{
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
		number tauSign;
		number sVelocity;
		MathVector<dim> velocityDirection;
	};
	enum DirectionMode
	{
		SHEAR_DIRECTION = 0,
		VELOCITY_ALIGNED_SHEAR_DIRECTION = 1,
		VELOCITY_DIRECTION = 2
	};


private:

	////////////////////////////////////////////////////////////////////////////
	// Flux evaluation
	////////////////////////////////////////////////////////////////////////////

	void compute_flux(MathVector<dim>& flux, const MathVector<dim>& gradC, const MathMatrix<dim,dim>& gradU, number nuT, const MathVector<dim>& velocity, const MathVector<dim>& normal) const
	{
		VecSet(flux, 0.0);

		if(!Inter) UG_THROW("SaltationFluxLinker: Interface pointer is null.");

		const number grav = std::fabs(Inter->gravity());

		if(grav <= 0.0) UG_THROW("SaltationFluxLinker: gravity magnitude must be positive.");

		SaltationData data;
		compute_saltation_data(data, gradC, gradU, nuT, velocity, normal);

		if(data.deltaGamma <= 0.0) return;
		if(data.uStar <= m_uStarThreshold) return;
		
		const number rho_bed = Inter->Density_s()*Inter->packing_factor();
		const number excess = data.uStar - m_uStarThreshold;
		const number factor = data.deltaGamma*(m_C/grav)*excess/rho_bed;

		if(m_directionMode == SHEAR_DIRECTION)
		{
			for(size_t d = 0; d < dim; ++d)
				flux[d] = factor*data.tau[d];
		}
		else if(m_directionMode == VELOCITY_ALIGNED_SHEAR_DIRECTION)
		{
			for(size_t d = 0; d < dim; ++d)
				flux[d] = factor*data.tauSign*data.tau[d];
		}
		else if(m_directionMode == VELOCITY_DIRECTION)
		{
			for(size_t d = 0; d < dim; ++d)
				flux[d] = factor*data.tauMag*data.velocityDirection[d];
		}
	}


	////////////////////////////////////////////////////////////////////////////
	// Analytical directional derivative
	////////////////////////////////////////////////////////////////////////////

	void compute_flux_derivative(MathVector<dim>& dFlux,
								 const MathVector<dim>& gradC,
								 const MathMatrix<dim,dim>& gradU,
								 number nuT,
								 const MathVector<dim>& velocity,
								 const MathVector<dim>& normal,
								 const MathVector<dim>& dGradC,
								 const MathMatrix<dim,dim>& dGradU,
								 const MathVector<dim>& dVelocity,
								 number dNuT = 0.0) const
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
		compute_saltation_data(data, gradC, gradU, nuT, velocity, normal);

		if(data.deltaGamma <= 0.0)
			return;

		if(data.tauMag <= 0.0)
			return;

		if(data.uStar <= m_uStarThreshold)
			return;

		number gradCDotDGradC = 0.0;

		for(size_t d = 0; d < dim; ++d)
			gradCDotDGradC += gradC[d]*dGradC[d];


		MathVector<dim> dn;
		VecSet(dn, 0.0);
		
		number dNormal2 = 0.0;


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

		number D2 = 0.0;
		number DDotDD = 0.0;

		for(size_t i = 0; i < dim; ++i)
		{
			for(size_t j = 0; j < dim; ++j)
			{
				D2 += data.D(i,j)*data.D(i,j);
				DDotDD += data.D(i,j)*dD(i,j);
			}
		}

		const number shearRate = std::sqrt(2.0*D2);

		number dShearRate = 0.0;

		if(shearRate > 0.0)
			dShearRate = 2.0*DDotDD/shearRate;

		const number dTauMag = dMuEff*shearRate + mu_eff*dShearRate;
		
		
		const number dUStar = dTauMag/(2.0*rho_a*data.uStar);
		const number excess = data.uStar - m_uStarThreshold;
		const number rho_bed = Inter->Density_s()*Inter->packing_factor();
		const number Cg = m_C/(grav*rho_bed);

		

		if(m_directionMode == SHEAR_DIRECTION)
		{
			for(size_t d = 0; d < dim; ++d)
				dFlux[d] = Cg*(dDeltaGamma*excess*data.tau[d] + data.deltaGamma*dUStar*data.tau[d] + data.deltaGamma*excess*dTau[d]);
		}
		else if(m_directionMode == VELOCITY_ALIGNED_SHEAR_DIRECTION)
		{
			for(size_t d = 0; d < dim; ++d)
				dFlux[d] = data.tauSign*Cg*(dDeltaGamma*excess*data.tau[d] + data.deltaGamma*dUStar*data.tau[d] + data.deltaGamma*excess*dTau[d]);
		}
		else if(m_directionMode == VELOCITY_DIRECTION)
		{
			number velocityDotDVelocity = 0.0;

			for(size_t d = 0; d < dim; ++d)
				velocityDotDVelocity += velocity[d]*dVelocity[d];

			const number sVelocity3 = data.sVelocity*data.sVelocity*data.sVelocity;

			MathVector<dim> dVelocityDirection;

			for(size_t d = 0; d < dim; ++d)
				dVelocityDirection[d] = dVelocity[d]/data.sVelocity - velocity[d]*velocityDotDVelocity/sVelocity3;

			for(size_t d = 0; d < dim; ++d)
			{
				dFlux[d] = Cg*(
					dDeltaGamma*excess*data.tauMag*data.velocityDirection[d]
					+ data.deltaGamma*dUStar*data.tauMag*data.velocityDirection[d]
					+ data.deltaGamma*excess*dTauMag*data.velocityDirection[d]
					+ data.deltaGamma*excess*data.tauMag*dVelocityDirection[d]);
			}
		}
	}
	
	void compute_saltation_data(SaltationData& data, const MathVector<dim>& gradC, const MathMatrix<dim,dim>& gradU, number nuT, const MathVector<dim>& velocity, const MathVector<dim>& normal) const
	{
		data.sDelta = 0.0;
		data.deltaGamma = 0.0;
		data.normal2 = 0.0;
		data.nDn = 0.0;
		data.tauMag = 0.0;
		data.uStar = 0.0;
		data.sVelocity = 0.0;
		VecSet(data.velocityDirection, 0.0);

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

		data.sDelta = std::sqrt(gradC2 + m_epsDelta*m_epsDelta);
		data.deltaGamma = data.sDelta - m_epsDelta;

		number normal2 = 0.0;

		for(size_t d = 0; d < dim; ++d)
			normal2 += normal[d]*normal[d];

		const number normalMag = std::sqrt(normal2 + m_epsNormal*m_epsNormal);

		for(size_t d = 0; d < dim; ++d)
			data.normal[d] = normal[d]/normalMag;

		data.normal2 = 0.0;

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

		number tauDotVelocity = 0.0;

		for(size_t d = 0; d < dim; ++d)
			tauDotVelocity += data.tau[d]*velocity[d];

		data.tauSign = tauDotVelocity < 0.0 ? -1.0 : 1.0;
		
		number D2 = 0.0;

		for(size_t i = 0; i < dim; ++i)
			for(size_t j = 0; j < dim; ++j)
				D2 += data.D(i,j)*data.D(i,j);

		data.tauMag = mu_eff*std::sqrt(2.0*D2);

		if(data.tauMag > 0.0)
			data.uStar = std::sqrt(data.tauMag/rho_a);
		
		number velocity2 = 0.0;

		for(size_t d = 0; d < dim; ++d)
			velocity2 += velocity[d]*velocity[d];

		data.sVelocity = std::sqrt(velocity2 + m_epsVelocity*m_epsVelocity);

		for(size_t d = 0; d < dim; ++d)
			data.velocityDirection[d] = velocity[d]/data.sVelocity;
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
		MathVector<dim> velocity;
		MathVector<dim> normal;
		
		(*m_spTurbulentKinViscosity)(nuT, globIP, time, si);
		(*m_spVolumeGrad)(gradC, globIP, time, si);
		(*m_spVelocityGrad)(gradU, globIP, time, si);
		(*m_spVelocity)(velocity, globIP, time, si);
		(*m_spNormal)(normal, globIP, time, si);

		compute_flux(value, gradC, gradU, nuT, velocity, normal);
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
		std::vector<MathVector<dim> > vVelocity(nip);
		std::vector<MathVector<dim> > vNormal(nip);
		
		
		(*m_spTurbulentKinViscosity)(&vNuT[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
		(*m_spVolumeGrad)(&vVolumeGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
		(*m_spVelocityGrad)(&vVelocityGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
		(*m_spVelocity)(&vVelocity[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
		(*m_spNormal)(&vNormal[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);

		for(size_t ip = 0; ip < nip; ++ip)
			compute_flux(vValue[ip], vVolumeGrad[ip], vVelocityGrad[ip], vNuT[ip], vVelocity[ip], vNormal[ip]);
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
		const int s_VEL_ = base_type::series_id(_VEL_, s);
		const int s_NORMAL_ = base_type::series_id(_NORMAL_, s);

		const MathVector<dim>* vVolumeGrad = m_spVolumeGrad->values(s_DC_);
		const MathMatrix<dim,dim>* vVelocityGrad = m_spVelocityGrad->values(s_DU_);
		const number* vNuT = m_spTurbulentKinViscosity->values(s_DNUT_);
		const MathVector<dim>* vVelocity = m_spVelocity->values(s_VEL_);
		const MathVector<dim>* vNormal = m_spNormal->values(s_NORMAL_);

		for(size_t ip = 0; ip < nip; ++ip)
			compute_flux(vValue[ip], vVolumeGrad[ip], vVelocityGrad[ip], vNuT[ip], vVelocity[ip], vNormal[ip]);

		if(!bDeriv || this->zero_derivative()) return;

		this->set_zero(vvvDeriv, nip);

		////////////////////////////////////////////////////////////////////////////
		// Derivatives with respect to grad(c)
		////////////////////////////////////////////////////////////////////////////

		if(m_spDVolumeGrad.valid() && !m_spDVolumeGrad->zero_derivative())
		{
			MathMatrix<dim,dim> zeroGradU;
			MathVector<dim> zeroVelocity;
			MathVector<dim> dFlux;
			
			for(size_t ip = 0; ip < nip; ++ip)
			{
				for(size_t fct = 0; fct < m_spDVolumeGrad->num_fct(); ++fct)
				{
					const MathVector<dim>* vDVolumeGrad = m_spDVolumeGrad->deriv(s_DC_, ip, fct);
					const size_t commonFct = this->input_common_fct(_DC_, fct);

					for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
					{

						MatSet(zeroGradU, 0.0);
						VecSet(zeroVelocity, 0.0);

						compute_flux_derivative(dFlux, vVolumeGrad[ip], vVelocityGrad[ip], vNuT[ip], vVelocity[ip], vNormal[ip], vDVolumeGrad[sh], zeroGradU, zeroVelocity);

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
			MathVector<dim> zeroGradC;
			MathVector<dim> zeroVelocity;
			MathVector<dim> dFlux;
			
			for(size_t ip = 0; ip < nip; ++ip)
			{
				for(size_t fct = 0; fct < m_spDVelocityGrad->num_fct(); ++fct)
				{
					const MathMatrix<dim,dim>* vDVelocityGrad = m_spDVelocityGrad->deriv(s_DU_, ip, fct);
					const size_t commonFct = this->input_common_fct(_DU_, fct);

					for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
					{

						VecSet(zeroGradC, 0.0);
						VecSet(zeroVelocity, 0.0);

						compute_flux_derivative(dFlux, vVolumeGrad[ip], vVelocityGrad[ip], vNuT[ip], vVelocity[ip], vNormal[ip], zeroGradC, vDVelocityGrad[sh], zeroVelocity);

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
			MathVector<dim> zeroVelocity;
			MathVector<dim> dFlux;

			for(size_t ip = 0; ip < nip; ++ip)
			{
				for(size_t fct = 0; fct < m_spDTurbulentKinViscosity->num_fct(); ++fct)
				{
					const number* vDNuT = m_spDTurbulentKinViscosity->deriv(s_DNUT_, ip, fct);
					const size_t commonFct = this->input_common_fct(_DNUT_, fct);

					for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
					{



						VecSet(zeroGradC, 0.0);
						VecSet(zeroVelocity, 0.0);
						MatSet(zeroGradU, 0.0);

						compute_flux_derivative(dFlux, vVolumeGrad[ip], vVelocityGrad[ip], vNuT[ip], vVelocity[ip], vNormal[ip], zeroGradC, zeroGradU, zeroVelocity, vDNuT[sh]);
						vvvDeriv[ip][commonFct][sh] += dFlux;
					}
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////
		// Derivatives with respect to velocity
		////////////////////////////////////////////////////////////////////////////

		if(m_directionMode == VELOCITY_DIRECTION && m_spDVelocity.valid() && !m_spDVelocity->zero_derivative())
		{
			MathVector<dim> zeroGradC;
			MathMatrix<dim,dim> zeroGradU;
			MathVector<dim> dFlux;

			

			for(size_t ip = 0; ip < nip; ++ip)
			{
				for(size_t fct = 0; fct < m_spDVelocity->num_fct(); ++fct)
				{
					const MathVector<dim>* vDVelocity = m_spDVelocity->deriv(s_VEL_, ip, fct);
					const size_t commonFct = this->input_common_fct(_VEL_, fct);

					for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
					{
						
						VecSet(zeroGradC, 0.0);
						MatSet(zeroGradU, 0.0);
						
						compute_flux_derivative(dFlux, vVolumeGrad[ip], vVelocityGrad[ip], vNuT[ip], vVelocity[ip], vNormal[ip], zeroGradC, zeroGradU, vDVelocity[sh]);

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
	void set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
	{
		if(!data.valid()) UG_THROW("SaltationFluxLinker: velocity input is null.");

		m_spVelocity = data;
		m_spDVelocity = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim> >();
		base_type::set_input(_VEL_, data, data);
	}
	void set_normal(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
	{
		if(!data.valid()) UG_THROW("SaltationFluxLinker: normal input is null.");

		m_spNormal = data;
		m_spDNormal = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim> >();
		base_type::set_input(_NORMAL_, data, data);
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
	
	void set_velocity_epsilon(number eps)
	{
		if(!(eps > 0.0) || !std::isfinite(eps))
			UG_THROW("SaltationFluxLinker: velocity epsilon must be finite and positive.");

		m_epsVelocity = eps;
	}
	
	void set_direction_mode(int mode)
	{
		if(mode < 0 || mode > 2)
			UG_THROW("SaltationFluxLinker: direction mode must be 0, 1, or 2.");

		m_directionMode = mode;
	}


protected:

	////////////////////////////////////////////////////////////////////////////
	// Input indices
	////////////////////////////////////////////////////////////////////////////

	static const size_t _DC_ = 0;
	static const size_t _DU_ = 1;
	static const size_t _DNUT_ = 2;
	static const size_t _VEL_ = 3;
	static const size_t _NORMAL_ = 4;


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
	// Vel u
	////////////////////////////////////////////////////////////////////////////
	
	SmartPtr<CplUserData<MathVector<dim>, dim> > m_spVelocity;
	SmartPtr<DependentUserData<MathVector<dim>, dim> > m_spDVelocity;
	
	////////////////////////////////////////////////////////////////////////////
	// normal n
	////////////////////////////////////////////////////////////////////////////
	
	SmartPtr<CplUserData<MathVector<dim>, dim> > m_spNormal;
	SmartPtr<DependentUserData<MathVector<dim>, dim> > m_spDNormal;


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
	number m_epsVelocity;
	int m_directionMode;
};


} // namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__SALTATION_FLUX_LINKER__ */
