/*
 * This file is part of the UG4 NavierStokes plugin.
 */

#include "rans_turbulence_fv1.h"

#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
RANSTurbulenceFV1<TDomain>::
RANSTurbulenceFV1(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset)
	: IElemDisc<TDomain>(vFct, vSubset)
{
	std::string functions;

	for(size_t i = 0; i < vFct.size(); ++i)
	{
		if(i > 0) functions.append(",");
		functions.append(vFct[i]);
	}

	init(functions);
}
template <typename TDomain>
RANSTurbulenceFV1<TDomain>::
RANSTurbulenceFV1(const char* functions, const char* subsets)
	: IElemDisc<TDomain>(functions, subsets)
{
	init(functions);
}

template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::
init(const std::string& functions)
{

	if(this->num_fct() != 2)
		UG_THROW("RANSTurbulenceFV1: Exactly two symbolic functions are required: k and omega.");

	m_exTurbulentKinViscosity = make_sp(new DataExport<number, dim>(functions.c_str()));
	m_exK = make_sp(new DataExport<number, dim>(functions.c_str()));
	m_exOmega = make_sp(new DataExport<number, dim>(functions.c_str()));
	
	
	m_imVelocity.set_comp_lin_defect(true);
	m_imVelocityGradientSCVF.set_comp_lin_defect(false);
	m_imVelocityGradientSCV.set_comp_lin_defect(false);
	m_imKinViscositySCVF.set_comp_lin_defect(false);
	m_imKinViscositySCV.set_comp_lin_defect(false);
	m_imWallDistanceSCVF.set_comp_lin_defect(false);
	m_imWallDistanceSCV.set_comp_lin_defect(false);

	this->register_import(m_imVelocity);
	this->register_import(m_imVelocityGradientSCVF);
	this->register_import(m_imVelocityGradientSCV);
	this->register_import(m_imKinViscositySCVF);
	this->register_import(m_imKinViscositySCV);
	this->register_import(m_imWallDistanceSCVF);
	this->register_import(m_imWallDistanceSCV);

	register_all_funcs(false);
}

////////////////////////////////////////////////////////////////////////////////
// Turbulence model
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::set_model(TurbulenceModel model)
{
	if(model != K_OMEGA_SST) UG_THROW("RANSTurbulenceFV1: Only SST-2003 is implemented.");
	m_model = model;
}
////////////////////////////////////////////////////////////////////////////////
// Input data
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::
set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
{
	m_imVelocity.set_data(data);
}
template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::
set_velocity(const std::vector<number>& velocity)
{
	if(velocity.size() != dim)
		UG_THROW("RANSTurbulenceFV1: Velocity vector must have "
				 << dim << " components.");

	SmartPtr<ConstUserVector<dim> > data(
		new ConstUserVector<dim>(velocity));

	set_velocity(data);
}

template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::
set_velocity_gradient(
	SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > data)
{
	m_imVelocityGradientSCVF.set_data(data);
	m_imVelocityGradientSCV.set_data(data);
	
	m_exTurbulentKinViscosity->add_needed_data(data);
}
template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::
set_velocity_gradient(const MathMatrix<dim, dim>& velocityGradient)
{
	SmartPtr<ConstUserMatrix<dim> > data(new ConstUserMatrix<dim>());

	for(size_t i = 0; i < dim; ++i)
		for(size_t j = 0; j < dim; ++j)
			data->set_entry(i, j, velocityGradient(i,j));

	set_velocity_gradient(data);
}

template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::
set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data)
{
	m_imKinViscositySCVF.set_data(data);
	m_imKinViscositySCV.set_data(data);
	
	m_exTurbulentKinViscosity->add_needed_data(data);
}
template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::
set_kinematic_viscosity(number viscosity)
{
	SmartPtr<ConstUserNumber<dim> > data(new ConstUserNumber<dim>(viscosity));
	set_kinematic_viscosity(data);
}

template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::
set_wall_distance(SmartPtr<CplUserData<number, dim> > data)
{
	m_imWallDistanceSCVF.set_data(data);
	m_imWallDistanceSCV.set_data(data);
	
	m_exTurbulentKinViscosity->add_needed_data(data);
}



////////////////////////////////////////////////////////////////////////////////
// Discretization setting
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	if(bNonRegularGrid)
		UG_THROW("RANSTurbulenceFV1: Hanging nodes are not implemented.");

	if(vLfeID.size() != 2)
		UG_THROW("RANSTurbulenceFV1: Exactly two functions are required.");

	for(size_t fct = 0; fct < vLfeID.size(); ++fct)
	{
		if(vLfeID[fct].type() != LFEID::LAGRANGE ||
		   vLfeID[fct].order() != 1)
		{
			UG_THROW("RANSTurbulenceFV1: FV1 expects first-order "
					 "Lagrange trial spaces for k and omega.");
		}
	}

	register_all_funcs(false);
}

////////////////////////////////////////////////////////////////////////////////
// Element preparation
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	// 	Only first order implementation
	if(!(TFVGeom::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
				 " Geometry set.");
	
	//	check, that convective upwinding has been set
	if(m_spConvUpwind.invalid())
		UG_THROW("RANSTurbulenceFV1: Upwind method has not been set.");
	
	m_spConvUpwind->template set_geometry_type<TFVGeom>();
	
	if(!m_imVelocity.data_given())
		UG_THROW("RANSTurbulenceFV1: Velocity field has not been set.");
	
	if(!m_imVelocityGradientSCVF.data_given() || !m_imVelocityGradientSCV.data_given())
		UG_THROW("RANSTurbulenceFV1: Velocity gradient has not been set.");

	if(!m_imKinViscositySCVF.data_given() || !m_imKinViscositySCV.data_given())
		UG_THROW("RANSTurbulenceFV1: Kinematic viscosity has not been set.");
	
	if(!m_imWallDistanceSCVF.data_given() || !m_imWallDistanceSCV.data_given())
		UG_THROW("RANSTurbulenceFV1: Wall distance has not been set.");
	
	

	

	
	
	
	//	set local positions for imports
	if(!TFVGeom::usesHangingNodes)
	{
		
		static const int refDim = TElem::dim;
		TFVGeom& geo = GeomProvider<TFVGeom>::get();
		const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();
		
		m_imVelocity.template set_local_ips<refDim>(vSCVFip, numSCVFip);
		m_imVelocityGradientSCVF.template set_local_ips<refDim>(vSCVFip, numSCVFip);
		m_imKinViscositySCVF.template set_local_ips<refDim>(vSCVFip, numSCVFip);
		m_imWallDistanceSCVF.template set_local_ips<refDim>(vSCVFip, numSCVFip);
		
		m_imVelocityGradientSCV.template set_local_ips<refDim>(vSCVip, numSCVip);
		m_imKinViscositySCV.template set_local_ips<refDim>(vSCVip, numSCVip);
		m_imWallDistanceSCV.template set_local_ips<refDim>(vSCVip, numSCVip);
		
	}


}

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::
prep_elem(const LocalVector& u,
		  GridObject* elem,
		  const ReferenceObjectID roid,
		  const MathVector<dim> vCornerCoords[])
{
	TFVGeom& geo = GeomProvider<TFVGeom>::get();

	try
	{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("RANSTurbulenceFV1::prep_elem: "
				   "Cannot update finite-volume geometry.");

	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();

	m_imVelocity.set_global_ips(vSCVFip, numSCVFip);
	m_imVelocityGradientSCVF.set_global_ips(vSCVFip, numSCVFip);
	m_imKinViscositySCVF.set_global_ips(vSCVFip, numSCVFip);
	m_imWallDistanceSCVF.set_global_ips(vSCVFip, numSCVFip);
	
	m_imVelocityGradientSCV.set_global_ips(vSCVip, numSCVip);
	m_imKinViscositySCV.set_global_ips(vSCVip, numSCVip);
	m_imWallDistanceSCV.set_global_ips(vSCVip, numSCVip);
}

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::
fsh_elem_loop()
{
}

////////////////////////////////////////////////////////////////////////////////
// Assembly placeholders
//
// These routines intentionally contain no turbulence-model contributions yet.
// They only establish the FV1 element-discretization interface.
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::
add_jac_A_elem(LocalMatrix& J,
			   const LocalVector& u,
			   GridObject* elem,
			   const MathVector<dim> vCornerCoords[])
{
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// Compute the same upwind interpolation weights used in the defect.
	m_spConvUpwind->update(&geo, m_imVelocity.values());

	const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;

	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		////////////////////////////////////////////////////////////
		// Interpolate k and omega at SCVF
		////////////////////////////////////////////////////////////
		
		number kIP = 0.0;
		number omegaIP = 0.0;
		
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			kIP += scvf.shape(sh) * u(_K_, sh);
			omegaIP += scvf.shape(sh) * u(_OMEGA_, sh);
		}
		////////////////////////////////////////////////////////////
		// Gradients of k and omega
		////////////////////////////////////////////////////////////
		
		MathVector<dim> gradK;
		MathVector<dim> gradOmega;

		VecSet(gradK, 0.0);
		VecSet(gradOmega, 0.0);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			for(int d1 = 0; d1 < dim; ++d1)
			{
				gradK[d1] += scvf.global_grad(sh)[d1] * u(_K_, sh);
				gradOmega[d1] += scvf.global_grad(sh)[d1] * u(_OMEGA_, sh);
			}
		}
		
		////////////////////////////////////////////////////////////
		// SST blending functions and turbulent viscosity
		////////////////////////////////////////////////////////////
		
		const number CDkw = cross_diffusion_CD(omegaIP, gradK, gradOmega);
		const number F1 = blending_function_F1(kIP, omegaIP, m_imKinViscositySCVF[ip], m_imWallDistanceSCVF[ip], CDkw);
		const number F2 = blending_function_F2(kIP, omegaIP, m_imKinViscositySCVF[ip], m_imWallDistanceSCVF[ip]);
		
		const number strainMag = strain_rate_magnitude(m_imVelocityGradientSCVF[ip]);
		const number nuT = turbulent_kinematic_viscosity(kIP, omegaIP, strainMag, F2);
		
	
		////////////////////////////////////////////////////////////
		// SST diffusion coefficients
		////////////////////////////////////////////////////////////
		
		const number sigmaK = blend_sst_coefficient(F1, m_sigmaK1, m_sigmaK2);
		const number sigmaOmega = blend_sst_coefficient(F1, m_sigmaOmega1, m_sigmaOmega2);
		
		const number nuEffK = m_imKinViscositySCVF[ip] + sigmaK * nuT;
		const number nuEffOmega = m_imKinViscositySCVF[ip] + sigmaOmega * nuT;
		
		number dNuT_dKIP = 0.0;
		number dNuT_dOmegaIP = 0.0;

		if(m_linearizeTurbulentViscosity) turbulent_viscosity_derivatives(dNuT_dKIP, dNuT_dOmegaIP, kIP, omegaIP, strainMag, F2, m_imKinViscositySCVF[ip], m_imWallDistanceSCVF[ip]);
		
		
		////////////////////////////////////////////////////////////
		// Convective volume flux
		////////////////////////////////////////////////////////////
		
		const number volFlux = VecDot(m_imVelocity[ip], scvf.normal());

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			// F1, sigmaK, and sigmaOmega are currently frozen in the Jacobian.
			// The F2 contribution to d(nu_t)/d(k,omega) can optionally be linearized.
			// All max() branches are frozen during the linearization.
			
			////////////////////////////////////////////////////////
			// Diffusion
			////////////////////////////////////////////////////////

			const number diffFluxShapeK = -nuEffK * VecDot(scvf.global_grad(sh), scvf.normal());
			const number diffFluxShapeOmega = -nuEffOmega * VecDot(scvf.global_grad(sh), scvf.normal());

			J(_K_, scvf.from(), _K_, sh) += diffFluxShapeK;
			J(_K_, scvf.to(), _K_, sh) -= diffFluxShapeK;
			
			J(_OMEGA_, scvf.from(), _OMEGA_, sh) += diffFluxShapeOmega;
			J(_OMEGA_, scvf.to(), _OMEGA_, sh) -= diffFluxShapeOmega;
			
			
			
			////////////////////////////////////////////////////////
			// Optional linearization of turbulent viscosity
			////////////////////////////////////////////////////////
			
			
			if(m_linearizeTurbulentViscosity)
			{
				const number dNuEffK_dK = sigmaK * dNuT_dKIP * scvf.shape(sh);
				const number diffFluxKK = -dNuEffK_dK * VecDot(gradK, scvf.normal());

				J(_K_, scvf.from(), _K_, sh) += diffFluxKK;
				J(_K_, scvf.to(), _K_, sh) -= diffFluxKK;


				const number dNuEffK_dOmega = sigmaK * dNuT_dOmegaIP * scvf.shape(sh);
				const number diffFluxKOmega = -dNuEffK_dOmega * VecDot(gradK, scvf.normal());

				J(_K_, scvf.from(), _OMEGA_, sh) += diffFluxKOmega;
				J(_K_, scvf.to(), _OMEGA_, sh) -= diffFluxKOmega;


				const number dNuEffOmega_dK = sigmaOmega * dNuT_dKIP * scvf.shape(sh);
				const number diffFluxOmegaK = -dNuEffOmega_dK * VecDot(gradOmega, scvf.normal());

				J(_OMEGA_, scvf.from(), _K_, sh) += diffFluxOmegaK;
				J(_OMEGA_, scvf.to(), _K_, sh) -= diffFluxOmegaK;


				const number dNuEffOmega_dOmega = sigmaOmega * dNuT_dOmegaIP * scvf.shape(sh);
				const number diffFluxOmegaOmegaCoeff = -dNuEffOmega_dOmega * VecDot(gradOmega, scvf.normal());

				J(_OMEGA_, scvf.from(), _OMEGA_, sh) += diffFluxOmegaOmegaCoeff;
				J(_OMEGA_, scvf.to(), _OMEGA_, sh) -= diffFluxOmegaOmegaCoeff;
			}
			


			////////////////////////////////////////////////////////
			// Convection
			////////////////////////////////////////////////////////

			const number convFluxShape = volFlux * upwind.upwind_shape_sh(ip, sh);

			J(_K_, scvf.from(), _K_, sh) += convFluxShape;
			J(_K_, scvf.to(), _K_, sh) -= convFluxShape;

			J(_OMEGA_, scvf.from(), _OMEGA_, sh) += convFluxShape;
			J(_OMEGA_, scvf.to(), _OMEGA_, sh) -= convFluxShape;
		}
	}
	
	////////////////////////////////////////////////////////////
	// k- and omega-equation source Jacobian
	////////////////////////////////////////////////////////////

	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		const typename TFVGeom::SCV& scv = geo.scv(ip);
		const size_t co = scv.node_id();
		const number volume = scv.volume();
		const number kRaw = u(_K_, co);
		const number omegaRaw = u(_OMEGA_, co);
		const number k = std::max(kRaw, 0.0);
		const number omega = std::max(omegaRaw, 1.0e-12);
		const bool omegaActive = omegaRaw > 1.0e-12;
		MathVector<dim> gradK;
		MathVector<dim> gradOmega;
		VecSet(gradK, 0.0);
		VecSet(gradOmega, 0.0);
		for(size_t sh = 0; sh < scv.num_sh(); ++sh)
		{
			for(int d1 = 0; d1 < dim; ++d1)
			{
				gradK[d1] += scv.global_grad(sh)[d1] * u(_K_, sh);
				gradOmega[d1] += scv.global_grad(sh)[d1] * u(_OMEGA_, sh);
			}
		}
		const number CDkw = cross_diffusion_CD(omega, gradK, gradOmega);
		const number F1 = blending_function_F1(k, omega, m_imKinViscositySCV[ip], m_imWallDistanceSCV[ip], CDkw);
		const number beta = blend_sst_coefficient(F1, m_beta1, m_beta2);
		J(_K_, co, _K_, co) += m_betaStar * omega * volume;
		if(m_linearizeDestructionCoupling && omegaActive) J(_K_, co, _OMEGA_, co) += m_betaStar * kRaw * volume;
		J(_OMEGA_, co, _OMEGA_, co) += beta * (omega + (omegaActive ? omegaRaw : 0.0)) * volume;
		const number crossDiffusionOmega = cross_diffusion_omega(F1, omega, gradK, gradOmega);
		const number crossCoeff = crossDiffusionOmega / omega;
		if(m_crossDiffusionLinearization == 0)
		{
			J(_OMEGA_, co, _OMEGA_, co) += std::max(-crossCoeff, 0.0) * volume;
		}
		else if(m_crossDiffusionLinearization == 1)
		{
			if(omegaActive) J(_OMEGA_, co, _OMEGA_, co) += std::max(crossCoeff, 0.0) * volume;
		}
		else
		{
			const number crossFactor = 2.0 * (1.0 - F1) * m_sigmaOmega2;
			for(size_t sh = 0; sh < scv.num_sh(); ++sh)
			{
				J(_OMEGA_, co, _K_, sh) -= crossFactor * VecDot(scv.global_grad(sh), gradOmega) * volume / omega;
				J(_OMEGA_, co, _OMEGA_, sh) -= crossFactor * VecDot(gradK, scv.global_grad(sh)) * volume / omega;
			}
			if(omegaActive) J(_OMEGA_, co, _OMEGA_, co) += crossCoeff * volume;
		}
	}
}

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::
add_jac_M_elem(LocalMatrix& J,
			   const LocalVector& u,
			   GridObject* elem,
			   const MathVector<dim> vCornerCoords[])
{
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		const size_t co = scv.node_id();
		const number volume = scv.volume();

		J(_K_, co, _K_, co) += volume;
		J(_OMEGA_, co, _OMEGA_, co) += volume;
	}
}

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::
add_def_A_elem(LocalVector& d,
			   const LocalVector& u,
			   GridObject* elem,
			   const MathVector<dim> vCornerCoords[])
{
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// Compute upwind interpolation weights using the prescribed
	// velocity evaluated at the SCVFs.
	m_spConvUpwind->update(&geo, m_imVelocity.values());

	const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;

	// Loop over sub-control-volume faces.
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		
		
		////////////////////////////////////////////////////////////
		// Interpolate k and omega at SCVF
		////////////////////////////////////////////////////////////
		
		number kIP = 0.0;
		number omegaIP = 0.0;
		
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			kIP += scvf.shape(sh) * u(_K_, sh);
			omegaIP += scvf.shape(sh) * u(_OMEGA_, sh);
		}
		////////////////////////////////////////////////////////////
		// Gradients of k and omega
		////////////////////////////////////////////////////////////
		
		MathVector<dim> gradK;
		MathVector<dim> gradOmega;

		VecSet(gradK, 0.0);
		VecSet(gradOmega, 0.0);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			for(int d1 = 0; d1 < dim; ++d1)
			{
				gradK[d1] += scvf.global_grad(sh)[d1] * u(_K_, sh);
				gradOmega[d1] += scvf.global_grad(sh)[d1] * u(_OMEGA_, sh);
			}
		}
		
		////////////////////////////////////////////////////////////
		// SST blending functions and turbulent viscosity
		////////////////////////////////////////////////////////////
		
		const number CDkw = cross_diffusion_CD(omegaIP, gradK, gradOmega);
		const number F1 = blending_function_F1(kIP, omegaIP, m_imKinViscositySCVF[ip], m_imWallDistanceSCVF[ip], CDkw);
		const number F2 = blending_function_F2(kIP, omegaIP, m_imKinViscositySCVF[ip], m_imWallDistanceSCVF[ip]);
		
		const number strainMag = strain_rate_magnitude(m_imVelocityGradientSCVF[ip]);
		const number nuT = turbulent_kinematic_viscosity(kIP, omegaIP, strainMag, F2);
		
		
		////////////////////////////////////////////////////////////
		// SST diffusion coefficients
		////////////////////////////////////////////////////////////
		
		const number sigmaK = blend_sst_coefficient(F1, m_sigmaK1, m_sigmaK2);
		const number sigmaOmega = blend_sst_coefficient(F1, m_sigmaOmega1, m_sigmaOmega2);
		
		const number nuEffK = m_imKinViscositySCVF[ip] + sigmaK * nuT;
		const number nuEffOmega = m_imKinViscositySCVF[ip] + sigmaOmega * nuT;
		
		////////////////////////////////////////////////////////////
		// Diffusion
		////////////////////////////////////////////////////////////


		const number diffFluxK = -nuEffK * VecDot(gradK, scvf.normal());
		const number diffFluxOmega = -nuEffOmega * VecDot(gradOmega, scvf.normal());

		d(_K_, scvf.from()) += diffFluxK;
		d(_K_, scvf.to()) -= diffFluxK;

		d(_OMEGA_, scvf.from()) += diffFluxOmega;
		d(_OMEGA_, scvf.to()) -= diffFluxOmega;

		////////////////////////////////////////////////////////////
		// Convection
		////////////////////////////////////////////////////////////

		const number volFlux = VecDot(m_imVelocity[ip], scvf.normal());

		const number kUp = upwind.upwind_value(ip, u, _K_);
		const number omegaUp = upwind.upwind_value(ip, u, _OMEGA_);

		const number fluxK = volFlux * kUp;
		const number fluxOmega = volFlux * omegaUp;

		d(_K_, scvf.from()) += fluxK;
		d(_K_, scvf.to()) -= fluxK;

		d(_OMEGA_, scvf.from()) += fluxOmega;
		d(_OMEGA_, scvf.to()) -= fluxOmega;
	}
	
	////////////////////////////////////////////////////////////
	// k- and omega-equation production and destruction
	////////////////////////////////////////////////////////////

	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		const size_t co = scv.node_id();
		const number volume = scv.volume();

		const number k = std::max(u(_K_, co), 0.0);
		const number omega = std::max(u(_OMEGA_, co), 1.0e-12);

		////////////////////////////////////////////////////////////
		// Gradients of k and omega at SCV
		////////////////////////////////////////////////////////////

		MathVector<dim> gradK;
		MathVector<dim> gradOmega;

		VecSet(gradK, 0.0);
		VecSet(gradOmega, 0.0);

		for(size_t sh = 0; sh < scv.num_sh(); ++sh)
		{
			for(int d1 = 0; d1 < dim; ++d1)
			{
				gradK[d1] += scv.global_grad(sh)[d1] * u(_K_, sh);
				gradOmega[d1] += scv.global_grad(sh)[d1] * u(_OMEGA_, sh);
			}
		}

		////////////////////////////////////////////////////////////
		// SST blending functions and turbulent viscosity
		////////////////////////////////////////////////////////////

		const number CDkw = cross_diffusion_CD(omega, gradK, gradOmega);
		const number F1 = blending_function_F1(k, omega, m_imKinViscositySCV[ip], m_imWallDistanceSCV[ip], CDkw);
		const number F2 = blending_function_F2(k, omega, m_imKinViscositySCV[ip], m_imWallDistanceSCV[ip]);

		const number strainMag = strain_rate_magnitude(m_imVelocityGradientSCV[ip]);
		const number nuT = turbulent_kinematic_viscosity(k, omega, strainMag, F2);

		////////////////////////////////////////////////////////////
		// SST coefficients
		////////////////////////////////////////////////////////////

		const number beta = blend_sst_coefficient(F1, m_beta1, m_beta2);
		const number gamma = blend_sst_coefficient(F1, m_gamma1, m_gamma2);

		////////////////////////////////////////////////////////////
		// k-equation production and destruction
		////////////////////////////////////////////////////////////

		const number production = nuT * strainMag * strainMag;
		const number productionLimit = m_productionLimiter * m_betaStar * k * omega;
		const number limitedProduction = std::min(production, productionLimit);

		const number destructionK = m_betaStar * u(_K_, co) * omega;

		d(_K_, co) += (destructionK - limitedProduction) * volume;

		////////////////////////////////////////////////////////////
		// omega-equation production, destruction and cross-diffusion
		////////////////////////////////////////////////////////////

		const number destructionOmega = beta * omega * u(_OMEGA_, co);
		const number crossDiffusionOmega = cross_diffusion_omega(F1, omega, gradK, gradOmega);
		
		const number nuTDenominator = std::max(m_a1 * omega, strainMag * F2);
		const number productionOmega = gamma * std::min(strainMag * strainMag, (m_productionLimiter * m_betaStar / m_a1) * omega * nuTDenominator);
		
		d(_OMEGA_, co) += (destructionOmega - productionOmega - crossDiffusionOmega) * volume;
	}
}

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::
add_def_M_elem(LocalVector& d,
			   const LocalVector& u,
			   GridObject* elem,
			   const MathVector<dim> vCornerCoords[])
{
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
		const typename TFVGeom::SCV& scv = geo.scv(ip);

		const size_t co = scv.node_id();
		const number volume = scv.volume();

		d(_K_, co) += u(_K_, co) * volume;
		d(_OMEGA_, co) += u(_OMEGA_, co) * volume;
	}
}

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::
add_rhs_elem(LocalVector& d,
			 GridObject* elem,
			 const MathVector<dim> vCornerCoords[])
{
}

template <typename TDomain>
number RANSTurbulenceFV1<TDomain>::
strain_rate_magnitude(const MathMatrix<dim, dim>& gradU) const
{
	number sum = 0.0;

	for(int i = 0; i < dim; ++i)
	{
		for(int j = 0; j < dim; ++j)
		{
			const number Sij =
				0.5 * (gradU(i,j) + gradU(j,i));

			sum += Sij * Sij;
		}
	}

	return std::sqrt(2.0 * sum);
}

template <typename TDomain>
number RANSTurbulenceFV1<TDomain>::
vorticity_magnitude(const MathMatrix<dim, dim>& gradU) const
{
	number sum = 0.0;

	for(int i = 0; i < dim; ++i)
	{
		for(int j = 0; j < dim; ++j)
		{
			const number Omegaij =
				0.5 * (gradU(i,j) - gradU(j,i));

			sum += Omegaij * Omegaij;
		}
	}

	return std::sqrt(2.0 * sum);
}
template <typename TDomain>
number RANSTurbulenceFV1<TDomain>::
turbulent_kinematic_viscosity(
	number k,
	number omega,
	number limiterMag,
	number F2) const
{

	const number kEff = std::max(k, 0.0);
	const number omegaEff = std::max(omega, 1.0e-12);

	const number denominator =
		std::max(m_a1 * omegaEff, limiterMag * F2);

	return m_a1 * kEff / denominator;
}

template <typename TDomain>
number RANSTurbulenceFV1<TDomain>::
blend_sst_coefficient(number F1, number innerValue, number outerValue) const
{
	return F1 * innerValue + (1.0 - F1) * outerValue;
}

template <typename TDomain>
number RANSTurbulenceFV1<TDomain>::
blending_function_F1(number k, number omega, number nu, number wallDistance, number CDkw) const
{
	if(!std::isfinite(k) || !std::isfinite(omega) || !std::isfinite(nu) || !std::isfinite(wallDistance) || nu < 0.0 || wallDistance < 0.0) UG_THROW("RANSTurbulenceFV1: Non-finite SST input or negative molecular viscosity/wall distance.");

	const number kEff = std::max(k, 0.0);
	const number omegaEff = std::max(omega, 1.0e-12);
	const number dEff = std::max(wallDistance, 1.0e-12);
	const number nuEff = std::max(nu, 0.0);
	const number CDEff = std::max(CDkw, 1.0e-10);

	const number arg1_1 = std::sqrt(kEff) / (m_betaStar * omegaEff * dEff);
	const number arg1_2 = 500.0 * nuEff / (dEff * dEff * omegaEff);
	const number arg1_3 = 4.0 * m_sigmaOmega2 * kEff / (CDEff * dEff * dEff);

	const number arg1 = std::min(std::min(std::max(arg1_1, arg1_2), arg1_3), 10.0);
	const number arg1Squared = arg1 * arg1;

	return std::tanh(arg1Squared * arg1Squared);
}

template <typename TDomain>
number RANSTurbulenceFV1<TDomain>::
blending_function_F2(number k, number omega, number nu, number wallDistance) const
{
	if(!std::isfinite(k) || !std::isfinite(omega) || !std::isfinite(nu) || !std::isfinite(wallDistance) || nu < 0.0 || wallDistance < 0.0) UG_THROW("RANSTurbulenceFV1: Non-finite SST input or negative molecular viscosity/wall distance.");

	const number kEff = std::max(k, 0.0);
	const number omegaEff = std::max(omega, 1.0e-12);
	const number dEff = std::max(wallDistance, 1.0e-12);

	const number nuEff = std::max(nu, 0.0);
	const number arg2_1 = 2.0 * std::sqrt(kEff) / (m_betaStar * omegaEff * dEff);
	const number arg2_2 = 500.0 * nuEff / (dEff * dEff * omegaEff);
	const number arg2 = std::min(std::max(arg2_1, arg2_2), 100.0);

	return std::tanh(arg2 * arg2);
}

template <typename TDomain>
number RANSTurbulenceFV1<TDomain>::
cross_diffusion_CD(number omega, const MathVector<dim>& gradK, const MathVector<dim>& gradOmega) const
{

	const number omegaEff = std::max(omega, 1.0e-12);
	const number crossDiffusion = 2.0 * m_sigmaOmega2 * VecDot(gradK, gradOmega) / omegaEff;

	return std::max(crossDiffusion, 1.0e-10);
}

template <typename TDomain>
number RANSTurbulenceFV1<TDomain>::
cross_diffusion_omega(number F1, number omega, const MathVector<dim>& gradK, const MathVector<dim>& gradOmega) const
{
	const number omegaEff = std::max(omega, 1.0e-12);

	return 2.0 * (1.0 - F1) * m_sigmaOmega2 * VecDot(gradK, gradOmega) / omegaEff;
}

template<typename TDomain>
void RANSTurbulenceFV1<TDomain>::turbulent_viscosity_derivatives(number& dNuT_dK,number& dNuT_dOmega,number k,number omega,number strainMag,
	number F2,number nu,number wallDist) const
{
	dNuT_dK = 0.0;
	dNuT_dOmega = 0.0;

	////////////////////////////////////////////////////////////
	// Regularized values
	////////////////////////////////////////////////////////////

	const number kEff = std::max(k, 0.0);
	const number omegaEff = std::max(omega, 1.0e-12);
	const number dEff = std::max(wallDist, 1.0e-12);
	const number nuEff = std::max(nu, 0.0);

	////////////////////////////////////////////////////////////
	// Active denominator branch of nu_t
	//
	// nu_t = a1*k / max(a1*omega, S*F2)
	////////////////////////////////////////////////////////////

	const number omegaBranch = m_a1 * omegaEff;
	const number strainBranch = strainMag * F2;
	const number denominator = std::max(omegaBranch, strainBranch);

	if(k <= 0.0)
		return;

	////////////////////////////////////////////////////////////
	// Basic derivatives
	//
	// F2 is frozen.
	////////////////////////////////////////////////////////////

	dNuT_dK = m_a1 / denominator;

	if(omega > 1.0e-12 && omegaBranch >= strainBranch)
	{
		dNuT_dOmega =
			-m_a1 * m_a1 * kEff /
			(denominator * denominator);
	}

	////////////////////////////////////////////////////////////
	// Optional F2 contribution
	//
	// Only contributes when the active nu_t denominator is S*F2.
	////////////////////////////////////////////////////////////

	if(!m_linearizeF2InTurbulentViscosity || strainBranch <= omegaBranch)
		return;

	////////////////////////////////////////////////////////////
	// F2 = tanh(arg2^2)
	//
	// arg2 = max(
	//     2*sqrt(k)/(betaStar*omega*d),
	//     500*nu/(d^2*omega)
	// )
	////////////////////////////////////////////////////////////

	const number arg2K = 2.0 * std::sqrt(kEff) / (m_betaStar * omegaEff * dEff);
	const number arg2Nu = 500.0 * nuEff / (dEff * dEff * omegaEff);
	const number arg2 = std::max(arg2K, arg2Nu);
	
	if(arg2 >= 100.0 || F2 >= 1.0) return;

	number dArg2_dK = 0.0;
	number dArg2_dOmega = 0.0;

	////////////////////////////////////////////////////////////
	// Freeze active max() branch of arg2
	////////////////////////////////////////////////////////////

	if(arg2K >= arg2Nu)
	{
		dArg2_dK = 1.0 / (m_betaStar * omegaEff * dEff * std::sqrt(kEff));
		if(omega > 1.0e-12)
			dArg2_dOmega = -arg2K / omegaEff;
	}
	else
	{
		if(omega > 1.0e-12)
			dArg2_dOmega = -arg2Nu / omegaEff;
	}

	////////////////////////////////////////////////////////////
	// Derivatives of F2
	////////////////////////////////////////////////////////////

	const number dF2_dArg2 = 2.0 * arg2 * (1.0 - F2 * F2);
	const number dF2_dK = dF2_dArg2 * dArg2_dK;
	const number dF2_dOmega = dF2_dArg2 * dArg2_dOmega;

	////////////////////////////////////////////////////////////
	// Include derivative of F2 in
	//
	// nu_t = a1*k / (S*F2)
	////////////////////////////////////////////////////////////

	const number denominatorSquared = denominator * denominator;
	dNuT_dK = m_a1 / denominator - m_a1 * kEff * strainMag * dF2_dK / denominatorSquared;
	dNuT_dOmega = -m_a1 * kEff * strainMag * dF2_dOmega / denominatorSquared;
}

template <typename TDomain>
number RANSTurbulenceFV1<TDomain>::
evaluate_turbulent_kinematic_viscosity(number k, number omega, const MathMatrix<dim, dim>& gradU, number nu, number wallDistance) const
{
	const number kEff = std::max(k, 0.0);
	const number omegaEff = std::max(omega, 1.0e-12);
	const number nuEff = std::max(nu, 0.0);
	const number wallDistanceEff = std::max(wallDistance, 1.0e-12);

	const number strainMag = strain_rate_magnitude(gradU);
	const number F2 = blending_function_F2(kEff, omegaEff, nuEff, wallDistanceEff);

	return turbulent_kinematic_viscosity(kEff, omegaEff, strainMag, F2);
}

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::ex_turbulent_kinematic_viscosity(number vValue[], const MathVector<dim> vGlobIP[], number time, int si, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const MathVector<TFVGeom::dim> vLocIP[], const size_t nip, bool bDeriv, std::vector<std::vector<number> > vvvDeriv[])
{
	if(nip == 0) return;
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t numSH = ref_elem_type::numCorners;
	const bool atSCVF = vLocIP == geo.scvf_local_ips();
	const bool atSCV = vLocIP == geo.scv_local_ips();
	UG_ASSERT(!atSCVF || nip == geo.num_scvf(), "SCVF export integration-point mismatch.");
	UG_ASSERT(!atSCV || nip == geo.num_scv(), "SCV export integration-point mismatch.");
	LocalVector* uAux = const_cast<LocalVector*>(&u);
	const FunctionIndexMapping& turbulenceMap = m_exTurbulentKinViscosity->map();
	struct MappingRestore
	{
		LocalVector& vector;
		const FunctionIndexMapping& map;
		MappingRestore(LocalVector& v, const FunctionIndexMapping& m) : vector(v), map(m) {}
		~MappingRestore() {vector.access_by_map(map);}
	} restore(*uAux, turbulenceMap);
	uAux->access_by_map(turbulenceMap);
	LagrangeP1<ref_elem_type>& trialSpace = Provider<LagrangeP1<ref_elem_type> >::get();
	std::vector<number> vK(nip, 0.0);
	std::vector<number> vOmega(nip, 0.0);
	std::vector<number> vNu(nip);
	std::vector<number> vWallDistance(nip);
	std::vector<MathMatrix<dim, dim> > vGradU(nip);
	std::vector<std::vector<number> > vShapes(nip, std::vector<number>(numSH, 0.0));
	for(size_t ip = 0; ip < nip; ++ip)
	{
		if(atSCVF)
		{
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
			for(size_t sh = 0; sh < numSH; ++sh) vShapes[ip][sh] = scvf.shape(sh);
			vGradU[ip] = m_imVelocityGradientSCVF[ip];
			vNu[ip] = m_imKinViscositySCVF[ip];
			vWallDistance[ip] = m_imWallDistanceSCVF[ip];
		}
		else if(atSCV)
		{
			const size_t co = geo.scv(ip).node_id();
			vShapes[ip][co] = 1.0;
			vGradU[ip] = m_imVelocityGradientSCV[ip];
			vNu[ip] = m_imKinViscositySCV[ip];
			vWallDistance[ip] = m_imWallDistanceSCV[ip];
		}
		else trialSpace.shapes(&vShapes[ip][0], vLocIP[ip]);
		for(size_t sh = 0; sh < numSH; ++sh)
		{
			vK[ip] += u(_K_, sh) * vShapes[ip][sh];
			vOmega[ip] += u(_OMEGA_, sh) * vShapes[ip][sh];
		}
	}
	if(!atSCVF && !atSCV)
	{
		(*(m_imVelocityGradientSCVF.user_data()))(&vGradU[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, uAux, NULL);
		uAux->access_by_map(turbulenceMap);
		(*(m_imKinViscositySCVF.user_data()))(&vNu[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, uAux, NULL);
		uAux->access_by_map(turbulenceMap);
		(*(m_imWallDistanceSCVF.user_data()))(&vWallDistance[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, uAux, NULL);
		uAux->access_by_map(turbulenceMap);
	}
	for(size_t ip = 0; ip < nip; ++ip)
	{
		vValue[ip] = evaluate_turbulent_kinematic_viscosity(vK[ip], vOmega[ip], vGradU[ip], vNu[ip], vWallDistance[ip]);
		if(!bDeriv) continue;
		for(size_t fct = 0; fct < vvvDeriv[ip].size(); ++fct)
		{
			for(size_t sh = 0; sh < vvvDeriv[ip][fct].size(); ++sh) vvvDeriv[ip][fct][sh] = 0.0;
		}
		if(!m_linearizeExportedViscosity) continue;
		const number strainMag = strain_rate_magnitude(vGradU[ip]);
		const number F2 = blending_function_F2(vK[ip], vOmega[ip], vNu[ip], vWallDistance[ip]);
		number dNuT_dK = 0.0;
		number dNuT_dOmega = 0.0;
		turbulent_viscosity_derivatives(dNuT_dK, dNuT_dOmega, vK[ip], vOmega[ip], strainMag, F2, vNu[ip], vWallDistance[ip]);
		for(size_t sh = 0; sh < numSH; ++sh)
		{
			vvvDeriv[ip][_K_][sh] = dNuT_dK * vShapes[ip][sh];
			vvvDeriv[ip][_OMEGA_][sh] = dNuT_dOmega * vShapes[ip][sh];
		}
	}
}
template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::
ex_turbulent_kinetic_energy(number vValue[],
							const MathVector<dim> vGlobIP[],
							number time, int si,
							const LocalVector& u,
							GridObject* elem,
							const MathVector<dim> vCornerCoords[],
							const MathVector<TFVGeom::dim> vLocIP[],
							const size_t nip,
							bool bDeriv,
							std::vector<std::vector<number> > vvvDeriv[])
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t numSH = ref_elem_type::numCorners;

	LagrangeP1<ref_elem_type>& trialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	number vShape[numSH];

	for(size_t ip = 0; ip < nip; ++ip)
	{
		trialSpace.shapes(vShape, vLocIP[ip]);

		vValue[ip] = 0.0;

		for(size_t sh = 0; sh < numSH; ++sh)
			vValue[ip] += u(_K_, sh) * vShape[sh];

		if(bDeriv)
		{
			for(size_t sh = 0; sh < numSH; ++sh)
				vvvDeriv[ip][_K_][sh] = vShape[sh];

			for(size_t sh = 0; sh < vvvDeriv[ip][_OMEGA_].size(); ++sh)
				vvvDeriv[ip][_OMEGA_][sh] = 0.0;
		}
	}
}
template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::
ex_specific_dissipation_rate(number vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si,
							 const LocalVector& u,
							 GridObject* elem,
							 const MathVector<dim> vCornerCoords[],
							 const MathVector<TFVGeom::dim> vLocIP[],
							 const size_t nip,
							 bool bDeriv,
							 std::vector<std::vector<number> > vvvDeriv[])
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t numSH = ref_elem_type::numCorners;

	LagrangeP1<ref_elem_type>& trialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	number vShape[numSH];

	for(size_t ip = 0; ip < nip; ++ip)
	{
		trialSpace.shapes(vShape, vLocIP[ip]);

		vValue[ip] = 0.0;

		for(size_t sh = 0; sh < numSH; ++sh)
			vValue[ip] += u(_OMEGA_, sh) * vShape[sh];

		if(bDeriv)
		{
			for(size_t sh = 0; sh < vvvDeriv[ip][_K_].size(); ++sh)
				vvvDeriv[ip][_K_][sh] = 0.0;

			for(size_t sh = 0; sh < numSH; ++sh)
				vvvDeriv[ip][_OMEGA_][sh] = vShape[sh];
		}
	}
}

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::
lin_def_velocity_convection(const LocalVector& u,
							std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
							const size_t nip)
{
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	UG_ASSERT(nip == geo.num_scvf(), "Number of velocity integration points does not match number of SCVFs.");

	m_spConvUpwind->update(&geo, m_imVelocity.values());

	const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;

	for(size_t ip = 0; ip < nip; ++ip)
	{
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		for(size_t fct = 0; fct < this->num_fct(); ++fct)
		{
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecSet(vvvLinDef[ip][fct][sh], 0.0);
		}

		const number kUp = upwind.upwind_value(ip, u, _K_);
		const number omegaUp = upwind.upwind_value(ip, u, _OMEGA_);

		MathVector<dim> dFluxK_dVelocity;
		MathVector<dim> dFluxOmega_dVelocity;

		VecScale(dFluxK_dVelocity, scvf.normal(), kUp);
		VecScale(dFluxOmega_dVelocity, scvf.normal(), omegaUp);

		VecAdd(vvvLinDef[ip][_K_][scvf.from()], vvvLinDef[ip][_K_][scvf.from()], dFluxK_dVelocity);
		VecSubtract(vvvLinDef[ip][_K_][scvf.to()], vvvLinDef[ip][_K_][scvf.to()], dFluxK_dVelocity);

		VecAdd(vvvLinDef[ip][_OMEGA_][scvf.from()], vvvLinDef[ip][_OMEGA_][scvf.from()], dFluxOmega_dVelocity);
		VecSubtract(vvvLinDef[ip][_OMEGA_][scvf.to()], vvvLinDef[ip][_OMEGA_][scvf.to()], dFluxOmega_dVelocity);
	}
}

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::lin_def_velocity_gradient_scv(const LocalVector& u, std::vector<std::vector<MathMatrix<dim, dim> > > vvvLinDef[], const size_t nip)
{
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	UG_ASSERT(nip == geo.num_scv(), "SCV integration-point mismatch.");
	for(size_t ip = 0; ip < nip; ++ip)
	{
		const typename TFVGeom::SCV& scv = geo.scv(ip);
		const size_t co = scv.node_id();
		const number volume = scv.volume();
		for(size_t fct = 0; fct < this->num_fct(); ++fct)
		{
			for(size_t sh = 0; sh < scv.num_sh(); ++sh) MatSet(vvvLinDef[ip][fct][sh], 0.0);
		}
		const number k = std::max(u(_K_, co), 0.0);
		const number omega = std::max(u(_OMEGA_, co), 1.0e-12);
		MathVector<dim> gradK;
		MathVector<dim> gradOmega;
		VecSet(gradK, 0.0);
		VecSet(gradOmega, 0.0);
		for(size_t sh = 0; sh < scv.num_sh(); ++sh)
		{
			for(int d1 = 0; d1 < dim; ++d1)
			{
				gradK[d1] += scv.global_grad(sh)[d1] * u(_K_, sh);
				gradOmega[d1] += scv.global_grad(sh)[d1] * u(_OMEGA_, sh);
			}
		}
		const number CDkw = cross_diffusion_CD(omega, gradK, gradOmega);
		const number F1 = blending_function_F1(k, omega, m_imKinViscositySCV[ip], m_imWallDistanceSCV[ip], CDkw);
		const number F2 = blending_function_F2(k, omega, m_imKinViscositySCV[ip], m_imWallDistanceSCV[ip]);
		const number gamma = blend_sst_coefficient(F1, m_gamma1, m_gamma2);
		const MathMatrix<dim, dim>& gradU = m_imVelocityGradientSCV[ip];
		const number strainMag = strain_rate_magnitude(gradU);
		const number strainSquared = strainMag * strainMag;
		const bool strainBranch = strainMag * F2 > m_a1 * omega;
		const number denominator = std::max(m_a1 * omega, strainMag * F2);
		const number nuT = turbulent_kinematic_viscosity(k, omega, strainMag, F2);
		const number productionK = nuT * strainSquared;
		const number limitK = m_productionLimiter * m_betaStar * k * omega;
		const number limitOmega = (m_productionLimiter * m_betaStar / m_a1) * omega * denominator;
		number coeffK = 0.0;
		number coeffOmega = 0.0;
		if(k > 0.0 && productionK <= limitK) coeffK = -(strainBranch ? 2.0 : 4.0) * nuT * volume;
		if(strainSquared <= limitOmega) coeffOmega = -4.0 * gamma * volume;
		else if(strainBranch && strainMag > 0.0) coeffOmega = -2.0 * gamma * (m_productionLimiter * m_betaStar / m_a1) * omega * F2 * volume / strainMag;
		for(int i = 0; i < dim; ++i)
		{
			for(int j = 0; j < dim; ++j)
			{
				const number Sij = 0.5 * (gradU(i,j) + gradU(j,i));
				vvvLinDef[ip][_K_][co](i,j) = coeffK * Sij;
				vvvLinDef[ip][_OMEGA_][co](i,j) = coeffOmega * Sij;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template <>
void RANSTurbulenceFV1<Domain1d>::
register_all_funcs(bool bHang)
{
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
	}
	else
	{
		UG_THROW("RANSTurbulenceFV1: Hanging nodes are not implemented.");
	}
}
#endif

#ifdef UG_DIM_2
template <>
void RANSTurbulenceFV1<Domain2d>::
register_all_funcs(bool bHang)
{
	if(!bHang)
	{
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
	}
	else
	{
		UG_THROW("RANSTurbulenceFV1: Hanging nodes are not implemented.");
	}
}
#endif

#ifdef UG_DIM_3
template <>
void RANSTurbulenceFV1<Domain3d>::
register_all_funcs(bool bHang)
{
	if(!bHang)
	{
		register_func<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, FV1Geometry<Prism, dim> >();
		register_func<Pyramid, FV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, FV1Geometry<Hexahedron, dim> >();
	}
	else
	{
		UG_THROW("RANSTurbulenceFV1: Hanging nodes are not implemented.");
	}
}
#endif

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RANSTurbulenceFV1<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;
	
	static const int refDim = reference_element_traits<TElem>::dim;

	this->clear_add_fct(id);

	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFVGeom>);
	
	m_imVelocity.set_fct(id, this, &T::template lin_def_velocity_convection<TElem, TFVGeom>);
	m_imVelocityGradientSCV.set_fct(id, this, &T::template lin_def_velocity_gradient_scv<TElem, TFVGeom>);
	
	m_exTurbulentKinViscosity->template set_fct<T, refDim>(id, this, &T::template ex_turbulent_kinematic_viscosity<TElem, TFVGeom>);
	m_exK->template set_fct<T, refDim>(id, this, &T::template ex_turbulent_kinetic_energy<TElem, TFVGeom>);
	m_exOmega->template set_fct<T, refDim>(id, this, &T::template ex_specific_dissipation_rate<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
// Explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class RANSTurbulenceFV1<Domain1d>;
#endif

#ifdef UG_DIM_2
template class RANSTurbulenceFV1<Domain2d>;
#endif

#ifdef UG_DIM_3
template class RANSTurbulenceFV1<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
