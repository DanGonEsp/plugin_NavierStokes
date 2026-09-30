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
RANSTurbulenceFV1(const char* functions, const char* subsets)
	: IElemDisc<TDomain>(functions, subsets),
	  m_model(K_OMEGA_SST)
{
	if(this->num_fct() != 2)
		UG_THROW("RANSTurbulenceFV1: Exactly two symbolic functions "
				 "are required: k and omega.");

	this->register_import(m_imVelocity);
	this->register_import(m_imVelocityGradient);
	this->register_import(m_imKinViscosity);
	this->register_import(m_imWallDistance);

	register_all_funcs(false);
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
	m_imVelocityGradient.set_data(data);
}

template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::
set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data)
{
	m_imKinViscosity.set_data(data);
}

template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::
set_wall_distance(SmartPtr<CplUserData<number, dim> > data)
{
	m_imWallDistance.set_data(data);
}

////////////////////////////////////////////////////////////////////////////////
// Turbulence model
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void RANSTurbulenceFV1<TDomain>::
set_model(TurbulenceModel model)
{
	m_model = model;
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
	static const int refDim = TElem::dim;

	TFVGeom& geo = GeomProvider<TFVGeom>::get();

	const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
	const size_t numSCVFip = geo.num_scvf_ips();

	m_imVelocity.template set_local_ips<refDim>(vSCVFip, numSCVFip, false);
	m_imVelocityGradient.template set_local_ips<refDim>(vSCVFip, numSCVFip, false);
	m_imKinViscosity.template set_local_ips<refDim>(vSCVFip, numSCVFip, false);
	m_imWallDistance.template set_local_ips<refDim>(vSCVFip, numSCVFip, false);

	if(m_spConvUpwind.valid())
		m_spConvUpwind->template set_geometry_type<TFVGeom>();
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

	m_imVelocity.set_global_ips(vSCVFip, numSCVFip);
	m_imVelocityGradient.set_global_ips(geo.scvf_global_ips(),geo.num_scvf_ips());
	m_imKinViscosity.set_global_ips(vSCVFip, numSCVFip);
	m_imWallDistance.set_global_ips(geo.scvf_global_ips(),geo.num_scvf_ips());
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

	if(!m_imVelocity.data_given())
		UG_THROW("RANSTurbulenceFV1: Velocity field has not been set.");

	if(m_spConvUpwind.invalid())
		UG_THROW("RANSTurbulenceFV1: Upwind method has not been set.");

	// Compute the same upwind interpolation weights used in the defect.
	m_spConvUpwind->update(&geo, m_imVelocity.values());

	const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;

	for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
	{
		// Molecular diffusion
		const number diffFluxShape = -m_imKinViscosity[ip] * VecDot(scvf.global_grad(sh), scvf.normal());

		J(_K_, scvf.from(), _K_, sh) += diffFluxShape;
		J(_K_, scvf.to(),   _K_, sh) -= diffFluxShape;

		J(_OMEGA_, scvf.from(), _OMEGA_, sh) += diffFluxShape;
		J(_OMEGA_, scvf.to(),   _OMEGA_, sh) -= diffFluxShape;

		// Convection
		const number fluxShape = volFlux * upwind.upwind_shape_sh(ip, sh);

		J(_K_, scvf.from(), _K_, sh) += fluxShape;
		J(_K_, scvf.to(),   _K_, sh) -= fluxShape;

		J(_OMEGA_, scvf.from(), _OMEGA_, sh) += fluxShape;
		J(_OMEGA_, scvf.to(),   _OMEGA_, sh) -= fluxShape;
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

	if(!m_imVelocity.data_given())
		UG_THROW("RANSTurbulenceFV1: Velocity field has not been set.");

	if(!m_imKinViscosity.data_given())
		UG_THROW("RANSTurbulenceFV1: Kinematic viscosity has not been set.");

	if(m_spConvUpwind.invalid())
		UG_THROW("RANSTurbulenceFV1: Upwind method has not been set.");

	// Compute upwind interpolation weights using the prescribed
	// velocity evaluated at the SCVFs.
	m_spConvUpwind->update(&geo, m_imVelocity.values());

	const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;

	// Loop over sub-control-volume faces.
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		////////////////////////////////////////////////////////////
		// Molecular diffusion
		////////////////////////////////////////////////////////////

		MathVector<dim> gradK;
		MathVector<dim> gradOmega;

		VecSet(gradK, 0.0);
		VecSet(gradOmega, 0.0);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			for(int d1 = 0; d1 < dim; ++d1)
			{
				gradK[d1] +=
					scvf.global_grad(sh)[d1] * u(_K_, sh);

				gradOmega[d1] +=
					scvf.global_grad(sh)[d1] * u(_OMEGA_, sh);
			}
		}

		const number diffFluxK =
			-m_imKinViscosity[ip] *
			VecDot(gradK, scvf.normal());

		const number diffFluxOmega =
			-m_imKinViscosity[ip] *
			VecDot(gradOmega, scvf.normal());

		d(_K_, scvf.from()) += diffFluxK;
		d(_K_, scvf.to())   -= diffFluxK;

		d(_OMEGA_, scvf.from()) += diffFluxOmega;
		d(_OMEGA_, scvf.to())   -= diffFluxOmega;

		////////////////////////////////////////////////////////////
		// Convection
		////////////////////////////////////////////////////////////

		// Conservative volume flux:
		//
		//      u . n
		//
		// The FV1 normal already contains the SCVF measure.
		const number volFlux =
			VecDot(m_imVelocity[ip], scvf.normal());

		const number kUp =
			upwind.upwind_value(ip, u, _K_);

		const number omegaUp =
			upwind.upwind_value(ip, u, _OMEGA_);

		const number fluxK =
			volFlux * kUp;

		const number fluxOmega =
			volFlux * omegaUp;

		// Conservative FV contribution:
		// what leaves "from" enters "to".
		d(_K_, scvf.from()) += fluxK;
		d(_K_, scvf.to())   -= fluxK;

		d(_OMEGA_, scvf.from()) += fluxOmega;
		d(_OMEGA_, scvf.to())   -= fluxOmega;
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
	const number a1 = 0.31;

	const number kEff = std::max(k, 0.0);
	const number omegaEff = std::max(omega, 1.0e-12);

	const number denominator =
		std::max(a1 * omegaEff, limiterMag * F2);

	return a1 * kEff / denominator;
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

	this->clear_add_fct(id);

	this->set_prep_elem_loop_fct(
		id, &T::template prep_elem_loop<TElem, TFVGeom>);

	this->set_prep_elem_fct(
		id, &T::template prep_elem<TElem, TFVGeom>);

	this->set_fsh_elem_loop_fct(
		id, &T::template fsh_elem_loop<TElem, TFVGeom>);

	this->set_add_jac_A_elem_fct(
		id, &T::template add_jac_A_elem<TElem, TFVGeom>);

	this->set_add_jac_M_elem_fct(
		id, &T::template add_jac_M_elem<TElem, TFVGeom>);

	this->set_add_def_A_elem_fct(
		id, &T::template add_def_A_elem<TElem, TFVGeom>);

	this->set_add_def_M_elem_fct(
		id, &T::template add_def_M_elem<TElem, TFVGeom>);

	this->set_add_rhs_elem_fct(
		id, &T::template add_rhs_elem<TElem, TFVGeom>);
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
