/*
 * This file is part of the UG4 NavierStokes plugin.
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__RANS_TURBULENCE_FV1__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__RANS_TURBULENCE_FV1__

// ug4
#include "common/common.h"
#include "lib_grid/lg_base.h"

#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"

#include "../../upwind_interface.h"

namespace ug{
namespace NavierStokes{

/// FV1 discretization for two-equation RANS turbulence models.
template <typename TDomain>
class RANSTurbulenceFV1
	: public IElemDisc<TDomain>
{
protected:
	/// Base class type
	typedef IElemDisc<TDomain> base_type;

	/// Own type
	typedef RANSTurbulenceFV1<TDomain> this_type;

public:
	/// World dimension
	static const int dim = base_type::dim;

	/// Available RANS closures
	enum TurbulenceModel
	{
		K_OMEGA,
		K_EPSILON,
		K_OMEGA_SST
	};

public:
	/// Constructor
	RANSTurbulenceFV1(const char* functions, const char* subsets);
	
	/// Wall distance
	void set_wall_distance(SmartPtr<CplUserData<number, dim> > data);

	/// Set velocity field
	void set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > data);
	void set_velocity(const std::vector<number>& velocity);
	
	/// Set velocity gradient
	void set_velocity_gradient(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > data);
	void set_velocity_gradient(const MathMatrix<dim, dim>& velocityGradient);

	/// Set molecular kinematic viscosity
	void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data);
	void set_kinematic_viscosity(number viscosity);
	
	/// Set upwind method for turbulent transport
	void set_upwind(SmartPtr<INavierStokesUpwind<dim> > spUpwind)
	{
		m_spConvUpwind = spUpwind;
	}

	/// Set upwind method by name
	void set_upwind(const std::string& name)
	{
		m_spConvUpwind = CreateNavierStokesUpwind<dim>(name);
	}

	/// Select turbulence model
	void set_model(TurbulenceModel model);

	/// Return selected turbulence model
	TurbulenceModel model() const {return m_model;}

	/// Check finite-element setting
	virtual void prepare_setting(const std::vector<LFEID>& vLfeID,
								 bool bNonRegularGrid);

	/// Discretization identifier
	virtual std::string disc_type() const {return "fv1";}

public:
	/// Prepare loop over elements of one type
	template <typename TElem, typename TFVGeom>
	void prep_elem_loop(const ReferenceObjectID roid, const int si);

	/// Prepare one element
	template <typename TElem, typename TFVGeom>
	void prep_elem(const LocalVector& u,
				   GridObject* elem,
				   const ReferenceObjectID roid,
				   const MathVector<dim> vCornerCoords[]);

	/// Finish loop over elements of one type
	template <typename TElem, typename TFVGeom>
	void fsh_elem_loop();

	/// Stiffness Jacobian
	template <typename TElem, typename TFVGeom>
	void add_jac_A_elem(LocalMatrix& J,
						const LocalVector& u,
						GridObject* elem,
						const MathVector<dim> vCornerCoords[]);

	/// Mass Jacobian
	template <typename TElem, typename TFVGeom>
	void add_jac_M_elem(LocalMatrix& J,
						const LocalVector& u,
						GridObject* elem,
						const MathVector<dim> vCornerCoords[]);

	/// Stiffness defect
	template <typename TElem, typename TFVGeom>
	void add_def_A_elem(LocalVector& d,
						const LocalVector& u,
						GridObject* elem,
						const MathVector<dim> vCornerCoords[]);

	/// Mass defect
	template <typename TElem, typename TFVGeom>
	void add_def_M_elem(LocalVector& d,
						const LocalVector& u,
						GridObject* elem,
						const MathVector<dim> vCornerCoords[]);

	/// Right-hand side
	template <typename TElem, typename TFVGeom>
	void add_rhs_elem(LocalVector& d,
					  GridObject* elem,
					  const MathVector<dim> vCornerCoords[]);

protected:
	void register_all_funcs(bool bHang);

	template <typename TElem, typename TFVGeom>
	void register_func();

	number strain_rate_magnitude( const MathMatrix<dim, dim>& gradU) const;

	number vorticity_magnitude( const MathMatrix<dim, dim>& gradU) const;
	
	number turbulent_kinematic_viscosity(number k,number omega,number limiterMag,number F2) const;
	
	number blending_function_F1(number k, number omega, number nu, number wallDistance, number CDkw) const;
	
	number blending_function_F2(number k, number omega, number nu, number wallDistance) const;
	
	number blend_sst_coefficient(number F1, number innerValue, number outerValue) const;
	
	number cross_diffusion_CD(number omega, const MathVector<dim>& gradK, const MathVector<dim>& gradOmega) const;
	

protected:
	/// Velocity field u
	DataImport<MathVector<dim>, dim> m_imVelocity;
	
	/// Velocity grad u
	DataImport<MathMatrix<dim, dim>, dim> m_imVelocityGradientSCVF;
	DataImport<MathMatrix<dim, dim>, dim> m_imVelocityGradientSCV;

	/// Molecular kinematic viscosity nu
	DataImport<number, dim> m_imKinViscositySCVF;
	DataImport<number, dim> m_imKinViscositySCV;
	
	/// WallDistance
	DataImport<number, dim> m_imWallDistanceSCVF;
	DataImport<number, dim> m_imWallDistanceSCV;
	
	/// Upwind method for convection of k and omega
	SmartPtr<INavierStokesUpwind<dim> > m_spConvUpwind;
	
	


	/// Selected turbulence model
	TurbulenceModel m_model;

private:
	/// Local function index of k
	static const size_t _K_ = 0;

	/// Local function index of omega
	static const size_t _OMEGA_ = 1;
};

} // namespace NavierStokes
} // namespace ug

#endif
