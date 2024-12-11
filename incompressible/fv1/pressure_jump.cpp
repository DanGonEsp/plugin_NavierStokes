/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include <string>
#include <locale>

#include "pressure_jump.h"

#include "common/math/math_vector_matrix/math_vector_functions.h"
#include "common/math/math_vector_matrix/math_matrix_functions.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

template <int dim>
SmartPtr<INavierStokesPressureJump<dim> > CreateNavierStokesPressureJump(const std::string& name)
{
    std::string n = TrimString(name);
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);

    if(n == "viscous") return SmartPtr<NavierStokesViscousPressureJump<dim> >(new NavierStokesViscousPressureJump<dim>());


    UG_THROW("NavierStokes: Pressure Jump type '"<<name<<"' not a valid name"
             " Options are: viscous");
}

/////////////////////////////////////////////////////////////////////////////
// Interface for Stabilization
/////////////////////////////////////////////////////////////////////////////

//    register a update function for a Geometry
template <int dim>

template <typename TFVGeom, typename TAssFunc>
void
INavierStokesPressureJump<dim>::
register_update_func(TAssFunc func)
{
//    get unique geometry id
    size_t id = GetUniqueFVGeomID<TFVGeom>();

//    make sure that there is enough space
    if((size_t)id >= m_vUpdateFunc.size())
        m_vUpdateFunc.resize(id+1, NULL);

//    set pointer
    m_vUpdateFunc[id] = (UpdateFunc)func;
}




/////////////////////////////////////////////////////////////////////////////
// Pressure Jump
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesViscousPressureJump<TDim>::
update(const FV1Geometry<TElem, dim>* geo,
       const LocalVector& vCornerValue,
       const DataImport<MathVector<dim>, dim>& n,
       const DataImport<number, dim>& kinViscoSCV,
       const DataImport<number, dim>& density,
       const DataImport<number, dim>& densitySCV,
       const DataImport<number, dim>& jump_shape,
       const DataImport<number, dim>& vol_fraction,
       const number mu_l,
       const number rho_l,
       const number mu_g,
       const number rho_g,
       const number interface_value)
{
//    abbreviation for pressure
    static const size_t _P_ = dim;

    static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
//    size of the system
    static const size_t N = numSh;
    
    
    
    MathVector<dim> x, DX;
    MathVector<dim> xRho_interface[numSh];
    number interN[numSh];
    size_t NumSCVF = geo->num_scvf();

    
    number theta_to, theta_from, c_to, c_from, DC;
    number rho;

    for(size_t ip = 0; ip < N; ++ip)
    {
        VecSet(xRho_interface[ip],0.0);
        interN[ip]=0.0;

    }
    for(size_t ip = 0; ip < NumSCVF; ++ip)
    {
        const typename FV1Geometry<TElem, dim>::SCVF scvf = geo->scvf(ip);

        const size_t from=scvf.from();
        const size_t to=scvf.to();
        
        if (jump_shape[from]*jump_shape[to] < 0.0)
        {

            c_from = vol_fraction[from];
            c_to = vol_fraction[to];

            
            DC=c_to-c_from;
            VecSubtract(DX,geo->scv_global_ips()[to],geo->scv_global_ips()[from]);
            
            theta_to=  (c_to   - interface_value)/DC;
            theta_from=(c_from - interface_value)/DC;
            
            VecScaleAppend(xRho_interface[to], theta_to,   DX);
            VecScaleAppend(xRho_interface[from],   theta_from, DX);
            
            interN[from] += 1.0;
            interN[to] += 1.0;
        }
    }
//    a fixed size matrix
    DenseMatrix< FixedArray2<number, N, N> > mat;
//    reset all values of the matrix to zero
    mat = 0.0;
    for(size_t ip = 0; ip < N; ++ip)
    {
        
        rho = (jump_shape[ip] > 0.0)? rho_g:rho_l;
        VecScale(xRho_interface[ip], xRho_interface[ip], (rho_l-rho_g) / (interN[ip] * rho) );
        
        
        mat(ip, ip) += 1.0;
        
        for(size_t ip2 = 0; ip2 < N; ++ip2)
        {
            const typename FV1Geometry<TElem, dim>::SCV& scv = geo->scv(ip2);
            if (jump_shape[ip] * jump_shape[ip2] > 0.0)
                mat(ip, ip2) += jump_shape[ip2] * VecProd(scv.global_grad(ip2),xRho_interface[ip]);
        }
        
    }



//    we now create a matrix, where we store the inverse matrix
    typename block_traits<DenseMatrix< FixedArray2<number, N, N> > >::inverse_type inv;
    
    if(!GetInverse(inv, mat))
        UG_THROW("Could not compute inverse.");
    
    
    

    DenseVector< FixedArray1<number, N> > SumInv;
    DenseVector< FixedArray1<number, N> > SumInv2;
    SumInv = 1.0;
    
    
    DenseVector< FixedArray1<number, N> > rhs;
    rhs = 0.0;
    MatMult(SumInv2, 1.0, inv, SumInv);    //// Remember to change this vartiable to sum the contribution of all ips in the actual one
    
    //MathMatrix<dim,dim> VelGrad[N];
    for(size_t ip = 0; ip < N; ++ip)
    {
        
        for(size_t k = 0; k < numSh; ++k)
        {
            const typename FV1Geometry<TElem, dim>::SCV& scv = geo->scv(k);

            for(size_t d1 = 0; d1 < dim; ++d1)
            
            {
                
                number sumVel = 2 * (mu_l-mu_g) * n[ip][d1] * VecProd(scv.global_grad(k), n[ip] ) * SumInv2[ip] ;
                
                shape_vel(ip, d1, k) = sumVel;

                
                rhs[ip] += sumVel * vCornerValue(d1, k);
                
                
            }
            
            
            
            number sumP = VecProd(xRho_interface[ip],scv.global_grad(k)) * SumInv2[ip];
            
            shape_p(ip, k) = sumP ;

            
            rhs[ip] += sumP * vCornerValue(_P_, k);
            
            
            
            
        }
    }
    bool f= true;
    for(size_t ip = 0; ip < N; ++ip)
    {
    
        if (!((geo->scv_global_ips()[ip][0] > 9.749) && (geo->scv_global_ips()[ip][0] < 10.001)))
        {
            f = f && false;
        }
    }
    
    DenseVector< FixedArray1<number, N> > P_jump;

    MatMult(P_jump, 1.0, inv, rhs);
    
    for(size_t ip = 0; ip < N; ++ip)
    {
        pressure_jump(ip) = P_jump[ip];
        if (f)
        {
            //const typename FV1Geometry<TElem, dim>::SCV& scv = geo->scv(ip);
            if(ip==0) printf("Pressure jump at model\n" );
            printf("Pressure[%zu] = %f\n",ip, P_jump[ip]);
            printf("rhs[%zu] = %f\n",ip, rhs[ip]);
            printf("SumInv2[%zu] = %f\n",ip, SumInv2[ip]);
            //printf("n[%zu] = %f\n",ip, n[ip][0]);
            //printf("n[%zu] = %f\n",ip, n[ip][1]);
            
            //printf("grad[%zu] = %f\n",ip, scv.global_grad(ip)[0]);
            //printf("grad[%zu] = %f\n",ip, scv.global_grad(ip)[1]);
            
            //printf("u[%zu] = %f\n",ip, vCornerValue(0, ip));
            //printf("v[%zu] = %f\n",ip, vCornerValue(1, ip));
            
            //printf(" Velocity grad[%zu]\n", ip);
            //printf("   %f     %f    \n", VelGrad[ip][0][0],VelGrad[ip][0][1]);
            //printf("   %f     %f    \n", VelGrad[ip][1][0],VelGrad[ip][1][1]);
            


            
            
        }
    }
}



template <>
void NavierStokesViscousPressureJump<1>::register_func()
{
    register_func<RegularEdge>();
}

template <>
void NavierStokesViscousPressureJump<2>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
}

template <>
void NavierStokesViscousPressureJump<3>::register_func()
{
    register_func<RegularEdge>();
    register_func<Triangle>();
    register_func<Quadrilateral>();
    register_func<Tetrahedron>();
    register_func<Pyramid>();
    register_func<Prism>();
    register_func<Hexahedron>();
}


////////////////////////////////////////////////////////////////////////////////
//    explicit instantiations
////////////////////////////////////////////////////////////////////////////////

/*#ifdef UG_DIM_1
template class INavierStokesPressureJump<1>;
template class NavierStokesViscousPressureJump<1>;
template class NavierStokesFLOWStabilization<1>;

template SmartPtr<INavierStokesPressureJump<1> >CreateNavierStokesPressureJump<1>(const std::string& name);
#endif*/
#ifdef UG_DIM_2
template class INavierStokesPressureJump<2>;
template class NavierStokesViscousPressureJump<2>;


template SmartPtr<INavierStokesPressureJump<2> >CreateNavierStokesPressureJump<2>(const std::string& name);
#endif
#ifdef UG_DIM_3
template class INavierStokesPressureJump<3>;
template class NavierStokesViscousPressureJump<3>;


template SmartPtr<INavierStokesPressureJump<3> >CreateNavierStokesPressureJump<3>(const std::string& name);
#endif

} // namespace NavierStokes
} // end namespace ug
