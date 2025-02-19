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
#include "diffusion_length.h"
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
// Common functions for the Schneider-Raw-type Stabilizations
/////////////////////////////////////////////////////////////////////////////

template <int dim>
void
INavierStokesPressureJump<dim>::
set_diffusion_length(std::string diffLength)
{
    std::string n = TrimString(diffLength);
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    
    if      (n == "raw")        m_diffLengthType = RAW;
    else if (n == "fivepoint")  m_diffLengthType = FIVEPOINT;
    else if (n == "cor")        m_diffLengthType = COR;
    else
        UG_THROW("Diffusion Length calculation method not found."
                 " Use one of [Raw, Fivepoint, Cor].");
}
template <int dim>
template <typename TFVGeom>
void
INavierStokesPressureJump<dim>::
compute_diff_length(const TFVGeom& geo)
{
    //     Compute Diffusion Length in corresponding IPs
    switch(m_diffLengthType)
    {
        case FIVEPOINT: NSDiffLengthFivePoint(m_vDiffLengthSqInv, geo); return;
        case RAW:       NSDiffLengthRaw(m_vDiffLengthSqInv, geo); return;
        case COR:       NSDiffLengthCor(m_vDiffLengthSqInv, geo); return;
        default: UG_THROW(" Diffusion Length type not found.");
    }
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
       const MathVector<dim> vStdVel[],
       const bool bStokes,
       const DataImport<MathVector<dim>, dim>& n,
       const DataImport<number, dim>& kinViscoSCV,
       const DataImport<number, dim>& density,
       const DataImport<number, dim>& densitySCV,
       const DataImport<number, dim>& jump_shape,
       const DataImport<number, dim>& vol_fraction,
       const DataImport<MathVector<dim>, dim>* pSource,
       const LocalVector* pvCornerValueOldTime, number dt,
       const number density_ref,
       const number mu_l,
       const number rho_l,
       const number mu_g,
       const number rho_g,
       const number interface_value)
{
    
//     Only first order implemented
    UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");
    
    if( non_zero_shape_ip())
    {
        UG_THROW("Not implemented for ip velocities depending on other ip.");
    }
//    abbreviation for pressure
    static const size_t _P_ = dim;

    static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
//    size of the system
    static const size_t N = numSh;
    
    //    compute upwind (no convective terms for the Stokes eq. => no upwind)
    if (! bStokes) this->compute_upwind(geo, vStdVel);
    
    //    compute diffusion length
    this->compute_diff_length(*geo);
    
    
    for(size_t ip = 0; ip < N; ++ip)
    {
        
        pressure_jump(ip) = 0.0;
        for(size_t k = 0; k < numSh; ++k)
        {
            pressure_shape_p(ip, k) = 0.0 ;
            for(size_t d1 =0; d1 < dim; ++d1)
                pressure_shape_vel(ip, d1, k) = 0.0;

        }
        
        
    }
    
    
    /////////////////////////////////////////////////////////////////////////////
    // Calculation X_interface
    /////////////////////////////////////////////////////////////////////////////
    MathVector<dim> x, DX;
    MathVector<dim> x_interface[numSh];
    number interN[numSh];
    size_t NumSCVF = geo->num_scvf();
    
    number theta_to, theta_from, c_to, c_from, DC;
    number rho;
    
    
    number Inv_DiffLenSq = 0.0;
    number vNormStdVelPerConvLen = 0.0;
    number count_interface=0;

    for(size_t ip = 0; ip < N; ++ip)
    {
        VecSet(x_interface[ip],0.0);
        interN[ip]=0.0;

    }
    for(size_t ip = 0; ip < NumSCVF; ++ip)
    {
        const typename FV1Geometry<TElem, dim>::SCVF scvf = geo->scvf(ip);

        const size_t from=scvf.from();
        const size_t to=scvf.to();
        
        if (jump_shape[from]*jump_shape[to] < 0.0)
        {
            diff_length_sq_inv(ip);
            c_from = vol_fraction[from];
            c_to = vol_fraction[to];

            
            DC=c_to-c_from;
            VecSubtract(DX,geo->scv_global_ips()[to],geo->scv_global_ips()[from]);
            
            theta_to=  (c_to   - interface_value)/DC;
            theta_from=(c_from - interface_value)/DC;
            
            VecScaleAppend(x_interface[to], theta_to,   DX);
            VecScaleAppend(x_interface[from],   theta_from, DX);
            
            interN[from] += 1.0;
            interN[to] += 1.0;
            
            
            
            Inv_DiffLenSq += diff_length_sq_inv(ip);
            count_interface += 1.0;
            
            
            if(!bStokes)
            {
                const number norm = VecTwoNorm(vStdVel[ip]);
                vNormStdVelPerConvLen += norm / upwind_conv_length(ip);
                //vNormStdVelPerDownLen[ip] = norm / (downwind_conv_length(ip) + upwind_conv_length(ip));
                
            }
        }
    }
    for(size_t ip = 0; ip < N; ++ip)
    {
        VecScale(x_interface[ip], x_interface[ip], 1.0 / interN[ip]  );

    }
    Inv_DiffLenSq *= 1.0/count_interface;
    
    if(!bStokes)
        vNormStdVelPerConvLen *= 1.0/count_interface;
    
    number diag2 = Inv_DiffLenSq * mu_l;
    number diag1 = Inv_DiffLenSq * mu_g;
    
    /*if(pvCornerValueOldTime != NULL)
    {
        diag2 += 1.0/dt;
        diag1 += 1.0/dt;
    }
    if(!bStokes)
    {
        diag2 += vNormStdVelPerConvLen;
        diag1 += vNormStdVelPerConvLen;
    }*/
    
    /////////////////////////////////////////////////////////////////////////////
    // SlipVelocity
    /////////////////////////////////////////////////////////////////////////////
    
    MathVector<dim> Tang, normal_vel, VelVel;
    VecSet(Tang,0.0);
    number VOL_t = 0.0;
    
    MathVector<dim> C_grad;
    VecSet(C_grad,0.0);
    
    for(size_t ip = 0; ip < N; ++ip)
    {
        const typename FV1Geometry<TElem, dim>::SCV& scv = geo->scv(ip);
        const number vol = scv.volume();
        VOL_t += vol;
        
        for(size_t d1 = 0; d1 < dim; ++d1)
        {
            Tang[d1] += vCornerValue(d1, ip) * vol;
            C_grad[d1] += scv.global_grad(ip)[d1] * vol_fraction[ip];
        }
    }
    VecScale(Tang, Tang, 1.0 / VOL_t);
    number C_grad_magnitud = sqrt(VecProd(C_grad,C_grad));
    
    VelVel= Tang;
    VecScale(normal_vel, n[0], VecProd(n[0],Tang));
    VecSubtract(Tang,Tang,normal_vel);
    number Tang_mag= sqrt(VecProd(Tang,Tang));
    if( Tang_mag < 1e-08)
    {
        VecSet(Tang,0.0);
        Tang[0] = 1.0;
        VecScale(normal_vel, n[0], VecProd(n[0],Tang));
        VecSubtract(Tang,Tang,normal_vel);
        VecScale(Tang, Tang, 1.0 / sqrt(VecProd(Tang,Tang)));
        
        
    }
    else
        VecScale(Tang, Tang, 1.0 / Tang_mag);
    
    
    number Visc_eff = mu_l;
    const typename FV1Geometry<TElem, dim>::SCV& scv = geo->scv(0);
    
    for(size_t sh = 0; sh < numSh; ++sh)
    {
        number Nl = ( jump_shape[sh]>0) ? 1.0 : 0.0;
        Visc_eff -= (mu_l - mu_g) * (vol_fraction[sh] - interface_value) * Nl * VecProd(scv.global_grad(sh), n[0]) / C_grad_magnitud;
        
    }
    for(size_t k = 0; k < numSh; ++k)
        VecSet(tang_vel(k),0.0);
    

    
    for(size_t k = 0; k < numSh; ++k)
    {
        
        for(size_t d1 =0; d1 < dim; ++d1)
        {
            number Deriv = 0.0;
            for(size_t d2 =0; d2 < dim; ++d2)
            {
                Deriv +=  (Tang[d1]*n[0][d2] +Tang[d2]*n[0][d1]) * scv.global_grad(k)[d2];
            }
            
            for(size_t d =0; d < dim; ++d)
            {
                for(size_t ip = 0; ip < numSh; ++ip)
                {
                    tang_vel_shape_vel(ip, d, d1, k) = -VecProd(n[ip],x_interface[ip]) * (mu_l - mu_g) * Deriv * Tang[d] / Visc_eff;
                    tang_vel(ip)[d] += tang_vel_shape_vel(ip, d, d1, k) * vCornerValue(d1, k);
                }
                
            }
            
        }
    }
    
    
    
    
    /////////////////////////////////////////////////////////////////////////////
    // Pressure Jump
    /////////////////////////////////////////////////////////////////////////////

//    a fixed size matrix
    DenseMatrix< FixedArray2<number, N, N> > mat;
//    reset all values of the matrix to zero
    mat = 0.0;
    for(size_t ip = 0; ip < N; ++ip)
    {
        
        //rho = (jump_shape[ip] > 0.0)? rho_g:rho_l;
        //VecScale(x_interface[ip], x_interface[ip], 1.0 );
        
        
        mat(ip, ip) += 1.0;
        
        for(size_t ip2 = 0; ip2 < N; ++ip2)
        {
            const typename FV1Geometry<TElem, dim>::SCV& scv = geo->scv(ip2);
            const number Ng = (jump_shape[ip]>0 && jump_shape[ip2] <0.0)? 1.0 : 0.0;
            const number Nl = (jump_shape[ip]<0 && jump_shape[ip2] >0.0)? 1.0 : 0.0;

            mat(ip, ip2) += -VecProd(scv.global_grad(ip2),x_interface[ip])*(Ng*mu_l/rho_l+Nl*mu_g/rho_g)/(mu_l/rho_l+mu_g/rho_g) ;
            
        }
        
    }
    /*printf("----------------------------------------------------------- Velocity grad[]\n");
    printf("-----------------------------------------------------------   %f     %f    %f\n", mat(0,0),mat(0,1),mat(0,2));
    printf("-----------------------------------------------------------   %f     %f    %f\n", mat(1,0),mat(1,1),mat(1,2));
    printf("-----------------------------------------------------------   %f     %f    %f\n", mat(2,0),mat(2,1),mat(2,2));

    
    printf("x    y     =    %f     %f    \n", geo->scv_global_ips()[0][0],geo->scv_global_ips()[0][1]);
    printf("x    y     =    %f     %f    \n", geo->scv_global_ips()[1][0],geo->scv_global_ips()[1][1]);
    printf("x    y     =    %f     %f    \n", geo->scv_global_ips()[2][0],geo->scv_global_ips()[2][1]);*/

//    we now create a matrix, where we store the inverse matrix
    typename block_traits<DenseMatrix< FixedArray2<number, N, N> > >::inverse_type inv;
    
    if(!GetInverse(inv, mat))
        UG_THROW("Could not compute inverse.");
    
    
    

    DenseVector< FixedArray1<number, N> > SumInv2;
    DenseVector< FixedArray1<number, N> > SumInv;
    SumInv2 = 1.0;

    DenseVector< FixedArray1<number, N> > rhs;
    rhs = 0.0;
    DenseVector< FixedArray1<number, N> > P_jump;
    P_jump = 0.0;
    MatMult(SumInv, 1.0, inv, SumInv2);    //// Remember to change this vartiable to sum the contribution of all ips in the actual one
    

    //MathMatrix<dim,dim> VelGrad[N];
    for(size_t ip = 0; ip < N; ++ip)
    {
        
        for(size_t k = 0; k < numSh; ++k)
        {
            const typename FV1Geometry<TElem, dim>::SCV& scv = geo->scv(k);

            for(size_t d1 = 0; d1 < dim; ++d1)
            
            {
                
                number sumVel =  2.0*(mu_l-mu_g) * n[ip][d1] * VecProd(scv.global_grad(k), n[ip] ) * SumInv[ip] ;
                
                pressure_shape_vel(ip, d1, k) = sumVel;

                
                P_jump[ip] += sumVel * vCornerValue(d1, k);
                
                
            }
            
            
            number sumP = 0.0;
            for(size_t ip2 = 0; ip2 < N; ++ip2)
            {
                sumP += mat(ip, ip2) * VecProd(x_interface[ip2],scv.global_grad(k)) * (mu_l/rho_l-mu_g/rho_g)  /(mu_l/rho_l+mu_g/rho_g) ;
            }
            
            pressure_shape_p(ip, k) = sumP ;

            
            P_jump[ip] += sumP * vCornerValue(_P_, k);
            
            
            
            
        }
    }
    bool f= true;
    number yy = 0.0;
    for(size_t ip = 0; ip < N; ++ip)
    {
    
        if (!((geo->scv_global_ips()[ip][0] > 11.99) && (geo->scv_global_ips()[ip][0] < 12.255)))
        {
            f = f && false;
            
        }
        yy += geo->scv_global_ips()[ip][1];
    }
    
    

    //MatMult(P_jump, 1.0, inv, rhs);
    
    for(size_t ip = 0; ip < N; ++ip)
    {
        pressure_jump(ip) = P_jump[ip];
        if (false)
        {
            //const typename FV1Geometry<TElem, dim>::SCV& scv = geo->scv(ip);
            if(ip==0) printf("Pressure jump at model\n" );
            printf("yy = %f\n", yy/N);
            printf("Pressure[%zu] = %f\n",ip, P_jump[ip]);

            printf("SumInv2[%zu] = %f\n",ip, SumInv[ip]);
            
            printf("Visc Effec = %f\n",Visc_eff);
            printf("mu_l  = %f\n",mu_l);
            printf("mu_g  = %f\n",mu_g);
            //printf("n[%zu] = %f\n",ip, n[ip][0]);
            //printf("n[%zu] = %f\n",ip, n[ip][1]);
            
            //printf("grad[%zu] = %f\n",ip, scv.global_grad(ip)[0]);
            //printf("grad[%zu] = %f\n",ip, scv.global_grad(ip)[1]);
            
            printf("Velocity  \n");
            
            printf("ut = %f\n", Tang[0]);
            printf("vt = %f\n", Tang[1]);
            
            //printf("u = %f\n", VEL_t[0]);
            //printf("v = %f\n", VEL_t[1]);
            
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
