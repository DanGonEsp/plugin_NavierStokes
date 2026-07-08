/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC_NV_CUT_ELEMENT_
#define __H__UG__LIB_DISC__SPATIAL_DISC_NV_CUT_ELEMENT_

// extern libraries
#include <cmath>
#include <map>
#include <vector>

// other ug4 modules
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"
#include "consistent_gravity_multiphase.h"

// library intern includes



namespace ug{

template <int TDim, int TWorldDim = TDim>
class Interface
{


	public:
		///	dimension of reference element
		static const int dim = TDim;

		///	dimension of world
		static const int worldDim = TWorldDim;

	public:
		/// construct object and initialize local values and sizes
		Interface():
			m_P0(0.0),
			rho_s(2400.0),
			rho_max(2400.0),
			rho_a(1.2),
			mu_a(1.48e-5),
			nu_s(1.48e-5),
			dp(0.001),
			Fr(0.05),
			B_phi(1.0),
			alpha_max(0.635),
			alpha_min(0.57),
			deltaGamma(1e-4),
			FricMu_1(0.38),
			FricMu_2(0.64),
			I_0(0.279),
			deltaI(1e-3),
			deltaPs(1.48e-04),
			m_limit(1e6),
			m_dt(1.0),
			drag_model(1),
			m_interface_value(0.5),
			m_packing_factor(0.6),
			m_gravitation(-9.81),
			m_RelVelError(2e-01),
			epsilon(1e-012),
			m_delta(0.05),
			m_bParticleGradientForce(false),
			m_bConsistentGravity(false),
			m_init(false)
		{
				
		}

		void vPs(number* ParticlePressure, number* DCParticlePressure, const number Gamma[], const LocalVector& u, const size_t _C_, const size_t numSH, const bool deriv)
		{
			bool DERIVATIVE = false;
			if(deriv || DCParticlePressure != NULL)
				DERIVATIVE = true;
				
			for(size_t sh = 0; sh < numSH; ++sh)
			{
				number gamma=Gamma[sh];
				//const number phi=fmin(1.0, fmax(u(_C_,sh),0));
				const number vol=u(_C_,sh);

				
				number Ps_val, DPs;
				Ps( Ps_val,  DPs, gamma, vol, DERIVATIVE);
				//Dynamic pressure
				ParticlePressure[sh] = Ps_val;
				if(DERIVATIVE)
					DCParticlePressure[sh] = DPs;
			}

		}
		void Ps(number& ParticlePressure, number& DCParticlePressure, const number gamma, const number vol, const bool DERIVATIVE)
		{
			const number delta = 1e-02;
			const number Co=alpha_max-delta;
			const number phi = m_packing_factor * vol;
			//Stokes number
			const number St=gamma*rho_s*pow(dp,2)/mu_a;
			
			//Permanent contact pressure
			number pff= 0.0;
			number dpff = 0.0;
			number dpa = 0.0;
			number pa = 0.0;
			if(phi >= Co)
			{
				const number p0 = Fr *pow(  Co-alpha_min,3) /pow(alpha_max-Co,5);
				const number dp0 = Fr * pow( Co-alpha_min,2) * (2.0*Co+3.0*alpha_max-5.0*alpha_min)/pow(alpha_max-Co,6.0);
				pff = dp0*(phi-Co)+p0;
				pa = 0.0;//*mu_a*(1.0+St)*pow(B_phi*phi/(alpha_max-phi),2)*gamma;
				if(DERIVATIVE)
				{
					dpff = dp0;
					//dpa = 2*mu_a*(1.0+St)*pow(B_phi,2.0)*phi*alpha_max*gamma/(pow(alpha_max-phi,3));
				}
			}
			else if(phi >= alpha_min)
			{
				pff = Fr *pow(  phi-alpha_min,3) /pow(alpha_max-phi,5);
				pa = 0.0;//*mu_a*(1.0+St)*pow(B_phi*phi/(alpha_max-phi),2)*gamma;
				if (DERIVATIVE)
				{
					dpff =  Fr * pow( phi-alpha_min,2) * (2.0*phi+3.0*alpha_max-5.0*alpha_min)/pow(alpha_max-phi,6.0);
					//dpa = 2*mu_a*(1.0+St)*pow(B_phi,2.0)*phi*alpha_max*gamma/(pow(alpha_max-phi,3));
				}
			}
			
			//Dynamic pressure
			ParticlePressure = pff+pa;
			if(DERIVATIVE)
				DCParticlePressure = dpff + dpa;

			if(std::isnan(pff) || std::isnan(pa) || pff < 0.0 || pa<0.0 ) UG_THROW("Error in PropertiesInterface: Export particlePressure: Ps = " << ParticlePressure << "   phi = "<< phi<< "   Pa = "<< pa<< "   Pf = "<< pff);
		}
	
		void MU_I_Viscosity(number& mu_s, number& Dmu_s, const number  gamma_nr, const number vol, const bool deriv)
		{
			
			const number gamma = sqrt(pow(deltaGamma,2) + pow(gamma_nr,2.0));
			const number phi = m_packing_factor * vol;
			number Ps_val, DPs, s, mu_friction;
			Ps( Ps_val,  DPs, gamma, vol, deriv);
			
			number St=gamma*rho_s*pow(dp,2)/mu_a;
			number I=mu_a*(1+St)*gamma/(Ps_val+deltaPs)+deltaI;//grad_vel_mag*Viscosity_fluid/((Pressure_s+0.001)*ParticleDensity);
			//DI=-I/(Ps+deltaPs);
			
			//number I=gamma*dp/pow((Ps_val+deltaPs)/rho_s,0.5)+deltaI;
			//DI=-0.5*I/(Ps+deltaPs);
			
			mu_friction=FricMu_1+(FricMu_2-FricMu_1)/(1.0+I_0/I);
			
			mu_s=mu_friction*Ps_val/gamma;
			Dmu_s = mu_friction*DPs/gamma;
			
			if(std::isnan(mu_s) || mu_s<0.0 ) UG_THROW("Error in MU(I) Viscosity: Value = NaN" <<"  Volume Fraction = "<<phi<<".");
			//if(phi > alpha_max) UG_LOG("Phi > phi_max in Einstein Viscosity\n");
		}
		void Einstein_viscosity(number& ss, number& Dss,const number vol, const bool deriv)
		{
			
			const number phi = m_packing_factor * vol;
			number C_r = alpha_max - 1e-03;
			number power = 2.5;
			if (phi<=C_r)
			{
				ss= mu_a*pow(1.0-phi/alpha_max,-power*alpha_max);
				if(deriv)
					Dss =  power*ss/(1.0-phi/alpha_max) ;
			}
			else
			{
				
				number ss_r = mu_a*pow(1.0-C_r/alpha_max,-power*alpha_max);
				number slope=power*ss_r/(1.0-C_r/alpha_max);
				ss = slope*(phi-C_r)+ss_r;
				if(deriv)
					Dss =  slope ;
			}
			
			if(std::isnan(ss) || ss<0.0 ) UG_THROW("Error in Einstein ViscosityLinker: Value = NaN" <<"  Volume Fraction = "<<phi<<"  m_packing_factor =  " << m_packing_factor << ".");
			//if(phi > alpha_max) UG_LOG("Phi > phi_max in Einstein Viscosity\n");
		}

		number RelVel_ext(const number vol, const number rho_a1, const number dp1, const number rho_s1, const number g1)
		{
			number Vel = 6.8598478663758;
			size_t iter;
			RelVel(Vel, iter,  vol,  rho_a1, rho_a1, dp1,  rho_s1,  fabs(g1));
			UG_LOG("Sediment Velocity Ws = " << Vel << "\n");
			//UG_LOG("Cd = " << cd  << "\n");
			//UG_LOG("Iter = " << iter << "\n");
			return Vel;
			
		}
		void RelVel(number& Rel, size_t& iter, const number vol, const number rho_mix, const number rho_a, const number dp1, const number rho_s, const number g1)
		{
			
			number mu_mix, Dss;
			mu_mix = 1.776e-05;//Einstein_viscosity(mu_mix , Dss, 0.0 ,false);
			size_t mod = drag_model;
			size_t i=0;
			number e=10.0;
			number w2=Rel;
			number w1=Rel;
			number re, c;
			
			
			while(e>m_RelVelError)
			{
				
				w1=w2;
				re=RE(mu_mix,rho_mix,dp1,w1);
				c = CD(re,mod);
				w2 = sqrt((4.0/3.0)*dp1*(rho_s-rho_a)*g1/ (c*rho_mix));
				e = 100.0 * fabs(w2-w1)/w2;
				i=i+1;
				if(i > 100) UG_THROW("Error in RelativeVelocity: Reached " << i <<" iterations in RelVel. Vel = "  << w2 <<"   Error = "<<e<< " Tol =  "<< m_RelVelError << "  Mu =  " << mu_mix << "  Gz =  "<< g1 <<"\n");
			}
			iter = i;
			Rel = w2;

			if(std::isnan(Rel) || Rel<0.0) UG_THROW("Error in RelativeVelocity, values out of range: Vel  = "<< Rel << "  Mu = "<<mu_mix<<".");
		}
		number RE(const number mu_a1, const number rho_a1, const number dp1, const number w1)
		{
			return rho_a1*dp1*w1/mu_a1;
		}
		number CD(const number RE, const size_t mod)
		{
			number c;
			if(mod == 0)
				c = 24.0 / RE;
			else if (mod == 1)
				c=pow(0.63+4.8/sqrt(RE),2.0);
			else if (mod == 2)
			{//Schiller-Naumann
				if(RE>1000.0)
					c = 0.44;
				else
					c = (24.0/RE)*(1.0+0.15*pow(RE,0.687));
				
			}
			else if (mod == 3)
			{//Turton and Levenspiel
				c=(24.0/RE)*(1.0+0.173*pow(RE,0.657))+0.413/(1.0+16300.0*pow(RE,-1.09));
			}
			else UG_THROW("Error in RelativeVelocityFunction: Drag Coefficient model not defined.");
			return c;
		}
		void RelativeVelDeriv(MathVector<dim>& RelVel_deriv, number Ws, number  uVal)
		{
			const number dc = 1e-02;
			size_t iter = 0;
			
			number Mu, dMu;
			Einstein_viscosity( Mu, dMu, uVal+dc,false);
			number Ws_2 =Ws;
			
			RelVel(Ws_2, iter, Mu, rho_a, rho_a, dp, rho_s, fabs(m_gravitation));
			VecSet(RelVel_deriv,0.0);
			RelVel_deriv[dim-1] = -(Ws_2 - Ws)/dc;
			//UG_LOG("Iterations = "<<iter<<"\n")
			
		}
		void GodunovCharacteristic_ip(MathVector<dim>& vVel, MathVector<dim>& vStdVel, const number RhoL, const number RhoR, const number MuL, const number MuR, const MathVector<dim> vUL, const MathVector<dim> vUR, const MathVector<dim> vNormal)
		{
			number power = 1.0;
			number FL, SL, FR, SR;
			number area = VecLength(vNormal);
			MathVector<dim> normal; VecScale(normal,vNormal, 1.0/area);
			MathVector<dim> vVel_aux = 0.0;
			MathVector<dim> vStdVel_aux = 0.0;
			
			SL=VecProd(vUL,normal);
			SR=VecProd(vUR,normal);
			
			FL = RhoL*SL;
			FR = RhoR*SR;
			
			VecScaleAdd(vStdVel_aux, RhoR/(RhoL + RhoR), vUR, RhoL/(RhoL + RhoR), vUL );
			VecScaleAdd(vVel_aux, 0.5 , vUR, 0.5, vUL);
			number dRho = fabs(RhoR- RhoL);
			number theta = pow( dRho/(1.0+dRho), power);
			
			number S_val = 0.0;
			
			if(fabs(RhoR-RhoL)< 1e-08)
			{
				//vStdVel = vVel_aux;
				//vVel = vVel_aux;
			}
			else
			{
				S_val = RankineHugoniot(FR, FL, RhoR,  RhoL);

				if(SL * SR >=0.0)
				{

					vStdVel = vStdVel_aux;
					vVel = vStdVel_aux;
					
				}else if(SL > 0.0 && SR <0.0)
				{
					vStdVel = vStdVel_aux;
					vVel = vStdVel_aux;
				}else if(SL < 0.0 && SR > 0.0)
				{
					vStdVel = vStdVel_aux;
					vVel = vStdVel_aux;
					
				}else
				{
					
					
					 UG_LOG("Shock or Rarefaction S_val   = "<< S_val   <<".\n");
					 UG_LOG("Shock or Rarefaction VL   = " << VecLength(vUL)  <<"  VR   = "<< VecLength(vUR)   <<".\n");
					 UG_LOG("Shock or Rarefaction SL   = " << SL  <<"  SR   = "<< SR   <<".\n");
					 UG_LOG("Shock or Rarefaction RhoL = "<< RhoL <<"  RhoR = "<< RhoR <<".\n");
					 UG_LOG("Shock or Rarefaction FL   = " << FL  <<"  FR   = "<< FR   <<".\n");
					UG_THROW("Starting point.\n");
				}
					
	
			
					
				
			}
			
			
			//VecScaleAdd(vStdVel, 0.5, vUR, 0.5, vUL );
			
			/*if(SL*SR >= 0 )
				//S_ip = vStdVel;
			else
			{
				//number power = 1.0;
				//number theta = 0.0;//pow(fabs(MuR- MuL) /(MuR + MuL), power);
				//VecScaleAdd(Vel, 0.5 , vUR, 0.5, vUL);
				
				//VecScaleAdd(S_ip, 1.0 - theta, vStdVel, theta, Vel);
				
				
			}*/
			
			//VecScaleAdd(vVel,   1.0 - theta, vStdVel, theta, vStdVel_aux);
				
	
		}
	
	

	
		number RankineHugoniot(const number wFR, const number wFL, const number UR, const number UL)
		{
			return (wFR-wFL)/(UR-UL);
		}
		number RankineHugoniotCharac(const number uL, const number uR, const number Vn, const number Wn)
		{
			return Vn + (1.0-(uL+uR))*Wn;
		}
		
		
		bool Compute_F_star(number& F_star, const number wFL, const number wFR, const number wSL, const number wSR, const number uL, const number uR, const number wL, const number wR, const number V)
		{

			number w_star = 0.5 * (wL + wR);
			number u_star = 0.5*(1.0 + V/w_star);
			F_star = u_star * V + u_star * (1.0 - u_star) * w_star;
			
			number F_L = uL * V + uL * (1.0 - uL) * wL;
			number F_R = uR * V + uR * (1.0 - uR) * wR;
			
			if(!(u_star < fmax(uL,uR) && u_star > fmin(uL,uR)))
			{
				
				F_star = 0.0;
				UG_LOG("Non valid number rarafection sonic V = " <<V<<"  uL = "<< uL << "  uR = "<< uR<<".\n");
				UG_LOG("Non valid number rarafection sonic u_star = " <<u_star<<"  F_star = "<< F_star << "  F_L = "<< F_L<<"  F_R = "<< F_R <<".\n");
				UG_THROW("RareFaction");
				
			}
			
			//Rusanov_flux(F_star, uL, uR,  SL,  SR,  wFL,  wFR);
			//printf("rare");
			
		}
		bool Compute_F_star_jac(number& F_star_jacL, number& F_star_jacR, const number wFL, const number wFR, const number wSL, const number wSR, const number uL, const number uR, const number WL, const number WR, const number Vel_ip)
		{
			F_star_jacL = 0.0;
			F_star_jacR = 0.0;
			/*size_t iter = 0;

			number u_star = uL-SL*(uR-uL)/(SR-SL);
			number w_star = WL + (WR-WL) * (u_star - uL) / (uR-uL);
			
			F_star_jacL = u_star*Vel_ip + u_star * (1.0-u_star)*w_star;*/
			
			//Rusanov_jac(F_star_jacL,F_star_jacR, uL, uR,  SL,  SR,  wFL,  wFR);
			
		}

		void Flux_ip(number& Flux_ip, const number UL, const number UR, const MathVector<dim> vWL, const MathVector<dim> vWR,  const MathVector<dim> vVel_ip, const MathVector<dim> vNormal, const int upwind_vol_method)
		{
			
			number wFL, wSL, wFR, wSR, WL, WR, Vn, Flux;
			bool GeoDegL, GeoDegR;
			bool MagDegL, MagDegR;
			number area = VecLength(vNormal);
			MathVector<dim> normal; VecScale(normal,vNormal, 1.0/area);
			
			Flux_function_interface(wFL,wSL, UL, WL, vWL, GeoDegL, MagDegL,
									wFR,wSR, UR, WR, vWR, GeoDegR, MagDegR, Vn, vVel_ip, normal, true,true);
			
				
			switch (upwind_vol_method) {
				case 1:
					Godunov_flux(Flux,UL,UR,wSL,wSR,wFL,wFR, WL,WR,Vn);
					break;
				case 2:
					Rusanov_flux(Flux, UL,UR,wSL,wSR,wFL,wFR,Vn);
					break;
				case 3:
					Roe_flux(Flux, UL,UR,wSL,wSR,wFL,wFR, WL,WR,Vn);
					break;
				default:
					UG_THROW("Wrong model selected for solving riemman problem: PropInterface has options"
							 " model= 0, 1 , 2, 3");
					break;
			}
			if(std::isnan(Flux))
			{
				UG_THROW("Non valid number for FLUX = " <<Flux<<".\n");
				
			}
			
			
			
			Flux_ip = area * Flux;
			
			
		}
		void Flux_Jac_ip(number& JacVL, number& JacVR, number& JacWL, number& JacWR, const number UL, const number UR, const MathVector<dim> vWL, const MathVector<dim> vWR,  const MathVector<dim> vVel_ip, const MathVector<dim> vNormal, const int upwind_vol_method)
		{
			
			number wFL, wSL, wFR, wSR, WL, WR, Vn, Flux;
			bool GeoDegL, GeoDegR;
			bool MagDegL, MagDegR;
			number area = VecLength(vNormal);
			MathVector<dim> normal; VecScale(normal,vNormal, 1.0/area);
			
			Flux_function_interface(wFL,wSL, UL, WL, vWL, GeoDegL, MagDegL,
									wFR,wSR, UR, WR, vWR, GeoDegR, MagDegR, Vn, vVel_ip, normal, true,true);
			
				
			switch (upwind_vol_method) {
				case 1:
					Godunov_jac(JacVL,JacVR, JacWL,JacWR,  UL, UR, wSL, wSR, wFL, wFR,  WL, WR,  Vn);
					break;
				case 2:
					Rusanov_jac(JacVL,JacVR,JacWL,JacWR, UL, UR,  wSL,  wSR,  wFL,  wFR,  Vn);
					break;
				case 3:
					Roe_jac(JacVL,JacVR,JacWL,JacWR, UL, UR,  wSL,  wSR,  wFL,  wFR, WL, WR,  Vn);
					break;
				default:
					UG_THROW("Wrong model selected for solving riemman problem: PropInterface has options"
								   " model= 0, 1 , 2")
					break;
			}
			
			if(std::isnan(JacWL + JacWR + JacVL + JacVR))
			{
				
				UG_THROW("Non valid number for Jacobian in NonLinearTransportEquation JacWL = " <<JacWL<<"  JacWR = "<< JacWR << "  JacVL = "<< JacVL<<"  JacVR = "<< JacVR <<".\n");
				
			}
			
			JacWL = area * JacWL;
			JacWR = area * JacWR;
			JacVL = area * JacVL;
			JacVR = area * JacVR;
			
			
		}
		void Flux_function_interface( number& wFL, number& wSL, const number uL, number& WL, const MathVector<dim> vWL, bool& boolGeomDegL, bool& boolMagDegL,
									  number& wFR, number& wSR, const number uR, number& WR, const MathVector<dim> vWR, bool& boolGeomDegR, bool& boolMagDegR,
						    number& Vn, const MathVector<dim> Vel_ip, const MathVector<dim> normal, const bool boolFlux, const bool boolFluxDeriv )
		{
			bool boolSameU = false;
			
			MathVector<dim> vDW;
			VecSubtract(vDW, vWR,vWL);
			number du = uR-uL;
			boolSameU = true;
			vDW = 0.0;
			/*if(fabs(du)> epsilon)
			{
				VecScale(vDW,vDW,1.0/du);
				boolSameU = false;
			}*/
	
			
			Flux_function( uL, wFL,  wSL, WL,  vWL, Vel_ip, vDW, normal, boolGeomDegL, boolMagDegL,  boolFlux, boolFluxDeriv );
			Flux_function( uR, wFR,  wSR, WR,  vWR, Vel_ip, vDW, normal, boolGeomDegR, boolMagDegR,  boolFlux, boolFluxDeriv );
			Vn = VecProd(Vel_ip,normal);
			
			
		}
	
		void Flux_function( const number u, number& Flux, number& Flux_prime, number& Wn, const MathVector<dim> vW, const MathVector<dim> Vel_ip, const MathVector<dim> vDW, const MathVector<dim> normal,
						   bool& boolGeomDeg, bool& boolMagDeg, const bool boolFlux, const bool boolFluxDeriv )
		{
			boolMagDeg = false;
			boolGeomDeg = false;
			
			MathVector<dim> vFlux;
			MathVector<dim> vFlux_prime;
			
			const number func = u*(1.0-u);
			const number func_prime = 1.0-2.0*u;
			
			VecScale(vFlux, vW,func);
			VecScaleAdd(vFlux_prime, func_prime, vW, func, vDW);
			
			
			number Smag = VecLength(vFlux_prime);

			Flux = VecProd(vFlux, normal);
			Flux_prime = VecProd(vFlux_prime, normal);
			
			
			if (Smag < epsilon)
				boolMagDeg = true;
			if(Flux_prime/Smag < epsilon &&  Smag>epsilon)
				boolGeomDeg = true;
				
			Wn = VecProd(vW, normal);
			if(std::isnan(Wn))
				UG_LOG("  U = "<< u << "  normal = " <<normal[0]<<"  "<<normal[1]<< "  vW = "<< vW[0]<<"  "<< vW[1] <<".\n");
				
			
		}

		void Godunov_flux(number& Flux, const number UL, const number UR, const number wSL, const number wSR, const number wFL, const number wFR, const number WL, const number WR, const number Vn)
		{
			
			number fL = Vn*UL + wFL;
			number fR = Vn*UR + wFR;
			
			number w_star = 0.5 * (WL + WR);
			number u_star = 0.5*(1.0 + Vn/w_star);
			number F_star = u_star * Vn + u_star * (1.0 - u_star) * w_star;
			
			const number SL = Vn + wSL;
			const number SR = Vn + wSR;
			const number Sw = RankineHugoniotCharac( UL, UR,  Vn,  0.5 *(WL + WR));
			const number S_abs = m_delta * fmax(fabs(SL), fabs(SR));

			if (UL <= UR) {
				// Rarefaction/Increasing case: Look for a Minimum
				if (UL <= u_star && u_star <= UR) {
					// Sonic point is inside the range
					Flux =  fmin(fmin(fL, fR),F_star);
					Flux += -0.5*S_abs*(UR-UL);
				}else{
					Flux = (Sw > 0)? fL:fR;
				}
			} else {
				// Shock/Decreasing case: Look for a Maximum
				if (UR <= u_star && u_star <= UL) {
					// Sonic point is inside the range
					Flux =  fmax(fmax(fL, fR),F_star);
					Flux += -0.5*S_abs*(UR-UL);
				} else {
					Flux = (Sw > 0)? fL:fR;
				}
			}
			/*
			number Sw = 0.0;
			number SR = Vn + wSR;
			number SL = Vn + wSL;
			int ShockCase; bool shock_elem = false; bool rarefaction_elem = false;
			const number Smax = fmax(fabs(SR), fabs(SL));
			const number Smin = fmin(fabs(SR), fabs(SL));
			number GodunovFlux = 0.0;
			
			if ((SR * SL > 0.0))
			{
				
				ShockCase = 0;

				
				Sw = RankineHugoniotCharac( UL, UR,  Vn,  0.5 *(WL + WR));
				
				GodunovFlux = (Sw > 0)? UL*Vn + wFL:UR*Vn +wFR;
				
				
			}
			else if ( SL >0.0 && SR<0.0  || fabs(UR-UL)<)
			{
				

				shock_elem = true;
				ShockCase= 1;

				
				Sw = RankineHugoniotCharac( UL, UR,  Vn,  0.5 *(WL + WR));
				
				GodunovFlux = (Sw > 0)? UL*Vn + wFL:UR*Vn +wFR;
				
				
			}
			else if((SL <=0.0 && SR>=0.0 ) )
			{
				
				number F_star = 0.0;
				Compute_F_star( F_star,  wFL,  wFR, SL, SR, UL,  UR, WL,  WR,  Vn);

				ShockCase= 2;
				rarefaction_elem = true;
				
				//if(!entropy)
				//UG_LOG("Rarefaction error \n");
				
				GodunovFlux = F_star;
				
			}
			else
			{
				UG_THROW("Non valid state in Godunov  .\n");
			}
			
			Flux = GodunovFlux;// - 0.5 * 0.01 *Smax * (UR-UL);*/
			
		}
		void Godunov_jac(number& JacVL, number& JacVR, number& JacWL, number& JacWR, const number UL, const number UR, const number wSL, const number wSR, const number wFL, const number wFR, const number WL, const number WR, const number Vn)
		{
			
			
			const number SL = Vn + wSL;
			const number SR = Vn + wSR;
			number Sw = RankineHugoniotCharac( UL, UR,  Vn,  0.5 *(WL + WR));
			const number S_abs = m_delta * fmax(fabs(SL), fabs(SR));
			
			int ShockCase; bool shock_elem; bool rarefaction_elem = false;
			const number Smax = fmax(fabs(SR), fabs(SL));
			const number Smin = fmin(fabs(SR), fabs(SL));
			JacVL = 0.0;
			JacVR = 0.0;
			JacWL = 0.0;
			JacWR = 0.0;
			
			number eps_diss = 1e-02;
			
			
			
			
			
			number fL = Vn*UL + wFL;
			number fR = Vn*UR + wFR;
			
			number w_star = 0.5 * (WL + WR);
			number u_star = 0.5*(1.0 + Vn/w_star);
			number F_star = u_star * Vn + u_star * (1.0 - u_star) * w_star;

			if (UL <= UR) {
				// Rarefaction/Increasing case: Look for a Minimum
				if (UL <= u_star && u_star <= UR) {
					// Sonic point is inside the range
					if(fmin(fL, fR) <= F_star)
					{
						JacWL = (Sw > 0)?  wSL: 0.0;
						JacWR = (Sw > 0)? 0.0 : wSR;
						
						JacVL = (Sw > 0)?  Vn : 0.0;
						JacVR = (Sw > 0)? 0.0 : Vn;
						
						
					}
					
					
					JacWL +=   0.5*S_abs;
					JacWR += - 0.5*S_abs;
					
				}
				else{
					JacWL = (Sw > 0)?  wSL: 0.0;
					JacWR = (Sw > 0)? 0.0 : wSR;
					
					JacVL = (Sw > 0)?  Vn : 0.0;
					JacVR = (Sw > 0)? 0.0 : Vn;
				}
			} else {
				// Shock/Decreasing case: Look for a Maximum
				if (UR <= u_star && u_star <= UL)
				{
					// Sonic point is inside the range
					if(fmax(fL, fR) >= F_star)
					{
						JacWL = (Sw > 0)?  wSL: 0.0;
						JacWR = (Sw > 0)? 0.0 : wSR;
						
						JacVL = (Sw > 0)?  Vn : 0.0;
						JacVR = (Sw > 0)? 0.0 : Vn;
						
						
					}
					
					JacWL +=   0.5*S_abs;
					JacWR += - 0.5*S_abs;
						
				}
				else {
					
					JacWL = (Sw > 0)?  wSL: 0.0;
					JacWR = (Sw > 0)? 0.0 : wSR;
					
					JacVL = (Sw > 0)?  Vn : 0.0;
					JacVR = (Sw > 0)? 0.0 : Vn;
					
				}
			}
			
			/*if(fabs(UL-UR)< epsilon)
			{
				
			}
			
			if(Smin < epsilon || fabs(uR-uL)< epsilon || fabs(FR-FL)< epsilon)
			 {
			 //VecScale( StdCharacteristicVel[ip],scvf.normal(), Smax/ face_norm);
			 //ShockCase[ip] = 3;
			 //chararcter[ip] = false;
			 //continue;
			 number Vel_n = VecProd(Vel_ip[ip],scvf.normal())/ face_norm;
			 number Ws_n = 0.0;
			 
			 number WL = 0.0;
			 number WR = 0.0;
			 
			 if(RelVelSCV.data_given())
			 {
			 Ws_n += 0.5 * (VecProd(RelVelSCV[shR],scvf.normal()) + VecProd(RelVelSCV[shL],scvf.normal()))/ face_norm;
			 WL +=  VecProd(RelVelSCV[shL],scvf.normal()) / face_norm;
			 WR +=  VecProd(RelVelSCV[shR],scvf.normal()) / face_norm;
			 }
			 if(SlipVelSCV.data_given())
			 {
			 Ws_n += 0.5 * (VecProd(SlipVelSCV[shR],scvf.normal()) + VecProd(SlipVelSCV[shL],scvf.normal()))/ face_norm;
			 WL +=  VecProd(SlipVelSCV[shL],scvf.normal()) / face_norm;
			 WR +=  VecProd(SlipVelSCV[shR],scvf.normal()) / face_norm;
			 
			 }
			 
			 UG_LOG("Ip = "<<ip<<" FL = "<< FL<<" FR = "<<FR<<" SL = "<< SL<<" SR = "<<SR<<" uL = "<<uL<< " uR =  "<<uR<<" WL = "<<WL<< " WR =  "<<WR<< " Vel = "<< Vel_n<< " Ws = "<<Ws_n << " \n");
			 
			 }*/
			
			/*if(fabs(uR-uL)< eps )
			 {
			 VecScaleAdd(StdCharacteristicVel[ip],  0.5, LocalCharact_total_L, 0.5,  LocalCharact_total_R );
			 ShockCase[ip] = 0;
			 chararcter[ip] = false;
			 }*/
			
			
			/*if ((SR * SL > 0.0) )
			{
				

				ShockCase = 0;
				
				Sw = RankineHugoniotCharac( UL, UR,  Vn,  0.5 *(WL + WR));
				
				JacWL = (Sw > 0)?  wSL: 0.0;
				JacWR = (Sw > 0)? 0.0 : wSR;
				
				JacVL = (Sw > 0)?  Vn : 0.0;
				JacVR = (Sw > 0)? 0.0 : Vn;
				
				
			}
			else if ( (SL >0.0 && SR<0.0 ) )
			{
				
				shock_elem = true;
				ShockCase= 1;
				
				Sw = RankineHugoniotCharac( UL, UR,  Vn,  0.5 *(WL + WR));
				
				JacWL = (Sw > 0)?  wSL: 0.0;
				JacWR = (Sw > 0)? 0.0 : wSR;
				
				JacVL = (Sw > 0)?  Vn : 0.0;
				JacVR = (Sw > 0)? 0.0 : Vn;
				
				//Roe_jac(JacVL,JacVR,JacWL,JacWR, UL, UR,  wSL,  wSR,  wFL,  wFR, WL, WR,  Vn);
				
				
			}
			else if((SL <=0.0 && SR>=0.0 ) )
			{
				
				//number F_star_jacL = 0.0;
				//number F_star_jacR = 0.0;
				//Compute_F_star_jac( F_star_jacL, F_star_jacR,  wFL,  wFR, SL, SR, UL,  UR, WL,  WR,  Vn);

				ShockCase= 2;
				rarefaction_elem = true;
				
				//if(!entropy)
				//UG_LOG("Rarefaction error \n");
				
				number w_star = 0.5*(WL + WR);
				number u_star = 0.5*(UL + UR);
				
				number wS_star =  (1- 2.0*u_star)*w_star;
				
				JacWL = 0.0;//0.5*wS_star;
				JacWR = 0.0;//0.5*wS_star;
				JacVL = 0.0;// 0.5 * Vn;
				JacVR = 0.0;//0.5 * Vn;
				
				//Roe_jac(JacVL,JacVR,JacWL,JacWR, UL, UR,  wSL,  wSR,  wFL,  wFR, WL, WR,  Vn);
				
			}
			else
			{
				UG_THROW("Non valid state in Godunov Jac .\n");
			}*/

			//JacVL +=  +0.5 * eps_diss * Smax;
			//JacVR +=  -0.5 * eps_diss * Smax;
		
			
			
		}
		void Rusanov_flux(number& Flux, const number UL, const number UR, const number wSL, const number wSR, const number wFL, const number wFR,  const number Vn )
		{

			const number S=SmoothMax(fabs(Vn + wSL),fabs(Vn + wSR), 1e-06);
			//const number S=fmax(fabs(wSL),fabs(wSR));
			
			Flux = 0.5*(UL+UR)*Vn;
			//Flux = (Vn > 0)? UL*Vn:UR*Vn;
			
			Flux +=  0.5 * (wFL + wFR) - 0.5*S*(UR-UL);
			
		}
		void Rusanov_jac(number& JacVL,number& JacVR,number& JacWL,number& JacWR, const number UL, const number UR, const number wSL, const number wSR, const number wFL, const number wFR, const number Vn)
		{

			const number S=SmoothMax(fabs(Vn + wSL),fabs(Vn + wSR), 1e-06);
			//const number S=fmax(fabs(wSL),fabs(wSR));
			
			JacWL = 0.5*(wSL) + 0.5*S;
			JacWR = 0.5*(wSR) - 0.5*S;
			
			JacVL = 0.5 * Vn;
			JacVR = 0.5 * Vn;
			
			//JacVL = (Vn > 0)?  Vn : 0.0;
			//JacVR = (Vn > 0)? 0.0 : Vn;
			
			
		}
		void Roe_flux(number& Flux, const number UL, const number UR, const number wSL, const number wSR, const number wFL, const number wFR, const number WL, const number WR, const number Vn)
		{

			number S;
			const number SR =Vn + wSR;
			const number SL =Vn + wSL;
			
			
			S = RankineHugoniotCharac( UL, UR,  Vn,  0.5 *(WL + WR));
			//Entropy fixing Harten-Hyman
			number S_abs = Harten_Hyman( S, SL, SR);
			
			
			//Flux = (Vn > 0)? UL*Vn:UR*Vn;
			
			
			Flux = 0.5 * (UL + UR)*Vn;
			Flux +=   0.5 * (wFL + wFR) - 0.5*S_abs*(UR-UL);
			
		}
		void Roe_jac(number& JacVL,number& JacVR, number& JacWL,number& JacWR, const number UL, const number UR, const number wSL, const number wSR, const number wFL, const number wFR, const number WL, const number WR, const number Vn)
		{
			number S;
			const number SR =Vn + wSR;
			const number SL =Vn + wSL;
			
			
			S = RankineHugoniotCharac( UL, UR,  Vn,  0.5 *(WL + WR));
			//Entropy fixing Harten-Hyman
			number S_abs = Harten_Hyman( S, SL, SR);
			
			JacWL = 0.5*(wSL) + 0.5*S_abs;
			JacWR = 0.5*(wSR) - 0.5*S_abs;
			
			
			//JacVL = (Vn > 0)?  Vn : 0.0;
			//JacVR = (Vn > 0)? 0.0 : Vn;
			JacVL =  0.5 * Vn;
			JacVR =  0.5 * Vn;
			
			
		}
		//Entropy fixing Harten-Hyman
		number Harten_Hyman(const number S, const number SL, const number SR)
		{
			// 1. Global Background Diffusion (The "Epsilon")
			// This ensures the matrix is never singular.
			//const number S_max = SmoothMax(fabs(SL), fabs(SR),1e-06);
			//const number epsilon = m_delta * S_max + 1e-12; // Small floor to avoid div by zero
	

			// 3. Rarefaction (Entropy) Fix
			/*if(SL < 0.0 && SR > 0.0)
			{
				// Use the physical speeds to determine the fan width
				number delta = SmoothMax(SmoothMax(0.0, S - SL,1e-06), SR - S, 1e-06);
				
				// Ensure delta is at least as large as our global epsilon
				delta = SmoothMax(delta, epsilon, 1e-06);

				if (fabs(S) < delta) {
					S_abs = (S*S + delta*delta) / (2.0 * delta);
				}
			}*/
			
			// Use the physical speeds to determine the fan width
			number epsilon = SmoothMax(SmoothMax(0.0, S - SL,1e-06), SR - S, 1e-06);
			
			// Ensure delta is at least as large as our global epsilon
			//delta = SmoothMax(delta, epsilon, 1e-06);
			
			// 2. Initial smoothed absolute value
			number S_abs = (fabs(S) <= epsilon) ? (S*S + epsilon*epsilon)/(2.0*epsilon) : fabs(S);
			
			return S_abs;
		}
	
		//Entropy fixing Harten-Hyman
		/*number Harten_Hyman(const number S, const number SL, const number SR)
		{

			number S_abs;
		
			
			// Use the physical speeds to determine the fan width
			number delta = fmax(fmax(0.0, S - SL), SR - S);
			S_abs = fabs(S);
			if (S_abs <= delta) {
				S_abs = (S*S + delta*delta) / (2.0 * delta);
			}


			
			return S_abs;
		}*/
	 

		/*
		 if(Inter->boolConsistenGravity())
		 {
			 StdLinConsistentGravity<refDim> RhoG;
			 MathVector<refDim> vConsGravity[numSh];
			 MathVector<dim> Gravity; VecSet(Gravity,0.0); Gravity[dim-1]=9.81;
			 RhoG.template prepare<dim>(vConsGravity, numSh, vCornerCoords, m_imDensitySCV.values(), Gravity);
			 
			 for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
			 {
				 //     get current SCVF
				 const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
				 
				 RhoG.template compute<dim>(vConsGravitySCVF[ip],scvf.local_ip(), scvf.JTInv(),scvf.local_grad_vector(), vConsGravity);
			 }
			 
			 for(size_t ip = 0; ip < geo.num_scv(); ++ip)
			 {
				 //     get current SCV
				 const typename TFVGeom::SCV& scv = geo.scv(ip);
				 
				 RhoG.template compute<dim>(vConsGravitySCV[ip],scv.local_ip(), scv.JTInv(),scv.local_grad_vector(), vConsGravity);
			 }
		 }
		 */
		template <typename TElem, typename TFVGeom>
		inline
		void ConsistentGravitySCV( MathVector<dim>* vConsGravitySCV, const TFVGeom& geo, const MathVector<dim> vCornerCoords[], size_t nip, const number* DensitySCV)
		{
			static const int refDim = TElem::dim;
			static const size_t numSh = reference_element_traits<TElem>::reference_element_type::numCorners;
			StdLinConsistentGravityX<refDim> RhoG;
			MathVector<refDim> vConsGravity[numSh];
			MathVector<dim> Gravity; VecSet(Gravity,0.0); Gravity[dim-1]=-9.81;
			RhoG.template prepare<dim>(vConsGravity, numSh, vCornerCoords, DensitySCV, Gravity);
			
			for(size_t ip = 0; ip < nip; ++ip)
			{
				//     get current SCV
				const typename TFVGeom::SCV& scv = geo.scv(ip);
				
				RhoG.template compute<dim>(vConsGravitySCV[ip],scv.local_ip(), scv.JTInv(),scv.local_grad_vector(), vConsGravity);
			}

		}

		template <typename TElem, typename TFVGeom>
		inline
		void ConsistentGravitySCVF( MathVector<dim>* vConsGravitySCVF, const TFVGeom& geo, const MathVector<dim> vCornerCoords[], size_t nip, const number* DensitySCV)
		{
			static const int refDim = TElem::dim;
			static const size_t numSh = reference_element_traits<TElem>::reference_element_type::numCorners;
			StdLinConsistentGravityX<refDim> RhoG;
			MathVector<refDim> vConsGravity[numSh];
			MathVector<dim> Gravity; VecSet(Gravity,0.0); Gravity[dim-1]=-9.81;
			RhoG.template prepare<dim>(vConsGravity, numSh, vCornerCoords, DensitySCV, Gravity);
			
			for(size_t ip = 0; ip < nip; ++ip)
			{
				//     get current SCVF
				const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
				
				RhoG.template compute<dim>(vConsGravitySCVF[ip],scvf.local_ip(), scvf.JTInv(),scvf.local_grad_vector(), vConsGravity);
			}
			

		}

		
		void integration_points( bool &b_scvf, bool &b_scv, GridObject* elem, const MathVector<dim> vCornerCoords[], const MathVector<dim> vGlobIP[], const size_t nip)
		{
			
			if(nip == 0) return;
			b_scvf=true;
			b_scv=true;
			
			DimFV1Geometry<dim> geo;
			geo.update(elem, vCornerCoords, NULL);
			
			const size_t numSCVFip = geo.num_scvf_ips();
			const size_t numSCVip = geo.num_scv_ips();
			
			
			if(geo.num_scvf_ips()!=nip)
			{
				b_scvf=false;
			}
			if(geo.num_scv_ips()!=nip || geo.num_sh() != nip)
			{
				b_scv=false;
			}
			
			for(size_t ip = 0; ip < nip; ++ip)
			{
				if (ip<numSCVFip)
				{
					b_scvf =(b_scvf&&(vGlobIP[ip] == geo.scvf_global_ips()[ip]));
				}
				if (ip<numSCVip)
				{
					b_scv  =(b_scv&&(vGlobIP[ip]  == geo.scv_global_ips()[ip]));
				}
			}

			/*if( !b_scv && !b_scvf )
			{
				printf("Other ips ----------------------------------------------\n");
				for(size_t ip = 0; ip < nip; ++ip)
				{
					printf("x1 vGlobIP  %f   vCornerCoords   %f      geo    %f\n" ,vGlobIP[ip][0], vCornerCoords[ip][0], geo.scvf_global_ips()[ip][0]   );
					printf("y1 vGlobIP  %f   vCornerCoords   %f      geo    %f\n" ,vGlobIP[ip][1], vCornerCoords[ip][1], geo.scvf_global_ips()[ip][1]   );
				   
				}
			}*/

		}
		bool scvf_interface(const number uL, const number uR)
		{
			
			int phaseL=(uL>m_interface_value)? 1:-1;
			int phaseR=(uR>m_interface_value)? 1:-1;
		
			
			int phase = phaseL * phaseR;
			//UG_LOG("SCVF uL uR = "<< uL<<"  "<<uR<<"  Phase = "<<phase<<"\n");
			bool interface = (phase>0)? false:true;
			return interface;
		}
		void cut_element(bool &cut_elem,bool& boolInside, LocalVector* u, const size_t _C_)
		{
			
			bool cut = false;
			(*u).access_all();

			const size_t numSH=(*u).num_all_dof(_C_);
			size_t inside=0;
			size_t outside=0;
			number c;
			for(size_t sh = 0; sh < numSH; ++sh)
			{
				c =(*u)(_C_,sh);
				if (c>m_interface_value)
					inside += 1;
				else
					outside +=1;
			}
			
			if (inside==numSH || outside == numSH)
			{
				cut = false;
			}
			else
			{
				cut = true;
			}
			boolInside = (inside==numSH)? true:false;
			cut_elem=cut;
		}
	
		void cut_element(number &value ,bool &cut_elem, bool &phase2, LocalVector* u, const number interface)
		{
			
			bool cut = false;
			size_t _C_=(*u).num_all_fct() -1 ;
			(*u).access_all();

			const size_t numSH=(*u).num_all_dof(_C_);
			size_t inside=0;
			size_t outside=0;
			number c;
			value = 1000;
			for(size_t sh = 0; sh < numSH; ++sh)
			{
				//c = fmin(1.0, fmax((*u)(_C_,sh),0));
				c = (*u)(_C_,sh);
				value = fmin(value, c);
				if (c>interface)
					inside += 1;
				else
					outside +=1;
			}
			
			if (inside==numSH || outside == numSH)
			{
				cut = false;
				phase2 = (inside==numSH)? true : false;
			}
			else
			{
				cut = true;
				phase2 = false;
			}
			
			cut_elem=cut;
		}
		/*void cut_elementSCV( bool &cut_scv , number* N_value, const number interface, size_t nip, LocalVector* u)
		{
			
			bool cut = false;

			size_t inside=0;
			number c;
			number N;
						
			size_t _C_=(*u).num_all_fct() -1 ;
			(*u).access_all();
			
			for(size_t ip = 0; ip < nip; ++ip)
			{
				c = fmin(1.0, fmax((*u)(_C_,ip),0));
				
				if (c>=interface)
				{
					inside += 1;
					
					 N = 1.0;
				}
				else
				{
					N = -1.0;
				}
				if (N_value != NULL)
					N_value[ip]=N;
				
			}
			
			if (inside != nip && inside != 0)
				cut = true;

			cut_scv=cut;
		}*/
		
		number SmoothMax( const number a, const number b, const number SmoothEps)
		{
			const number diff = a-b;
			return (a + b + sqrt(diff*diff + SmoothEps))/2.0;
		}


		void MatAddTraspose( MathMatrix<dim,dim>& mOut, const MathMatrix<dim,dim> m1)
		{
			for(size_t i = 0; i < dim; ++i)
				for(size_t j = 0; j < dim; ++j)
				{
					mOut[i][j] = m1[i][j] + m1[j][i];
				}
		}
		number MatMultiplyElment( const MathMatrix<dim,dim> m1, const MathMatrix<dim,dim> m2)
		{
			number Sum=0.0;
			for(size_t i = 0; i < dim; ++i)
				for(size_t j = 0; j < dim; ++j)
				{
					Sum += m1[i][j] * m2[i][j];
				}
			return Sum;
		}

		number MatDotTraspose( MathMatrix<dim,dim> m1, const MathVector<dim> m2)
		{
			number s = 0.0;
			for(size_t i = 0; i < dim; ++i)
				for(size_t j = 0; j < dim; ++j)
					s += m2[i] * m1[i][j] * m2[j];
				
			
			return s;
		}

		MathVector<2> MaxMin( number* N_value, const number* density, const size_t nip)
		{
			MathVector<2> DensityJump;
			DensityJump[1] = -100000;
			DensityJump[0] = 100000000;
			for(size_t ip = 0; ip < nip; ++ip)
			{
				if (N_value[ip]>0)  DensityJump[1] = fmax( DensityJump[1], density[ip] );
				else DensityJump[0] = fmin( DensityJump[0], density[ip] );
			}
			return DensityJump;
		}
		MathVector<2> MaxMin( number* N_value, std::vector<number> density, const size_t nip)
		{
			MathVector<2> DensityJump;
			DensityJump[1] = -100000;
			DensityJump[0] = 100000000;
			for(size_t ip = 0; ip < nip; ++ip)
			{
				if (N_value[ip]>0)  DensityJump[1] = fmax( DensityJump[1], density[ip] );
				else DensityJump[0] = fmin( DensityJump[0], density[ip] );
			}
			return DensityJump;
		}
		MathVector<2> MaxMin( number* N_value,const number* viscosity, const number* density, const size_t nip)
		{
			MathVector<2> ViscJump;
			ViscJump[1] = -100000;
			ViscJump[0] = 100000000;
			for(size_t ip = 0; ip < nip; ++ip)
			{
				if (N_value[ip]>0)  ViscJump[1] = fmax( ViscJump[1], density[ip] * viscosity[ip]);
				else ViscJump[0] = fmin( ViscJump[0], density[ip] * viscosity[ip]);
			}
			return ViscJump;
		}
		MathVector<2> MaxMin( number* N_value, std::vector<number> viscosity, std::vector<number> density, const size_t nip)
		{
			MathVector<2> ViscJump;
			ViscJump[1] = -100000;
			ViscJump[0] = 100000000;
			for(size_t ip = 0; ip < nip; ++ip)
			{
				if (N_value[ip]>0)  ViscJump[1] = fmax( ViscJump[1], density[ip] * viscosity[ip]);
				else ViscJump[0] = fmin( ViscJump[0], density[ip] * viscosity[ip]);
			}
			return ViscJump;
		}
		void interface_position( MathVector<dim>* x_interface, const MathVector<dim> vCornerCoords[], const number* N_value, const number rho_l, const number rho_g, const number interface ,DimFV1Geometry<dim> geo,const size_t nip, LocalVector*  u)
		{
						
			size_t _C_=(*u).num_all_fct() -1 ;
			(*u).access_all();
			
			MathVector<dim> x, DX;
			size_t NumSCVF = geo.num_scvf();
			number n[nip];
			number theta_to, theta_from, c_to, c_from, DC;
			number rho;
			for(size_t ip = 0; ip < nip; ++ip)
			{
				n[ip] = 0.0;
				VecSet(x_interface[ip],0.0);
			}
			for(size_t ip = 0; ip < NumSCVF; ++ip)
			{
				const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
				size_t from=scvf.from();
				size_t to=scvf.to();
				
				if (N_value[from]*N_value[to] < 0.0)
				{

					//c_from = fmin(1.0, fmax((*u)(_C_,from),0));
					//c_to = fmin(1.0, fmax((*u)(_C_,to),0));
					
					c_from = (*u)(_C_,from);
					c_to = (*u)(_C_,to);

					
					DC=c_to-c_from;
					VecSubtract(DX,vCornerCoords[to],vCornerCoords[from]);
					
					theta_to=  (c_to   - interface)/DC;
					theta_from=(c_from - interface)/DC;
					
					VecScaleAppend(x_interface[to], theta_to,   DX);
					VecScaleAppend(x_interface[from],   theta_from, DX);
					
					n[from] += 1.0;
					n[to] += 1.0;
				}
			}
			for(size_t ip = 0; ip < nip; ++ip)
			{

				rho = (N_value[ip] > 0.0)? rho_g:rho_l;
				VecScale(x_interface[ip], x_interface[ip], (rho_l-rho_g) / (n[ip] * rho) );
				
			}
			
		}
		number interface_density(number VolFrac_surf, LocalVector* u, GridObject* elem, const MathVector<dim> vCornerCoords[], const size_t ip, const number rho_s, const number rho_a, const bool m_scvf, const bool m_scv)
		{
			
			size_t _C_=(*u).num_all_fct() -1 ;
			size_t numSH=(*u).num_all_dof(_C_);
			if (m_scvf)
			{
				number c1, c2, rho1, rho2;
				
				DimFV1Geometry<dim> geo;
				geo.update(elem, vCornerCoords, NULL);
				
				
				const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
				
				//c1 = fmin(1.0, fmax((*u)(_C_,scvf.from()),0));
				//c2 = fmin(1.0, fmax((*u)(_C_,scvf.to()),0));
				
				c1 = (*u)(_C_,scvf.from());
				c2 = (*u)(_C_,scvf.to());
				
				rho1=(rho_s-rho_a)*c1+rho_a;
				rho2=(rho_s-rho_a)*c2+rho_a;
				
				return fmin(rho1,rho2);
			}
			else
			{
				number rho=1e10, c;
				number rho_aux;
				for(size_t sh = 0; sh < numSH; ++sh)
				{
					//c = fmin(1.0, fmax((*u)(_C_,sh),0));
					c = (*u)(_C_,sh);
					rho_aux=(rho_s-rho_a)*c+rho_a;
					rho = fmin (rho_aux,rho);
					
				}
				return rho;
				
				
			}


		}

		template <typename TElem, typename TFVGeom>
		inline
		void PropertiesJump(LocalVector* u, const size_t _C_, const TFVGeom& geo, int* JumpShape, const DataImport<number, dim>& DensitySCV, const DataImport<number, dim>& KinViscSCV, DataImport<MathVector<dim>, dim>& SourceSCV,size_t numSh, const bool interface, bool* Phase2, number& mu_l, number& mu_g, number& rho_l, number& rho_g, MathVector<dim>& Source_l, MathVector<dim>& Source_g, bool* GFM)
		{
			
			/*UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");
			
			interface = cut_interface(JumpShape, numSh);*/
			
			
			
			if(!interface)
			{
				return;
			}
			
			Jump_Shape( u, _C_, numSh, JumpShape);
			/*if (true)
			{
				
				for(size_t sh = 0; sh < geo.num_scvf(); ++sh)
					printf("Vol[%zu] = %f\n", sh, VolFraction[sh]);
				
				printf("Pos[0] = %f     %f\n", geo.scv_global_ips()[0][0],geo.scv_global_ips()[0][1]);
				printf("Pos[1] = %f     %f\n", geo.scv_global_ips()[1][0],geo.scv_global_ips()[1][1]);
				printf("Pos[2] = %f     %f\n", geo.scv_global_ips()[2][0],geo.scv_global_ips()[2][1]);
				

			}*/
			
			RhoMuSource(mu_l,  mu_g,  rho_l, rho_g, Source_l, Source_g, DensitySCV, KinViscSCV, SourceSCV, JumpShape, numSh);

			
			/*number Cval=0.0;
			number vol = 0.0;
			//     loop Sub Control Volumes (SCV)
			for(size_t ip = 0; ip < geo.num_scv(); ++ip)
			{
				//     get current SCV
				const typename TFVGeom::SCV& scv = geo.scv(ip);
				
				//     get associated node
				const int sh = scv.node_id();
				
				Cval += VolFraction[sh]*scv.volume();
				vol += scv.volume();

			}
			Cval *= 1.0/vol;*/
			this->template phase<TFVGeom>(geo,JumpShape, Phase2, GFM);
			
		}
	
		void Jump_Shape( LocalVector* u, const size_t _C_, const size_t numSH,  int* JumpShape)
		{
			
			for(size_t sh = 0; sh < numSH; ++sh)
			{
				number c =(*u)(_C_,sh);
				if (c>m_interface_value)
					JumpShape[sh]= 1;
				else
					JumpShape[sh]= -1;
			}
			
		}
	
		bool cut_interface(const DataImport<number, dim>& JumpShape, const size_t numSh)
		{
			
			number max=0;
			number min=0;
			
			bool interface=false;
			for(size_t sh = 0; sh < numSh; ++sh)
			{
				max=fmax(JumpShape[sh],max);
				min=fmin(JumpShape[sh],min);
				
			}
			if (-0.5<(max+min) && (max+min)<0.5)
			{
				interface=true;
				
			}
			return interface;


		}

		template <typename TElem, typename TFVGeom>
		inline
		void phase(const TFVGeom& geo,int* JumpShape, bool* Phase2, bool* GFM)
		{
			
			for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
			{
				//     get current SCV
				const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
				GFM[ip]=true;
				if(JumpShape[scvf.to()]*JumpShape[scvf.from()]<0.0)
				{
					Phase2[ip]=true;
				}
				else
				{
					
					
					if(JumpShape[scvf.to()]>0.0)
					{
						Phase2[ip]=true;
					}
					else
					{
						Phase2[ip]=false;
					}
				}
				
			}


		}
		template <typename TElem, typename TFVGeom>
		inline
		bool cut_interface_concentration( const TFVGeom& geo, const LocalVector& u, const size_t _C_, const size_t numSh, number interface_value, size_t& num_inside, size_t& num_outside, number& value_in, number& value_out, number* JumpShape)
		{
			bool interface = false;

			size_t inside=0;
			size_t outside=0;
			value_in = 0.0;
			value_out = 0.0;
			
			for(size_t sh = 0; sh < numSh; ++sh)
			{

				if (u(_C_, sh)>interface_value)
				{
					inside += 1;
					JumpShape[sh] = 1;
					value_in += u(_C_, sh);
				}
				else
				{
					outside +=1;
					JumpShape[sh] = -1;
					value_out += u(_C_, sh);
				}
			}
			
			if (inside==numSh || outside == numSh)
			{
				interface = false;
				value_in = (value_in + value_out) / (inside+outside);
				value_out = value_in;
				inside = inside+outside;
				outside = inside;
			}
			else
			{
				interface = true;
				value_in = value_in / inside;
				value_out = value_out / outside;
			}
			
			num_inside = inside;
			num_outside = outside;
			
			
			return interface;


		}

		void RhoMuSource( number& mu_l, number& mu_g, number& rho_l, number& rho_g, MathVector<dim>& vSource_l, MathVector<dim>& vSource_g, const DataImport<number, dim>& DensitySCV, const DataImport<number, dim>& KinViscSCV, const DataImport<MathVector<dim>, dim>& SourceSCV, int* JumpShape, const size_t numSh)
		{

			bool boolSource = (SourceSCV.data_given()) ? true : false;
			
			number mu_2 = 0.0, mu_1 = 0.0;
			number rho_2 = 0.0, rho_1 = 0.0;
			

			MathVector<dim> Source_1;
			MathVector<dim> Source_2;
			
			VecSet(Source_1,0.0);
			VecSet(Source_2,0.0);
			
			
			int Count_1 = 0;
			int Count_2 = 0;
			
			for(size_t sh = 0; sh < numSh; ++sh)
			{
				if (JumpShape[sh]>0)
				{
					mu_2  +=  1.0/(DensitySCV[sh] * KinViscSCV[sh]);
					rho_2 +=  (DensitySCV[sh]) ;

					if(boolSource) {VecScaleAppend(Source_2 , 1.0  ,SourceSCV[sh] );}
					
					Count_2 +=1;
					
				}
				else
				{
					mu_1  +=  1.0/(DensitySCV[sh] * KinViscSCV[sh]);
					rho_1 +=  DensitySCV[sh] ;
					if(boolSource) {VecScaleAppend(Source_1   , 1.0  ,SourceSCV[sh] );}
					Count_1 +=1;
					
				}
				
			}
			mu_1 =  Count_1 / mu_1  ;
			mu_2 =  Count_2 / mu_2  ;
			rho_1 = rho_1 / Count_1;
			rho_2 = rho_2 / Count_2;
			
			if(boolSource)
			{
				VecScale(Source_1,Source_1,1.0/Count_1);
				VecScale(Source_2,Source_2,1.0/Count_2);
			}
			
			if ((mu_2 < mu_1)||(mu_2<=0.0)||(mu_1<=0.0))
			{
				printf("Mu2 = %f    Mu1 = %f\n",mu_2,mu_1);
				printf("Rho2 = %f    Rho1 = %f\n",rho_2,rho_1);
				
				for(size_t sh = 0; sh < numSh; ++sh)
					printf("mu[%zu] = %f\n", sh, DensitySCV[sh] * KinViscSCV[sh]);
				for(size_t sh = 0; sh < numSh; ++sh)
					printf("rho[%zu] = %f\n", sh, DensitySCV[sh] );
				for(size_t sh = 0; sh < numSh; ++sh)
					printf("JumpShape[%zu] = %f\n", sh, JumpShape[sh] );
				
				
				
				UG_THROW("Viscosity in phase 1 is lower that phase 2");
			}
			if ((rho_2 < rho_1)||(rho_2<=0.0)||(rho_1<=0.0))
				UG_THROW("Density in phase 1 is lower that phase 2");
			
			
			
			
			mu_l = mu_2;
			mu_g = mu_1;
			rho_l = rho_2;
			rho_g = rho_1;

			vSource_l = Source_2;
			vSource_g = Source_1;

		}
		void InterMu( number& mu_l, number& mu_g, const number Mu[], int JumpShape[], const size_t numSh)
		{

	
			
			number mu_2 = 0.0, mu_1 = 0.0;
			
			
			int Count_1 = 0;
			int Count_2 = 0;
			
			for(size_t sh = 0; sh < numSh; ++sh)
			{
				if (JumpShape[sh]>0)
				{
					mu_2  +=  1.0/Mu[sh];
					Count_2 +=1;
					
				}
				else
				{
					mu_1  +=  1.0/Mu[sh];
					Count_1 +=1;
					
				}
				
			}
			mu_1 =  Count_1 / mu_1  ;
			mu_2 =  Count_2 / mu_2  ;

			
			if ((mu_2 < mu_1)||(mu_2<=0.0)||(mu_1<=0.0))
			{
				printf("Mu2 = %f    Mu1 = %f\n",mu_2,mu_1);
				
				for(size_t sh = 0; sh < numSh; ++sh)
					printf("mu[%zu] = %f\n", sh, Mu[sh]);
				for(size_t sh = 0; sh < numSh; ++sh)
					printf("JumpShape[%zu] = %f\n", sh, JumpShape[sh] );
				
				
				
				UG_THROW("Viscosity in phase 1 is lower that phase 2");
			}
			
			
			
			
			mu_l = mu_2;
			mu_g = mu_1;

		}

		template <typename TElem, typename TFVGeom>
		inline
		void InterfaceSCVFShapes( number** SCVFinterShape, const TFVGeom& geo, const LocalVector& u, const DataImport<number, dim>& JumpShape, size_t numSh, const size_t _C_ )
		{
			/////////////////////////////////////////////////////////////////////////////
			// Calculation X_interface
			/////////////////////////////////////////////////////////////////////////////
			
			/*number** SCVFinterShape = new number*[numSCVF];
			for(size_t count=0; count<numSCVF; count++)
				SCVFinterShape[count] = new number[numSh];*/
			
			UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");
			
			size_t numSCVF = geo.num_scvf();
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			static const int refDim = ref_elem_type::dim;
			
			MathVector<refDim> vLocIP_inter;
			VecSet(vLocIP_inter,0.0);
			
			MathVector<refDim> vLocIP_SCVF[numSCVF];
			
			number theta_to, c_to, c_from, DC;
			

			number count_interface=0;

			for(size_t ip = 0; ip < numSCVF; ++ip)
			{
				const typename FV1Geometry<TElem, dim>::SCVF scvf = geo.scvf(ip);

				const size_t from=scvf.from();
				const size_t to=scvf.to();
				
				VecSet(vLocIP_SCVF[ip],0.0);
				
				if (JumpShape[from]*JumpShape[to] < 0.0)
				{
					c_from = u(_C_,from);
					c_to = u(_C_,to);

					DC=c_to-c_from;
					
					theta_to=  (c_to   - m_interface_value)/DC;
					
					count_interface += 1.0;
					VecScaleAppend(vLocIP_SCVF[ip], 1.0  ,geo.scv_local_ips()[to],-1.0 * theta_to, geo.scv_local_ips()[to],theta_to,geo.scv_local_ips()[from]);
					VecScaleAppend(vLocIP_inter   , 1.0  ,vLocIP_SCVF[ip]);
					
				}
				else
				{
					VecScaleAppend(vLocIP_SCVF[ip], 1.0  ,geo.scv_local_ips()[to],-1.0 * 0.5, geo.scv_local_ips()[to],0.5,geo.scv_local_ips()[from]);
				}

			}
			

			VecScale(vLocIP_inter,vLocIP_inter,1.0/count_interface);
			
			LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();


		//    Reference Mapping

			ReferenceMapping<ref_elem_type, dim> mapping(geo.scv_global_ips());
			
			
			

			
			for(size_t ip = 0; ip < numSCVF; ++ip)
			{
				
				VecScale(vLocIP_SCVF[ip],vLocIP_SCVF[ip],0.5);
				VecScaleAppend(vLocIP_SCVF[ip], 0.5  ,vLocIP_inter);
							   
				rTrialSpace.shapes(SCVFinterShape[ip], vLocIP_SCVF[ip]);
				
				/*const typename FV1Geometry<TElem, dim>::SCVF scvf = geo.scvf(ip);
				const size_t from=scvf.from();
				const size_t to=scvf.to();
				printf("vGlobIP_SCVF[from][0] =  %f \n", geo.scv_global_ips()[from][0]);
				printf("vGlobIP_SCVF[from][1] =  %f \n", geo.scv_global_ips()[from][1]);
				printf("vGlobIP_SCVF[to][0] =  %f \n",geo.scv_global_ips()[to][0]);
				printf("vGlobIP_SCVF[to][1] =  %f \n",geo.scv_global_ips()[to][1]);
				printf(" Shapes  ------------------------------------");
				for(size_t sh = 0; sh < numSh; ++sh)
					printf("SCVF shape [%zu][%zu] = %f\n",ip,sh,interShape[ip][sh]);*/

				
			}
		}

		template <typename TElem, typename TFVGeom>
		inline
		void InterfaceShape(number* InterfaceShape, const TFVGeom& geo, const DataImport<number, dim>& VolFraction, const DataImport<number, dim>& JumpShape, size_t numSh)
		{
			/////////////////////////////////////////////////////////////////////////////
			// Calculation X_interface
			/////////////////////////////////////////////////////////////////////////////
			
			/*number** SCVFinterShape = new number*[numSCVF];
			for(size_t count=0; count<numSCVF; count++)
				SCVFinterShape[count] = new number[numSh];*/
			
			UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");
			
			size_t numSCVF = geo.num_scvf();
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			static const int refDim = ref_elem_type::dim;
			
			MathVector<refDim> vLocIP_inter;
			VecSet(vLocIP_inter,0.0);

			
			number theta_to, c_to, c_from, DC;
			

			number count_interface=0;

			for(size_t ip = 0; ip < numSCVF; ++ip)
			{
				const typename FV1Geometry<TElem, dim>::SCVF scvf = geo.scvf(ip);

				const size_t from=scvf.from();
				const size_t to=scvf.to();
				
				
				if (JumpShape[from]*JumpShape[to] < 0.0)
				{
					c_from = VolFraction[from];
					c_to = VolFraction[to];

					DC=c_to-c_from;
					
					theta_to=  (c_to   - m_interface_value)/DC;
					
					count_interface += 1.0;
					VecScaleAppend(vLocIP_inter, 1.0  ,geo.scv_local_ips()[to],-1.0 * theta_to, geo.scv_local_ips()[to],theta_to,geo.scv_local_ips()[from]);
					
				}

			}
			

			VecScale(vLocIP_inter,vLocIP_inter,1.0/count_interface);
			
			LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();


		//    Reference Mapping

			ReferenceMapping<ref_elem_type, dim> mapping(geo.scv_global_ips());
			
			rTrialSpace.shapes(InterfaceShape, vLocIP_inter);
			

		}
	public:
		void set_particle_density(float R) {
			rho_s = R;
		}
		void set_air_density(float R) {
			rho_a = R;
		}
		void set_fluid_Visc(float R) {
			mu_a = R;
		}
		void set_particle_kinVisc(float R) {
			nu_s = R;
		}
		void set_particle_diameter(float R) {
			dp = R;
		}
		void set_alpha_max(float R) {
			alpha_max = R;
		}
		void set_alpha_min(float R) {
			alpha_min = R;
		}
		void set_packing_factor(float R) {
			m_packing_factor = R;
		}
		void set_FR(float R) {
			Fr = R;
		}
		void set_B_phi(float R) {
			B_phi = R;
		}
		void set_deltaGamma(float R) {
			deltaGamma = R;
		}
		void set_limit(float R) {
			m_limit = R;
		}
		void set_time_step_factor(float R) {
			m_dt = R;
		}
		void set_interface_volume_fraction(float R) {
			m_interface_value = R;
		}
		void set_FrictionMu_1(float R) {
			FricMu_1 = R;
		}
		void set_FrictionMu_2(float R) {
			FricMu_2 = R;
		}
		void set_I_0(float R) {
			I_0 = R;
		}
		void set_deltaI(float R) {
			deltaI = R;
		}
		void set_deltaPs(float R) {
			deltaPs = R;
		}
		void set_reference_pressure(float R) {
			m_P0 = R;
		}
		void set_gravity(float R) {
			m_gravitation = R;
		}
	
		void set_relative_vel_error(float R) {
			m_RelVelError = R;
		}
	
		void set_bool_particle_pressure_force(bool R) {
			m_bParticleGradientForce = R;
		}
		void set_bool_consistent_gravity(bool R) {
			m_bConsistentGravity = R;
		}
		void set_drag_model(size_t R) {
			drag_model = R;
		}
		void set_bool_initialized(bool R) {
			m_init = R;
			rho_max = (rho_s-rho_a) * m_packing_factor + rho_a;
			
		}

		number Density_s(){ return rho_s;}
		number Density_a(){ return rho_a;}
		number Density_max(){ return rho_max;}
		number Viscosity_a(){ return mu_a;}
		number KinViscosity_s(){ return nu_s;}
		number Alpha_max(){ return alpha_max;}
		number packing_factor(){ return m_packing_factor;}
		number DT(){ return m_dt;}
		number ReferencePressure(){ return m_P0;}
		number gravity(){ return m_gravitation;}
		number diameter(){ return dp;}
		number interface_value(){ return m_interface_value;}
		number FrictionMu_1(){ return FricMu_1;}
		number FrictionMu_2(){ return FricMu_2;}
		number param_I_0(){ return I_0;}
		number param_deltaI(){ return deltaI;}
		number param_deltaPs(){ return deltaPs;}
		number param_deltaGamma(){ return deltaGamma;}
		number Epsilon(){return epsilon;}
		bool ParticleGradientForce(){ return m_bParticleGradientForce;}
		bool boolConsistentGravity(){ return m_bConsistentGravity;}
		
		size_t DragModel(){ return drag_model;}
		bool valid(){ return m_init;}


	protected:
		float m_P0, rho_max, rho_s, rho_a, mu_a, nu_s, dp, Fr, B_phi, alpha_max, alpha_min, deltaGamma, FricMu_1, FricMu_2, I_0, deltaI, deltaPs, m_limit, m_dt, m_interface_value, m_packing_factor, m_gravitation, m_RelVelError, epsilon, m_delta;
		size_t drag_model;
		bool m_bParticleGradientForce, m_bConsistentGravity, m_init;
        
};


template <int TDim, int TWorldDim = TDim>
class LineWriter
{


	public:
	///	dimension of reference element
		static const int dim = TDim;

	///	dimension of world
		static const int worldDim = TWorldDim;

	public:
	/// construct object and initialize local values and sizes
		LineWriter():
			num_columns(4),
			id_width(5),
			double_width(20),
			double_precision(10),
			separator(","),
			newline("\n")
		{
			LINE_SIZE =  id_width + num_columns * double_width + (num_columns ) * separator.size()  + newline.size();
		}
		
	

	
		
	public:


		std::string format_line(int rank, const std::vector<double>& values)
		{
			
			UG_ASSERT(values.size() == (size_t)num_columns, "PropertiesInterface Writer: Values number mismatch colums number!");
			
			std::ostringstream oss;
			oss << std::setw(id_width) << std::setfill('0') << (rank+1);
			
			
			// Append each double column
			for (int i = 0; i < num_columns; ++i) {
				oss << separator;
				
				std::ostringstream tmp;
				tmp << std::showpos << std::fixed << std::setprecision(double_precision) << values[i];
				std::string num_str = tmp.str();
				char sign = num_str[0];
				num_str = num_str.substr(1);
				
				// Pad left with zeros if too short
				if ((int)num_str.size() < double_width-1)
				{
					num_str.insert(0, double_width -1 - num_str.size(), '0');
				}
				// Truncate if too long
				else if ((int)num_str.size() > double_width-1)
				{
					num_str = num_str.substr(0, double_width-1);
				}
				oss <<	sign	<<num_str;

			}
			
			// Append newline
			oss << newline;

			std::string line = oss.str();
	
		   // Safety check: line length
			UG_ASSERT((int)(line.size()) == LINE_SIZE, "PropertiesInterface Writer: Line length mismatch!");

			return line;
			
		}
 
		void write_line(const std::string filename, int rank, const std::string Headers, double value1, double value2, double value3, int boolSolution)
		{
			int local_rank = 0;
			#ifdef UG_PARALLEL
				local_rank = pcl::ProcRank();
			#endif
			if(local_rank != 0) return;

			std::vector<double> values(num_columns);

			values[0] = value1;
			values[1] = value2;
			values[2] = value3;
			values[3] = boolSolution;
		
			std::string line = format_line(rank, values);
			//std::string empty_line(LINE_SIZE, ' '); // or '\0'
			// Compute offset using plain int
			int offset =(rank == 0)? 0 : Headers.size() + rank * LINE_SIZE;
			
			
			// Open file in binary read/write mode
			std::fstream file(filename, std::ios::in | std::ios::out | std::ios::binary);

			// If file does not exist, create it
			if (!file.is_open()) {
				std::ofstream create_file(filename, std::ios::out | std::ios::binary);
				create_file.close();
				file.open(filename, std::ios::in | std::ios::out | std::ios::binary);
			}

			// Seek to the offset
			file.seekp(offset); // just plain int, no pos_type needed
			if (!file) {
				std::cerr << "Error seeking to offset " << offset << std::endl;
				return;
			}

			// Write the line
			if(rank == 0) file.write(Headers.c_str(), Headers.size());
			file.write(line.c_str(), LINE_SIZE);
			if (!file) {
				std::cerr << "Error writing line at offset " << offset << std::endl;
			}

			file.close();
		}

	

	protected:
		int num_columns, id_width, double_width, double_precision, LINE_SIZE;
		std::string separator,newline;


		


};




}

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC_NV_CUT_ELEMENT_ */


