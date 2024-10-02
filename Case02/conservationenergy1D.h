#ifndef CONSERVATION_ENERGY_1D_H
#define CONSERVATION_ENERGY_1D_H

#include "sparselizard.h"

using namespace sl;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ELECTRON ENERGY CONSERVATION: ∇[(5/2)Tnv + (-K.∇T)] = Q_R + Q_E  [J/(m³s)]

First two LHS terms: Qconv = ∇[(0.5*m_i*v² + 5T)nv]
Third LHS term     : Qcond = ∇.(-K.∇T)

heat conduction coefficient: K = (k0e+koi) * T^(5/2)
k0e ≃ 2000
k0i ≃ 60

Q_R = R(ve - vi) = 0; Due to quasineutrality and ambipolarity
Q_E = 0, Assuming cold neutrals
-----------------------------------------------------------------*/ 

// double a=1.0, b=0.03695; // for qconv only
double a=1.0, b=1.0; // qconv + qcond
parameter alpha, beta;
alpha.setvalue(selectall(), a);
 beta.setvalue(selectall(), b);

// Heat convection:
plasma_flow += integral(SOL, alpha * e*(n_o*c_st) * (5.0) * (
      dof(T)*N*grad(M) + dof(T)*grad(N)*M + grad(dof(T))*N*M
    + T*dof(N)*grad(M) + T*grad(dof(N))*M + grad(T)*dof(N)*M
    + T*N*grad(dof(M)) + T*grad(N)*dof(M) + grad(T)*N*dof(M)
    - 2*(T*N*grad(M) + T*grad(N)*M + grad(T)*N*M)
) * tf(T)); // Linearized

// Vdp term from the pressure equation
plasma_flow += integral(SOL, -alpha*e*2*n_o*c_st*M*(N*grad(dof(T))+dof(T)*grad(N)) * tf(T));

// Heat conduction: block number = 2
expression k0e = 2000, k0i = 60;
expression K = k0e*pow(T, 5.0/2.0);  // Thermal diffusivity
plasma_flow += integral(SOL, beta * (k0e)*(         pow(T, 5.0/2.0)        *grad(dof(T)) 
                                            + 2.5 * pow(T, 3.0/2.0)*dof(T) *grad(T) 
                                            - 2.5 * pow(T, 5.0/2.0)        *grad(T) 
) * grad(tf(T)), 0, 2); // Newton linearization

// Heat source
// plasma_flow += integral(upstream, -q_u * tf(T));    // NBC AT UPSTREAM (or)
plasma_flow += integral(SOL, -Spn * tf(T)); // Volumetric net source

plasma_flow += integral(target, on(SOL, +Spn) * tf(T)); // Make Source zero at target

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
STABILIZED FORMULATION:
----------------------*/

parameter epsT;
epsT.setvalue(selectall(), 1e-20);
expression diffusivityT = K;

//-----------------------------------------------------------------------------------------
// A) Stabilization parameter (tau)
parameter tauT;

//  A.1) by SANJAY
const int problemdimensionT = mymesh.getdimension();
expression hT = pow(meshsize(2), 1.0 / problemdimensionT);    // characteristic element length
if (SUPG_TYPE=="SANJAY"){
    double mk = 1.0/3.0;
    tauT.setvalue(SOL, 1.0 / sqrt(pow(2*norm(M)/hT, 2) + pow(8*diffusivityT/(mk*hT*hT), 2))); 
}
// A.2) by CODINA
expression invtaujacaT = norm(invjac() * M);
expression invtaujacdT = doubledotproduct(invjac(), invjac());
expression denomT  = sqrt((invtaujacaT) * (invtaujacaT) + (invtaujacdT * diffusivityT) * (invtaujacdT * diffusivityT));
if (SUPG_TYPE=="CODINA"){
    tauT.setvalue(SOL, ifpositive(denomT-epsT, 1.0 / denomT, 0.0));
}

//-----------------------------------------------------------------------------------------
// B) RESIDUALS (R)
expression Rt_supg = alpha * e*(n_o*c_st) * (5.0) * (
      dof(T)*N*grad(M) + dof(T)*grad(N)*M + grad(dof(T))*N*M
); // Linearized  // T is dof. Linearized terms whereever necessary

//-----------------------------------------------------------------------------------------
// C) MODIFIED TEST FUNCTIONS (P)
expression Pwt_supg = grad(tf(M))*M;

// if (SUPGT == true){ plasma_flow += integral(SOL, Pwt_supg * tauT * Rt_supg); }


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#endif