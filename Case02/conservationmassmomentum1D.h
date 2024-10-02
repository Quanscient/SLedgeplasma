#ifndef CONSERVATION_MASSMOMENTUM_1D_H
#define CONSERVATION_MASSMOMENTUM_1D_H

#include "sparselizard.h"

using namespace sl;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MASS CONSERVATION    :  ∇.(nv) - Si = 0
Normalization factor :  n_o * c_st
Normalized PDE       :  ∇.(NM) - Sin = 0 -> (N∇M + M∇N) - Sin = 0
--------------------------------------------------------------*/
plasma_flow += integral(SOL, (dof(N)*grad(M) + M*grad(dof(N))
                            +  N*grad(dof(M)) + dof(M)*grad(N)
                            -  N*grad(M) - M*grad(N)) * tf(N) );    // 1st term: Linearization of divergence term
plasma_flow += integral(SOL, -Sin * tf(N));                         // 2nd term: Source term


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MOMENTUM CONSERVATION :  ∇(mnv² + 2nT) = -mv*Rcx -> ∇(mnv² + 2nT) + mv*Rcx = 0
Normalization factor  :  n_o * c_st²
Normalized PDE        :  (2NM*dM/dx + M²dN/dx) + (2/(m_i*c_s²))*(TdN/dx + NdT/dx) + M*Scxn = 0
with: 2/(m_i*c_s²) = 1/Tt
---------------------------------------------------------------------------------------*/
plasma_flow += integral(SOL, 2*(dof(N)*grad(M)*M + N*grad(M)*dof(M) + N*grad(dof(M))*M - 2*N*grad(M)*M) * tf(M));    // Linearization of 1st term
plasma_flow += integral(SOL, (M*M*grad(dof(N)) + 2*M*dof(M)*grad(N) - 2*M*M*grad(N) ) * tf(M) );    // Linearization of 2nd term
// plasma_flow += integral(SOL, (1.0/Tt) * (T*grad(dof(N)) + dof(N)*grad(T)) * tf(M) );    // 3rd term
plasma_flow += integral(SOL, (1.0/Tt) * (
      dof(T)*grad(N) + N*grad(dof(T))  +  T*grad(dof(N)) + dof(N)*grad(T)  -  T*grad(N) - N*grad(T)
) * tf(M) );    // 3rd term Linearized
plasma_flow += integral(SOL, -dofSmn * tf(M));     // 4th term: zero for 2point model

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
STABILIZED FORMULATION:
----------------------*/

parameter eps;
eps.setvalue(selectall(), 1e-20);
expression diffusivity = 0.0;

//-----------------------------------------------------------------------------------------
// A) Stabilization parameter (tau)
parameter tau;

//  A.1) by SANJAY
const int problemdimension = mymesh.getdimension();
expression h = pow(meshsize(2), 1.0 / problemdimension);    // characteristic element length
if (SUPG_TYPE=="SANJAY"){
    double mk = 1.0/3.0;
    tau.setvalue(SOL, 1.0 / sqrt(pow(2*norm(M)/h, 2) + pow(8*diffusivity/(mk*h*h), 2))); 
}
// A.2) by CODINA
expression invtaujaca = norm(invjac() * M);
expression invtaujacd = doubledotproduct(invjac(), invjac());
expression denom  = sqrt((invtaujaca) * (invtaujaca) + (invtaujacd * diffusivity) * (invtaujacd * diffusivity));
if (SUPG_TYPE=="CODINA"){
    tau.setvalue(SOL, ifpositive(denom-eps, 1.0 / denom, 0.0));
}

//-----------------------------------------------------------------------------------------
// B) RESIDUALS (R)
// works better for Rcx=0, but not for Rcx!=0: Use this for 2-point model Validation
expression Rm_supg = ( 2*N*(grad(M)*dof(M) + grad(dof(M))*M - grad(M)*M) 
                    + (2*M*dof(M) - M*M) * grad(N)
                    + (1.0/Tt) * (T*grad(N) + N*grad(T))
                    + (-dofSmn)
);  // M is dof. Linearized terms whereever necessary

expression Rm_pspg = ( 2*dof(N)*grad(M)*M
                    + M*M*grad(dof(N))
                    + (1.0/Tt) * (T*grad(dof(N)) + dof(N)*grad(T)) 
                    + (-Smn)
);  // N is dof. Linearized terms whereever necessary
expression Rc_lsic = (N*div(dof(M)) + dof(M)*grad(N) - Sin);

//-----------------------------------------------------------------------------------------
// C) MODIFIED TEST FUNCTIONS (P)
expression Pw_lsic = div(tf(M));
expression Pw_supg = grad(tf(M))*M;
expression Pq_pspg = grad(tf(N));

if (LSIC == true){ plasma_flow += integral(SOL, Pw_lsic * tau * Rc_lsic); }
if (SUPG == true){ plasma_flow += integral(SOL, Pw_supg * tau * Rm_supg); }
if (PSPG == true){ plasma_flow += integral(SOL, Pq_pspg * tau * Rm_pspg); }


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#endif