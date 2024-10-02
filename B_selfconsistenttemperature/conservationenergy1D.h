#ifndef CONSERVATION_ENERGY_1D_H
#define CONSERVATION_ENERGY_1D_H

#include "sparselizard.h"

using namespace sl;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ENERGY CONSERVATION: ∇[(0.5*m_i*v² + 5T)nv + (-K.∇T)] = Q_R + Q_E

First two LHS terms: Qconv = ∇[(0.5*m_i*v² + 5T)nv]
Third LHS term     : Qcond = ∇.(-K.∇T)

heat conduction coefficient: K = (k0e+koi) * T^(5/2)
k0e ≃ 2000
k0i ≃ 60

Q_R = R(ve - vi) = 0; Due to quasineutrality and ambipolarity
Q_E = 0, Assuming cold neutrals
-----------------------------------------------------------------*/ 

// Heat convection: block number = 1
// plasma_flow += integral(SOL, (n_o*c_st*c_st*c_st) 
//                            * (m_i/2.0) * (3*N*M*M*gradM + gradN*M*M*M) * tfT()), 0, 1);
// plasma_flow += integral(SOL, (n_o*c_st) * 5.0 * (T*N*gradM + T*gradN*M + gradT*N*M) * tfT, 0, 1);

// Heat conduction: block number = 2
expression k0e = 2000, k0i = 60;
expression K = k0e*pow(T, 5.0/2.0);  // Thermal diffusivity
plasma_flow += integral(SOL, (k0e)*(        pow(T, 5.0/2.0)        *grad(dof(T)) 
                                    + 2.5 * pow(T, 3.0/2.0)*dof(T) *grad(T) 
                                    - 2.5 * pow(T, 5.0/2.0)        *grad(T) 
) * grad(tf(T)), 0, 2); // Netwon linearization

// NBC AT UPSTREAM
plasma_flow += integral(upstream, -q_u * tf(T) );

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#endif