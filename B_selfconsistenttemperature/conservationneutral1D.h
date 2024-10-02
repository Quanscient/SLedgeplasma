#ifndef CONSERVATION_NEUTRAL_1D_H
#define CONSERVATION_NEUTRAL_1D_H

#include "sparselizard.h"

using namespace sl;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DIFFUSION OF NEUTRALS:
∇.(-D.∇nn) - Sn = 0
NBC at target ->  -D.∇nn = Γ_recycling = -R*(nv)

Normalized with Sp_o=n_o * c_st:
∇.((-D/c_st).∇Nn) - Snn = 0
NBC at target ->  (-D/c_st).∇Nn = -R*(NM)
----------------------------------------------*/

// DIFFUSION COEFFICIENT, D [m²/s] : As described in UEDGE paper   : https://doi.org/10.1063/1.873488
expression Diz   = T/(m_i*n_o*N*Riz);       // Considering only electron impact on H-atom and ignoring charge exchange.
expression Dizcx = T/(m_i*n_o*N*(Riz+Rcx)); // Considering electron impact on H-atom and charge exchange.
expression D = Dizcx/c_st;  // Effective diffusion coefficient

// WEAK FORMULATION
plasma_flow += integral(SOL,     D*grad(dof(Nn))*grad(tf(Nn)));    // integ. by parts: term 1
plasma_flow += integral(target,  R*N*(-M) * tf(Nn) );            // integ. by parts: term 2 -> NBC target: source at divertor target
plasma_flow += integral(SOL,    -dofSnn   * tf(Nn) );            // -(-n*nn*Riz): sink --> ionization of neutrals

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#endif