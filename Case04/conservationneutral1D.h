#ifndef CONSERVATION_NEUTRAL_1D_H
#define CONSERVATION_NEUTRAL_1D_H

#include "sparselizard.h"

using namespace sl;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DIFFUSION OF NEUTRALS:
∇.(-D.∇nn) - Sn = 0
NBC at target ->  -D.∇nn = Γ_recycling = -R*(nv)

Normalized with NVo=n_o * c_st:
∇.((-D/c_st).∇Nn) - Snn = 0
NBC at target ->  (-D/c_st).∇Nn = -R*(NM)

As described in UEDGE paper   : https://doi.org/10.1063/1.873488
----------------------------------------------*/

// charge exchange frequency, nucx [1/s]
expression nucx = (n_o*N)  * sigVcx;

// ionization frequency, nuiz [1/s]
expression nuiz = (n_o*N)  * sigViz;

// neutral-neutral collision frequency, nunn [1/s]
expression vth_n = sqrt(T/m_i);                             // thermal velocity of neutral atom [m/s]
expression a0 = getpi() * pow(5.29e-11, 2);                 // cross-sectional area of a neutral hydrogen atom [m²]
expression Lmax = 0.1;                                      // Max neutral-neutral mean free path
expression lambda_nn = ifpositive(1.0/(a0*(n_o*Nn)) - Lmax, 
                            Lmax, 1.0/(a0*(n_o*Nn)));       // Neutral-neutral mean free path
expression nunn = vth_n / lambda_nn;

// Total neutral collision frequency
expression nu;

if(diffcoeff_type == "UEDGE")
    nu = nucx + nuiz;

if(diffcoeff_type == "SD1D")
    nu = nucx + nuiz + nunn;

// Diffusion coefficient, D [m²/s]
expression dneut = 10.0;
expression D = dneut * pow(vth_n,2) / nu;

// Normalized diffusion coefficient
expression Dn = D/c_st;

// WEAK FORMULATION 
plasma_flow += integral(SOL,     Dn*grad(dof(Nn))*grad(tf(Nn)), 0, 3);    // integ. by parts: term 1
plasma_flow += integral(target,  Frecycle*N*(-M) * tf(Nn), 0, 3);            // integ. by parts: term 2 -> NBC target: source at divertor target
plasma_flow += integral(SOL,    -dofSnn   * tf(Nn), 0, 3);            // -(-n*nn*sigViz): sink --> ionization of neutrals


// Nn.setconstraint(SOL, 1e-20);
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#endif