#ifndef PREPROCESS_1D_SOL_H
#define PREPROCESS_1D_SOL_H

/* ------------------------------------------------------------------------------------------------------------
ASSUMPTIONS:


------------------------------------------------------------------------------------------------------------ */

#include "sparselizard.h"

using namespace sl;

/*-------------------------------------------------CONSTANTS-------------------------------------------------*/
// SPEED OF LIGHT [m/s]
const double c = 3e+8;

// MASS ENERGY EQUIVALENT (E = mc²) [eV]
const double Emc2_e =   0.511 * 1e6;            // m_e*c_l*c_l =   0.511 [MeV]
const double Emc2_i = 938.272 * 1e6;            // m_i*c_l*c_l = 938.272 [MeV]
/* -------------------------------------------------------------------------------------------------------- */

// LENGTH of SOL [m] & corresponding FILENAME
const double L = 100.0;

// MESH
std::string meshfile = "1D_SOL_geometry_L100.msh";
int SOL=1, upstream=2, target=3;    // physical regions
mesh mymesh(meshfile);              // loading the mesh

// CONTROL PARAMETERS
double n_o  = 3.0e+19;   // upstream plasma density. Also the normalizing quantity for n, nn
double qsol = 2.5e+7;       // Heat flux W/m² entering the SOL at upstream (NBC)

// MASS [eV/c²]
const double m_e =   0.511e+6 / (c*c);          // Mass of electron                    [eV/(m/s)^2]
const double m_i = 2 * 938.272e+6 / (c*c);      // Mass of ion (proton in this case)   [eV/(m/s)^2]
const double m   = m_i;                         // Mass of plasma = m_i + m_e ~= m_i
/* -------------------------------------------------------------------------------------------------------- */

// INITIAL TEMPERATURES (depends on qsol)
double Ti  = 60.0;     // initial SOL temperature
double Tti = 20.0;     // initial target temperatuer. Do not set this to zero.

// STABILIZATION METHODS
const bool SUPG = true;
const bool PSPG = true;
const bool LSIC = true;

std::string SUPG_TYPE = "CODINA";   // CODINA is slower. With SANJAY runs till q_sol 3.5e+7

const bool verify2ptmodel= true;    // Basic 2-point model verification
const bool useAMR = true;           // Adaptive Mesh Refinement (Recommended)

// RECYCLING COEFFICIENT
const double R = 0.99999;                        // [-]

// FITTING COEFFICIENTS <σv> [m²/s] (from eirene.de/amjuel.pdf)
// Electron ionization: Section 2.17:  Reaction 2.1.5FJ: e + H --> H⁺ + 2e
std::vector<double> b_ionization = {-0.317385000000e+02,  0.114381800000e+02, -0.383399800000e+01,
                                     0.704669200000e+00, -0.743148620000e-01,  0.415374900000e-02,
                                    -0.948696700000e-04,  0.000000000000e-00,  0.000000000000e+00
};
// Charge exchange: Section 2.19: Reaction 3.1.8: H⁰ + p --> p + H⁰
std::vector<double> b_chargeexchange = {-1.850280000000E+01,  3.708409000000E-01,  7.949876000000E-03,
                                        -6.143769000000E-04, -4.698969000000E-04, -4.096807000000E-04,
                                         1.440382000000E-04, -1.514243000000E-05,  5.122435000000E-07
};

// RAMP-UP
const double   rampup_ratio = 1.1;      // M_{k+1} / M_{k} : geometric progression ramp-up
const double rampdown_ratio = 1.0/1.1;  // q_{k+1} / q_{k} : geometric progression ramp-down 
double M_start=0.01, M_next=M_start, M_end=1.0;    // M_end = 1.0 --> v = c_s at target   # seems to help in convergence
double N_start=0.01, N_next=N_start, N_end=1.0;    // N_end = 1.0 --> n = n_o at upstream # seems to help in convergence
double q_start=qsol, q_next=q_start, q_end=qsol;
/* ------------------------------------------------------------------------------------------------------ */

/*------------------------------------------ User Input Ends ------------------------------------------*/

/* FIELDS */
field N("h1"), M("h1");   // normalized particle density and velocity(Mach number)
field Nn("h1");           // normalized neutral density
field T("h1");            // Temperature
field x("x");             // x-coordinates of the calculation domain

/* ELEMENT INTERPOLATION ORDER */
const int N_order  = 1;
const int M_order  = 1;
const int Nn_order = 1;
const int T_order  = 1; // with 2, better energy balance but hgh computational time. 1 is better.

// SET THE INTERPOLATION ORDER
N.setorder(SOL,   N_order);
M.setorder(SOL,   M_order);
Nn.setorder(SOL, Nn_order);
T.setorder(SOL,   T_order);

parameter Tt, q_u;          // target temperature and q_upstream parameter
Tt.setvalue(SOL, Tti);      // Do not set this zero
q_u.setvalue(SOL, q_next);

// INITIAL CONDITIONS
N.setvalue(SOL, N_start);
M.setvalue(SOL, M_start);
Nn.setvalue(SOL, 1e-20);
T.setvalue(SOL, Ti);

// BOUNDARY CONDITIONS
N.setconstraint(upstream, N_start);
M.setconstraint(target, M_start);
Nn.setconstraint(upstream, 1e-20);
T.setconstraint(target, Tt);  // DBC
        
// MAX ACOUSTIC SPEED [m/s]
parameter c_st;
c_st.setvalue(SOL, sqrt(2.0*Tt/m_i));   // c_s² = k*(Te + Ti) / (m_e + m_i) ≃ 2*T / m_i , where T[eV] = k*T[K]

// AMR: Adaptation criterion for hp-adaptivity:
expression adaptcriterion = norm(grad(M));
if(useAMR){
    // Mesh adaptivity up to 5 refinement levels:
    mymesh.setadaptivity(adaptcriterion, 0, 5);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* RATE COEFFICIENTS <σv> [m³/s] --> Riz, Rcx*/
expression ln_T = sl::log(T);
// 1. ionization rate coefficient, Riz
expression ln_Riz = 0.0;
for (int n=0; n<=b_ionization.size(); n++){
    ln_Riz = ln_Riz + (b_ionization[n] * pow(ln_T, n));             // eirene.de/amjuel.pdf: Example of Use of Fits
}
// 2. charge exchange rate coefficient, Rcx
expression ln_Rcx = 0.0;
for (int n=0; n<=b_chargeexchange.size(); n++){
    ln_Rcx = ln_Rcx + (b_chargeexchange[n] * pow(ln_T, n));         // eirene.de/amjuel.pdf: Example of Use of Fits
}

expression Riz = exp(ln_Riz) * 1e-6;    // [cm³/s] * 1e-6 --> [m³/s]
expression Rcx = exp(ln_Rcx) * 1e-6;    // [cm³/s] * 1e-6 --> [m³/s]

// ln_Riz.print();
// Riz.print();
// std::exit(-1);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// NORMALZIATION OF PARTICLE SOURCE
expression Sp_o = n_o*c_st;

// LOCAL IONIZATION OF NEUTRALS: S = nn*ni*<σv>
expression Siz    = n_o*n_o *     Nn *N*Riz;
expression dofSiz = n_o*n_o * dof(Nn)*N*Riz;
expression Scx    = (verify2ptmodel) ? 0.0 : n_o*n_o*Nn*N*Rcx;

// SOURCE [normalized]:
expression Spn =  Siz / Sp_o;     // ionization of neutrals are a source for plasma, hence positive source
expression Snn = -Siz / Sp_o;     // ionization of neutrals are a sink for neutrals, hence negative source
expression dofSnn = -dofSiz / Sp_o;
expression Scxn = Scx / Sp_o;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// FORMULATIONS
formulation plasma_flow;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#endif