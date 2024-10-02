#ifndef PREPROCESS_1D_SOL_H
#define PREPROCESS_1D_SOL_H

#include "sparselizard.h"

using namespace sl;

/*-------------------------------------------------CONSTANTS-------------------------------------------------*/
// SPEED OF LIGHT [m/s]
const double c = 3e+8;

// MASS ENERGY EQUIVALENT (E = mc²) [eV]
const double Emc2_e =   0.511 * 1e6;            // m_e*c_l*c_l =   0.511 [MeV]
const double Emc2_i = 938.272 * 1e6;            // m_i*c_l*c_l = 938.272 [MeV]

double e = 1.6022e-19;

double gamma = 1.0;
/* -------------------------------------------------------------------------------------------------------- */

// LENGTH of SOL [m] & corresponding FILENAME
std::string filename = "1D_SOL_geometry_L25.msh";
// std::string filename = "1D_SOL_geometry_L25_nonuniform.msh";

// MESH    
std::string meshpath = "/home/rahulkn/AllSimcore/AllCADgeometries/Plasma1D_mesh/";
std::string meshfile = meshpath + filename;
int SOL=1, upstream=2, target=3;    // physical regions
mesh mymesh(meshfile);      // loading the mesh

// LENGTH of SOL [m]
const double L = mymesh.getdimensions()[0];     // [m] Domain length or Length of SOL
const double Xpt = L;                           // [m] Length of X-point from midplane(upstream)
std::cout << "\nLength of SOL = " << L << " m" << std::endl;

// CONTROL PARAMETERS
double n_o  = 1.0e+20;          // upstream plasma density. Also the normalizing quantity for n
double qsol = 2e7;  // 1e+9     // Heat flux W/m² entering the SOL at upstream (NBC)

// MASS [eV/c²]
const double m_e =       0.511e+6 / (c*c);      // Mass of electron                    [eV/(m/s)^2]
const double m_i = 2 * 938.272e+6 / (c*c);      // Mass of ion (proton in this case)   [eV/(m/s)^2]
const double m   = m_i;                         // Mass of plasma = m_i + m_e ~= m_i
/* -------------------------------------------------------------------------------------------------------- */

// INITIAL TEMPERATURES (depends on qsol)
double Ti  = 100.0;     // initial SOL temperature
double Tti = 75.0;     // initial target temperatuer. Do not set this to zero.

// STABILIZATION METHODS
const bool SUPG = false;
const bool PSPG = true;
const bool LSIC = true;
const bool SUPGT = true;

std::string SUPG_TYPE = "CODINA";   // CODINA is slower. With SANJAY runs till q_sol 3.5e+7

const bool verify2ptmodel= false;    // Basic 2-point model verification
const bool useAMR = true;           // Adaptive Mesh Refinement (Recommended)
const bool excitation = false;  // Togo et al.
const bool updatedradiatedpower = true;
const bool noNeutrals = true;
const std::string diffcoeff_type = "SD1D"; // UEDGE or SD1D

// RECYCLING COEFFICIENT
const double Frecycle = 0.2;                        // [-]

// RAMP-UP
const double   rampup_ratio = 1.1;      // M_{k+1} / M_{k} : geometric progression ramp-up
const double rampdown_ratio = 1.0/1.1;  // q_{k+1} / q_{k} : geometric progression ramp-down 
double M_start=0.01, M_next=M_start, M_end=1.0;     // M_end = 1.0 --> v = c_s at target   # seems to help in convergence
double N_start=0.01, N_next=N_start, N_end=0.15514151716039512;    // N_end = 1.0 --> n = n_o at upstream # seems to help in convergence
double q_start=qsol, q_next=q_start, q_end=qsol;


/*
Nend=0.15202712705826016 --> diffusiveneutral, Xpt = L/2
Nend=0.15514151716039512 --> diffusiveneutral, Xpt = L
*/
/* ------------------------------------------------------------------------------------------------------ */

/*------------------------------------------ User Input Ends ------------------------------------------*/

/* FIELDS */
field N("h1"), M("h1");   // normalized particle density and velocity(Mach number)
field T("h1");            // Temperature
field x("x");             // x-coordinates of the calculation domain
field Nn("h1");

/* ELEMENT INTERPOLATION ORDER */
const int N_order = 2;
const int M_order = 2;
const int T_order = 2; // with 2, better energy balance but hgh computational time. 1 is better.
const int Nn_order = 2;

// SET THE INTERPOLATION ORDER
N.setorder(SOL, N_order);
M.setorder(SOL, M_order);
T.setorder(SOL, T_order);
Nn.setorder(SOL, Nn_order);

parameter Tt, q_u;          // Target temperature and q_upstream parameter
Tt.setvalue(SOL, Tti);      // Do not set this zero
q_u.setvalue(SOL, q_next);

// INITIAL CONDITIONS
N.setvalue(SOL, N_start);
M.setvalue(SOL, M_start);
T.setvalue(SOL, Tti);
Nn.setvalue(SOL, 1e-20); //??

// BOUNDARY CONDITIONS
N.setconstraint(upstream, N_start);
M.setconstraint(target, M_start);
T.setconstraint(target, Tt); // intial constraint
Nn.setconstraint(upstream, 1e-20);

// MAX ACOUSTIC SPEED [m/s]
parameter c_st;
c_st.setvalue(SOL, sqrt(gamma * 2.0*Tt/m_i));   // c_s² = k*(Te + Ti) / (m_e + m_i) ≃ 2*T / m_i , where T[eV] = k*T[K]


// AMR: Adaptation criterion for hp-adaptivity:
expression adaptcriterion = norm(grad(M));
if(useAMR){
    // Mesh adaptivity up to 5 refinement levels:
    mymesh.setadaptivity(adaptcriterion, 0, 5);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/* RATE COEFFICIENTS <σv> [m³/s] --> sigV*/
#include "ratecoefficients.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// NORMALZIATION OF PARTICLE SOURCE
expression NVo = n_o*c_st;
std::cout << "NVo = " << NVo.evaluate() << std::endl;

// LOCAL IONIZATION OF NEUTRALS: S = n_a*n_b*<σv>
expression Rrcn = n_o*n_o*N*N*sigVrc / NVo;
expression Rizn = n_o*n_o*N*Nn*sigViz / NVo;
expression Rcxn = n_o*n_o*N*Nn*sigVcx / NVo;
expression Rexn = n_o*n_o*N*Nn*sigVex / NVo;

expression dofRizn = n_o*n_o*N*dof(Nn)*sigViz / NVo;

// Energy lost per ionization [eV]
expression Eiz = expression(excitation-1, 13.6, 30.0);

//++++++++++++
const double particleflux = 1.6e+22 * L; //3.2e+22*L/2=4e23; 1.6e+22*L=4e23;  // [1/(m²s)]
expression particlesource = (particleflux / Xpt) * H(Xpt - x);  // [1/(m³s)]

const double energyflux = qsol;     // [W/m²]
expression   energysource = (energyflux / Xpt) * H(Xpt - x);;  // [J/(m³s)] or [W/m³]
// expression pressuresource = (2.0/3.0)*energysource / e;       // [eV/(m³s)]
//++++++++++++

// SOURCE [normalized]:
expression Sin_internal = Rizn - Rrcn;     // ionization of neutrals are a source for plasma, hence positive source
expression Sin_external = particlesource / NVo;
expression Sin = Sin_internal + Sin_external;

expression Snn_internal = -Rizn + Rrcn;
expression dofSnn_internal = -dofRizn + Rrcn;
expression Snn_external = 0.0;
expression Snn = Snn_internal + Snn_external;
expression dofSnn = dofSnn_internal + Snn_external;

expression    Smn_internal = -M*Rrcn - M*Rcxn;
expression dofSmn_internal = -dof(M)*Rrcn - dof(M)*Rrcn;
expression    Smn = Smn_internal;
expression dofSmn = dofSmn_internal;

expression En, Rn; 
En = (3.0/2.0)*T*Rrcn + (3.0/2.0)*T*Rizn + (3.0/2.0)*T*Rcxn;          // transfer of energy to neutrals
Rn = (1.09*T - 13.6)*Rrcn + Eiz*Rizn + Rexn;   // plasma cooling due to radiation
expression Spn_internal = (-En - Rn) * e;   // [J/(m³s)] or [W/m³]
expression Spn_external = energysource;     // [J/(m³s)]
expression Spn = Spn_internal + Spn_external;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// FORMULATIONS
formulation plasma_flow;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#endif