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
const double Xpt = L/2.0;                           // [m] Length of X-point from midplane(upstream)
std::cout << "\nLength of SOL = " << L << " m" << std::endl;

// CONTROL PARAMETERS
double n_o  = 1.0e+20;          // upstream plasma density. Also the normalizing quantity for plasma density
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
const bool SUPG = true;
const bool PSPG = true;
const bool LSIC = true;

std::string SUPG_TYPE = "CODINA";   // CODINA is slower. With SANJAY runs till q_sol 3.5e+7

const bool useAMR = true;           // Adaptive Mesh Refinement (Recommended)
const bool sourceProjection = true;
const bool noNeutrals = true;
const bool excitation = false;  // Togo et al.
const bool heat_conduction = false; 

// RECYCLING COEFFICIENT
const double R = 0.9;                        // [-]

// TOTAL SHEATH HEAT TRANSMISSION COEFFICIENT
const double gamma_se = 6.0;

// RAMP-UP
const double   rampup_ratio = 1.1;      // M_{k+1} / M_{k} : geometric progression ramp-up
const double rampdown_ratio = 1.0/1.1;  // q_{k+1} / q_{k} : geometric progression ramp-down 
double M_start=0.005, M_next=M_start, M_end=1.0;     // M_end = 1.0 --> v = c_s at target   # seems to help in convergence
double N_start=0.005, N_next=N_start, N_end=0.09236905993134277;    // N_end = 1.0 --> n = n_o at upstream # seems to help in convergence
double q_start=qsol, q_next=q_start, q_end=qsol;


/*
Nend=0.11377807271204562, fluxsource=3.2e+22*Xpt=4e23 for qcond
Nend=0.09236905993134277, fluxsource=3.2e+22*Xpt=4e23 for qconv
*/
/* ------------------------------------------------------------------------------------------------------ */

/*------------------------------------------ User Input Ends ------------------------------------------*/

/* FIELDS */
field N("h1"), M("h1");   // normalized particle density and velocity(Mach number)
field T("h1");            // Temperature
field x("x");             // x-coordinates of the calculation domain

/* ELEMENT INTERPOLATION ORDER */
const int N_order = 2;
const int M_order = 2;
const int T_order = 2; // with 2, better energy balance but hgh computational time. 1 is better.

// SET THE INTERPOLATION ORDER
N.setorder(SOL, N_order);
M.setorder(SOL, M_order);
T.setorder(SOL, T_order);

parameter Tt, q_u;          // Target temperature and q_upstream parameter
Tt.setvalue(SOL, Tti);      // Do not set this zero
q_u.setvalue(SOL, q_next);

// INITIAL CONDITIONS
N.setvalue(SOL, N_start);
M.setvalue(SOL, M_start);
T.setvalue(SOL, Tti);

// BOUNDARY CONDITIONS
N.setconstraint(upstream, N_start);
M.setconstraint(target, M_start);
T.setconstraint(target, Tt); // intial constraint

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

// 1. recombination rate coefficient, sigVrec
expression Tcr(T-1.0, T, 1.0);    // if (T<1eV) --> T = 1eV;
expression E = n_o*N * 1e-14;
double reccoeffs[9][9]= {
    {-2.855728479302E+01, -7.664042607917E-01, -4.930424003280E-03, -5.386830982777E-03, -1.626039237665E-04,  6.080907650243E-06,  2.101102051942E-05, -2.770717597683E-06,  1.038235939800E-07,}, 
    { 3.488563234375E-02, -3.583233366133E-03, -3.620245352252E-03, -9.532840484460E-04,  1.888048628708E-04, -1.014890683861E-05,  2.245676563601E-05, -4.695982369246E-06,  2.523166611507E-07,}, 
    {-2.799644392058E-02, -7.452514292790E-03,  6.958711963182E-03,  4.631753807534E-04,  1.288577690147E-04, -1.145028889459E-04, -2.245624273814E-06,  3.250878872873E-06, -2.145390398476E-07,},
    { 1.209545317879E-02,  2.709299760454E-03, -2.139257298118E-03, -5.371179699661E-04, -1.634580516353E-05,  5.942193980802E-05, -2.944873763540E-06, -9.387290785993E-07,  7.381435237585E-08,},
    {-2.436630799820E-03, -7.745129766167E-04,  4.603883706734E-04,  1.543350502150E-04, -9.601036952725E-06, -1.211851723717E-05,  1.002105099354E-06,  1.392391630459E-07, -1.299713684966E-08,},
    { 2.837893719800E-04,  1.142444698207E-04, -5.991636837395E-05, -2.257565836876E-05,  3.425262385387E-06,  1.118965496365E-06, -1.291320799814E-07, -1.139093288575E-08,  1.265189576423E-09,},
    {-1.886511169084E-05, -9.382783518064E-06,  4.729262545726E-06,  1.730782954588E-06, -4.077019941998E-07, -4.275321573501E-08,  7.786155463269E-09,  5.178505597480E-10, -6.854203970018E-11,},
    { 6.752155602894E-07,  3.902800099653E-07, -1.993485395689E-07, -6.618240780594E-08,  2.042041097083E-08,  3.708616111085E-10, -2.441127783437E-10, -9.452402157390E-12,  1.836615031798E-12,},
    {-1.005893858779E-08, -6.387411585521E-09,  3.352589865190E-09,  1.013364275013E-09, -3.707977721109E-10,  7.068450112690E-12,  3.773208484020E-12, -4.672724022059E-14, -1.640492364811E-14,},
};
expression ln_sigVrec = 0.0;
for (int n=0; n<=8; n++){
    for (int m=0; m<=8; m++){
        ln_sigVrec = ln_sigVrec + (reccoeffs[n][m] * pow(log(Tcr),n) * pow(log(E),m));   // eirene.de/amjuel.pdf: Example of Use of Fits
    }
}

// 2. excitation rate coefficient, sigVex
expression Yex = 10.2/Tcr;
expression sigVex = (49.0e-14/(0.28+Yex)) * exp(-Yex) * sqrt(Yex*(1.0+Yex));

// [cm³/s] * 1e-6 --> [m³/s]
expression sigVrec = ifpositive(N*n_o-1e3, exp(ln_sigVrec) * 1e-6, 0.0);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// NORMALZIATION OF PARTICLE SOURCE
expression NVo = n_o*c_st;
std::cout << "NVo = " << NVo.evaluate() << std::endl;

// LOCAL IONIZATION OF NEUTRALS: S = n_a*n_b*<σv>
expression Rrc = n_o*n_o*N*N*sigVrec;
expression Rrcn = Rrc / NVo;

expression Rex = n_o*n_o*N*N*sigVex;
expression Rexn = Rex / NVo;

// Energy lost per ionization [eV]
expression Eiz = expression(excitation-1, 13.6, 30.0);

//++++++++++++
const double particleflux = 3.2e+22 * Xpt; //3.2e+22*Xpt=4e23;   // [1/(m²s)]
expression particlesource = (particleflux / Xpt) * H(Xpt - x);  // [1/(m³s)]

const double energyflux = qsol;     // [W/m²]
expression   energysource = (energyflux / Xpt) * H(Xpt - x);    // [J/(m³s)] or [W/m³]
// expression pressuresource = (2.0/3.0)*energysource / e;       // [eV/(m³s)]

// projection
field NeSource("h1"); NeSource.setorder(SOL, 2);
field EeSource("h1"); EeSource.setorder(SOL, 2);

if (sourceProjection){
    formulation projection;
    projection += integral(SOL, (dof(NeSource) - particlesource)*tf(NeSource));
    projection += integral(SOL, (dof(EeSource) - energysource)*tf(EeSource));
    projection.allsolve(1e-6, 1000, 1e-6, -1);
}
else{
    NeSource.setvalue(SOL, particlesource);
    EeSource.setvalue(SOL, energysource);
}
//++++++++++++

// SOURCE [normalized]:
expression Sin_internal = -Rrcn;     // ionization of neutrals are a source for plasma, hence positive source
expression Sin_external = NeSource / NVo;
expression Sin = Sin_internal + Sin_external;
// expression Sin = Sin_external;

expression    Smn_internal = -M*Rrcn;
expression dofSmn_internal = -dof(M)*Rrcn;
expression    Smn = Smn_internal;
expression dofSmn = dofSmn_internal;

expression En, Rn; 
En = (3.0/2.0) * T  *Rrcn;          // transfer of energy to neutrals
Rn = (1.09*T - 13.6)*Rrcn + Rexn;   // plasma cooling due to radiation
expression Spn_internal = (-En - Rn) * e;       // [J/(m³s)] or [W/m³]
expression Spn_external = EeSource;             // [J/(m³s)] or [W/m³]
expression Spn = Spn_internal + Spn_external;   // [J/(m³s)] or [W/m³]
// expression Spn = Spn_external;

// NeSource.write(SOL, "NeSource.vtu", 2);
// EeSource.write(SOL, "EeSource.vtu", 2);
// std::exit(-1);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// FORMULATIONS
formulation plasma_flow;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#endif