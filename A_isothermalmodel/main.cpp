/*
    The following script implements the 1D ISOTHERMAL FLUID EQUATIONS.

    ---------------------------------------------
                PLASMA PARAMETERS               
    ---------------------------------------------
    <SYMBOL>    <TERMINOLOGY>       <UNIT>              <DESCRIPTION>
    ---------------------------------------------
    n           number density      particles/m3        ions/electrons density
    v           velocity            m/s                 plasma velocity
    T           Temperature         eV or K             T[eV] = Kb * T[K]
    p           pressure            N/m2                plasma pressure
    nv          flux density        particles/m2s       particles per unit area per unit time
    ---------------------------------------------

    p = n*T     , where T is in [eV]
      = n*Kb*T  , where T is in [Kelvin] and Kb is Boltzmann constant

    For isothermal compression,
    grad(p) = grad(nT) = T*grad(n)
*/

#include "sparselizard.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace sl;

/*--------------------- GLOBAL CONSTANTS ---------------------*/

// Stabilization methods
const bool SUPG = true;
const bool PSPG = true;
const bool LSIC = true;

// Source model: Not self-consistent: C is a constraint determined analytically.
const int Sp_A_LINEAR    = 1;  // Sp = C*n, where C = (c_s/L)*(pi/2 - 1)
const int Sp_B_CONSTANT  = 2;  // Sp = C
const int Sp_C_NONLINEAR = 3;  // Sp = C*(1-M²)/(2M(1+M²)) // NOT IMPLEMENTED YET

/* ################################################################################################################ */

int main(void)
{
    /*------------------------------------------ User Input Starts ------------------------------------------*/
    const int Sp_TYPE = Sp_B_CONSTANT;

    /* MESH */
    const std::string meshfile = "1D_SOL_geometry_2500Els.msh";   // mesh file
    int SOL=1, upstream=2, downstream=3;                          // physical regions
    mesh mymesh(meshfile);                                        // loading the mesh

    /* INPUT PARAMETERS */
    // Length of SOL [m]
    const double L = 40.0;  

    // ISOTHERMAL TEMPERATURE of SOL [eV]. T = Te = Ti
    const double T = 25.0;                        

    // SPEED OF LIGHT [m/s]
    const double c = 3e+8;

    // MASS ENERGY EQUIVALENT (E = mc²) [eV]
    const double Emc2_e =   0.511 * 1e6;            // m_e*c_l*c_l =   0.511 [MeV]
    const double Emc2_i = 938.272 * 1e6;            // m_i*c_l*c_l = 938.272 [MeV]

    // MASS [eV/c²]
    const double m_e =   0.511e+6 / (c*c);          // Mass of electron                    [eV/(m/s)^2]
    const double m_i = 938.272e+6 / (c*c);          // Mass of ion (proton in this case)   [eV/(m/s)^2]
    const double m   = m_i;                         // Mass of plasma = m_i + m_e ~= m_i

    // ISOTHERMAL ACOUSTIC SPEED [m/s]
    const double c_s = sqrt(2*T/m_i);               // c_s² = k*(Te + Ti) / (m_e + m_i) ≃ 2*T / m_i , where T[eV] = k*T[K]

    // PLASMA DENSITY UPSTREAM of SOL (at x=0)
    const double n_o = 1e+19;

    // VELOCITIES
    double M_start = 0.01;                // Initial ramp-up velocity
    const double rampup_ratio = 1.1;      // M_{k+1} / M_{k} : geometric progression ramp-up
    double M_next = M_start;

    /*------------------------------------------ User Input Ends ------------------------------------------*/

    /* FIELDS */
    field N("h1"), M("h1");   // normalized particle density and velocity (Mach number)

    /* ELEMENT INTERPOLATION ORDER */
    int N_order = 2;
    int M_order = 2;

    // SET THE INTERPOLATION ORDER
    N.setorder(SOL, N_order);
    M.setorder(SOL, M_order);

    // INITIAL CONDITIONS
    N.setvalue(SOL, 1.0);
    M.setvalue(SOL, M_start);

    // BOUNDARY CONDITIONS
    N.setconstraint(upstream, 1.0);
    M.setconstraint(downstream, M_start);

    // NORMALIZATION OF PARTICLE SOURCE
    const double Sp_o = (n_o*c_s/L);
  
    // PARTICLE SOURCE Spn [Normalized]: 
    // See Jupyter notebook for the constrained equation 'C'
    expression C, Sp, Spn;
    if (Sp_TYPE == Sp_A_LINEAR){
        C   = (c_s/L) * (getpi()/2 - 1);
        Sp  = C * (n_o*N); // Sp = C*n = C*(n_o*N)
        Spn = Sp / Sp_o;
    }
    else if (Sp_TYPE == Sp_B_CONSTANT){
        C   = (1.0/2.0) * Sp_o;
        Sp  = C;
        Spn = Sp / Sp_o;    // = 1.0/2.0 =  0.5
    }
    else if (Sp_TYPE == Sp_C_NONLINEAR){
        C   = getpi()/4 * Sp_o;
        Sp  = C * (1-M*M)/(1+M*M);
        Spn = Sp / Sp_o;
    }

    /*-------------------------------------- NORMALIZED FORMULATIONS ------------------------------------*/
    formulation plasma_flow;

    /* MASS CONSERVATION : d(NM)/dX - Spn = 0 */
    plasma_flow += integral(SOL, (dof(N)*grad(M) + M*grad(dof(N))
                               +  N*grad(dof(M)) + dof(M)*grad(N)
                               -  N*grad(M) - M*grad(N)) * tf(N) );  // Linearization of 1st term: Divergence term
    plasma_flow += integral(SOL, -Spn * tf(N));                      // 2nd term: Source term

    /* MOMENTUM CONSERVATION : NM*dM/dX + dN/dX + M*Spn = 0 */
    plasma_flow += integral(SOL, (  dof(N) * grad(M)      * M
                                    + N      * grad(M)      * dof(M)
                                    + N      * grad(dof(M)) * M
                                    - 2*N    * grad(M)      * M     ) * tf(M)    
    ); // Linearization of 1st term
    plasma_flow += integral(SOL, grad(dof(N)) * tf(M));   // 2nd term
    plasma_flow += integral(SOL, dof(M)*Spn * tf(M));     // 3rd term


    /*------------------------------------------ STABILIZATION ------------------------------------------*/

    const int problemdimension = mymesh.getdimension();
    expression h = pow(meshsize(2), 1.0 / problemdimension);    // characteristic element length

    // stabilization parameter
    expression tau = 1.0 / sqrt(pow(2*norm(M)/h, 2));

    expression Rm_supg = N*(grad(M)*dof(M) + grad(dof(M))*M - grad(M)*M) + grad(N) + dof(M)*Spn;
    expression Rm_pspg = dof(N)*grad(M)*M + grad(dof(N)) + M*Spn;
    expression Rc_lsic = N*div(dof(M)) + dof(M)*grad(N) - Spn;

    expression Pw_supg = grad(tf(M))*M;
    expression Pq_pspg = grad(tf(N));
    expression Pw_lsic = div(tf(M));

    if (SUPG == true){ plasma_flow += integral(SOL, Pw_supg * tau * Rm_supg); }
    if (LSIC == true){ plasma_flow += integral(SOL, Pw_lsic * tau * Rc_lsic); }
    if (PSPG == true){ plasma_flow += integral(SOL, Pq_pspg * tau * Rm_pspg); }
    
    /*------------------------------------------ NEWTON ITERATION ------------------------------------------*/
    
    /* SOLVER CONFIGURATION */

    // relaxation factor for coupled velocity-density
    const double w = 1.0;

    const int max_iterations = 250;

    const double tol_relres = 1e-5;

    // initial configuration
    int iter = 1;               // iteration 
    bool is_converged = false;
    bool ramped_up = false;     // status of velocity ramp-up

    // solution vectors
    vec sol_curr(plasma_flow), sol_prev(plasma_flow), sol_pred(plasma_flow);

    // variables for convergence check
    double relres = 0.0;

    // mass conservation variables
    double N_upstream    = 0.0, N_downstream    = 0.0;
    double M_upstream    = 0.0, M_downstream    = 0.0;
    double flux_upstream = 0.0, flux_downstream = 0.0, flux_balance = 0.0;

    // log file for convergence evolution data
    ofstream log_file;
    log_file.open("convergence_log.csv", ios::out | ios::trunc);
    if (log_file.is_open()){
        log_file << "#iter"           << "," \
                 << "N_upstream"      << "," \
                 << "M_upstream"      << "," \
                 << "N_downstream"    << "," \
                 << "M_downstream"    << "," \
                 << "flux_upstream"   << "," \
                 << "flux_downstream" << "," \
                 << "flux_balance"    << "," \
                 << "relres"          << "," \
                 << "is_converged"    << std::endl;
    }
    else{
        std::cout << "Unable to open the log file" << std::endl;
    }

    std::cout << "Solution loop started" << std::endl;

    /* NEWTON SOLVER LOOP */
    while ( is_converged==false && iter<max_iterations )
    {
        // Solve plasma_flow formulation:
        plasma_flow.solve();
        sol_pred.setdata();                             // prediction
        sol_curr = sol_pred*w + sol_prev*(1-w);         // relaxation
        N.setdata(SOL, sol_curr);                       // set field 'N'
        M.setdata(SOL, sol_curr);                       // set field 'M'

        // calculate residual
        plasma_flow.generate();
        mat A = plasma_flow.A();
        vec b = plasma_flow.b();
        relres = (A*sol_curr - b).norm() / b.norm();


        // Check convergence only after rampup is complete
        if (ramped_up == true && iter>100)
        {
            is_converged = (relres < tol_relres);
        }

        // Output the upstream and downstream flux at each iteration
        N_upstream      = (N).integrate(upstream, 2);
        M_upstream      = (M).integrate(upstream, 2);
        flux_upstream   = (N*M * abs(normal(SOL)) ).integrate(upstream, 2);
        N_downstream    = (N).integrate(downstream, 2);
        M_downstream    = (M).integrate(downstream, 2);
        flux_downstream = (N*M * abs(normal(SOL)) ).integrate(downstream, 2);

        flux_balance    = flux_upstream + Spn.integrate(SOL,2);   // (nv)_downstream = (nv)_upstream + Sp*L;

        // Append the convergence evolution data to the log file
        if (log_file.is_open()){
            log_file << std::fixed      << std::setprecision(0) << iter            << "," \
                     << std::fixed      << std::setprecision(1) << N_upstream      << "," \
                     << std::scientific << std::setprecision(6) << M_upstream      << "," \
                     << std::fixed      << std::setprecision(6) << N_downstream    << "," \
                     << std::fixed      << std::setprecision(1) << M_downstream    << "," \
                     << std::scientific << std::setprecision(6) << flux_upstream   << "," \
                     << std::fixed      << std::setprecision(6) << flux_downstream << "," \
                     << std::fixed      << std::setprecision(6) << flux_balance    << "," \
                     << std::scientific << std::setprecision(6) << relres          << "," \
                     << std::boolalpha  << is_converged << std::endl;
        }
        else {
            std::cout << "Unable to open the log file" << std::endl;
        }

        // Update for next iteration        
        iter++; 
        if (ramped_up == false) {
            M_next *= rampup_ratio;
            
            if (M_next < 1.0) {
                M.setconstraint(downstream, M_next);
                is_converged = false; // avoid convergence before complete_rampup
            }
            else {
                M.setconstraint(downstream, 1.0);
                ramped_up = true;
            }
        }
        sol_prev = sol_curr.copy();

    } // solver loop

    // OUTPUT THE FIELDS
    N.write(SOL, "N.vtu", 1);
    M.write(SOL, "M.vtu", 1);

    if (is_converged){
        std::cout << "Solution converged at " << iter-1 << " iterations to complete." << std::endl;
    }
    else {
        std::cout << "Solution NOT converged. Solver stopped as maximum iterations was reached." << std::endl;
    }

    log_file.close();


    /* 
    VALIDATION: Expected values 
    if Sp_TYPE == Sp_A_LINEAR
            M_upstream = 7.53803e-05
            N_downstream = 0.500038

    if Sp_TYPE == Sp_B_CONSTANT
            M_upstream = 5.8534e-05
            N_downstream = 0.500059

    if Sp_TYPE == Sp_C_NONLINEAR
            M_upstream = -0.000103198
            N_downstream = 0.5
    */

    std::cout << "Validation: " << std::endl;    
    std::cout << ( (abs(M_upstream)<2e-4) && (N_downstream > 0.4999) && (N_downstream < 0.5001) ) << std::endl;

    return 0;    

} // main