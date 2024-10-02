#include "sparselizard.h"
#include "utilityfunctions.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace sl;

// function alias
const auto relerrnorm = sl::allmeasuredistance;


int main(void){
    wallclock myclock;        // initialize the wall clock object


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #include "preprocess_1D_SOL_case-02.h"

    #include "conservationmassmomentum1D.h"
    #include "conservationenergy1D.h"

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /* NEWTON SOLVER CONFIGURATION: CS=coupled_solver
    --------------------------------------------------*/

    // rampup/rampdown flags
    bool rampedup=false, M_rampedup=false, N_rampedup=false;
    bool q_rampeddown=false;

    // relaxation factor
    double w_CS = 0.50;
    double w_Tt = 0.99; // 0.50 results in NAN/INF

    // Max iterations
    int maxiters_CS = 1200;          // max iterations for coupled solver for a specified Tt

    // Convergence criteria
    double tolrelerr_CS = 1e-6; // 1e-5 for qconv only
    double tolrelerr_Tt = 1e-6; // 1e-5 for qconv only

    // Variables used in convergence check
    double relerr_CS = 1.0;
    double relerr_Tt = 1.0;
    bool isconverged_CS = false;
    bool isconverged_Tt = false;

    // Solver iterations
    int iter_total = 1;
    int iter_CS    = 0;

    // Variables for upstream/target field values
    double  N_upstream=0.0,  N_target=0.0;
    double  M_upstream=0.0,  M_target=0.0;
    double  T_upstream=0.0,  T_target=0.0;

    // Flux balance variables
    double massflux_upstream = 0.0, momnetumflux_upstream = 0.0, neutralflux_upstream = 0.0, energyflux_upstream = 0.0;
    double massflux_target   = 0.0, momnetumflux_target   = 0.0, neutralflux_target   = 0.0, energyflux_target   = 0.0;
    double massflux_balance  = 0.0, momnetumflux_balance  = 0.0, neutralflux_balance  = 0.0, energyflux_balance  = 0.0;

    // log file for convergence evolution data
    ofstream log_file;
    log_file.open("convergence_log.csv", ios::out | ios::trunc);
    if (log_file.is_open()){
        log_file << "q_sol"          << "," \
                 << "CS_converged@"  << "," \
                 << "relerr_CS"      << "," \
                 << "N_upstream"     << "," \
                 << "M_upstream"     << "," \
                 << "T_upstream"     << "," \
                 << "N_target"       << "," \
                 << "M_target"       << "," \
                 << "T_target"       << std::endl;
    }
    else{
        std::cout << "Unable to open the log file" << std::endl;
    }

    // solution vectors (used in coupled sovler loop)
    vec sol_prev(plasma_flow), sol_pred(plasma_flow), sol_curr(plasma_flow);
    double Tt_prev=0.0, Tt_pred=0.0, Tt_curr=0.0;


    bool maxbeta_reached = false;

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /* NETWON SOVLER
    ----------------------*/

    // HEAT FLUX RAMPEDDOWN LOOP
    while (q_rampeddown==false)
    {
        // reset for each new q_sol
        iter_CS = 0;
        isconverged_CS = false;

        // Log output
        std::cout << std::endl << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "Current q_sol = " << q_next << std::endl;
    
        std::cout << "\n\n------------------------------------" << std::endl;
        std::cout << "iter_total:" << iter_total << std::endl;
        std::cout << "------------" << std::endl;

        /*-------------------------------------- COUPLED SOLVER (NONLINEAR) -------------------------------------- */
        // COUPLED SOLVER LOOP (Solves for <N, M, T> for the current target temperature)
        while( isconverged_CS==false && iter_CS<maxiters_CS ){

            plasma_flow.solve("lu", false); 

            // if(iter_CS < 81){
            //     plasma_flow.solve("lu", false, {0,2});                         // solve the coupled plasma formulation
            // }
            // else{
            //     plasma_flow.solve("lu", false, {0,1});  
            // }
            
            /*Floor field values */ // Allow the to evolve freely
            if(rampedup==false){
                M.setvalue(selectall(), max(1e-06, M));
                N.setvalue(selectall(), max(1e-10, N));
                T.setvalue(selectall(), max(5.0, T));
            }
            
            sol_pred.setdata();                                     // prediction
            sol_curr = sol_pred*w_CS + sol_prev*(1-w_CS);           // relaxation
            sl::setdata(sol_curr);                                  // update the fields with the relaxed solution vector
            relerr_CS = relerrnorm(sol_curr, sol_prev, sol_curr);   // relative solution error

            // Stop the simulation if NAN or INF occurs
            if ( std::isnan(relerr_CS) || std::isinf(relerr_CS) ){
                    std::cout << std::endl << "ENCOUNTERED NAN or INF in COUPLED SOLVER. SIMULATION TERMINATED!" << std::endl;
                    std::exit(-1);
            }

            isconverged_CS = relerr_CS < tolrelerr_CS;              // check for convergence of coupled solver
            sol_prev = sol_curr.copy();                             // update for next coupled solver iteration
            iter_CS++;                                              // move to the next coupled solver iteration
        
            std::cout << "------------------------------------" << std::endl;
            std::cout <<   "iter_CS: " <<   iter_CS << ";\t";
            std::cout << "relerr_CS: " << relerr_CS << ";\t";
            std::cout << "isconverged_CS: " << std::boolalpha << isconverged_CS << std::endl;

            /*------------------------------------------ POST-PROCESSING FOR CS------------------------------------------ */
            // Get the field values at upstream and target at each iteration
            M_upstream  =  M.integrate(upstream, 5);
            N_upstream  =  N.integrate(upstream, 5);
            T_upstream  =  T.integrate(upstream, 5);
            M_target    =  M.integrate(target, 5);
            N_target    =  N.integrate(target, 5);
            T_target    =  T.integrate(target, 5);

            double tau_upstream  =  expression(tau).integrate(upstream, 5);
            double tau_target    =  expression(tau).integrate(target, 5);
            double tauT_upstream =  expression(tauT).integrate(upstream, 5);
            double tauT_target   =  expression(tauT).integrate(target, 5);

            std::cout << "M  = [" <<  M_upstream << ", " <<  M_target << "]" << std::endl;
            std::cout << "N  = [" <<  N_upstream << ", " <<  N_target << "]" << std::endl;
            std::cout << "T  = [" <<  T_upstream << ", " <<  T_target << "]" << std::endl;
            std::cout << "tau= [" <<  tau_upstream << ", " <<  tau_target << "]" << std::endl;
            std::cout << "tauT= [" <<  tauT_upstream << ", " <<  tauT_target << "]" << std::endl;

            // Append the convergence evolution data to the log file
            if (log_file.is_open()){
                log_file << std::scientific << std::setprecision(1) << q_next         << "," \
                            << std::scientific << std::setprecision(5) << iter_CS     << "," \
                            << std::scientific << std::setprecision(6) << relerr_CS   << "," \
                            << std::scientific << std::setprecision(1) << N_upstream  << "," \
                            << std::scientific << std::setprecision(5) << M_upstream  << "," \
                            << std::scientific << std::setprecision(6) << T_upstream  << "," \
                            << std::scientific << std::setprecision(5) << N_target    << "," \
                            << std::scientific << std::setprecision(1) << M_target    << "," \
                            << std::scientific << std::setprecision(5) << T_target    << std::endl;
            }
            else {
                std::cout << "Unable to open the log file" << std::endl;
            }


            expression massflux = N*M;
            massflux_upstream = massflux.integrate(upstream, 5);
            massflux_target   = massflux.integrate(target, 5);

            expression momnetumflux = N*M*M + N*T*(1.0/Tt);
            momnetumflux_upstream = momnetumflux.integrate(upstream, 5);
            momnetumflux_target   = momnetumflux.integrate(target, 5);

            // expression energyflux = beta*on(SOL, -K*grad(T)) + alpha*(5*T*N*M*n_o*c_st)*e + alpha*(0.5*m_i*N*M*M*M*n_o*c_st*c_st*c_st)*e;
            expression energyflux = beta*on(SOL, -K*grad(T)) + alpha*((5.0/2.0)*2*T*N*M*n_o*c_st)*e;
            energyflux_upstream = energyflux.integrate(upstream, 5); // expression(q_u).evaluate(); // q_next // imposed //energyflux.integrate(upstream, 5);
            energyflux_target   = energyflux.integrate(target, 5);

            // Flux balances
            massflux_balance     = massflux_upstream     + (Sin).integrate(SOL, 5);
            momnetumflux_balance = momnetumflux_upstream + (Smn).integrate(SOL, 5);
            expression Vdp = e*2*n_o*c_st* M*(N*grad(T) + T*grad(N));
            energyflux_balance   = energyflux_upstream   + (Vdp + Spn).integrate(SOL, 5);

            std::cout << "massbalance     : " << massflux_target     << "\t\t" << massflux_balance     << std::endl;
            std::cout << "momentumbalance : " << momnetumflux_target << "\t\t" << momnetumflux_balance << std::endl;
            std::cout << "energybalance   : " << energyflux_target   << "\t\t" << energyflux_balance   << std::endl;

            // Pressure balance
            double P_upstream = (N*T).integrate(upstream, 5);
            double P_target   = (N*T).integrate(target, 5);
            std::cout << "pressurebalance : " <<  2*P_target << "\t\t" << P_upstream << std::endl;


            /*------------------------------------------ UPDATE FOR NEXT Tt Loop ------------------------------------------ */
            // Update the stabilization parameter (1e-3 helps. You are in trouble if it is changed to 1.0: Revisit again)
            if (SUPG == true){
                double factor = 1e-3;
                if (q_next<2e+8) {factor = 1.0;}
                if (SUPG_TYPE=="SANJAY"){ tau.setvalue(SOL, factor / sqrt(pow(2*norm(M)/h, 2))); }
                if (SUPG_TYPE=="CODINA"){ tau.setvalue(SOL, factor * ifpositive(denom-eps, 1.0 / denom, 0.0));}
            }

            // Ramp-up M & N 
            if (rampedup == false) {
                M_next *= rampup_ratio;
                N_next *= rampup_ratio;

                if (M_next < 1.0) {
                    M.setconstraint(target, M_next); 
                }
                else { 
                    M.setconstraint(target, M_end);
                    M_rampedup = true;
                }

                if (N_next < 1.0) {
                    N.setconstraint(upstream, N_next); 
                }
                else {
                    N.setconstraint(upstream, N_end);
                    N_rampedup = true;
                }

                if(M_rampedup==true && N_rampedup==true){ 
                    rampedup = true;
                    std::cout << "M & N ramp up complete" << std::endl;
                }
                else {
                    isconverged_CS = false;
                    std::cout << "M & N ramping up ..." << std::endl;
                }
            }

            /*-------------------------------------- TARGET TEMPERATURE SOLVER -------------------------------------- */
            // Solve for new target temperature with the latest N and M (only after N and M)
            if(rampedup==true){
                Tt_pred = T.integrate(target, 5);
                Tt_curr = Tt_pred*w_Tt + Tt_prev*(1-w_Tt);

                relerr_Tt = abs( (Tt_curr-Tt_prev)/Tt_curr );
                isconverged_Tt = relerr_Tt < tolrelerr_Tt;
                std::cout << "relerr_Tt: " << relerr_Tt << ";\t";
                std::cout << "isconverged_Tt: " << std::boolalpha << isconverged_Tt << std::endl;

                isconverged_CS = isconverged_CS && isconverged_Tt;
                // if(expression(beta).evaluate()>max_beta) isconverged_CS = false;

                if(isconverged_CS==false){
                    // double gamma_se = 6.0;  // total sheath heat transmission coefficient, γ ≃ 7–8
                    Tt.setvalue(SOL, Tt_curr); // Tt.setvalue(SOL, Tt_curr);
                    c_st.setvalue(SOL, sqrt(gamma *2.0*Tt/m_i));
                    expression GAMMA_se = n_o*c_st * N*M;

                    T.setconstraint(target, qsol/(gamma_se*GAMMA_se*e));
                    std::cout << "New target temperature = " << expression(Tt).evaluate() << std::endl;
                }

                Tt_prev = Tt_curr;

                if (useAMR){
                    std::cout << "meshadaptivity: " << std::boolalpha << adapt() << std::endl; // important for maintaining flux balances.
                }

                // Only if heat conduction is false
                if (heat_conduction==false){
                    double max_beta = 9e-6; // lowest:9e-6; move this to preprocess file
                    if(isconverged_Tt && !maxbeta_reached){
                        // a = a + 0.05;
                        // alpha.setvalue(selectall(), min(a, 1.0));
                        b = b/10.0; if(b-max_beta<1e-6){ b=max_beta; maxbeta_reached=true; } // b = b - 0.1
                        // b = b/7.0; if(b < max_beta){ maxbeta_reached=true; }
                        beta.setvalue(selectall(), max(b, max_beta)); // allows for artificial diffusion, increased stability
                        isconverged_CS = false;
                    }
                    std::cout << "beta = " << expression(beta).evaluate() << std::endl;
                }
            } // new target temperature       

            myclock.print("Total time elapsed so far: ");
        }

        /*------------------------------------------ POST-PROCESSING FOR Tt------------------------------------------ */
        ostringstream qsol_ss;
        qsol_ss << std::scientific << std::setprecision(2)  << q_next;

        //  T.write(SOL,  "T_q" + qsol_ss.str() + ".vtu", 1);
        //  M.write(SOL,  "M_q" + qsol_ss.str() + ".vtu", 1);
        //  N.write(SOL,  "N_q" + qsol_ss.str() + ".vtu", 1);

        /*------------------------------------------ UPDATE FOR NEXT q_sol Loop ------------------------------------------ */
        if (q_rampeddown==false){
            q_next = expression(q_u).evaluate() * rampdown_ratio;
            if(q_next > q_end){
                q_u.setvalue(SOL, q_next);
                q_rampeddown = true;
                break;
            }
            else{
                q_next = q_end;
                q_u.setvalue(SOL, q_next);
                q_rampeddown = true;
            }
        }
        else{
            break;
        }

    } // q_u rampdown

    // OUTPUT THE FIELDS
     T.write(SOL, "T.vtu", 1);
     N.write(SOL, "N.vtu", 1);
     M.write(SOL, "M.vtu", 1);

    expression qcond = beta*on(SOL, -K*grad(T));
    expression qconv = alpha*(5*T*N*M * (n_o*c_st) + 0.5*m_i*N*M*M*M * (n_o*c_st*c_st*c_st))*e;
    expression qtotal = qcond + qconv;
    qcond.write(selectall(), "qcond.vtu", 2);
    qconv.write(selectall(), "qconv.vtu", 2);
    qtotal.write(selectall(), "qtotal.vtu", 2);

    sigVrec.write(selectall(), "sigVrec.vtu", 2);
    sigVex.write(selectall(), "sigVex.vtu", 2);

    Rrcn.write(selectall(), "Rrcn.vtu", 2);
    Rexn.write(selectall(), "Rexn.vtu", 2);

    Sin.write(selectall(), "Sin.vtu", 2);
    Sin_internal.write(selectall(), "SinInt.vtu", 2);
    Sin_external.write(selectall(), "SinExt.vtu", 2);
    Smn.write(selectall(), "Smn.vtu", 2);
    Spn.write(selectall(), "Spn.vtu", 2);
    Spn_internal.write(selectall(), "SpnInt.vtu", 2);
    Spn_external.write(selectall(), "SpnExt.vtu", 2);

    expression(tau).write(selectall(), "tau.vtu", 2);

    log_file.close();

    std::cout << "\n\n------------------------------------" << std::endl;
    std::cout << "c_st = " << expression(c_st).evaluate() << " m/s" << std::endl;
    myclock.print("\nTime elapsed for completion of simulation is ");

    return 0;

} // main