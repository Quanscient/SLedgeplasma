#include "sparselizard.h"

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
    #include "preprocess_1D_SOL.h"

    #include "conservationmassmomentum1D.h"
    #include "conservationneutral1D.h"
    #include "conservationenergy1D.h"

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /* NEWTON SOLVER CONFIGURATION: CS=coupled_solver
    --------------------------------------------------*/

    // rampup/rampdown flags
    bool rampedup=false, M_rampedup=false, N_rampedup=false;
    bool q_rampeddown=false;

    // relaxation factor
    double w_Tt  = 1.00;
    double w_CS  = 1.00;

    // Max iterations
    int maxiters_Tt = 100;          // max iterations considering all q_sol
    int maxiters_CS = 300;          // max iterations for coupled solver for a specified Tt

    // Convergence criteria
    double tolrelerr_qt = 1e-6;
    double tolrelerr_Tt = 1e-6;
    double tolrelerr_CS = 1e-10;

    // Variables used in convergence check
    double relerr_qt = 1.0;
    double relerr_Tt = 1.0;
    double relerr_CS = 1.0;
    bool isconverged_Tt = false;
    bool isconverged_CS = false;

    // Solver iterations
    int iter_total = 1;
    int iter_CS    = 0;
    int iter_Tt    = 0;

    // Variables for upstream/target field values
    double  N_upstream=0.0,  N_target=0.0;
    double  M_upstream=0.0,  M_target=0.0;
    double Nn_upstream=0.0, Nn_target=0.0;
    double  T_upstream=0.0,  T_target=0.0;

    // Flux balance variables
    double massflux_upstream = 0.0, momnetumflux_upstream = 0.0, neutralflux_upstream = 0.0, energyflux_upstream = 0.0;
    double massflux_target   = 0.0, momnetumflux_target   = 0.0, neutralflux_target   = 0.0, energyflux_target   = 0.0;
    double massflux_balance  = 0.0, momnetumflux_balance  = 0.0, neutralflux_balance  = 0.0, energyflux_balance  = 0.0;

    // log file for convergence evolution data
    ofstream log_file;
    log_file.open("convergence_log.csv", ios::out | ios::trunc);
    if (log_file.is_open()){
        log_file << "#iter_Tt"       << "," \
                 << "q_sol"          << "," \
                 << "CS_converged@"  << "," \
                 << "relerr_CS"      << "," \
                 << "N_upstream"     << "," \
                 << "M_upstream"     << "," \
                 << "Nn_upstream"    << "," \
                 << "T_upstream"     << "," \
                 << "N_target"       << "," \
                 << "M_target"       << "," \
                 << "Nn_target"      << "," \
                 << "T_target"       << "," \
                 << "relerr_Tt"      << "," \
                 << "relerr_qt"      << "," \
                 << "isconverged_Tt" << std::endl;
    }
    else{
        std::cout << "Unable to open the log file" << std::endl;
    }

    // solution vectors (used in coupled sovler loop)
    vec sol_prev(plasma_flow), sol_pred(plasma_flow), sol_curr(plasma_flow);

    // Target temperature variables (used in target temperature loop)
    double Tt_curr=0.0, Tt_pred=0.0, Tt_next=0.0;
    Tt_curr = Tti;

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /* NETWON SOVLER
    ----------------------*/

    // HEAT FLUX RAMPEDDOWN LOOP
    while (true)
    {
        // reset for each new q_sol
        iter_Tt = 1;
        isconverged_Tt = false;

        // Log output
        std::cout << std::endl << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "Current q_sol = " << q_next << std::endl;
        
        // TARGET TEMPERATURE LOOP (Solves for target temperature for the current q_sol)
        while( isconverged_Tt==false && iter_Tt<maxiters_Tt )
        {
            // reset for each new Tt
            iter_CS = 1;
            isconverged_CS = false;
            
            // Log output
            std::cout << "\n\n------------------------------------" << std::endl;
            std::cout << "iter_total:" << iter_total << "\niter_Tt: " << iter_Tt << std::endl;
            std::cout << "------------" << std::endl;

            /*-------------------------------------- COUPLED SOLVER (NONLINEAR) -------------------------------------- */
            // COUPLED SOLVER LOOP (Solves for <N, M, Nn, T> for the current target temperature)
            while( isconverged_CS==false && iter_CS<maxiters_CS ){
                plasma_flow.solve("lu", false);                         // solve the coupled plasma formulation                
                Nn.setvalue(selectall(), max(1e-30, Nn));               // floor Nn (because negative Nn are unphysical)
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
            }
            std::cout <<   "iter_CS: " <<   iter_CS << ";\t";
            std::cout << "relerr_CS: " << relerr_CS << ";\t";
            std::cout << "isconverged_CS: " << std::boolalpha << isconverged_CS << std::endl;

            /*------------------------------------------ POST-PROCESSING FOR CS------------------------------------------ */
            // Get the field values at upstream and target at each iteration
            M_upstream  =  M.integrate(upstream, 5);
            N_upstream  =  N.integrate(upstream, 5);
            Nn_upstream = Nn.integrate(upstream, 5);
            T_upstream  =  T.integrate(upstream, 5);
            M_target    =  M.integrate(target, 5);
            N_target    =  N.integrate(target, 5);
            Nn_target   = Nn.integrate(target, 5);
            T_target    =  T.integrate(target, 5);
    
            std::cout << "M  = [" <<  M_upstream << ", " <<  M_target << "]" << std::endl;
            std::cout << "N  = [" <<  N_upstream << ", " <<  N_target << "]" << std::endl;
            std::cout << "Nn = [" << Nn_upstream << ", " << Nn_target << "]" << std::endl;
            std::cout << "T  = [" <<  T_upstream << ", " <<  T_target << "]" << std::endl;

            // Append the convergence evolution data to the log file
            if (log_file.is_open()){
                log_file << std::scientific << std::setprecision(0) << iter_Tt     << "," \
                         << std::scientific << std::setprecision(1) << q_next         << "," \
                         << std::scientific << std::setprecision(5) << iter_CS     << "," \
                         << std::scientific << std::setprecision(6) << relerr_CS   << "," \
                         << std::scientific << std::setprecision(1) << N_upstream  << "," \
                         << std::scientific << std::setprecision(5) << M_upstream  << "," \
                         << std::scientific << std::setprecision(6) << Nn_upstream << "," \
                         << std::scientific << std::setprecision(6) << T_upstream  << "," \
                         << std::scientific << std::setprecision(5) << N_target    << "," \
                         << std::scientific << std::setprecision(1) << M_target    << "," \
                         << std::scientific << std::setprecision(5) << Nn_target   << "," \
                         << std::scientific << std::setprecision(5) << T_target    << "," \
                         << std::scientific << std::setprecision(1) << relerr_Tt   << "," \
                         << std::scientific << std::setprecision(5) << relerr_qt   << "," \
                         << std::boolalpha  << isconverged_Tt << std::endl;
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

            expression neutralflux = on(SOL, -D*grad(Nn));
            neutralflux_upstream = neutralflux.integrate(upstream, 5);
            neutralflux_target   = -R*massflux_target;  // imposed // neutralflux.integrate(target, 5);

            expression energyflux = on(SOL, -K*grad(T));
            energyflux_upstream = expression(q_u).evaluate(); // q_next // imposed //energyflux.integrate(upstream, 5);
            energyflux_target   = energyflux.integrate(target, 5);

            // Flux balances
            massflux_balance     = massflux_upstream     +     (Spn).integrate(SOL, 5);
            momnetumflux_balance = momnetumflux_upstream + (-M*Scxn).integrate(SOL, 5);
            neutralflux_balance  = neutralflux_upstream  +     (Snn).integrate(SOL, 5);
            energyflux_balance   = energyflux_upstream; // no volumetric energy source/sink considered

            std::cout << "massbalance     : " << massflux_target     << "\t\t" << massflux_balance     << std::endl;
            std::cout << "momentumbalance : " << momnetumflux_target << "\t\t" << momnetumflux_balance << std::endl;
            std::cout << "neutralbalance  : " << neutralflux_target  << "\t\t" << neutralflux_balance  << std::endl;
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
                    isconverged_Tt = false;
                    std::cout << "M & N ramping up ..." << std::endl;
                }
            }

            /*-------------------------------------- TARGET TEMPERATURE SOLVER -------------------------------------- */
            // Solve for new target temperature with the latest N and M (only after N and M)
            if(rampedup==true){
                expression gradT = on(SOL, grad(T));        // calculate temperature gradient
                double unit_eV = 1.6022e-19;                // unit charge: 1eV = 1.6022e-19 Joules
                double gamma_se = 7.0;                      // total sheath heat transmission coefficient, γ ≃ 7–8
                expression GAMMA_se = n_o*c_st * N*M;       // plasma flux at sheath edge

                // heat flux
                expression q_sim = -K * gradT;         // q = -K .∇T
                double qt_sim = q_sim.integrate(target, 5);
                double qt_exp = (verify2ptmodel) ? q_next : (gamma_se*T*GAMMA_se*unit_eV).integrate(target, 5); // expected
                relerr_qt = sl::abs( (qt_sim - qt_exp) / q_u).evaluate();

                std::cout << "------------" << std::endl;
                std::cout << "(qt_sim,    qt_exp)    = (" << qt_sim << ", " << qt_exp << ")" << std::endl;

                // Update target temperature if relerr_qt is above tolerance
                if(relerr_qt > tolrelerr_qt){
                    // calculate target temperature from --> q_se = γ eT (nv)_target
                    if(verify2ptmodel){
                        Tt_pred = (q_next/(gamma_se*unit_eV*GAMMA_se)).integrate(target, 5); // Onion skin model
                    }
                    else{
                        Tt_pred = (q_sim/(gamma_se*unit_eV*GAMMA_se)).integrate(target, 5);
                    }
                    
                    Tt_next = (q_next<1e-8) ? Tt_pred*w_Tt + Tt_curr*(1-w_Tt) : Tt_pred;
                    Tt.setvalue(SOL, Tt_next);              // update the target temperature parameter
                    c_st.setvalue(SOL, sqrt(2.0*Tt/m_i));   // update the c_st parameter
                    T.setconstraint(target, Tt);            // update the target temperature constraint
                                
                    relerr_Tt = sl::abs(Tt_next - Tt_curr).evaluate() / Tt_next;
                    isconverged_Tt = (relerr_Tt < tolrelerr_Tt) ? true : false;

                    /* relerr_qt is much better for energybalance than relerr_Tt, but also slower convergence. 
                        May not be necessary.*/
                    
                    std::cout << "(relerr_qt, relerr_Tt) = (" << relerr_qt << ", " << relerr_Tt << ")" << std::endl;
                    std::cout << "Tt = " << Tt_next;

                    // Update for next iteration of target temperature
                    Tt_curr = Tt_next;
                }
                else {
                    isconverged_Tt=true;
                }
                std::cout << "; isconverged_Tt: " << std::boolalpha << isconverged_Tt <<std::endl;
                std::cout << "------------" << std::endl;

                std::cout << "meshadaptivity: " << std::boolalpha << adapt() << std::endl; // important for maintaining flux balances.

                if(isconverged_Tt){
                    std::cout << "Simulation completed for q_sol=" << q_next << std::endl;
                }
            } // new target temperature       

            // move to the next iteration in Tt loop
            iter_Tt++;
            iter_total++;

            myclock.print("Total time elapsed so far: ");

        } // target temperature loop

        /*------------------------------------------ POST-PROCESSING FOR Tt------------------------------------------ */
        ostringstream qsol_ss;
        qsol_ss << std::scientific << std::setprecision(2)  << q_next;

         T.write(SOL,  "T_q" + qsol_ss.str() + ".vtu", 1);
         M.write(SOL,  "M_q" + qsol_ss.str() + ".vtu", 1);
         N.write(SOL,  "N_q" + qsol_ss.str() + ".vtu", 1);
        Nn.write(SOL, "Nn_q" + qsol_ss.str() + ".vtu", 1);

        /*------------------------------------------ UPDATE FOR NEXT q_sol Loop ------------------------------------------ */
        if (q_rampeddown==false){
            q_next = expression(q_u).evaluate() * rampdown_ratio;
            if(q_next > q_end){
                q_u.setvalue(SOL, q_next);
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
    Nn.write(SOL,"Nn.vtu", 1);

    log_file.close();

    myclock.print("\nTime elapsed for completion of simulation is ");

    return 0;

} // main