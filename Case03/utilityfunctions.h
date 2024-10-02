#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include <iostream>
#include "sparselizard.h"

using namespace sl;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
// UTILITY FUNCTIONS

// spitzer heat conductivity
// K = k0e T^(5/2) [W/m/eV]

double spitzerheatconductivitycoefficient(void)
{
    /* Calculates the value of k0e */

    const double    c = 3e8;                     // [m/s]       speed of light
    const double   pi = 3.141592653589793238463; // [-]         value of pi
    const double    e = 1.6022e-19;              // [C]         unit charge of electron
    const double eps0 = 55.26349406e6 * e*e;     // [C²/(eV.m)] vaccum permittivity
    const double  m_e = 0.511e6 / (c*c);         // [eV/c²]     mass of electron
    const double  lnA = 15.0;                    // [-]         coulomb factor or plasma parameter
        
    // electron collision time coefficient
    double coeff_taue = 6.0 * pow(2.0,0.5) * pow(pi,3.0/2.0) * pow(eps0,2) * sqrt(m_e) / pow(e, 4); // ~= 3.4384e+11
    
    // Spitzer heat conductivity coefficient
    double coeff_k0e = 3.2 * e * coeff_taue / m_e; //  ~= 31048.8

    // Spitzer heat conductivity coefficient
    double k0e = coeff_k0e / lnA; //  ~= 2069.92

    std::cout << "\nSPITZER HEAT CONDUCTIVITY COEFFICIENT:"     << std::endl;
    std::cout << "---------------------------------------"     << std::endl;
    std::cout << "k0e = " << k0e << std::endl;

    return k0e;
}

double H(double x)
{
    /* Returns a Heavside step function*/

    return (x >= 0) ? 1.0 : 0.0;
}

expression H(expression x)
{
    /* Returns a Heavside step function*/

    return expression(x, 1.0, 0.0);
}


#endif