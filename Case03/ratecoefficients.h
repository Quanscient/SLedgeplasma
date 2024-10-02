#include "sparselizard.h"
#include <iostream>

/* RATE COEFFICIENTS <σv> [m³/s] --> sigV*/

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// HYDROGEN RADIATED POWER

// // 1. recombination rate coefficient, sigVrec
// expression Tcr(T-1.0, T, 1.0);    // if (T<1eV) --> T = 1eV;
// expression E = n_o*N * 1e-14;
// double reccoeffs[9][9] = {
//     {-2.855728479302E+01, -7.664042607917E-01, -4.930424003280E-03, -5.386830982777E-03, -1.626039237665E-04,  6.080907650243E-06,  2.101102051942E-05, -2.770717597683E-06,  1.038235939800E-07,},
//     { 3.488563234375E-02, -3.583233366133E-03, -3.620245352252E-03, -9.532840484460E-04,  1.888048628708E-04, -1.014890683861E-05,  2.245676563601E-05, -4.695982369246E-06,  2.523166611507E-07,},
//     {-2.799644392058E-02, -7.452514292790E-03,  6.958711963182E-03,  4.631753807534E-04,  1.288577690147E-04, -1.145028889459E-04, -2.245624273814E-06,  3.250878872873E-06, -2.145390398476E-07,},
//     { 1.209545317879E-02,  2.709299760454E-03, -2.139257298118E-03, -5.371179699661E-04, -1.634580516353E-05,  5.942193980802E-05, -2.944873763540E-06, -9.387290785993E-07,  7.381435237585E-08,},
//     {-2.436630799820E-03, -7.745129766167E-04,  4.603883706734E-04,  1.543350502150E-04, -9.601036952725E-06, -1.211851723717E-05,  1.002105099354E-06,  1.392391630459E-07, -1.299713684966E-08,},
//     { 2.837893719800E-04,  1.142444698207E-04, -5.991636837395E-05, -2.257565836876E-05,  3.425262385387E-06,  1.118965496365E-06, -1.291320799814E-07, -1.139093288575E-08,  1.265189576423E-09,},
//     {-1.886511169084E-05, -9.382783518064E-06,  4.729262545726E-06,  1.730782954588E-06, -4.077019941998E-07, -4.275321573501E-08,  7.786155463269E-09,  5.178505597480E-10, -6.854203970018E-11,},
//     { 6.752155602894E-07,  3.902800099653E-07, -1.993485395689E-07, -6.618240780594E-08,  2.042041097083E-08,  3.708616111085E-10, -2.441127783437E-10, -9.452402157390E-12,  1.836615031798E-12,},
//     {-1.005893858779E-08, -6.387411585521E-09,  3.352589865190E-09,  1.013364275013E-09, -3.707977721109E-10,  7.068450112690E-12,  3.773208484020E-12, -4.672724022059E-14, -1.640492364811E-14,},
// };
// expression ln_sigVrec = 0.0;
// for (int n=0; n<=8; n++){
//     for (int m=0; m<=8; m++){
//         ln_sigVrec = ln_sigVrec + (reccoeffs[n][m] * pow(log(E),m) * pow(log(Tcr),n));   // eirene.de/amjuel.pdf: Example of Use of Fits
//     }
// }
// // [cm³/s] * 1e-6 --> [m³/s]
// expression sigVrec = ifpositive(N*n_o-1e3, exp(ln_sigVrec) * 1e-6, 0.0);

// // 2. excitation rate coefficient, sigVex
// expression Yex = 10.2/Tcr;
// expression sigVex = (49.0e-14/(0.28+Yex)) * exp(-Yex) * sqrt(Yex*(1.0+Yex));


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// UPDATED RADIATED POWER

// 1. recombination rate coefficient, sigVrec
expression Tcr(T-0.025, T, 0.025);    // if (T<1eV) --> T = 1eV;
expression E = n_o*N * 1e-14;
double reccoeffs[9][9] = {
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
expression ln_sigHAV = 0.0;
for (int n=0; n<=8; n++){
    for (int m=0; m<=8; m++){
        ln_sigHAV = ln_sigHAV + (reccoeffs[n][m] * pow(log(E),m) * pow(log(Tcr),n));   // eirene.de/amjuel.pdf: Example of Use of Fits
    }
}
// [cm³/s] * 1e-6 --> [m³/s]
expression sigHAV = ifpositive(N*n_o-1e3, exp(ln_sigHAV) * 1e-6, 0.0);

double A = 3.92E-20;
double B = 3.0E-124 * pow(1.6E-19, -4.5);
double Ry = 13.60569;
double chi = 0.35;
expression sigRAD = A * pow(Ry, 1.5) * 1.0 / (sqrt(T) * (Ry + chi * T));

expression sigVrec = sigHAV + sigRAD + B*(n_o*N)*pow(T, -5.0);

// 2. excitation rate coefficient, sigVex
expression Yex = 10.2/Tcr;
expression sigVex = (49.0e-14/(0.28+Yex)) * exp(-Yex) * sqrt(Yex*(1.0+Yex));
