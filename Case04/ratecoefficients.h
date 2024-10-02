#ifndef RATE_COEFFICIENTS_H
#define RATE_COEFFICIENTS_H

/*
| recombination     | Rrc | ion-electron     | n*n*sigVrc  |
| ionization        | Riz | electron-neutral | n*nn*sigViz |
| chargeexchange    | Rcx | ion-neutral      | n*nn*sigVcx |
| elasticscattering | Rel | ion-neutral      | n*nn*sigVel |

| excitation        | Rex | electron-neutral | n*nn*sigVex |
| impurityradiation | Rrd |                  |             |


Sp = Sp_ext + Sp_int
Sn = Sn_ext + Sn_int
Sm = Sm_ext + Sn_int
Se = Se_ext + Se_int

                                        in SD1D
Sp_int = -Src + Siz                 -->   -S
Sn_int =  Src - Siz                 -->    S
Sm_int = -Frc + Fiz - Fcx - Fel     -->   -F
Se_int = -Erc + Eiz - Ecx - Eel     -->   -E
         -Crc - Ciz - Cex - Crd     -->   -R
*/




/////////////////////////////////////////////////////////////////////////////
// HYDROGEN RATE COEFFICIENTS <σv> [m³/s] --> sigV
// (Taken from SD1D source code)

// 1. RECOMBINATION, sigVrc
expression sigVrc, Trc;

if (!updatedradiatedpower)
    Trc = expression(T-1.0, T, 1.0);
else
    Trc = expression(T-0.025, T, 0.025);

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

expression ln_sigVrc = 0.0;
for (int n=0; n<=8; n++){
    for (int m=0; m<=8; m++){
        ln_sigVrc = ln_sigVrc + (reccoeffs[n][m] * pow(log(E),m) * pow(log(Trc),n));   // eirene.de/amjuel.pdf: Example of Use of Fits
    }
}
// convert [cm³/s] * 1e-6 --> [m³/s]
sigVrc = ifpositive(N*n_o-1e3, exp(ln_sigVrc) * 1e-6, 0.0);

if (updatedradiatedpower){
    double A = 3.92E-20;
    double B = 3.0E-124 * pow(1.6E-19, -4.5);
    double Ry = 13.60569;
    double chi = 0.35;
    expression sigRAD = A * pow(Ry, 1.5) * 1.0 / (sqrt(T) * (Ry + chi * T));

    sigVrc = sigVrc + sigRAD + B*(n_o*N)*pow(T, -5.0);
}

/////////////////////////////////////////////////////////////////////////////
// 2. IONIZATION, sigViz
expression sigViz, Tiz;

if (!updatedradiatedpower)
{

    Tiz = expression(T-1.0, T, 1.0);
    expression X = log(Tiz)/log(10); 

    expression S, S1, S2;
    S1 = -0.5151*X -2.563/X - 5.231 - 6.0;
    S2 = -3.054*X -15.72*exp(-X) + 1.603*exp(-X*X) - 6.0;
    S = expression(Tiz-20, S1, S2);
  
    sigViz = pow(10.0, S);   

}
else
{
    expression Tiz(T-0.025, T, 0.025);

    double ioncoeffs[9] = { -3.271397E1,   1.353656E1,  -5.739329,    \
                             1.563155,    -2.877056E-1,  3.482560e-2, \
                            -2.631976E-3,  1.119544E-4, -2.039150E-6};

    expression ln_sigViz = 0.0;
    for (int n=0; n<=8; n++)
    {
        ln_sigViz = ln_sigViz + ioncoeffs[n] * pow(log(Tiz),n);
    }
    
    // convert [cm³/s] * 1e-6 --> [m³/s]
    sigViz = exp(ln_sigViz) * 1e-6;

}

/////////////////////////////////////////////////////////////////////////////
// 3. CHARGE EXCHANGE, sigVcx
expression sigVcx, Tcx; 

if (!updatedradiatedpower)
{
    Tcx = expression(T-1.0, T, 1.0);    // if (T<1eV) --> T = 1eV;
    expression S = -14.0 + (log(Tcx)/log(10))/3.0;
    sigVcx = pow(10.0,S); 
}
else{
    Tcx = expression(T-0.025, T, 0.025);    // if (T<0.025eV) --> T = 0.025eV;
    expression E = 10.0;

    double cxcoeffs[9][9] = {
        {-1.829079582E+01,  1.640252721E-01,  3.364564509E-02,  9.530225559E-03, -8.519413900E-04, -1.247583861E-03,  3.014307546E-04, -2.499323170E-05,  6.932627238E-07,},
        { 2.169137616E-01, -1.106722014E-01, -1.382158680E-03,  7.348786287E-03, -6.343059502E-04, -1.919569450E-04,  4.075019352E-05, -2.850044983E-06,  6.966822400E-08,},
        { 4.307131244E-02,  8.948693625E-03, -1.209480567E-02, -3.675019470E-04,  1.039643391E-03, -1.553840718E-04,  2.670827249E-06,  7.695300598E-07, -3.783302282E-08,},
        {-5.754895093E-04,  6.062141761E-03,  1.075907882E-03, -8.119301728E-04,  8.911036876E-06,  3.175388950E-05, -4.515123642E-06,  2.187439284E-07, -2.911233952E-09,},
        {-1.552077120E-03, -1.210431588E-03,  8.297212634E-04,  1.361661817E-04, -1.008928628E-04,  1.080693990E-05,  5.106059414E-07, -1.299275586E-07,  5.117133050E-09,},
        {-1.876800283E-04, -4.052878752E-05, -1.907025663E-04,  1.141663042E-05,  1.775681984E-05, -3.149286924E-06,  3.105491555E-08,  2.274394089E-08, -1.130988251E-09,},
        { 1.125490271E-04,  2.875900436E-05,  1.338839629E-05, -4.340802793E-06, -7.003521917E-07,  2.318308730E-07, -6.030983538E-09, -1.755944926E-09,  1.005189187E-10,},
        {-1.238982763E-05, -2.616998140E-06, -1.171762874E-07,  3.517971869E-07, -4.928692833E-08,  1.756388999E-10, -1.446756796E-10,  7.143183138E-11, -3.989884106E-12,},
        { 4.163596197E-07,  7.558092849E-08, -1.328404104E-08, -9.170850254E-09,  3.208853884E-09, -3.952740759E-10,  2.739558476E-11, -1.693040209E-12,  6.388219930E-14,},
    };

    expression ln_sigVcx = 0.0;
    for (int n=0; n<=8; n++){
        for (int m=0; m<=8; m++){
            ln_sigVcx = ln_sigVcx + (reccoeffs[n][m] * pow(log(E),m) * pow(log(Tcx),n));   // eirene.de/amjuel.pdf: Example of Use of Fits
        }
    }
    // convert [cm³/s] * 1e-6 --> [m³/s]
    sigVcx = exp(ln_sigVcx) * 1e-6;

}


/////////////////////////////////////////////////////////////////////////////
// 4. EXCITATION, sigVex
expression sigVex, Tex;

Tex = expression(T-1.0, T, 1.0);    // if (T<1eV) --> T = 1eV;
expression Yex = 10.2/Tex;
sigVex = (49.0e-14/(0.28+Yex)) * exp(-Yex) * sqrt(Yex*(1.0+Yex));


/////////////////////////////////////////////////////////////////////////////
// HYDROGEN RATE OF EVENTS
// expression Rrcn = n_o*n_o*N*N*sigVrc  / NVo;
// expression Rizn = n_o*n_o*N*Nn*sigViz / NVo;
// expression Rexn = n_o*n_o*N*N*sigVex  / NVo;

#endif
