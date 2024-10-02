From SD1D:

Case 2 : Fluid model 
====================
*Like case-01 but with a source only in the first half of the domain*

Constant and uniform source functions, only in the first half of the domain.
No heat conduction
No neutrals --> Mn, Nn = 0. Tn=Te for Diffusion coefficient calculation and 0 otherwise.


**My notes:**
The source is constant and uniform but only till the first half of the domain, i.e until the X-point.

Projection of source:
The results seem good with projection of NeSource and EeSource.
Results remained the same without projection as well. 

with heat conduction, with AMR and with Projection (SUPG was required) - RESULTS ARE GOOD

with only heat convection: ramping down really helps
(For presentation: progressive graph with beta starting from 1 and reducing it to as low as possible)

To try:
inequality boundary condition for Mach number at the target.

The uniform difference in temperature profiles between SD1D and Sparselizard could be due to Energy source and sink. 
Recheck the corresponding rate coefficients.
