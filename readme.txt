///////////////////////////////////////////////////// 
// This code is part of the diagfrob MAGMA library //
//                                                 //
// copyright (c) 2017 Jan Tuitman                  //
/////////////////////////////////////////////////////


This is the diagfrob (Diagonal Frobenius) Magma library v1.1 for computing p-adic cohomology and L-factor of a diagonal hypersurface in projective space over a finite field F_p with p prime.

The code uses the algorithm developed in section 4 of:

Tuitman and Pancratz, "Improvements to the deformation method for counting points on smooth projective hypersurfaces".

but now also includes an improvement lowering the complexity in p from quasilinear in p to quasilinear in p^(1/2). This improvement uses the algorithm developed in:

Bostan, Gaudry and Schost, "Linear recurrences with polynomial coefficients and application to integer factorization and Cartier-Manin operator",

which was implemented in Magma by Minzlaff and Fontein.

An example of how to use the code can be found in the file t_diagfrob.m.

Jan Tuitman,

February 2017. 
