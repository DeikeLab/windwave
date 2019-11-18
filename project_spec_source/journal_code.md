## Major change to the source code should be logged

-----------
Sometime around 11.10.19 or before
1. Changed the definition of c in the definition of ustar so that now it is `Ustar = sqrt(g_/k_+f.sigma*k_)*UstarRATIO`;
2. However, the definition of Re is not changed (still 1/\nu_w). It relies on outside conversion. Also the output time nondimensionalization is also not changed.

-----------
11.13.19
1. Add new files where there is a initial profile in water. 

-----------
11.18.19
1. Change the minmal level of refinement to see if that changes any result.
