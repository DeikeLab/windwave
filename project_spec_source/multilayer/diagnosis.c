/**
   Diagnosis of multilayer breaking wave field.
*/

#include "grid/multigrid.h"
#include "view.h"
#include "layered/hydro.h"
//#include "layered/hydro_test.h"
#include "layered/nh.h"
#include "layered/remap.h"
//#include "remap_test.h"
#include "layered/perfs.h"
#include "input.h"
#include "output_mpi.h"

#define g_ 9.8
double h_ = 10; // depth of the water
double gpe_base = 0; // gauge of potential energy
/* double ETAE = 0.1; // refinement criteria for eta */
/* int MAXLEVEL = 8; // max level of refinement in adapt_wavelet function */
/* int MINLEVEL = 6; // min level */
double TEND = 50.; // t end
int NLAYER = 10; // number of layers
int LEVEL_data = 7;


int main(int argc, char * argv[])
{
  if (argc > 1)
    NLAYER = atoi(argv[1]);
  if (argc > 2)
    LEVEL_data = atoi(argv[2]);
  if (argc > 3)
    L0 = atof(argv[3]);
  else
    L0 = 50.;
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  N = 1 << LEVEL_data; // start with a grid of 128
  nl = NLAYER;
  G = g_;
  h_ = L0/5.; // set the water depth to be 1/5 the field size (ad hoc)
  run();
}

event init (i = 0)
{
  if (!restore ("dump")) {
    fprintf (ferr, "Not restored!\n");
    return 1;
  }
  char s[80];
  char filename1[50], filename2[50], filename3[50], filename4[50];
  for (int j=0; j<nl; ++j) {
    sprintf (filename1, "field/ux_matrix_l%d", j);
    sprintf (filename2, "field/uy_matrix_l%d", j);  
    sprintf (filename3, "field/uz_matrix_l%d", j);  
    sprintf (filename4, "field/h_matrix_l%d", j);  
    sprintf (s, "u%d", j);
    vector u_temp = lookup_vector (s);
    FILE * fux = fopen (filename1, "w");
    output_matrix_mpi (u_temp.x, fux, N, linear = true);
    fclose (fux);
    FILE * fuy = fopen (filename2, "w");
    output_matrix_mpi (u_temp.y, fuy, N, linear = true);
    fclose (fuy);
    FILE * fuz = fopen (filename3, "w");
    sprintf (s, "w%d", j);
    scalar w_temp = lookup_field (s);
    output_matrix_mpi (w_temp, fuz, N, linear = true);
    fclose (fuz);
    FILE * fh = fopen (filename4, "w");
    sprintf (s, "h%d", j);
    scalar h_temp = lookup_field (s);
    output_matrix_mpi (h_temp, fh, N, linear = true);
    fclose (fh);    
  }
}

