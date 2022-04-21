/**
   # Field scale wave breaking (multilayer solver)

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

/**
   The initial condition is a externally imported wave field. Some controlling parameters. */

#define g_ 9.8
double h_ = 10; // depth of the water
double gpe_base = 0; // gauge of potential energy
/* double ETAE = 0.1; // refinement criteria for eta */
/* int MAXLEVEL = 8; // max level of refinement in adapt_wavelet function */
/* int MINLEVEL = 6; // min level */
double TEND = 50.; // t end
int NLAYER = 10; // number of layers
int LEVEL_data = 7;

/* double kp_ = 2.*pi/10.; */
/* double P_ = 0.01; */
#define N_mode_ 32 // Corresponds to input of 32 modes in kx and 33 modes in ky
double F_kxky_[N_mode_*(N_mode_+1)], omega[N_mode_*(N_mode_+1)], phase[N_mode_*(N_mode_+1)];
double kx_[N_mode_], ky_[N_mode_+1];
double dkx_, dky_;

int RANDOM; // integer to seed random number generator

double randInRange(int min, int max)
{
  return min + (rand() / (double) (RAND_MAX) * (max - min + 1));
}

void power_input() {
  /* for (int i=0; i<N_mode_; i++) { */
  /*   kx_[i] = 2.*pi/L0*(i+1); */
  /*   ky_[i] = 2.*pi/L0*(i-N_mode_/2); */
  /* } */
  /* ky_[N_mode_] = 2.*pi/L0*N_mode_/2; */
/** 
A MPI compatible function that reads in kx_, ky_, and F_kxky_. Next step is to generate F_kxky_ inside basilisk too. For now kx_, ky_ are 1D arrays while F_kxky_ is a 2D array. The root process reads in and then broadcast to all other processes. It looks for files named F_kxky, 
*/
#if _MPI 
  int length1D, length2D; // The length of the array to be read
  char message[20];
  int i, rank, size;
  MPI_Status status;
  int root = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == root) {
    // First read in F_kxky_
    length2D = N_mode_*(N_mode_+1);
    float * a = (float*) malloc (sizeof(float)*length2D);
    char filename[100];
    sprintf (filename, "F_kxky");
    FILE * fp = fopen (filename, "rb");
    fread (a, sizeof(float), length2D, fp);
    for (int i=0;i<length2D;i++) {
      F_kxky_[i] = (double)a[i];
    }
    fclose (fp);

    // Then read in kx_, ky_
    length1D = N_mode_;
    float * b1 = (float*) malloc (sizeof(float)*length1D);
    sprintf (filename, "kx");
    FILE *fp1 = fopen (filename, "rb");
    fread (b1, sizeof(float), length1D, fp1);
    for (int i=0;i<length1D;i++) {
      kx_[i] = (double)b1[i];
    }
    fclose (fp1);

    // One more mode in ky
    float * b2 = (float*) malloc (sizeof(float)*(length1D+1));
    sprintf (filename, "ky");
    FILE *fp2 = fopen (filename, "rb");
    fread (b2, sizeof(float), length1D+1, fp2);
    for (int i=0;i<length1D+1;i++) {
      ky_[i] = (double)b2[i];
    }
    fclose (fp2);

    // Wave frequency omega, and randomly generated phase
    double kmod = 0;
    int index = 0;
    srand(RANDOM); // We can seed it differently for different runs
    for (int i=0; i<N_mode_; i++) {
      for (int j=0; j<N_mode_+1; j++) {
	index = j*N_mode_ + i;
	kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
	omega[index] = sqrt(g_*kmod);
	phase[index] = randInRange (0, 2.*pi);
      }
    }
  }
  MPI_Bcast(&kx_, length1D, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&ky_, length1D+1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&F_kxky_, length2D, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&omega, length2D, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&phase, length2D, MPI_DOUBLE, root, MPI_COMM_WORLD);
  // Make sure that the inputs are correct 
  char checkout[100];
  sprintf (checkout, "F-%d", pid());
  FILE * fout = fopen (checkout, "w");
  for (int i=0; i<length2D; i++)
    fprintf (fout, "%g ", F_kxky_[i]);
  fclose (fout);
  sprintf (checkout, "ky-%d", pid());
  fout = fopen (checkout, "w");
  for (int i=0; i<length1D+1; i++)
    fprintf (fout, "%g ", ky_[i]);
  fclose (fout);
#else
  int length2D = N_mode_*(N_mode_+1);
  float * a = (float*) malloc (sizeof(float)*length2D);
  char filename[100];
  sprintf (filename, "F_kxky");
  FILE * fp = fopen (filename, "rb");
  fread (a, sizeof(float), length2D, fp);
  for (int i=0;i<length2D;i++) {
    F_kxky_[i] = (double)a[i];
  }
  fclose (fp);
  // Phase and omega, next focusing phase
  double kmod = 0;
  int index = 0;
  srand(0); 
  for (int i=0; i<N_mode_;i++) {
    for (int j=0;j<N_mode_+1;j++) {
      index = j*N_mode_ + i;
      kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      omega[index] = sqrt(g_*kmod);
      phase[index] = randInRange (0, 2.*pi);
    }
  }
#endif
}

double wave (double x, double y)
{
  double eta = 0;
  double ampl = 0, a = 0;
  int index = 0;
  for (int i=0; i<N_mode_; i++) {
    for (int j=0; j<N_mode_+1; j++) {
      index = j*N_mode_ + i;
      ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      a = (kx_[i]*x + ky_[j]*y + phase[index]);
      eta += ampl*cos(a);
    }
  }
  return eta;
}
double u_x (double x, double y, double z) {
  int index = 0;
  double u_x = 0;
  double ampl = 0, a = 0;
  double z_actual = 0, kmod = 0, theta = 0;
  for (int i=0; i<N_mode_; i++) {
    for (int j=0; j<N_mode_+1; j++) {
      index = j*N_mode_ + i;
      ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      z_actual = (z < ampl ? (z) : ampl);
      // fprintf(stderr, "z = %g, ampl = %g, z_actual = %g\n", z, ampl, z_actual);
      kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      theta = atan(ky_[j]/kx_[i]);
      a = (kx_[i]*x + ky_[j]*y + phase[index]);
      u_x += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*cos(a)*cos(theta);
    }
  }
  return u_x;
}
// #if dimension = 2
double u_y (double x, double y, double z) {
  int index = 0;
  double u_y = 0;
  double ampl = 0, a = 0;
  double z_actual = 0, kmod = 0, theta = 0;
  for (int i=0; i<N_mode_; i++) {
    for (int j=0; j<N_mode_+1; j++) {
      index = j*N_mode_ + i;
      ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      z_actual = (z < ampl ? (z) : ampl);
      kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      theta = atan(ky_[j]/kx_[i]);
      a = (kx_[i]*x + ky_[j]*y + phase[index]);
      u_y += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*cos(a)*sin(theta);
    }
  }
  return u_y;
}
// #endif
double u_z (double x, double y, double z) {
  int index = 0;
  double u_z = 0;
  double ampl = 0, a = 0;
  double z_actual = 0, kmod = 0;
  for (int i=0; i<N_mode_; i++) {
    for (int j=0; j<N_mode_+1; j++) {
      index = j*N_mode_ + i;
      ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      z_actual = (z < ampl ? (z) : ampl);
      kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      a = (kx_[i]*x + ky_[j]*y + phase[index]);
      u_z += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*sin(a);
    }
  }
  return u_z;
}

/** 
## Parameters
*/

int main(int argc, char * argv[])
{
  if (argc > 1)
    NLAYER = atoi(argv[1]); // # of layers
  if (argc > 2)
    LEVEL_data = atoi(argv[2]); // Horizontal resolution
  if (argc > 3)
    TEND = atof(argv[3]); // End time
  if (argc > 4)
    nu = atof(argv[4]); // Viscosity for vertical diffusion
  else
    nu = 0.;
  if (argc > 5) 
    RANDOM = atoi(argv[5]); // An integer to seed the random number generator (don't use 0 or 1)
  if (argc > 6)
    L0 = atof(argv[6]); // Box size (not necessarily related to peak wave number)
  else
    L0 = 50.;
  if (argc > 7)
    theta_H = atof(argv[7]); // Numerical parameter to dump fast barotropic modes
  else
    theta_H = 0.5;
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  N = 1 << LEVEL_data; 
  nl = NLAYER;
  G = g_;
  h_ = L0/5.; // set the water depth to be 1/5 the field size (ad hoc)
  /** Use the already written remapping function. 
      coeff is the one defined with remap_test.h but is now replaced by the geometric beta function. 
      See example/breaking.c for details. 
      theta_H also follows the one defined in the breaking example. */
  // coeff = 0.05;
  // geometric_beta (1./3., true);
  // theta_H = 0.51;

#if dimension == 2
  gpe_base = -0.5*sq(h_)*sq(L0)*g_;
#else
  gpe_base = -0.5*sq(h_)*L0*g_;
#endif
  CFL_H = 1; // Smaller time step
  // max_slope = 0.8;
  run();
}

/**
## Initialization
Read in the power spectrum and initialize the wave field. 
 */
event init (i = 0)
{
  if (!restore ("restart")) {
    power_input();
    dkx_ = kx_[1] - kx_[0];
    dky_ = ky_[1] - ky_[0];
    fprintf (stderr, "dkx = %g, dky = %g\n", dkx_, dky_);
    geometric_beta (1./3., true); // Varying layer thickness
    foreach() {
      zb[] = -h_;
      eta[] = wave(x, y);
      double H = wave(x, y) - zb[];
      foreach_layer() {
	h[] = H/nl;
      }
    }
    // remap?
    vertical_remapping (h, tracers);
    foreach() {
      double z = zb[];
      foreach_layer() {
	z += h[]/2.;
	u.x[] = u_x(x, y, z);
	u.y[] = u_y(x, y, z);
	w[] = u_z(x, y, z);
	z += h[]/2.;
      }
    }
    fprintf (stderr,"Done initialization!\n");
    dump("initial");
  }
}

event energy_before_remap (i++, last)
{
  if (i==10) {
    fprintf(stderr, "energy output before remap!\n");
    fflush(stderr);
  }
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double zc = zb[];
    foreach_layer () {
      double norm2 = sq(w[]);
      foreach_dimension()
	norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
    }
  }
  static FILE * fp = fopen("energy_before_remap.dat","w");
  fprintf (fp, "%g %g %g\n", t, ke/2., g_*gpe - gpe_base);
  fflush (fp);
}

event energy_after_remap (i++, last)
{
  if (i==10) {
    fprintf(stderr, "energy output after remap!\n");
    fflush(stderr);
  }
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double zc = zb[];
    foreach_layer () {
      double norm2 = sq(w[]);
      foreach_dimension()
	norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
    }
  }
  static FILE * fp = fopen("energy_after_remap.dat","w");
  fprintf (fp, "%g %g %g\n", t, ke/2., g_*gpe - gpe_base);
  fflush (fp);
}

/**
   Note that the movie generation below is very expensive. */
#  define POPEN(name, mode) fopen (name ".ppm", mode)
#if 1
event movie (t += 0.1; t <= TEND)
{
  char s[80];
  view (fov = 20, quat = {0.475152,0.161235,0.235565,0.832313}, width = 800, height = 600);
  sprintf (s, "t = %.2f", t);
  draw_string (s, size = 30);
  sprintf (s, "u%d.x", nl-1);
  squares (s, linear = true, z = "eta", min = -2./7.*sqrt(L0), max = 2./7.*sqrt(L0));
  {
  static FILE * fp = POPEN ("ux", "a");
  save (fp = fp);
  }
  scalar slope[];
  foreach () {
    slope[] = (eta[1]-eta[-1])/(2.*Delta);
  }
  clear();
  squares ("slope", linear = true, z = "eta", min = -1./50.*L0, max = 1./50.*L0);
  sprintf (s, "t = %.2f", t);
  draw_string (s, size = 30);
  {
  static FILE * fp = POPEN ("slope", "a");
  save (fp = fp);
  }
  char filename1[50], filename2[50], filename3[50];
  sprintf (filename1, "surface/eta_matrix_%g", t);
  sprintf (filename2, "surface/ux_matrix_%g", t);
  sprintf (filename3, "surface/uy_matrix_%g", t);  
  FILE * feta = fopen (filename1, "w");
  // Might need to change to mpi function later
  output_matrix_mpi (eta, feta, N, linear = true);
  fclose (feta);
  sprintf (s, "u%d", nl-1);
  vector u_temp = lookup_vector (s);
  FILE * fux = fopen (filename2, "w");
  output_matrix_mpi (u_temp.x, fux, N, linear = true);
  fclose (fux);
  FILE * fuy = fopen (filename3, "w");
  output_matrix_mpi (u_temp.y, fuy, N, linear = true);
  fclose (fuy);  
}

#endif

#if QUADTREE
event adapt (i++) {
  /* fprintf(stderr, "Adapting start!\n"); */
  /* fflush(stderr); */
  my_adapt();
}
#endif


event endrun (t = TEND) {
  dump();
}
