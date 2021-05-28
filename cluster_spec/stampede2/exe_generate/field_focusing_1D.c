/**
   # Field scale wave breaking (multilayer solver)

*/

#include "grid/multigrid.h"
#include "view.h"
#include "layered/hydro.h"
#include "layered/nh.h"
//#include "layered/remap.h"
#include "remap_test.h"
#include "layered/perfs.h"
#include "input.h"
#include "output_mpi.h"

/**
   The initial condition is a externally imported wave field. Some controlling parameters. */

#define g_ 9.8
#define h_ 10
double gpe_base = 0;
/* double ETAE = 0.1; // refinement criteria for eta */
/* int MAXLEVEL = 8; // max level of refinement in adapt_wavelet function */
/* int MINLEVEL = 6; // min level */
double TEND = 50.; // t end
int NLAYER = 10; // number of layers
int LEVEL_data = 7;

double kp_ = 2.*pi/10.;
double P_ = 0.01;
int N_power_ = 5;
#define N_mode_ 16
double F_kx_[N_mode_], omega[N_mode_], phase[N_mode_]; // One more node for ky_
double kx_[N_mode_];
double dkx_;
// Focusing time and place
double xb_ = 0, tb_ = 40;

void power_input () {
  for (int i=0; i<N_mode_; i++) {
    kx_[i] = 2.*pi/L0*(i+1);
  }
  dkx_ = kx_[1]-kx_[0];
  // A function that reads in F_kxky_. Next step is to generate F_kxky_ all inside
#if _MPI
  int length = N_mode_;
  char message[20];
  int  i, rank, size;
  MPI_Status status;
  int root = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == root) {
    int length = N_mode_;
    float * a = (float*) malloc (sizeof(float)*length);
    char filename[100];
    sprintf (filename, "F_kx");
    FILE * fp = fopen (filename, "rb");
    fread (a, sizeof(float), length, fp);
    for (int i=0;i<length;i++) {
      F_kx_[i] = (double)a[i];
      //fprintf(stderr, "%g \n", F_kxky_[i]);
    }
    fclose (fp);
    // Phase and omega, next focusing phase
    double kmod = 0;
    int index = 0;
    for (int i=0; i<N_mode_;i++) {
    	index = i;
    	kmod = kx_[i];
    	omega[index] = sqrt(g_*kmod);
    	phase[index] = - kx_[i]*xb_ + omega[index]*tb_;
    }
    /* for (int j=0; j<(N_mode_+1); j++) { */
    /*   for (int i=0; i<N_mode_; i++) { */
    /* 	index = i*N_mode_ + j; */
    /* 	kmod = sqrt(sq(kx_[i]) + sq(ky_[j])); */
    /* 	omega[index] = sqrt(g_*kmod); */
    /* 	phase[index] = - kx_[i]*xb_ - ky_[j]*yb_ + omega[index]*tb_; */
    /*   } */
    /* } */
  }
  MPI_Bcast(&F_kx_, length, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&omega, length, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&phase, length, MPI_DOUBLE, root, MPI_COMM_WORLD);
  //printf( "Message from process %d : %s\n", rank, message);
  char checkout[100];
  sprintf (checkout, "F-%d", pid());
  FILE * fout = fopen (checkout, "w");
  for (int i=0; i<length; i++)
    fprintf (fout, "%g ", F_kx_[i]);
  fclose (fout);
  sprintf (checkout, "phase-%d", pid());
  fout = fopen (checkout, "w");
  for (int i=0; i<length; i++)
    fprintf (fout, "%g ", phase[i]);
  fclose (fout);
#else
  int length = N_mode_;
  float * a = (float*) malloc (sizeof(float)*length);
  char filename[100];
  sprintf (filename, "F_kx");
  FILE * fp = fopen (filename, "rb");
  fread (a, sizeof(float), length, fp);
  for (int i=0;i<length;i++) {
    F_kx_[i] = (double)a[i];
    //fprintf(stderr, "%g \n", F_kxky_[i]);
  }
  fclose (fp);
  char checkout[100];
  sprintf (checkout, "F");
  FILE * fout = fopen (checkout, "w");
  for (int i=0; i<length; i++)
    fprintf (fout, "%g ", F_kx_[i]);
  fclose (fout);
  // Phase and omega, next focusing phase
  double kmod = 0;
  int index = 0; 
  for (int i=0; i<N_mode_;i++) {
      index = i;
      kmod = kx_[i];
      omega[index] = sqrt(g_*kmod);
      phase[index] = -kx_[i]*xb_ + omega[index]*tb_;
  }
#endif
}

double wave (double x, double y)
{
  double eta = 0;
  double ampl = 0, a = 0;
  int index = 0;
  for (int i=0; i<N_mode_;i++) {
      index = i;
      ampl = sqrt(2.*F_kx_[index]*dkx_);
      a = kx_[i]*x + phase[index];
      eta += ampl*cos(a);
  }
  return eta;
}
double u_x (double x, double y, double z) {
  int index = 0;
  double u_x = 0;
  double ampl = 0, a = 0;
  double z_actual = 0, kmod = 0, theta = 0;
  for (int i=0; i<N_mode_;i++) {
      index = i;
      ampl = sqrt(2.*F_kx_[index]*dkx_);
      z_actual = (z < ampl ? (z) : ampl);
      //fprintf(stderr, "z = %g, ampl = %g, z_actual = %g\n", z, ampl, z_actual);
      kmod = kx_[i];
      a = kx_[i]*x + phase[index];
      u_x += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*cos(a)*cos(theta);
  }
  return u_x;
}
// #if dimension = 2
double u_y (double x, double y, double z) {
  double u_y = 0;
  return u_y;
}
// #endif
double u_z (double x, double y, double z) {
  int index = 0;
  double u_z = 0;
  double ampl = 0, a = 0;
  double z_actual = 0, kmod = 0;
  for (int i=0; i<N_mode_;i++) {
      index = i;
      ampl = sqrt(2.*F_kx_[index]*dkx_);
      z_actual = (z < ampl ? (z) : ampl);
      kmod = kx_[i];
      a = kx_[i]*x + phase[index];
      u_z += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*sin(a);
  }
  return u_z;
}


int main(int argc, char * argv[])
{
  if (argc > 1)
    NLAYER = atoi(argv[1]);
  if (argc > 2)
    LEVEL_data = atoi(argv[2]);
  if (argc > 3)
    TEND = atof(argv[3]);
  if (argc > 4)
    nu = atof(argv[4]);
  else
    nu = 0.;
  L0 = 50.;
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  N = 1 << LEVEL_data; // start with a grid of 128
  nl = NLAYER;
  G = g_;
  coeff = 0.05;
#if dimension == 2
  gpe_base = -0.5*sq(h_)*sq(L0)*g_;
#else
  gpe_base = -0.5*sq(h_)*L0*g_;
#endif
  CFL_H = 1; // Smaller time step
  //  max_slope = 1.;
  run();
}

event init (i = 0)
{
  if (!restore ("restart")) {
    power_input();
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
  sprintf (s, "u%d", nl-1);
  vector u_temp = lookup_vector (s);
  view (fov = 20, quat = {0.475152,0.161235,0.235565,0.832313}, width = 800, height = 600);
  sprintf (s, "t = %.2f", t);
  draw_string (s, size = 30);
  sprintf (s, "u%d.x", nl-1);
  squares (s, linear = true, z = "eta", min = -2, max = 2);
  {
  static FILE * fp = POPEN ("ux", "a");
  save (fp = fp);
  }
  scalar slope[];
  foreach () {
    slope[] = (eta[1]-eta[-1])/(2.*Delta);
  }
  clear();
  squares ("slope", linear = true, z = "eta", min = -1, max = 1);
  sprintf (s, "t = %.2f", t);
  draw_string (s, size = 30);
  {
  static FILE * fp = POPEN ("slope", "a");
  save (fp = fp);
  }
  char filename1[50], filename2[50];
  sprintf (filename1, "surface/eta_matrix_%g", t);
  sprintf (filename2, "surface/ux_matrix_%g", t);
  FILE * feta = fopen (filename1, "w");
  // Might need to change to mpi function later
  output_matrix_mpi (eta, feta, N, linear = true);
  fclose (feta);
  FILE * fux = fopen (filename2, "w");
  output_matrix_mpi (u_temp.x, fux, N, linear = true);
  fclose (fux);
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
