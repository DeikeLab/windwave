
/**
   # Field scale wave breaking (multilayer solver)

*/

#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "input.h"
#include "output_mpi.h"

/**
   The initial condition is a externally imported wave field. Some controlling parameters. */

#define g_ 9.8
double ETAE = 0.1; // refinement criteria for eta
int MAXLEVEL = 8; // max level of refinement in adapt_wavelet function
int MINLEVEL = 6; // min level
double TEND = 50.; // t end
int NLAYER = 10; // number of layers
int LEVEL_data = 7;

int read_xy_float(char * fname, scalar s, int dlev){
  unsigned long long int size = (1 << (dimension*dlev));
#if _MPI 
  MPI_Win win;
  float * a;
  if (pid() == 0){
    MPI_Win_allocate_shared (size*sizeof(float), sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &a, &win);
  }
  else{ // Slaves obtain the location of the pid()=0 allocated array    
    int disp_unit;
    MPI_Aint  ssize; 
    MPI_Win_allocate_shared (0, sizeof(float), MPI_INFO_NULL,MPI_COMM_WORLD, &a, &win);
    MPI_Win_shared_query (win, 0, &ssize, &disp_unit, &a);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  if (pid() == 0){ 
    MPI_Win_lock (MPI_LOCK_EXCLUSIVE, 0, MPI_MODE_NOCHECK, win);
    FILE * fp = fopen (fname, "rb");
    fread (a, sizeof(float), size,fp);
    MPI_Win_unlock (0, win);
  }
  MPI_Barrier (MPI_COMM_WORLD);
	
  /**
     In serial, life is a bit easier.
  */
#else 
  float * a = (float*) malloc (sizeof(float)*size);
  FILE * fp = fopen (fname, "rb");
  fread (a, sizeof(float), size, fp);
#endif
#if MULTIGRID_MPI  
  int ny = pid() % (mpi_dims[1]);
  int nx = pid() / (mpi_dims[1]);
  int nxoffset = ((1 << depth())*nx);
  int nyoffset = ((1 << depth())*ny);
  fprintf (stderr, "pid = %d, mpi_dims[1] = %d, mpi_dim[0] = %d, nx = %d, ny = %d \n", pid(), mpi_dims[1], mpi_dims[0], nx, ny);
  fflush (stderr);
#else // Non MG-MPI, no offset 
  int nxoffset = 0;
  int nyoffset = 0;
#endif
  unsigned long long int CR = (1 << dlev);
  // o is found to be -3 which is causing the trouble in Antoon's code
  // int o = -BGHOSTS - 1; 
  int o = -BGHOSTS-1;
  /* fprintf(stderr, "o = %d\n", o); */
  unsigned long long int index;
  /** Loading the data itself is now straightforward*/
  // int iternumber = 0;
  foreach(){
    index = ((nxoffset + point.i + o) + (CR*(nyoffset + point.j + o)));
    /* fprintf(ferr, "point.i = %d, point.j = %d, index = %d\n", point.i, point.j, index); */
    s[] = (double)a[index];
    // iternumber += 1;
  }
  // fprintf(ferr, "iteration number = %d", iternumber);  
  return 1;
}

/**
   The adaptive function. */
int my_adapt() {
#if QUADTREE
  /* scalar eta[]; */
  /* foreach() */
  /*   eta[] = h[] > dry ? h[] + zb[] : 0; */
  /* boundary ({eta}); */
  astats s = adapt_wavelet ({eta}, (double[]){ETAE},
			    MAXLEVEL, MINLEVEL);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
#else // Cartesian
  return 0;
#endif
}

/**
   The domain is periodic in $x$ and resolved using 256$^2$
   points and 60 layers. */

int main(int argc, char * argv[])
{
  if (argc > 1)
    NLAYER = atoi(argv[1]);
  if (argc > 2)
    LEVEL_data = atoi(argv[2]);
  if (argc > 3)
    MAXLEVEL = atoi(argv[3]);
  if (argc > 4)
    MINLEVEL = atoi(argv[4]);
  if (argc > 5)
    ETAE = atof(argv[5]);
  if (argc > 6)
    TEND = atof(argv[6]);
  if (argc > 7)
    breaking = atof(argv[7]);
  L0 = 500.;
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  N = 1 << LEVEL_data; // start with a grid of 128
  nl = NLAYER;
  G = g_;
  nu = 1/40000.;
  CFL_H = 0.1; // Smaller time step
  run();
}

/**
   The intial conditions for the free-surface and velocity are given by
   the input field. */


event init (i = 0)
{
  if (!restore ("restart")) {
    /** Import surface position. h of each layer is (H-zb)/nl. */
    /* FILE * feta = fopen ("./pre/eta", "r"); */
    // H is not defined previously
    read_xy_float ("pre/eta", eta, LEVEL_data);
    /* fclose (feta); */
    fprintf(stderr, "Read in eta!\n");
    foreach() {
      zb[] = -5.; // zb[] is a scalar field already defined in header file
      foreach_layer() { 
	h[] = (eta[]-zb[])/nl;
      }
    }
    /** Import velocity. */
    char filename_u[50], filename_v[50], filename_w[50];
    for (int ii=0; ii<nl; ii++) {
      sprintf (filename_u, "pre/u_layer%d", ii);
      sprintf (filename_v, "pre/v_layer%d", ii);
      sprintf (filename_w, "pre/w_layer%d", ii);
      char s[100];
      sprintf (s, "u%d", ii+1);
      vector u_temp = lookup_vector (s);
      read_xy_float (filename_u, u_temp.x, LEVEL_data);
      read_xy_float (filename_v, u_temp.y, LEVEL_data);
      sprintf (s, "w%d", ii+1);
      scalar w_temp = lookup_field (s);
      read_xy_float (filename_w, w_temp, LEVEL_data);
      fprintf(stderr, "Read in velocity, layer index = %d!\n", ii);
    }
    foreach() {
      foreach_layer () 
	zb[] = -5.;
    }
    dump("initial");
  }
}

/** 
    We log the evolution of kinetic and potential energy.
*/

event logfile (i++)
{
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
  static FILE * fp = fopen("energy.dat","a");
  fprintf (fp, "%g %g %g\n", t, ke/2., g_*gpe);
  fflush (fp);
}

/**
   Note that the movie generation below is very expensive. */

int first_frame = 0;
#if 1
event movie (t += 0.1; t <= TEND)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (first_frame == 0) 
    fprintf (fp, "set term pngcairo font ',9' size 1024,768;"
	     "unset key\n"
	     "set pm3d interpolate 4,4 lighting specular 0.6\n"
	     "set zrange [-50:50]\n"
	     "set cbrange [-6:6]\n"
	     "set xlabel 'x'\n"
	     "set ylabel 'y'\n"
	     );
  fprintf (fp,
	   "set output 'plot%05d.png'\n"
	   "set title 't = %.2f'\n"
	   "splot '-' u 1:2:3:4 w pm3d\n",
	   i, t);
  char s[10];
  sprintf (s, "u%d", nl-1);
  vector u_temp = lookup_vector (s);
  output_field ({eta,u_temp.x}, fp, linear = true);
  fprintf (fp, "e\n\n");
  fflush (fp);
  char filename1[50], filename2[50], filename3[50];
  sprintf (filename1, "surface/eta_matrix_%g", t);
  sprintf (filename2, "surface/ux_matrix_%g", t);
  sprintf (filename3, "surface/uy_matrix_%g", t);
  FILE * feta = fopen (filename1, "w");
  // Might need to change to mpi function later
  output_matrix_mpi (eta, feta, N, linear = true);
  fclose (feta);
  FILE * fux = fopen (filename2, "w");
  output_matrix_mpi (u_temp.x, fux, N, linear = true);
  fclose (fux);
  FILE * fuy = fopen (filename3, "w");
  output_matrix_mpi (u_temp.y, fuy, N, linear = true);
  fclose (fuy);
}

/* event closeupmovie (t += 0.1; t <= TEND) */
/* { */
/*   view (fov = 17.3106, quat = {0.475152,0.161235,0.235565,0.832313}, */
/* 	tx = -0.0221727, ty = -0.0140227, width = 1200, height = 768); */
/*   char s[80]; */
/*   sprintf (s, "t = %.2f T0", t/T0); */
/*   draw_string (s, size = 80); */
/*   for (double x = -1; x <= 1; x++) */
/*     translate (x) */
/*       squares ("u59.x", linear = true, z = "eta", min = -2., max = 2.); */
/*   save ("movie.mp4"); */
/* } */

event moviemaker (t = end)
{
  system ("for f in plot*.png; do convert $f ppm:-; done | ppm2mp4 movie.mp4");
}
#endif

#if QUADTREE
event adapt (i++) {
  /* fprintf(stderr, "Adapting start!\n"); */
  /* fflush(stderr); */
  my_adapt();
}
#endif


event dumpstep (t += 1) {
  char dname[100];
  sprintf (dname, "dump%g", t);
  dump (dname);
}


