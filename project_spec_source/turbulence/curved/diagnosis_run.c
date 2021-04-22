#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "lambda2.h"

#define POPEN(name, mode) fopen (name ".ppm", mode)
double snapshot_time = 0;
int OUTLEVEL=7;

double RELEASETIME = 200; 
int MAXLEVEL = dimension == 2 ? 10 : 5; // max level if not use limited refinement

double RE_tau = 180.; // Friction Reynolds number that we don't easily change 
double BO = 200; // Default Bond number
double RE = 100.; // Wave Reynolds number that is dependent on c
double k_ = 1;
double g_ = 1;
double h_ = 1;
double ak = 0.1;
double c_;
double UstarRATIO = 0.5;
double Ustar;

double uemax = 0.01;
double femax = 0.0001;
double uemaxRATIO = 0.01;

#define RATIO 1.225/1000. //density ratio, air to water
#define MURATIO 18.31e-6/10.0e-4 //dynamic viscosity ratio, air to water
// kinematic viscosity air = 16*water
double alter_MU = 1.; //not matching MURATIO extra factor. Default is 1, which 
                      //means MURATIO is used.

vector uwater[];
scalar pair[];

/**
    We need to store the variable forcing term. */
double amp_force = 0.1; //amplitude of the forcing
face vector av[];

void sliceXZ(char * fname,scalar s,double yp, int maxlevel){
  FILE *fpver =fopen (fname,"w"); 
  int nn = (1<<maxlevel);
  double ** field = matrix_new (nn, nn, sizeof(double));
  double stp = L0/(double)nn;
  for (int i = 0; i < nn; i++)
    {
      double xp = stp*i + X0 + stp/2.;
      for (int j = 0; j < nn; j++) 
	{
	  double zp = stp*j + Y0 + stp/2.;
	  Point point = locate (xp, yp,zp);
	  field[i][j] = point.level >= 0 ? s[] : nodata;
	}
    }
  if (pid() == 0){ // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
	fprintf (fpver, "%.9g\t", field[i][j]);
      }
      fputc ('\n', fpver);
    }
    fflush (fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, nn*nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
}

void sliceXY(char * fname,scalar s,double zp, int maxlevel){
  FILE *fpver =fopen (fname,"w"); 
  int nn = (1<<maxlevel);
  double ** field = matrix_new (nn, nn, sizeof(double));
  double stp = L0/(double)nn;
  for (int i = 0; i < nn; i++)
    {
      double xp = stp*i + X0 + stp/2.;
      for (int j = 0; j < nn; j++) 
	{
	  double yp = stp*j + Y0 + stp/2.;
	  Point point = locate (xp, yp,zp);
	  field[i][j] = point.level >= 0 ? s[] : nodata;
	}
    }
  if (pid() == 0){ // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
	fprintf (fpver, "%g\t", field[i][j]);
      }
      fputc ('\n', fpver);
    }
    fflush (fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, nn*nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
}

void output_slice ()
{
  char filename[100];
  int Nslice = 256;
  double L0 = 2*pi;
  double zslice = -L0/2+L0/2./Nslice;
  for (int i=0; i<Nslice; i++) {
    zslice += L0/Nslice;
    sprintf (filename, "./field/ux_t%g_slice%d", snapshot_time, i);
    sliceXY (filename,u.x,zslice,OUTLEVEL);
    sprintf (filename, "./field/uy_t%g_slice%d", snapshot_time, i);
    sliceXY (filename,u.y,zslice,OUTLEVEL);
    sprintf (filename, "./field/f_t%g_slice%d", snapshot_time, i);
    sliceXY (filename,f,zslice,OUTLEVEL);
  } 
  /* double L0 = 2*pi; */
  /* double yslice = pi/Nslice; */
  /* for (int i=0; i<Nslice; i++) { */
  /*   yslice += L0/Nslice; */
  /*   sprintf (filename, "./field/ux_t%g_slice%d", snapshot_time, i); */
  /*   sliceXZ (filename,u.x,yslice,OUTLEVEL); */
  /*   sprintf (filename, "./field/uy_t%g_slice%d", snapshot_time, i); */
  /*   sliceXZ (filename,u.y,yslice,OUTLEVEL); */
  /*   sprintf (filename, "./field/f_t%g_slice%d", snapshot_time, i); */
  /*   sliceXZ (filename,f,yslice,OUTLEVEL); */
  /* }   */
}

void output_twophase () {
  scalar pos[];
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
  position (f, pos, G, Z);
  char etaname[100];
  sprintf (etaname, "./eta/eta_t%g_%d", snapshot_time, pid());
  FILE * feta = fopen (etaname, "w");
  fprintf(feta, "x,z,pos,p,p_p1,p_p2\n");
  // printing out quantities: p_p1 for p at plus 1, p_m1 for p at minus 1 etc.
  foreach(){
    if (interfacial (point, f)){
      fprintf (feta, "%g,%g,%g,%g,%g,%g\n", 
	       x, z, pos[], p[], p[0,1], p[0,2]);
      /* tau.x[], tau.y[], u.x[], u.y[], n.x/norm_2, n.y/norm_2); */
    }
  }
  fclose (feta);
} 

int main(int argc, char *argv[]) {
  if (argc > 1)
    RE_tau  = atof (argv[1]);
  if (argc > 2)
    BO = atof(argv[2]);
  if (argc > 3)
    MAXLEVEL = atoi(argv[3]);
  if (argc > 4)
    g_ = atof(argv[4]);
  if (argc > 5)
    ak = atof(argv[5]);
  if (argc > 6)
    RELEASETIME = atof(argv[6]);
  if (argc > 7)
    uemaxRATIO = atof(argv[7]);
  if (argc > 8)
    alter_MU = atof(argv[8]);
  if (argc > 9)
    snapshot_time = atof(argv[9]);
  if (argc > 10)
    OUTLEVEL = atoi(argv[10]);

  L0 = 2.*pi;
  h_ = 1; // Water depth
  k_ = 4; // Four waves per box
  origin (-L0/2., 0, -L0/2.);
  // According to http://basilisk.fr/Basilisk%20C#boundary-conditions
  // for top, u.n = u.y, u.t = u.z, u.r = u.x
  u.r[top] = neumann(0);
  u.r[bottom] = neumann(0);
  u.n[top] = dirichlet(0); // This is supposed to be neumann 
  u.n[bottom] = neumann(0);
  u.t[top] = neumann(0);
  u.t[bottom] = neumann(0);
  // TO-DO: Test if setting to neumann change 
  periodic (right);
  periodic (front);
  rho1 = 1.;
  rho2 = RATIO;
  f.sigma = g_*rho1/(BO*sq(k_));
  G.y = -g_;
  c_ = sqrt(g_/k_+f.sigma/rho1*k_);
  // Ustar = c_*UstarRATIO; // Depleted
  Ustar = 0.25; // Pick a fixed value
  mu2 = Ustar*rho2*(L0-h_)/RE_tau;
  // mu1 = (2.*pi/k_)*c_*rho1/RE; // Depleted: using wavelength as length scale
  mu1 = mu2/(MURATIO)/alter_MU;
  RE = rho1*c_*(2*pi/k_)/mu1; // RE now becomes a dependent Non-dim number on c 
  // Give the address of av to a so that acceleration can be changed
  a = av;
  init_grid (1 << 7);
  // Refine according to 
  uemax = uemaxRATIO*Ustar;
  //output_twophase ();
  //output_slice();
  //output();
  //movies();
  run();
}

int phony = 1;
int j = 0;

event init (i=0) {
  char targetname[100];
  sprintf (targetname, "dump%g", snapshot_time);
  if (!restore (targetname)) {
    fprintf(ferr, "Not restored!\n");
    return 1;
  }
  output_slice();
  output_twophase();
  phony = 0;
  fprintf (stderr, "Output from dump done! \n");
}

double WaveProfile (double x, double y)
{
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + eta2 + eta3 + h_;
}

/**
  Set the wave velocity 0. */
event set_wave(i=0; i++; t<RELEASETIME) {
  fraction (f, WaveProfile(x,z)-y);
  foreach(){
    u.x[] = (1.0 - f[])*u.x[];
    u.y[] = (1.0 - f[])*u.y[];
    u.z[] = (1.0 - f[])*u.z[];
  }
  boundary ((scalar *){u});
}

event acceleration (i++) {
  double ampl = sq(Ustar)/(L0-h_);
  foreach_face(x)
    av.x[] += ampl*(1.-f[]);
}

event pair_out (i++){
  if (i == j+5) {
    fprintf (stderr, "i = %d, start figure! \n", i);
    foreach () {
      pair[] = p[]*(1.-f[]);
    }
    foreach () {
      uwater.x[] = u.x[]*f[];
    }
    view (ty = -0.5);
    squares ("uwater.x", n = {0,0,1}, alpha = 0.1);
    {
      static FILE * fp = POPEN ("uwater", "a");
      save (fp = fp);
    }
    squares ("p", n = {0,0,1}, alpha = -0.1);
    {
      static FILE * fp = POPEN ("p", "a");
      save (fp = fp);
    }
    squares ("pair", n = {0,0,1}, alpha = -0.1);
    {
      static FILE * fp = POPEN ("pair", "a");
      save (fp = fp);
    }
    clear ();
    view (ty = -0.2, fov = 8);
    cells (n = {0,0,1}, alpha = -0.1);
    squares ("level", n = {0,0,1}, alpha = -0.1, max = 10, min = 5);
    {
      static FILE * fp = POPEN ("level", "a");
      save (fp = fp);
    }
    fprintf (stderr, "i = %d, end figure! \n", i);
    char filename[100];
    int Nslice = 256;
    double L0 = 2*pi;
    double zslice = -L0/2+L0/2./Nslice;
    for (int i=0; i<Nslice; i++) {
      zslice += L0/Nslice;
      sprintf (filename, "./field/p_run_t%g_slice%d", snapshot_time, i);
      sliceXY (filename,p,zslice,OUTLEVEL);
      sprintf (filename, "./field/pair_run_t%g_slice%d", snapshot_time, i);
      sliceXY (filename,pair,zslice,OUTLEVEL);
    }
    return 1;
  }
}

#if TREE
event adapt (i++) {
  if (i == 5)
    fprintf(stderr, "uemaxRATIO = %g\n", uemaxRATIO);
  if (t < RELEASETIME)
    adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL);
  if (t >= RELEASETIME) {
    foreach () {
      foreach_dimension ()
	uwater.x[] = u.x[]*f[];
    }
    adapt_wavelet ({f,u,uwater}, (double[]){femax,uemax,uemax,uemax,0.001,0.001,0.001}, MAXLEVEL);
  }
  /**  Ad hoc way of going around resetting the index i */
  if (phony == 0) {
    j = i;
    phony = 1;
    fprintf (stderr, "New indexing set to j = %d!\n", j);
  }
}
#endif

event end (t = 1000) {
  fprintf (stderr, "i = %d t = %g\n", i, t);
  dump ("end");
}
