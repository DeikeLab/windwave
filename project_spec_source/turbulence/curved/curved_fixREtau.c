//#include "grid/multigrid3D.h"
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "tag.h"
#include "lambda2.h"
#include "navier-stokes/perfs.h"

//#include "adapt_wavelet_limited.h"
#include "sandbox/frac-dist.h"
#include "sandbox/profile6.h"

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

vector u_water[];

/**
    We need to store the variable forcing term. */
double amp_force = 0.1; //amplitude of the forcing
face vector av[];

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
    
  L0 = 2*pi;
  h_ = 1;
  k_ = 4;
  origin (-L0/2., 0, -L0/2.);
  // According to http://basilisk.fr/Basilisk%20C#boundary-conditions
  // for top, u.n = u.y, u.t = u.z, u.r = u.x
  u.r[top] = neumann(0);
  u.r[bottom] = neumann(0);
  u.n[top] = dirichlet(0);
  u.n[bottom] = neumann(0);
  u.t[top] = neumann(0);
  u.t[bottom] = neumann(0);
  // Test if setting to neumann change 
  periodic (right);
  periodic (front);
  rho1 = 1.;
  rho2 = RATIO;
  f.sigma = rho1/(BO*sq(k_));
  G.y = -g_;
  c_ = sqrt(g_/k_+f.sigma/rho1*k_);
  // Ustar = c_*UstarRATIO; // Depleted
  Ustar = 0.25; // Pick a fixed value
  mu2 = Ustar*rho2*(L0-h_)/RE_tau;
  // mu1 = (2.*pi/k_)*c_*rho1/RE; // Depleted: using wavelength as length scale
  mu1 = mu2/(MURATIO);
  RE = rho1*c_*(2*pi/k_)/mu1; // RE now becomes a dependent Non-dim number on c 
  fprintf (stderr, "g = %g, c = %g, Ustar = %g, MURATIO = %g, mu_w = %g, rho_w = %g, mu_a = %g, rho_a = %g, sigma = %g, Bo = %g, RE = %g, Re_tau = %g\n", g_, c_, Ustar, MURATIO, mu1, rho1, mu2, rho2, f.sigma, BO, RE, RE_tau);
  // Give the address of av to a so that acceleration can be changed
  a = av;
  init_grid (1 << 7);
  // Refine according to 
  uemax = uemaxRATIO*Ustar;
  fprintf (stderr, "RELEASETIME = %g, uemax = %g \n", RELEASETIME, uemax);
  run();
}


/** 
    Specify the interface shape. */
double WaveProfile(double x, double z) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  return eta1 + h_;
}

/**
   A function to generate a random field initially. */
double randInRange(int min, int max)
{
  return min + (rand() / (double) (RAND_MAX) * (max - min + 1));
}

event init (i = 0) {
  if (!restore("restart")){
    double rand = 0;
    double ytau = (mu2/rho2)/Ustar;
    do {
      fraction (f, WaveProfile(x,z)-y);
      foreach() {
	// rand = randInRange (0,0.1);
	// Initialize with a more accurate profile
	if ((y-WaveProfile(x,z)) > 0.05)
	  u.x[] = (1-f[])*(log((y-WaveProfile(x,z))/ytau)*Ustar/0.41);
	else
	  u.x[] = 0.;
	// Random noise gets killed by adaptation anyway. We wait for the instability to naturally develop instead.
	// u.x[] = (2+rand)*(1.-f[]); 
	u.y[] = 0.;
	u.z[] = 0.;
      }
      boundary ((scalar *){u});
    }
#if TREE
    // No need for adaptation when starting 
    while (0);
#else
    while (0);
#endif
  }
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

/**
   Start the wave at RELEASETIME. We don't do any adaptation at this step. 
   And we use linear wave instead of stokes. */
// #include "test/stokes.h"
double u_x (double x, double y) {
  return ak*c_*cos(x*k_)*exp(y*k_);
}
double u_y (double x, double y) {
  return ak*c_*sin(x*k_)*exp(y*k_);
}
event start(t = RELEASETIME) {
  // A slightly changed version of stokes wave as y = 0 at the bottom now so y+h -> y
  fraction (f, WaveProfile(x,z)-y);
  foreach () {
    u.x[] += u_x(x, y-h_)*f[];
    u.y[] += u_y(x, y-h_)*f[];
  }
  boundary ((scalar *){u});
}

/**
   Forcing term equivalent to pressure gradient in x. */
event acceleration (i++) {
  double ampl = sq(Ustar)/(L0-h_);
  foreach_face(x)
    av.x[] += ampl*(1.-f[]);
}

/** 
    Output video and field. */
#  define POPEN(name, mode) fopen (name ".ppm", mode)

event movies (t += 0.1) {

  /**
     We first do simple movies of the volume fraction, level of
     refinement fields. In 3D, these are in a $z=0$ cross-section. */

  scalar omega[];
  vorticity (u, omega);
  clear();
  view (fov = 40, camera = "iso", ty = -0.25,
  width = 600, height = 600, bg = {1,1,1}, samples = 4);
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -3.1415);
  squares ("omega", linear = true, n = {1,0,0}, alpha = -3.1415);
  cells (n = {1,0,0}, alpha = -3.1415);
  draw_vof ("f", color = "u.x");
  char s[80];
  sprintf (s, "t = %0.1f", t);
  draw_string (s, size = 30);
  // scalar l2[];
  // lambda2 (u, l2);
  // isosurface ("l2", -1);
  {
    static FILE * fp = POPEN ("3D", "a");
    save (fp = fp);
  }
}

/**
   Generate averaged profile in y direction o nthe fly. */
event profile_output (t += 1) {
  char file[99];
  sprintf (file, "prof_%g", t);
  scalar uxuy[],uxux[],uyuy[],uzuz[];
  foreach () {
    uxuy[] = u.x[]*u.y[];
    uxux[] = u.x[]*u.x[];
    uyuy[] = u.y[]*u.y[];
    uzuz[] = u.z[]*u.z[];
  }
  vertex scalar phi[];
  foreach_vertex ()
    phi[] = y;
  // default phi is y
  profiles ({u.x, u.y, u.z, uxuy, uxux, uyuy, uzuz}, phi, rf = 0.5, fname = file, min = 0.8, max = 2.*pi);
}


/**
   ## End 
   The wave period is `k_/sqrt(g_*k_)`. We want to run up to 2
   (alternatively 4) periods. */

event end (t = 1000.) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}

event dumpstep (t += 5) {
  char dname[100];
  p.nodump = false;
  sprintf (dname, "dump%g", t);
  dump (dname);
}

/** 
    Adaptive function. uemax is tuned. We need a more strict criteria for water speed once the waves starts moving. */ 
#if TREE
event adapt (i++) {
  if (i == 5)
    fprintf(stderr, "uemaxRATIO = %g\n", uemaxRATIO);
  if (t < RELEASETIME)
    adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL);
  if (t >= RELEASETIME) {
    foreach () {
      foreach_dimension ()
	u_water.x[] = u.x[]*f[];
    }
    adapt_wavelet ({f,u,u_water}, (double[]){femax,uemax,uemax,uemax,0.001,0.001,0.001}, MAXLEVEL);
  }
}
#endif

