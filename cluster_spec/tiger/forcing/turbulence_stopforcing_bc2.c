#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "tag.h"
#include "lambda2.h"


int MAXLEVEL = dimension == 2 ? 10 : 7;

double BO = 200.;
double RE = 1000.;
double k_ = 2.*pi;
double h_ = 0.5;
double g_ = 1.;

double uemax = 0.001;
double femax = 0.0001;

#define BO 1000.0 //Bond number
#define RE 40000.0 //Reynolds number
#define RATIO 1.0/850.0 //density ratio, water to air
#define MURATIO 17.4e-6/8.9e-4 //dynamic viscosity ratio, water to air

/**
    We need to store the variable forcing term. */
double amp_force = 0.1; //amplitude of the forcing
face vector av[];

int main() {
  L0 = 2 ;
  origin (-L0/2., -L0/2., -L0/2.);
  // According to http://basilisk.fr/Basilisk%20C#boundary-conditions
  // for top, u.n = u.y, u.t = u.z, u.r = u.x
  u.r[top] = neumann (1);
  u.r[bottom] = dirichlet (y);
  u.n[top] = dirichlet(0);
  u.n[bottom] = dirichlet(0);
  u.t[top] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
  // Test if setting to neumann change 
  periodic (right);
  periodic (front);
  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0/RE; //using wavelength as length scale
  mu2 = 1.0/RE*MURATIO;
  f.sigma = 1./(BO*sq(k_));
  G.y = -g_;

  a = av;

  init_grid (1 << (MAXLEVEL-2));
  run();
}


/** 
    We use the moving cylinder example to inform the implementation
    of the beach.*/
double WaveProfile(double x, double z) {
//Define beach profile H(x)
  double H = 0.0;
  // Write a small amplitude first order wave function
  return H;
}

/**
  Set the wave velocity 0. */
event set_wave(i=0;i++) {
  /** 
      And set the corresponding momentum B.C.s within the beach region:*/
  foreach(){
    // foreach_dimension()
    //   u.x[] = (1.0 - f[])*u.x[];
    u.x[] = (1.0 - f[])*u.x[];
    u.y[] = (1.0 - f[])*u.y[];
    u.z[] = (1.0 - f[])*u.z[];
  }
  boundary ((scalar *){u});
}

event init (i = 0) {
  if (!restore("restart")){
    double femax = 2e-4;
    double uemax = 2e-4;
    do {
  	  fraction (f, WaveProfile(x,z)-y);
  	  foreach(){
        u.x[] = y;  
        u.y[] = 0.;
        u.z[] = 0.;
      }
  	  boundary ((scalar *){u});
    }
    while (adapt_wavelet ({f,u}, (double[]){femax, uemax, uemax, uemax}, MAXLEVEL, 5).nf);
  }
}



/** Output video and field. */
#  define POPEN(name, mode) fopen (name ".ppm", mode)

event movies (t += k_/sqrt(g_*k_)/32) {

  /**
     We first do simple movies of the volume fraction, level of
     refinement fields. In 3D, these are in a $z=0$ cross-section. */

  {
    static FILE * fp = POPEN ("f", "w");
    output_ppm (f, fp, min = 0, max = 1, n = 512);
  }

#if TREE
  {
    scalar l[];
    foreach()
      l[] = level;
    static FILE * fp = POPEN ("level", "w");
    output_ppm (l, fp, min = 5, max = MAXLEVEL, n = 512);
  }
#endif

  scalar omega[];
  vorticity (u, omega);
  view (width = 1200, height = 1200);
  clear();
  squares ("u.x", linear = true, n = {0,0,1}, alpha = 0);
  {
    static FILE * fp = POPEN ("Ux", "w");
    save (fp = fp);
  }
  squares ("u.y", linear = true, n = {0,0,1}, alpha = 0);
  {
    static FILE * fp = POPEN ("Uy", "w");
    save (fp = fp);
  }
  squares ("u.z", linear = true, n = {1,0,0}, alpha = 0);
  {
    static FILE * fp = POPEN ("Uz", "w");
    save (fp = fp);
  }
  squares ("omega", linear = true, n = {0,0,1}, alpha = 0);
  {
    static FILE * fp = POPEN ("omega", "w");
    save (fp = fp);
  }
  char dname[100];
  sprintf (dname, "field%g", t/(k_/sqrt(g_*k_)));
  FILE * fp = fopen(dname, "w");
  output_field({u.x,u.y,u.z,omega,f}, fp, n = 128, linear = true);
  fclose (fp);
}

event vortex (t += 0.25)
{
  view (fov = 44, camera = "iso", ty = .2,
  width = 600, height = 600, bg = {1,1,1}, samples = 4);
  clear();
  squares ("u.y", linear = true);
  squares ("u.x", linear = true, n = {1,0,0});
  scalar omega[];
  vorticity (u, omega);
  squares ("omega", linear = true, n = {0,1,0});
  scalar l2[];
  lambda2 (u, l2);
  isosurface ("l2", -1);
  {
    static FILE * fp = POPEN ("3Dvortex", "w");
    save (fp = fp);
  }
}

event logfile (i++) {
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }

  double ke = 0., vd = 0., vol = 0.;
  foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
    vol += dv();
    foreach_dimension() {
      // mean fluctuating kinetic energy
      ke += dv()*sq(u.x[] - ubar.x);
      // viscous dissipation
      vd += dv()*(sq(u.x[1] - u.x[-1]) +
                  sq(u.x[0,1] - u.x[0,-1]) +
                  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
  }
  ke /= 2.*vol;
  vd *= mu1/vol;
  static FILE * fd = fopen("stats.dat","a");//"w" before
  if (i == 0) {
    fprintf (ferr, "t dissipation energy Reynolds\n");          //screen
    fprintf (fd, "t dissipation energy Reynolds\n");
  }
  fprintf (ferr, "%g %g %g %g\n",
           t, vd, ke, 2./3.*ke/mu1*sqrt(15.*mu1/vd));
  fprintf (fd, "%g %g %g %g\n",
           t, vd, ke, 2./3.*ke/mu1*sqrt(15.*mu1/vd));
}

/**
   ## End 
   The wave period is `k_/sqrt(g_*k_)`. We want to run up to 2
   (alternatively 4) periods. */

event end (t = 200.*k_/sqrt(g_*k_)) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}

event dumpstep (t += k_/sqrt(g_*k_)/32) {
  char dname[100];
  sprintf (dname, "dump%g", t/(k_/sqrt(g_*k_)));
  dump (dname);
}


event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){femax, uemax, uemax, uemax}, MAXLEVEL, 5);
}
