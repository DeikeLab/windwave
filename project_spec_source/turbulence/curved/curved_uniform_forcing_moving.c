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

bool limitedAdaptation = 0;
double RELEASETIME = 18;


int MAXLEVEL = dimension == 2 ? 10 : 5; // max level if not use limited refinement
int MAXLEVEL1 = dimension == 2 ? 10 : 7; // max level for finer part
int MAXLEVEL2 = dimension == 2 ? 10 : 7; // max level for coarser part

double BO = 0.845;
double RE = 100.;
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
   Set the regional refinement criteria. */
int refRegion(double x,double y, double z){
  int lev;
  if( y < 1.2 && y > 0)
    lev = MAXLEVEL1;
  /* else if( y > 2 && y < 2*pi) */
  /*  lev = MAXLEVEL1; */
  else if (y >= 1.2 && y < 4)
    lev = MAXLEVEL2;
  else 
    lev = MAXLEVEL2-1;
  return lev;
}

// We might need a diffferent refinement for waves
int refRegion_wave(double x,double y, double z){
  int lev;
  if( y < 1.1 && y > 0)
    lev = MAXLEVEL1+1;
  /* else if( y > 2 && y < 2*pi) */
  /*  lev = MAXLEVEL1; */
  else if (y >= 1.1 && y < 2)
    lev = MAXLEVEL1;
  else if (y >= 2 && y < 4)
    lev = MAXLEVEL2;
  else
    lev = MAXLEVEL2-1;
  return lev;
}

/**
    We need to store the variable forcing term. */
double amp_force = 0.1; //amplitude of the forcing
face vector av[];

int main(int argc, char *argv[]) {
  if (argc > 1)
    RE = atof (argv[1]);
  if (argc > 2)
    MAXLEVEL1 = atoi(argv[2]);
  if (argc > 3)
    MAXLEVEL2 = atoi(argv[3]);
  if (argc > 4)
    UstarRATIO = atof(argv[4]);
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
  u.r[bottom] = dirichlet(0);
  u.n[top] = dirichlet(0);
  u.n[bottom] = dirichlet(0);
  u.t[top] = neumann(0);
  u.t[bottom] = dirichlet(0);
  // Test if setting to neumann change 
  periodic (right);
  periodic (front);
  rho1 = 1.;
  rho2 = RATIO;
  f.sigma = rho1/(BO*sq(k_));
  G.y = -g_;
  c_ = sqrt(g_/k_+f.sigma/rho1*k_);
  Ustar = c_*UstarRATIO;
  mu1 = (2.*pi/k_)*c_*rho1/RE; //using wavelength as length scale
  mu2 = mu1*MURATIO;
  fprintf (ferr, "c = %g, Ustar = %g, nu_w = %g, nu_a = %g\n", c_, Ustar, mu1/rho1, mu2/rho2);
  // Give the address of av to a so that acceleration can be changed
  a = av;
  init_grid (1 << 5);
  // Refine according to 
  uemax = uemaxRATIO*Ustar;
  run();
}


/** 
    Specify the interface shape. */
double WaveProfile(double x, double z) {
  /* double H = 1 + 0.1*sin(4.*x); */
  /* // Write a small amplitude first order wave function */
  /* return H; */
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3 + h_;
}

/**
A function to generate a random field initially. 
*/
double randInRange(int min, int max)
{
  return min + (rand() / (double) (RAND_MAX) * (max - min + 1));
}

event init (i = 0) {
  if (!restore("restart")){
    double rand = 0;
    double ytau = (mu1/rho1)/Ustar;
    do {
      fraction (f, WaveProfile(x,z)-y);
      foreach() {
	rand = randInRange (0,0.1);
	// Initialize with a more accurate profile
	if ((y-WaveProfile(x,z)) > 0.001)
	  u.x[] = (1-f[])*(log((y-WaveProfile(x,z))/ytau)*Ustar/0.41+rand);
	else
	  u.x[] = 0.;
	//u.x[] = (2+rand)*(1.-f[]);  
	u.y[] = 0.;
	u.z[] = 0.;
      }
      boundary ((scalar *){u});
    }
#if TREE  
    // while (adapt_wavelet_limited ({f,u}, (double[]){femax,uemax*10,uemax*10,uemax*10}, refRegion, 5).nf); //if not adapting anymore, return zero
    //while (adapt_wavelet ({f,u}, (double[]){femax,uemax*1,uemax*1,uemax*1}, MAXLEVEL1-1).nf);
    while (0);
#else
    while (0);
#endif

  }
  /* do { */
  /*   fraction (f, WaveProfile(x,z)-y); */
  /*   foreach(){ */
  /*     u.x[] = y; */
  /*     u.y[] = 0.; */
  /*     u.z[] = 0.; */
  /*   } */
  /*   boundary ((scalar *){u}); */
  /* } */
  /* while (adapt_wavelet ({f,u}, (double[]){femax, uemax, uemax, uemax}, MAXLEVEL, MAXLEVEL-2).nf); */
}

/**
  Set the wave velocity 0. */
event set_wave(i=0; i++; t<RELEASETIME) {
  fraction (f, WaveProfile(x,z)-y);
  foreach(){
    // foreach_dimension()
    //   u.x[] = (1.0 - f[])*u.x[];
    u.x[] = (1.0 - f[])*u.x[];
    u.y[] = (1.0 - f[])*u.y[];
    u.z[] = (1.0 - f[])*u.z[];
  }
  boundary ((scalar *){u});
}

/** 
    Start the wave at t=50. */
/* event start(t=55) { */
/*   // A slightly changed version of stokes wave as y=0 at the bottom now so y+h -> y */
/*   fraction (f, WaveProfile(x,z)-y); */
/*   do{ */
/*     scalar Phi[]; */
/*     foreach() { */
/*       double alpa = 1./tanh(k_*h_); */
/*       double a_ = ak/k_; */
/*       double sgma = sqrt(g_*k_*tanh(k_*h_)* */
/* 			 (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)* */
/* 					    (sq(alpa) - 1.) + sq(alpa)))); */
/*       double A_ = a_*g_/sgma; */
/*       double phi1 = A_*cosh(k_*(y))/cosh(k_*h_)*sin(k_*x); */
/*       double phi2 = 3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)* */
/* 	cosh(2.0*k_*(y))*sin(2.0*k_*x)/cosh(2.0*k_*h_); */
/*       double phi3 = 1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)* */
/* 	(9.*sq(alpa) - 13.)*cosh(3.*k_*(y))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*k_*x); */
/*       Phi[] = phi1 + ak*phi2 + ak*ak*phi3; */
/*     }  */
/*     boundary ({Phi}); */
/*     foreach(){ */
/*       foreach_dimension() */
/* 	u.x[] += (Phi[1] - Phi[-1])/(2.0*Delta) * f[]; // f[] is not strictly 0 or 1 I suppose */
/*     } */
/*     boundary ((scalar *){u}); */
/*     /\* do { *\/ */
/*     /\*   fraction (f, WaveProfile(x,z)-y);  *\/ */
/*     /\*   foreach(){  *\/ */
/*     /\* 	u.x[] = y;    *\/ */
/*     /\*     u.y[] = 0.;  *\/ */
/*     /\*     u.z[] = 0.; *\/ */
/*     /\*   }  *\/ */
/*     /\*   boundary ((scalar *){u}); *\/ */
/*     /\* }  *\/ */
/*     /\* while (adapt_wavelet ({f,u}, (double[]){femax, uemax, uemax, uemax}, MAXLEVEL, MAXLEVEL-2).nf);  *\/ */
/*   } */
/* #if TREE   */
/*   while (adapt_wavelet_limited ({f,u}, (double[]){femax,uemax,uemax,uemax}, refRegion, 5).nf); //if not adapting anymore, return zero */
/*   // while (adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL1, 5).nf); */
/* #else */
/*   while (0); */
/* #endif */

/* } */

#include "test/stokes.h"
event start(t = RELEASETIME) {
  // A slightly changed version of stokes wave as y=0 at the bottom now so y+h -> y
  fraction (f, WaveProfile(x,z)-y);
  foreach () {
    u.x[] += u_x(x, y-h_)*f[];
    u.y[] += u_y(x, y-h_)*f[];
  }
  boundary ((scalar *){u});
}

event acceleration (i++) {
  /**
     Forcing term equivalent to pressure gradient in x. */
  double ampl = sq(Ustar)/(L0-h_);
  foreach_face(x)
    av.x[] += ampl*(1.-f[]);
}

/** Output video and field. */
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
    Log turbulent statistics. Since there is a mean flow, Re and ke should be calculated 
    differently from the isotropics.c */

event logfile (i+=100) {
  /* coord ubar; */
  /* foreach_dimension() { */
  /*   stats s = statsf(u.x); */
  /*   ubar.x = s.sum/s.volume; */
  /* } */
  double ke = 0., vd = 0., vol = 0.;
  foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
    vol += dv();
    foreach_dimension() {
      // mean fluctuating kinetic energy
      // ke += dv()*sq(u.x[] - ubar.x);
      // viscous dissipation
      vd += dv()*(sq(u.x[1] - u.x[-1]) +
                  sq(u.x[0,1] - u.x[0,-1]) +
                  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
  }
  // ke /= 2.*vol;
  vd *= mu1/vol;
  static FILE * fd = fopen("stats.dat","a");//"w" before
  if (i == 0) {
    fprintf (ferr, "t dissipation energy Reynolds\n");          //screen
    fprintf (fd, "t dissipation energy Reynolds\n");
  }
  /* fprintf (ferr, "%g %g %g %g\n", */
  /*          t, vd, ke, 2./3.*ke/mu1*sqrt(15.*mu1/vd)); */
  fprintf (fd, "%g %g %g %g\n",
           t, vd, 0, 0);
  fflush(fd);
}

/** Profiling function provided by Antoon. */
/* event profile_maker (t += 5) { */
/*   char file[99]; */
/*   sprintf (file, "prof%g", t); */
/*   vertex scalar phi[]; */
/*   distance_to_surface (f, phi = phi); */
/*   profiles ({u.x, u.y, u.z}, phi, rf = 0.5, */
/* 	    fname = file, min = 0.05, max = 5); */
/* } */


/**
   ## End 
   The wave period is `k_/sqrt(g_*k_)`. We want to run up to 2
   (alternatively 4) periods. */

event end (t = 1000.) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}

event dumpstep (t += 0.1) {
  char dname[100];
  p.nodump = false;
  sprintf (dname, "dump%g", t);
  dump (dname);
}

#if TREE
event adapt (i++) {
  //adapt_wavelet ({f, u}, (double[]){femax, uemax, uemax, uemax}, MAXLEVEL, 5);
  /* if(limitedAdaptation) { */
  /*   if (i > 5) */
  /*     adapt_wavelet_limited ({f,u}, (double[]){femax,uemax,uemax,uemax}, refRegion, 5); */
  /* } */
  /* else */
  if (i == 5)
    fprintf(stderr, "uemaxRATIO = %g\n", uemaxRATIO);
  if (t < RELEASETIME)
    adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL1);
  if (t >= RELEASETIME) {
    foreach () {
      foreach_dimension ()
	u_water.x[] = u.x[]*f[];
    }
    adapt_wavelet ({f,u,u_water}, (double[]){femax,uemax,uemax,uemax,0.001,0.001,0.001}, MAXLEVEL1);
  }
}
/* event adapt2 (i++;t>=55) { */
/*   //adapt_wavelet ({f, u}, (double[]){femax, uemax, uemax, uemax}, MAXLEVEL, 5); */
/*   if(limitedAdaptation) */
/*     adapt_wavelet_limited ({f,u}, (double[]){femax,uemax,uemax,uemax}, refRegion_wave, 5); */
/*   else */
/*     adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL, 5); */
/*  } */
#endif

