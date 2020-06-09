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

#include "adapt_wavelet_limited.h"

bool limitedAdaptation = 1;

int MAXLEVEL = dimension == 2 ? 10 :7; // max level if not use limited refinement
int MAXLEVEL1 = dimension == 2 ? 10 : 7; // max level for finer part
int MAXLEVEL2 = dimension == 2 ? 10 : 7; // max level for coarser part

double BO = 200.;
double RE = 100.;
double k_ = 1;
double g_ = 1.;

double uemax = 0.0001;
double femax = 0.0001;

#define RATIO 1.0/850.0 //density ratio, air to water
#define MURATIO 17.4e-6/8.9e-4 //dynamic viscosity ratio, air to water
// kinematic viscosity air = 16*water

/**
   Set the regional refinement criteria. */
int refRegion(double x,double y, double z){
int lev;
if( y < 2 && y > 0)
   lev = MAXLEVEL1;
  /* else if( y > 2 && y < 2*pi) */
  /*  lev = MAXLEVEL1; */
 else
   lev = MAXLEVEL2;
return lev;
}

/**
    We need to store the variable forcing term. */
double amp_force = 0.001; //amplitude of the forcing
face vector av[];

int main(int argc, char *argv[]) {
  if (argc > 1)
    RE = atof (argv[1]);
  if (argc > 2)
    MAXLEVEL1 = atoi(argv[2]);
  if (argc > 3)
    MAXLEVEL2 = atoi(argv[3]);

  L0 = 2*pi ;
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
  mu1 = 2*pi/RE; //using wavelength as length scale
  mu2 = mu1*MURATIO;
  f.sigma = 1./(BO*sq(k_));
  G.y = -g_;
  // Give the address of av to a so that acceleration can be changed
  a = av;
  init_grid (1 << (MAXLEVEL+1));
  run();
}


/** 
    Specify the interface shape. */
double WaveProfile(double x, double z) {
  double H = 1.0;
  // Write a small amplitude first order wave function
  return H;
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
    fraction (f, WaveProfile(x,z)-y);
    foreach() {
        rand = randInRange (0,0.1);
        u.x[] = 1+rand;  
        u.y[] = 0.;
        u.z[] = 0.;
      }
    boundary ((scalar *){u});
    /* do { */
    /*   fraction (f, WaveProfile(x,z)-y);  */
    /*   foreach(){  */
    /* 	u.x[] = y;    */
    /*     u.y[] = 0.;  */
    /*     u.z[] = 0.; */
    /*   }  */
    /*   boundary ((scalar *){u}); */
    /* }  */
    /* while (adapt_wavelet ({f,u}, (double[]){femax, uemax, uemax, uemax}, MAXLEVEL, MAXLEVEL-2).nf);  */
  }
}

/**
  Set the wave velocity 0. */
event set_wave(i=0;i++) {
  foreach(){
    // foreach_dimension()
    //   u.x[] = (1.0 - f[])*u.x[];
    u.x[] = (1.0 - f[])*u.x[];
    u.y[] = (1.0 - f[])*u.y[];
    u.z[] = (1.0 - f[])*u.z[];
  }
  boundary ((scalar *){u});
}

event acceleration (i++) {
  /**
     Forcing term equivalent to pressure gradient in x. */
  foreach_face(x)
    av.x[] += 0.001*(1.-f[]);
}

/** Output video and field. */
#  define POPEN(name, mode) fopen (name ".ppm", mode)

event movies (t += 1.) {

  /**
     We first do simple movies of the volume fraction, level of
     refinement fields. In 3D, these are in a $z=0$ cross-section. */

  /* { */
  /*   static FILE * fp = POPEN ("f", "w"); */
  /*   output_ppm (f, fp, min = 0, max = 1, n = 512); */
  /* } */

/* #if TREE */
/*   { */
/*     scalar l[]; */
/*     foreach() */
/*       l[] = level; */
/*     static FILE * fp = POPEN ("level", "w"); */
/*     output_ppm (l, fp, min = 5, max = MAXLEVEL, n = 512); */
/*   } */
/* #endif */

  scalar omega[];
  vorticity (u, omega);
  /* clear(); */
  /* view (width = 1200, height = 1200); */  
  /* squares ("u.x", linear = true, n = {0,0,1}, alpha = 0); */
  /* { */
  /*   static FILE * fp = POPEN ("Ux", "w"); */
  /*   save (fp = fp); */
  /* } */
  /* squares ("u.y", linear = true, n = {0,0,1}, alpha = 0); */
  /* { */
  /*   static FILE * fp = POPEN ("Uy", "w"); */
  /*   save (fp = fp); */
  /* } */
  /* squares ("u.z", linear = true, n = {1,0,0}, alpha = 0); */
  /* { */
  /*   static FILE * fp = POPEN ("Uz", "w"); */
  /*   save (fp = fp); */
  /* } */
  /* squares ("omega", linear = true, n = {0,0,1}, alpha = 0); */
  /* { */
  /*   static FILE * fp = POPEN ("omega", "w"); */
  /*   save (fp = fp); */
  /* } */
  /* clear(); */
  /* squares("level"); */
  /* cells(); */
  /* { */
  /*   static FILE * fp = POPEN("level", "w"); */
  /*   save (fp = fp); */
  /* } */
  clear();
  view (fov = 44, camera = "iso", ty = .2,
  width = 600, height = 600, bg = {1,1,1}, samples = 4);
  squares ("u.y", linear = true);
  squares ("u.x", linear = true, n = {1,0,0});
  squares ("omega", linear = true, n = {0,1,0});
  scalar l2[];
  lambda2 (u, l2);
  isosurface ("l2", -1);
  {
    static FILE * fp = POPEN ("3Dvortex", "w");
    save (fp = fp);
  }
}


event logfile (i+=100) {
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
  fflush(fd);
}

/**
   ## End 
   The wave period is `k_/sqrt(g_*k_)`. We want to run up to 2
   (alternatively 4) periods. */

event end (t = 1000.) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}

event dumpstep (t += 10) {
  char dname[100];
  sprintf (dname, "dump%g", t);
  dump (dname);
}


/* event adapt (i++) { */
/*   //adapt_wavelet ({f, u}, (double[]){femax, uemax, uemax, uemax}, MAXLEVEL, 5);  */
/*   if(limitedAdaptation) */
/*     adapt_wavelet_limited ({f,u}, (double[]){femax,uemax,uemax,uemax}, refRegion, 5); */
/*   else */
/*     adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL, 5); */
/*  } */
