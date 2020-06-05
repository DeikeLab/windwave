//#include "grid/multigrid3D.h"
#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
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
double g_ = 9.8;

double uemax = 0.0001;
double femax = 0.0001;

#define RATIO 1.225/1000. //density ratio, air to water
#define MURATIO 18.31e-6/10.0e-4 //dynamic viscosity ratio, air to water
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

face vector av[];
#define MU 0.0004
double ue = 0.05, cse = 0.01; //Refinement criteria 

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
  const face vector muc[] = {MU,MU,MU};
  mu = muc;
  // Give the address of av to a so that acceleration can be changed
  a = av;
  init_grid (1 << (MAXLEVEL+1));
  run();
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
    vertex scalar phi[];
    foreach_vertex()
      phi[] = -1. - 0.1*sin(4*x) + y;
    fractions (phi, cs, fs);
    foreach() {
        rand = randInRange (0,0.1);
        u.x[] = (1+rand)*(cs[]);  
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

event acceleration (i++) {
  /**
     Forcing term equivalent to pressure gradient in x. */
  foreach_face(x)
    av.x[] += 0.01;
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

event logfile2 (i+=100) {
  static FILE * fp1 = fopen("DNS.dat", "a");
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
  double vd;
  double ke = 0., vdm = 0., vol = 0.;
  foreach(reduction(+:ke) reduction(+:vdm) reduction(+:vol)) {
    vol += dv();
    foreach_dimension() {
      // mean fluctuating kinetic energy    
      ke += dv()*sq(u.x[] - ubar.x);
      // viscous dissipation                                                        
      vd = dv()*(sq(u.x[1] - u.x[-1]) + sq(u.x[0,1] - u.x[0,-1]) + sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
      vdm += MU*vd;
    }
  }
  ke /= 2.*vol;
  vdm /= vol;
  if (i == 0)
    fprintf (fp1, "t dissipationt energy perf.t perf.speed\n");
  fprintf (fp1, "%g %g %g %g %g\n", t, vdm, ke, perf.t, perf.speed);
  fflush (fp1);
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
  p.nodump = false;
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
