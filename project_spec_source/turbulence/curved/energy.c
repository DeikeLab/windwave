/**
   Compute various quantities related to energy budget.
   Looping over dump files is achieved with a python script.  */

// #include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "tag.h"
#include "curvature.h"
#include "vof.h"

#  define POPEN(name, mode) fopen (name ".ppm", mode)

double uemax = 0.01;
double femax = 0.00001;
int LEVEL = dimension == 2 ? 11 : 7;

scalar omega[];

double snapshot_time; // Current snapshot time
double omega0; // The initial omega value

// The code needs to know a few parameters to 1. Compute Utop; 2. Compute omega0
double BO = 200.;
double RE = 40000.;
double RATIO = (1./850.);
double MURATIO = (17.4e-6/8.9e-4);
double k_ = 2.*pi;
double h_ = 0.5;
double g_ = 1.;
double UstarRATIO = 1;   // Ratio between Ustar and c                                                                        
double Ustar;
double Utop;
double omega_wave;  // wave frequency, for normalization

void dissipation () {
  /* coord ubar; */
  /* foreach_dimension() { */
  /*   stats s = statsf(u.x); */
  /*   ubar.x = s.sum/s.volume; */
  /* } */
  char dissname[100] = "dissipation.dat";
  FILE * fp = fopen (dissname, "a");
  if (snapshot_time == 0) {
    fprintf(fp, "t,diss\n");
  }
  double ke = 0., vd = 0., vol = 0.;
  foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
    vol += dv()*f[]; // Only count the water phase
    foreach_dimension() {
      // mean fluctuating kinetic energy
      // ke += dv()*sq(u.x[] - ubar.x);
      // viscous dissipation
      vd += sqrt(f[])*f[]*dv()*(sq(u.x[1] - u.x[-1]) +
		      sq(u.x[0,1] - u.x[0,-1]))/sq(2.*Delta);
    }
  }
  // ke /= 2.*vol;
  vd *= mu1;

  /* if (i == 0) */
  /*   fprintf (stderr, "t dissipation energy Reynolds\n"); */
  fprintf (fp, "%g,%g \n", snapshot_time, vd);
}

/**
   Output the interface related quantities. For serial only. */

void output_twophase_locate () {
  scalar pos[];
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
  position (f, pos, G, Z);
  char etaname[100];
  sprintf (etaname, "./field/eta_loc_t%g", snapshot_time);
  FILE * feta = fopen (etaname, "w");
  fprintf(feta, "x,pos,ux,ux_p1,ux_p2,uy,epsilon,delta,p,p_p1,dudy1,dudy2,dvdx1,dvdx2,dudx1,dudx2,dvdy1,dvdy2,uxw,uyw\n");
  // printing out quantities: p_p1 for p at plus 1, p_m1 for p at minus 1 etc.
  double stp = L0/1024.;
  foreach(){
    if (interfacial (point, f)){
      if (point.level == 11) {
	coord n = mycs (point, f);
	double norm_2 = sqrt(sq(n.x) + sq(n.y));
	fprintf (feta, "%g,%g,%g,%g,%g,%g,", x, pos[], u.x[],u.x[0,1],u.x[0,2],u.y[]);
	/** Find the airside stress. */
 	double yp = y + stp*2;
	point = locate (x, yp);
	if (point.level > 0) {
	  POINT_VARIABLES;
	  double dudy1 = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
	  double dudy2 = (u.x[0,2] - u.x[0,0])/(2.*Delta);
	  // double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
	  double dvdx1 = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
	  double dvdx2 = (u.y[1,1] - u.y[1,-1])/(2.*Delta);
	  double dudx1 = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
	  double dudx2 = (u.x[1,1] - u.x[1,-1])/(2.*Delta);
	  double dvdy1 = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
	  double dvdy2 = (u.y[0,2] - u.y[0,0])/(2.*Delta);
	  /*  tau.x[] = 2*mu2*(SDeform.x.x[]*n.x + SDeform.y.x[]*n.y)/norm_2; */
	  /*  tau.y[] = 2*mu2*(SDeform.x.y[]*n.x + SDeform.y.y[]*n.y)/norm_2;  */
	  /* double tau1 = 2*mu2*(dudy1 + dvdx); */
	  /* double tau2 = 2*mu2*(dudy2 + dvdx); */
	  fprintf (feta, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,", \
		   -n.x/n.y, Delta, p[0,0], p[0,1], dudy1, dudy2, dvdx1, dvdx2, \
		   dudx1, dudx2, dvdy1, dvdy2);
	/* tau.x[], tau.y[], u.x[], u.y[], n.x/norm_2, n.y/norm_2); */
	}

	/** Find the water side velocity. */
 	yp = y - stp;
	point = locate (x, yp);
	if (point.level > 0) {
	  fprintf (feta, "%g,%g\n", u.x[], u.y[]);
	}
      }
    }
  }
  fflush (feta);
  fclose (feta);
}

int phony = 1;
int j = 0;

int main (int argc, char * argv[])
{
  if (argc > 1)
    snapshot_time = atof(argv[1]);
  if (argc > 2)
    BO = atof(argv[2]);
  if (argc > 3)
    RE = atof(argv[3]);
  if (argc > 4)
    UstarRATIO = atof(argv[4]);
  if (argc > 5)
    TOLERANCE = atof(argv[5]); // Tolerance for Poisson solver
  
  // A few field information, copied from the simulation code
  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0/RE;                                                                      
  mu2 = 1.0/RE*MURATIO;
  f.sigma = 1./(BO*sq(k_));
  G.y = -g_;
  Ustar = sqrt(g_/k_+f.sigma*k_)*UstarRATIO;
  omega0 = sq(Ustar)/(mu2/rho2);
  Utop = omega0*L0/2.;
  // fprintf(ferr, "omega0 = %g, Utop = %g\n", omega0, Utop);  
  omega_wave = sqrt(g_*k_ + f.sigma*k_*k_*k_/rho1);
  fprintf(ferr, "omega_wave = %g, Utop = %g\n", omega_wave, Utop);
  
  // boundary() so that ghost cells have correct values
  // origin (-L0/2, -L0/2, -L0/2);
  periodic (right);
  u.n[top] = dirichlet(0);
  u.t[top] = dirichlet(Utop);

  // TOLERANCE = 1e-7;
  // restore (targetname);
  // boundary((scalar *){u});
  run ();
}

event init (i=0) {
  char targetname[1000];
  sprintf (targetname, "dump%.7g", snapshot_time);
  fprintf(ferr, "%s ", targetname);
  if (!restore (targetname)) {
    fprintf(ferr, "Not restored!\n");
    return 1;
  }
  phony = 0;
}

event adapt (i++) {
  if (i == j+5) {
    fprintf(stderr, "Output running!\n");
    dissipation();
    output_twophase_locate();
    return 1;
  }
  /**  Ad hoc way of going around resetting the index i */
  if (phony == 0) {
    j = i;
    phony = 1;
    fprintf (stderr, "New indexing set to j = %d!\n", j);
  }
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, LEVEL, 5);
}

event end_simu (t = 100) {
  dump ("end");
}


