/**
   # Breaking wave

   We solve the two-phase Navier--Stokes equations with surface tension
   and using a momentum-conserving transport of each phase. Gravity is
   taken into account using the "reduced gravity approach" and the
   results are visualised using Basilisk view. */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "tag.h"
#include "lambda2.h"

/**
   We log some profiling information. */

// #include "navier-stokes/perfs.h"
// #include "profiling.h"

/**
   The primary parameters are the wave steepness $ak$, the Bond and
   Reynolds numbers. */

double ak = 0.05;
double BO = 200.;
double RE = 40000.;

/**
   The default maximum level of refinement depends on the dimension. */

int LEVEL = dimension == 2 ? 10 : 7;
/**
   The error on the components of the velocity field used for adaptive
   refinement. */

double uemax = 0.001;
double femax = 0.00001;

/**
   The density and viscosity ratios are those of air and water. */

double RATIO = (1./850.);
double MURATIO = (17.4e-6/8.9e-4);

/**
   Define if we want to use a Dirac viscous layer initialization. */
int DIRAC = 0;

/**
   Define the velocity output file number. */
int countvelocity = 0;

/**
   The wave number, fluid depth and acceleration of gravity are set to
   these values. */
double k_ = 2.*pi;
double h_ = 0.5;
double g_ = 1.;

/** define the wind profile related parameters. */

double m = 5.;  // vary between 5 and 8
double B = 0.;
double Karman = 0.41;   // Karman universal turbulence constant
double UstarRATIO = 1;   // Ratio between Ustar and c



/**
   The program takes optional arguments which are the level of
   refinement, steepness, Bond and Reynolds numbers. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  if (argc > 2)
    ak = atof(argv[2]);
  if (argc > 3)
    BO = atof(argv[3]);
  if (argc > 4)
    RE = atof(argv[4]);   
  if (argc > 5)
    m = atof(argv[5]);
  if (argc > 6)
    B = atof(argv[6]);
  if (argc > 7)
    UstarRATIO = atof(argv[7]);
  if (argc > 8)
    DIRAC = atof(argv[8]);
  

  /**
     The domain is a cubic box centered on the origin and of length
     $L0=1$, periodic in the x- and z-directions. */
   
  origin (-L0/2, -L0/2, -L0/2);
  periodic (right);
  u.n[top] = dirichlet(0);
  u.t[top] = neumann(0);
#if dimension > 2
  periodic (front);
#endif

  /**
     Here we set the densities and viscosities corresponding to the
     parameters above. Note that these variables are defined in two-phase.h already.*/
 
  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0/RE; //using wavelength as length scale
  mu2 = 1.0/RE*MURATIO;
  f.sigma = 1./(BO*sq(k_));
  G.y = -g_;

  /**
     When we use adaptive refinement, we start with a coarse mesh which
     will be refined as required when initialising the wave. */
  
#if TREE  
  N = 32; // If coarsened or not
#else
  N = 1 << LEVEL;
#endif
  run();
}

/**
   ## Initial conditions

   These functions return the shape of a third-order Stokes wave with the
   wavenumber and steepness given by the parameters above ($ak$ and
   $_k_$). */

double wave (double x, double y) {

  return 0 - y;
}

double eta (double x, double y) {

  return 0;
}
/**

   We also calculate an approximation to a Dirac distribution on the wave surface.
   This allows us to calculate a vortex sheet on the surface to provide a boundary
   layer in the air above the water surface. */

double gaus (double y, double yc, double T){
  double deltaw = sqrt(2.0/RE)/k_;
  double deltaa = sqrt(2.0/RE*MURATIO/RATIO)/k_;
  double r = y - yc;
  return 2.0/(sqrt(2.0*pi*sq(deltaa)) + sqrt(2.0*pi*sq(deltaw))) *
    (T*exp(-sq(r)/(2.0*sq(deltaw))) + (1.0 - T)*exp(-sq(r)/(2.0*sq(deltaa))));
}


/**
   We either restart (if a "restart" file exists), or initialise the wave
   using the third-order Stokes wave solution. */
event init (i = 0)
{

  // calculate profile related info
  double Ustar = sqrt(g_/k_)*UstarRATIO;
  double y1 = m*mu2/rho2/Ustar;
  double Udrift = B*Ustar;   
  fprintf(stderr, "UstarRATIO=%g B=%g\n m=%g ak=%g", UstarRATIO, B, m, ak);

  if (!restore ("restart")) {
    do {
      fraction (f, wave(x,y));

      /**
	 To initialise the velocity field, we first define the potential. */
      
	foreach(){
	  foreach_dimension()
	    u.x[] = 0; // f[] is not strictly 0 or 1 I suppose
	}
	foreach(){
	  if ((y-eta(x,y))<y1){
	    u.x[] += Udrift + sq(Ustar)/(mu2/rho2)* (y-eta(x,y)) * (1-f[]);
	  }
	  else{
	    double beta = 2*Karman*Ustar/mu2*rho2*((y-eta(x,y))-y1);
	    double alpha = log(beta+sqrt(sq(beta)+1));
	    double tanhtemp = (exp(alpha/2)-exp(-alpha/2))/(exp(alpha/2)+exp(-alpha/2));
	    u.x[] += (Udrift + m*Ustar + Ustar/Karman*(alpha-tanhtemp)) * (1-f[]);
	  }
	}
        //fprintf(stderr, "Added line running for each cell!");
        //fprintf(stderr, "Added line running!");
      boundary ((scalar *){u});  // type casting   
    }
    /**
       On trees, we repeat this initialisation until mesh adaptation does
       not refine the mesh anymore. */

#if TREE  
    while (adapt_wavelet ({f,u},
			  (double[]){0.01,0.05,0.05,0.05}, LEVEL, 5).nf); //if not adapting anymore, return zero
#else
    while (0);
#endif
    //fprintf(stderr, "break3!");
  }
}

/**
   ## Outputs

   We are interested in the viscous dissipation rate. */

/**
   ## Outputs

   We are interested in the viscous dissipation rate in both water and air. */


int dissipation_rate (double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			     sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			     sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
    rateWater += mu1/rho[]*f[]*sqterm; //water
    rateAir   += mu2/rho[]*(1. - f[])*sqterm; //air
  }
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}


/**
   We log the evolution of the kinetic and potential energies and
   dissipation rate as functions of the non-dimensional time. */

event graphs (i++) {
  static FILE * fpwater = fopen("budgetWaterwind.dat", "a");
  static FILE * fpair = fopen("budgetAirwind.dat", "a");
  double ke = 0., gpe = 0.;
  double keAir = 0., gpeAir = 0.;
  foreach(reduction(+:ke) reduction(+:gpe) 
	  reduction(+:keAir) reduction(+:gpeAir)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    ke += rho[]*norm2*f[]*dv();
    keAir += rho[]*norm2*(1.0-f[])*dv();
    gpe += rho1*g_*y*f[]*dv();
    gpeAir += rho2*g_*y*(1.0-f[])*dv();
  }
  double rates[2];
  dissipation_rate(rates);
  double dissWater = rates[0];
  double dissAir   = rates[1];
  if (i == 0) {
    fprintf (fpwater, "t ke gpe dissipation\n");
    fprintf (fpair, "t ke gpe dissipation\n");
  }
  fprintf (fpwater, "%g %g %g %g\n",
	   t/(k_/sqrt(g_*k_)), ke/2., gpe + 0.125, dissWater);
  fprintf (fpair, "%g %g %g %g\n",
	   t/(k_/sqrt(g_*k_)), keAir/2., gpeAir + 0.125, dissAir);
  fprintf (ferr, "%g %g %g %g\n",
	   t/(k_/sqrt(g_*k_)), ke/2., gpe + 0.125, dissWater);
}

/**
   Visualisation

   We use Basilisk view (and output_ppm()) to display animations of the
   results.

   On some parallel systems, pipes tend to cause problems, so we switch
   to simple uncompressed PPM outputs when running with MPI. Otherwise,
   we use MPEG-4 file compression (which requires working pipes and
   ffmpeg). */

#  define POPEN(name, mode) fopen (name ".ppm", mode)

event movies (t += 0.1) {

  /**
     We first do simple movies of the volume fraction, level of
     refinement fields. In 3D, these are in a $z=0$ cross-section. */

  {
    static FILE * fp = POPEN ("f", "a");
    output_ppm (f, fp, min = 0, max = 1, n = 512);
  }

#if TREE
  {
    scalar l[];
    foreach()
      l[] = level;
    static FILE * fp = POPEN ("level", "a");
    output_ppm (l, fp, min = 5, max = LEVEL, n = 512);
  }
#endif

  /**
     <p><center>
     <video width="512" height="512" controls>
     <source src="wave/level.mp4" type="video/mp4">
     Your browser does not support the video tag.
     </video><br>
     Wave breaking. Animation of the level of refinement.
     </center></p>

     We use Basilisk view differently in 2D and 3D. */
  
  scalar omega[];
  vorticity (u, omega);

#if dimension == 2
  view (width = 800, height = 600, fov = 18.8);
  clear();

  /**
     We repeat the drawing periodically in the x-direction. */
  
  for (double x = -L0; x <= L0; x += L0)
    translate (x) {
      draw_vof ("f");
      squares ("omega", linear = true);
    }

  /**
     This gives the following movie.
     <p><center>
     <video width="800" height="600" controls>
     <source src="wave/movie.mp4" type="video/mp4">
     Your browser does not support the video tag.
     </video></center></p>
  */
#else // dimension == 3
  /**
     In 3D, we generate a first movie seen from below. */
  
  view (width = 1600, height = 1200, theta = pi/4, phi = -pi/6, fov = 20);
  clear();
  for (double x = -2*L0; x <= L0; x += L0)
    translate (x) {
      squares ("omega", linear = true, n = {0,0,1}, alpha = -L0/2);
      for (double z = -3*L0; z <= L0; z += L0)
	translate (z = z)
	  draw_vof ("f");
    }
  {
    static FILE * fp = POPEN ("below", "a");
    save (fp = fp);
  }

  /**
     And a second movie, seen from above. */
  
  view (width = 1600, height = 1200, theta = pi/4, phi = pi/6, fov = 20);
  clear();

  /**
     In 3D, we are doubly-periodic (along x and z). */
  
  for (double x = -2*L0; x <= L0; x += L0)
    translate (x) {
      squares ("omega", linear = true, n = {0,0,1}, alpha = -L0/2);
      for (double z = -3*L0; z <= L0; z += L0)
	translate (z = z)
	  draw_vof ("f");
    }
#endif 
  {
    static FILE * fp = POPEN ("movie", "a");
    save (fp = fp);
  }
}

/**
   ## Dump/restore

   To be able to restart, we dump the entire simulation at regular
   intervals. */

event snapshot (i += 200) {
  dump ("dump");
}

/**
   ## End 

   The wave period is `k_/sqrt(g_*k_)`. We want to run up to 2
   (alternatively 4) periods. */

event end (t = 40.*k_/sqrt(g_*k_)) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}


/** 
    ## Dump every 1/32 period.  */

event dumpstep (t += k_/sqrt(g_*k_)/32) {
  char dname[100];
  sprintf (dname, "dump%g", t/(k_/sqrt(g_*k_)));
  dump (dname);
}

/**
   ## Mesh adaptation

   On trees, we adapt the mesh according to the error on volume fraction
   and velocity. */

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, LEVEL, 5);
}
#endif

/**
   ## Running in parallel

   This file will work in 2D or 3D, either with parallel multigrid
   (without adaptivity), using for example:

   ~~~bash
   qcc -source -D_MPI=1 -grid=multigrid3D wave.c
   scp _wave.c occigen.cines.fr:
   ~~~

   and then following a recipe similar to that of the
   [atomisation](/src/examples/atomisation.c#on-occigen) example.

   To use adaptivity, just do something like:

   ~~~bash
   qcc -source -D_MPI=1 -grid=octree wave.c
   scp _wave.c occigen.cines.fr:
   ~~~
*/
