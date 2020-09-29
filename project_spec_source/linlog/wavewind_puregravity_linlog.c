/** 
    This is the wave wind interaction simulation with linear wind profile and adaptive grid.
    Pressure is output on the run by output_matrix_mpi function and also dumped into dump file.
    For velocity initialization, velocity at the top is initialized with slope*height and kept 
    the same across x position. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "tag.h"
#include "navier-stokes/perfs.h"
#include "adapt_wavelet_limited.h"

bool limitedAdaptation = 1;

/**
    Output pressure matrix on the run. */

void output_matrix_mpi (struct OutputMatrix p)
{
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  float fn = p.n, Delta = L0/fn;
  float ** field = matrix_new (p.n, p.n, sizeof(float));

  for (int i = 0; i < p.n; i++) {
    float xp = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2.;
      if (p.linear) {
        field[i][j] = interpolate (p.f, xp, yp);
      }
      else {
        Point point = locate (xp, yp);
        field[i][j] = point.level >= 0 ? val(p.f) : nodata;
      }
    }
  }

  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif


    fwrite (&fn, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2.;
      fwrite (&yp, sizeof(float), 1, p.fp);
    }

    for (int i = 0; i < p.n; i++){
      float xp = Delta*i + X0 + Delta/2.;
      fwrite (&xp, sizeof(float), 1, p.fp);
      for (int j = 0; j < p.n; j++) {
        fwrite (&field[i][j], sizeof(float), 1, p.fp);
      }
    }
    fflush (p.fp);
  }
@if _MPI
  else // slave
  MPI_Reduce (field[0], NULL, p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif

  matrix_free (field);
}


/**
   We log some profiling information. */

// #include "navier-stokes/perfs.h"
// #include "profiling.h"

/**
   The primary parameters are the wave steepness $ak$, the Bond and
   Reynolds numbers. */

double ak = 0.05;
double BO = 200.;
double REw = 40000., REa = 0;

/**
   The default maximum level of refinement depends on the dimension. */

int LEVEL = dimension == 2 ? 10 : 6;

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
   The wave number, fluid depth and acceleration of gravity are set to
   these values. */
double k_ = 2.*pi;
double h_ = 0.5;
double g_ = 9.8;
double T0;
double c_;
int Nwave = 1;
double gpe_base = 0;

/** define the wind profile related parameters. */

double m = 5.;  // vary between 5 and 8
double B = 0.;
double Karman = 0.41;   // Karman universal turbulence constant
double UstarRATIO = 1;   // Ratio between Ustar and c
double Ustar;
double Utop;
double y_1;
double Udrift;   
double TEND = 10;

scalar omega[];
scalar uwater[];
scalar pair[], pairdiff[];


double REGION = 0.2;
int refRegion(double x,double y, double z){
int lev;
if( y < REGION*L0 )
   lev = LEVEL;
 else
   lev = LEVEL-4;
return lev;
}

/**
    We need to store the variable forcing term. */
double amp_force = 0.1; //amplitude of the forcing
face vector av[];

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
    UstarRATIO = atof(argv[3]);
  if (argc > 4)
    Nwave = atoi(argv[4]);
  if (argc > 5)
    L0 = atof(argv[5]);
  if (argc > 6)
    TEND = atof(argv[6]);
  if (argc > 7) 
    REGION = atof(argv[7]);

  /**
     Here we set the densities and viscosities corresponding to the
     parameters above. Note that these variables are defined in two-phase.h already.*/
  h_ = L0/2.;
  k_ = 2.*pi*(Nwave/L0);  
  rho1 = 1000.;
  rho2 = rho1*RATIO;
  mu1 = 8.9e-4;
  mu2 = mu1*MURATIO;
  f.sigma = 0;
  BO = (rho1-rho2)*g_/sq(k_)/f.sigma; // Defined with k instead of lambda
  c_ = sqrt(g_/k_+f.sigma/rho1*k_);
  T0 = 2.*pi/sqrt(g_*k_+f.sigma/rho1*sq(k_)*k_);
  REw = rho1*c_*(2.*pi/k_)/mu1;
  Ustar = c_*UstarRATIO;
  Utop = sq(Ustar)/(mu2/rho2)*L0/2.;
  y_1 = m*mu2/rho2/Ustar;
  Udrift = B*Ustar;
  amp_force = sq(Ustar)/(L0/2.);
#if dimension == 2
  gpe_base = -0.5*sq(h_)*L0*g_;
#else
  gpe_base = -0.5*sq(h_)*sq(L0)*g_;
#endif
  REa = rho2*Ustar*h_/mu2;
  fprintf(stderr, "k = %g, BO = inf, REw = %g, REa = %g \n", k_, REw, REa);
  fprintf(stderr, "a = %g \n", amp_force);
  G.y = -g_;  
  a = av;  
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
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3 - y;
}

double eta (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3;
}

/**
   We either restart (if a "restart" file exists), or initialise the wave
   using the third-order Stokes wave solution. */
event init (i = 0)
{
  fprintf(stderr, "UstarRATIO=%g B=%g\n m=%g ak=%g\n", UstarRATIO, B, m, ak);

  if (!restore ("restart")) {
    do {
      fraction (f, wave(x,y));

      /**
	 To initialise the velocity field, we first define the potential. */
      
      scalar Phi[];
      foreach() {
      	double alpa = 1./tanh(k_*h_);
      	double a_ = ak/k_;
      	double sgma = sqrt(g_*k_*tanh(k_*h_)*
      			   (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
      					      (sq(alpa) - 1.) + sq(alpa))));
      	double A_ = a_*g_/sgma;
      	double phi1 = A_*cosh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x);
      	double phi2 = 3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
      	  cosh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_);
      	double phi3 = 1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
      	  (9.*sq(alpa) - 13.)*cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*k_*x);
      	Phi[] = phi1 + ak*phi2 + ak*ak*phi3;
      } 
      boundary ({Phi});
      foreach(){
	foreach_dimension()
	  u.x[] = (Phi[1] - Phi[-1])/(2.0*Delta) * f[]; // f[] is not strictly 0 or 1 I suppose
      }
      /**
	 Superimpose the air velocity on top. */
      /* double h_local = 0.5; */
      /* foreach(){ */
      /* 	h_local = L0/2. - eta(x,0); // This or just 0.5? */
      /* 	u.x[] += (1-f[])*amp_force*sq(h_local)/(2.*(mu2/rho2))*(1-sq((y-L0/2.)/h_local)); */
      /* } */
      foreach(){
	if ((y-eta(x,y))<y_1){
	  u.x[] += Udrift + sq(Ustar)/(mu2/rho2)* (y-eta(x,y)) * (1-f[]);
	}
	else{
	  double beta = 2*Karman*Ustar/mu2*rho2*((y-eta(x,y))-y_1);
	  double alpha = log(beta+sqrt(sq(beta)+1));
	  double tanhtemp = (exp(alpha/2)-exp(-alpha/2))/(exp(alpha/2)+exp(-alpha/2));
	  u.x[] += (Udrift + m*Ustar + Ustar/Karman*(alpha-tanhtemp)) * (1-f[]);
	}
      }
      boundary ((scalar *){u});  // type casting   
    }
    /**
       On trees, we repeat this initialisation until mesh adaptation does
       not refine the mesh anymore. */

#if TREE  
    while (adapt_wavelet ({f, u},
			  (double[]){0.01,0.05,0.05,0.05}, LEVEL, 5).nf); //if not adapting anymore, return zero
#else
    while (0);
#endif
  }
}

event acceleration (i++) {
  /**
     Forcing term equivalent to pressure gradient in x. */
  foreach_face(x)
    av.x[] += amp_force*(1-f[]);
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

event graphs (i+=10) {
  static FILE * fpwater = fopen("budgetWaterwind.dat", "a");
  static FILE * fpair = fopen("budgetAirwind.dat", "a");
  double ke = 0., gpe = 0.;
  double keAir = 0., gpeAir = 0.;
  foreach(reduction(+:ke) reduction(+:gpe) 
	  reduction(+:keAir) reduction(+:gpeAir)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    /* ke += rho[]*norm2*f[]*dv(); */
    /* keAir += rho[]*norm2*(1.0-f[])*dv(); */
    /* gpe += rho[]*g_*y*f[]*dv(); */
    /* gpeAir += rho[]*g_*y*(1.0-f[])*dv(); */
    ke += rho1*norm2*f[]*dv();
    keAir += rho2*norm2*(1.0-f[])*dv();
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
	   t, ke/2., gpe - gpe_base*rho1, dissWater);
  fprintf (fpair, "%g %g %g %g\n",
	   t, keAir/2., gpeAir - gpe_base*rho2, dissAir);
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
     Draw velocity only in waer. */
  char s[80];
  sprintf (s, "t = %0.2f", t);
  vorticity (u, omega);
  view (width = 600, height = 600, fov = 18.8);
  clear();
  draw_string (s, size = 30);
  squares("omega", linear=true);
  draw_vof("f"); 
  {
    static FILE * fp = POPEN ("omega", "a");
    save (fp = fp);
  }
  clear();
  draw_string (s, size = 30);
  squares("u.x", linear=true, max=5, min=-5);
  draw_vof("f");
  {
    static FILE * fp = POPEN ("u_air", "a");
    save (fp = fp);
  }
  clear();
  foreach ()
    uwater[] = u.x[]*f[];
  draw_string (s, size = 30);
  squares("uwater", linear=true);
  draw_vof("f");
  {
    static FILE * fp = POPEN ("u_water", "a");
    save (fp = fp);
  }
  foreach()
    pair[] = p[]*(1-f[]);
  norm n1 = normf(pair);
  // norm n2 = normf(f);
  foreach()
    pairdiff[] = (pair[] - n1.avg/0.5)*(1-f[]);
  fprintf(stderr, "pair average = %g\n", n1.avg);
  clear();
  draw_string (s, size = 30);
  squares("pairdiff", linear=true);
  draw_vof("f");
  {
    static FILE * fp = POPEN ("p_air", "a");
    save (fp = fp);
  }
}

/**
   ## Dump/restore

   To be able to restart, we dump the entire simulation at regular
   intervals. */

/* event snapshot (i += 200) { */
/*   dump ("dump"); */
/* } */

/**
   ## End 

   The wave period is `k_/sqrt(g_*k_)`. We want to run up to 2
   (alternatively 4) periods. */

event end (t = TEND*T0) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}

event pseudo_initial (i = 10) {
  dump ("initial");
}

/** 
    ## Dump every 1/32 period.  */

event dumpstep (t = 0; t += 0.1) {
  fprintf (stderr, "running dumpstep!\n");
  char dname[100];
  sprintf (dname, "dump%g", t);
  p.nodump = false;
  dump (dname);
  // Output pressure
  char pressurename1[100], pressurename2[100], pressurename3[100];
  sprintf (pressurename1, "./pressure/p_matrix%g.dat", t);
  sprintf (pressurename2, "./pressure/p_air_matrix%g.dat", t);
  sprintf (pressurename3, "./pressure/p_air_diff_matrix%g.dat", t);
  FILE * fpressure1 = fopen (pressurename1, "w");
  output_matrix_mpi (p, fpressure1, n = 512);
  FILE * fpressure2 = fopen (pressurename2, "w");
  output_matrix_mpi (pair, fpressure2, n = 512);
  FILE * fpressure3 = fopen (pressurename3, "w");
  output_matrix_mpi (pairdiff, fpressure3, n = 512);
  char fname[100];
  sprintf (fname, "./pressure/f_matrix%g.dat", t);
  FILE * ff = fopen (fname, "w");
  output_matrix_mpi (f, ff, n = 512);
}

/**
   ## Mesh adaptation

   On trees, we adapt the mesh according to the error on volume fraction
   and velocity. */

#if TREE
event adapt (i++) {
  if(limitedAdaptation)
    adapt_wavelet_limited ({f,u}, (double[]){femax,uemax,uemax,uemax}, refRegion, 5);
  else
    adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, LEVEL, 5);
  //adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, LEVEL, 5);
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
