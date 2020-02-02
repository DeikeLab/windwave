 /**
   Run from dump with one or two level more refined grid. Direct output for more accurate interpolation.
   Looping over dump files is achieved with a python script.
   The output files are field_direct+time, field_interp+time, eta+time  */

#include "grid/multigrid.h"
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


/* Direct output on each grid point */
void outfield (int index, double time) {
  
  char fieldname[100], etaname[100], pressurename[100];
  scalar pos[], p_air[], omega[], omega_air[], omega_water[];
  vector tau[];
  tensor SDeform = new tensor;

#if dimension == 2 
  coord G = {0.,1.}, Z = {0.,0.};
#else
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
#endif
  position (f, pos, G, Z);    
  vorticity (u, omega);
  foreach()
  {
    omega_air[] = omega[]*(1-f[]);
    omega_water[] = omega_water[]*f[];
    p_air[] = p[]*(1-f[]);
  }

  // output_ppm (omega_air, // min = -1e-2, max = 1e-2, 
  //       n = 512, linear = true, file = "vort.mp4");   
  sprintf (pressurename, "./diagnostics/pressure_matrix%d.dat", index);
  FILE * fpressure = fopen (pressurename, "w");  
  output_matrix_mpi (p_air, fpressure, n = 512);
  sprintf (fieldname, "./diagnostics/field_direct%d_%d.dat", index, pid());
  FILE * ffield = fopen (fieldname, "w");
  // fprintf(ffield, "x,y,u.x,u.y,f,p,p_air,tau.x,tau.y, omega, omega_air\n");
  foreach()
  {
    p_air[] = p[]*(1-f[]);
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    // double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    // double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    // double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    // double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    // double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    SDeform.x.x[] = dudx;
    SDeform.x.y[] = 0.5*(dudy + dvdx);
    // double SDeformxz = 0.5*(dudz + dwdx);
    SDeform.y.x[] = SDeform.x.y[];
    SDeform.y.y[] = dvdy;
    double mu_eff = mu1/rho[]*f[] + mu2/rho[]*(1. - f[]); // compute effective viscosity
    tau.x[] = 2*mu_eff*(SDeform.x.x[]*0 + SDeform.y.x[]*1);
    tau.y[] = 2*mu_eff*(SDeform.x.y[]*0 + SDeform.y.y[]*1);
    fprintf (ffield, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", 
      x, y, u.x[], u.y[], f[], p[], p_air[], tau.x[], tau.y[], omega[], omega_air[]);
  }
  fclose (ffield);

  // direct output of the interface (x, eta)
  // also direct output of pressure, tau in cartesian coord, normal vector at the interface 
  sprintf (etaname, "./diagnostics/eta%d_%d.dat", index, pid());
  FILE * feta = fopen (etaname, "w");

  // fprintf(feta, "x,pos,f,p,p_p1,p_p2,p_m1,p_m2,tau_x,tau_y,u_x,u_y,n_x,n_y\n");
  // printing out quantities: p_p1 for p at plus 1, p_m1 for p at minus 1 etc.
  foreach(){
    if (interfacial (point, f)){
      // Getting the local normal vector
      coord n = mycs (point, f);
      double norm_2 = sqrt(sq(n.x) + sq(n.y));
      // coord n2 = interface_normal (point, f);  
      // #define interface_normal(point, c) interface_normal (point, c) in src/contact.h
      // n is norm 1 and has to be normalized to norm 2
      double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
      double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
      // double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
      double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
      double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
      // double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
      // double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
      // double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
      // double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
      SDeform.x.x[] = dudx;
      SDeform.x.y[] = 0.5*(dudy + dvdx);
      // double SDeformxz = 0.5*(dudz + dwdx);
      SDeform.y.x[] = SDeform.x.y[];
      SDeform.y.y[] = dvdy;
      double mu_eff = mu2;  // compute effective viscosity
      // double mu_eff = mu1/rho[]*f[] + mu2/rho[]*(1. - f[]); 
      tau.x[] = 2*mu_eff*(SDeform.x.x[]*n.x + SDeform.y.x[]*n.y)/norm_2;
      tau.y[] = 2*mu_eff*(SDeform.x.y[]*n.x + SDeform.y.y[]*n.y)/norm_2;
      // double SDeformyz = 0.5*(dvdz + dwdy);
      // double SDeformzx = SDeformxz;
      // double SDeformzy = SDeformyz;
      // double SDeformzz = dwdz; 
      fprintf (feta, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", 
        x, pos[], f[], p[], p[0,1], p[0,2], p[0,-1], p[0,-2], 
        tau.x[], tau.y[], u.x[], u.y[], n.x/norm_2, n.y/norm_2);
    }
  }
  fclose (feta);
  // return 0;

  dump("dump", list = all);
}

/* Some info about the snapshot time and the refining criteria */
int NUMBER = 32; // Time gap in taking the snapshots
double TIME = 2; // Total number of run time
double snapshot_time; // Current snapshot time
int snapshot_number;
double uemax = 0.0000001;
double femax = 0.0000001;
int LEVEL = 12;
int counting = 0;

/* All the constant that comes with the case, will be changed after input */
double ak = 0.05;
double BO = 200.;
double RE = 40000.;
// double PRESSURE = 0;
double m = 5.;  // vary between 5 and 8
double B = 0.;
double Karman = 0.41;   // Karman universal turbulence constant
double UstarRATIO = 1; // Ratio between Ustar and c

double k_ = 2.*pi;
double h_ = 0.5;
double g_ = 1.;

double RATIO = (1./850.);
double MURATIO = (17.4e-6/8.9e-4); // dynamic viscosity: air/water



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
  if (argc > 8){
    snapshot_number = atoi(argv[8]);
    snapshot_time = snapshot_number*1./32.;
  }
  
  // if (argc > 8)
  //   PRESSURE = atof(argv[8]);
  
  origin (-L0/2, -L0/2, -L0/2);
  periodic (right);
  // We temporally give up the pressure idea
  // p[left] = dirichlet(PRESSURE);
  // p[right] = dirichlet(0); 
  // Notice that the following should be changed according to the specific condition in the case
  u.n[top] = dirichlet(0);
  u.t[top] = neumann(0);

  rho1 = 1.; // 1 stands for water and 2 for air
  rho2 = RATIO;
  mu1 = 1.0/RE; //using wavelength as length scale
  mu2 = 1.0/RE*MURATIO;
  f.sigma = 1./(BO*sq(k_));
  G.y = -g_;

  run();
  return 0;
}

int phony = 1;  // initialize phony as 1, just a trick to play around the way i is treated in basilisk
int j = 0;  // j substitute i as the indexing number

event init (i = 0)
{
  char targetname[100], imagename[100];
  sprintf (targetname, "dump%g", snapshot_time);
  if (!restore (targetname)) 
  	fprintf(ferr, "Not restored!\n");
  restore (targetname);
  phony = 0;
}

event trick (i++) {
  if (phony == 0) {
    j = i;
    phony = 1;
    fprintf(ferr, "New indexing set to j = %d!\n", j);
  }
}

event end (t = 1000) {
	fprintf(ferr, "end\n");
  // return 1;
}

/*event outfield_end (t=end) {
  fprintf(ferr, "outfield_end running!\n");
  outfield(snapshot_time);
  }*/

event outfield_mid (i++) {
if(i == j+2) {
  /* char imagename[100]; */
  fprintf(ferr, "outfield_mid running!\n");
  outfield(snapshot_number, snapshot_time);
  /* scalar omega[], omega_air[]; */
  /* vorticity (u, omega); */
  /* foreach() */
  /* { */
  /*   omega_air[] = omega[]*(1-f[]); */
  /* } */
  /* squares("omega_air", min = -150, max = 150); */
  /* draw_vof("f"); */
  /* sprintf (imagename, "omega%d", snapshot_number); */
  /* { */
  /*   static FILE * fp = fopen (imagename, "w"); */
  /*   save (fp = fp); */
  /* } */
  /* sprintf (imagename, "ux%d", snapshot_number); */
  /* { */
  /*   static FILE * fp = fopen (imagename, "w"); */
  /*   output_ppm (u.x, fp, n = 512); */
  /* } */
  return 1;
}
}

event stop (t = snapshot_time + 200) {
	fprintf(ferr, "end\n");
  return 1;
}

/*event outfield_end (t=end) {
  fprintf(ferr, "outfield_end running!\n");
  outfield(snapshot_time);
  }*/



