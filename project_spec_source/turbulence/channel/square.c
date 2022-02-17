#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "./sandbox/profile6.h"
#include "navier-stokes/perfs.h"

face vector muv[];
face vector av[];
double amp_force = 0.1;
double Ustar = 0.1;
double delta = 1.;
int MAXLEVEL = dimension == 2 ? 10 : 5; // max level if not use limited refinement
double RE_tau = 180.;
double uemax = 0.01;
double uemaxRATIO = 0.1;

int main(int argc, char *argv[]) {
  // If we want to do slender channel we still need embedded
  if (argc > 1)
    MAXLEVEL = atoi(argv[1]);
  if (argc > 2)
    RE_tau = atof(argv[2]);
  if (argc > 3)
    uemaxRATIO = atof(argv[3]);
  L0 = 2.;
  origin (-L0/2., -L0/2., -L0/2.);
  N = 128;
  mu = muv;
  a = av;
  u.r[top] = dirichlet(0);
  u.r[bottom] = dirichlet(0);
  u.n[top] = dirichlet(0);
  u.n[bottom] = dirichlet(0);
  u.t[top] = dirichlet(0);
  u.t[bottom] = dirichlet(0); 
  periodic (right);
  periodic (front);
  uemax = uemaxRATIO*Ustar;
  fprintf (stderr, "RE_tau = %g, MAXLEVEL = %d, uemax = %g \n", RE_tau, MAXLEVEL, uemax);
  run ();
}

event properties (i++)
{
  foreach_face()
    muv.x[] = delta*Ustar/RE_tau;  // Use channel half height and u_tao 
}

//u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
//u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);

double randInRange(int min, int max)
{
  return min + (rand() / (double) (RAND_MAX) * (max - min + 1));
}

event init (t = 0) {
  /* vertex scalar phi[]; */
  /* foreach_vertex() { */
  /*   phi[] = intersection (delta-y, delta+y);  // Slender channel using embedded */
  /* } */
  /* boundary ({phi}); */
  /* fractions (phi,cs,fs); */
  if (!restore("restart")) {
    double rand = 0;
    double ytau = delta/RE_tau;
    foreach() {
      rand = randInRange (0,0.1);
      // Initialize with a more accurate profile
      if (y < 0.)
	u.x[] = (log((y+L0/2.)/ytau)*Ustar/0.41);
      else
	u.x[] = (log((L0/2.-y)/ytau)*Ustar/0.41);;
      //u.x[] = (2+rand)*(1.-f[]);  
      u.y[] = 0.;
      u.z[] = 0.;
    }
    boundary ((scalar *){u});
  }
}

event acceleration (i++) {
  /**
     Forcing term equivalent to pressure gradient in x. */
  double ampl = sq(Ustar)/delta;
  foreach_face(x)
    av.x[] += ampl;
}

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
  int npoint = 1<<MAXLEVEL;
  profiles ({u.x, u.y, u.z, uxuy, uxux, uyuy, uzuz}, phi, rf = 0.5, fname = file, n = npoint, min = -1., max = 1.);
}

#  define POPEN(name, mode) fopen (name ".ppm", mode)
event movies (t += 1) {

  /**
     Movie generation. */

  scalar omega[];
  vorticity (u, omega);
  clear();
  view (fov = 40, camera = "iso", 
  width = 600, height = 600, bg = {1,1,1}, samples = 4);
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -0.9);
  squares ("u.x", linear = true, n = {0,1,0}, alpha = -0.9);
  squares ("omega", linear = true, n = {1,0,0}, alpha = -0.9);
  cells (n = {1,0,0}, alpha = -0.9);
  // draw_vof ("f", color = "u.x");
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

event adapt (i++) {
    adapt_wavelet ({u}, (double[]){uemax,uemax,uemax}, MAXLEVEL, 5);
}

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
	  double zp = stp*j + Z0 + stp/2.;
	  Point point = locate (xp, yp, zp);
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
	  Point point = locate (xp, yp, zp);
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

/* event output_slice (t +=1) */
/* { */
/*   char filename[100]; */
/*   int Nslice = 256; */
/*   double zslice = -L0/2. + L0/2./Nslice;; */
/*   for (int i=0; i<Nslice; i++) { */
/*     sprintf (filename, "./field/ux_t%g_slice%d", t, i); */
/*     sliceXY (filename,u.x,zslice,9); */
/*     sprintf (filename, "./field/uy_t%g_slice%d", t, i); */
/*     sliceXY (filename,u.y,zslice,9); */
/*     sprintf (filename, "./field/uz_t%g_slice%d", t, i); */
/*     sliceXY (filename,u.z,zslice,9); */
/*     zslice += L0/Nslice; */
/*   }  */
/*   /\* double yslice = -L0/2. + L0/2./Nslice; *\/ */
/*   /\* for (int i=0; i<Nslice; i++) { *\/ */
/*   /\*   sprintf (filename, "./field/ux_t%g_slice%d", t, i); *\/ */
/*   /\*   sliceXZ (filename,u.x,yslice,MAXLEVEL); *\/ */
/*   /\*   sprintf (filename, "./field/uy_t%g_slice%d", t, i); *\/ */
/*   /\*   sliceXZ (filename,u.y,yslice,MAXLEVEL); *\/ */
/*   /\*   sprintf (filename, "./field/uz_t%g_slice%d", t, i); *\/ */
/*   /\*   sliceXZ (filename,u.z,yslice,MAXLEVEL); *\/ */
/*   /\*   yslice += L0/Nslice; *\/ */
/*   /\* }   *\/ */
/* } */

event dumpstep (t += 5) {
  char dname[100];
  p.nodump = false;
  sprintf (dname, "dump%g", t);
  dump (dname);
}

event dump_restart (t += 1) {
  dump ("restart");
  //output_slice ();
}

event end (t = 1000.) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}
