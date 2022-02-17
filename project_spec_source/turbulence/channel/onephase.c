#include "grid/octree.h"
#include "embed.h"
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
int phony = 10000;

int main(int argc, char *argv[]) {
  // If we want to do slender channel we still need embedded
  if (argc > 1)
    MAXLEVEL = atoi(argv[1]);
  if (argc > 2)
    RE_tau = atof(argv[2]);
  L0 = 4*pi;
  origin (0, -L0/2., 0);
  N = 1 << (MAXLEVEL-2);
  mu = muv;
  a = av;
  /* u.r[top] = dirichlet(0); */
  /* u.r[bottom] = dirichlet(0); */
  /* u.n[top] = dirichlet(0); */
  /* u.n[bottom] = dirichlet(0); */
  /* u.t[top] = dirichlet(0); */
  /* u.t[bottom] = dirichlet(0);  */
  /* TOLERANCE = 1e-4; */
  /* NITERMAX = 20; */
  periodic (right);
  periodic (front);
  run ();
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*delta*Ustar/180.;  // Use channel half height and u_tao 
}

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

double randInRange(int min, int max)
{
  return min + (rand() / (double) (RAND_MAX) * (max - min + 1));
}

event init (t = 0) {
  if (!restore("restart")) {
    // restore from level 9 
    phony = 0;
    vertex scalar phi[];
    foreach_vertex() {
      phi[] = intersection (delta-y, delta+y);  // Slender channel using embedded
    }
    boundary ({phi});
    fractions (phi,cs,fs);
    double ytau = delta/180.;
    foreach() {
      if ((y < 0.) && (y > -delta))
	u.x[] = cs[] ? (log((y+delta)/ytau)*Ustar/0.41) : 0.;
      else if ((y > 0.) && (y < delta))
	u.x[] = cs[] ? (log((delta-y)/ytau)*Ustar/0.41) : 0.;
      else 
	u.x[] = 0.;  
      u.y[] = 0.;
      u.z[] = 0.;
    } 
    boundary((scalar *){u});
  }
  else {
    phony = 0;
    for (scalar s in {u, g})
      s.prolongation = refine_injection;
    restore ("restart");
    for (scalar s in {u, g})
      s.prolongation = refine_embed_linear;
    vertex scalar phi[];
    foreach_vertex() {
      phi[] = intersection (delta-y, delta+y);  // Slender channel using embedded
    }
    boundary ({phi});
    fractions (phi,cs,fs);
    boundary (all);
  }
}

event acceleration (i++) {
  /**
     Forcing term equivalent to pressure gradient in x. */
  double ampl = sq(Ustar)/delta;
  foreach_face(x){
      av.x[] += ampl*fm.x[]*cs[];
  }
}

event profile_output (t += 10) {
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
  profiles ({u.x, u.y, u.z, uxuy, uxux, uyuy, uzuz}, phi, rf = 0.5,  fname = file, min = -1.1, max = 1.1);
}

#  define POPEN(name, mode) fopen (name ".ppm", mode)
event movies_diagnosis (i += 1) {
  clear();
  view (fov = 20, camera = "front", tx = -0.5, ty = 0.05);
  squares("u.x", n = {0,0,1}, alpha = 0.1);
  cells(n = {0,0,1}, alpha = 0.1);
  char s[80];
  sprintf (s, "i = %d", i);
  draw_string (s, size = 30);
  {
    static FILE * fp1 = POPEN ("ux", "a");
    save (fp = fp1);
  }
  clear();
  view (fov = 20, camera = "front", tx = -0.5, ty = 0.05);
  squares("cs", n = {0,0,1}, alpha = 0.1);
  cells(n = {0,0,1}, alpha = 0.1);
  sprintf (s, "i = %d", i);
  draw_string (s, size = 30);
  {
    static FILE * fp2 = POPEN ("cs", "a");
    save (fp = fp2);
  }
}

event movies (t += 1) {

  /**
     Movie generation. */

  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  clear();
  view (fov = 40, camera = "iso", ty = 0.5, tx = 0.05,
  width = 600, height = 600, bg = {1,1,1}, samples = 4);
  squares ("u.x", linear = true, n = {0,1,0}, alpha = -0.1);
  squares ("u.x", linear = true, n = {0,0,1}, alpha = 0.1);
  squares ("omega", linear = true, n = {1,0,0}, alpha = 0.1);
  cells (n = {1,0,0}, alpha = 0.1);
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
  /* output_ppm (u.x, file = "2D.ppm", box = {{0,4*pi},{-1,1}}, */
  /* 	      min = -4, max = 4, linear = true, mask = m); */
}

event adapt (i++) {
  if (phony <= 10)
    adapt_wavelet ({cs,u}, (double[]){1e-4,3e-2,3e-2,3e-2}, 9, 6);
  if (phony > 10)
    adapt_wavelet ({cs,u}, (double[]){1e-4,1e-1,1e-1,1e-1}, MAXLEVEL, MAXLEVEL-3);
  phony += 1;
}

event dumpstep (t += 10.) {
  char dname[100];
  //p.nodump = false;
  //phi.nodump = true;
  sprintf (dname, "dump%g", t);
  scalar pid[];
  foreach ()
    pid[] = fmod(pid()*(npe()+37),npe());
  boundary ({pid});
  dump (dname);
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
	  Point point = locate (xp, yp,zp);
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

event turbulence_stat (t += 1) {
  char filename[100];
  int Nslice = 256;
  double zslice = -L0/2+L0/2./Nslice;
  for (int i=0; i<Nslice; i++) {
    sprintf (filename, "./field/ux_t%g_slice%d", t, i);
    sliceXY (filename,u.x,zslice,MAXLEVEL-1);
    sprintf (filename, "./field/uy_t%g_slice%d", t, i);
    sliceXY (filename,u.y,zslice,MAXLEVEL-1);
    sprintf (filename, "./field/uz_t%g_slice%d", t, i);
    sliceXY (filename,u.z,zslice,MAXLEVEL-1);
    zslice += L0/Nslice;
  }
}

event dump_restart (t += 1) {
  dump ("restart");
}

event end (t = 1000.) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}
