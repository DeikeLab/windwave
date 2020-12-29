#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "view.h"

face vector muv[];
face vector av[];
double amp_force = 0.1;
double Ustar = 0.1;
double delta = 1.;

int main() {
  // If we want to do slender channel we still need embedded
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
  run ();
}

event properties (i++)
{
  foreach_face()
    muv.x[] = delta*Ustar/180.;  // Use channel half height and u_tao 
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
    double ytau = delta/180.;
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

#  define POPEN(name, mode) fopen (name ".ppm", mode)
event movies (t += 1) {

  /**
     Movie generation. */

  scalar omega[];
  vorticity (u, omega);
  clear();
  view (fov = 40, camera = "iso", 
  width = 600, height = 600, bg = {1,1,1}, samples = 4);
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -0.9, max = 4, min = -4);
  squares ("u.x", linear = true, n = {0,1,0}, alpha = -0.9, max = 4, min = -4);
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
  if (i <= 20)
    adapt_wavelet ({u}, (double[]){3e-2,3e-2,3e-2}, 9, 5); 
  if (i > 20)
    adapt_wavelet ({u}, (double[]){3e-2,3e-2,3e-2}, 10, 5);
}

event dumpstep (t += 10.) {
  char dname[100];
  p.nodump = false;
  sprintf (dname, "dump%g", t);
  dump (dname);
}

event end (t = 500.) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}
