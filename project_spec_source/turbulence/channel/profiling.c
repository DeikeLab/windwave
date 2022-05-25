#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "./sandbox/profile6.h"

#define POPEN(name, mode) fopen (name ".ppm", mode)
double snapshot_time = 0;
int MAXLEVEL=7;
/* double X0, Z0, Y0; */


void profile_output () {
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
  profiles ({u.x, u.y, u.z, uxuy, uxux, uyuy, uzuz}, phi, rf = 0.5, fname = file, n = 1 << MAXLEVEL, min = -1., max = 1.);
}


int main (int argc, char * argv[] )
{
  if (argc > 1)
    snapshot_time = atof(argv[1]);
  if (argc > 2)
    MAXLEVEL = atoi(argv[2]);
  char targetname[100];
  sprintf (targetname, "dump%g", snapshot_time);
  if (!restore (targetname)) {
    fprintf(ferr, "Not restored!\n");
    return 1;
  }
  L0 = 2.;
  /* X0 = 0.; */
  /* Z0 = 0.; */
  /* Y0 = -2.*pi; */
  restore (targetname);
  profile_output();
  run();
}
