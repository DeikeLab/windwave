#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "lambda2.h"

#define POPEN(name, mode) fopen (name ".ppm", mode)

double snapshot_time = 0;
int MAXLEVEL=7;

void movies () {
  scalar omega[];
  vorticity (u, omega);
  clear();
  view (fov = 32, camera = "iso", ty = -0.25,
  width = 600, height = 600, bg = {1,1,1}, samples = 4);
  squares ("u.x", linear = false, n = {0,0,1}, alpha = -3.1415, max = 8, min = -2, map = cool_warm);
  squares ("u.x", linear = false, n = {1,0,0}, alpha = -3.1415, max = 8, min = -2, map = cool_warm);
  cells (n = {1,0,0}, alpha = -3.1415);
  draw_vof ("f", color = "u.x");
  char s[80];
  sprintf (s, "t = %0.1f", t);
  draw_string (s, size = 30);
  // scalar l2[];
  // lambda2 (u, l2);
  // isosurface ("l2", -1);
  {
    FILE * fp = fopen ("plot.ppm", "a");
    save (fp = fp);
  }
}

int main (int argc, char * argv[] )
{
  double START=0, DT=0.1, FRAMES=1;
  if (argc > 1)
    START = atof (argv[1]);
  if (argc > 2)
    DT = atof (argv[2]);
  if (argc > 3)
    FRAMES = atof (argv[3]);
  char targetname[100];
  L0 = 2.*pi;
  int i;
  for (i=0;i<FRAMES;i++) {
    sprintf (targetname, "dump%g", START+i*DT);
    if (!restore (targetname)) {
      fprintf(ferr, "Not restored!\n");
      return 1;
    }
    restore (targetname);
    movies();
  }
}
