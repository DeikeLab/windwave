/**
   Main file for processing the dump files */


#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "tag.h"
#include "curvature.h"
#include "vof.h"

int outfiled (int i);
int direct (double X, int i, char *name);
int shear (int i);

int direct (double X, int i, char *name) {
  int maxlevel = 10; //define the maximum refinement level (used in refine funcetion)
  scalar posx[];
  scalar posy[];
  scalar kap[];
  scalar tau_xy[]; // shear stress
  double RE = 40000;
  double MURATIO = 17.4e-6/8.9e-4; // don't know if dump keeps global parameters like mu so have to hard code
  double MU1 = 1.0/RE;
  double MU2 = MU1*MURATIO;
  char filename[100];
  position (f, posy, {0,1});
  position (f, posx, {1,0});
  curvature (f, kap);
  double Y = 0;
  sprintf (filename, "direct_t%d_%s.dat", i, name);
  FILE *fp = fopen (filename, "w");

  // find the y coordinate of the surface at xpos
  // compute the shear stress in a mu averaged way
  foreach () {
    if ((x > X - Delta) && (x < X + Delta) && (posy[] != nodata))
      Y = posy[];
    double dudy = (u.x[0,1] - u.x[0,-1])/(2*Delta);
    tau_xy[] = (MU1*f[] + MU2*(1-f[]))*dudy;
  }

  // refine the column
  // refine ((x > (X - Delta)) && (x < (X + Delta)) && (level < maxlevel));
  // output velocity above the surface
  foreach () {
    if ((x > (X - Delta)) && (x < (X + Delta)) && (y > Y))
      fprintf (fp, "%g %g %g %g\n", x, y, u.x[], u.y[]);
  }
  fclose (fp);

  // the shear gets written multiple times but it's a temporal solution
  FILE *fp_shear = fopen ("shear.dat", "w");
  foreach () {
    if (posy[] != nodata)
      fprintf (fp_shear, "%g %g %g\n", x, y, tau_xy[]);
  }
  fclose (fp_shear);
  return 0;
}

int outfield (int i) {
  
  char fieldname[100], etaname[100];
  scalar pos[];
#if dimension == 2 
  coord G = {0.,1.}, Z = {0.,0.};
#else
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
#endif
  position (f, pos, G, Z);
  scalar m[];  //define the mask for plotting
    
  /**
     foreach() {
     if (pos[]!=nodata)
     fprintf (fpeta, "%g,%g\n", x, pos[]);
     }
  */

  sprintf (fieldname, "fieldnl%d", i);
  FILE * ffield = fopen (fieldname, "w");
  output_field ({f, u.x, pos}, ffield, n = 128); // linear = true
  fclose (ffield);
  sprintf (etaname, "eta%d", i);
  FILE * feta = fopen (etaname, "w");
  foreach(){
    if (interfacial (point, f)) 
      fprintf (feta, "%g,%g\n", x, pos[]);
  }
  fclose (feta);
  return 0;
}

int main ()
{
  for (int i = 0; i<16; i++) {
    char targetname[100];
    sprintf (targetname, "dump%d", i);
    if (restore(targetname)){
      restore (targetname); //This brings back all the scalar fields f, u, etc
      // Find the X position with the largest eta
      scalar posy[];
      position (f, posy, {0,1});
      double Xmax = 0; // Xmax is the position with the largest eta
      double etamax = -100;
      double Xmin = 0; // Xmin is the position with the minimum eta
      double etamin = 100;
      double Xmid1 = 0;
      double Xmid2 = 0;
      // search for the minimum and maximum position
      foreach(){
	if (etamax < posy[]){
	  etamax = posy[];
	  Xmax = x;
	}
	if (etamin > posy[]){
	  etamin = posy[];
	  Xmin = x;
	}   
      }
      Xmid1 = 0.5*(Xmax + Xmin);
      //getting the velocity profile at X
      direct (Xmax, i, "max");
      direct (Xmin, i, "min");
      direct (Xmid1, i, "mid1");
      //getting the amplitude evolution
      outfield (i);
    }
  }    
  return 0;
}

