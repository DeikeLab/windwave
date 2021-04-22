#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "lambda2.h"


double snapshot_time = 0;
int OUTLEVEL=7;

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
	  double zp = stp*j + Y0 + stp/2.;
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
	fprintf (fpver, "%.9g\t", field[i][j]);
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

void output_slice ()
{
  char filename[100];
  int Nslice = 256;
  double L0 = 2*pi;
  double zslice = -L0/2+L0/2./Nslice;
  for (int i=0; i<Nslice; i++) {
    zslice += L0/Nslice;
    sprintf (filename, "./field/ux_t%g_slice%d", snapshot_time, i);
    sliceXY (filename,u.x,zslice,OUTLEVEL);
    sprintf (filename, "./field/uy_t%g_slice%d", snapshot_time, i);
    sliceXY (filename,u.y,zslice,OUTLEVEL);
    sprintf (filename, "./field/f_t%g_slice%d", snapshot_time, i);
    sliceXY (filename,f,zslice,OUTLEVEL);
    /* sprintf (filename, "./field/p_t%g_slice%d", snapshot_time, i); */
    /* sliceXY (filename,p,zslice,OUTLEVEL); */
  } 
  /* double L0 = 2*pi; */
  /* double yslice = pi/Nslice; */
  /* for (int i=0; i<Nslice; i++) { */
  /*   yslice += L0/Nslice; */
  /*   sprintf (filename, "./field/ux_t%g_slice%d", snapshot_time, i); */
  /*   sliceXZ (filename,u.x,yslice,OUTLEVEL); */
  /*   sprintf (filename, "./field/uy_t%g_slice%d", snapshot_time, i); */
  /*   sliceXZ (filename,u.y,yslice,OUTLEVEL); */
  /*   sprintf (filename, "./field/f_t%g_slice%d", snapshot_time, i); */
  /*   sliceXZ (filename,f,yslice,OUTLEVEL); */
  /* }   */
}

void output ()
{
  char filename[100];
  sprintf (filename, "field_%g", snapshot_time);
  FILE * fp = fopen (filename, "w");
  fprintf(fp, "x,y,z,u.x,u.y,u.z,dudy,p,g.x\n");
  // direct output of 5 fields
  // %g in C gives the shortest representation possible
  // The short answer to print floating point numbers losslessly 
  // (such that they can be read back in to exactly the same number, except NaN and Infinity):
  // If your type is float: use printf("%.9g", number).
  // If your type is double: use printf("%.17g", number).
  int count = 0;
  foreach()
    {
      double dudy = (u.x[0,1]-u.x[0,-1])/(2.*Delta);
      fprintf (fp, "%g,%g,%g,%g,%g,%g,%g,%g,%g\n", x, y, z, u.x[], u.y[], u.z[], dudy,p[],g.x[]);
      count += 1;
    }
  fprintf(stderr, "%i\n", count);
  fflush (fp);
  fclose (fp);
}

scalar wateru[];
void output_twophase () {
  foreach () {
    wateru[] = u.x[]*f[];
  }
  squares("wateru", n = {0, 0, 1}, alpha = 0);
  save ("uwater.ppm");
  scalar pos[];
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
  position (f, pos, G, Z);
  char etaname[100];
  sprintf (etaname, "./eta/eta_t%g_%d", snapshot_time, pid());
  FILE * feta = fopen (etaname, "w");
  fprintf(feta, "x,z,pos,p,p_p1,p_p2\n");
  // printing out quantities: p_p1 for p at plus 1, p_m1 for p at minus 1 etc.
  foreach(){
    if (interfacial (point, f)){
      fprintf (feta, "%g,%g,%g,%g,%g,%g\n", 
	       x, z, pos[], p[], p[0,1], p[0,2]);
      /* tau.x[], tau.y[], u.x[], u.y[], n.x/norm_2, n.y/norm_2); */
    }
  }
  fclose (feta);
} 

int main(int argc, char *argv[]) {
  if (argc > 1)
    snapshot_time = atof(argv[1]);
  if (argc > 2)
    OUTLEVEL = atoi(argv[2]);

  run();
}

event init(i=0)
{
  char targetname[100];
  sprintf (targetname, "dump%g", snapshot_time);
  if (!restore (targetname)) {
    fprintf(ferr, "Not restored!\n");
    return 1;
  }
  //output_twophase();
  output_slice();
}
/* event acceleration (i++) { */
/*   double ampl = sq(Ustar)/(L0-h_); */
/*   foreach_face(x) */
/*     av.x[] += ampl*(1.-f[]); */
/* } */

/* event p_air (i=5) */
/* { */
/*   scalar pair[]; */
/*   foreach () { */
/*     pair[] = p[]*(1.-f[]); */
/*   } */
/*   squares ("pair", n = {0,0,1}, alpha = -0.1); */
/*   save ("pair.ppm"); */
/*   cells (n = {0,0,1}, alpha = -0.1); */
/*   squares ("level", n = {0,0,1}, alpha = -0.1, max = 10, min = 5); */
/*   save ("LEVEL.ppm"); */

/*   char filename[100]; */
/*   int Nslice = 256; */
/*   double L0 = 2*pi; */
/*   double zslice = -L0/2+L0/2./Nslice; */
/*   for (int i=0; i<Nslice; i++) { */
/*     zslice += L0/Nslice; */
/*     sprintf (filename, "./field/p_run_t%g_slice%d", snapshot_time, i); */
/*     sliceXY (filename,p,zslice,OUTLEVEL); */
/*     sprintf (filename, "./field/pair_run_t%g_slice%d", snapshot_time, i); */
/*     sliceXY (filename,pair,zslice,OUTLEVEL); */
/*   }    */
/* } */

/* #if TREE */
/* event adapt (i++) { */
/*   if (i == 5) */
/*     fprintf(stderr, "uemaxRATIO = %g\n", uemaxRATIO); */
/*   if (t < RELEASETIME) */
/*     adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL); */
/*   if (t >= RELEASETIME) { */
/*     foreach () { */
/*       foreach_dimension () */
/* 	u_water.x[] = u.x[]*f[]; */
/*     } */
/*     adapt_wavelet ({f,u,u_water}, (double[]){femax,uemax,uemax,uemax,0.001,0.001,0.001}, MAXLEVEL); */
/*   } */
/* } */
/* #endif */
