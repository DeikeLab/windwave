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
/* double X0, Z0, Y0; */

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

void output_slice ()
{
  char filename[100];
  int Nslice = 256;
  /* double zslice = -L0/2; */
  /* double L0 = 2*pi; */
  /* for (int i=0; i<Nslice; i++) { */
  /*   zslice += L0/Nslice; */
  /*   sprintf (filename, "./field/ux_t%g_slice%d", snapshot_time, i); */
  /*   sliceXY (filename,u.x,zslice,MAXLEVEL); */
  /*   sprintf (filename, "./field/uy_t%g_slice%d", snapshot_time, i); */
  /*   sliceXY (filename,u.y,zslice,MAXLEVEL); */
  /*   sprintf (filename, "./field/f_t%g_slice%d", snapshot_time, i); */
  /*   sliceXY (filename,f,zslice,MAXLEVEL); */
  /* } */ 
  double L0 = 4*pi;
  double yslice = 2*pi/Nslice;
  for (int i=0; i<25; i++) {
    yslice += L0/Nslice;
    sprintf (filename, "./field/ux_t%g_slice%d", snapshot_time, i);
    sliceXZ (filename,u.x,yslice,MAXLEVEL);
    sprintf (filename, "./field/uy_t%g_slice%d", snapshot_time, i);
    sliceXZ (filename,u.y,yslice,MAXLEVEL);
    sprintf (filename, "./field/uz_t%g_slice%d", snapshot_time, i);
    sliceXZ (filename,u.z,yslice,MAXLEVEL);
    sprintf (filename, "./field/ux_t%g_slice%d", snapshot_time, -i);
    sliceXZ (filename,u.x,-yslice,MAXLEVEL);
    sprintf (filename, "./field/uy_t%g_slice%d", snapshot_time, -i);
    sliceXZ (filename,u.y,-yslice,MAXLEVEL);
    sprintf (filename, "./field/uz_t%g_slice%d", snapshot_time, -i);
    sliceXZ (filename,u.z,-yslice,MAXLEVEL);
  }  
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
  L0 = 4.*pi;
  /* X0 = 0.; */
  /* Z0 = 0.; */
  /* Y0 = -2.*pi; */
  restore (targetname);
  output_slice();
  run();
}

/* event run_output (i=0) */
/* { */
/* 	/\* refine (level < 8); *\/ */
/* 	output (); */
/* 	cells (n = {0, 0, 1}, alpha = -pi); */
/* 	squares ("u.x", n = {0, 0, 1}, alpha = -pi); */
/* 	save ("ux.ppm"); */
/* 	return 1; */
/* } */
