#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "lambda2.h"
#include "save_data.h"

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

/** Not used direct otput of fields. */
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

/** Output interface info. Deprecated because of added runtime output. */
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

/** Output dissipation in one text file. May not be in the time order. */
void dissipation () {
  /* coord ubar; */
  /* foreach_dimension() { */
  /*   stats s = statsf(u.x); */
  /*   ubar.x = s.sum/s.volume; */
  /* } */
  char dissname[100] = "dissipation.dat";
  FILE * fp = fopen (dissname, "a");
  if (snapshot_time == 0) {
    fprintf(fp, "t,diss\n");
  }
  double ke = 0., vd = 0., vol = 0.;
  foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
    vol += dv()*f[]; // Only count the water phase
    foreach_dimension() {
      // mean fluctuating kinetic energy
      // ke += dv()*sq(u.x[] - ubar.x);
      // viscous dissipation
      vd += sqrt(f[])*f[]*dv()*(sq(u.x[1] - u.x[-1]) +
		      sq(u.x[0,1] - u.x[0,-1]))/sq(2.*Delta);
    }
  }
  // ke /= 2.*vol;
  vd *= mu1;

  /* if (i == 0) */
  /*   fprintf (stderr, "t dissipation energy Reynolds\n"); */
  fprintf (fp, "%g,%g \n", snapshot_time, vd);
}

/**
   Output the underwater velocity as supplementary info for runtime eta file. */

void output_twophase_locate () {
  scalar pos[];
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
  position (f, pos, G, Z);
  char etaname[100];
  sprintf (etaname, "./eta/eta_supplement_t%g", snapshot_time);
  FILE * feta = fopen (etaname, "w");
  fprintf(feta, "x,z,uxw,uyw\n");
  // printing out quantities: p_p1 for p at plus 1, p_m1 for p at minus 1 etc.
  double stp = L0/256.;
  foreach(){
    if (interfacial (point, f)){
      if (point.level == 10) {
	fprintf (feta, "%g,%g,", x, z);
	/** Find the water side velocity. */
 	double yp = y - stp;
	point = locate (x, yp);
	if (point.level > 0) {
	  fprintf (feta, "%g,%g\n", u.x[], u.y[]);
	}
      }
    }
  }
  fflush (feta);
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
  output_twophase_locate ();
  output_slice ();
  char name[80];
  save_data ((scalar*){f,p}, (vector*){u}, i, t, "pvta");
}
