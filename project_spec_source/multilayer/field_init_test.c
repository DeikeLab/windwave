/**
# Field scale wave breaking (multilayer solver)

 */

#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
//#include "layered/remap.h"
#include "remap_test.h"
#include "layered/perfs.h"
#include "input.h"
#include "output_mpi.h"

/**
The initial condition is a externally imported wave field. Some controlling parameters. */

#define g_ 9.8
#define h_ 25
/* double ETAE = 0.1; // refinement criteria for eta */
/* int MAXLEVEL = 8; // max level of refinement in adapt_wavelet function */
/* int MINLEVEL = 6; // min level */
double TEND = 50.; // t end
int NLAYER = 10; // number of layers
int LEVEL_data = 7;

double kp_ = 2.*pi/10.;
double P_ = 0.01;
int N_power_ = 5;
#define N_mode_ 16
double F_kxky_[N_mode_*N_mode_], omega[N_mode_*N_mode_], phase[N_mode_*N_mode_];
double kx_[N_mode_], ky_[N_mode_];
double dkx_, dky_;

double randInRange(int min, int max)
{
  return min + (rand() / (double) (RAND_MAX) * (max - min + 1));
}

void power_input () {
	for (int i=0; i<N_mode_; i++) {
		kx_[i] = 2.*pi/L0*(i+1);
		ky_[i] = 2.*pi/L0*(i-N_mode_/2);
	}
	dkx_ = kx_[1]-kx_[0];
	dky_ = ky_[1]-ky_[0];
	// A function that reads in F_kxky_. Next step is to generate F_kxky_ all inside
	double size = N_mode_*N_mode_;
	float * a = (float*) malloc (sizeof(float)*size);
	char filename[100];
	sprintf (filename, "F_kxky");
  FILE * fp = fopen (filename, "rb");
  fread (a, sizeof(float), size, fp);
	for (int i=0;i<size;i++) {
		F_kxky_[i] = (double)a[i];
		//fprintf(stderr, "%g \n", F_kxky_[i]);
	}
	// Phase and omega, next focusing phase
	double kmod = 0;
	int index = 0;
	for (int i=0; i<N_mode_;i++) {
		for (int j=0;j<N_mode_;j++) {
			index = i*N_mode_ + j;
			kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
			omega[index] = sqrt(g_*kmod);
			phase[index] = randInRange (0, 2.*pi);
		}
	}
}
double wave (double x, double y)
{
	double eta = 0;
	double ampl = 0, a = 0;
	int index = 0;
	for (int i=0; i<N_mode_;i++) {
		for (int j=0;j<N_mode_;j++) {
			index = i*N_mode_ + j;
			ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
			a = (kx_[i]*x + ky_[j]*y + phase[index]);
			eta += ampl*cos(a);
		}
	}
	return eta;
}
double u_x (double x, double y, double z) {
	int index = 0;
	double u_x = 0;
	double ampl = 0, a = 0;
	double z_actual = 0, kmod = 0, theta = 0;
	for (int i=0; i<N_mode_;i++) {
		for (int j=0;j<N_mode_;j++) {
			index = i*N_mode_ + j;
			ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
			z_actual = (z < ampl ? (z) : ampl);
			//fprintf(stderr, "z = %g, ampl = %g, z_actual = %g\n", z, ampl, z_actual);
			kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
			theta = atan(ky_[j]/kx_[i]);
			a = (kx_[i]*x + ky_[j]*y + phase[index]);
			u_x += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*cos(a)*cos(theta);
		}
	}
	return u_x;
}
// #if dimension = 2
double u_y (double x, double y, double z) {
	int index = 0;
	double u_y = 0;
	double ampl = 0, a = 0;
	double z_actual = 0, kmod = 0, theta = 0;
	for (int i=0; i<N_mode_;i++) {
		for (int j=0;j<N_mode_;j++) {
			index = i*N_mode_ + j;
			ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
			z_actual = (z < ampl ? (z) : ampl);
			kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
			theta = atan(ky_[j]/kx_[i]);
			a = (kx_[i]*x + ky_[j]*y + phase[index]);
			u_y += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*cos(a)*sin(theta);
		}
	}
	return u_y;
}
// #endif
double u_z (double x, double y, double z) {
	int index = 0;
	double u_z = 0;
	double ampl = 0, a = 0;
	double z_actual = 0, kmod = 0;
	for (int i=0; i<N_mode_;i++) {
		for (int j=0;j<N_mode_;j++) {
			index = i*N_mode_ + j;
			ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
			z_actual = (z < ampl ? (z) : ampl);
			kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
			a = (kx_[i]*x + ky_[j]*y + phase[index]);
			u_z += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*sin(a);
		}
	}
	return u_z;
}


int main(int argc, char * argv[])
{
	if (argc > 1)
		NLAYER = atoi(argv[1]);
	if (argc > 2)
		LEVEL_data = atoi(argv[2]);
	if (argc > 3)
    TEND = atof(argv[3]);
	L0 = 50.;
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  N = 1 << LEVEL_data; // start with a grid of 128
  nl = NLAYER;
  G = g_;
  nu = 1/40000.;
  CFL_H = 1; // Smaller time step
  run();
}

event init (i = 0)
{
	power_input();
  foreach() {
    zb[] = -h_;
    double H = wave(x, y) - zb[];
    foreach_layer() {
      h[] = H/nl;
    }
	}
		// remap?
	vertical_remapping (h, tracers);
	foreach() {
		double z = zb[];
		foreach_layer() {
			z += h[]/2.;
      u.x[] = u_x(x, y, z);
			u.y[] = u_y(x, y, z);
      w[] = u_z(x, y, z);
      z += h[]/2.;
		}
  }
	fprintf (stderr,"Done initialization!\n");
	dump("initial");
}

event energy_before_remap (i++, last)
{
	if (i==10) {
		fprintf(stderr, "energy output before remap!\n");
		fflush(stderr);
	}
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
		double zc = zb[];
		foreach_layer () {
		  double norm2 = sq(w[]);
      foreach_dimension()
			  norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
		}
  }
	static FILE * fp = fopen("energy_before_remap.dat","w");
  fprintf (fp, "%g %g %g\n", t, ke/2., g_*gpe+7656250.);
  fflush (fp);
}

event energy_after_remap (i++, last)
{
	if (i==10) {
		fprintf(stderr, "energy output after remap!\n");
		fflush(stderr);
	}
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
		double zc = zb[];
		foreach_layer () {
		  double norm2 = sq(w[]);
      foreach_dimension()
			  norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
		}
  }
	static FILE * fp = fopen("energy_after_remap.dat","w");
  fprintf (fp, "%g %g %g\n", t, ke/2., g_*gpe+7656250.);
  fflush (fp);
}

/**
Note that the movie generation below is very expensive. */

#if 1
event movie (t += 0.1; t <= TEND)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0) 
    fprintf (fp, "set term pngcairo font ',9' size 1024,768;"
	     "unset key\n"
	     "set pm3d interpolate 4,4 lighting specular 0.6\n"
	     "set zrange [-5:5]\n"
	     "set cbrange [-1:1]\n"
	     "set xlabel 'x'\n"
	     "set ylabel 'y'\n"
	     );
  fprintf (fp,
	   "set output 'plot%05d.png'\n"
	   "set title 't = %.2f'\n"
	   "splot '-' u 1:2:3:4 w pm3d\n",
	   i, t);
	char s[10];
	sprintf (s, "u%d", nl-1);
  vector u_temp = lookup_vector (s);
  output_field ({eta,u_temp.x}, fp, linear = true);
  fprintf (fp, "e\n\n");
  fflush (fp);
        char filename1[50], filename2[50];
	sprintf (filename1, "surface/eta_matrix_%g", t);
	sprintf (filename2, "surface/ux_matrix_%g", t);
	FILE * feta = fopen (filename1, "w");
	// Might need to change to mpi function later
	output_matrix_mpi (eta, feta, N, linear = true);
	fclose (feta);
	FILE * fux = fopen (filename2, "w");
	output_matrix_mpi (u_temp.x, fux, N, linear = true);
	fclose (fux);
}

/* event closeupmovie (t += 0.1; t <= TEND) */
/* { */
/*   view (fov = 17.3106, quat = {0.475152,0.161235,0.235565,0.832313}, */
/* 	tx = -0.0221727, ty = -0.0140227, width = 1200, height = 768); */
/*   char s[80]; */
/*   sprintf (s, "t = %.2f T0", t/T0); */
/*   draw_string (s, size = 80); */
/*   for (double x = -1; x <= 1; x++) */
/*     translate (x) */
/*       squares ("u59.x", linear = true, z = "eta", min = -2., max = 2.); */
/*   save ("movie.mp4"); */
/* } */

event moviemaker (t = end)
{
  system ("for f in plot*.png; do convert $f ppm:-; done | ppm2mp4 movie.mp4");
}
#endif

#if QUADTREE
event adapt (i++) {
	/* fprintf(stderr, "Adapting start!\n"); */
	/* fflush(stderr); */
	my_adapt();
}
#endif


event endrun (t = TEND) {
	dump();
}
