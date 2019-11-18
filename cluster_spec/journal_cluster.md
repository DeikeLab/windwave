x102.11 Trouble-shooting:

On adroit:
1. When running slurm script, run through relative path. Use srun when running parallel, otherwise it causes convergence problem. 
2. Added compiling flag -ldl to link dynamic libraries.
3. Reorganized the Makefiles. 
To-do: 
wave.slurm script fails to run from /scratch/jiarongw/... [FIXED]
Set PATH variables in ./bashrc and add "module load mpi" as well. [DONE]
Explore the possibility to code output directory into the source code. The question is how to change the function of basilisk like save().

On tiger:
1. -I{BASILISK} flag should be placed behind the {SRC_FILE}. It includes the basic basilisk header files. For graphic libraries, it is more important to add the linking path than the including path. {OPENGLINC} name should be changed to avoid confusion. 

---------------
02.12 Trouble-shooting

1. NOTICE: /scratch/netid only exists on the particular head node you are using.
Must use /scratch/network/netid on adroit and /scratch/gpfs/netid on tiger.
2. Better not to use _ when naming the directory. Still haven't figured out a good way to put the strings together.


--------------
04.17

Modified contains the result from the new source code where the vertical kinetic energy is taken into account. Linear contains where the whole profile is set to a linear shape.

-------------
08.13

A repeated problem of not being able to restart, either "empty flag" or "error mapping":
The dump file was not fully written because quota was running out.
Use command checkquota to check.
Also look at the dump file size to determine if it's the cause.

--------------
08.29

1. Found that the way time unit is calculated now cannot be generated to cases where k is not equal to 2pi. Should use T = 2pi/sqrt(kg) instead. Changed that under /tiger/linear

-------------
09.04

1. Edited .bashrc to load python 2.7 automatically. Installed python package Pillow with pip install --user Pillow, so that bview client can be used on cluster.
"you can install package without root privilege. use pip with '--user' option. Then pip install package in your '~/.local/' folder."



-------------
11.05.19

1. Installed Basilisk on Della. Error when doing `cd $BASILISK/gl; make`, which also appeared on stampede2.

```
cc -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2   -c -o trackball.o trackball.c
cc -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -Igl2ps -c gl2ps/gl2ps.c -o gl2ps.o
cc -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2   -c -o utils.o utils.c
cc -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2   -c -o polygonize.o polygonize.c
cc -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2   -c -o og_font.o og_font.c
cc -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2   -c -o og_stroke_mono_roman.o og_stroke_mono_roman.c
ar cr libglutils.a trackball.o gl2ps.o utils.o polygonize.o og_font.o og_stroke_mono_roman.o
cc -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2   -c -o fb_osmesa.o fb_osmesa.c
ar cr libfb_osmesa.a fb_osmesa.o
cc -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2   -c -o fb_glx.o fb_glx.c
cc -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2   -c -o OffscreenContextGLX.o OffscreenContextGLX.c
In file included from OffscreenContextGLX.c:53:0:
system-gl.h:2:21: fatal error: GL/glew.h: No such file or directory
 #include <GL/glew.h>
                     ^
compilation terminated.
make: *** [OffscreenContextGLX.o] Error 1
```

2. Della and Tiger have two different versions of openmpi.

------------
11.18.19

1. `pip install --user Pillow` on both adroit and della so that bview server is now good to run.
2. Add to ./bashrc on adroit to load python 2.7 