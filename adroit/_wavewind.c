#line 1 "wavewind-cpp.c"
#line 1 "<built-in>"
#line 1 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 1 "<command-line>"
#line 1 "wavewind-cpp.c"
#if _XOPEN_SOURCE < 700
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#if _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif



#line 1 "/home/jiarongw/basilisk/src/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#if _OPENMP
# include <omp.h>
#elif 1
# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe
#endif

#if _CADNA
# include <cadna.h>
#endif

#if __cplusplus
# define delete delete_qcc
# define right right_qcc
# define left left_qcc
# define norm norm_qcc
# define new new_qcc
#endif

#define pi 3.14159265358979
#undef HUGE
#define HUGE ((double)1e30)
#define nodata HUGE
#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : nodata)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))
#define swap(type,a,b) { type tmp = a; a = b; b = tmp; }
#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout

#if 1
static FILE * qstderr (void);
static FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 72

# define system(command) (pid() == 0 ? system(command) : 0)
#else
# define qstderr() stderr
# define qstdout() stdout
# define ferr stderr
# define fout stdout
# define not_mpi_compatible()
#endif



#if MTRACE

struct {
  FILE * fp;
  size_t total, max;
  size_t overhead, maxoverhead;
  size_t nr;
  size_t startrss, maxrss;
  char * fname;
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
      c, f->id, pmtrace.nr, pmtrace.total, f->total);
#if _GNU_SOURCE
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
#endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
       const char * func, const char * file, int line,
       char c)
{
  assert (d != NULL);
  d->id = pmfunc_index(func, file, line);
  d->size = size;
  pmfunc * f = &pmfuncs[d->id - 1];
  f->total += size;
  if (f->total > f->max)
    f->max = f->total;
  pmtrace.total += size;
  pmtrace.overhead += sizeof(pmdata);
  if (pmtrace.total > pmtrace.max) {
    pmtrace.max = pmtrace.total;
    pmtrace.maxoverhead = pmtrace.overhead;
  }
  pmfunc_trace (f, c);
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", qstderr());
    if (d->size == 0)
      fputs (", possible double free()", qstderr());
    else
      fputs (", not traced?", qstderr());
    fputs (", aborting...\n", qstderr());
    abort();
    return ptr;
  }
  else {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (qstderr(), "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (qstderr(), "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
    return d;
  }
}

static void * pmalloc (size_t size,
         const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),
         size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
         const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
   const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),
           sizeof(pmdata) + size),
         size, func, file, line, '>');
}

static void pfree (void * ptr,
     const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
         const char * func, const char * file, int line)
{
  char * d = (char *) pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

#if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
#endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = (const pmfunc *) a, * p2 = (const pmfunc *) b;
  return p1->total < p2->total;
}

static void pmfuncs_free()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
#if MTRACE < 3
  fprintf (qstderr(),
    "*** MTRACE: max resident  set size: %10ld bytes\n"
    "*** MTRACE: max traced memory size: %10ld bytes"
    " (tracing overhead %.1g%%)\n"
    "%10s    %20s   %s\n",
    pmtrace.maxrss*1024,
    pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
    "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (qstderr(), "%10ld    %20s   %s:%d\n",
      p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
      "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
      "total(\"%s\") w l t 'total'",
      fname,
      pmtrace.startrss*1024.,
      pmtrace.startrss*1024.,
      fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
        ",func(\"%s\",%d) w l t '%s'",
        fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (qstderr(),
      "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
      fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
#endif

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (qstderr(), "%s:%d: error: %ld bytes leaked here\n",
        p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
#if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", qstderr());
#endif
    pmfuncs_free();
  }
}

#else
# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)
#endif







typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,__LINE__));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
  return p;
}



#if TRACE == 1
#include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = pstrdup (func,__func__,__FILE__,__LINE__);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

static void trace_pop (Trace * t, const char * func)
{
  assert (t->stack.len > 0);
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    pfree (*func,__func__,__FILE__,__LINE__);
  pfree (t->index.p,__func__,__FILE__,__LINE__);
  pfree (t->stack.p,__func__,__FILE__,__LINE__);
}

static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}






# define trace(func, file, line) trace_push (&trace_func, func)
# define end_trace(func, file, line) trace_pop (&trace_func, func)

#elif TRACE

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
#if 1
  double min, max;
#endif
} TraceIndex;

struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
         double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {pstrdup(func,__func__,__FILE__,__LINE__), pstrdup(file,__func__,__FILE__,__LINE__), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));




}

static void end_trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  assert (Trace.stack.len >= 2*sizeof(double));
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];




  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

#if 1
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
#endif

void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
#if 1
  qsort (index->p, len, sizeof(TraceIndex), compar_func);
  double tot[len], self[len], min[len], max[len];
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    tot[i] = t->total, self[i] = t->self;
  MPI_Reduce (self, min, len, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (self, max, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? self : MPI_IN_PLACE,
       self, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? tot : MPI_IN_PLACE,
       tot, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total = 0.;
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    t->total = tot[i]/npe(), t->self = self[i]/npe(),
      t->max = max[i], t->min = min[i], total += t->self;
#endif
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
        t->calls, t->total, t->self, t->self*100./total);
#if 1
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
#endif
      fprintf (fp, "   %s():%s:%d\n", t->func, t->file, t->line);
    }
  fflush (fp);
  array_free (index);
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    t->calls = t->total = t->self = 0.;
}

static void trace_off()
{
  trace_print (fout, 0.);

  int i, len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    pfree (t->func,__func__,__FILE__,__LINE__), pfree (t->file,__func__,__FILE__,__LINE__);

  pfree (Trace.index.p,__func__,__FILE__,__LINE__);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,__LINE__);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else
# define trace(...)
# define end_trace(...)
#endif



#if _OPENMP

#define OMP(x) _Pragma(#x)
#define tid() omp_get_thread_num()
#define pid() 0
#define npe() omp_get_num_threads()
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)

#elif 1

#define OMP(x)

static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  assert (!in_prof); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 645

#define prof_stop()\
  assert (in_prof); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 650


#if FAKE_MPI
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)
#else

int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{ trace ("mpi_all_reduce0", "/home/jiarongw/basilisk/src/common.h", 659);
  { int _ret =  MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm); end_trace("mpi_all_reduce0", "/home/jiarongw/basilisk/src/common.h", 660);  return _ret; }
 end_trace("mpi_all_reduce0", "/home/jiarongw/basilisk/src/common.h", 661); }
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 669

#define mpi_all_reduce_double(v,op) {\
  prof_start ("mpi_all_reduce");\
  double global, tmp = v;\
  mpi_all_reduce0 (&tmp, &global, 1, MPI_DOUBLE, op, MPI_COMM_WORLD);\
  v = global;\
  prof_stop();\
}\

#line 677


#endif

#define QFILE FILE

static FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

static FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL){
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = qstderr();
      fout = qstdout();
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else

#define OMP(x)
#define tid() 0
#define pid() 0
#define npe() 1
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)

#endif

void init_solver()
{
#if _CADNA
  cadna_init (-1);
#endif
#if 1
  mpi_init();
#elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
#endif
}

#define OMP_PARALLEL() OMP(omp parallel)

#define NOT_UNUSED(x) (void)(x)

#define VARIABLES ;
#define val(a,k,l,m) data(k,l,m)[a.i]

double _val_higher_dimension = 0.;
#define _val_higher_dimension(x,a,b,c) _val_higher_dimension
#line 803 "/home/jiarongw/basilisk/src/common.h"
#if (_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
# define enable_fpe(flags) feenableexcept (flags)
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  assert (sizeof (int64_t) == sizeof (double));
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined ((double) DBL_MAX)
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;


int N = 64;




typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;




} vector;

typedef struct {
  vector x;

  vector y;




} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;
#line 882 "/home/jiarongw/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  {
#line 885

    norm += sq(n->x);
#line 885

    norm += sq(n->y);}
  norm = sqrt(norm);
  {
#line 888

    n->x /= norm;
#line 888

    n->y /= norm;}
}

struct _origin { double x, y, z; };

void origin (struct _origin p) {
  X0 = p.x; Y0 = p.y; Z0 = p.z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }






  enum { right, left, top, bottom };



int nboundary = 2*2;



#define dirichlet(x) (2.*(x) - val(_s,0,0,0))
#define dirichlet_homogeneous() (- val(_s,0,0,0))
#define neumann(x) (Delta*(x) + val(_s,0,0,0))
#define neumann_homogeneous() (val(_s,0,0,0))

double * _constant = NULL;
extern size_t datasize;
typedef struct _Point Point;

#line 1 "/home/jiarongw/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,__LINE__);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 47 "/home/jiarongw/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 927 "/home/jiarongw/basilisk/src/common.h"



typedef struct {

#line 932 "/home/jiarongw/basilisk/src/common.h"

  double (** boundary) (Point, Point, scalar);
  double (** boundary_homogeneous) (Point, Point, scalar);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;




  } d;
  vector v;
  bool face, nodump;

#line 17 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);

#line 8 "/home/jiarongw/basilisk/src/grid/tree-common.h"

  void (* refine) (Point, scalar);

#line 94 "/home/jiarongw/basilisk/src/grid/tree-common.h"

  void (* coarsen) (Point, scalar);

#line 26 "/home/jiarongw/basilisk/src/vof.h"

  scalar * tracers;
  bool inverse;

#line 76 "/home/jiarongw/basilisk/src/fractions.h"

  vector n;

#line 20 "/home/jiarongw/basilisk/src/iforce.h"

  scalar phi;

#line 21 "/home/jiarongw/basilisk/src/tension.h"

  double sigma;

} _Attributes;
_Attributes * _attribute;
#line 930 "/home/jiarongw/basilisk/src/common.h"























int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  if (list) for (scalar s = *list, *_i0 = list; ((scalar *)&s)->i >= 0; s = *++_i0) ns++;
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  if (list) for (scalar t = *list, *_i1 = list; ((scalar *)&t)->i >= 0; t = *++_i1)
    if (t.i == s.i)
      return list;
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    if (l) for (scalar s1 = *l, *_i2 = l; ((scalar *)&s1)->i >= 0; s1 = *++_i2)
      if (s1.i == s.i)
 return true;
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    if (l) for (scalar s = *l, *_i3 = l; ((scalar *)&s)->i >= 0; s = *++_i3)
      list = list_append (list, s);
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  if (l2) for (scalar s = *l2, *_i4 = l2; ((scalar *)&s)->i >= 0; s = *++_i4)
    l3 = list_append (l3, s);
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  if (l) for (scalar s = *l, *_i5 = l; ((scalar *)&s)->i >= 0; s = *++_i5)
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  if (list) for (vector v = *list, *_i6 = list; ((scalar *)&v)->i >= 0; v = *++_i6) nv++;
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,__LINE__);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  if (list) for (vector w = *list, *_i7 = list; ((scalar *)&w)->i >= 0; w = *++_i7) {
    bool id = true;
    {
#line 1033

      if (w.x.i != v.x.i)
 id = false;
#line 1033

      if (w.y.i != v.y.i)
 id = false;}
    if (id)
      return list;
  }
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    if (l) for (vector v = *l, *_i8 = l; ((scalar *)&v)->i >= 0; v = *++_i8)
      list = vectors_append (list, v);
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
    {
#line 1056
 {
      assert (s->i >= 0);
      v.x = *s++;
    }
#line 1056
 {
      assert (s->i >= 0);
      v.y = *s++;
    }}
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  if (list) for (tensor t = *list, *_i9 = list; ((scalar *)&t)->i >= 0; t = *++_i9) nt++;
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,__LINE__);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
    {
#line 1087
 {
      assert (v->x.i >= 0);
      t.x = *v++;
    }
#line 1087
 {
      assert (v->y.i >= 0);
      t.y = *v++;
    }}
    list = tensors_append (list, t);
  }
  return list;
}

scalar * all = NULL;



scalar (* init_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
vector (* init_face_vector) (vector, const char *);





typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
void _init_solver (void);



#if 1
static double mpi_time = 0.;
#endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
#if 1
  t.tm = mpi_time;
#endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



vector zerof= {{_NVARMAX + 0},{_NVARMAX + 1}};
vector unityf= {{_NVARMAX + 2},{_NVARMAX + 3}};
scalar unity= {_NVARMAX + 4};
scalar zeroc= {_NVARMAX + 5};



 vector fm = {{_NVARMAX + 2},{_NVARMAX + 3}};
 scalar cm = {(_NVARMAX + 4)};



static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}



void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = HUGE;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 swap (double, m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 swap (double, m[k][indxr[l]], m[k][indxc[l]]);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}
#line 14 "wavewind-cpp.c"
#line 1 "grid/quadtree.h"
#line 1 "/home/jiarongw/basilisk/src/grid/quadtree.h"


#line 1 "grid/tree.h"
#line 1 "/home/jiarongw/basilisk/src/grid/tree.h"
#line 1 "grid/mempool.h"
#line 1 "/home/jiarongw/basilisk/src/grid/mempool.h"





typedef struct _Pool Pool;

struct _Pool {
  Pool * next;
};

typedef struct {
  char * first, * lastb;
  size_t size;
  size_t poolsize;
  Pool * pool, * last;
} Mempool;

typedef struct {
  char * next;
} FreeBlock;

Mempool * mempool_new (size_t poolsize, size_t size)
{

  assert (poolsize % 8 == 0);
  assert (size >= sizeof(FreeBlock));


  poolsize = min(1 << 20, poolsize + sizeof(Pool));
  Mempool * m = ((Mempool *) pcalloc (1, sizeof(Mempool),__func__,__FILE__,__LINE__));
  m->poolsize = poolsize;
  m->size = size;
  return m;
}

void mempool_destroy (Mempool * m)
{
  Pool * p = m->pool;
  while (p) {
    Pool * next = p->next;
    pfree (p,__func__,__FILE__,__LINE__);
    p = next;
  }
  pfree (m,__func__,__FILE__,__LINE__);
}

void * mempool_alloc (Mempool * m)
{
  if (!m->first) {

    Pool * p = (Pool *) pmalloc (m->poolsize,__func__,__FILE__,__LINE__);
    p->next = NULL;
    if (m->last)
      m->last->next = p;
    else
      m->pool = p;
    m->last = p;
    m->first = m->lastb = ((char *)m->last) + sizeof(Pool);
    FreeBlock * b = (FreeBlock *) m->first;
    b->next = NULL;
  }
  void * ret = m->first;
  FreeBlock * b = (FreeBlock *) ret;
  char * next = b->next;
  if (!next) {
    m->lastb += m->size;
    next = m->lastb;
    if (next + m->size > ((char *) m->last) + m->poolsize)
      next = NULL;
    else {
      FreeBlock * b = (FreeBlock *) next;
      b->next = NULL;
    }
  }
  m->first = next;
#if TRASH
  double * v = (double *) ret;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  return ret;
}

void * mempool_alloc0 (Mempool * m)
{
  void * ret = mempool_alloc (m);
  memset (ret, 0, m->size);
  return ret;
}

void mempool_free (Mempool * m, void * p)
{
#if TRASH
  double * v = (double *) p;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  FreeBlock * b = (FreeBlock *) p;
  b->next = m->first;
  m->first = (char *) p;
}
#line 2 "/home/jiarongw/basilisk/src/grid/tree.h"
#line 22 "/home/jiarongw/basilisk/src/grid/tree.h"
typedef struct {
  unsigned short flags;

  unsigned short neighbors;
  int pid;
} Cell;

enum {
  active = 1 << 0,
  leaf = 1 << 1,
  border = 1 << 2,
  vertex = 1 << 3,
  user = 4,

  face_x = 1 << 0

  , face_y = 1 << 1




};

#define is_active(cell) ((cell).flags & active)
#define is_leaf(cell) ((cell).flags & leaf)
#define is_coarse() ((cell).neighbors > 0)
#define is_border(cell) ((cell).flags & border)
#define is_local(cell) ((cell).pid == pid())



typedef struct {
  int i;

  int j;




} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i;

  int j;




  int level, flags;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;




static char * new_refarray (size_t len, size_t size) {
  return (char *) pcalloc (len + 1, size,__func__,__FILE__,__LINE__);
}

static void refarray (void * p, size_t len, size_t size) {
  int * refcount = (int *)(((char *)p) + len*size);
  (*refcount)++;
}

static bool unrefarray (void * p, size_t len, size_t size) {
  int * refcount = (int *)(((char *)p) + len*size);
  (*refcount)--;
  if (*refcount == 0) {
    pfree (p,__func__,__FILE__,__LINE__);
    return true;
  }
  return false;
}




typedef struct {



  char *** m;



  Mempool * pool;
  int nc;
  int len;
} Layer;

static size_t _size (size_t depth)
{
  return (1 << depth) + 2*2;
}

static size_t poolsize (size_t depth, size_t size)
{




  return sq(_size(depth))*size;



}


#line 139

static inline
void assign_periodic_x (void ** m, int i, int nl, void * b)
{
  m[i] = b;
  if (Period.x) {
    for (int j = i; j < nl + 2*2; j += nl)
      m[j] = b;
    for (int j = i - nl; j >= 0; j -= nl)
      m[j] = b;
  }
}
#line 139

static inline
void assign_periodic_y (void ** m, int i, int nl, void * b)
{
  m[i] = b;
  if (Period.y) {
    for (int j = i; j < nl + 2*2; j += nl)
      m[j] = b;
    for (int j = i - nl; j >= 0; j -= nl)
      m[j] = b;
  }
}

static Layer * new_layer (int depth)
{
  Layer * l = ((Layer *) pmalloc ((1)*sizeof(Layer),__func__,__FILE__,__LINE__));
  l->len = _size (depth);
  if (depth == 0)
    l->pool = NULL;
  else {
    size_t size = sizeof(Cell) + datasize;


    l->pool = mempool_new (poolsize (depth, size), (1 << 2)*size);
  }



  l->m = ((char ** *) pcalloc (l->len, sizeof(char **),__func__,__FILE__,__LINE__));



  l->nc = 0;
  return l;
}

static void destroy_layer (Layer * l)
{
  if (l->pool)
    mempool_destroy (l->pool);
  pfree (l->m,__func__,__FILE__,__LINE__);
  pfree (l,__func__,__FILE__,__LINE__);
}
#line 199 "/home/jiarongw/basilisk/src/grid/tree.h"
static void layer_add_row (Layer * l, int i, int j)
{
  if (!l->m[i]) {
    assign_periodic_x ((void **) l->m, i, l->len - 2*2,
         (void *) new_refarray (l->len, sizeof (char *)));
    l->nc++;
  }
  refarray (l->m[i], l->len, sizeof(char *));






}

static bool layer_remove_row (Layer * l, int i, int j)
{




  if (unrefarray (l->m[i], l->len, sizeof (char *))) {
    assign_periodic_x ((void **) l->m, i, l->len - 2*2, NULL);
    if (--l->nc == 0) {
      destroy_layer (l);
      return true;
    }
    assert (l->nc >= 0);
  }
  return false;
}




typedef struct {
  Grid g;
  Layer ** L;

  Cache leaves;
  Cache faces;
  Cache vertices;
  Cache refined;
  CacheLevel * active;
  CacheLevel * prolongation;
  CacheLevel * boundary;

  CacheLevel * restriction;

  bool dirty;
} Tree;



struct _Point {

  int i;

  int j;




  int level;
};
static Point last_point;



static void cache_level_append (CacheLevel * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (IndexLevel *) prealloc (c->p, (c->nm)*sizeof(IndexLevel),__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;




  c->n++;
}

static void cache_level_shrink (CacheLevel * c)
{
  if (c->nm > (c->n/128 + 1)*128) {
    c->nm = (c->n/128 + 1)*128;
    assert (c->nm > c->n);
    c->p = (IndexLevel *) prealloc (c->p, sizeof (Index)*c->nm,__func__,__FILE__,__LINE__);
  }
}

static void cache_append (Cache * c, Point p, unsigned short flags)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (Index *) prealloc (c->p, (c->nm)*sizeof(Index),__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;




  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

void cache_shrink (Cache * c)
{
  cache_level_shrink ((CacheLevel *)c);
}
#line 342 "/home/jiarongw/basilisk/src/grid/tree.h"
#define allocated(k,l,n) (point.i+k >= 0 &&\
         point.i+k < (1 << point.level) + 2*2 &&\
         ((Tree *)grid)->L[point.level]->m[point.i+k] &&\
         point.j+l >= 0 &&\
         point.j+l < (1 << point.level) + 2*2 &&\
         ((Tree *)grid)->L[point.level]->m[point.i+k][point.j+l])\

#line 348


#define NEIGHBOR(k,l,n) (((Tree *)grid)->L[point.level]->m[point.i+k][point.j+l])
#define PARENT(k,l,n) (((Tree *)grid)->L[point.level-1]->m[(point.i+2)/2+k]\
    [(point.j+2)/2+l])\

#line 353

#define allocated_child(k,l,n) (level < depth() &&\
         point.i > 0 && point.i <= (1 << level) + 2 &&\
         point.j > 0 && point.j <= (1 << level) + 2 &&\
         ((Tree *)grid)->L[point.level+1]->m[2*point.i-2 +k]\
   && ((Tree *)grid)->L[point.level+1]->m[2*point.i-2 +k][2*point.j-2 +l])\

#line 359

#define CHILD(k,l,n) (((Tree *)grid)->L[point.level+1]->m[2*point.i-2 +k]\
   [2*point.j-2 +l])\

#line 362

#line 402 "/home/jiarongw/basilisk/src/grid/tree.h"
#define CELL(m) (*((Cell *)(m)))


#define depth() (grid->depth)
#define aparent(k,l,n) CELL(PARENT(k,l,n))
#define child(k,l,n) CELL(CHILD(k,l,n))


#define cell CELL(NEIGHBOR(0,0,0))
#define neighbor(k,l,n) CELL(NEIGHBOR(k,l,n))
#define neighborp(l,m,n) (Point) {\
    point.i + l,\
\
    point.j + m,\
\
\
\
\
    point.level }\

#line 421



#define data(k,l,n) ((double *) (NEIGHBOR(k,l,n) + sizeof(Cell)))
#define fine(a,k,l,n) ((double *) (CHILD(k,l,n) + sizeof(Cell)))[a.i]
#define coarse(a,k,l,n) ((double *) (PARENT(k,l,n) + sizeof(Cell)))[a.i]

#define POINT_VARIABLES\
  VARIABLES\
  int level = point.level; NOT_UNUSED(level);\
\
\
\
  struct { int x, y; } child = {\
    2*((point.i+2)%2)-1, 2*((point.j+2)%2)-1\
  };\
\
\
\
\
\
  NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2;\
\
  parent.j = (point.j + 2)/2;\
\

#line 457


#line 1 "grid/foreach_cell.h"
#line 1 "/home/jiarongw/basilisk/src/grid/foreach_cell.h"
#line 66 "/home/jiarongw/basilisk/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point;\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
  Point root = {2,2,0};\
\
\
\
  foreach_cell_root (root)\

#line 134

#define end_foreach_cell() end_foreach_cell_root() }

#define foreach_cell_all() {\
  Point root = { .level = 0 };\
  for (root.i = 2*Period.x; root.i <= 2*(2 - Period.x); root.i++)\
\
    for (root.j = 2*Period.y; root.j <= 2*(2 - Period.y); root.j++)\
\
\
\
\
 foreach_cell_root (root)\

#line 147

#define end_foreach_cell_all() end_foreach_cell_root() }

#define foreach_cell_post_root(condition, root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point;\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
\
\
\
\
\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
\
      default: {\
        POINT_VARIABLES;\
\

#line 244

#define end_foreach_cell_post_root()\
      }\
      }\
    }\
  }\

#line 250


#define foreach_cell_post(condition)\
  {\
\
\
\
    Point root = {2,2,0};\
\
\
\
    foreach_cell_post_root(condition, root)\

#line 262

#define end_foreach_cell_post() end_foreach_cell_post_root() }

#define foreach_cell_post_all(condition) {\
  Point root = { .level = 0 };\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 281

#define end_foreach_leaf() } continue; } end_foreach_cell()
#line 460 "/home/jiarongw/basilisk/src/grid/tree.h"
#line 477 "/home/jiarongw/basilisk/src/grid/tree.h"
#define foreach_child() {\
  int _i = 2*point.i - 2, _j = 2*point.j - 2;\
  point.level++;\
  for (int _k = 0; _k < 2; _k++) {\
    point.i = _i + _k;\
    for (int _l = 0; _l < 2; _l++) {\
      point.j = _j + _l;\
      POINT_VARIABLES;\

#line 485

#define end_foreach_child()\
    }\
  }\
  point.i = (_i + 2)/2; point.j = (_j + 2)/2;\
  point.level--;\
}\

#line 492

#define foreach_child_break() _k = _l = 2
#line 523 "/home/jiarongw/basilisk/src/grid/tree.h"
#define is_refined_check() ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) &&\
    point.i > 0 && point.i < (1 << level) + 2*2 - 1\
\
    && point.j > 0 && point.j < (1 << level) + 2*2 - 1\
\
\
\
\
    )\

#line 532


#define foreach_cache(_cache) {\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
\
\
\
  Point point = {2,2,0};\
\
\
\
  int _k; unsigned short _flags; NOT_UNUSED(_flags);\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
\
\
    point.level = _cache.p[_k].level;\
    _flags = _cache.p[_k].flags;\
    POINT_VARIABLES;\

#line 557

#define end_foreach_cache() } } }

#define foreach_cache_level(_cache,_l) {\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
\
\
\
  Point point = {2,2,0};\
\
\
\
  point.level = _l;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
\
\
    POINT_VARIABLES;\

#line 582

#define end_foreach_cache_level() } } }

#define foreach_boundary_level(_l) {\
  if (_l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _boundary = ((Tree *)grid)->boundary[_l];\
    foreach_cache_level (_boundary,_l)\

#line 590

#define end_foreach_boundary_level() end_foreach_cache_level(); }}



#define foreach_boundary(_b) {\
  for (int _l = depth(); _l >= 0; _l--)\
    foreach_boundary_level(_l) {\
      if ((- cell.pid - 1) == _b)\
 for (int _d = 0; _d < 2; _d++) {\
   for (int _i = -1; _i <= 1; _i += 2) {\
     if (_d == 0) ig = _i; else if (_d == 1) jg = _i; else kg = _i;\
     if (allocated(-ig,-jg,-kg) &&\
  is_leaf (neighbor(-ig,-jg,-kg)) &&\
  !(neighbor(-ig,-jg,-kg).pid < 0) &&\
  is_local(neighbor(-ig,-jg,-kg))) {\
       point.i -= ig; x -= ig*Delta/2.;\
\
       point.j -= jg; y -= jg*Delta/2.;\
\
\
\
\

#line 613

#define end_foreach_boundary()\
       point.i += ig; x += ig*Delta/2.;\
\
       point.j += jg; y += jg*Delta/2.;\
\
\
\
\
            }\
   }\
   ig = jg = kg = 0;\
 }\
    } end_foreach_boundary_level(); }\

#line 627


#define foreach_halo(_name,_l) {\
  if (_l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _cache = ((Tree *)grid)->_name[_l];\
    foreach_cache_level (_cache, _l)\

#line 634

#define end_foreach_halo() end_foreach_cache_level(); }}

#line 1 "grid/neighbors.h"
#line 1 "/home/jiarongw/basilisk/src/grid/neighbors.h"
#line 17 "/home/jiarongw/basilisk/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j;\
  for (int _k = - _nn; _k <= _nn; _k++) {\
    point.i = _i + _k;\
    for (int _l = - _nn; _l <= _nn; _l++) {\
      point.j = _j + _l;\
      POINT_VARIABLES;\

#line 25

#define end_foreach_neighbor()\
    }\
  }\
  point.i = _i; point.j = _j;\
}\

#line 31

#define foreach_neighbor_break() _k = _l = _nn + 1
#line 638 "/home/jiarongw/basilisk/src/grid/tree.h"

static inline bool has_local_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 640 "/home/jiarongw/basilisk/src/grid/tree.h"

   { foreach_child()
    if (is_local(cell))
      return true; end_foreach_child(); }
  return false;
}

static inline void cache_append_face (Point point, unsigned short flags)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 648 "/home/jiarongw/basilisk/src/grid/tree.h"

  Tree * q = ((Tree *)grid);
  cache_append (&q->faces, point, flags);

  if (!(cell.flags & vertex)) {
    cache_append (&q->vertices, point, 0);
    cell.flags |= vertex;
  }
  {
#line 656

    if ((flags & face_y) && !(neighbor(1,0,).flags & vertex)) {
      cache_append (&q->vertices, neighborp(1,0,), 0);
      neighbor(1,0,).flags |= vertex;
    }
#line 656

    if ((flags & face_x) && !(neighbor(0,1,).flags & vertex)) {
      cache_append (&q->vertices, neighborp(0,1,), 0);
      neighbor(0,1,).flags |= vertex;
    }}
#line 671 "/home/jiarongw/basilisk/src/grid/tree.h"
}



void check_periodic (Tree * q)
{
#line 702 "/home/jiarongw/basilisk/src/grid/tree.h"
}

static void update_cache_f (void)
{
  Tree * q = ((Tree *)grid);

  check_periodic (q);

   { foreach_cache (q->vertices){

#line 710 "/home/jiarongw/basilisk/src/grid/tree.h"

    if (level <= depth() && allocated(0,0,0))
      cell.flags &= ~vertex; } end_foreach_cache(); }


  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n =
      q->boundary[l].n = q->restriction[l].n = 0;

  const unsigned short fboundary = 1 << user;
   { foreach_cell(){

#line 721 "/home/jiarongw/basilisk/src/grid/tree.h"
 {



    if (is_local(cell) && is_active(cell)) {


      cache_level_append (&q->active[level], point);
    }
#line 745 "/home/jiarongw/basilisk/src/grid/tree.h"
    if (!(cell.pid < 0)) {

       { foreach_neighbor (2)
 if (allocated(0,0,0) && (cell.pid < 0) && !(cell.flags & fboundary)) {
   cache_level_append (&q->boundary[level], point);
   cell.flags |= fboundary;
 } end_foreach_neighbor(); }
    }

    else if (level > 0 && is_local(aparent(0,0,0)))
      cache_level_append (&q->restriction[level], point);

    if (is_leaf (cell)) {
      if (is_local(cell)) {
 cache_append (&q->leaves, point, 0);

 unsigned short flags = 0;
 {
#line 762

   if ((neighbor(-1,0,).pid < 0) || (!is_leaf(neighbor(-1,0,)) && !neighbor(-1,0,).neighbors && neighbor(-1,0,).pid >= 0) ||
       is_leaf(neighbor(-1,0,)))
     flags |= face_x;
#line 762

   if ((neighbor(0,-1,).pid < 0) || (!is_leaf(neighbor(0,-1,)) && !neighbor(0,-1,).neighbors && neighbor(0,-1,).pid >= 0) ||
       is_leaf(neighbor(0,-1,)))
     flags |= face_y;}
 if (flags)
   cache_append (&q->faces, point, flags);
 {
#line 768

   if ((neighbor(1,0,).pid < 0) || (!is_leaf(neighbor(1,0,)) && !neighbor(1,0,).neighbors && neighbor(1,0,).pid >= 0) ||
       (!is_local(neighbor(1,0,)) && is_leaf(neighbor(1,0,))))
     cache_append (&q->faces, neighborp(1,0,), face_x);
#line 768

   if ((neighbor(0,1,).pid < 0) || (!is_leaf(neighbor(0,1,)) && !neighbor(0,1,).neighbors && neighbor(0,1,).pid >= 0) ||
       (!is_local(neighbor(0,1,)) && is_leaf(neighbor(0,1,))))
     cache_append (&q->faces, neighborp(0,1,), face_y);}

 for (int i = 0; i <= 1; i++)

   for (int j = 0; j <= 1; j++)




       if (!(neighbor(i,j,k).flags & vertex)) {
  cache_append (&q->vertices, neighborp(i,j,k), 0);
  neighbor(i,j,k).flags |= vertex;
       }

        if (cell.neighbors > 0)
   cache_level_append (&q->prolongation[level], point);
      }
      else if (!(cell.pid < 0) || is_local(aparent(0,0,0))) {

 unsigned short flags = 0;
 {
#line 791

   if (allocated(-1,0,) &&
       is_local(neighbor(-1,0,)) && (!is_leaf(neighbor(-1,0,)) && !neighbor(-1,0,).neighbors && neighbor(-1,0,).pid >= 0))
     flags |= face_x;
#line 791

   if (allocated(0,-1,) &&
       is_local(neighbor(0,-1,)) && (!is_leaf(neighbor(0,-1,)) && !neighbor(0,-1,).neighbors && neighbor(0,-1,).pid >= 0))
     flags |= face_y;}
 if (flags)
   cache_append_face (point, flags);
 {
#line 797

   if (allocated(1,0,) && is_local(neighbor(1,0,)) &&
       (!is_leaf(neighbor(1,0,)) && !neighbor(1,0,).neighbors && neighbor(1,0,).pid >= 0))
     cache_append_face (neighborp(1,0,), face_x);
#line 797

   if (allocated(0,1,) && is_local(neighbor(0,1,)) &&
       (!is_leaf(neighbor(0,1,)) && !neighbor(0,1,).neighbors && neighbor(0,1,).pid >= 0))
     cache_append_face (neighborp(0,1,), face_y);}
      }

      continue;

    }
  } } end_foreach_cell(); }


  cache_shrink (&q->leaves);
  cache_shrink (&q->faces);
  cache_shrink (&q->vertices);
  for (int l = 0; l <= depth(); l++) {
    cache_level_shrink (&q->active[l]);
    cache_level_shrink (&q->prolongation[l]);
    cache_level_shrink (&q->boundary[l]);
    cache_level_shrink (&q->restriction[l]);
}

  q->dirty = false;


  for (int l = depth(); l >= 0; l--)
     { foreach_boundary_level (l){

#line 823 "/home/jiarongw/basilisk/src/grid/tree.h"

      cell.flags &= ~fboundary; } end_foreach_boundary_level(); }



  grid->n = q->leaves.n;

#if !1
  grid->tn = grid->n;
  grid->maxdepth = grid->depth;
#endif
}

#define foreach() { if (((Tree *)grid)->dirty) update_cache_f(); }; foreach_cache(((Tree *)grid)->leaves)
#define end_foreach() end_foreach_cache()

#define foreach_face_generic()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->faces) 
#line 840

#define end_foreach_face_generic() end_foreach_cache()

#define is_face_x() (_flags & face_x)

#define is_face_y() (_flags & face_y)





#define foreach_vertex()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->vertices) {\
    x -= Delta/2.;\
\
    y -= Delta/2.;\
\
\
\
\

#line 862

#define end_foreach_vertex() } end_foreach_cache()
#line 874 "/home/jiarongw/basilisk/src/grid/tree.h"
#define foreach_level(l) {\
  if (l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _active = ((Tree *)grid)->active[l];\
    foreach_cache_level (_active,l)\

#line 879

#define end_foreach_level() end_foreach_cache_level(); }}

#define foreach_coarse_level(l) foreach_level(l) if (!is_leaf(cell)) {
#define end_foreach_coarse_level() } end_foreach_level()

#define foreach_level_or_leaf(l) {\
  for (int _l1 = l; _l1 >= 0; _l1--)\
    foreach_level(_l1)\
      if (_l1 == l || is_leaf (cell)) {\

#line 889

#define end_foreach_level_or_leaf() } end_foreach_level(); }

#if TRASH
# undef trash
# define trash(list) reset(list, undefined)
#endif

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  Tree * q = ((Tree *)grid);

  for (int l = 0; l <= depth(); l++) {
    Layer * L = q->L[l];
    for (int i = 0; i < L->len; i++)
      if (L->m[i])





 for (int j = 0; j < L->len; j++)
   if (L->m[i][j])

     if (list) for (scalar s = *list, *_i10 = list; ((scalar *)&s)->i >= 0; s = *++_i10)
       if (!is_constant(s))
  ((double *)(L->m[i][j] + sizeof(Cell)))[s.i] = val;
#line 925 "/home/jiarongw/basilisk/src/grid/tree.h"
  }
}

#define cache_level_resize(name, a)\
{\
  for (int i = 0; i <= depth() - a; i++)\
    pfree (q->name[i].p,__func__,__FILE__,__LINE__);\
  pfree (q->name,__func__,__FILE__,__LINE__);\
  q->name = ((CacheLevel *) pcalloc (depth() + 1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));\
}\

#line 935


static void update_depth (int inc)
{
  Tree * q = ((Tree *)grid);
  grid->depth += inc;
  q->L = &(q->L[-1]);
  q->L = (Layer * *) prealloc (q->L, (grid->depth + 2)*sizeof(Layer *),__func__,__FILE__,__LINE__);
  q->L = &(q->L[1]);
  if (inc > 0)
    q->L[grid->depth] = new_layer (grid->depth);
  cache_level_resize (active, inc);
  cache_level_resize (prolongation, inc);
  cache_level_resize (boundary, inc);
  cache_level_resize (restriction, inc);
}

static void alloc_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 953 "/home/jiarongw/basilisk/src/grid/tree.h"

  if (point.level == grid->depth)
    update_depth (+1);
  else if (allocated_child(0,0,0))
    return;


  Layer * L = ((Tree *)grid)->L[point.level + 1];
  size_t len = sizeof(Cell) + datasize;
  char * b = (char *) mempool_alloc0 (L->pool);
  int nl = L->len - 2*2;
  int i = 2*point.i - 2;
  for (int k = 0; k < 2; k++, i++) {






    layer_add_row (L, i, 0);
    int j = 2*point.j - 2;
    for (int l = 0; l < 2; l++, j++) {
      assert (!L->m[i][j]);
      assign_periodic_y ((void **) L->m[i], j, nl, (void *) b);
      b += len;
    }
#line 991 "/home/jiarongw/basilisk/src/grid/tree.h"
  }

  int pid = cell.pid;
   { foreach_child() {
    cell.pid = pid;
#if TRASH
    if (all) for (scalar s = *all, *_i11 = all; ((scalar *)&s)->i >= 0; s = *++_i11)
      val(s,0,0,0) = undefined;
#endif
  } end_foreach_child(); }
}
#line 1020 "/home/jiarongw/basilisk/src/grid/tree.h"
static void free_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1021 "/home/jiarongw/basilisk/src/grid/tree.h"


  Layer * L = ((Tree *)grid)->L[point.level + 1];
  int i = 2*point.i - 2, nl = L->len - 2*2;
  assert (L->m[i][2*point.j - 2]);
  mempool_free (L->pool, L->m[i][2*point.j - 2]);
  for (int k = 0; k < 2; k++, i++) {
    int j = 2*point.j - 2;
    for (int l = 0; l < 2; l++, j++)
      assign_periodic_y ((void **) L->m[i], j, nl, NULL);
    if (layer_remove_row (L, i, j)) {
      assert (point.level + 1 == grid->depth);
      update_depth (-1);
    }
  }
}
#line 1060 "/home/jiarongw/basilisk/src/grid/tree.h"
void increment_neighbors (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1061 "/home/jiarongw/basilisk/src/grid/tree.h"

  ((Tree *)grid)->dirty = true;
  if (cell.neighbors++ == 0)
    alloc_children (point);
   { foreach_neighbor (2/2)
    if (cell.neighbors++ == 0)
      alloc_children (point); end_foreach_neighbor(); }
  cell.neighbors--;
}

void decrement_neighbors (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1072 "/home/jiarongw/basilisk/src/grid/tree.h"

  ((Tree *)grid)->dirty = true;
   { foreach_neighbor (2/2)
    if (allocated(0,0,0)) {
      cell.neighbors--;
      if (cell.neighbors == 0)
 free_children (point);
    } end_foreach_neighbor(); }
  if (cell.neighbors) {
    int pid = cell.pid;
     { foreach_child() {
      cell.flags = 0;
      cell.pid = pid;
    } end_foreach_child(); }
  }
}

static void apply_periodic_elem (char ** m, int len)
{
  if (m) {
    int end = len - 2;
    for (int k = 0; k < 2; k++) {
      m[k] = m[k + end - 2];
      m[end + k] = m[k + 2];
    }
  }
}

static void apply_periodic (Tree * q)
{
#line 1110 "/home/jiarongw/basilisk/src/grid/tree.h"
  if (Period.y) {
    for (int i = 0; i < q->L[0]->len; i++)
      for (int j = 0; j < q->L[0]->len; j++)
 q->L[0]->m[i][j] = q->L[0]->m[i][2];
    for (int l = 1; l <= depth(); l++) {
      Layer * L = q->L[l];
      for (int i = 0; i < L->len; i++)
 apply_periodic_elem (L->m[i], L->len);
    }
  }
#line 1135 "/home/jiarongw/basilisk/src/grid/tree.h"
}

void realloc_scalar (void)
{

  Tree * q = ((Tree *)grid);
  size_t newlen = sizeof(Cell) + datasize;
  size_t oldlen = newlen - sizeof(double);

  Layer * L = q->L[0];
  int len = L->len;
  for (int i = Period.x*2; i < len - Period.x*2; i++) {



    for (int j = Period.y*2; j < len - Period.y*2; j++) {

      L->m[i][j] = (char *) prealloc (L->m[i][j], (newlen)*sizeof(char),__func__,__FILE__,__LINE__);




    }

  }

  for (int l = 1; l <= depth(); l++) {
    Layer * L = q->L[l];
    int len = L->len;
    Mempool * oldpool = L->pool;
    L->pool = mempool_new (poolsize (l, newlen), (1 << 2)*newlen);
    for (int i = Period.x*2; i < len - Period.x*2; i += 2)
      if (L->m[i]) {
#line 1176 "/home/jiarongw/basilisk/src/grid/tree.h"
 for (int j = Period.y*2; j < len - Period.y*2; j += 2)
   if (L->m[i][j]) {

     char * new = (char *) mempool_alloc (L->pool);
     for (int k = 0; k < 2; k++)
       for (int o = 0; o < 2; o++) {
  memcpy (new, L->m[i+k][j+o], oldlen);
  L->m[i+k][j+o] = new;
  new += newlen;
       }
#line 1199 "/home/jiarongw/basilisk/src/grid/tree.h"
   }

      }
    mempool_destroy (oldpool);
  }
  apply_periodic (q);
  check_periodic (q);
}



#define VN v.x
#define VT v.y
#define VR v.z




#if 1
# define disable_fpe_for_mpi() disable_fpe (FE_DIVBYZERO|FE_INVALID)
# define enable_fpe_for_mpi() enable_fpe (FE_DIVBYZERO|FE_INVALID)
#else
# define disable_fpe_for_mpi()
# define enable_fpe_for_mpi()
#endif

static inline void no_restriction (Point point, scalar s);

static bool normal_neighbor (Point point, scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1228 "/home/jiarongw/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    {
#line 1230

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(i,0,) && !(neighbor(i,0,).pid < 0))) {
   Point neighbor = neighborp(i,0,);
   int id = (- cell.pid - 1);
   if (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
     val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s);
   if (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13) {
     scalar vn = VN;
     val(v.x,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x);

     scalar vt = VT;
     val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y);





   }
   return true;
 }
#line 1230

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(0,i,) && !(neighbor(0,i,).pid < 0))) {
   Point neighbor = neighborp(0,i,);
   int id = (- cell.pid - 1);
   if (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
     val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s);
   if (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13) {
     scalar vn = VN;
     val(v.y,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y);

     scalar vt = VT;
     val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x);





   }
   return true;
 }}
  return false;
}

static bool diagonal_neighbor_2D (Point point,
      scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1256 "/home/jiarongw/basilisk/src/grid/tree.h"


  for (int k = 1; k <= 2; k++)



      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(i,j,0) && (allocated(i,j,0) && !(neighbor(i,j,0).pid < 0)) &&
       allocated(i,0,0) && (neighbor(i,0,0).pid < 0) &&
       allocated(0,j,0) && (neighbor(0,j,0).pid < 0)) {
     Point n = neighborp(i,j,0),
       n1 = neighborp(i,0,0), n2 = neighborp(0,j,0);
     int id1 = (- neighbor(i,0,0).pid - 1), id2 = (- neighbor(0,j,0).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i14 = scalars; ((scalar *)&s)->i >= 0; s = *++_i14)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s) + _attribute[s.i].boundary[id2](n,n2,s) -
       val(s,i,j,0));
     if (vectors) for (vector v = *vectors, *_i15 = vectors; ((scalar *)&v)->i >= 0; v = *++_i15) {
       scalar vt = VT, vn = VN;
       val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.x) +
         _attribute[vn.i].boundary[id2](n,n2,v.x) -
         val(v.x,i,j,0));
       val(v.y,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.y) +
         _attribute[vt.i].boundary[id2](n,n2,v.y) -
         val(v.y,i,j,0));






     }
     return true;
   }

  return false;
}

static bool diagonal_neighbor_3D (Point point,
      scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1296 "/home/jiarongw/basilisk/src/grid/tree.h"

#line 1338 "/home/jiarongw/basilisk/src/grid/tree.h"
  return false;
}



#line 1342

static Point tangential_neighbor_x (Point point, bool * zn)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1344 "/home/jiarongw/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,j,) && !(neighbor(0,j,).pid < 0)) || (allocated(-1,j,) && !(neighbor(-1,j,).pid < 0))) {
 *zn = false;
 return neighborp(0,j,);
      }







    }
  return (Point){.level = -1};
}
#line 1342

static Point tangential_neighbor_y (Point point, bool * zn)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1344 "/home/jiarongw/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(j,0,) && !(neighbor(j,0,).pid < 0)) || (allocated(j,-1,) && !(neighbor(j,-1,).pid < 0))) {
 *zn = false;
 return neighborp(j,0,);
      }







    }
  return (Point){.level = -1};
}


static inline bool is_boundary_point (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1363 "/home/jiarongw/basilisk/src/grid/tree.h"

  return (cell.pid < 0);
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  disable_fpe_for_mpi();
  scalar * scalars = NULL;
  vector * vectors = NULL, * faces = NULL;
  if (list) for (scalar s = *list, *_i16 = list; ((scalar *)&s)->i >= 0; s = *++_i16)
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i) {
 if (_attribute[s.i].face)
   faces = vectors_add (faces, _attribute[s.i].v);
 else
   vectors = vectors_add (vectors, _attribute[s.i].v);
      }
      else if (_attribute[s.i].v.x.i < 0)
 scalars = list_add (scalars, s);
    }

   { foreach_boundary_level (l){

#line 1384 "/home/jiarongw/basilisk/src/grid/tree.h"
 {
    if (!normal_neighbor (point, scalars, vectors) &&
 !diagonal_neighbor_2D (point, scalars, vectors) &&
 !diagonal_neighbor_3D (point, scalars, vectors)) {

      if (scalars) for (scalar s = *scalars, *_i17 = scalars; ((scalar *)&s)->i >= 0; s = *++_i17)
 val(s,0,0,0) = undefined;
      if (vectors) for (vector v = *vectors, *_i18 = vectors; ((scalar *)&v)->i >= 0; v = *++_i18)
 {
#line 1392

   val(v.x,0,0,0) = undefined;
#line 1392

   val(v.y,0,0,0) = undefined;}
    }
    if (faces) {
      int id = (- cell.pid - 1);
      {
#line 1397

 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(i,0,) && !(neighbor(i,0,).pid < 0))) {
     Point neighbor = neighborp(i,0,);
     if (faces) for (vector v = *faces, *_i19 = faces; ((scalar *)&v)->i >= 0; v = *++_i19) {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
  val(v.x,(i + 1)/2,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x);
     }
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_x (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(-1,0,).pid - 1) : (- cell.pid - 1);
       if (faces) for (vector v = *faces, *_i20 = faces; ((scalar *)&v)->i >= 0; v = *++_i20) {

  scalar vt = VT;



  val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x);
       }
     }
     else

       if (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21)
  val(v.x,0,0,0) = 0.;
   }

 }
#line 1397

 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(0,i,) && !(neighbor(0,i,).pid < 0))) {
     Point neighbor = neighborp(0,i,);
     if (faces) for (vector v = *faces, *_i19 = faces; ((scalar *)&v)->i >= 0; v = *++_i19) {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
  val(v.y,0,(i + 1)/2,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y);
     }
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_y (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(0,-1,).pid - 1) : (- cell.pid - 1);
       if (faces) for (vector v = *faces, *_i20 = faces; ((scalar *)&v)->i >= 0; v = *++_i20) {

  scalar vt = VT;



  val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y);
       }
     }
     else

       if (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21)
  val(v.y,0,0,0) = 0.;
   }

 }}
    }
  } } end_foreach_boundary_level(); }

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (vectors,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
  enable_fpe_for_mpi();
}



#undef VN
#undef VT
#define VN _attribute[s.i].v.x
#define VT _attribute[s.i].v.y

static double masked_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1449 "/home/jiarongw/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (!(cell.pid < 0) && val(s,0,0,0) != nodata)
      sum += val(s,0,0,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;
}


#line 1457

static double masked_average_x (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1459 "/home/jiarongw/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (child.x < 0 && (!(cell.pid < 0) || !(neighbor(1,0,).pid < 0)) &&
 val(s,1,0,0) != nodata)
      sum += val(s,1,0,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;
}
#line 1457

static double masked_average_y (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1459 "/home/jiarongw/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (child.y < 0 && (!(cell.pid < 0) || !(neighbor(0,1,).pid < 0)) &&
 val(s,0,1,0) != nodata)
      sum += val(s,0,1,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;
}

static void masked_boundary_restriction (const Boundary * b,
      scalar * list, int l)
{
  scalar * scalars = NULL;
  vector * faces = NULL;
  if (list) for (scalar s = *list, *_i22 = list; ((scalar *)&s)->i >= 0; s = *++_i22)
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i && _attribute[s.i].face)
 faces = vectors_add (faces, _attribute[s.i].v);
      else
 scalars = list_add (scalars, s);
    }

   { foreach_halo (restriction, l){

#line 1481 "/home/jiarongw/basilisk/src/grid/tree.h"
 {
    if (scalars) for (scalar s = *scalars, *_i23 = scalars; ((scalar *)&s)->i >= 0; s = *++_i23)
      val(s,0,0,0) = masked_average (parent, s);
    if (faces) for (vector v = *faces, *_i24 = faces; ((scalar *)&v)->i >= 0; v = *++_i24)
      {
#line 1485
 {
 double average = masked_average_x (parent, v.x);
 if ((neighbor(-1,0,).pid < 0))
   val(v.x,0,0,0) = average;
 if ((neighbor(1,0,).pid < 0))
   val(v.x,1,0,0) = average;
      }
#line 1485
 {
 double average = masked_average_y (parent, v.y);
 if ((neighbor(0,-1,).pid < 0))
   val(v.y,0,0,0) = average;
 if ((neighbor(0,1,).pid < 0))
   val(v.y,0,1,0) = average;
      }}
  } } end_foreach_halo(); }

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
}
#line 1521 "/home/jiarongw/basilisk/src/grid/tree.h"
static void free_cache (CacheLevel * c)
{
  for (int l = 0; l <= depth(); l++)
    pfree (c[l].p,__func__,__FILE__,__LINE__);
  pfree (c,__func__,__FILE__,__LINE__);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Tree * q = ((Tree *)grid);
  pfree (q->leaves.p,__func__,__FILE__,__LINE__);
  pfree (q->faces.p,__func__,__FILE__,__LINE__);
  pfree (q->vertices.p,__func__,__FILE__,__LINE__);
  pfree (q->refined.p,__func__,__FILE__,__LINE__);


  Layer * L = q->L[0];




  for (int i = Period.x*2; i < L->len - Period.x*2; i++) {
    for (int j = Period.y*2; j < L->len - Period.y*2; j++)
      pfree (L->m[i][j],__func__,__FILE__,__LINE__);
    pfree (L->m[i],__func__,__FILE__,__LINE__);
  }

  for (int l = 1; l <= depth(); l++) {
    Layer * L = q->L[l];
    for (int i = Period.x*2; i < L->len - Period.x*2; i++)
      pfree (L->m[i],__func__,__FILE__,__LINE__);
  }
#line 1577 "/home/jiarongw/basilisk/src/grid/tree.h"
  for (int l = 0; l <= depth(); l++)
    destroy_layer (q->L[l]);
  q->L = &(q->L[-1]);
  pfree (q->L,__func__,__FILE__,__LINE__);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->boundary);
  free_cache (q->restriction);
  pfree (q,__func__,__FILE__,__LINE__);
  grid = NULL;
}

static void refine_level (int depth);


void init_grid (int n)
{ trace ("init_grid", "/home/jiarongw/basilisk/src/grid/tree.h", 1593);

  assert (sizeof(Cell) % 8 == 0);

  free_grid();
  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (qstderr(), "tree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  Tree * q = ((Tree *) pcalloc (1, sizeof(Tree),__func__,__FILE__,__LINE__));
  grid = (Grid *) q;
  grid->depth = 0;


  q->L = ((Layer * *) pmalloc ((2)*sizeof(Layer *),__func__,__FILE__,__LINE__));

  q->L[0] = NULL; q->L = &(q->L[1]);

  Layer * L = new_layer (0);
  q->L[0] = L;
#line 1633 "/home/jiarongw/basilisk/src/grid/tree.h"
  for (int i = Period.x*2; i < L->len - Period.x*2; i++) {
    layer_add_row (L, i, 0);
    for (int j = Period.y*2; j < L->len - Period.y*2; j++)
      L->m[i][j] = (char *) pcalloc (1, sizeof(Cell) + datasize,__func__,__FILE__,__LINE__);
  }
  apply_periodic (q);
  CELL(L->m[2][2]).flags |= leaf;
  if (pid() == 0)
    CELL(L->m[2][2]).flags |= active;
  for (int k = - 2*(1 - Period.x); k <= 2*(1 - Period.x); k++)
    for (int l = -2*(1 - Period.y); l <= 2*(1 - Period.y); l++)
      CELL(L->m[2 +k][2 +l]).pid = (k < 0 ? -1 - left :
         k > 0 ? -1 - right :
         l > 0 ? -1 - top :
         l < 0 ? -1 - bottom :
         0);
  CELL(L->m[2][2]).pid = 0;
#line 1673 "/home/jiarongw/basilisk/src/grid/tree.h"
  q->active = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->prolongation = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->boundary = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->restriction = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->dirty = true;
  N = 1 << depth;
#if 1
  void mpi_boundary_new();
  mpi_boundary_new();
#endif

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
  b->level = box_boundary_level;
  b->restriction = masked_boundary_restriction;
  add_boundary (b);
  refine_level (depth);
  reset (all, 0.);
  { if (((Tree *)grid)->dirty) update_cache_f(); };
 end_trace("init_grid", "/home/jiarongw/basilisk/src/grid/tree.h", 1691); }


void check_two_one (void)
{
   { foreach_leaf(){

#line 1696 "/home/jiarongw/basilisk/src/grid/tree.h"

    if (level > 0)
      for (int k = -1; k <= 1; k++)
 for (int l = -1; l <= 1; l++) {

   int i = (point.i + 2)/2 + k;
   int j = (point.j + 2)/2 + l;
   double x = ((i - 2 + 0.5)*(1./(1 << point.level))*2. - 0.5);
   double y = ((j - 2 + 0.5)*(1./(1 << point.level))*2. - 0.5);
   if (x > -0.5 && x < 0.5 && y > -0.5 && y < 0.5 &&
       !(aparent(k,l,0).flags & active)) {
     FILE * fp = fopen("check_two_one_loc", "w");
     fprintf (fp,
       "# %d %d\n"
       "%g %g\n%g %g\n",
       k, l,
       (((point.i - 2) + 0.5)*(1./(1 << point.level)) - 0.5),
       (((point.j - 2) + 0.5)*(1./(1 << point.level)) - 0.5),
       x, y);
     fclose (fp);





     assert (false);
   }
 } } end_foreach_leaf(); }
}


struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = { .level = l };
    int n = 1 << point.level;
    point.i = (p.x - X0)/L0*n + 2;

    point.j = (p.y - Y0)/L0*n + 2;




    if (point.i >= 0 && point.i < n + 2*2

 && point.j >= 0 && point.j < n + 2*2




 ) {
      if (allocated(0,0,0) && is_local(cell) && is_leaf(cell))
 return point;
    }
    else
      break;
  }
  Point point = { .level = -1 };
  return point;
}



bool tree_is_full()
{
  { if (((Tree *)grid)->dirty) update_cache_f(); };
  return (grid->tn == 1L << grid->maxdepth*2);
}

#line 1 "grid/tree-common.h"
#line 1 "/home/jiarongw/basilisk/src/grid/tree-common.h"



#line 1 "grid/multigrid-common.h"
#line 1 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/jiarongw/basilisk/src/grid/events.h"







static void event_error (Event * ev, const char * s)
{
  fprintf (qstderr(), "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = ev->t = -1;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = i;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == -123456) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = ev->t = -1;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == 1234567890 || ev->t == 1234567890) {
 ev->i = 1234567890; ev->t = -1;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != -1)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

void event_register (Event event) {
  assert (Events);
  assert (!event.last);
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      assert (parent < 0);
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,__LINE__);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,__LINE__));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 131 "/home/jiarongw/basilisk/src/grid/events.h"
static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (iter == ev->i || fabs (t - ev->t) <= 1e-9) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == 1234567890 && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = 1234567890; tnext = HUGE;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if ((!cond || cond1) && (tnext != HUGE || inext != 1234567890)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != HUGE && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    assert (n < INT_MAX);
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt + 1e-9)
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
#line 2 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])

#undef VARIABLES
#define VARIABLES\
  double Delta = L0*(1./(1 << point.level));\
  \
    double Delta_x = Delta;\
    double Delta_y = Delta;\
\
  double x = (ig/2. + (point.i - 2) + 0.5)*Delta + X0; NOT_UNUSED(x);\
\
  double y = (jg/2. + (point.j - 2) + 0.5)*Delta + Y0;\
\
\
\
 NOT_UNUSED(y);\
\
\
\
  double z = 0.;\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  \
    NOT_UNUSED(Delta_x);\
    NOT_UNUSED(Delta_y);\
\
  ;\

#line 32


#line 1 "grid/fpe.h"
#line 1 "/home/jiarongw/basilisk/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static void gdb()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', qstderr());
    fflush (qstderr());
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  system (command);
}

static void caught_abort (int sig)
{
  fprintf (qstderr(), "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (qstderr(), "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (qstderr(), "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 35 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

#define end_foreach_face()

scalar new_scalar (const char * name)
{
  int nvar = datasize/sizeof(double);
  scalar s;
  for (s.i = 0; s.i < nvar; s.i++)
    if (!list_lookup (all, s)) {
      init_scalar (s, name);
      trash (((scalar []){s, {-1}}));
      all = list_append (all, s);
      return s;
    }


  assert (nvar < _NVARMAX);
  datasize += sizeof(double); nvar++;
  _attribute = (_Attributes *) prealloc (_attribute, (nvar)*sizeof(_Attributes),__func__,__FILE__,__LINE__);
  memset (&_attribute[nvar-1], 0, sizeof (_Attributes));
  s = (scalar){nvar - 1};
  realloc_scalar();
  init_scalar (s, name);
  trash (((scalar []){s, {-1}}));
  all = list_append (all, s);
  return s;
}

scalar new_vertex_scalar (const char * name)
{
  scalar s = new_scalar (name);
  {
#line 66

    _attribute[s.i].d.x = -1;
#line 66

    _attribute[s.i].d.y = -1;}
  return s;
}

static vector alloc_vector (const char * name)
{
  vector v;
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  {
#line 76
 {
    sprintf (cname, ext.x, name);
    v.x = new_scalar (cname);
  }
#line 76
 {
    sprintf (cname, ext.y, name);
    v.y = new_scalar (cname);
  }}
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_vector (name);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_vector (name);
  init_face_vector (v, NULL);
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
  {
#line 102
 {
    sprintf (cname, ext.x, name);
    t.x = new_vector (cname);
  }
#line 102
 {
    sprintf (cname, ext.y, name);
    t.y = new_vector (cname);
  }}
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  char cname[strlen(name) + 5];
  struct { char * x, * y, * z; } ext = {"%s.x.x", "%s.y.y", "%s.z.z"};
  tensor t;
  {
#line 115
 {
    sprintf (cname, ext.x, name);
    t.x.x = new_scalar(cname);
  }
#line 115
 {
    sprintf (cname, ext.y, name);
    t.y.y = new_scalar(cname);
  }}

    sprintf (cname, "%s.x.y", name);
    t.x.y = new_scalar(cname);
    t.y.x = t.x.y;
#line 135 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = (double *) prealloc (_constant, (nconst)*sizeof(double),__func__,__FILE__,__LINE__);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  {
#line 159

    init_const_scalar (v.x, name, *val++);
#line 159

    init_const_scalar (v.y, name, *val++);}
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  {
#line 166

    v.x.i = _NVARMAX + i++;
#line 166

    v.y.i = _NVARMAX + i++;}
  init_const_vector (v, name, val);
  return v;
}

void scalar_clone (scalar a, scalar b)
{
  char * name = _attribute[a.i].name;
  double (** boundary) (Point, Point, scalar) = _attribute[a.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar) =
    _attribute[a.i].boundary_homogeneous;
  _attribute[a.i] = _attribute[b.i];
  _attribute[a.i].name = name;
  _attribute[a.i].boundary = boundary;
  _attribute[a.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[a.i].boundary[i] = _attribute[b.i].boundary[i];
    _attribute[a.i].boundary_homogeneous[i] = _attribute[b.i].boundary_homogeneous[i];
  }
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  if (l) for (scalar s = *l, *_i25 = l; ((scalar *)&s)->i >= 0; s = *++_i25) {
    scalar c = new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }
  if (list) for (scalar s = *list, *_i26 = list; ((scalar *)&s)->i >= 0; s = *++_i26)
    {
#line 201

      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
#line 201

      if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  if (list) for (scalar f = *list, *_i27 = list; ((scalar *)&f)->i >= 0; f = *++_i27) {
    if (_attribute[f.i].delete)
      _attribute[f.i].delete (f);
    pfree (_attribute[f.i].name,__func__,__FILE__,__LINE__); _attribute[f.i].name = NULL;
    pfree (_attribute[f.i].boundary,__func__,__FILE__,__LINE__); _attribute[f.i].boundary = NULL;
    pfree (_attribute[f.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[f.i].boundary_homogeneous = NULL;
  }

  if (list == all) {
    all[0].i = -1;
    return;
  }

  trash (list);
  if (list) for (scalar f = *list, *_i28 = list; ((scalar *)&f)->i >= 0; f = *++_i28) {
    scalar * s = all;
    for (; s->i >= 0 && s->i != f.i; s++);
    if (s->i == f.i)
      for (; s->i >= 0; s++)
 s[0] = s[1];
  }
}

typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}

void free_solver()
{
  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,__LINE__);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree (_attribute,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;
  free_grid();
  qpclose_all();
#if TRACE
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
#if _CADNA
  cadna_end();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_flux) (vector *);


void boundary (scalar * list)
{ trace ("boundary", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 289);
  if (list == NULL)
    { ; end_trace("boundary", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 291);  return; }
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i29 = list; ((scalar *)&s)->i >= 0; s = *++_i29)
    if (!is_constant(s) && _attribute[s.i].face)
      listf = vectors_add (listf, _attribute[s.i].v);
  if (listf) {
    boundary_flux (listf);
    pfree (listf,__func__,__FILE__,__LINE__);
  }
  boundary_level (list, -1);
 end_trace("boundary", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 301); }

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_flux (vector * list)
{

}

static double symmetry (Point point, Point neighbor, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 314 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 319 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  return -val(s,0,0,0);
}

double (* default_scalar_bc[]) (Point, Point, scalar) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  pfree (_attribute[s.i].boundary,__func__,__FILE__,__LINE__);
  pfree (_attribute[s.i].boundary_homogeneous,__func__,__FILE__,__LINE__);

  _attribute[s.i] = (const _Attributes){0};
  _attribute[s.i].name = pname;

  _attribute[s.i].boundary = (double (**)(Point, Point, scalar))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*2 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
  {
#line 351
 {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  }
#line 351
 {
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  }}
  _attribute[s.i].face = false;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 368
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_scalar (v.x, cname);
    }
    else
      init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  }
#line 368
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_scalar (v.y, cname);
    }
    else
      init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  }}

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*2 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  {
#line 388
 {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  }
#line 388
 {
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  }}
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 400
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_vector (t.x, cname);
    }
    else
      init_vector (t.x, NULL);
  }
#line 400
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_vector (t.y, cname);
    }
    else
      init_vector (t.y, NULL);
  }}






    for (int b = 0; b < nboundary; b++) {
      _attribute[t.x.x.i].boundary[b] = _attribute[t.y.x.i].boundary[b] =
 _attribute[t.x.x.i].boundary_homogeneous[b] = _attribute[t.y.y.i].boundary_homogeneous[b] =
 b < 2*2 ? default_scalar_bc[b] : symmetry;
      _attribute[t.x.y.i].boundary[b] = _attribute[t.y.y.i].boundary[b] =
 _attribute[t.x.y.i].boundary_homogeneous[b] = _attribute[t.y.x.i].boundary_homogeneous[b] =
 b < 2*2 ? default_vector_bc[b] : antisymmetry;
    }



  return t;
}

void output_cells (FILE * fp)
{
   { foreach(){

#line 431 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
 {
    Delta /= 2.;



      fprintf (fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
        x - Delta, y - Delta,
        x - Delta, y + Delta,
        x + Delta, y + Delta,
        x + Delta, y - Delta,
        x - Delta, y - Delta);
#line 456 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
  } } end_foreach(); }
  fflush (fp);
}

static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,__LINE__), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"




    "  plot '%s' w l lc 0, "
    "'%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 1 title columnhead(3+3*v)",





    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,__LINE__);
}

void cartesian_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 494 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells (fp);
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  if (all) for (scalar v = *all, *_i30 = all; ((scalar *)&v)->i >= 0; v = *++_i30)



    fprintf (fp, "x y %s ", _attribute[v.i].name);



  fputc ('\n', fp);
#line 527 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++) {
 if (all) for (scalar v = *all, *_i31 = all; ((scalar *)&v)->i >= 0; v = *++_i31) {
   fprintf (fp, "%g %g ",
     x + k*Delta + _attribute[v.i].d.x*Delta/2.,
     y + l*Delta + _attribute[v.i].d.y*Delta/2.);
   if (allocated(k,l,0))
     fprintf (fp, "%g ", val(v,k,l,0));
   else
     fputs ("n/a ", fp);
 }
 fputc ('\n', fp);
      }
#line 557 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  if (all) for (scalar s = *all, *_i32 = all; ((scalar *)&s)->i >= 0; s = *++_i32) {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }
  fclose (fp);

  fprintf (qstderr(), "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (qstderr(), _attribute[0].name, name, stencil);
  fflush (qstderr());

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar = cartesian_init_scalar;
  init_vector = cartesian_init_vector;
  init_tensor = cartesian_init_tensor;
  init_face_vector = cartesian_init_face_vector;
  boundary_level = cartesian_boundary_level;
  boundary_flux = cartesian_boundary_flux;
  debug = cartesian_debug;
}

struct _interpolate {
  scalar v;
  double x, y, z;
};

static double interpolate_linear (Point point, struct _interpolate p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 597 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  scalar v = p.v;







  x = (p.x - x)/Delta - _attribute[v.i].d.x/2.;
  y = (p.y - y)/Delta - _attribute[v.i].d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);

  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);
#line 625 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
}


double interpolate (struct _interpolate p)
{ trace ("interpolate", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 629);
  Point point = locate ((struct _locate){p.x, p.y, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 630 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  if (point.level < 0)
    { double _ret =  nodata; end_trace("interpolate", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 632);  return _ret; }
  { double _ret =  interpolate_linear (point, p); end_trace("interpolate", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 633);  return _ret; }
 end_trace("interpolate", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 634); }


void interpolate_array (scalar * list, coord * a, int n, double * v, bool linear)
{ trace ("interpolate_array", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 638);
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point point = locate ((struct _locate){a[i].x, a[i].y, a[i].z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 641 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

    if (point.level >= 0) {
      if (list) for (scalar s = *list, *_i33 = list; ((scalar *)&s)->i >= 0; s = *++_i33)
 v[j++] = !linear ? val(s,0,0,0) :
   interpolate_linear (point,
         (struct _interpolate){s, a[i].x, a[i].y, a[i].z});
    }
    else
      if (list) for (scalar s = *list, *_i34 = list; ((scalar *)&s)->i >= 0; s = *++_i34)
 v[j++] = nodata;
  }
#if 1
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (v, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
#endif
 end_trace("interpolate_array", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 660); }



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  if (all) for (scalar s = *all, *_i35 = all; ((scalar *)&s)->i >= 0; s = *++_i35) {
    _attribute[s.i].boundary = (double (**)(Point, Point, scalar))
      prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar))
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  }
  if (all) for (scalar s = *all, *_i36 = all; ((scalar *)&s)->i >= 0; s = *++_i36) {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      {
#line 680

 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
#line 680

 _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;}
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 692 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  return nodata;
}

static void periodic_boundary (int d)
{

  if (all) for (scalar s = *all, *_i37 = all; ((scalar *)&s)->i >= 0; s = *++_i37)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;

  if (all) for (scalar s = *all, *_i38 = all; ((scalar *)&s)->i >= 0; s = *++_i38)
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{



    assert (dir <= bottom);




  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}
#line 4 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf foreach_level
# define end_foreach_level_or_leaf end_foreach_level
#endif

#ifndef foreach_coarse_level
# define foreach_coarse_level foreach_level
# define end_foreach_coarse_level end_foreach_level
#endif










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 27 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  double sum = 0.;
   { foreach_child()
    sum += val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 2);
}

static inline void restriction_volume_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 35 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 2)/val_cm(cm,0,0,0);
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 2)/val_cm(cm,0,0,0);
 }}

static inline void face_average (Point point, vector v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 43 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  {
#line 44
 {




      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0))/2.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0))/2.;






  }
#line 44
 {




      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,1,0,0))/2.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,1,2,0))/2.;






  }}
}

static inline void restriction_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 61 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  face_average (point, _attribute[s.i].v);
}

static inline void no_restriction (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 64 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
}

static inline void no_data (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 67 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = nodata; end_foreach_child(); }
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar []){s,{-1}}));
  for (int l = depth() - 1; l >= 0; l--)
     { foreach_coarse_level (l){

#line 76 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
 {
      double sc[1 << 2];
      int c = 0;
       { foreach_child()
 sc[c++] = val(s,0,0,0); end_foreach_child(); }
      _attribute[s.i].prolongation (point, s);
      c = 0;
       { foreach_child() {

 val(w,0,0,0) = sc[c] - val(s,0,0,0);
 val(s,0,0,0) = sc[c++];
      } end_foreach_child(); }
    } } end_foreach_coarse_level(); }

   { foreach_level(0){

#line 90 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
 val(w,0,0,0) = 0.; } end_foreach_level(); }
}

static inline double bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 94 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"




    return (9.*coarse(s,0,0,0) +
     3.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0)) +
     coarse(s,child.x,child.y,0))/16.;
#line 109 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
}

static inline void refine_bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 112 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = bilinear (point, s); end_foreach_child(); }
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 123 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"




  return
    quadratic (quadratic (coarse(s,0,0,0),
     coarse(s,child.x,0,0),
     coarse(s,-child.x,0,0)),
        quadratic (coarse(s,0,child.y,0),
     coarse(s,child.x,child.y,0),
     coarse(s,-child.x,child.y,0)),
        quadratic (coarse(s,0,-child.y,0),
     coarse(s,child.x,-child.y,0),
     coarse(s,-child.x,-child.y,0)));




}

static inline double biquadratic_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 144 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"




  return (36.*val(s,0,0,0) + 18.*(val(s,-1,0,0) + val(s,0,-1,0)) - 6.*(val(s,1,0,0) + val(s,0,1,0)) +
   9.*val(s,-1,-1,0) - 3.*(val(s,1,-1,0) + val(s,-1,1,0)) + val(s,1,1,0))/64.;




}

static inline void refine_biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 157 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = biquadratic (point, s); end_foreach_child(); }
}

static inline void refine_linear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 163 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 163

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 166

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 166

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
#line 169

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 169

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 2);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 175

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 175

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  assert (fabs(sum) < 1e-10);
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 163

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 166

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 166

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
#line 169

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 169

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 2);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 175

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 175

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  assert (fabs(sum) < 1e-10);
 }}

static inline void refine_reset (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 183 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(v,0,0,0) = 0.; end_foreach_child(); }
}

static inline void refine_injection (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 189 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  double val = val(v,0,0,0);
   { foreach_child()
    val(v,0,0,0) = val; end_foreach_child(); }
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  {
#line 206

    _attribute[v.y.i].restriction = no_restriction;
#line 206

    _attribute[v.x.i].restriction = no_restriction;}
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

void multigrid_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 213 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 233 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++) {
   if (all) for (scalar v = *all, *_i39 = all; ((scalar *)&v)->i >= 0; v = *++_i39)
     fprintf (fp, "%g %g %g ",
       xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
       yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
       coarse(v,k*child.x,l*child.y,0));
   fputc ('\n', fp);
 }
      fprintf (qstderr(), ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
#line 264 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 286 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++) {
   if (all) for (scalar v = *all, *_i40 = all; ((scalar *)&v)->i >= 0; v = *++_i40) {
     fprintf (fp, "%g %g ",
       xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
       yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.);
     if (allocated_child(k,l,0))
       fprintf (fp, "%g ", fine(v,k,l,0));
     else
       fputs ("n/a ", fp);
   }
   fputc ('\n', fp);
 }
      fprintf (qstderr(), ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
#line 324 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }
  fflush (qstderr());
  fclose (plot);
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  if (list) for (scalar s = *list, *_i41 = list; ((scalar *)&s)->i >= 0; s = *++_i41)
    if (!is_constant (s)) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
#line 342

     list2 = list_add (list2, _attribute[s.i].v.x);
#line 342

     list2 = list_add (list2, _attribute[s.i].v.y);}
 else
   list2 = list_add (list2, s);
      }
    }

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
       { foreach_coarse_level(l){

#line 351 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
 {
 if (listdef) for (scalar s = *listdef, *_i42 = listdef; ((scalar *)&s)->i >= 0; s = *++_i42)
   restriction_average (point, s);
 if (listc) for (scalar s = *listc, *_i43 = listc; ((scalar *)&s)->i >= 0; s = *++_i43)
   _attribute[s.i].restriction (point, s);
      } } end_foreach_coarse_level(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  debug = multigrid_debug;
  init_scalar = multigrid_init_scalar;
  init_face_vector = multigrid_init_face_vector;
  restriction = multigrid_restriction;
}







void subtree_size (scalar size, bool leaves)
{




   { foreach(){

#line 386 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

    val(size,0,0,0) = 1; } end_foreach(); }





  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
     { foreach_coarse_level(l){

#line 395 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
 {
      double sum = !leaves;
       { foreach_child()
 sum += val(size,0,0,0); end_foreach_child(); }
      val(size,0,0,0) = sum;
    } } end_foreach_coarse_level(); }
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), l); };
  }
}
#line 5 "/home/jiarongw/basilisk/src/grid/tree-common.h"






#line 21 "/home/jiarongw/basilisk/src/grid/tree-common.h"
int refine_cell (Point point, scalar * list, int flag, Cache * refined)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 22 "/home/jiarongw/basilisk/src/grid/tree-common.h"

  int nr = 0;


  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)

      for (int l = 0; l != 2*child.y; l += child.y)




   if (aparent(k,l,m).pid >= 0 && is_leaf(aparent(k,l,m))) {
     Point p = point;


     p.level = point.level - 1;
     p.i = (point.i + 2)/2 + k;
     do { if (p.i < 2) p.i += 1 << p.level; else if (p.i >= 2 + (1 << p.level)) p.i -= 1 << p.level; } while(0);

       p.j = (point.j + 2)/2 + l;
       do { if (p.j < 2) p.j += 1 << p.level; else if (p.j >= 2 + (1 << p.level)) p.j -= 1 << p.level; } while(0);





     nr += refine_cell (p, list, flag, refined);
     aparent(k,l,m).flags |= flag;
   }



  cell.flags &= ~leaf;


  increment_neighbors (point);

  int cflag = is_active(cell) ? (active|leaf) : leaf;
   { foreach_child()
    cell.flags |= cflag; end_foreach_child(); }


  if (list) for (scalar s = *list, *_i44 = list; ((scalar *)&s)->i >= 0; s = *++_i44)
    if (is_local(cell) || _attribute[s.i].face)
      _attribute[s.i].refine (point, s);

#if 1
  if (is_border(cell)) {
     { foreach_child() {
      bool bord = false;
       { foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0))))
   bord = true, foreach_neighbor_break();
 if (is_refined_check())
    { foreach_child()
     if (!is_local(cell))
       bord = true, foreach_child_break(); end_foreach_child(); }
 if (bord)
   foreach_neighbor_break();
      } end_foreach_neighbor(); }
      if (bord)
 cell.flags |= border;
    } end_foreach_child(); }
    if (refined)
      cache_append (refined, point, cell.flags);
    nr++;
  }
#endif
  return nr;
}





bool coarsen_cell (Point point, scalar * list)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 99 "/home/jiarongw/basilisk/src/grid/tree-common.h"




  int pid = cell.pid;
   { foreach_child()
    if (cell.neighbors || (cell.pid < 0 && cell.pid != pid))
      return false; end_foreach_child(); }



  if (list) for (scalar s = *list, *_i45 = list; ((scalar *)&s)->i >= 0; s = *++_i45) {
    _attribute[s.i].restriction (point, s);
    if (_attribute[s.i].coarsen)
      _attribute[s.i].coarsen (point, s);
  }


  cell.flags |= leaf;


  decrement_neighbors (point);

#if 1
  if (!is_local(cell)) {
    cell.flags &= ~(active|border);
     { foreach_neighbor(1)
      if (cell.neighbors)
  { foreach_child()
   if (allocated(0,0,0) && is_local(cell) && !is_border(cell))
     cell.flags |= border; end_foreach_child(); } end_foreach_neighbor(); }
  }
#endif

  return true;
}

void coarsen_cell_recursive (Point point, scalar * list)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 137 "/home/jiarongw/basilisk/src/grid/tree-common.h"



   { foreach_child()
    if (cell.neighbors)
       { foreach_neighbor(1)
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   coarsen_cell_recursive (point, list); end_foreach_neighbor(); } end_foreach_child(); }

  assert (coarsen_cell (point, list));
}

void mpi_boundary_refine (scalar *);
void mpi_boundary_coarsen (int, int);
void mpi_boundary_update (scalar *);

typedef struct {
  int nc, nf;
} astats;

struct Adapt {
  scalar * slist;
  double * max;
  int maxlevel;
  int minlevel;
  scalar * list;
};


astats adapt_wavelet (struct Adapt p)
{ trace ("adapt_wavelet", "/home/jiarongw/basilisk/src/grid/tree-common.h", 167);
  if (p.list == NULL)
    p.list = all;
  if (is_constant(cm))
    restriction (p.slist);
  else {
    scalar * listr = list_concat (((scalar []){cm,{-1}}), p.slist);
    restriction (listr);
    pfree (listr,__func__,__FILE__,__LINE__);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  if (p.list) for (scalar s = *p.list, *_i46 = p.list; ((scalar *)&s)->i >= 0; s = *++_i46)
    if (!is_constant(s) && _attribute[s.i].restriction != no_restriction)
      listc = list_add (listc, s);


  if (p.minlevel < 1)
    p.minlevel = 1;
  ((Tree *)grid)->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
   { foreach_cell(){

#line 189 "/home/jiarongw/basilisk/src/grid/tree-common.h"
 {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
 if (cell.flags & too_coarse) {
   cell.flags &= ~too_coarse;
   refine_cell (point, listc, refined, &((Tree *)grid)->refined);
   st.nf++;
 }
 continue;
      }
      else {
 if (cell.flags & refined) {

   cell.flags &= ~too_coarse;
   continue;
 }

 bool local = is_local(cell);
 if (!local)
    { foreach_child()
     if (is_local(cell))
       local = true, foreach_child_break(); end_foreach_child(); }
 if (local) {
   int i = 0;
   static const int just_fine = 1 << (user + 3);
   if (p.slist) for (scalar s = *p.slist, *_i47 = p.slist; ((scalar *)&s)->i >= 0; s = *++_i47) {
     double max = p.max[i++], sc[1 << 2];
     int c = 0;
      { foreach_child()
       sc[c++] = val(s,0,0,0); end_foreach_child(); }
     _attribute[s.i].prolongation (point, s);
     c = 0;
      { foreach_child() {
       double e = fabs(sc[c] - val(s,0,0,0));
       if (e > max && level < p.maxlevel) {
  cell.flags &= ~too_fine;
  cell.flags |= too_coarse;
       }
       else if ((e <= max/1.5 || level > p.maxlevel) &&
         !(cell.flags & (too_coarse|just_fine))) {
  if (level >= p.minlevel)
    cell.flags |= too_fine;
       }
       else if (!(cell.flags & too_coarse)) {
  cell.flags &= ~too_fine;
  cell.flags |= just_fine;
       }
       val(s,0,0,0) = sc[c++];
     } end_foreach_child(); }
   }
    { foreach_child() {
     cell.flags &= ~just_fine;
     if (!is_leaf(cell)) {
       cell.flags &= ~too_coarse;
       if (level >= p.maxlevel)
  cell.flags |= too_fine;
     }
     else if (!is_active(cell))
       cell.flags &= ~too_coarse;
   } end_foreach_child(); }
 }
      }
    }
    else
      continue;
  } } end_foreach_cell(); }
  mpi_boundary_refine (listc);



  for (int l = depth(); l >= 0; l--) {
     { foreach_cell(){

#line 261 "/home/jiarongw/basilisk/src/grid/tree-common.h"

      if (!(cell.pid < 0)) {
 if (level == l) {
   if (!is_leaf(cell)) {
     if (cell.flags & refined)

       cell.flags &= ~(refined|too_fine);
     else if (cell.flags & too_fine) {
       if (is_local(cell) && coarsen_cell (point, listc))
  st.nc++;
       cell.flags &= ~too_fine;
     }
   }
   if (cell.flags & too_fine)
     cell.flags &= ~too_fine;
   else if (level > 0 && (aparent(0,0,0).flags & too_fine))
     aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
 else if (is_leaf(cell))
   continue;
      } } end_foreach_cell(); }
    mpi_boundary_coarsen (l, too_fine);
  }
  pfree (listc,__func__,__FILE__,__LINE__);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (p.list);

  { astats _ret =  st; end_trace("adapt_wavelet", "/home/jiarongw/basilisk/src/grid/tree-common.h", 292);  return _ret; }
 end_trace("adapt_wavelet", "/home/jiarongw/basilisk/src/grid/tree-common.h", 293); }
#line 314 "/home/jiarongw/basilisk/src/grid/tree-common.h"
static void refine_level (int depth)
{
  int refined;
  do {
    refined = 0;
    ((Tree *)grid)->refined.n = 0;
     { foreach_leaf(){

#line 320 "/home/jiarongw/basilisk/src/grid/tree-common.h"

      if (level < depth) {
 refine_cell (point, NULL, 0, &((Tree *)grid)->refined);
 refined++;
 continue;
      } } end_foreach_leaf(); }
    mpi_all_reduce (refined, MPI_INT, MPI_SUM);
    if (refined) {
      mpi_boundary_refine (NULL);
      mpi_boundary_update (NULL);
    }
  } while (refined);
}
#line 359 "/home/jiarongw/basilisk/src/grid/tree-common.h"
static void halo_flux (vector * list)
{
  vector * listv = NULL;
  if (list) for (vector v = *list, *_i48 = list; ((scalar *)&v)->i >= 0; v = *++_i48)
    if (!is_constant(v.x))
      listv = vectors_add (listv, v);

  if (listv) {
    for (int l = depth() - 1; l >= 0; l--)
       { foreach_halo (prolongation, l){

#line 368 "/home/jiarongw/basilisk/src/grid/tree-common.h"

 {
#line 369
 {
#line 378 "/home/jiarongw/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(-1,0,)) && neighbor(-1,0,).neighbors && neighbor(-1,0,).pid >= 0))
     if (listv) for (vector f = *listv, *_i49 = listv; ((scalar *)&f)->i >= 0; f = *++_i49)
       val(f.x,0,0,0) = (fine(f.x,0,0,0) + fine(f.x,0,1,0))/2.;
   if ((!is_leaf (neighbor(1,0,)) && neighbor(1,0,).neighbors && neighbor(1,0,).pid >= 0))
     if (listv) for (vector f = *listv, *_i50 = listv; ((scalar *)&f)->i >= 0; f = *++_i50)
       val(f.x,1,0,0) = (fine(f.x,2,0,0) + fine(f.x,2,1,0))/2.;
#line 394 "/home/jiarongw/basilisk/src/grid/tree-common.h"
      }
#line 369
 {
#line 378 "/home/jiarongw/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(0,-1,)) && neighbor(0,-1,).neighbors && neighbor(0,-1,).pid >= 0))
     if (listv) for (vector f = *listv, *_i49 = listv; ((scalar *)&f)->i >= 0; f = *++_i49)
       val(f.y,0,0,0) = (fine(f.y,0,0,0) + fine(f.y,1,0,0))/2.;
   if ((!is_leaf (neighbor(0,1,)) && neighbor(0,1,).neighbors && neighbor(0,1,).pid >= 0))
     if (listv) for (vector f = *listv, *_i50 = listv; ((scalar *)&f)->i >= 0; f = *++_i50)
       val(f.y,0,1,0) = (fine(f.y,0,2,0) + fine(f.y,1,2,0))/2.;
#line 394 "/home/jiarongw/basilisk/src/grid/tree-common.h"
      }} } end_foreach_halo(); }
    pfree (listv,__func__,__FILE__,__LINE__);
  }
}



static scalar tree_init_scalar (scalar s, const char * name)
{
  s = multigrid_init_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation;
  return s;
}


#line 408

static void refine_face_x (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 410 "/home/jiarongw/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;
  if (!(!is_leaf (neighbor(-1,0,)) && neighbor(-1,0,).neighbors && neighbor(-1,0,).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(-1,0,)))) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0))/8.;
    double g2 = (val(v.x,0,0,+1) - val(v.x,0,0,-1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,0,j,k) = val(v.x,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(1,0,)) && neighbor(1,0,).neighbors && neighbor(1,0,).pid >= 0) && neighbor(1,0,).neighbors &&
      (is_local(cell) || is_local(neighbor(1,0,)))) {
    double g1 = (val(v.x,1,+1,0) - val(v.x,1,-1,0))/8.;
    double g2 = (val(v.x,1,0,+1) - val(v.x,1,0,-1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,2,j,k) = val(v.x,1,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.x,0,+1,0) + val(v.x,1,+1,0) - val(v.x,0,-1,0) - val(v.x,1,-1,0))/16.;
    double g2 = (val(v.x,0,0,+1) + val(v.x,1,0,+1) - val(v.x,0,0,-1) - val(v.x,1,0,-1))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,1,j,k) = (val(v.x,0,0,0) + val(v.x,1,0,0))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
}
#line 408

static void refine_face_y (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 410 "/home/jiarongw/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;
  if (!(!is_leaf (neighbor(0,-1,)) && neighbor(0,-1,).neighbors && neighbor(0,-1,).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,-1,)))) {
    double g1 = (val(v.y,+1,0,0) - val(v.y,-1,0,0))/8.;
    double g2 = (val(v.y,0,0,+1) - val(v.y,0,0,-1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,j,0,k) = val(v.y,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(0,1,)) && neighbor(0,1,).neighbors && neighbor(0,1,).pid >= 0) && neighbor(0,1,).neighbors &&
      (is_local(cell) || is_local(neighbor(0,1,)))) {
    double g1 = (val(v.y,+1,1,0) - val(v.y,-1,1,0))/8.;
    double g2 = (val(v.y,0,1,+1) - val(v.y,0,1,-1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,j,2,k) = val(v.y,0,1,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.y,+1,0,0) + val(v.y,+1,1,0) - val(v.y,-1,0,0) - val(v.y,-1,1,0))/16.;
    double g2 = (val(v.y,0,0,+1) + val(v.y,0,1,+1) - val(v.y,0,0,-1) - val(v.y,0,1,-1))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,j,1,k) = (val(v.y,0,0,0) + val(v.y,0,1,0))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
}

void refine_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 438 "/home/jiarongw/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;
  {
#line 440

    _attribute[v.x.i].prolongation (point, v.x);
#line 440

    _attribute[v.y.i].prolongation (point, v.y);}
}

void refine_face_solenoidal (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 445 "/home/jiarongw/basilisk/src/grid/tree-common.h"

  refine_face (point, s);

  if (is_local(cell)) {

    vector v = _attribute[s.i].v;
    double d[1 << 2], p[1 << 2];
    int i = 0;
     { foreach_child() {
      d[i] = 0.;
      {
#line 455

 d[i] += val(v.x,1,0,0) - val(v.x,0,0,0);
#line 455

 d[i] += val(v.y,0,1,0) - val(v.y,0,0,0);}
      i++;
    } end_foreach_child(); }

    p[0] = 0.;
    p[1] = (3.*d[3] + d[0])/4. + d[2]/2.;
    p[2] = (d[3] + 3.*d[0])/4. + d[2]/2.;
    p[3] = (d[3] + d[0])/2. + d[2];
    fine(v.x,1,1,0) += p[1] - p[0];
    fine(v.x,1,0,0) += p[3] - p[2];
    fine(v.y,0,1,0) += p[0] - p[2];
    fine(v.y,1,1,0) += p[1] - p[3];
#line 495 "/home/jiarongw/basilisk/src/grid/tree-common.h"
  }

}

vector tree_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  {
#line 502

    _attribute[v.x.i].restriction = _attribute[v.x.i].refine = no_restriction;
#line 502

    _attribute[v.y.i].restriction = _attribute[v.y.i].refine = no_restriction;}
  _attribute[v.x.i].restriction = restriction_face;
  _attribute[v.x.i].refine = refine_face;
  {
#line 506

    _attribute[v.x.i].prolongation = refine_face_x;
#line 506

    _attribute[v.y.i].prolongation = refine_face_y;}
  return v;
}

static void tree_boundary_level (scalar * list, int l)
{
  int depth = l < 0 ? depth() : l;

  if (tree_is_full()) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, depth); };
    return;
  }

  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  if (list) for (scalar s = *list, *_i51 = list; ((scalar *)&s)->i >= 0; s = *++_i51)
    if (!is_constant (s)) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
#line 530

     list2 = list_add (list2, _attribute[s.i].v.x);
#line 530

     list2 = list_add (list2, _attribute[s.i].v.y);}
 else
   list2 = list_add (list2, s);
      }
    }

  if (listdef || listc) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, depth); };
    for (int l = depth - 1; l >= 0; l--) {
       { foreach_coarse_level(l){

#line 540 "/home/jiarongw/basilisk/src/grid/tree-common.h"
 {
 if (listdef) for (scalar s = *listdef, *_i52 = listdef; ((scalar *)&s)->i >= 0; s = *++_i52)
   restriction_average (point, s);
 if (listc) for (scalar s = *listc, *_i53 = listc; ((scalar *)&s)->i >= 0; s = *++_i53)
   _attribute[s.i].restriction (point, s);
      } } end_foreach_coarse_level(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }

  scalar * listr = NULL;
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i54 = list; ((scalar *)&s)->i >= 0; s = *++_i54)
    if (!is_constant (s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listr = list_add (listr, s);
    }

  if (listr || listf) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, 0); };
    for (int i = 0; i < depth; i++) {
       { foreach_halo (prolongation, i){

#line 566 "/home/jiarongw/basilisk/src/grid/tree-common.h"
 {
 if (listr) for (scalar s = *listr, *_i55 = listr; ((scalar *)&s)->i >= 0; s = *++_i55)
          _attribute[s.i].prolongation (point, s);
 if (listf) for (vector v = *listf, *_i56 = listf; ((scalar *)&v)->i >= 0; v = *++_i56)
   {
#line 570

     _attribute[v.x.i].prolongation (point, v.x);
#line 570

     _attribute[v.y.i].prolongation (point, v.y);}
      } } end_foreach_halo(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, i + 1); };
    }
    pfree (listr,__func__,__FILE__,__LINE__);
    pfree (listf,__func__,__FILE__,__LINE__);
  }
}

double treex (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 580 "/home/jiarongw/basilisk/src/grid/tree-common.h"

  if (level == 0)
    return 0;

  double i = 2*child.x - child.y;
  if (i <= 1 && i >= -1) i = -i;




  return treex(parent) + i/(1 << 2*(level - 1));
}

double treey (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 593 "/home/jiarongw/basilisk/src/grid/tree-common.h"

  if (level == 0)
    return 0;
  return treey(parent) + 4./(1 << 2*(level - 1));
}

void output_tree (FILE * fp)
{
   { foreach_cell(){

#line 601 "/home/jiarongw/basilisk/src/grid/tree-common.h"

    if (cell.neighbors)
       { foreach_child()
 if (is_local(cell))
   fprintf (fp, "%g %g\n%g %g\n\n",
     treex(parent), treey(parent), treex(point), treey(point)); end_foreach_child(); }; } end_foreach_cell(); }
}

void tree_check()
{


  long nleaves = 0, nactive = 0;
   { foreach_cell_all(){

#line 614 "/home/jiarongw/basilisk/src/grid/tree-common.h"
 {
    if (is_leaf(cell)) {
      assert (cell.pid >= 0);
      nleaves++;
    }
    if (is_local(cell))
      assert (is_active(cell) || (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0));
    if (is_active(cell))
      nactive++;

    int neighbors = 0;
     { foreach_neighbor(1)
      if (allocated(0,0,0) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 neighbors++; end_foreach_neighbor(); }
    assert (cell.neighbors == neighbors);


    if (!cell.neighbors)
      assert (!allocated_child(0,0,0));
  } } end_foreach_cell_all(); }


  long reachable = 0;
   { foreach_cell(){

#line 637 "/home/jiarongw/basilisk/src/grid/tree-common.h"
 {
    if (is_active(cell))
      reachable++;
    else
      continue;
  } } end_foreach_cell(); }
  assert (nactive == reachable);


  reachable = 0;
   { foreach_cell(){

#line 647 "/home/jiarongw/basilisk/src/grid/tree-common.h"

    if (is_leaf(cell)) {
      reachable++;
      continue;
    } } end_foreach_cell(); }
  assert (nleaves == reachable);
}

static void tree_restriction (scalar * list) {
  if (tree_is_full())
    multigrid_restriction (list);

}

void tree_methods()
{
  multigrid_methods();
  init_scalar = tree_init_scalar;
  init_face_vector = tree_init_face_vector;
  boundary_level = tree_boundary_level;
  boundary_flux = halo_flux;
  restriction = tree_restriction;
}
#line 1768 "/home/jiarongw/basilisk/src/grid/tree.h"


void tree_periodic (int dir)
{
  int depth = grid ? depth() : -1;
  if (grid)
    free_grid();
  periodic (dir);
  if (depth >= 0)
    init_grid (1 << depth);
}


#if 1
#line 1 "grid/tree-mpi.h"
#line 1 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

int debug_iteration = -1;

void debug_mpi (FILE * fp1);

typedef struct {
  CacheLevel * halo;
  void * buf;
  MPI_Request r;
  int depth;
  int pid;
  int maxdepth;
} Rcv;

typedef struct {
  Rcv * rcv;
  char * name;
  int npid;
} RcvPid;

typedef struct {
  RcvPid * rcv, * snd;
} SndRcv;

typedef struct {
  Boundary parent;

  SndRcv mpi_level, mpi_level_root, restriction;
  Array * send, * receive;
} MpiBoundary;

static void cache_level_init (CacheLevel * c)
{
  c->p = NULL;
  c->n = c->nm = 0;
}

static void rcv_append (Point point, Rcv * rcv)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 39 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

  if (level > rcv->depth) {
    rcv->halo = (CacheLevel *) prealloc (rcv->halo, (level + 1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__);
    for (int j = rcv->depth + 1; j <= level; j++)
      cache_level_init (&rcv->halo[j]);
    rcv->depth = level;
  }
  cache_level_append (&rcv->halo[level], point);
  if (level > rcv->maxdepth)
    rcv->maxdepth = level;
}

void rcv_print (Rcv * rcv, FILE * fp, const char * prefix)
{
  for (int l = 0; l <= rcv->depth; l++)
    if (rcv->halo[l].n > 0)
       { foreach_cache_level(rcv->halo[l], l){

#line 55 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

 fprintf (fp, "%s%g %g %g %d %d\n", prefix, x, y, z, rcv->pid, level); } end_foreach_cache_level(); }
}

static void rcv_free_buf (Rcv * rcv)
{
  if (rcv->buf) {
    prof_start ("rcv_pid_receive");
    MPI_Wait (&rcv->r, MPI_STATUS_IGNORE);
    pfree (rcv->buf,__func__,__FILE__,__LINE__);
    rcv->buf = NULL;
    prof_stop();
  }
}

static void rcv_destroy (Rcv * rcv)
{
  rcv_free_buf (rcv);
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      pfree (rcv->halo[i].p,__func__,__FILE__,__LINE__);
  pfree (rcv->halo,__func__,__FILE__,__LINE__);
}

static RcvPid * rcv_pid_new (const char * name)
{
  RcvPid * r = ((RcvPid *) pcalloc (1, sizeof(RcvPid),__func__,__FILE__,__LINE__));
  r->name = pstrdup (name,__func__,__FILE__,__LINE__);
  return r;
}

static Rcv * rcv_pid_pointer (RcvPid * p, int pid)
{
  assert (pid >= 0 && pid < npe());

  int i;
  for (i = 0; i < p->npid; i++)
    if (pid == p->rcv[i].pid)
      break;

  if (i == p->npid) {
    p->rcv = (Rcv *) prealloc (p->rcv, (++p->npid)*sizeof(Rcv),__func__,__FILE__,__LINE__);
    Rcv * rcv = &p->rcv[p->npid-1];
    rcv->pid = pid;
    rcv->depth = rcv->maxdepth = 0;
    rcv->halo = ((CacheLevel *) pmalloc ((1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__));
    rcv->buf = NULL;
    cache_level_init (&rcv->halo[0]);
  }
  return &p->rcv[i];
}

static void rcv_pid_append (RcvPid * p, int pid, Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 108 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

  rcv_append (point, rcv_pid_pointer (p, pid));
}

static void rcv_pid_append_pids (RcvPid * p, Array * pids)
{

  for (int i = 0; i < p->npid; i++) {
    int pid = p->rcv[i].pid, j, * a;
    for (j = 0, a = pids->p; j < pids->len/sizeof(int); j++,a++)
      if (*a == pid)
 break;
    if (j == pids->len/sizeof(int))
      array_append (pids, &pid, sizeof(int));
  }
}

void rcv_pid_write (RcvPid * p, const char * name)
{
  for (int i = 0; i < p->npid; i++) {
    Rcv * rcv = &p->rcv[i];
    char fname[80];
    sprintf (fname, "%s-%d-%d", name, pid(), rcv->pid);
    FILE * fp = fopen (fname, "w");
    rcv_print (rcv, fp, "");
    fclose (fp);
  }
}

static void rcv_pid_print (RcvPid * p, FILE * fp, const char * prefix)
{
  for (int i = 0; i < p->npid; i++)
    rcv_print (&p->rcv[i], fp, prefix);
}

static void rcv_pid_destroy (RcvPid * p)
{
  for (int i = 0; i < p->npid; i++)
    rcv_destroy (&p->rcv[i]);
  pfree (p->rcv,__func__,__FILE__,__LINE__);
  pfree (p->name,__func__,__FILE__,__LINE__);
  pfree (p,__func__,__FILE__,__LINE__);
}

static Boundary * mpi_boundary = NULL;






void debug_mpi (FILE * fp1);

static void apply_bc (Rcv * rcv, scalar * list, vector * listf, int l,
        MPI_Status s)
{
  double * b = rcv->buf;
   { foreach_cache_level(rcv->halo[l], l){

#line 165 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
    if (list) for (scalar s = *list, *_i57 = list; ((scalar *)&s)->i >= 0; s = *++_i57)
      val(s,0,0,0) = *b++;
    if (listf) for (vector v = *listf, *_i58 = listf; ((scalar *)&v)->i >= 0; v = *++_i58) {
      {
#line 169
 {
 val(v.x,0,0,0) = *b++;
 if (allocated(1,0,))
   val(v.x,1,0,0) = *b;
 b++;
      }
#line 169
 {
 val(v.y,0,0,0) = *b++;
 if (allocated(0,1,))
   val(v.y,0,1,0) = *b;
 b++;
      }}
    }
  } } end_foreach_cache_level(); }
  size_t size = b - (double *) rcv->buf;
  pfree (rcv->buf,__func__,__FILE__,__LINE__);
  rcv->buf = NULL;

  int rlen;
  MPI_Get_count (&s, MPI_DOUBLE, &rlen);
  if (rlen != size) {
    fprintf (qstderr(),
      "rlen (%d) != size (%ld), %d receiving from %d at level %d\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      rlen, size, pid(), rcv->pid, l);
    fflush (qstderr());
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -2);
  }
}
#line 215 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
static void mpi_recv_check (void * buf, int count, MPI_Datatype datatype,
       int source, int tag,
       MPI_Comm comm, MPI_Status * status,
       const char * name)
{
#line 250 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
  int errorcode = MPI_Recv (buf, count, datatype, source, tag, comm, status);
  if (errorcode != MPI_SUCCESS) {
    char string[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string (errorcode, string, &resultlen);
    fprintf (qstderr(),
      "ERROR MPI_Recv \"%s\" (count = %d, source = %d, tag = %d):\n%s\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      name, count, source, tag, string);
    fflush (qstderr());
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -1);
  }





}


static int mpi_waitany (int count, MPI_Request array_of_requests[], int *indx,
   MPI_Status *status)
{ trace ("mpi_waitany", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 274);
  { int _ret =  MPI_Waitany (count, array_of_requests, indx, status); end_trace("mpi_waitany", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 275);  return _ret; }
 end_trace("mpi_waitany", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 276); }

static void rcv_pid_receive (RcvPid * m, scalar * list, vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_receive");

  int len = list_len (list) + 2*2*vectors_len (listf);

  MPI_Request r[m->npid];
  Rcv * rrcv[m->npid];
  int nr = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      assert (!rcv->buf);
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);






      MPI_Irecv (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid,
   (l), MPI_COMM_WORLD, &r[nr]);
      rrcv[nr++] = rcv;






    }
  }


  if (nr > 0) {
    int i;
    MPI_Status s;
    mpi_waitany (nr, r, &i, &s);
    while (i != MPI_UNDEFINED) {
      Rcv * rcv = rrcv[i];
      assert (l <= rcv->depth && rcv->halo[l].n > 0);
      assert (rcv->buf);
      apply_bc (rcv, list, listf, l, s);
      mpi_waitany (nr, r, &i, &s);
    }
  }

  prof_stop();
}


static void rcv_pid_wait (RcvPid * m)
{ trace ("rcv_pid_wait", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 332);

  for (int i = 0; i < m->npid; i++)
    rcv_free_buf (&m->rcv[i]);
 end_trace("rcv_pid_wait", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 336); }

static void rcv_pid_send (RcvPid * m, scalar * list, vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_send");

  int len = list_len (list) + 2*2*vectors_len (listf);


  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      assert (!rcv->buf);
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);
      double * b = rcv->buf;
       { foreach_cache_level(rcv->halo[l], l){

#line 354 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
 if (list) for (scalar s = *list, *_i59 = list; ((scalar *)&s)->i >= 0; s = *++_i59)
   *b++ = val(s,0,0,0);
 if (listf) for (vector v = *listf, *_i60 = listf; ((scalar *)&v)->i >= 0; v = *++_i60)
   {
#line 358
 {
     *b++ = val(v.x,0,0,0);
     *b++ = allocated(1,0,) ? val(v.x,1,0,0) : undefined;
   }
#line 358
 {
     *b++ = val(v.y,0,0,0);
     *b++ = allocated(0,1,) ? val(v.y,0,1,0) : undefined;
   }}
      } } end_foreach_cache_level(); }





      MPI_Isend (rcv->buf, (b - (double *) rcv->buf),
   MPI_DOUBLE, rcv->pid, (l), MPI_COMM_WORLD,
   &rcv->r);
    }
  }

  prof_stop();
}

static void rcv_pid_sync (SndRcv * m, scalar * list, int l)
{
  scalar * listr = NULL;
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i61 = list; ((scalar *)&s)->i >= 0; s = *++_i61)
    if (!is_constant(s)) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listr = list_add (listr, s);
    }
  rcv_pid_send (m->snd, listr, listf, l);
  rcv_pid_receive (m->rcv, listr, listf, l);
  rcv_pid_wait (m->snd);
  pfree (listr,__func__,__FILE__,__LINE__);
  pfree (listf,__func__,__FILE__,__LINE__);
}

static void snd_rcv_destroy (SndRcv * m)
{
  rcv_pid_destroy (m->rcv);
  rcv_pid_destroy (m->snd);
}

static void snd_rcv_init (SndRcv * m, const char * name)
{
  char s[strlen(name) + 5];
  strcpy (s, name);
  strcat (s, ".rcv");
  m->rcv = rcv_pid_new (s);
  strcpy (s, name);
  strcat (s, ".snd");
  m->snd = rcv_pid_new (s);
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  snd_rcv_destroy (&m->mpi_level);
  snd_rcv_destroy (&m->mpi_level_root);
  snd_rcv_destroy (&m->restriction);
  array_free (m->send);
  array_free (m->receive);
  pfree (m,__func__,__FILE__,__LINE__);
}


static void mpi_boundary_level (const Boundary * b, scalar * list, int l)
{ trace ("mpi_boundary_level", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 425);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->mpi_level, list, l);
  rcv_pid_sync (&m->mpi_level_root, list, l);
 end_trace("mpi_boundary_level", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 429); }


static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{ trace ("mpi_boundary_restriction", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 433);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l);
 end_trace("mpi_boundary_restriction", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 436); }

void mpi_boundary_new()
{
  mpi_boundary = (Boundary *) ((MpiBoundary *) pcalloc (1, sizeof(MpiBoundary),__func__,__FILE__,__LINE__));
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->level = mpi_boundary_level;
  mpi_boundary->restriction = mpi_boundary_restriction;
  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  snd_rcv_init (&mpi->mpi_level, "mpi_level");
  snd_rcv_init (&mpi->mpi_level_root, "mpi_level_root");
  snd_rcv_init (&mpi->restriction, "restriction");
  mpi->send = array_new();
  mpi->receive = array_new();
  add_boundary (mpi_boundary);
}

static FILE * fopen_prefix (FILE * fp, const char * name, char * prefix)
{
  if (fp) {
    sprintf (prefix, "%s-%d ", name, pid());
    return fp;
  }
  else {
    strcpy (prefix, "");
    char fname[80];
    if (debug_iteration >= 0)
      sprintf (fname, "%s-%d-%d", name, debug_iteration, pid());
    else
      sprintf (fname, "%s-%d", name, pid());
    return fopen (fname, "w");
  }
}

void debug_mpi (FILE * fp1)
{
  void output_cells (FILE * fp);

  char prefix[80];
  FILE * fp;


  if (fp1 == NULL) {
    char name[80];
    sprintf (name, "halo-%d", pid()); remove (name);
    sprintf (name, "cells-%d", pid()); remove (name);
    sprintf (name, "faces-%d", pid()); remove (name);
    sprintf (name, "neighbors-%d", pid()); remove (name);
    sprintf (name, "mpi-level-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-border-%d", pid()); remove (name);
    sprintf (name, "exterior-%d", pid()); remove (name);
    sprintf (name, "depth-%d", pid()); remove (name);
    sprintf (name, "refined-%d", pid()); remove (name);
  }


  fp = fopen_prefix (fp1, "halo", prefix);
  for (int l = 0; l < depth(); l++)
     { foreach_halo (prolongation, l){

#line 499 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

       { foreach_child()
      fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); end_foreach_child(); }; } end_foreach_halo(); }
  if (!fp1)
    fclose (fp);

  if (!fp1) {
    fp = fopen_prefix (fp1, "cells", prefix);
    output_cells (fp);
    fclose (fp);
  }

  fp = fopen_prefix (fp1, "faces", prefix);
   { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 512
{

#line 512 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 512
{

#line 512 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  end_foreach_face_generic()
#line 513
 end_foreach_face(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "neighbors", prefix);
   { foreach(){

#line 518 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
    int n = 0;
     { foreach_neighbor(1)
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 n++; end_foreach_neighbor(); }
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, cell.neighbors);
    assert (cell.neighbors == n);
  } } end_foreach(); }
  if (!fp1)
    fclose (fp);

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;

  fp = fopen_prefix (fp1, "mpi-level-rcv", prefix);
  rcv_pid_print (mpi->mpi_level.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-rcv", prefix);
  rcv_pid_print (mpi->mpi_level_root.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-rcv", prefix);
  rcv_pid_print (mpi->restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-snd", prefix);
  rcv_pid_print (mpi->mpi_level.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-snd", prefix);
  rcv_pid_print (mpi->mpi_level_root.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-snd", prefix);
  rcv_pid_print (mpi->restriction.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-border", prefix);
   { foreach_cell(){

#line 562 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
    if (is_border(cell))
      fprintf (fp, "%s%g %g %g %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors, cell.pid);
    else
      continue;
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "exterior", prefix);
   { foreach_cell(){

#line 575 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
    if (!is_local(cell))
      fprintf (fp, "%s%g %g %g %d %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors,
        cell.pid, cell.flags & leaf);






  } } end_foreach_cell(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "depth", prefix);
  fprintf (fp, "depth: %d %d\n", pid(), depth());
  fprintf (fp, "======= mpi_level.snd ======\n");
  RcvPid * snd = mpi->mpi_level.snd;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  fprintf (fp, "======= mpi_level.rcv ======\n");
  snd = mpi->mpi_level.rcv;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "refined", prefix);
   { foreach_cache (((Tree *)grid)->refined){

#line 604 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); } end_foreach_cache(); }
  if (!fp1)
    fclose (fp);
}

static void snd_rcv_free (SndRcv * p)
{
  char name[strlen(p->rcv->name) + 1];
  strcpy (name, p->rcv->name);
  rcv_pid_destroy (p->rcv);
  p->rcv = rcv_pid_new (name);
  strcpy (name, p->snd->name);
  rcv_pid_destroy (p->snd);
  p->snd = rcv_pid_new (name);
}

static bool is_root (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 622 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

  if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
     { foreach_child()
      if (is_local(cell))
 return true; end_foreach_child(); }
  return false;
}


static bool is_local_prolongation (Point point, Point p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 632 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"


  struct { int x, y; } dp = {p.i - point.i, p.j - point.j};



  {
#line 638
 {
    if (dp.x == 0 && ((!is_leaf (neighbor(-1,0,)) && neighbor(-1,0,).neighbors && neighbor(-1,0,).pid >= 0) || (!is_leaf (neighbor(1,0,)) && neighbor(1,0,).neighbors && neighbor(1,0,).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(dp.x,0,)) && neighbor(dp.x,0,).neighbors && neighbor(dp.x,0,).pid >= 0))
      return true;
  }
#line 638
 {
    if (dp.y == 0 && ((!is_leaf (neighbor(0,-1,)) && neighbor(0,-1,).neighbors && neighbor(0,-1,).pid >= 0) || (!is_leaf (neighbor(0,1,)) && neighbor(0,1,).neighbors && neighbor(0,1,).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(0,dp.y,)) && neighbor(0,dp.y,).neighbors && neighbor(0,dp.y,).pid >= 0))
      return true;
  }}
  return false;
}



static void append_pid (Array * pids, int pid)
{
  for (int i = 0, * p = (int *) pids->p; i < pids->len/sizeof(int); i++, p++)
    if (*p == pid)
      return;
  array_append (pids, &pid, sizeof(int));
}

static int locals_pids (Point point, Array * pids)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 658 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

  if (is_leaf(cell)) {
    if (is_local(cell)) {
      Point p = point;
       { foreach_neighbor(1) {
 if ((cell.pid >= 0 && cell.pid != pid()) &&
     ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p)))
   append_pid (pids, cell.pid);
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    { foreach_child()
     if ((cell.pid >= 0 && cell.pid != pid()))
       append_pid (pids, cell.pid); end_foreach_child(); }
      } end_foreach_neighbor(); }
    }
  }
  else
     { foreach_neighbor(1) {
      if ((cell.pid >= 0 && cell.pid != pid()))
 append_pid (pids, cell.pid);
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
  { foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
     append_pid (pids, cell.pid); end_foreach_child(); }
    } end_foreach_neighbor(); }
  return pids->len/sizeof(int);
}

static int root_pids (Point point, Array * pids)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 686 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

   { foreach_child()
    if ((cell.pid >= 0 && cell.pid != pid()))
      append_pid (pids, cell.pid); end_foreach_child(); }
  return pids->len/sizeof(int);
}







static void rcv_pid_row (RcvPid * m, int l, int * row)
{
  for (int i = 0; i < npe(); i++)
    row[i] = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0)
      row[rcv->pid] = rcv->halo[l].n;
  }
}

void check_snd_rcv_matrix (SndRcv * sndrcv, const char * name)
{
  int maxlevel = depth();
  mpi_all_reduce (maxlevel, MPI_INT, MPI_MAX);
  int * row = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,__LINE__));
  for (int l = 0; l <= maxlevel; l++) {
    int status = 0;
    if (pid() == 0) {


      int ** send = matrix_new (npe(), npe(), sizeof(int));
      int ** receive = matrix_new (npe(), npe(), sizeof(int));
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, &send[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, &receive[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);

      int * astatus = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,__LINE__));
      for (int i = 0; i < npe(); i++)
 astatus[i] = 0;
      for (int i = 0; i < npe(); i++)
 for (int j = 0; j < npe(); j++)
   if (send[i][j] != receive[j][i]) {
     fprintf (qstderr(), "%s: %d sends    %d to   %d at level %d\n",
       name, i, send[i][j], j, l);
     fprintf (qstderr(), "%s: %d receives %d from %d at level %d\n",
       name, j, receive[j][i], i, l);
     fflush (qstderr());
     for (int k = i - 2; k <= i + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
     for (int k = j - 2; k <= j + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
   }
      MPI_Scatter (astatus, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      pfree (astatus,__func__,__FILE__,__LINE__);

      matrix_free (send);
      matrix_free (receive);
    }
    else {
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Scatter (NULL, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (status) {
      fprintf (qstderr(),
        "check_snd_rcv_matrix \"%s\" failed\n"
        "Calling debug_mpi(NULL)...\n"
        "Aborting...\n",
        name);
      fflush (qstderr());
      debug_mpi (NULL);
      MPI_Abort (MPI_COMM_WORLD, -3);
    }
  }
  pfree (row,__func__,__FILE__,__LINE__);
}

static bool has_local_child (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 775 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

   { foreach_child()
    if (is_local(cell))
      return true; end_foreach_child(); }
  return false;
}


void mpi_boundary_update_buffers()
{ trace ("mpi_boundary_update_buffers", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 784);
  if (npe() == 1)
    { ; end_trace("mpi_boundary_update_buffers", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 786);  return; }

  prof_start ("mpi_boundary_update_buffers");

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  SndRcv * mpi_level = &m->mpi_level;
  SndRcv * mpi_level_root = &m->mpi_level_root;
  SndRcv * restriction = &m->restriction;

  snd_rcv_free (mpi_level);
  snd_rcv_free (mpi_level_root);
  snd_rcv_free (restriction);

  static const unsigned short used = 1 << user;
   { foreach_cell(){

#line 800 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
    if (is_active(cell) && !is_border(cell))



      continue;

    if (cell.neighbors) {

      Array pids = {NULL, 0, 0};
      int n = locals_pids (point, &pids);
      if (n) {
  { foreach_child()
   if (is_local(cell))
     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (mpi_level->snd, *p, point); end_foreach_child(); }
 pfree (pids.p,__func__,__FILE__,__LINE__);
      }

      bool locals = false;
      if (is_leaf(cell)) {
 if ((cell.pid >= 0 && cell.pid != pid())) {
   Point p = point;
    { foreach_neighbor(1)
     if ((is_local(cell) &&
   ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p))) ||
  is_root(point))
       locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
 }
      }
      else
  { foreach_neighbor(1)
   if (is_local(cell) || is_root(point))
     locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
      if (locals)
  { foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
            rcv_pid_append (mpi_level->rcv, cell.pid, point),
       cell.flags |= used; end_foreach_child(); }


      if (!is_leaf(cell)) {

 if (is_local(cell)) {
   Array pids = {NULL, 0, 0};

   int n = root_pids (point, &pids);
   if (n) {
      { foreach_neighbor()
       for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
  if (cell.pid >= 0 && cell.pid != *p)
    rcv_pid_append (mpi_level_root->snd, *p, point); end_foreach_neighbor(); }

     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (restriction->snd, *p, point);
     pfree (pids.p,__func__,__FILE__,__LINE__);
   }
 }

 else if ((cell.pid >= 0 && cell.pid != pid())) {
   bool root = false;
    { foreach_child()
     if (is_local(cell))
       root = true, foreach_child_break(); end_foreach_child(); }
   if (root) {
     int pid = cell.pid;
      { foreach_neighbor()
       if ((cell.pid >= 0 && cell.pid != pid()))
  rcv_pid_append (mpi_level_root->rcv, pid, point),
    cell.flags |= used; end_foreach_neighbor(); }

     rcv_pid_append (restriction->rcv, pid, point);
   }
 }
      }
    }


    if (level > 0) {
      if (is_local(cell)) {

 Array pids = {NULL, 0, 0};
 if ((aparent(0,0,0).pid >= 0 && aparent(0,0,0).pid != pid()))
   append_pid (&pids, aparent(0,0,0).pid);
 int n = root_pids (parent, &pids);
 if (n) {
   for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
     rcv_pid_append (restriction->snd, *p, point);
   pfree (pids.p,__func__,__FILE__,__LINE__);
 }
      }
      else if ((cell.pid >= 0 && cell.pid != pid())) {

 if (is_local(aparent(0,0,0)) || has_local_child (parent))
   rcv_pid_append (restriction->rcv, cell.pid, point);
      }
    }
  } } end_foreach_cell(); }





  static const unsigned short keep = 1 << (user + 1);
  for (int l = depth(); l >= 0; l--)
     { foreach_cell(){

#line 905 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

      if (level == l) {
 if (level > 0 && (cell.pid < 0 || is_local(cell) || (cell.flags & used)))
   aparent(0,0,0).flags |= keep;
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !(cell.flags & keep))
   coarsen_cell (point, NULL);
 cell.flags &= ~(used|keep);
 continue;
      } } end_foreach_cell(); }


  m->send->len = m->receive->len = 0;
  rcv_pid_append_pids (mpi_level->snd, m->send);
  rcv_pid_append_pids (mpi_level_root->snd, m->send);
  rcv_pid_append_pids (mpi_level->rcv, m->receive);
  rcv_pid_append_pids (mpi_level_root->rcv, m->receive);

  prof_stop();
#line 937 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 end_trace("mpi_boundary_update_buffers", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 937); }


void mpi_boundary_refine (scalar * list)
{ trace ("mpi_boundary_refine", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 941);
  prof_start ("mpi_boundary_refine");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;


  Array * snd = mpi->send;
  MPI_Request r[2*snd->len/sizeof(int)];
  int nr = 0;
  for (int i = 0, * dest = snd->p; i < snd->len/sizeof(int); i++,dest++) {
    int len = ((Tree *)grid)->refined.n;
    MPI_Isend (&((Tree *)grid)->refined.n, 1, MPI_INT, *dest,
        (128), MPI_COMM_WORLD, &r[nr++]);
    if (len > 0)
      MPI_Isend (((Tree *)grid)->refined.p, sizeof(Index)/sizeof(int)*len,
   MPI_INT, *dest, (128), MPI_COMM_WORLD, &r[nr++]);
  }



  Array * rcv = mpi->receive;
  Cache rerefined = {NULL, 0, 0};
  for (int i = 0, * source = rcv->p; i < rcv->len/sizeof(int); i++,source++) {
    int len;
    mpi_recv_check (&len, 1, MPI_INT, *source, (128),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE,
      "mpi_boundary_refine (len)");
    if (len > 0) {
      Index p[len];
      mpi_recv_check (p, sizeof(Index)/sizeof(int)*len,
        MPI_INT, *source, (128),
        MPI_COMM_WORLD, MPI_STATUS_IGNORE,
        "mpi_boundary_refine (p)");
      Cache refined = {p, len, len};
       { foreach_cache (refined){

#line 975 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

 if (level <= depth() && allocated(0,0,0)) {
   if (is_leaf(cell)) {
     bool neighbors = false;
      { foreach_neighbor()
       if (allocated(0,0,0) && (is_active(cell) || is_local(aparent(0,0,0))))
  neighbors = true, foreach_neighbor_break(); end_foreach_neighbor(); }

     if (neighbors)
       refine_cell (point, list, 0, &rerefined);
   }
 } } end_foreach_cache(); }
    }
  }


  if (nr)
    MPI_Waitall (nr, r, MPI_STATUSES_IGNORE);


  pfree (((Tree *)grid)->refined.p,__func__,__FILE__,__LINE__);
  ((Tree *)grid)->refined = rerefined;

  prof_stop();



  mpi_all_reduce (rerefined.n, MPI_INT, MPI_SUM);
  if (rerefined.n)
    mpi_boundary_refine (list);
 end_trace("mpi_boundary_refine", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 1005); }

static void check_depth()
{
#line 1040 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
}

typedef struct {
  int refined, leaf;
} Remote;




void mpi_boundary_coarsen (int l, int too_fine)
{ trace ("mpi_boundary_coarsen", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 1050);
  if (npe() == 1)
    { ; end_trace("mpi_boundary_coarsen", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 1052);  return; }

  check_depth();

  assert (sizeof(Remote) == sizeof(double));

  scalar remote= new_scalar("remote");
   { foreach_cell(){

#line 1059 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
    if (level == l) {
      if (is_local(cell)) {
 ((Remote *)&val(remote,0,0,0))->refined = (!is_leaf (cell) && cell.neighbors && cell.pid >= 0);
 ((Remote *)&val(remote,0,0,0))->leaf = is_leaf(cell);
      }
      else {
 ((Remote *)&val(remote,0,0,0))->refined = true;
 ((Remote *)&val(remote,0,0,0))->leaf = false;
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  mpi_boundary_level (mpi_boundary, ((scalar []){remote,{-1}}), l);

   { foreach_cell(){

#line 1076 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
    if (level == l) {
      if (!is_local(cell)) {
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !((Remote *)&val(remote,0,0,0))->refined)
   coarsen_cell_recursive (point, NULL);
 else if (is_leaf(cell) && cell.neighbors && ((Remote *)&val(remote,0,0,0))->leaf) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  check_depth();

  if (l > 0) {
     { foreach_cell(){

#line 1096 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
      if (level == l) {
 val(remote,0,0,0) = is_local(cell) ? cell.neighbors : 0;
 continue;
      }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
    mpi_boundary_level (mpi_boundary, ((scalar []){remote,{-1}}), l);
     { foreach_cell(){

#line 1105 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
      if (level == l)
 if (!is_local(cell) && is_local(aparent(0,0,0)) && val(remote,0,0,0)) {
   aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
  }
 delete (((scalar []){remote,{-1}}));  end_trace("mpi_boundary_coarsen", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 1115); }

static void flag_border_cells()
{
   { foreach_cell(){

#line 1119 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
    if (is_active(cell)) {
      short flags = cell.flags & ~border;
       { foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0))))
   flags |= border, foreach_neighbor_break();

 if (is_refined_check())
    { foreach_child()
     if (!is_local(cell))
       flags |= border, foreach_child_break(); end_foreach_child(); }
 if (flags & border)
   foreach_neighbor_break();
      } end_foreach_neighbor(); }
      cell.flags = flags;
    }
    else {
      cell.flags &= ~border;

    }
    if (is_leaf(cell)) {
      if (cell.neighbors) {
  { foreach_child()
   cell.flags &= ~border; end_foreach_child(); }
 if (is_border(cell)) {
   bool remote = false;
    { foreach_neighbor (2/2)
     if (!is_local(cell))
       remote = true, foreach_neighbor_break(); end_foreach_neighbor(); }
   if (remote)
      { foreach_child()
       cell.flags |= border; end_foreach_child(); }
 }
      }
      continue;
    }
  } } end_foreach_cell(); }
}

static int balanced_pid (long index, long nt, int nproc)
{
  long ne = max(1, nt/nproc), nr = nt % nproc;
  int pid = index < nr*(ne + 1) ?
    index/(ne + 1) :
    nr + (index - nr*(ne + 1))/ne;
  return min(nproc - 1, pid);
}



void mpi_partitioning()
{ trace ("mpi_partitioning", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 1170);
  prof_start ("mpi_partitioning");

  long nt = 0;
   { foreach(){

#line 1174 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

    nt++; } end_foreach(); }


  long i = 0;
  ((Tree *)grid)->dirty = true;
   { foreach_cell_post (is_active (cell)){

#line 1180 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

    if (is_active (cell)) {
      if (is_leaf (cell)) {
 cell.pid = balanced_pid (i++, nt, npe());
 if (cell.neighbors > 0) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else {
 cell.pid = child(0,0,0).pid;
 bool inactive = true;
  { foreach_child()
   if (is_active(cell))
     inactive = false, foreach_child_break(); end_foreach_child(); }
 if (inactive)
   cell.flags &= ~active;
      }
    } } end_foreach_cell_post(); }

  flag_border_cells();

  prof_stop();

  mpi_boundary_update_buffers();
 end_trace("mpi_partitioning", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 1208); }

void restore_mpi (FILE * fp, scalar * list1)
{
  long index = 0, nt = 0, start = ftell (fp);
  scalar size= new_scalar("size"), * list = list_concat (((scalar []){size,{-1}}), list1);;
  long offset = sizeof(double)*list_len(list);


  static const unsigned short set = 1 << user;
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector []){{fm.x,fm.y},{{-1},{-1}}});
   { foreach_cell(){

#line 1219 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

    if (balanced_pid (index, nt, npe()) <= pid()) {
      unsigned flags;
      if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
 fprintf (qstderr(), "restore(): error: expecting 'flags'\n");
 exit (1);
      }
      if (list) for (scalar s = *list, *_i62 = list; ((scalar *)&s)->i >= 0; s = *++_i62) {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (qstderr(), "restore(): error: expecting scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }
      if (level == 0)
 nt = val(size,0,0,0);
      cell.pid = balanced_pid (index, nt, npe());
      cell.flags |= set;
      if (!(flags & leaf) && is_leaf(cell)) {
 if (balanced_pid (index + val(size,0,0,0) - 1, nt, npe()) < pid()) {
   fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
   index += val(size,0,0,0);
   continue;
 }
 refine_cell (point, listm, 0, NULL);
      }
      index++;
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }


  fseek (fp, start, SEEK_SET);
  index = 0;
   { foreach_cell(){

#line 1255 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (qstderr(), "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    if (cell.flags & set)
      fseek (fp, offset, SEEK_CUR);
    else {
      if (list) for (scalar s = *list, *_i63 = list; ((scalar *)&s)->i >= 0; s = *++_i63) {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (qstderr(), "restore(): error: expecting a scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }
      cell.pid = balanced_pid (index, nt, npe());
      if (is_leaf(cell) && cell.neighbors) {
 int pid = cell.pid;
  { foreach_child()
   cell.pid = pid; end_foreach_child(); }
      }
    }
    if (!(flags & leaf) && is_leaf(cell)) {
      bool locals = false;
       { foreach_neighbor(1)
 if ((cell.flags & set) && (is_local(cell) || is_root(point)))
   locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
      if (locals)
 refine_cell (point, listm, 0, NULL);
      else {
 fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
 index += val(size,0,0,0);
 continue;
      }
    }
    index++;
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }


   { foreach_cell_post (is_active (cell)){

#line 1299 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
    cell.flags &= ~set;
    if (is_active (cell)) {
      if (is_leaf (cell)) {
 if (cell.neighbors > 0) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else if (!is_local(cell)) {
 bool inactive = true;
  { foreach_child()
   if (is_active(cell))
     inactive = false, foreach_child_break(); end_foreach_child(); }
 if (inactive)
   cell.flags &= ~active;
      }
    }
  } } end_foreach_cell_post(); }

  flag_border_cells();

  mpi_boundary_update (list);
  pfree (list,__func__,__FILE__,__LINE__);
 delete (((scalar []){size,{-1}})); }
#line 1348 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

double z_indexing (scalar index, bool leaves)
{ trace ("z_indexing", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 1350);



  scalar size= new_scalar("size");
  subtree_size (size, leaves);






  double maxi = -1.;
  if (pid() == 0)
     { foreach_level(0){

#line 1364 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

      maxi = val(size,0,0,0) - 1.; } end_foreach_level(); }




   { foreach_level(0){

#line 1370 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"

    val(index,0,0,0) = 0; } end_foreach_level(); }
  for (int l = 0; l < depth(); l++) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){index,{-1}}), l); };
     { foreach_cell(){

#line 1374 "/home/jiarongw/basilisk/src/grid/tree-mpi.h"
 {
      if (level == l) {
 if (is_leaf(cell)) {
   if (is_local(cell) && cell.neighbors) {
     int i = val(index,0,0,0);
      { foreach_child()
       val(index,0,0,0) = i; end_foreach_child(); }
   }
 }
 else {
   bool loc = is_local(cell);
   if (!loc)
      { foreach_child()
       if (is_local(cell))
  loc = true, foreach_child_break(); end_foreach_child(); }
   if (loc) {
     int i = val(index,0,0,0) + !leaves;
      { foreach_child() {
       val(index,0,0,0) = i;
       i += val(size,0,0,0);
     } end_foreach_child(); }
   }
 }
 continue;
      }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
  }
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){index,{-1}}), depth()); };

  { double _ret =  maxi; delete (((scalar []){size,{-1}}));  end_trace("z_indexing", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 1405);  return _ret; }
 delete (((scalar []){size,{-1}}));  end_trace("z_indexing", "/home/jiarongw/basilisk/src/grid/tree-mpi.h", 1406); }
#line 1783 "/home/jiarongw/basilisk/src/grid/tree.h"
#line 1 "grid/balance.h"
#line 1 "/home/jiarongw/basilisk/src/grid/balance.h"


typedef struct {
  short leaf, prolongation;
  int pid;
} NewPid;



#if TRASH
# define is_newpid() (!isnan(val(newpid,0,0,0)) && ((NewPid *)&val(newpid,0,0,0))->pid > 0)
#else
# define is_newpid() (((NewPid *)&val(newpid,0,0,0))->pid > 0)
#endif

Array * linear_tree (size_t size, scalar newpid)
{
  const unsigned short sent = 1 << user, next = 1 << (user + 1);
  Array * a = array_new();

   { foreach_cell_post_all (true){

#line 21 "/home/jiarongw/basilisk/src/grid/balance.h"

    if (level > 0 && (cell.flags & (sent|next)))
      aparent(0,0,0).flags |= next; } end_foreach_cell_post_all(); }

  bool empty = true;
   { foreach_cell_all(){

#line 26 "/home/jiarongw/basilisk/src/grid/balance.h"
 {
    if (cell.flags & sent) {
      array_append (a, &cell, size);
      cell.flags &= ~sent;
      empty = false;
    }
    else {
      if (cell.pid >= 0 && ((NewPid *)&val(newpid,0,0,0))->leaf)
 assert (is_leaf(cell));
      if (is_refined_check()) {


 bool prolo = false;
  { foreach_child()
   if (((NewPid *)&val(newpid,0,0,0))->prolongation)
     prolo = true; end_foreach_child(); }
 if (prolo) {

   cell.flags |= leaf;
   array_append (a, &cell, sizeof(Cell));
   cell.flags &= ~leaf;
 }
 else
   array_append (a, &cell, sizeof(Cell));
      }
      else
 array_append (a, &cell, sizeof(Cell));
    }
    if (cell.flags & next)
      cell.flags &= ~next;
    else
      continue;
  } } end_foreach_cell_all(); }

  if (empty)
    a->len = 0;
  return a;
}

#define foreach_tree(t, size, list)\
{\
  const unsigned short _sent = 1 << user, _next = 1 << (user + 1);\
  scalar * _list = list;\
  char * _i = (char *) (t)->p;\
  foreach_cell_all() {\
    Cell * c = (Cell *) _i;\
    if (c->flags & _sent) {\
      _i += size;\

#line 74


#define end_foreach_tree()\
    }\
    else\
      _i += sizeof(Cell);\
    if (c->flags & _next) {\
      assert (c->neighbors);\
      if (!(c->flags & leaf) && is_leaf(cell) &&\
   (!is_newpid() || !((NewPid *)&val(newpid,0,0,0))->leaf))\
\
 refine_cell (point, _list, 0, NULL);\
      else if (!cell.neighbors)\
\
 alloc_children (point);\
    }\
    else\
      continue;\
  } end_foreach_cell_all();\
}\

#line 94


Array * neighborhood (scalar newpid, int nextpid, FILE * fp)
{
  const unsigned short sent = 1 << user;
   { foreach_cell(){

#line 99 "/home/jiarongw/basilisk/src/grid/balance.h"
 {

    bool root = false;
    if ((!is_local(cell) || ((NewPid *)&val(newpid,0,0,0))->pid - 1 != nextpid) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0)) {
       { foreach_child()
 if (is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid)
   root = true, foreach_child_break(); end_foreach_child(); }
      if (root && cell.pid != nextpid) {
  { foreach_neighbor()
   if (cell.pid != nextpid && is_newpid()) {
     if (fp)
       fprintf (fp, "%g %g %g %d %d root\n",
         x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
     cell.flags |= sent;
   } end_foreach_neighbor(); }
      }
    }

    if ((is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid) || root) {
       { foreach_neighbor(1)
 if (cell.neighbors && cell.pid != nextpid)
    { foreach_child()
     if (cell.pid != nextpid && is_newpid()) {
       if (fp)
  fprintf (fp, "%g %g %g %d %d nextpid\n",
    x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
       cell.flags |= sent;
     } end_foreach_child(); } end_foreach_neighbor(); }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  return linear_tree (sizeof(Cell) + datasize, newpid);
}

static void send_tree (Array * a, int to, MPI_Request * r)
{
  MPI_Isend (&a->len, 1, MPI_LONG, to, (256), MPI_COMM_WORLD, &r[0]);
  if (a->len > 0) {
    MPI_Isend (a->p, a->len, MPI_BYTE, to, (256), MPI_COMM_WORLD, &r[1]);
    ((Tree *)grid)->dirty = true;
  }
}

static void receive_tree (int from, scalar newpid, FILE * fp)
{
  Array a;
  mpi_recv_check (&a.len, 1, MPI_LONG, from, (256),
    MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (len)");
  if (a.len > 0) {
    a.p = pmalloc (a.len,__func__,__FILE__,__LINE__);
    if (fp)
      fprintf (fp, "receiving %ld from %d\n", a.len, from);
    mpi_recv_check (a.p, a.len, MPI_BYTE, from, (256),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (p)");

     { foreach_tree (&a, sizeof(Cell) + datasize, NULL){

#line 156 "/home/jiarongw/basilisk/src/grid/balance.h"
 {
      memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
       datasize);
      assert (((NewPid *)&val(newpid,0,0,0))->pid > 0);
      if (fp)
 fprintf (fp, "%g %g %g %d %d %d %d %d %d recv\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
   c->flags & leaf,
   cell.flags & leaf, from, ((NewPid *)&val(newpid,0,0,0))->leaf);
    } } end_foreach_tree(); }
    pfree (a.p,__func__,__FILE__,__LINE__);
    ((Tree *)grid)->dirty = true;
  }
}

static void wait_tree (Array * a, MPI_Request * r)
{
  MPI_Wait (&r[0], MPI_STATUS_IGNORE);
  if (a->len > 0)
    MPI_Wait (&r[1], MPI_STATUS_IGNORE);
}

static void check_flags()
{







}

struct {
  int min;
  bool leaves;

  int npe;
} mpi = {
  1,
  true
};


bool balance()
{ trace ("balance", "/home/jiarongw/basilisk/src/grid/balance.h", 201);
  if (npe() == 1)
    { bool _ret =  false; end_trace("balance", "/home/jiarongw/basilisk/src/grid/balance.h", 203);  return _ret; }

  assert (sizeof(NewPid) == sizeof(double));

  check_flags();

  long nl = 0, nt = 0;
   { foreach_cell(){

#line 210 "/home/jiarongw/basilisk/src/grid/balance.h"
 {
    if (is_local(cell)) {
      nt++;
      if (is_leaf(cell))
 nl++;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  grid->n = grid->tn = nl;
  grid->maxdepth = depth();
  long nmin = nl, nmax = nl;

  mpi_all_reduce (nmax, MPI_LONG, MPI_MAX);
  mpi_all_reduce (nmin, MPI_LONG, MPI_MIN);
  mpi_all_reduce (grid->tn, MPI_LONG, MPI_SUM);
  mpi_all_reduce (grid->maxdepth, MPI_INT, MPI_MAX);
  if (mpi.leaves)
    nt = grid->tn;
  else
    mpi_all_reduce (nt, MPI_LONG, MPI_SUM);

  long ne = max(1, nt/npe());

  if (ne < mpi.min) {
    mpi.npe = max(1, nt/mpi.min);
    ne = max(1, nt/mpi.npe);
  }
  else
    mpi.npe = npe();

  if (nmax - nmin <= 1)
    { bool _ret =  false; end_trace("balance", "/home/jiarongw/basilisk/src/grid/balance.h", 243);  return _ret; }

  scalar newpid= new_scalar("newpid");
  double zn = z_indexing (newpid, mpi.leaves);
  if (pid() == 0)
    assert (zn + 1 == nt);

  FILE * fp = NULL;
#line 260 "/home/jiarongw/basilisk/src/grid/balance.h"
  bool next = false, prev = false;
   { foreach_cell_all(){

#line 261 "/home/jiarongw/basilisk/src/grid/balance.h"
 {
    if (is_local(cell)) {
      int pid = balanced_pid (val(newpid,0,0,0), nt, mpi.npe);
      pid = clamp (pid, cell.pid - 1, cell.pid + 1);
      if (pid == pid() + 1)
 next = true;
      else if (pid == pid() - 1)
 prev = true;
      ((NewPid *)&val(newpid,0,0,0))->pid = pid + 1;
      ((NewPid *)&val(newpid,0,0,0))->leaf = is_leaf(cell);
      ((NewPid *)&val(newpid,0,0,0))->prolongation = (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0);
      if (fp)
 fprintf (fp, "%g %g %d %d newpid\n", x, y, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
    }
    else
      val(newpid,0,0,0) = 0;
  } } end_foreach_cell_all(); }
  for (int l = 0; l <= depth(); l++)
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, ((scalar []){newpid,{-1}}), l); };
#line 304 "/home/jiarongw/basilisk/src/grid/balance.h"
  Array * anext = next ? neighborhood (newpid, pid() + 1, fp) : array_new();
  Array * aprev = prev ? neighborhood (newpid, pid() - 1, fp) : array_new();

  if (fp)
    fflush (fp);

  check_flags();


  MPI_Request rprev[2], rnext[2];
  if (pid() > 0)
    send_tree (aprev, pid() - 1, rprev);
  if (pid() < npe() - 1)
    send_tree (anext, pid() + 1, rnext);


  if (pid() < npe() - 1)
    receive_tree (pid() + 1, newpid, fp);
  if (pid() > 0)
    receive_tree (pid() - 1, newpid, fp);


  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);

  if (fp)
    fflush (fp);


  int pid_changed = false;
   { foreach_cell_all(){

#line 338 "/home/jiarongw/basilisk/src/grid/balance.h"
 {
    if (cell.pid >= 0) {
      if (is_newpid()) {
 if (fp)
   fprintf (fp, "%g %g %g %d %d %d %d %d new\n",
     x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
     is_leaf(cell), cell.neighbors, ((NewPid *)&val(newpid,0,0,0))->leaf);
 if (cell.pid != ((NewPid *)&val(newpid,0,0,0))->pid - 1) {
   cell.pid = ((NewPid *)&val(newpid,0,0,0))->pid - 1;
   cell.flags &= ~(active|border);
   if (is_local(cell))
     cell.flags |= active;
   pid_changed = true;
 }
 if (((NewPid *)&val(newpid,0,0,0))->leaf && !is_leaf(cell) && cell.neighbors)
   coarsen_cell_recursive (point, NULL);
      }
      else if (level > 0 && ((NewPid *)&coarse(newpid,0,0,0))->leaf)
 cell.pid = aparent(0,0,0).pid;
    }

    if (!cell.neighbors && allocated_child(0,0,0)) {
      if (fp)
 fprintf (fp, "%g %g %g %d %d freechildren\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
      free_children (point);
    }
  } } end_foreach_cell_all(); }

  if (((Tree *)grid)->dirty || pid_changed) {


     { foreach_cell_post (!is_leaf (cell)){

#line 370 "/home/jiarongw/basilisk/src/grid/balance.h"

      if (!is_leaf(cell) && !is_local(cell)) {
 unsigned short flags = cell.flags & ~active;
  { foreach_child()
   if (is_active(cell))
     flags |= active, foreach_child_break(); end_foreach_child(); }
 cell.flags = flags;
      } } end_foreach_cell_post(); }

    flag_border_cells();
    pid_changed = true;
  }

  if (fp)
    fclose (fp);

  mpi_all_reduce (pid_changed, MPI_INT, MPI_MAX);
  if (pid_changed)
    mpi_boundary_update_buffers();

  { bool _ret =  pid_changed; delete (((scalar []){newpid,{-1}}));  end_trace("balance", "/home/jiarongw/basilisk/src/grid/balance.h", 390);  return _ret; }
 delete (((scalar []){newpid,{-1}}));  end_trace("balance", "/home/jiarongw/basilisk/src/grid/balance.h", 391); }

void mpi_boundary_update (scalar * list)
{
  mpi_boundary_update_buffers();
  grid->tn = 0;
  boundary (list);
  while (balance());
}
#line 1784 "/home/jiarongw/basilisk/src/grid/tree.h"
#else
void mpi_boundary_refine (scalar * list){}
void mpi_boundary_coarsen (int a, int b){}
void mpi_boundary_update (scalar * list) {
  boundary (list);
}
#endif
#line 4 "/home/jiarongw/basilisk/src/grid/quadtree.h"

void quadtree_methods() {
  tree_methods();
}
#line 15 "wavewind-cpp.c"
static double _boundary0 (Point point, Point neighbor, scalar _s);
static double _boundary0_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary1 (Point point, Point neighbor, scalar _s);
static double _boundary1_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary2 (Point point, Point neighbor, scalar _s);
static double _boundary2_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary3 (Point point, Point neighbor, scalar _s);
static double _boundary3_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary4 (Point point, Point neighbor, scalar _s);
static double _boundary4_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary5 (Point point, Point neighbor, scalar _s);
static double _boundary5_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary6 (Point point, Point neighbor, scalar _s);
static double _boundary6_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary7 (Point point, Point neighbor, scalar _s);
static double _boundary7_homogeneous (Point point, Point neighbor, scalar _s);
#line 1 "wavewind.c"
#line 9 "wavewind.c"
#line 1 "navier-stokes/centered.h"
#line 1 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
#line 26 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
#line 1 "./run.h"
#line 1 "/home/jiarongw/basilisk/src/run.h"
#line 9 "/home/jiarongw/basilisk/src/run.h"
double dt = 1.;

#line 1 "./utils.h"
#line 1 "/home/jiarongw/basilisk/src/utils.h"







double DT = 1e10, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf;





void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}






typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;






timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
#if 1
  s.avg = mpi_time - t.tm;
#endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _n = n; 
#line 69
foreach(){

#line 69 "/home/jiarongw/basilisk/src/utils.h"
 _n++; } end_foreach();OMP(omp critical) n += _n;
mpi_all_reduce_double (n, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 69
 }
    s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
#if _GNU_SOURCE
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
#else
  s.mem = 0;
#endif
#if 1
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
#else
  s.min = s.max = s.avg = 0.;
#endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}




void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Quadtree"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if 1
  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);
#endif
}







typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _max = max; double _avg = avg; double _rms = rms; double _volume = volume; 
#line 135

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 135
foreach(){

#line 136 "/home/jiarongw/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      double v = fabs(val(f,0,0,0));
      if (v > _max) _max = v;
      _volume += (sq(Delta)*val_cm(cm,0,0,0));
      _avg += (sq(Delta)*val_cm(cm,0,0,0))*v;
      _rms += (sq(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 135
foreach(){

#line 136 "/home/jiarongw/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      double v = fabs(val(f,0,0,0));
      if (v > _max) _max = v;
      _volume += (sq(Delta)*val_cm(cm,0,0,0));
      _avg += (sq(Delta)*val_cm(cm,0,0,0))*v;
      _rms += (sq(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }OMP(omp critical) if (_max > max) max = _max;
mpi_all_reduce_double (max, MPI_MAX);
OMP(omp critical) avg += _avg;
mpi_all_reduce_double (avg, MPI_SUM);
OMP(omp critical) rms += _rms;
mpi_all_reduce_double (rms, MPI_SUM);
OMP(omp critical) volume += _volume;
mpi_all_reduce_double (volume, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
 }
  norm n;
  n.avg = avg/volume;
  n.rms = sqrt(rms/volume);
  n.max = max;
  n.volume = volume;
  return n;
}





typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _sum = sum; double _sum2 = sum2; double _volume = volume; double _max = max; double _min = min; 
#line 163

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 163
foreach(){

#line 164 "/home/jiarongw/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      _volume += (sq(Delta)*val_cm(cm,0,0,0));
      _sum += (sq(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      _sum2 += (sq(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > _max) _max = val(f,0,0,0);
      if (val(f,0,0,0) < _min) _min = val(f,0,0,0);
    } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 163
foreach(){

#line 164 "/home/jiarongw/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      _volume += (sq(Delta)*val_cm(cm,0,0,0));
      _sum += (sq(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      _sum2 += (sq(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > _max) _max = val(f,0,0,0);
      if (val(f,0,0,0) < _min) _min = val(f,0,0,0);
    } } end_foreach(); }OMP(omp critical) sum += _sum;
mpi_all_reduce_double (sum, MPI_SUM);
OMP(omp critical) sum2 += _sum2;
mpi_all_reduce_double (sum2, MPI_SUM);
OMP(omp critical) volume += _volume;
mpi_all_reduce_double (volume, MPI_SUM);
OMP(omp critical) if (_max > max) max = _max;
mpi_all_reduce_double (max, MPI_MAX);
OMP(omp critical) if (_min < min) min = _min;
mpi_all_reduce_double (min, MPI_MIN);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
 }
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
#line 186 "/home/jiarongw/basilisk/src/utils.h"
static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}
#line 212 "/home/jiarongw/basilisk/src/utils.h"
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}
#line 236 "/home/jiarongw/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
   { foreach(){

#line 239 "/home/jiarongw/basilisk/src/utils.h"
 {
    scalar s; vector v;
    scalar * _i0 = f; vector * _i1 = g; if (f) for (s = *f, v = *g; ((scalar *)&s)->i >= 0; s = *++_i0, v = *++_i1) {
      if (_attribute[s.i].gradient)
 {
#line 243

   val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
#line 243

   val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;}
      else
 {
#line 246

   val(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
#line 246

   val(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);}
    }
  } } end_foreach(); }
  boundary ((scalar *) g);
}
#line 262 "/home/jiarongw/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{

     { foreach(){

#line 265 "/home/jiarongw/basilisk/src/utils.h"

      val(omega,0,0,0) = (val(u.y,1,0,0) - val(u.y,-1,0,0) + val(u.x,0,-1,0) - val(u.x,0,1,0))/(2.*Delta); } end_foreach(); }
    boundary (((scalar []){omega,{-1}}));

}





double change (scalar v, scalar vn)
{
  double max = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _max = max; 
#line 278
foreach(){

#line 278 "/home/jiarongw/basilisk/src/utils.h"
 {
    double dv = fabs (val(v,0,0,0) - val(vn,0,0,0));
    if (dv > _max)
      _max = dv;
    val(vn,0,0,0) = val(v,0,0,0);
  } } end_foreach();OMP(omp critical) if (_max > max) max = _max;
mpi_all_reduce_double (max, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 283
 }
  return max;
}

#line 1 "./output.h"
#line 1 "/home/jiarongw/basilisk/src/output.h"
#line 33 "/home/jiarongw/basilisk/src/output.h"
struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
};


void output_field (struct OutputField p)
{ trace ("output_field", "/home/jiarongw/basilisk/src/output.h", 42);
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = qstdout();
  p.n++;

  int len = list_len(p.list);
  double ** field = (double **) matrix_new (p.n, p.n, len*sizeof(double));

  double Delta = 0.999999*L0/(p.n - 1);
  for (int i = 0; i < p.n; i++) {
    double x = Delta*i + X0;
    for (int j = 0; j < p.n; j++) {
      double y = Delta*j + Y0;
      if (p.linear) {
 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i64 = p.list; ((scalar *)&s)->i >= 0; s = *++_i64)
   field[i][len*j + k++] = interpolate ((struct _interpolate){s, x, y});
      }
      else {
 Point point = locate ((struct _locate){x, y});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 62 "/home/jiarongw/basilisk/src/output.h"

 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i65 = p.list; ((scalar *)&s)->i >= 0; s = *++_i65)
   field[i][len*j + k++] = point.level >= 0 ? val(s,0,0,0) : nodata;
      }
    }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    if (p.list) for (scalar s = *p.list, *_i66 = p.list; ((scalar *)&s)->i >= 0; s = *++_i66)
      fprintf (p.fp, " %d:%s", i++, _attribute[s.i].name);
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double x = Delta*i + X0;
      for (int j = 0; j < p.n; j++) {
 double y = Delta*j + Y0;

 fprintf (p.fp, "%g %g", x, y);
 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i67 = p.list; ((scalar *)&s)->i >= 0; s = *++_i67)
   fprintf (p.fp, " %g", field[i][len*j + k++]);
 fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
#if 1
  else
    MPI_Reduce (field[0], NULL, len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (field);
 end_trace("output_field", "/home/jiarongw/basilisk/src/output.h", 102); }
#line 130 "/home/jiarongw/basilisk/src/output.h"
struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};


void output_matrix (struct OutputMatrix p)
{ trace ("output_matrix", "/home/jiarongw/basilisk/src/output.h", 139);
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = qstdout();
  float fn = p.n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 157 "/home/jiarongw/basilisk/src/output.h"

 assert (point.level >= 0);
 v = val(p.f,0,0,0);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
 end_trace("output_matrix", "/home/jiarongw/basilisk/src/output.h", 165); }
#line 174 "/home/jiarongw/basilisk/src/output.h"
typedef void (* colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}





typedef struct {
  unsigned char r, g, b;
} color;

color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  color c;
  if (val == nodata) {
    c.r = c.g = c.b = 0;
    return c;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = 127 - 2, coef = 1.;
  else {
    i = val*(127 - 1);
    coef = val*(127 - 1) - i;
  }
  assert (i < 127 - 1);
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 311 "/home/jiarongw/basilisk/src/output.h"
static const char * extension (const char * file, const char * ext) {
  int len = strlen(file);
  return len > 4 && !strcmp (file + len - 4, ext) ? file + len - 4 : NULL;
}

static const char * is_animation (const char * file) {
  const char * ext;
  if ((ext = extension (file, ".mp4")) ||
      (ext = extension (file, ".ogv")) ||
      (ext = extension (file, ".gif")))
    return ext;
  return NULL;
}

static struct {
  FILE ** fp;
  char ** names;
  int n;
} open_image_data = {NULL, NULL, 0};

static void open_image_cleanup()
{
  for (int i = 0; i < open_image_data.n; i++) {
    qpclose (open_image_data.fp[i]);
    pfree (open_image_data.names[i],__func__,__FILE__,__LINE__);
  }
  pfree (open_image_data.fp,__func__,__FILE__,__LINE__);
  pfree (open_image_data.names,__func__,__FILE__,__LINE__);
  open_image_data.fp = NULL;
  open_image_data.names = NULL;
  open_image_data.n = 0;
}

static FILE * open_image_lookup (const char * file)
{
  for (int i = 0; i < open_image_data.n; i++)
    if (!strcmp (file, open_image_data.names[i]))
      return open_image_data.fp[i];
  return NULL;
}

static bool which (const char * command)
{
  char * s = getenv ("PATH");
  if (!s)
    return false;
  char path[strlen(s) + 1];
  strcpy (path, s);
  s = strtok (path, ":");
  while (s) {
    char f[strlen(s) + strlen(command) + 2];
    strcpy (f, s);
    strcat (f, "/");
    strcat (f, command);
    FILE * fp = fopen (f, "r");
    if (fp) {
      fclose (fp);
      return true;
    }
    s = strtok (NULL, ":");
  }
  return false;
}

static FILE * ppm_fallback (const char * file, const char * mode)
{
  char filename[strlen(file) + 5];
  strcpy (filename, file);
  strcat (filename, ".ppm");
  FILE * fp = fopen (filename, mode);
  if (!fp) {
    perror (file);

    MPI_Abort (MPI_COMM_WORLD, 1);

    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  assert (pid() == 0);
  const char * ext;
  if ((ext = is_animation (file))) {
    FILE * fp = open_image_lookup (file);
    if (fp)
      return fp;

    int len = strlen ("ppm2???    ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "ppm2"); strcat (command, ext + 1);

    static int has_ffmpeg = -1;
    if (has_ffmpeg < 0) {
      if (which (command) && (which ("ffmpeg") || which ("avconv")))
 has_ffmpeg = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find '%s' or 'ffmpeg'/'avconv'\n"
   "  falling back to raw PPM outputs\n", command);
 has_ffmpeg = false;
      }
    }
    if (!has_ffmpeg)
      return ppm_fallback (file, "a");

    static bool added = false;
    if (!added) {
      free_solver_func_add (open_image_cleanup);
      added = true;
    }
    open_image_data.n++;
    open_image_data.names = (char * *) prealloc (open_image_data.names, (open_image_data.n)*sizeof(char *),__func__,__FILE__,__LINE__);
    open_image_data.names[open_image_data.n - 1] = pstrdup (file,__func__,__FILE__,__LINE__);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    open_image_data.fp = (FILE * *) prealloc (open_image_data.fp, (open_image_data.n)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    return open_image_data.fp[open_image_data.n - 1] = qpopen (command, "w");
  }
  else {
    static int has_convert = -1;
    if (has_convert < 0) {
      if (which ("convert"))
 has_convert = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find 'convert'\n"
   "  falling back to raw PPM outputs\n");
 has_convert = false;
      }
    }
    if (!has_convert)
      return ppm_fallback (file, "w");

    int len = strlen ("convert ppm:-   ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "convert ppm:- ");
    if (options) {
      strcat (command, options);
      strcat (command, " ");
    }
    strcat (command, file);
    return qpopen (command, "w");
  }
}

void close_image (const char * file, FILE * fp)
{
  assert (pid() == 0);
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    qpclose (fp);
  else
    fclose (fp);
}
#line 541 "/home/jiarongw/basilisk/src/output.h"
struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
  char * opt;
};


void output_ppm (struct OutputPPM p)
{ trace ("output_ppm", "/home/jiarongw/basilisk/src/output.h", 556);

  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
    p.min = avg - spread; p.max = avg + spread;
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;

  color ** ppm = (color **) matrix_new (ny, p.n, sizeof(color));
  double cmap[127][3];
  p.map (cmap);
  OMP_PARALLEL() {
    OMP(omp for schedule(static))
      for (int j = 0; j < ny; j++) {
 double yp = Delta*j + p.box[0][1] + Delta/2.;
 for (int i = 0; i < p.n; i++) {
   double xp = Delta*i + p.box[0][0] + Delta/2., v;
   if (p.mask.i) {
     if (p.linear) {
       double m = interpolate ((struct _interpolate){p.mask, xp, yp, p.z});
       if (m < 0.)
  v = nodata;
       else
  v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
     }
     else {
       Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 594 "/home/jiarongw/basilisk/src/output.h"

       if (point.level < 0 || val(p.mask,0,0,0) < 0.)
  v = nodata;
       else
  v = val(p.f,0,0,0);
     }
   }
   else if (p.linear)
     v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
   else {
     Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 604 "/home/jiarongw/basilisk/src/output.h"

     v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
   }
   ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
 }
      }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif
    if (!p.fp) p.fp = qstdout();
    if (p.file)
      p.fp = open_image (p.file, p.opt);

    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);

    if (p.file)
      close_image (p.file, p.fp);
    else
      fflush (p.fp);
  }
#if 1
  else
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (ppm);
 end_trace("output_ppm", "/home/jiarongw/basilisk/src/output.h", 636); }
#line 668 "/home/jiarongw/basilisk/src/output.h"
struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};


void output_grd (struct OutputGRD p)
{ trace ("output_grd", "/home/jiarongw/basilisk/src/output.h", 679);

  if (!p.fp) p.fp = qstdout();
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;


  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) {
 if (p.linear) {
   double m = interpolate ((struct _interpolate){p.mask, xp, yp});
   if (m < 0.)
     v = nodata;
   else
     v = interpolate ((struct _interpolate){p.f, xp, yp});
 }
 else {
   Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 715 "/home/jiarongw/basilisk/src/output.h"

   if (point.level < 0 || val(p.mask,0,0,0) < 0.)
     v = nodata;
   else
     v = val(p.f,0,0,0);
 }
      }
      else if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 725 "/home/jiarongw/basilisk/src/output.h"

 v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
      }
      if (v == nodata)
 fprintf (p.fp, "-9999 ");
      else
 fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
 end_trace("output_grd", "/home/jiarongw/basilisk/src/output.h", 737); }
#line 764 "/home/jiarongw/basilisk/src/output.h"
struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,__LINE__);
  }
  char * name = pstrdup (input,__func__,__FILE__,__LINE__), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}


void output_gfs (struct OutputGfs p)
{ trace ("output_gfs", "/home/jiarongw/basilisk/src/output.h", 794);
  char * fname = p.file;

#if 1



  FILE * fp = p.fp;
  if (p.file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,__LINE__));
    snprintf (fname, 80, ".output-%ld", pid);
    p.fp = NULL;
  }
#endif

  bool opened = false;
  if (p.fp == NULL) {
    if (fname == NULL)
      p.fp = qstdout();
    else if (!(p.fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * list = p.list ? p.list : list_copy (all);

  restriction (list);
  fprintf (p.fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);




  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,__LINE__);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
 fprintf (p.fp, ",%s", name);
 pfree (name,__func__,__FILE__,__LINE__);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  fprintf (p.fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (p.fp, "  VariableTracerVOF f\n");
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

#if 1
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  if (list) for (scalar s = *list, *_i68 = list; ((scalar *)&s)->i >= 0; s = *++_i68)
    if (_attribute[s.i].name)
      cell_size += sizeof(double);
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
#endif



   { foreach_cell(){

#line 872 "/home/jiarongw/basilisk/src/output.h"
 {
#if 1
    if (is_local(cell))
#endif
    {
#if 1
      if (fseek (p.fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags =
 level == 0 ? 0 :



      child.x == -1 && child.y == -1 ? 0 :
 child.x == -1 && child.y == 1 ? 1 :
 child.x == 1 && child.y == -1 ? 2 :
 3;
#line 902 "/home/jiarongw/basilisk/src/output.h"
      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      if (list) for (scalar s = *list, *_i69 = list; ((scalar *)&s)->i >= 0; s = *++_i69)
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && val(s,0,0,0) != nodata ? val(s,0,0,0) : (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && val(s,0,0,0) != nodata ? - val(s,0,0,0) : (double) DBL_MAX;
     }





   }
   else
     a = is_local(cell) && val(s,0,0,0) != nodata ? val(s,0,0,0) : (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, p.fp);
 }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

#if 1
  delete (((scalar []){index,{-1}}));
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
#endif
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    pfree (list,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);

#if 1
  if (p.file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (fp == NULL)
 fp = qstdout();
      p.fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, p.fp)) > 0)
 fwrite (buffer, 1, l, fp);
      fflush (fp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,__LINE__);
  }
#endif
 end_trace("output_gfs", "/home/jiarongw/basilisk/src/output.h", 970); }
#line 990 "/home/jiarongw/basilisk/src/output.h"
struct Dump {
  char * file;
  scalar * list;
  FILE * fp;
};

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =

  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat (((scalar []){cm,{-1}}), NULL);
  if (lista) for (scalar s = *lista, *_i70 = lista; ((scalar *)&s)->i >= 0; s = *++_i70)
    if (!_attribute[s.i].face && !_attribute[s.i].nodump && s.i != cm.i)
      list = list_add (list, s);
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  if (list) for (scalar s = *list, *_i71 = list; ((scalar *)&s)->i >= 0; s = *++_i71) {
    unsigned len = strlen(_attribute[s.i].name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (_attribute[s.i].name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

#if !1

void dump (struct Dump p)
{ trace ("dump", "/home/jiarongw/basilisk/src/output.h", 1043);
  FILE * fp = p.fp;
  char * file = p.file;

  if (file && (fp = fopen (file, "w")) == NULL) {
    perror (file);
    exit (1);
  }
  assert (fp);

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size= new_scalar("size");
  scalar * list = list_concat (((scalar []){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };
  dump_header (fp, &header, list);

  subtree_size (size, false);

   { foreach_cell(){

#line 1062 "/home/jiarongw/basilisk/src/output.h"
 {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    if (list) for (scalar s = *list, *_i72 = list; ((scalar *)&s)->i >= 0; s = *++_i72)
      if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
 perror ("dump(): error while writing scalars");
 exit (1);
      }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  pfree (list,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);
 delete (((scalar []){size,{-1}}));  end_trace("dump", "/home/jiarongw/basilisk/src/output.h", 1080); }
#else

void dump (struct Dump p)
{ trace ("dump", "/home/jiarongw/basilisk/src/output.h", 1084);
  FILE * fp = p.fp;
  char * file = p.file;

  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  FILE * fh = fopen (file, "w");

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size= new_scalar("size");
  scalar * list = list_concat (((scalar []){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };







  if (pid() == 0)
    dump_header (fh, &header, list);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  if (list) for (scalar s = *list, *_i73 = list; ((scalar *)&s)->i >= 0; s = *++_i73)
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(_attribute[s.i].name);

  subtree_size (size, false);

   { foreach_cell(){

#line 1121 "/home/jiarongw/basilisk/src/output.h"
 {

    if (is_local(cell)) {
      long offset = sizeofheader + val(index,0,0,0)*cell_size;
      fseek (fh, offset, SEEK_SET);
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      if (list) for (scalar s = *list, *_i74 = list; ((scalar *)&s)->i >= 0; s = *++_i74)
 fwrite (&val(s,0,0,0), 1, sizeof(double), fh);
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  delete (((scalar []){index,{-1}}));

  pfree (list,__func__,__FILE__,__LINE__);
  fclose (fh);
 delete (((scalar []){size,{-1}}));  end_trace("dump", "/home/jiarongw/basilisk/src/output.h", 1139); }
#endif


bool restore (struct Dump p)
{ trace ("restore", "/home/jiarongw/basilisk/src/output.h", 1144);
  FILE * fp = p.fp;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    { bool _ret =  false; end_trace("restore", "/home/jiarongw/basilisk/src/output.h", 1148);  return _ret; }
  assert (fp);

  struct DumpHeader header;
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }


  init_grid (1);
   { foreach_cell(){

#line 1159 "/home/jiarongw/basilisk/src/output.h"
 {
    cell.pid = pid();
    cell.flags |= active;
  } } end_foreach_cell(); }
  ((Tree *)grid)->dirty = true;
#line 1184 "/home/jiarongw/basilisk/src/output.h"
  bool restore_all = (p.list == all);
  scalar * list = dump_list (p.list ? p.list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (list)) {
      fprintf (ferr,
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (list));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (ferr,
        "restore(): error: file version mismatch: "
        "%d (file) != %d (code)\n",
        header.version, dump_version);
      exit (1);
    }

    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting len\n");
 exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting s.name\n");
 exit (1);
      }
      name[len] = '\0';

      if (i > 0) {
 bool found = false;
 if (list) for (scalar s = *list, *_i75 = list; ((scalar *)&s)->i >= 0; s = *++_i75)
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (list,__func__,__FILE__,__LINE__);
    list = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin ((struct _origin){o[0], o[1], o[2]});
    size (o[3]);
  }
#line 1259 "/home/jiarongw/basilisk/src/output.h"
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector []){{fm.x,fm.y},{{-1},{-1}}});

  restore_mpi (fp, list);
#line 1287 "/home/jiarongw/basilisk/src/output.h"
  boundary (listm);

  scalar * other = NULL;
  if (all) for (scalar s = *all, *_i76 = all; ((scalar *)&s)->i >= 0; s = *++_i76)
    if (!list_lookup (list, s) && !list_lookup (listm, s))
      other = list_append (other, s);
  reset (other, 0.);
  pfree (other,__func__,__FILE__,__LINE__);

  pfree (list,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  { bool _ret =  true; end_trace("restore", "/home/jiarongw/basilisk/src/output.h", 1309);  return _ret; }
 end_trace("restore", "/home/jiarongw/basilisk/src/output.h", 1310); }
#line 288 "/home/jiarongw/basilisk/src/utils.h"
#line 12 "/home/jiarongw/basilisk/src/run.h"


void run (void)
{ trace ("run", "/home/jiarongw/basilisk/src/run.h", 15);
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {





    update_perf();
    iter = inext, t = tnext;
  }




  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
 end_trace("run", "/home/jiarongw/basilisk/src/run.h", 37); }
#line 27 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
#line 1 "./timestep.h"
#line 1 "/home/jiarongw/basilisk/src/timestep.h"

double timestep (const vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _dtmax = dtmax; 
#line 6

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 6
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 6
{

#line 6 "/home/jiarongw/basilisk/src/timestep.h"

    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.x,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 6
{

#line 6 "/home/jiarongw/basilisk/src/timestep.h"

    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.y,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  end_foreach_face_generic()
#line 10
 end_foreach_face(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 6
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 6
{

#line 6 "/home/jiarongw/basilisk/src/timestep.h"

    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.x,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 6
{

#line 6 "/home/jiarongw/basilisk/src/timestep.h"

    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.y,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  end_foreach_face_generic()
#line 10
 end_foreach_face(); }OMP(omp critical) if (_dtmax < dtmax) dtmax = _dtmax;
mpi_all_reduce_double (dtmax, MPI_MIN);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 10
 }
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
#line 28 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
#line 1 "./bcg.h"
#line 1 "/home/jiarongw/basilisk/src/bcg.h"
#line 11 "/home/jiarongw/basilisk/src/bcg.h"
void tracer_fluxes (scalar f,
      vector uf,
      vector flux,
      double dt,
       scalar src)
{





  vector g= new_vector("g");
  gradients (((scalar []){f,{-1}}), ((vector []){{g.x,g.y},{{-1},{-1}}}));




   { 
if (!is_constant(fm.x) && !is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_src
#define val_src(a,i,j,k) val(a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) fine(a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) coarse(a,i,j,k)
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/jiarongw/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);







    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/jiarongw/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(4.*Delta);







    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 55
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(src)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_src
#define val_src(a,i,j,k) val(a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) fine(a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) coarse(a,i,j,k)
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/jiarongw/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);







    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/jiarongw/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(4.*Delta);







    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 55
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/jiarongw/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);







    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/jiarongw/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(4.*Delta);







    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 55
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(src)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/jiarongw/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);







    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/jiarongw/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(4.*Delta);







    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 55
 end_foreach_face(); } }





  boundary_flux (((vector []){{flux.x,flux.y},{{-1},{-1}}}));
 delete (((scalar []){g.x,g.y,{-1}})); }






struct Advection {
  scalar * tracers;
  vector u;
  double dt;
  scalar * src;
};

void advection (struct Advection p)
{




  scalar * lsrc = p.src;
  if (!lsrc) {
    scalar zero= new_const_scalar("zero", 6,  0.);
    if (p.tracers) for (scalar s = *p.tracers, *_i77 = p.tracers; ((scalar *)&s)->i >= 0; s = *++_i77)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, src;
  scalar * _i2 = p.tracers; scalar * _i3 = lsrc; if (p.tracers) for (f = *p.tracers, src = *lsrc; ((scalar *)&f)->i >= 0; f = *++_i2, src = *++_i3) {
    vector flux= new_face_vector("flux");
    tracer_fluxes (f, p.u, flux, p.dt, src);
     { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 94
foreach(){

#line 94 "/home/jiarongw/basilisk/src/bcg.h"

      {
#line 95

        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 95

        val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}; } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 94
foreach(){

#line 94 "/home/jiarongw/basilisk/src/bcg.h"

      {
#line 95

        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 95

        val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}; } end_foreach(); } }
   delete (((scalar []){flux.x,flux.y,{-1}})); }
  boundary (p.tracers);

  if (!p.src)
    pfree (lsrc,__func__,__FILE__,__LINE__);
}
#line 29 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
#line 1 "./viscosity.h"
#line 1 "/home/jiarongw/basilisk/src/viscosity.h"
#line 1 "./poisson.h"
#line 1 "/home/jiarongw/basilisk/src/poisson.h"
#line 32 "/home/jiarongw/basilisk/src/poisson.h"
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{




  restriction (res);





  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
       { foreach_level_or_leaf (l){

#line 54 "/home/jiarongw/basilisk/src/poisson.h"

 if (da) for (scalar s = *da, *_i78 = da; ((scalar *)&s)->i >= 0; s = *++_i78)
   val(s,0,0,0) = 0.; } end_foreach_level_or_leaf(); }





    else
       { foreach_level (l){

#line 63 "/home/jiarongw/basilisk/src/poisson.h"

 if (da) for (scalar s = *da, *_i79 = da; ((scalar *)&s)->i >= 0; s = *++_i79)
   val(s,0,0,0) = bilinear (point, s); } end_foreach_level(); }





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




   { foreach(){

#line 81 "/home/jiarongw/basilisk/src/poisson.h"
 {
    scalar s, ds;
    scalar * _i4 = a; scalar * _i5 = da; if (a) for (s = *a, ds = *da; ((scalar *)&s)->i >= 0; s = *++_i4, ds = *++_i5)
      val(s,0,0,0) += val(ds,0,0,0);
  } } end_foreach(); }
  boundary (a);
}
#line 99 "/home/jiarongw/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
} mgstats;
#line 120 "/home/jiarongw/basilisk/src/poisson.h"
struct MGSolve {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
         void * data);
  void (* relax) (scalar * da, scalar * res, int depth,
    void * data);
  void * data;

  int nrelax;
  scalar * res;
};

mgstats mg_solve (struct MGSolve p)
{





  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    if (p.a) for (scalar s = *p.a, *_i80 = p.a; ((scalar *)&s)->i >= 0; s = *++_i80) {
      scalar r = new_scalar("r");
      res = list_append (res, r);
    }






  for (int b = 0; b < nboundary; b++)
    if (da) for (scalar s = *da, *_i81 = da; ((scalar *)&s)->i >= 0; s = *++_i81)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];




  mgstats s = {0};
  double sum = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _sum = sum; 
#line 160
foreach (){

#line 160 "/home/jiarongw/basilisk/src/poisson.h"

    if (p.b) for (scalar s = *p.b, *_i82 = p.b; ((scalar *)&s)->i >= 0; s = *++_i82)
      _sum += val(s,0,0,0); } end_foreach();OMP(omp critical) sum += _sum;
mpi_all_reduce_double (sum, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 162
 }
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;




  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);






  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > TOLERANCE);
       s.i++) {
    mg_cycle (p.a, res, da, p.relax, p.data, s.nrelax, 0, grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);







    if (s.resa > TOLERANCE) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
 s.nrelax--;
    }
    resb = s.resa;
  }




  if (s.resa > TOLERANCE)
    fprintf (ferr,
      "WARNING: convergence not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d\n",
      s.i, s.resa, s.sum, s.nrelax), fflush (ferr);




  if (!p.res)
    delete (res), pfree (res,__func__,__FILE__,__LINE__);
  delete (da), pfree (da,__func__,__FILE__,__LINE__);

  return s;
}
#line 237 "/home/jiarongw/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
   vector alpha;
   scalar lambda;
  double tolerance;
  int nrelax;
  scalar * res;
};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;
#line 272 "/home/jiarongw/basilisk/src/poisson.h"
  scalar c = a;






   { 
if (!is_constant(lambda) && !is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 279
foreach_level_or_leaf (l){

#line 279 "/home/jiarongw/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 281
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 281
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (is_constant(lambda) && !is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 279
foreach_level_or_leaf (l){

#line 279 "/home/jiarongw/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 281
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 281
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (!is_constant(lambda) && is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 279
foreach_level_or_leaf (l){

#line 279 "/home/jiarongw/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 281
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 281
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (is_constant(lambda) && is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 279
foreach_level_or_leaf (l){

#line 279 "/home/jiarongw/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 281
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 281
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); } }
#line 304 "/home/jiarongw/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;
  double maxres = 0.;


  vector g= new_face_vector("g");
   { 
if (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 321
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 321
{

#line 321 "/home/jiarongw/basilisk/src/poisson.h"

    val(g.x,0,0,0) = val_alpha_x(alpha.x,0,0,0)*(val(a,0,0,0) - val(a,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 321
{

#line 321 "/home/jiarongw/basilisk/src/poisson.h"

    val(g.y,0,0,0) = val_alpha_y(alpha.y,0,0,0)*(val(a,0,0,0) - val(a,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 322
 end_foreach_face(); }
if (is_constant(alpha.x)) {
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 321
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 321
{

#line 321 "/home/jiarongw/basilisk/src/poisson.h"

    val(g.x,0,0,0) = val_alpha_x(alpha.x,0,0,0)*(val(a,0,0,0) - val(a,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 321
{

#line 321 "/home/jiarongw/basilisk/src/poisson.h"

    val(g.y,0,0,0) = val_alpha_y(alpha.y,0,0,0)*(val(a,0,0,0) - val(a,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 322
 end_foreach_face(); } }
  boundary_flux (((vector []){{g.x,g.y},{{-1},{-1}}}));
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 324

if (!is_constant(lambda)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#line 324
foreach (){

#line 324 "/home/jiarongw/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 326

      val(res,0,0,0) += (val(g.x,0,0,0) - val(g.x,1,0,0))/Delta;
#line 326

      val(res,0,0,0) += (val(g.y,0,0,0) - val(g.y,0,1,0))/Delta;}
    if (fabs (val(res,0,0,0)) > _maxres)
      _maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }
if (is_constant(lambda)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#line 324
foreach (){

#line 324 "/home/jiarongw/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 326

      val(res,0,0,0) += (val(g.x,0,0,0) - val(g.x,1,0,0))/Delta;
#line 326

      val(res,0,0,0) += (val(g.y,0,0,0) - val(g.y,0,1,0))/Delta;}
    if (fabs (val(res,0,0,0)) > _maxres)
      _maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 330
 }
#line 342 "/home/jiarongw/basilisk/src/poisson.h"
  boundary (resl);
  { double _ret =  maxres; delete (((scalar []){g.x,g.y,{-1}}));  return _ret; }
 delete (((scalar []){g.x,g.y,{-1}})); }
#line 355 "/home/jiarongw/basilisk/src/poisson.h"
mgstats poisson (struct Poisson p)
{






  if (!p.alpha.x.i) {
    vector alpha= new_const_vector("alpha", 7, (double []) {1.,1.,1.});
    p.alpha = alpha;
  }
  if (!p.lambda.i) {
    scalar lambda= new_const_scalar("lambda", 9,  0.);
    p.lambda = lambda;
  }




  vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction (((scalar []){alpha.x,alpha.y,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve ((struct MGSolve){((scalar []){a,{-1}}), ((scalar []){b,{-1}}), residual, relax, &p, p.nrelax, p.res});




  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}
#line 416 "/home/jiarongw/basilisk/src/poisson.h"
struct Project {
  vector u;
  scalar p;
  vector alpha;
  double dt;
  int nrelax;
};


mgstats project (struct Project q)
{ trace ("project", "/home/jiarongw/basilisk/src/poisson.h", 426);
  vector u = q.u;
  scalar p = q.p;
   vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;






  scalar div= new_scalar("div");
   { foreach(){

#line 439 "/home/jiarongw/basilisk/src/poisson.h"
 {
    val(div,0,0,0) = 0.;
    {
#line 441

      val(div,0,0,0) += val(u.x,1,0,0) - val(u.x,0,0,0);
#line 441

      val(div,0,0,0) += val(u.y,0,1,0) - val(u.y,0,0,0);}
    val(div,0,0,0) /= dt*Delta;
  } } end_foreach(); }
#line 455 "/home/jiarongw/basilisk/src/poisson.h"
  mgstats mgp = poisson ((struct Poisson){p, div, alpha,
    .tolerance = TOLERANCE/sq(dt), .nrelax = nrelax});




   { 
if (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 461
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 461
{

#line 461 "/home/jiarongw/basilisk/src/poisson.h"

    val(u.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 461
{

#line 461 "/home/jiarongw/basilisk/src/poisson.h"

    val(u.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 462
 end_foreach_face(); }
if (is_constant(alpha.x)) {
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 461
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 461
{

#line 461 "/home/jiarongw/basilisk/src/poisson.h"

    val(u.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 461
{

#line 461 "/home/jiarongw/basilisk/src/poisson.h"

    val(u.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 462
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}));

  { mgstats _ret =  mgp; delete (((scalar []){div,{-1}}));  end_trace("project", "/home/jiarongw/basilisk/src/poisson.h", 465);  return _ret; }
 delete (((scalar []){div,{-1}}));  end_trace("project", "/home/jiarongw/basilisk/src/poisson.h", 466); }
#line 2 "/home/jiarongw/basilisk/src/viscosity.h"

struct Viscosity {
  vector u;
  vector mu;
  scalar rho;
  double dt;
  int nrelax;
  scalar * res;
};
#line 25 "/home/jiarongw/basilisk/src/viscosity.h"
static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
   vector mu = p->mu;
   scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0])));




  vector w = u;


   { 
if (!is_constant(rho) && !is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/jiarongw/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
#line 58 "/home/jiarongw/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)




        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
#line 58 "/home/jiarongw/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)




        ));}
  } } end_foreach_level_or_leaf(); }
if (is_constant(rho) && !is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/jiarongw/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
#line 58 "/home/jiarongw/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)




        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
#line 58 "/home/jiarongw/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)




        ));}
  } } end_foreach_level_or_leaf(); }
if (!is_constant(rho) && is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/jiarongw/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
#line 58 "/home/jiarongw/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)




        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
#line 58 "/home/jiarongw/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)




        ));}
  } } end_foreach_level_or_leaf(); }
if (is_constant(rho) && is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
const struct { double x, y; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/jiarongw/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
#line 58 "/home/jiarongw/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)




        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
#line 58 "/home/jiarongw/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)




        ));}
  } } end_foreach_level_or_leaf(); } }
#line 85 "/home/jiarongw/basilisk/src/viscosity.h"
}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl,
      void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
   vector mu = p->mu;
   scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0]))), res = (*((vector *)&(resl[0])));
  double maxres = 0.;


  {
#line 98
 {
    vector taux= new_face_vector("taux");
     { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 100
{

#line 100 "/home/jiarongw/basilisk/src/viscosity.h"

      val(taux.x,0,0,0) = 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); }
if (is_constant(mu.x)) {
const struct { double x, y; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 100
{

#line 100 "/home/jiarongw/basilisk/src/viscosity.h"

      val(taux.x,0,0,0) = 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); } }

       { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 103
{

#line 103 "/home/jiarongw/basilisk/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); }
if (is_constant(mu.x)) {
const struct { double x, y; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 103
{

#line 103 "/home/jiarongw/basilisk/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); } }







    boundary_flux (((vector []){{taux.x,taux.y},{{-1},{-1}}}));
     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 115

if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#line 115
foreach (){

#line 115 "/home/jiarongw/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);}
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.}).x*val(u.x,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > _maxres)
 _maxres = fabs (val(res.x,0,0,0));
    } } end_foreach(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 115
foreach (){

#line 115 "/home/jiarongw/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);}
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.}).x*val(u.x,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > _maxres)
 _maxres = fabs (val(res.x,0,0,0));
    } } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 122
 }
   delete (((scalar []){taux.x,taux.y,{-1}})); }
#line 98
 {
    vector taux= new_face_vector("taux");
     { 
if (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,j,i,k) val(a,j,i,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,i,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,i,k)
#undef val_mu_x
#define val_mu_x(a,j,i,k) val(a,j,i,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,i,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,i,k)
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 100
{

#line 100 "/home/jiarongw/basilisk/src/viscosity.h"

      val(taux.y,0,0,0) = 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); }
if (is_constant(mu.y)) {
const struct { double x, y; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,j,i,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_x
#define val_mu_x(a,j,i,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 100
{

#line 100 "/home/jiarongw/basilisk/src/viscosity.h"

      val(taux.y,0,0,0) = 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); } }

       { 
if (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,j,i,k) val(a,j,i,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,i,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,i,k)
#undef val_mu_x
#define val_mu_x(a,j,i,k) val(a,j,i,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,i,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,i,k)
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_x()) {
#line 103
{

#line 103 "/home/jiarongw/basilisk/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); }
if (is_constant(mu.y)) {
const struct { double x, y; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,j,i,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_x
#define val_mu_x(a,j,i,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_x()) {
#line 103
{

#line 103 "/home/jiarongw/basilisk/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); } }







    boundary_flux (((vector []){{taux.x,taux.y},{{-1},{-1}}}));
     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 115

if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,j,i,k) val(a,j,i,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,j,i,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,j,i,k)
#line 115
foreach (){

#line 115 "/home/jiarongw/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);}
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.}).y*val(u.y,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > _maxres)
 _maxres = fabs (val(res.y,0,0,0));
    } } end_foreach(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,j,i,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 115
foreach (){

#line 115 "/home/jiarongw/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);}
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.}).y*val(u.y,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > _maxres)
 _maxres = fabs (val(res.y,0,0,0));
    } } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 122
 }
   delete (((scalar []){taux.x,taux.y,{-1}})); }}
  boundary (resl);
#line 153 "/home/jiarongw/basilisk/src/viscosity.h"
  return maxres;
}



mgstats viscosity (struct Viscosity p)
{
  vector u = p.u, r= new_vector("r");
   { foreach(){

#line 161 "/home/jiarongw/basilisk/src/viscosity.h"

    {
#line 162

      val(r.x,0,0,0) = val(u.x,0,0,0);
#line 162

      val(r.y,0,0,0) = val(u.y,0,0,0);}; } end_foreach(); }

  vector mu = p.mu;
  scalar rho = p.rho;
  restriction (((scalar []){mu.x,mu.y,rho,{-1}}));

  { mgstats _ret =  mg_solve ((struct MGSolve){(scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}), (scalar *)((vector []){{r.x,r.y},{{-1},{-1}}}),
     residual_viscosity, relax_viscosity, &p, p.nrelax, p.res}); delete (((scalar []){r.x,r.y,{-1}}));  return _ret; }
 delete (((scalar []){r.x,r.y,{-1}})); }

mgstats viscosity_explicit (struct Viscosity p)
{
  vector u = p.u, r= new_vector("r");
  mgstats mg = {0};
  mg.resb = residual_viscosity ((scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}), (scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}), (scalar *)((vector []){{r.x,r.y},{{-1},{-1}}}), &p);
   { foreach(){

#line 178 "/home/jiarongw/basilisk/src/viscosity.h"

    {
#line 179

      val(u.x,0,0,0) += val(r.x,0,0,0);
#line 179

      val(u.y,0,0,0) += val(r.y,0,0,0);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}));
  { mgstats _ret =  mg; delete (((scalar []){r.x,r.y,{-1}}));  return _ret; }
 delete (((scalar []){r.x,r.y,{-1}})); }
#line 30 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
#line 39 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
scalar p= {0};
vector u= {{1},{2}}, g= {{3},{4}};
scalar pf= {5};
vector uf= {{6},{7}};
#line 65 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
 vector mu = {{_NVARMAX + 0},{_NVARMAX + 1}}, a = {{_NVARMAX + 0},{_NVARMAX + 1}}, alpha = {{_NVARMAX + 2},{_NVARMAX + 3}};
 scalar rho = {(_NVARMAX + 4)};
mgstats mgp, mgpf, mgu;
bool stokes = false;
#line 79 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
static void _set_boundary0 (void) { _attribute[p.i].boundary[right] = _boundary0; _attribute[p.i].boundary_homogeneous[right] = _boundary0_homogeneous; } 
static void _set_boundary1 (void) { _attribute[p.i].boundary[left] = _boundary1; _attribute[p.i].boundary_homogeneous[left] = _boundary1_homogeneous; } 







static void _set_boundary2 (void) { _attribute[p.i].boundary[top] = _boundary2; _attribute[p.i].boundary_homogeneous[top] = _boundary2_homogeneous; } 
static void _set_boundary3 (void) { _attribute[p.i].boundary[bottom] = _boundary3; _attribute[p.i].boundary_homogeneous[bottom] = _boundary3_homogeneous; } 
#line 100 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
static int defaults_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults (const int i, const double t, Event * _ev) { trace ("defaults", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 100); 
{

  CFL = 0.8;




  _attribute[p.i].nodump = _attribute[pf.i].nodump = true;




  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    vector alphav = alpha;
     { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 119
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 119
{

#line 119 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

      val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 119
{

#line 119 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

      val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0); }  }}  end_foreach_face_generic()
#line 120
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 119
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 119
{

#line 119 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

      val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 119
{

#line 119 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

      val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0); }  }}  end_foreach_face_generic()
#line 120
 end_foreach_face(); } }
    boundary ((scalar *)((vector []){{alpha.x,alpha.y},{{-1},{-1}}}));
  }






  _attribute[uf.x.i].refine = refine_face_solenoidal;

 end_trace("defaults", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 131); } return 0; } 





double dtmax;

static int init_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init (const int i, const double t, Event * _ev) { trace ("init", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 139); 
{
  boundary ((scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}));
  trash (((vector []){{uf.x,uf.y},{{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 143
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 143
{

#line 143 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 143
{

#line 143 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.; }  }}  end_foreach_face_generic()
#line 144
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 143
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 143
{

#line 143 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 143
{

#line 143 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.; }  }}  end_foreach_face_generic()
#line 144
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y},{{-1},{-1}}}));




  event ("properties");





  dtmax = DT;
  event ("stability");
 end_trace("init", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 158); } return 0; } 
#line 167 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int set_dtmax (const int i, const double t, Event * _ev) { trace ("set_dtmax", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 167);  dtmax = DT; end_trace("set_dtmax", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 167);  return 0; } 

static int stability_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability (const int i, const double t, Event * _ev) { trace ("stability", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 169);  {
  dt = dtnext (timestep (uf, dtmax));
 end_trace("stability", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 171); } return 0; } 







static int vof_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof (const int i, const double t, Event * _ev) { trace ("vof", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 179); ; end_trace("vof", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 179);  return 0; } 
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_advection (const int i, const double t, Event * _ev) { trace ("tracer_advection", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 180); ; end_trace("tracer_advection", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 180);  return 0; } 
static int tracer_diffusion_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_diffusion (const int i, const double t, Event * _ev) { trace ("tracer_diffusion", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 181); ; end_trace("tracer_diffusion", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 181);  return 0; } 






static int properties_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int properties (const int i, const double t, Event * _ev) { trace ("properties", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 188);  {
  boundary (((scalar []){alpha.x,alpha.y,mu.x,mu.y,rho,{-1}}));
 end_trace("properties", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 190); } return 0; } 
#line 202 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
void prediction()
{
  vector du;
  {
#line 205
 {
    scalar s = new_scalar("s");
    du.x = s;
  }
#line 205
 {
    scalar s = new_scalar("s");
    du.y = s;
  }}

  if (_attribute[u.x.i].gradient)
     { foreach(){

#line 211 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

      {
#line 212

        val(du.x,0,0,0) = _attribute[u.x.i].gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;
#line 212

        val(du.y,0,0,0) = _attribute[u.y.i].gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;}; } end_foreach(); }
  else
     { foreach(){

#line 215 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

      {
#line 216

        val(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);
#line 216

        val(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{du.x,du.y},{{-1},{-1}}}));

  trash (((vector []){{uf.x,uf.y},{{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 221
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 221
{

#line 221 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);





    val(uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 221
{

#line 221 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

      double fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);





    val(uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 234
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 221
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 221
{

#line 221 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);





    val(uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 221
{

#line 221 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

      double fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);





    val(uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 234
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y},{{-1},{-1}}}));

  delete ((scalar *)((vector []){{du.x,du.y},{{-1},{-1}}}));
}
#line 249 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
static int advection_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int advection_term (const int i, const double t, Event * _ev) { trace ("advection_term", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 249); 
{
  if (!stokes) {
    prediction();
    mgpf = project ((struct Project){uf, pf, alpha, dt/2., mgpf.nrelax});
    advection ((struct Advection){(scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}), uf, dt, (scalar *)((vector []){{g.x,g.y},{{-1},{-1}}})});
  }
 end_trace("advection_term", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 256); } return 0; } 







static void correction (double dt)
{
   { foreach(){

#line 266 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    {
#line 267

      val(u.x,0,0,0) += dt*val(g.x,0,0,0);
#line 267

      val(u.y,0,0,0) += dt*val(g.y,0,0,0);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}));
}
#line 279 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int viscous_term (const int i, const double t, Event * _ev) { trace ("viscous_term", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 279); 
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity ((struct Viscosity){u, mu, rho, dt, mgu.nrelax});
    correction (-dt);
  }






  vector af = a;
  trash (((vector []){{uf.x,uf.y},{af.x,af.y},{{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 294
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 294
{

#line 294 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.;
    if (!is_constant(af.x))
      val(af.x,0,0,0) = 0.;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 294
{

#line 294 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.;
    if (!is_constant(af.y))
      val(af.y,0,0,0) = 0.;
  } }  }}  end_foreach_face_generic()
#line 298
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 294
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 294
{

#line 294 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.;
    if (!is_constant(af.x))
      val(af.x,0,0,0) = 0.;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 294
{

#line 294 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.;
    if (!is_constant(af.y))
      val(af.y,0,0,0) = 0.;
  } }  }}  end_foreach_face_generic()
#line 298
 end_foreach_face(); } }
 end_trace("viscous_term", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 299); } return 0; } 
#line 314 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
static int acceleration_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration (const int i, const double t, Event * _ev) { trace ("acceleration", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 314); 
{
  boundary ((scalar *)((vector []){{a.x,a.y},{{-1},{-1}}}));
   { 
if (!is_constant(fm.x) && !is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#line 317
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 317
{

#line 317 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 317
{

#line 317 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  end_foreach_face_generic()
#line 318
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(a.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#line 317
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 317
{

#line 317 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 317
{

#line 317 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  end_foreach_face_generic()
#line 318
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#line 317
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 317
{

#line 317 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 317
{

#line 317 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  end_foreach_face_generic()
#line 318
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(a.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#line 317
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 317
{

#line 317 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 317
{

#line 317 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  end_foreach_face_generic()
#line 318
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y},{{-1},{-1}}}));
 end_trace("acceleration", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 320); } return 0; } 
#line 329 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
static int projection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int projection (const int i, const double t, Event * _ev) { trace ("projection", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 329); 
{
  mgp = project ((struct Project){uf, p, alpha, dt, mgp.nrelax});





  vector gf= new_face_vector("gf");
   { 
if (!is_constant(a.x) && !is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (is_constant(a.x) && !is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (!is_constant(a.x) && is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (is_constant(a.x) && is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (!is_constant(a.x) && !is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (is_constant(a.x) && !is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (!is_constant(a.x) && is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (is_constant(a.x) && is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); } }
  boundary_flux (((vector []){{gf.x,gf.y},{{-1},{-1}}}));





  trash (((vector []){{g.x,g.y},{{-1},{-1}}}));
   { foreach(){

#line 347 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

    {
#line 348

      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/2.;
#line 348

      val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/2.;}; } end_foreach(); }
  boundary ((scalar *)((vector []){{g.x,g.y},{{-1},{-1}}}));




  correction (dt);
 delete (((scalar []){gf.x,gf.y,{-1}}));  end_trace("projection", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 356); } return 0; } 







static int adapt_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int adapt (const int i, const double t, Event * _ev) { trace ("adapt", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 364);  {
  event ("properties");
 end_trace("adapt", "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 366); } return 0; } 
#line 10 "wavewind.c"
#line 1 "two-phase.h"
#line 1 "/home/jiarongw/basilisk/src/two-phase.h"
#line 13 "/home/jiarongw/basilisk/src/two-phase.h"
#line 1 "vof.h"
#line 1 "/home/jiarongw/basilisk/src/vof.h"
#line 26 "/home/jiarongw/basilisk/src/vof.h"








#line 1 "fractions.h"
#line 1 "/home/jiarongw/basilisk/src/fractions.h"
#line 11 "/home/jiarongw/basilisk/src/fractions.h"
#line 1 "geometry.h"
#line 1 "/home/jiarongw/basilisk/src/geometry.h"
#line 28 "/home/jiarongw/basilisk/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    swap (double, n1, n2);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}
#line 133 "/home/jiarongw/basilisk/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}
#line 237 "/home/jiarongw/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
  {
#line 240
 {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  }
#line 240
 {
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  }}
  return line_area(n1.x, n1.y, alpha);
}
#line 262 "/home/jiarongw/basilisk/src/geometry.h"
int facets (coord n, double alpha, coord p[2])
{
  int i = 0;
  for (double s = -0.5; s <= 0.5; s += 1.)
    {
#line 266

      if (fabs (n.y) > 1e-4 && i < 2) {
 double a = (alpha - s*n.x)/n.y;
 if (a >= -0.5 && a <= 0.5) {
   p[i].x = s;
   p[i++].y = a;
 }
      }
#line 266

      if (fabs (n.x) > 1e-4 && i < 2) {
 double a = (alpha - s*n.y)/n.x;
 if (a >= -0.5 && a <= 0.5) {
   p[i].y = s;
   p[i++].x = a;
 }
      }}
  return i;
}
#line 352 "/home/jiarongw/basilisk/src/geometry.h"
double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p->x = p->y = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  {
#line 371

    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
#line 371

    if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }}

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

  {
#line 397
 {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  }
#line 397
 {
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }}

  return sqrt (ax*ax + ay*ay);
}
#line 12 "/home/jiarongw/basilisk/src/fractions.h"






#line 1 "myc2d.h"
#line 1 "/home/jiarongw/basilisk/src/myc2d.h"





coord mycs (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 7 "/home/jiarongw/basilisk/src/myc2d.h"

  int ix;
  double c_t,c_b,c_r,c_l;
  double mx0,my0,mx1,my1,mm1,mm2;


  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);
  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);
  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);
  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);



  mx0 = 0.5*(c_l-c_r);
  my0 = 0.5*(c_b-c_t);


  if (fabs(mx0) <= fabs(my0)) {
    my0 = my0 > 0. ? 1. : -1.;
    ix = 1;
  }
  else {
    mx0 = mx0 > 0. ? 1. : -1.;
    ix = 0;
  }


  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);
  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);
  mx1 = mm1 - mm2 + 1.e-30;
  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);
  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);
  my1 = mm1 - mm2 + 1.e-30;


  if (ix) {
    mm1 = fabs(my1);
    mm1 = fabs(mx1)/mm1;
    if (mm1 > fabs(mx0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
  else {
    mm1 = fabs(mx1);
    mm1 = fabs(my1)/mm1;
    if (mm1 > fabs(my0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }



  mm1 = fabs(mx0) + fabs(my0);
  coord n = {mx0/mm1, my0/mm1};

  return n;
}
#line 19 "/home/jiarongw/basilisk/src/fractions.h"
#line 32 "/home/jiarongw/basilisk/src/fractions.h"
void fraction_refine (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 33 "/home/jiarongw/basilisk/src/fractions.h"






  double cc = val(c,0,0,0);
  if (cc <= 0. || cc >= 1.)
     { foreach_child()
      val(c,0,0,0) = cc; end_foreach_child(); }
  else {




    coord n = mycs (point, c);
    double alpha = line_alpha (cc, n);






    coord a, b;
    {
#line 57
 {
      a.x = 0.; b.x = 0.5;
    }
#line 57
 {
      a.y = 0.; b.y = 0.5;
    }}

     { foreach_child() {
      coord nc;
      {
#line 63

 nc.x = child.x*n.x;
#line 63

 nc.y = child.y*n.y;}
      val(c,0,0,0) = rectangle_fraction (nc, alpha, a, b);
    } end_foreach_child(); }
  }
}











static void alpha_refine (Point point, scalar alpha)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 81 "/home/jiarongw/basilisk/src/fractions.h"

  vector n = _attribute[alpha.i].n;
  double alphac = 2.*val(alpha,0,0,0);
  coord m;
  {
#line 85

    m.x = val(n.x,0,0,0);
#line 85

    m.y = val(n.y,0,0,0);}
   { foreach_child() {
    val(alpha,0,0,0) = alphac;
    {
#line 89

      val(alpha,0,0,0) -= child.x*m.x/2.;
#line 89

      val(alpha,0,0,0) -= child.y*m.y/2.;}
  } end_foreach_child(); }
}
#line 116 "/home/jiarongw/basilisk/src/fractions.h"
struct Fractions {
  scalar Phi;
  scalar c;
  vector s;
};


void fractions (struct Fractions a)
{ trace ("fractions", "/home/jiarongw/basilisk/src/fractions.h", 124);
  scalar Phi = a.Phi;
  scalar c = a.c;
  vector s = (a.s).x.i ? (a.s) : new_face_vector("s");
#line 138 "/home/jiarongw/basilisk/src/fractions.h"
  vector p;
  p.x = s.y; p.y = s.x;
#line 148 "/home/jiarongw/basilisk/src/fractions.h"
   { foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 148
{

#line 148 "/home/jiarongw/basilisk/src/fractions.h"
 {





    if (val(Phi,0,0,0)*val(Phi,1,0,0) < 0.) {






      val(p.x,0,0,0) = val(Phi,0,0,0)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < 0.)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 173 "/home/jiarongw/basilisk/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > 0. || val(Phi,1,0,0) > 0.);
  } }  }}  { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 148
{

#line 148 "/home/jiarongw/basilisk/src/fractions.h"
 {





    if (val(Phi,0,0,0)*val(Phi,0,1,0) < 0.) {






      val(p.y,0,0,0) = val(Phi,0,0,0)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < 0.)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 173 "/home/jiarongw/basilisk/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > 0. || val(Phi,0,1,0) > 0.);
  } }  }}  end_foreach_face_generic()
#line 175
 end_foreach_face(); }
#line 190 "/home/jiarongw/basilisk/src/fractions.h"
  boundary_flux (((vector []){{s.x,s.y},{{-1},{-1}}}));
  scalar s_z = c;
   { foreach(){

#line 192 "/home/jiarongw/basilisk/src/fractions.h"


  {
#line 226 "/home/jiarongw/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 228
 {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    }
#line 228
 {
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }}





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      {
#line 245

 n.x /= nn;
#line 245

 n.y /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 255

   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0))*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
#line 255

   if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0))*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 269 "/home/jiarongw/basilisk/src/fractions.h"
      val(s_z,0,0,0) = ni ? line_area (n.x, n.y, alpha/ni) : max (val(p.x,0,0,0), val(p.y,0,0,0));
    }
  } } end_foreach(); }
#line 323 "/home/jiarongw/basilisk/src/fractions.h"
  boundary (((scalar []){c,{-1}}));
 { if (!(a.s).x.i) delete (((scalar []){s.x,s.y,{-1}})); }  end_trace("fractions", "/home/jiarongw/basilisk/src/fractions.h", 324); }
#line 349 "/home/jiarongw/basilisk/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 350 "/home/jiarongw/basilisk/src/fractions.h"

  coord n;
  double nn = 0.;
  assert (2 == 2);
  {
#line 354
 {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  }
#line 354
 {
    n.y = (val(c,1,-1,0) + 2.*val(c,0,-1,0) + val(c,-1,-1,0) -
    val(c,1,+1,0) - 2.*val(c,0,+1,0) - val(c,-1,+1,0));
    nn += fabs(n.y);
  }}

  if (nn > 0.)
    {
#line 361

      n.x /= nn;
#line 361

      n.y /= nn;}
  else
    n.x = 1.;
  return n;
}





coord facet_normal (Point point, scalar c, vector s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 373 "/home/jiarongw/basilisk/src/fractions.h"

  coord n;
  if (s.x.i < 0)
    n = mycs (point, c);
  else {
    double nn = 0.;
    {
#line 379
 {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    }
#line 379
 {
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }}
    assert (nn > 0.);
    {
#line 384

      n.x /= nn;
#line 384

      n.y /= nn;}
  }
  return n;
}
#line 397 "/home/jiarongw/basilisk/src/fractions.h"

void reconstruction (const scalar c, vector n, scalar alpha)
{ trace ("reconstruction", "/home/jiarongw/basilisk/src/fractions.h", 399);
   { foreach(){

#line 400 "/home/jiarongw/basilisk/src/fractions.h"
 {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      {
#line 408

 val(n.x,0,0,0) = 0.;
#line 408

 val(n.y,0,0,0) = 0.;}
    }
    else {






      coord m = mycs (point, c);

      {
#line 420

 val(n.x,0,0,0) = m.x;
#line 420

 val(n.y,0,0,0) = m.y;}
      val(alpha,0,0,0) = line_alpha (val(c,0,0,0), m);
    }
  } } end_foreach(); }
#line 433 "/home/jiarongw/basilisk/src/fractions.h"
  {
#line 433

    _attribute[n.x.i].refine = _attribute[n.x.i].prolongation = refine_injection;
#line 433

    _attribute[n.y.i].refine = _attribute[n.y.i].prolongation = refine_injection;}




  _attribute[alpha.i].n = n;
  _attribute[alpha.i].refine = _attribute[alpha.i].prolongation = alpha_refine;







  boundary (((scalar []){n.x,n.y,alpha,{-1}}));
 end_trace("reconstruction", "/home/jiarongw/basilisk/src/fractions.h", 449); }
#line 469 "/home/jiarongw/basilisk/src/fractions.h"
struct OutputFacets {
  scalar c;
  FILE * fp;
  vector s;
};


void output_facets (struct OutputFacets p)
{ trace ("output_facets", "/home/jiarongw/basilisk/src/fractions.h", 477);
  scalar c = p.c;
  vector s = p.s;
  if (!p.fp) p.fp = qstdout();
  if (!s.x.i) s.x.i = -1;

   { foreach(){

#line 483 "/home/jiarongw/basilisk/src/fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = line_alpha (val(c,0,0,0), n);

      coord segment[2];
      if (facets (n, alpha, segment) == 2)
 fprintf (p.fp, "%g %g\n%g %g\n\n",
   x + segment[0].x*Delta, y + segment[0].y*Delta,
   x + segment[1].x*Delta, y + segment[1].y*Delta);
#line 502 "/home/jiarongw/basilisk/src/fractions.h"
    } } end_foreach(); }

  fflush (p.fp);
 end_trace("output_facets", "/home/jiarongw/basilisk/src/fractions.h", 505); }








double interface_area (scalar c)
{ trace ("interface_area", "/home/jiarongw/basilisk/src/fractions.h", 515);
  double area = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _area = area; 
#line 517
foreach (){

#line 517 "/home/jiarongw/basilisk/src/fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = line_alpha (val(c,0,0,0), n);
      _area += pow(Delta, 2 - 1)*line_length_center(n,alpha,&p);
    } } end_foreach();OMP(omp critical) area += _area;
mpi_all_reduce_double (area, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 522
 }
  { double _ret =  area; end_trace("interface_area", "/home/jiarongw/basilisk/src/fractions.h", 523);  return _ret; }
 end_trace("interface_area", "/home/jiarongw/basilisk/src/fractions.h", 524); }
#line 35 "/home/jiarongw/basilisk/src/vof.h"
#line 43 "/home/jiarongw/basilisk/src/vof.h"
extern scalar * interfaces;
extern vector uf;
extern double dt;





static int defaults_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_0 (const int i, const double t, Event * _ev) { trace ("defaults_0", "/home/jiarongw/basilisk/src/vof.h", 51); 
{

  if (interfaces) for (scalar c = *interfaces, *_i83 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i83)
    _attribute[c.i].refine = _attribute[c.i].prolongation = fraction_refine;

 end_trace("defaults_0", "/home/jiarongw/basilisk/src/vof.h", 57); } return 0; } 





static int stability_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability_0 (const int i, const double t, Event * _ev) { trace ("stability_0", "/home/jiarongw/basilisk/src/vof.h", 63);  {
  if (CFL > 0.5)
    CFL = 0.5;
 end_trace("stability_0", "/home/jiarongw/basilisk/src/vof.h", 66); } return 0; } 
#line 80 "/home/jiarongw/basilisk/src/vof.h"

#line 80

static void sweep_x (scalar c, scalar cc)
{
  vector n= new_vector("n");
  scalar alpha= new_scalar("alpha"), flux= new_scalar("flux");
  double cfl = 0.;
#line 94 "/home/jiarongw/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    if (tracers) for (scalar t = *tracers, *_i84 = tracers; ((scalar *)&t)->i >= 0; t = *++_i84) {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }







     { foreach(){

#line 108 "/home/jiarongw/basilisk/src/vof.h"
 {
      scalar t, gf;
      scalar * _i6 = tracers; scalar * _i7 = gfl; if (tracers) for (t = *tracers, gf = *gfl; ((scalar *)&t)->i >= 0; t = *++_i6, gf = *++_i7) {
 double cl = val(c,-1,0,0), cc = val(c,0,0,0), cr = val(c,1,0,0);
 if (_attribute[t.i].inverse)
   cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
 val(gf,0,0,0) = 0.;
 static const double cmin = 0.5;
 if (cc >= cmin && _attribute[t.i].gradient != zero) {
   if (cr >= cmin) {
     if (cl >= cmin) {
       if (_attribute[t.i].gradient)
  val(gf,0,0,0) = _attribute[t.i].gradient (val(t,-1,0,0)/cl, val(t,0,0,0)/cc, val(t,1,0,0)/cr)/Delta;
       else
  val(gf,0,0,0) = (val(t,1,0,0)/cr - val(t,-1,0,0)/cl)/(2.*Delta);
     }
     else
        val(gf,0,0,0) = (val(t,1,0,0)/cr - val(t,0,0,0)/cc)/Delta;
   }
   else if (cl >= cmin)
     val(gf,0,0,0) = (val(t,0,0,0)/cc - val(t,-1,0,0)/cl)/Delta;
 }
      }
    } } end_foreach(); }
    boundary (gfl);
  }






  reconstruction (c, n, alpha);

   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _cfl = cfl; 
#line 142

if (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 142
{

#line 142 "/home/jiarongw/basilisk/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/jiarongw/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension(n.x,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 142
{

#line 142 "/home/jiarongw/basilisk/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/jiarongw/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension(n.x,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 142
{

#line 142 "/home/jiarongw/basilisk/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/jiarongw/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension(n.x,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 142
{

#line 142 "/home/jiarongw/basilisk/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/jiarongw/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension(n.x,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }OMP(omp critical) if (_cfl > cfl) cfl = _cfl;
mpi_all_reduce_double (cfl, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 197
 }
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);
#line 208 "/home/jiarongw/basilisk/src/vof.h"
  scalar * fluxl = list_concat (NULL, tfluxl);
  fluxl = list_append (fluxl, flux);
  for (int l = depth() - 1; l >= 0; l--)
     { foreach_halo (prolongation, l){

#line 211 "/home/jiarongw/basilisk/src/vof.h"
 {
#line 220 "/home/jiarongw/basilisk/src/vof.h"
      if ((!is_leaf (neighbor(-1,0,)) && neighbor(-1,0,).neighbors && neighbor(-1,0,).pid >= 0))
 if (fluxl) for (scalar fl = *fluxl, *_i85 = fluxl; ((scalar *)&fl)->i >= 0; fl = *++_i85)
   val(fl,0,0,0) = (fine(fl,0,0,0) + fine(fl,0,1,0))/2.;
      if ((!is_leaf (neighbor(1,0,)) && neighbor(1,0,).neighbors && neighbor(1,0,).pid >= 0))
 if (fluxl) for (scalar fl = *fluxl, *_i86 = fluxl; ((scalar *)&fl)->i >= 0; fl = *++_i86)
   val(fl,1,0,0) = (fine(fl,2,0,0) + fine(fl,2,1,0))/2.;
#line 236 "/home/jiarongw/basilisk/src/vof.h"
    } } end_foreach_halo(); }
  pfree (fluxl,__func__,__FILE__,__LINE__);





  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 265 "/home/jiarongw/basilisk/src/vof.h"
   { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 265
foreach(){

#line 265 "/home/jiarongw/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tfluxl; if (tracers) for (t = *tracers, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tflux = *++_i12)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,1,0,0))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 265
foreach(){

#line 265 "/home/jiarongw/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tfluxl; if (tracers) for (t = *tracers, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tflux = *++_i12)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,1,0,0))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); } }
  boundary (((scalar []){c,{-1}}));
  boundary (tracers);

  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);
 delete (((scalar []){flux,alpha,n.x,n.y,{-1}})); }
#line 80

static void sweep_y (scalar c, scalar cc)
{
  vector n= new_vector("n");
  scalar alpha= new_scalar("alpha"), flux= new_scalar("flux");
  double cfl = 0.;
#line 94 "/home/jiarongw/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    if (tracers) for (scalar t = *tracers, *_i84 = tracers; ((scalar *)&t)->i >= 0; t = *++_i84) {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }







     { foreach(){

#line 108 "/home/jiarongw/basilisk/src/vof.h"
 {
      scalar t, gf;
      scalar * _i6 = tracers; scalar * _i7 = gfl; if (tracers) for (t = *tracers, gf = *gfl; ((scalar *)&t)->i >= 0; t = *++_i6, gf = *++_i7) {
 double cl = val(c,0,-1,0), cc = val(c,0,0,0), cr = val(c,0,1,0);
 if (_attribute[t.i].inverse)
   cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
 val(gf,0,0,0) = 0.;
 static const double cmin = 0.5;
 if (cc >= cmin && _attribute[t.i].gradient != zero) {
   if (cr >= cmin) {
     if (cl >= cmin) {
       if (_attribute[t.i].gradient)
  val(gf,0,0,0) = _attribute[t.i].gradient (val(t,0,-1,0)/cl, val(t,0,0,0)/cc, val(t,0,1,0)/cr)/Delta;
       else
  val(gf,0,0,0) = (val(t,0,1,0)/cr - val(t,0,-1,0)/cl)/(2.*Delta);
     }
     else
        val(gf,0,0,0) = (val(t,0,1,0)/cr - val(t,0,0,0)/cc)/Delta;
   }
   else if (cl >= cmin)
     val(gf,0,0,0) = (val(t,0,0,0)/cc - val(t,0,-1,0)/cl)/Delta;
 }
      }
    } } end_foreach(); }
    boundary (gfl);
  }






  reconstruction (c, n, alpha);

   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _cfl = cfl; 
#line 142

if (!is_constant(fm.y) && !is_constant(cm)) {
#undef val_fm_y
#define val_fm_y(a,j,i,k) val(a,j,i,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,j,i,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,j,i,k)
#undef val_fm_x
#define val_fm_x(a,j,i,k) val(a,j,i,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,j,i,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,j,i,k)
#undef val_cm
#define val_cm(a,j,i,k) val(a,j,i,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,j,i,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,j,i,k)
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 142
{

#line 142 "/home/jiarongw/basilisk/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/jiarongw/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension(n.y,i,0,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (is_constant(fm.y) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.y.i -_NVARMAX], _constant[fm.x.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_y
#define val_fm_y(a,j,i,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_x
#define val_fm_x(a,j,i,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_cm
#define val_cm(a,j,i,k) val(a,j,i,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,j,i,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,j,i,k)
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 142
{

#line 142 "/home/jiarongw/basilisk/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/jiarongw/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension(n.y,i,0,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (!is_constant(fm.y) && is_constant(cm)) {
#undef val_fm_y
#define val_fm_y(a,j,i,k) val(a,j,i,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,j,i,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,j,i,k)
#undef val_fm_x
#define val_fm_x(a,j,i,k) val(a,j,i,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,j,i,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,j,i,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,i,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 142
{

#line 142 "/home/jiarongw/basilisk/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/jiarongw/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension(n.y,i,0,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (is_constant(fm.y) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.y.i -_NVARMAX], _constant[fm.x.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_y
#define val_fm_y(a,j,i,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_x
#define val_fm_x(a,j,i,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,i,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 142
{

#line 142 "/home/jiarongw/basilisk/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/jiarongw/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension(n.y,i,0,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }OMP(omp critical) if (_cfl > cfl) cfl = _cfl;
mpi_all_reduce_double (cfl, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 197
 }
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);
#line 208 "/home/jiarongw/basilisk/src/vof.h"
  scalar * fluxl = list_concat (NULL, tfluxl);
  fluxl = list_append (fluxl, flux);
  for (int l = depth() - 1; l >= 0; l--)
     { foreach_halo (prolongation, l){

#line 211 "/home/jiarongw/basilisk/src/vof.h"
 {
#line 220 "/home/jiarongw/basilisk/src/vof.h"
      if ((!is_leaf (neighbor(0,-1,)) && neighbor(0,-1,).neighbors && neighbor(0,-1,).pid >= 0))
 if (fluxl) for (scalar fl = *fluxl, *_i85 = fluxl; ((scalar *)&fl)->i >= 0; fl = *++_i85)
   val(fl,0,0,0) = (fine(fl,0,0,0) + fine(fl,1,0,0))/2.;
      if ((!is_leaf (neighbor(0,1,)) && neighbor(0,1,).neighbors && neighbor(0,1,).pid >= 0))
 if (fluxl) for (scalar fl = *fluxl, *_i86 = fluxl; ((scalar *)&fl)->i >= 0; fl = *++_i86)
   val(fl,0,1,0) = (fine(fl,0,2,0) + fine(fl,1,2,0))/2.;
#line 236 "/home/jiarongw/basilisk/src/vof.h"
    } } end_foreach_halo(); }
  pfree (fluxl,__func__,__FILE__,__LINE__);





  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 265 "/home/jiarongw/basilisk/src/vof.h"
   { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,j,i,k) val(a,j,i,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,j,i,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,j,i,k)
#line 265
foreach(){

#line 265 "/home/jiarongw/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tfluxl; if (tracers) for (t = *tracers, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tflux = *++_i12)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,1,0))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,i,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 265
foreach(){

#line 265 "/home/jiarongw/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tfluxl; if (tracers) for (t = *tracers, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tflux = *++_i12)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,1,0))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); } }
  boundary (((scalar []){c,{-1}}));
  boundary (tracers);

  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);
 delete (((scalar []){flux,alpha,n.x,n.y,{-1}})); }






void vof_advection (scalar * interfaces, int i)
{
  if (interfaces) for (scalar c = *interfaces, *_i87 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i87) {
#line 294 "/home/jiarongw/basilisk/src/vof.h"
    scalar cc= new_scalar("cc");
     { foreach(){

#line 295 "/home/jiarongw/basilisk/src/vof.h"

      val(cc,0,0,0) = (val(c,0,0,0) > 0.5); } end_foreach(); }






    void (* sweep[2]) (scalar, scalar);
    int d = 0;
    {
#line 305

      sweep[d++] = sweep_x;
#line 305

      sweep[d++] = sweep_y;}
    boundary (((scalar []){c,{-1}}));
    for (d = 0; d < 2; d++)
      sweep[(i + d) % 2] (c, cc);
   delete (((scalar []){cc,{-1}})); }
}

static int vof_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof_0 (const int i, const double t, Event * _ev) { trace ("vof_0", "/home/jiarongw/basilisk/src/vof.h", 313); 
  vof_advection (interfaces, i); end_trace("vof_0", "/home/jiarongw/basilisk/src/vof.h", 314);  return 0; } 
#line 14 "/home/jiarongw/basilisk/src/two-phase.h"

scalar f= {8}, * interfaces = ((scalar []){{8},{-1}});
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;





vector alphav= {{9},{10}};
scalar rhov= {11};

static int defaults_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_1 (const int i, const double t, Event * _ev) { trace ("defaults_1", "/home/jiarongw/basilisk/src/two-phase.h", 25);  {
  alpha = alphav;
  rho = rhov;





  if (mu1 || mu2)
    mu = new_face_vector("mu");
 end_trace("defaults_1", "/home/jiarongw/basilisk/src/two-phase.h", 35); } return 0; } 
#line 59 "/home/jiarongw/basilisk/src/two-phase.h"
static int properties_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int properties_0 (const int i, const double t, Event * _ev) { trace ("properties_0", "/home/jiarongw/basilisk/src/two-phase.h", 59);  {
#line 84 "/home/jiarongw/basilisk/src/two-phase.h"
  _attribute[f.i].prolongation = refine_bilinear;
  boundary (((scalar []){f,{-1}}));


   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 88
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 88
{

#line 88 "/home/jiarongw/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 88
{

#line 88 "/home/jiarongw/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  end_foreach_face_generic()
#line 95
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 88
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 88
{

#line 88 "/home/jiarongw/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 88
{

#line 88 "/home/jiarongw/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  end_foreach_face_generic()
#line 95
 end_foreach_face(); } }
   { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 96
foreach(){

#line 96 "/home/jiarongw/basilisk/src/two-phase.h"

    val(rhov,0,0,0) = val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2); } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 96
foreach(){

#line 96 "/home/jiarongw/basilisk/src/two-phase.h"

    val(rhov,0,0,0) = val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2); } end_foreach(); } }


  _attribute[f.i].prolongation = fraction_refine;
  boundary (((scalar []){f,{-1}}));

 end_trace("properties_0", "/home/jiarongw/basilisk/src/two-phase.h", 103); } return 0; } 
#line 11 "wavewind.c"
#line 1 "navier-stokes/conserving.h"
#line 1 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"
#line 14 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"
static void momentum_refine (Point point, scalar u) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 14 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 14

  refine_bilinear (point, u);
  double rhou = 0.;
   { foreach_child()
    rhou += val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2)*val(u,0,0,0); end_foreach_child(); }
  double du = val(u,0,0,0) - rhou/((1 << 2)*val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2));
   { foreach_child()
    val(u,0,0,0) += du; end_foreach_child(); }
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 14

  refine_bilinear (point, u);
  double rhou = 0.;
   { foreach_child()
    rhou += val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2)*val(u,0,0,0); end_foreach_child(); }
  double du = val(u,0,0,0) - rhou/((1 << 2)*val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2));
   { foreach_child()
    val(u,0,0,0) += du; end_foreach_child(); }
 }}

static void momentum_restriction (Point point, scalar u)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 25 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 25

  double rhou = 0.;
   { foreach_child()
    rhou += val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2)*val(u,0,0,0); end_foreach_child(); }
  val(u,0,0,0) = rhou/((1 << 2)*val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2));
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 25

  double rhou = 0.;
   { foreach_child()
    rhou += val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2)*val(u,0,0,0); end_foreach_child(); }
  val(u,0,0,0) = rhou/((1 << 2)*val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2));
 }}






static int defaults_2_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_2 (const int i, const double t, Event * _ev) { trace ("defaults_2", "/home/jiarongw/basilisk/src/navier-stokes/conserving.h", 37);  {
  stokes = true;
#line 48 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"
  int i = 0;
  while (all[i].i != f.i) i++;
  while (all[i].i)
    all[i] = all[i-1], i--;
  all[i] = f;





  {
#line 58
 {
    _attribute[u.x.i].refine = _attribute[u.x.i].prolongation = momentum_refine;
    _attribute[u.x.i].restriction = momentum_restriction;
  }
#line 58
 {
    _attribute[u.y.i].refine = _attribute[u.y.i].prolongation = momentum_refine;
    _attribute[u.y.i].restriction = momentum_restriction;
  }}

 end_trace("defaults_2", "/home/jiarongw/basilisk/src/navier-stokes/conserving.h", 63); } return 0; } 








#line 71

static double boundary_q1_x (Point neighbor, Point point, scalar q1)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 73 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"

  return clamp(val(f,0,0,0),0.,1.)*rho1*val(u.x,0,0,0);
}
#line 71

static double boundary_q1_y (Point neighbor, Point point, scalar q1)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 73 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"

  return clamp(val(f,0,0,0),0.,1.)*rho1*val(u.y,0,0,0);
}


#line 77

static double boundary_q2_x (Point neighbor, Point point, scalar q2)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 79 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"

  return (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.x,0,0,0);
}
#line 77

static double boundary_q2_y (Point neighbor, Point point, scalar q2)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 79 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"

  return (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.y,0,0,0);
}







#line 88

static void prolongation_q1_x (Point point, scalar q1) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 89 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"

   { foreach_child()
    val(q1,0,0,0) = clamp(val(f,0,0,0),0.,1.)*rho1*val(u.x,0,0,0); end_foreach_child(); }
}
#line 88

static void prolongation_q1_y (Point point, scalar q1) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 89 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"

   { foreach_child()
    val(q1,0,0,0) = clamp(val(f,0,0,0),0.,1.)*rho1*val(u.y,0,0,0); end_foreach_child(); }
}


#line 94

static void prolongation_q2_x (Point point, scalar q2) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 95 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"

   { foreach_child()
    val(q2,0,0,0) = (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.x,0,0,0); end_foreach_child(); }
}
#line 94

static void prolongation_q2_y (Point point, scalar q2) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 95 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"

   { foreach_child()
    val(q2,0,0,0) = (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.y,0,0,0); end_foreach_child(); }
}






static scalar * interfaces1 = NULL;

static int vof_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof_1 (const int i, const double t, Event * _ev) { trace ("vof_1", "/home/jiarongw/basilisk/src/navier-stokes/conserving.h", 107);  {






  vector q1= new_vector("q1"), q2= new_vector("q2");
  for (scalar s = *((scalar []){q1.x,q1.y,q2.x,q2.y,{-1}}), *_i88 = ((scalar []){q1.x,q1.y,q2.x,q2.y,{-1}}); ((scalar *)&s)->i >= 0; s = *++_i88)
    {
#line 116

      _attribute[s.i].v.x.i = -1;
#line 116

      _attribute[s.i].v.y.i = -1;}
  for (int i = 0; i < nboundary; i++)
    {
#line 119
 {
      _attribute[q1.x.i].boundary[i] = boundary_q1_x;
      _attribute[q2.x.i].boundary[i] = boundary_q2_x;
    }
#line 119
 {
      _attribute[q1.y.i].boundary[i] = boundary_q1_y;
      _attribute[q2.y.i].boundary[i] = boundary_q2_y;
    }}

  {
#line 124
 {
    _attribute[q1.x.i].prolongation = prolongation_q1_x;
    _attribute[q2.x.i].prolongation = prolongation_q2_x;
  }
#line 124
 {
    _attribute[q1.y.i].prolongation = prolongation_q1_y;
    _attribute[q2.y.i].prolongation = prolongation_q2_y;
  }}






   { foreach(){

#line 134 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"

    {
#line 135
 {
      double fc = clamp(val(f,0,0,0),0,1);
      val(q1.x,0,0,0) = fc*rho1*val(u.x,0,0,0);
      val(q2.x,0,0,0) = (1. - fc)*rho2*val(u.x,0,0,0);
    }
#line 135
 {
      double fc = clamp(val(f,0,0,0),0,1);
      val(q1.y,0,0,0) = fc*rho1*val(u.y,0,0,0);
      val(q2.y,0,0,0) = (1. - fc)*rho2*val(u.y,0,0,0);
    }} } end_foreach(); }
  boundary ((scalar *)((vector []){{q1.x,q1.y},{q2.x,q2.y},{{-1},{-1}}}));






  {
#line 147
 {
    _attribute[q2.x.i].inverse = true;
    _attribute[q1.x.i].gradient = _attribute[q2.x.i].gradient = _attribute[u.x.i].gradient;
  }
#line 147
 {
    _attribute[q2.y.i].inverse = true;
    _attribute[q1.y.i].gradient = _attribute[q2.y.i].gradient = _attribute[u.y.i].gradient;
  }}





  _attribute[f.i].tracers = (scalar *)((vector []){{q1.x,q1.y},{q2.x,q2.y},{{-1},{-1}}});
  vof_advection (((scalar []){f,{-1}}), i);





   { foreach(){

#line 163 "/home/jiarongw/basilisk/src/navier-stokes/conserving.h"

    {
#line 164

      val(u.x,0,0,0) = (val(q1.x,0,0,0) + val(q2.x,0,0,0))/(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);
#line 164

      val(u.y,0,0,0) = (val(q1.y,0,0,0) + val(q2.y,0,0,0))/(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}));





  interfaces1 = interfaces, interfaces = NULL;
 delete (((scalar []){q2.x,q2.y,q1.x,q1.y,{-1}}));  end_trace("vof_1", "/home/jiarongw/basilisk/src/navier-stokes/conserving.h", 173); } return 0; } 




static int tracer_advection_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_advection_0 (const int i, const double t, Event * _ev) { trace ("tracer_advection_0", "/home/jiarongw/basilisk/src/navier-stokes/conserving.h", 178);  {
  interfaces = interfaces1;
 end_trace("tracer_advection_0", "/home/jiarongw/basilisk/src/navier-stokes/conserving.h", 180); } return 0; } 
#line 12 "wavewind.c"
#line 1 "tension.h"
#line 1 "/home/jiarongw/basilisk/src/tension.h"
#line 15 "/home/jiarongw/basilisk/src/tension.h"
#line 1 "iforce.h"
#line 1 "/home/jiarongw/basilisk/src/iforce.h"
#line 20 "/home/jiarongw/basilisk/src/iforce.h"










static int defaults_3_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_3 (const int i, const double t, Event * _ev) { trace ("defaults_3", "/home/jiarongw/basilisk/src/iforce.h", 30);  {
  if (is_constant(a.x)) {
    a = new_face_vector("a");
     { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 33
{

#line 33 "/home/jiarongw/basilisk/src/iforce.h"

      val(a.x,0,0,0) = 0.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 33
{

#line 33 "/home/jiarongw/basilisk/src/iforce.h"

      val(a.y,0,0,0) = 0.; }  }}  end_foreach_face_generic()
#line 34
 end_foreach_face(); }
    boundary ((scalar *)((vector []){{a.x,a.y},{{-1},{-1}}}));
  }
 end_trace("defaults_3", "/home/jiarongw/basilisk/src/iforce.h", 37); } return 0; } 






static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_0 (const int i, const double t, Event * _ev) { trace ("acceleration_0", "/home/jiarongw/basilisk/src/iforce.h", 44); 
{





  scalar * list = NULL;
  if (interfaces) for (scalar f = *interfaces, *_i89 = interfaces; ((scalar *)&f)->i >= 0; f = *++_i89)
    if (_attribute[f.i].phi.i) {
      list = list_add (list, f);






       { foreach(){

#line 61 "/home/jiarongw/basilisk/src/iforce.h"

 val(f,0,0,0) = clamp (val(f,0,0,0), 0., 1.); } end_foreach(); }
      boundary (((scalar []){f,{-1}}));
    }
#line 74 "/home/jiarongw/basilisk/src/iforce.h"
  if (list) for (scalar f = *list, *_i90 = list; ((scalar *)&f)->i >= 0; f = *++_i90)
    _attribute[f.i].prolongation = _attribute[p.i].prolongation;
  boundary (list);
#line 87 "/home/jiarongw/basilisk/src/iforce.h"
  vector ia = a;
   { 
if (!is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 88
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 88
{

#line 88 "/home/jiarongw/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i91 = list; ((scalar *)&f)->i >= 0; f = *++_i91)
      if (val(f,0,0,0) != val(f,-1,0,0)) {
#line 100 "/home/jiarongw/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,-1,0,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,-1,0,0) < nodata ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 88
{

#line 88 "/home/jiarongw/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i91 = list; ((scalar *)&f)->i >= 0; f = *++_i91)
      if (val(f,0,0,0) != val(f,0,-1,0)) {
#line 100 "/home/jiarongw/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,-1,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,-1,0) < nodata ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      } }  }}  end_foreach_face_generic()
#line 109
 end_foreach_face(); }
if (is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 88
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 88
{

#line 88 "/home/jiarongw/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i91 = list; ((scalar *)&f)->i >= 0; f = *++_i91)
      if (val(f,0,0,0) != val(f,-1,0,0)) {
#line 100 "/home/jiarongw/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,-1,0,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,-1,0,0) < nodata ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 88
{

#line 88 "/home/jiarongw/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i91 = list; ((scalar *)&f)->i >= 0; f = *++_i91)
      if (val(f,0,0,0) != val(f,0,-1,0)) {
#line 100 "/home/jiarongw/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,-1,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,-1,0) < nodata ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      } }  }}  end_foreach_face_generic()
#line 109
 end_foreach_face(); }
if (!is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 88
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 88
{

#line 88 "/home/jiarongw/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i91 = list; ((scalar *)&f)->i >= 0; f = *++_i91)
      if (val(f,0,0,0) != val(f,-1,0,0)) {
#line 100 "/home/jiarongw/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,-1,0,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,-1,0,0) < nodata ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 88
{

#line 88 "/home/jiarongw/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i91 = list; ((scalar *)&f)->i >= 0; f = *++_i91)
      if (val(f,0,0,0) != val(f,0,-1,0)) {
#line 100 "/home/jiarongw/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,-1,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,-1,0) < nodata ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      } }  }}  end_foreach_face_generic()
#line 109
 end_foreach_face(); }
if (is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 88
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 88
{

#line 88 "/home/jiarongw/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i91 = list; ((scalar *)&f)->i >= 0; f = *++_i91)
      if (val(f,0,0,0) != val(f,-1,0,0)) {
#line 100 "/home/jiarongw/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,-1,0,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,-1,0,0) < nodata ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 88
{

#line 88 "/home/jiarongw/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i91 = list; ((scalar *)&f)->i >= 0; f = *++_i91)
      if (val(f,0,0,0) != val(f,0,-1,0)) {
#line 100 "/home/jiarongw/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,-1,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,-1,0) < nodata ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      } }  }}  end_foreach_face_generic()
#line 109
 end_foreach_face(); } }






  if (list) for (scalar f = *list, *_i92 = list; ((scalar *)&f)->i >= 0; f = *++_i92)
    _attribute[f.i].prolongation = fraction_refine;
  boundary (list);






  if (list) for (scalar f = *list, *_i93 = list; ((scalar *)&f)->i >= 0; f = *++_i93) {
    scalar phi = _attribute[f.i].phi;
    delete (((scalar []){phi,{-1}}));
    _attribute[f.i].phi.i = 0;
  }
  pfree (list,__func__,__FILE__,__LINE__);
 end_trace("acceleration_0", "/home/jiarongw/basilisk/src/iforce.h", 131); } return 0; } 
#line 16 "/home/jiarongw/basilisk/src/tension.h"
#line 1 "curvature.h"
#line 1 "/home/jiarongw/basilisk/src/curvature.h"
#line 12 "/home/jiarongw/basilisk/src/curvature.h"
static void curvature_restriction (Point point, scalar kappa)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 13 "/home/jiarongw/basilisk/src/curvature.h"

  double k = 0., s = 0.;
   { foreach_child()
    if (val(kappa,0,0,0) != nodata)
      k += val(kappa,0,0,0), s++; end_foreach_child(); }
  val(kappa,0,0,0) = s ? k/s : nodata;
}







static void curvature_prolongation (Point point, scalar kappa)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 28 "/home/jiarongw/basilisk/src/curvature.h"

   { foreach_child() {
    double sk = 0., s = 0.;
    for (int i = 0; i <= 1; i++)

      for (int j = 0; j <= 1; j++)




   if (coarse(kappa,child.x*i,child.y*j,child.z*k) != nodata)
     sk += coarse(kappa,child.x*i,child.y*j,child.z*k), s++;
    val(kappa,0,0,0) = s ? sk/s : nodata;
  } end_foreach_child(); }
}
#line 62 "/home/jiarongw/basilisk/src/curvature.h"
#line 1 "heights.h"
#line 1 "/home/jiarongw/basilisk/src/heights.h"
#line 29 "/home/jiarongw/basilisk/src/heights.h"
static inline double height (double H) {
  return H > 20./2. ? H - 20. : H < -20./2. ? H + 20. : H;
}

static inline int orientation (double H) {
  return fabs(H) > 20./2.;
}
#line 49 "/home/jiarongw/basilisk/src/heights.h"
static void half_column (Point point, scalar c, vector h, vector cs, int j)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 50 "/home/jiarongw/basilisk/src/heights.h"







  const int complete = -1;

  {
#line 59
 {







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.x,0,0,0) == 300.)
 state.s = complete, state.h = nodata;




      else {
 int s = (val(h.x,0,0,0) + 20./2.)/100.;
 state.h = val(h.x,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/jiarongw/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,i*j,0,0) : val(cs.x,(i - 2)*j,0,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/jiarongw/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/jiarongw/basilisk/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.x,0,0,0) = 300.;
      else if (S == complete)
 val(h.x,0,0,0) = H;
      else





 val(h.x,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/jiarongw/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.x,0,0,0) = nodata;
      else
 val(h.x,0,0,0) = (state.h > 1e10 ? nodata : state.h);
    }
  }
#line 59
 {







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.y,0,0,0) == 300.)
 state.s = complete, state.h = nodata;




      else {
 int s = (val(h.y,0,0,0) + 20./2.)/100.;
 state.h = val(h.y,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/jiarongw/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,0,i*j,0) : val(cs.y,0,(i - 2)*j,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/jiarongw/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/jiarongw/basilisk/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.y,0,0,0) = 300.;
      else if (S == complete)
 val(h.y,0,0,0) = H;
      else





 val(h.y,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/jiarongw/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.y,0,0,0) = nodata;
      else
 val(h.y,0,0,0) = (state.h > 1e10 ? nodata : state.h);
    }
  }}
}
#line 222 "/home/jiarongw/basilisk/src/heights.h"
static void column_propagation (vector h)
{
   { foreach(){

#line 224 "/home/jiarongw/basilisk/src/heights.h"

    for (int i = -2; i <= 2; i++)
      {
#line 226

 if (fabs(height(val(h.x,i,0,0))) <= 3.5 &&
     fabs(height(val(h.x,i,0,0)) + i) < fabs(height(val(h.x,0,0,0))))
   val(h.x,0,0,0) = val(h.x,i,0,0) + i;
#line 226

 if (fabs(height(val(h.y,0,i,0))) <= 3.5 &&
     fabs(height(val(h.y,0,i,0)) + i) < fabs(height(val(h.y,0,0,0))))
   val(h.y,0,0,0) = val(h.y,0,i,0) + i;}; } end_foreach(); }
  boundary ((scalar *)((vector []){{h.x,h.y},{{-1},{-1}}}));
}
#line 291 "/home/jiarongw/basilisk/src/heights.h"

#line 291

static void refine_h_x (Point point, scalar h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 293 "/home/jiarongw/basilisk/src/heights.h"





  bool complete = true;
   { foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(i,0,) &&
   !(!is_leaf(neighbor(i,0,)) && !neighbor(i,0,).neighbors && neighbor(i,0,).pid >= 0) && !(neighbor(i,0,).pid < 0) &&
   fabs(height(val(h,i,0,0))) <= 3.5 &&
   fabs(height(val(h,i,0,0)) + i) < fabs(height(val(h,0,0,0))))
 val(h,0,0,0) = val(h,i,0,0) + i;
    if (val(h,0,0,0) == nodata)
      complete = false;
  } end_foreach_child(); }
  if (complete)
    return;
#line 319 "/home/jiarongw/basilisk/src/heights.h"
  int ori = orientation(val(h,0,0,0));

  for (int i = -1; i <= 1; i++)
    if (val(h,0,i,0) == nodata || orientation(val(h,0,i,0)) != ori)
      return;

  double h0 = (30.*height(val(h,0,0,0)) + height(val(h,0,1,0)) + height(val(h,0,-1,0)))/16.
    + 20.*ori;
  double dh = (height(val(h,0,1,0)) - height(val(h,0,-1,0)))/4.;
   { foreach_child()
    if (val(h,0,0,0) == nodata)
      val(h,0,0,0) = h0 + dh*child.y - child.x/2.; end_foreach_child(); }
#line 353 "/home/jiarongw/basilisk/src/heights.h"
}
#line 291

static void refine_h_y (Point point, scalar h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 293 "/home/jiarongw/basilisk/src/heights.h"





  bool complete = true;
   { foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(0,i,) &&
   !(!is_leaf(neighbor(0,i,)) && !neighbor(0,i,).neighbors && neighbor(0,i,).pid >= 0) && !(neighbor(0,i,).pid < 0) &&
   fabs(height(val(h,0,i,0))) <= 3.5 &&
   fabs(height(val(h,0,i,0)) + i) < fabs(height(val(h,0,0,0))))
 val(h,0,0,0) = val(h,0,i,0) + i;
    if (val(h,0,0,0) == nodata)
      complete = false;
  } end_foreach_child(); }
  if (complete)
    return;
#line 319 "/home/jiarongw/basilisk/src/heights.h"
  int ori = orientation(val(h,0,0,0));

  for (int i = -1; i <= 1; i++)
    if (val(h,i,0,0) == nodata || orientation(val(h,i,0,0)) != ori)
      return;

  double h0 = (30.*height(val(h,0,0,0)) + height(val(h,1,0,0)) + height(val(h,-1,0,0)))/16.
    + 20.*ori;
  double dh = (height(val(h,1,0,0)) - height(val(h,-1,0,0)))/4.;
   { foreach_child()
    if (val(h,0,0,0) == nodata)
      val(h,0,0,0) = h0 + dh*child.x - child.y/2.; end_foreach_child(); }
#line 353 "/home/jiarongw/basilisk/src/heights.h"
}







void heights (scalar c, vector h)
{ trace ("heights", "/home/jiarongw/basilisk/src/heights.h", 362);
  vector cs= new_vector("cs");
  {
#line 364

    for (int i = 0; i < nboundary; i++)
      _attribute[cs.x.i].boundary[i] = _attribute[c.i].boundary[i];
#line 364

    for (int i = 0; i < nboundary; i++)
      _attribute[cs.y.i].boundary[i] = _attribute[c.i].boundary[i];}





  restriction (((scalar []){c,{-1}}));
  for (int j = -1; j <= 1; j += 2) {





     { foreach_level(0){

#line 379 "/home/jiarongw/basilisk/src/heights.h"

      {
#line 380

        val(h.x,0,0,0) = nodata;
#line 380

        val(h.y,0,0,0) = nodata;}; } end_foreach_level(); }

    for (int l = 1; l <= depth(); l++) {




       { foreach_level (l){

#line 388 "/home/jiarongw/basilisk/src/heights.h"

 {
#line 389

   val(cs.x,0,0,0) = val(c,2*j,0,0);
#line 389

   val(cs.y,0,0,0) = val(c,0,2*j,0);}; } end_foreach_level(); }
#line 400 "/home/jiarongw/basilisk/src/heights.h"
       { foreach_level (l - 1){

#line 400 "/home/jiarongw/basilisk/src/heights.h"

 {
#line 401
 {
   val(cs.x,0,0,0) = val(c,j,0,0);
   val(cs.x,j,0,0) = val(c,2*j,0,0);
        }
#line 401
 {
   val(cs.y,0,0,0) = val(c,0,j,0);
   val(cs.y,0,j,0) = val(c,0,2*j,0);
        }} } end_foreach_level(); }






       { foreach_halo (prolongation, l - 1){

#line 411 "/home/jiarongw/basilisk/src/heights.h"

 {
#line 412

   _attribute[c.i].prolongation (point, cs.x);
#line 412

   _attribute[c.i].prolongation (point, cs.y);}; } end_foreach_halo(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, (scalar *)((vector []){{cs.x,cs.y},{{-1},{-1}}}), l); };





       { foreach_level (l){

#line 420 "/home/jiarongw/basilisk/src/heights.h"

        half_column (point, c, h, cs, j); } end_foreach_level(); }
    }
  }






  {
#line 430
 {
    _attribute[h.x.i].prolongation = no_data;
    _attribute[h.x.i].restriction = no_restriction;
  }
#line 430
 {
    _attribute[h.y.i].prolongation = no_data;
    _attribute[h.y.i].restriction = no_restriction;
  }}
  boundary ((scalar *)((vector []){{h.x,h.y},{{-1},{-1}}}));






  {
#line 441

    _attribute[h.x.i].prolongation = refine_h_x;
#line 441

    _attribute[h.y.i].prolongation = refine_h_y;}




  column_propagation (h);
 delete (((scalar []){cs.x,cs.y,{-1}}));  end_trace("heights", "/home/jiarongw/basilisk/src/heights.h", 448); }
#line 63 "/home/jiarongw/basilisk/src/curvature.h"



#line 65

static double kappa_y (Point point, vector h) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 66 "/home/jiarongw/basilisk/src/curvature.h"

  int ori = orientation(val(h.y,0,0,0));
  for (int i = -1; i <= 1; i++)
    if (val(h.y,i,0,0) == nodata || orientation(val(h.y,i,0,0)) != ori)
      return nodata;
  double hx = (val(h.y,1,0,0) - val(h.y,-1,0,0))/2.;
  double hxx = (val(h.y,1,0,0) + val(h.y,-1,0,0) - 2.*val(h.y,0,0,0))/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);
}
#line 65

static double kappa_x (Point point, vector h) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 66 "/home/jiarongw/basilisk/src/curvature.h"

  int ori = orientation(val(h.x,0,0,0));
  for (int i = -1; i <= 1; i++)
    if (val(h.x,0,i,0) == nodata || orientation(val(h.x,0,i,0)) != ori)
      return nodata;
  double hx = (val(h.x,0,1,0) - val(h.x,0,-1,0))/2.;
  double hxx = (val(h.x,0,1,0) + val(h.x,0,-1,0) - 2.*val(h.x,0,0,0))/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);
}
#line 99 "/home/jiarongw/basilisk/src/curvature.h"
static double height_curvature (Point point, scalar c, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 100 "/home/jiarongw/basilisk/src/curvature.h"







  typedef struct {
    double n;
    double (* kappa) (Point, vector);
  } NormKappa;
  struct { NormKappa x, y, z; } n;
  {
#line 112

    n.x.n = val(c,1,0,0) - val(c,-1,0,0), n.x.kappa = kappa_x;
#line 112

    n.y.n = val(c,0,1,0) - val(c,0,-1,0), n.y.kappa = kappa_y;}
  double (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);




  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormKappa, n.x, n.y);
#line 131 "/home/jiarongw/basilisk/src/curvature.h"
  double kappa = nodata;
  {
#line 132

    if (kappa == nodata) {
      kappa = n.x.kappa (point, h);
      if (kappa != nodata) {
 kappaf = n.x.kappa;
 if (n.x.n < 0.)
   kappa = - kappa;
      }
    }
#line 132

    if (kappa == nodata) {
      kappa = n.y.kappa (point, h);
      if (kappa != nodata) {
 kappaf = n.y.kappa;
 if (n.y.n < 0.)
   kappa = - kappa;
      }
    }}

  if (kappa != nodata) {




    if (fabs(kappa) > 1./Delta)
      kappa = sign(kappa)/Delta;
#line 167 "/home/jiarongw/basilisk/src/curvature.h"
  }

  return kappa;
}
#line 184 "/home/jiarongw/basilisk/src/curvature.h"
#line 1 "parabola.h"
#line 1 "/home/jiarongw/basilisk/src/parabola.h"
#line 1 "utils.h"
#line 2 "/home/jiarongw/basilisk/src/parabola.h"






typedef struct {
  coord o;

  coord m;
  double ** M, rhs[3], a[3];
#line 21 "/home/jiarongw/basilisk/src/parabola.h"
} ParabolaFit;

static void parabola_fit_init (ParabolaFit * p, coord o, coord m)
{
  {
#line 25

    p->o.x = o.x;
#line 25

    p->o.y = o.y;}

  {
#line 28

    p->m.x = m.x;
#line 28

    p->m.y = m.y;}
  normalize (&p->m);
  int n = 3;
#line 65 "/home/jiarongw/basilisk/src/parabola.h"
  p->M = (double **) matrix_new (n, n, sizeof(double));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      p->M[i][j] = 0.;
    p->rhs[i] = 0.;
  }
}

static void parabola_fit_add (ParabolaFit * p, coord m, double w)
{

  double x1 = m.x - p->o.x, y1 = m.y - p->o.y;
  double x = p->m.y*x1 - p->m.x*y1;
  double y = p->m.x*x1 + p->m.y*y1;
  double x2 = w*x*x, x3 = x2*x, x4 = x3*x;
  p->M[0][0] += x4;
  p->M[1][0] += x3; p->M[1][1] += x2;
  p->M[2][1] += w*x; p->M[2][2] += w;
  p->rhs[0] += x2*y; p->rhs[1] += w*x*y; p->rhs[2] += w*y;
#line 111 "/home/jiarongw/basilisk/src/parabola.h"
}

static double parabola_fit_solve (ParabolaFit * p)
{

  p->M[0][1] = p->M[1][0];
  p->M[0][2] = p->M[2][0] = p->M[1][1];
  p->M[1][2] = p->M[2][1];
  double pivmin = matrix_inverse (p->M, 3, 1e-10);
  if (pivmin) {
    p->a[0] = p->M[0][0]*p->rhs[0] + p->M[0][1]*p->rhs[1] + p->M[0][2]*p->rhs[2];
    p->a[1] = p->M[1][0]*p->rhs[0] + p->M[1][1]*p->rhs[1] + p->M[1][2]*p->rhs[2];
  }
  else
    p->a[0] = p->a[1] = 0.;
#line 158 "/home/jiarongw/basilisk/src/parabola.h"
  matrix_free (p->M);
  return pivmin;
}

static double parabola_fit_curvature (ParabolaFit * p,
          double kappamax, double * kmax)
{
  double kappa;

  double dnm = 1. + sq(p->a[1]);
  kappa = - 2.*p->a[0]/pow(dnm, 3/2.);
  if (kmax)
    *kmax = fabs (kappa);
#line 190 "/home/jiarongw/basilisk/src/parabola.h"
  if (fabs (kappa) > kappamax) {
    if (kmax)
      *kmax = kappamax;
    return kappa > 0. ? kappamax : - kappamax;
  }
  return kappa;
}
#line 185 "/home/jiarongw/basilisk/src/curvature.h"






static int independents (coord * p, int n)
{
  if (n < 2)
    return n;
  int ni = 1;
  for (int j = 1; j < n; j++) {
    bool depends = false;
    for (int i = 0; i < j && !depends; i++) {
      double d2 = 0.;
      {
#line 200

 d2 += sq(p[i].x - p[j].x);
#line 200

 d2 += sq(p[i].y - p[j].y);}
      depends = (d2 < sq(0.5));
    }
    ni += !depends;
  }
  return ni;
}






static double height_curvature_fit (Point point, scalar c, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 215 "/home/jiarongw/basilisk/src/curvature.h"






  coord ip[2 == 2 ? 6 : 27];
  int n = 0;




  {
#line 227
 {





    int n1 = 0, n2 = 0;

    for (int i = -1; i <= 1; i++)
      if (val(h.y,i,0,0) != nodata) {
 if (orientation(val(h.y,i,0,0))) n1++; else n2++;
      }







    int ori = (n1 > n2);







    for (int i = -1; i <= 1; i++)
      if (val(h.y,i,0,0) != nodata && orientation(val(h.y,i,0,0)) == ori)
 ip[n].x = i, ip[n++].y = height(val(h.y,i,0,0));






  }
#line 227
 {





    int n1 = 0, n2 = 0;

    for (int i = -1; i <= 1; i++)
      if (val(h.x,0,i,0) != nodata) {
 if (orientation(val(h.x,0,i,0))) n1++; else n2++;
      }







    int ori = (n1 > n2);







    for (int i = -1; i <= 1; i++)
      if (val(h.x,0,i,0) != nodata && orientation(val(h.x,0,i,0)) == ori)
 ip[n].y = i, ip[n++].x = height(val(h.x,0,i,0));






  }}





  if (independents (ip, n) < (2 == 2 ? 3 : 9))
    return nodata;





  coord m = mycs (point, c), fc;
  double alpha = line_alpha (val(c,0,0,0), m);
  double area = line_length_center(m,alpha,&fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);

  NOT_UNUSED(area);
  parabola_fit_add (&fit, fc, .1);
#line 292 "/home/jiarongw/basilisk/src/curvature.h"
  for (int i = 0; i < n; i++)
    parabola_fit_add (&fit, ip[i], 1.);
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;
}






static double centroids_curvature_fit (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 308 "/home/jiarongw/basilisk/src/curvature.h"






  coord m = mycs (point, c), fc;
  double alpha = line_alpha (val(c,0,0,0), m);
  line_length_center(m,alpha,&fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);





  coord r = {x,y,z};
   { foreach_neighbor(1)
    if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
      coord m = mycs (point, c), fc;
      double alpha = line_alpha (val(c,0,0,0), m);
      double area = line_length_center(m,alpha,&fc);
      coord rn = {x,y,z};
      {
#line 331

 fc.x += (rn.x - r.x)/Delta;
#line 331

 fc.y += (rn.y - r.y)/Delta;}
      parabola_fit_add (&fit, fc, area);
    } end_foreach_neighbor(); }
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;
}
#line 354 "/home/jiarongw/basilisk/src/curvature.h"
static inline bool interfacial (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 355 "/home/jiarongw/basilisk/src/curvature.h"

  if (val(c,0,0,0) >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      {
#line 358

 if (val(c,i,0,0) <= 0.)
   return true;
#line 358

 if (val(c,0,i,0) <= 0.)
   return true;}
  }
  else if (val(c,0,0,0) <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      {
#line 364

 if (val(c,i,0,0) >= 1.)
   return true;
#line 364

 if (val(c,0,i,0) >= 1.)
   return true;}
  }
  else
    return true;
  return false;
}
#line 384 "/home/jiarongw/basilisk/src/curvature.h"
typedef struct {
  int h;
  int f;
  int a;
  int c;
} cstats;

struct Curvature {
  scalar c, kappa;
  double sigma;
  bool add;
};


cstats curvature (struct Curvature p)
{ trace ("curvature", "/home/jiarongw/basilisk/src/curvature.h", 399);
  scalar c = p.c, kappa = p.kappa;
  double sigma = p.sigma ? p.sigma : 1.;
  cstats s = {0,0,0,0};
  vector h= new_vector("h");
  heights (c, h);






  _attribute[kappa.i].refine = _attribute[kappa.i].prolongation = curvature_prolongation;
  _attribute[kappa.i].restriction = curvature_restriction;






  scalar k= new_scalar("k");
  scalar_clone (k, kappa);

   { foreach(){

#line 422 "/home/jiarongw/basilisk/src/curvature.h"
 {




    if (!interfacial (point, c))
      val(k,0,0,0) = nodata;





    else if ((val(k,0,0,0) = height_curvature (point, c, h)) != nodata)
      s.h++;
    else if ((val(k,0,0,0) = height_curvature_fit (point, c, h)) != nodata)
      s.f++;
  } } end_foreach(); }
  boundary (((scalar []){k,{-1}}));

   { foreach(){

#line 441 "/home/jiarongw/basilisk/src/curvature.h"
 {





    double kf;
    if (val(k,0,0,0) < nodata)
      kf = val(k,0,0,0);
    else if (interfacial (point, c)) {





      double sk = 0., a = 0.;
       { foreach_neighbor(1)
 if (val(k,0,0,0) < nodata)
   sk += val(k,0,0,0), a++; end_foreach_neighbor(); }
      if (a > 0.)
 kf = sigma*sk/a, s.a++;
      else




 kf = sigma*centroids_curvature_fit (point, c), s.c++;
    }
    else
      kf = nodata;




    if (kf == nodata)
      val(kappa,0,0,0) = nodata;
    else if (p.add)
      val(kappa,0,0,0) += sigma*kf;
    else
      val(kappa,0,0,0) = sigma*kf;
  } } end_foreach(); }
  boundary (((scalar []){kappa,{-1}}));

  { cstats _ret =  s; delete (((scalar []){k,h.x,h.y,{-1}}));  end_trace("curvature", "/home/jiarongw/basilisk/src/curvature.h", 484);  return _ret; }
 delete (((scalar []){k,h.x,h.y,{-1}}));  end_trace("curvature", "/home/jiarongw/basilisk/src/curvature.h", 485); }
#line 504 "/home/jiarongw/basilisk/src/curvature.h"

#line 504

static double pos_x (Point point, vector h, coord * G, coord * Z)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 506 "/home/jiarongw/basilisk/src/curvature.h"

  if (fabs(height(val(h.x,0,0,0))) > 1.)
    return nodata;
  coord o = {x, y, z};
  o.x += height(val(h.x,0,0,0))*Delta;
  double pos = 0.;
  {
#line 512

    pos += (o.x - Z->x)*G->x;
#line 512

    pos += (o.y - Z->y)*G->y;}
  return pos;
}
#line 504

static double pos_y (Point point, vector h, coord * G, coord * Z)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 506 "/home/jiarongw/basilisk/src/curvature.h"

  if (fabs(height(val(h.y,0,0,0))) > 1.)
    return nodata;
  coord o = {x, y, z};
  o.y += height(val(h.y,0,0,0))*Delta;
  double pos = 0.;
  {
#line 512

    pos += (o.y - Z->y)*G->y;
#line 512

    pos += (o.x - Z->x)*G->x;}
  return pos;
}







static double height_position (Point point, scalar f, vector h,
          coord * G, coord * Z)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 525 "/home/jiarongw/basilisk/src/curvature.h"







  typedef struct {
    double n;
    double (* pos) (Point, vector, coord *, coord *);
  } NormPos;
  struct { NormPos x, y, z; } n;
  {
#line 537

    n.x.n = val(f,1,0,0) - val(f,-1,0,0), n.x.pos = pos_x;
#line 537

    n.y.n = val(f,0,1,0) - val(f,0,-1,0), n.y.pos = pos_y;}




  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormPos, n.x, n.y);
#line 555 "/home/jiarongw/basilisk/src/curvature.h"
  double pos = nodata;
  {
#line 556

    if (pos == nodata)
      pos = n.x.pos (point, h, G, Z);
#line 556

    if (pos == nodata)
      pos = n.y.pos (point, h, G, Z);}

  return pos;
}
#line 572 "/home/jiarongw/basilisk/src/curvature.h"
struct Position {
  scalar f, pos;
  coord G, Z;
  bool add;
};

void position (struct Position p)
{
  scalar f = p.f, pos = p.pos;
  coord * G = &p.G, * Z = &p.Z;






  _attribute[pos.i].refine = _attribute[pos.i].prolongation = curvature_prolongation;
  _attribute[pos.i].restriction = curvature_restriction;


  vector h= new_vector("h");
  heights (f, h);
   { foreach(){

#line 594 "/home/jiarongw/basilisk/src/curvature.h"
 {
    if (interfacial (point, f)) {
      double hp = height_position (point, f, h, G, Z);
      if (hp == nodata) {





 coord n = mycs (point, f), o = {x,y,z}, c;
 double alpha = line_alpha (val(f,0,0,0), n);
 line_length_center(n,alpha,&c);
 hp = 0.;
 {
#line 607

   hp += (o.x + Delta*c.x - Z->x)*G->x;
#line 607

   hp += (o.y + Delta*c.y - Z->y)*G->y;}
      }
      if (p.add)
 val(pos,0,0,0) += hp;
      else
 val(pos,0,0,0) = hp;
    }
    else
      val(pos,0,0,0) = nodata;
  } } end_foreach(); }
  boundary (((scalar []){pos,{-1}}));
 delete (((scalar []){h.x,h.y,{-1}})); }
#line 17 "/home/jiarongw/basilisk/src/tension.h"







#line 36 "/home/jiarongw/basilisk/src/tension.h"
static int stability_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability_1 (const int i, const double t, Event * _ev) { trace ("stability_1", "/home/jiarongw/basilisk/src/tension.h", 36);  {





  double amin = HUGE, amax = -HUGE, dmin = HUGE;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _amin = amin; double _amax = amax; double _dmin = dmin; 
#line 43

if (!is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 43
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 43
{

#line 43 "/home/jiarongw/basilisk/src/tension.h"
 {
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > _amax) _amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < _amin) _amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 43
{

#line 43 "/home/jiarongw/basilisk/src/tension.h"
 {
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > _amax) _amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < _amin) _amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  end_foreach_face_generic()
#line 47
 end_foreach_face(); }
if (is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 43
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 43
{

#line 43 "/home/jiarongw/basilisk/src/tension.h"
 {
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > _amax) _amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < _amin) _amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 43
{

#line 43 "/home/jiarongw/basilisk/src/tension.h"
 {
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > _amax) _amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < _amin) _amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  end_foreach_face_generic()
#line 47
 end_foreach_face(); }
if (!is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 43
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 43
{

#line 43 "/home/jiarongw/basilisk/src/tension.h"
 {
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > _amax) _amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < _amin) _amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 43
{

#line 43 "/home/jiarongw/basilisk/src/tension.h"
 {
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > _amax) _amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < _amin) _amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  end_foreach_face_generic()
#line 47
 end_foreach_face(); }
if (is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 43
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 43
{

#line 43 "/home/jiarongw/basilisk/src/tension.h"
 {
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > _amax) _amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < _amin) _amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 43
{

#line 43 "/home/jiarongw/basilisk/src/tension.h"
 {
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > _amax) _amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < _amin) _amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  end_foreach_face_generic()
#line 47
 end_foreach_face(); }OMP(omp critical) if (_amin < amin) amin = _amin;
mpi_all_reduce_double (amin, MPI_MIN);
OMP(omp critical) if (_amax > amax) amax = _amax;
mpi_all_reduce_double (amax, MPI_MAX);
OMP(omp critical) if (_dmin < dmin) dmin = _dmin;
mpi_all_reduce_double (dmin, MPI_MIN);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 47
 }
  double rhom = (1./amin + 1./amax)/2.;





  double sigma = 0.;
  if (interfaces) for (scalar c = *interfaces, *_i94 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i94)
    sigma += _attribute[c.i].sigma;
  if (sigma) {
    double dt = sqrt (rhom*cube(dmin)/(pi*sigma));
    if (dt < dtmax)
      dtmax = dt;
  }
 end_trace("stability_1", "/home/jiarongw/basilisk/src/tension.h", 62); } return 0; } 







static int acceleration_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_1 (const int i, const double t, Event * _ev) { trace ("acceleration_1", "/home/jiarongw/basilisk/src/tension.h", 70); 
{




  if (interfaces) for (scalar f = *interfaces, *_i95 = interfaces; ((scalar *)&f)->i >= 0; f = *++_i95)
    if (_attribute[f.i].sigma) {





      scalar phi = _attribute[f.i].phi;
      if (phi.i)
 curvature ((struct Curvature){f, phi, _attribute[f.i].sigma, .add = true});
      else {
 phi = new_scalar("phi");
 curvature ((struct Curvature){f, phi, _attribute[f.i].sigma, .add = false});
 _attribute[f.i].phi = phi;
      }
    }
 end_trace("acceleration_1", "/home/jiarongw/basilisk/src/tension.h", 92); } return 0; } 
#line 13 "wavewind.c"
#line 1 "reduced.h"
#line 1 "/home/jiarongw/basilisk/src/reduced.h"
#line 19 "/home/jiarongw/basilisk/src/reduced.h"
coord G = {0.,0.,0.}, Z = {0.,0.,0.};
#line 36 "/home/jiarongw/basilisk/src/reduced.h"
static int acceleration_2_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_2 (const int i, const double t, Event * _ev) { trace ("acceleration_2", "/home/jiarongw/basilisk/src/reduced.h", 36); 
{
  scalar phi = _attribute[f.i].phi;
  coord G1;
  {
#line 40

    G1.x = (rho2 - rho1)*G.x;
#line 40

    G1.y = (rho2 - rho1)*G.y;}

  if (phi.i)
    position ((struct Position){f, phi, G1, Z, .add = true});
  else {
    phi = new_scalar("phi");
    position ((struct Position){f, phi, G1, Z, .add = false});
    _attribute[f.i].phi = phi;
  }
 end_trace("acceleration_2", "/home/jiarongw/basilisk/src/reduced.h", 50); } return 0; } 
#line 14 "wavewind.c"
#line 1 "view.h"
#line 1 "/home/jiarongw/basilisk/src/view.h"
#line 64 "/home/jiarongw/basilisk/src/view.h"
# include <GL/gl.h>
# include <GL/glu.h>


#include <gl/framebuffer.h>
#include <gl/gl2ps/gl2ps.h>
#include <gl/trackball.h>
#include <gl/utils.h>



#line 1 "input.h"
#line 1 "/home/jiarongw/basilisk/src/input.h"
#line 16 "/home/jiarongw/basilisk/src/input.h"
struct InputPGM {

  scalar s;
  FILE * fp;

  double ox, oy, width;
};

void input_pgm (struct InputPGM p)
{
  scalar s = p.s;
  if (p.width == 0.) p.width = L0;

  char line[81];
  if (!fgets (line, 81, p.fp)) {
    fprintf (qstderr(), "input_pgm: could not read magic number\n");
    exit (1);
  }
  if (strcmp (line, "P2\n") && strcmp (line, "P5\n")) {
    fprintf (qstderr(), "input_pgm: magic number '%s' does not match PGM\n",
      line);
    exit (1);
  }
  int binary = !strcmp (line, "P5\n");
  if (!fgets (line, 81, p.fp)) {
    fprintf (qstderr(), "input_pgm: could not read width and height\n");
    exit (1);
  }
  int width, height;
  while (line[0] == '#' && fgets (line, 81, p.fp));
  if (line[0] == '#' || sscanf (line, "%d %d", &width, &height) != 2) {
    fprintf (qstderr(), "input_pgm: could not read width and height\n");
    exit (1);
  }
  if (!fgets (line, 81, p.fp)) {
    fprintf (qstderr(), "input_pgm: could not read maxval\n");
    exit (1);
  }
  int maxval;
  if (sscanf (line, "%d", &maxval) != 1) {
    fprintf (qstderr(), "input_pgm: could not read maxval\n");
    exit (1);
  }
  if (maxval < 256) {
    unsigned char * a = ((unsigned char *) pmalloc ((width*height)*sizeof(unsigned char),__func__,__FILE__,__LINE__));
    size_t n = 0;
    if (binary)
      n = fread (a, 1, width*height, p.fp);
    else {
      int v;
      while (n < width*height && fscanf (p.fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != width*height) {
      fprintf (qstderr(), "input_pgm: read only %ld values\n", n);
      exit (1);
    }
     { foreach(){

#line 73 "/home/jiarongw/basilisk/src/input.h"
 {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < height)
 val(s,0,0,0) = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    } } end_foreach(); }
    pfree (a,__func__,__FILE__,__LINE__);
  }
  else {
    unsigned short * a = ((unsigned short *) pmalloc ((width*height)*sizeof(unsigned short),__func__,__FILE__,__LINE__));
    size_t n = 0;
    if (binary)
      n = fread (a, 2, width*height, p.fp);
    else {
      int v;
      while (n < width*height && fscanf (p.fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != width*height) {
      fprintf (qstderr(), "input_pgm: read only %ld values\n", n);
      exit (1);
    }
     { foreach(){

#line 96 "/home/jiarongw/basilisk/src/input.h"
 {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < height)
 val(s,0,0,0) = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    } } end_foreach(); }
    pfree (a,__func__,__FILE__,__LINE__);
  }
}

static void next_char (FILE * fp, int target)
{
  int c = fgetc(fp), para = 0;
  while (c != EOF && (c != target || para > 0)) {
    if (c == '{') para++;
    if (c == '}') para--;
    c = fgetc(fp);
  }
  if (c != target) {
    fprintf (qstderr(), "input_gfs(): error: expecting '%c'\n", target);
    exit (1);
  }
}

static int next_string (FILE * fp, const char * target)
{
  int slen = strlen (target), para = 0;
  char s[slen + 1];
  s[slen] = '\0';
  int len = 0, c = fgetc (fp);
  while (c != EOF && len < slen) {
    if (c == '{') para++;
    if (c == '}') para--;
    s[len++] = c;
    c = fgetc (fp);
  }
  while (c != EOF && para >= 0) {
    if (!strcmp (s, target) && para == 0)
      break;
    if (c == '{') para++;
    if (c == '}') para--;
    for (int i = 0; i < slen - 1; i++)
      s[i] = s[i+1];
    s[slen - 1] = c;
    c = fgetc (fp);
  }
  if (strcmp (s, target))
    c = -1;
  return c;
}
#line 166 "/home/jiarongw/basilisk/src/input.h"

void input_gfs (struct OutputGfs p)
{ trace ("input_gfs", "/home/jiarongw/basilisk/src/input.h", 168);
  not_mpi_compatible();

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdin;
    else if (!(p.fp = fopen (p.file, "r"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }
  bool input_all = (p.list == all);
  if (p.list == NULL) p.list = all;


  init_grid (1);


  next_char (p.fp, '{');

  char * s = ((char *) pmalloc ((1)*sizeof(char),__func__,__FILE__,__LINE__));
  int len = 0;
  int c = fgetc(p.fp);
  while (c != EOF && c != '}') {
    s[len++] = c;
    s = (char *) prealloc (s, (len + 1)*sizeof(char),__func__,__FILE__,__LINE__);
    s[len] = '\0';
    c = fgetc(p.fp);
  }
  if (c != '}') {
    fprintf (qstderr(), "input_gfs(): error: expecting '}'\n");
    exit (1);
  }

  char * s1 = strstr (s, "variables");
  if (!s1) {
    fprintf (qstderr(), "input_gfs(): error: expecting 'variables'\n");
    exit (1);
  }

  s1 = strstr (s1, "=");
  if (!s1) {
    fprintf (qstderr(), "input_gfs(): error: expecting '='\n");
    exit (1);
  }
  s1++;

  while (strchr (" \t", *s1))
    s1++;

  scalar * input = NULL;
  s1 = strtok (s1, ", \t");
  while (s1) {
    char * name = replace (s1, '_', '.', false);
    bool found = false;
    if (p.list) for (scalar s = *p.list, *_i96 = p.list; ((scalar *)&s)->i >= 0; s = *++_i96)
      if (!is_constant(s) && _attribute[s.i].name && !strcmp (_attribute[s.i].name, name)) {
 input = list_append (input, s);
 found = true; break;
      }
    if (!found) {
      if (input_all) {
 scalar s = new_scalar("s");
 pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
 _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
 input = list_append (input, s);
      }
      else
 input = list_append (input, (scalar){INT_MAX});
    }
    pfree (name,__func__,__FILE__,__LINE__);
    s1 = strtok (NULL, ", \t");
  }
  pfree (s,__func__,__FILE__,__LINE__);

  next_char (p.fp, '{');
  double t1 = 0.;
  if (next_string (p.fp, "Time") >= 0) {
    next_char (p.fp, '{');
    next_char (p.fp, 't');
    next_char (p.fp, '=');
    if (fscanf (p.fp, "%lf", &t1) != 1) {
      fprintf (qstderr(), "input_gfs(): error: expecting 't'\n");
      exit (1);
    }
    next_char (p.fp, '}');
    next_char (p.fp, '}');
  }

  if (next_string (p.fp, "Box") < 0) {
    fprintf (qstderr(), "input_gfs(): error: expecting 'GfsBox'\n");
    exit (1);
  }

  next_char (p.fp, '{');
  next_char (p.fp, '{');
  next_char (p.fp, '\n');

  scalar * listm = ((scalar []){cm,fm.x,fm.y,{-1}});
  scalar * listr = !is_constant(cm) ? listm : NULL;
  NOT_UNUSED (listr);

   { foreach_cell(){

#line 273 "/home/jiarongw/basilisk/src/input.h"
 {
    unsigned flags;
    if (fread (&flags, sizeof (unsigned), 1, p.fp) != 1) {
      fprintf (qstderr(), "input_gfs(): error: expecting 'flags'\n");
      exit (1);
    }
    if (!(flags & (1 << 4)) && is_leaf(cell))
      refine_cell (point, listr, 0, NULL);
    double a;
    if (fread (&a, sizeof (double), 1, p.fp) != 1 || a != -1) {
      fprintf (qstderr(), "input_gfs(): error: expecting '-1'\n");
      exit (1);
    }
    if (input) for (scalar s = *input, *_i97 = input; ((scalar *)&s)->i >= 0; s = *++_i97) {
      if (fread (&a, sizeof (double), 1, p.fp) != 1) {
 fprintf (qstderr(), "input_gfs(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX) {
 if (_attribute[s.i].v.x.i >= 0) {



   if (_attribute[s.i].v.x.i == s.i) {
     s = _attribute[s.i].v.y;
     val(s,0,0,0) = a;
   }
   else if (_attribute[s.i].v.y.i == s.i) {
     s = _attribute[s.i].v.x;
     val(s,0,0,0) = - a;
   }





 }
 else
   val(s,0,0,0) = a;
      }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  boundary (listm);
  boundary (input);

  pfree (input,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);


  while (t < t1 && events (false))
    t = tnext;
  events (false);
 end_trace("input_gfs", "/home/jiarongw/basilisk/src/input.h", 328); }
#line 357 "/home/jiarongw/basilisk/src/input.h"
struct InputGRD {
  scalar s;
  FILE * fp;
  char * file;
  double nodatavalue;
  bool linear;
};

void input_grd (struct InputGRD p)
{
  scalar input = p.s;

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdin;
    else if (!(p.fp = fopen (p.file, "r"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }


  double DeltaGRD;
  int nx, ny;
  double XG0, YG0, ndv;


  char waste[100];
  fscanf (p.fp, "%s %d", waste, &nx);
  fscanf (p.fp, "%s %d", waste, &ny);
  fscanf (p.fp, "%s %lf", waste, &XG0);
  fscanf (p.fp, "%s %lf", waste, &YG0);
  fscanf (p.fp, "%s %lf", waste, &DeltaGRD);
  fscanf (p.fp, "%s %lf", waste, &ndv);


  if (!p.nodatavalue)
    p.nodatavalue = ndv;


  double * value = ((double *) pmalloc ((nx*ny)*sizeof(double),__func__,__FILE__,__LINE__));
  for (int i = ny - 1; i >= 0; i--)
    for (int j = 0 ; j < nx; j++)
      fscanf (p.fp, "%lf ", &value[j + i*nx]);

  double LGx0 = nx*DeltaGRD;
  double LGy0 = ny*DeltaGRD;
  bool warning = false;

   { foreach(){

#line 409 "/home/jiarongw/basilisk/src/input.h"
 {

    bool incl = (x >= XG0 - DeltaGRD/2. &&
   x <= XG0 + LGx0 + DeltaGRD/2. &&
   y >= YG0 - DeltaGRD/2. &&
   y <= YG0 + LGy0 + DeltaGRD/2.);
    if (incl) {
      double val;

      bool ring = (x >= XG0 &&
     x <= XG0 + LGx0 &&
     y >= YG0 &&
     y <= YG0 + LGy0);
      if (p.linear && ring ) {
 int j = (x - XG0)/DeltaGRD; int i = (y - YG0)/DeltaGRD;
 double dx = x -(j*DeltaGRD + XG0); double dy = y - (i*DeltaGRD + YG0);
 val = value[j + i*nx]
   + dx*(value[j + 1 + i*nx] - value[j + i*nx])/DeltaGRD
   + dy*(value[j + (i + 1)*nx] - value[j + i*nx])/DeltaGRD
   + dx*dy*(value[j + i*nx] + value[j +1 + (i+1)*nx]
     - value[j + (i + 1)*nx]-value[j + 1 + i*nx])/sq(DeltaGRD);
      }
      else {
 int j = (x - XG0 + DeltaGRD/2.)/DeltaGRD;
 int i = (y - YG0 + DeltaGRD/2.)/DeltaGRD;
 val = value[j + i*nx];
      }
      val(input,0,0,0) = val;
      if (val == ndv)
 val(input,0,0,0) = p.nodatavalue;
    }
    else {
      val(input,0,0,0) = p.nodatavalue;
      warning = true;
    }
  } } end_foreach(); }
  pfree (value,__func__,__FILE__,__LINE__);

  if (warning)
    fprintf (qstderr(),
      "input_grd(): Warning: Raster data is not covering all"
      " the simulation area\n");

  if (opened)
    fclose (p.fp);
}
#line 76 "/home/jiarongw/basilisk/src/view.h"






struct _bview {
  float tx, ty, sx, sy, sz;
  float quat[4];
  float fov;

  bool gfsview;
  bool vector;

  float bg[3];
  float lc;
  float res;

  unsigned width, height, samples;

  framebuffer * fb;
  Frustum frustum;

  void (* map) (coord *);

  int ni;

  bool active;
};

typedef struct _bview bview;




bview * bview_new() {
  bview * p = ((bview *) pcalloc (1, sizeof(bview),__func__,__FILE__,__LINE__));

  p->tx = p->ty = 0;
  p->sx = p->sy = p->sz = 1.;
  p->quat[0] = p->quat[1] = p->quat[2] = 0; p->quat[3] = 1;
  p->fov = 24.;
  gl_trackball (p->quat, 0.0, 0.0, 0.0, 0.0);


  p->bg[0] = 1; p->bg[1] = 1; p->bg[2] = 1;



  p->res = 1.;
  p->lc = 0.001;

  p->vector = false;

  p->samples = 4;
  p->width = 600*p->samples, p->height = 600*p->samples;

  p->fb = framebuffer_new (p->width, p->height);

  init_gl();
  p->active = false;

  return p;
}




void bview_destroy (bview * p)
{
  framebuffer_destroy (p->fb);
  pfree (p,__func__,__FILE__,__LINE__);
}




static bview * _view = NULL;






static void destroy_view()
{
  assert (_view);
  bview_destroy (_view);
}

bview * get_view() {
  if (!_view) {
    _view = bview_new();
    free_solver_func_add (destroy_view);
  }
  return _view;
}




static void redraw() {
  bview * view = get_view();


  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  double max = 2.;





  gluPerspective (view->fov, view->width/(float)view->height, 1., 1. + 2.*max);
  glMatrixMode (GL_MODELVIEW);

  glLoadIdentity ();
  glTranslatef (view->tx, view->ty, - (1. + max));

  GLfloat m[4][4];
  gl_build_rotmatrix (m, view->quat);
  glMultMatrixf (&m[0][0]);

  if (view->gfsview) {
    m[0][0] = 0., m[0][1] = 0., m[0][2] = -1.;
    m[1][0] = 0., m[1][1] = -1., m[1][2] = 0.;
    m[2][0] = 1., m[2][1] = 0., m[2][2] = 0.;
    glMultMatrixf (&m[0][0]);
  }

  glScalef (view->sx/L0, view->sy/L0, view->sz/L0);

  glClearColor (view->bg[0], view->bg[1], view->bg[2], 0.);
  glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  gl_get_frustum (&view->frustum);

  view->active = true;
  view->ni = 0;
}




bview * draw() {
  bview * view = get_view();
  if (!view->active)
    redraw();
  return view;
}







typedef void * pointer;
#line 244 "/home/jiarongw/basilisk/src/view.h"
typedef struct {
  GLubyte a[4];
} RGBA;

static void compose_image_op (void * pin, void * pout, int * len,
         MPI_Datatype * dptr)
{
  RGBA * in = pin, * out = pout;
  for (int i = 0; i < *len; i++,in++,out++)
    if (out->a[3] == 0)
      *out = *in;
}


static pointer compose_image (bview * view)
{ trace ("compose_image", "/home/jiarongw/basilisk/src/view.h", 259);
  unsigned char * image = framebuffer_image (view->fb);
  if (npe() > 1) {
    MPI_Op op;
    MPI_Op_create (compose_image_op, true, &op);
    MPI_Datatype rgba;
    MPI_Type_contiguous (4, MPI_BYTE, &rgba);
    MPI_Type_commit (&rgba);
    int size = view->width*view->height;
    if (pid() == 0)
      MPI_Reduce (MPI_IN_PLACE, image, size, rgba, op, 0, MPI_COMM_WORLD);
    else
      MPI_Reduce (image, image, size, rgba, op, 0, MPI_COMM_WORLD);
    MPI_Op_free (&op);
    MPI_Type_free (&rgba);
  }
  { pointer _ret =  image; end_trace("compose_image", "/home/jiarongw/basilisk/src/view.h", 275);  return _ret; }
 end_trace("compose_image", "/home/jiarongw/basilisk/src/view.h", 276); }
#line 339 "/home/jiarongw/basilisk/src/view.h"
#line 1 "draw.h"
#line 1 "/home/jiarongw/basilisk/src/draw.h"





#line 1 "gl/font.h"
#line 1 "/home/jiarongw/basilisk/src/gl/font.h"
void gl_StrokeCharacter( int character );
void gl_StrokeString( const char *string );
float gl_StrokeWidth( int character );
float gl_StrokeLength( const char *string );
GLfloat gl_StrokeHeight( );
#line 7 "/home/jiarongw/basilisk/src/draw.h"




void clear()
{
  bview * view = get_view();
  if (view->active)
    view->active = false;
  draw();
}
#line 46 "/home/jiarongw/basilisk/src/draw.h"
struct _view_set {
  float tx, ty;
  float fov;
  float quat[4];
  float sx, sy, sz;
  unsigned width, height, samples;
  float bg[3];
  float theta, phi, psi;
  bool relative;
  float res;
  char * camera;
  void (* map) (coord *);
  float p1x, p1y, p2x, p2y;
  bview * view;
};

void view (struct _view_set p)
{
  bview * v = p.view ? p.view : get_view();
  if (p.fov) {
    if (p.relative)
      v->fov += (0.1 + 3.*v->fov)*p.fov;
    else
      v->fov = p.fov;
    v->fov = clamp(v->fov,0.01,100.);
  }
  for (int i = 0; i < 4; i++)
    if (p.quat[i]) {
      for (int j = 0; j < 4; j++)
 v->quat[j] = p.quat[j];
      break;
    }
  if (p.tx) v->tx = p.relative ? v->tx + p.tx*0.02*(0.01 + 3.*v->fov) : p.tx;
  if (p.ty) v->ty = p.relative ? v->ty + p.ty*0.02*(0.01 + 3.*v->fov) : p.ty;
  if (p.sx) v->sx = p.sx;
  if (p.sy) v->sy = p.sy;
  if (p.sz) v->sz = p.sz;
  if (p.bg[0] || p.bg[1] || p.bg[2])
    for (int i = 0; i < 3; i++)
      v->bg[i] = p.bg[i];

  if (p.camera) {
    v->gfsview = false;
    if (strlen(p.camera) >= 4 &&
 !strcmp (&p.camera[strlen(p.camera) - 4], ".gfv")) {
      FILE * fp = fopen (p.camera, "r");
      if (!fp) {
 perror (p.camera);
 exit (1);
      }
      char s[81];
      float q[4], fov;
      int nq = 0, nf = 0;
      while (fgets (s, 81, fp) && (!nq || !nf)) {
 if (!nq)
   nq = sscanf (s, "  q0 = %f q1 = %f q2 = %f q3 = %f",
         &q[0], &q[1], &q[2], &q[3]);
 if (!nf)
   nf = sscanf (s, "  fov = %f", &fov);
      }
      if (nq != 4 || nf != 1) {
 fprintf (qstderr(), "%s: not a valid gfv file\n", p.camera);
 exit (1);
      }
      for (int j = 0; j < 4; j++)
 v->quat[j] = q[j];
      v->fov = fov;
      v->gfsview = true;
    }
    else if (!strcmp (p.camera, "left"))
      gl_axis_to_quat ((float[]){0,1,0}, - pi/2., v->quat);
    else if (!strcmp (p.camera, "right"))
      gl_axis_to_quat ((float[]){0,1,0}, pi/2., v->quat);
    else if (!strcmp (p.camera, "top"))
      gl_axis_to_quat ((float[]){1,0,0}, - pi/2., v->quat);
    else if (!strcmp (p.camera, "bottom"))
      gl_axis_to_quat ((float[]){1,0,0}, pi/2., v->quat);
    else if (!strcmp (p.camera, "front"))
      gl_axis_to_quat ((float[]){0,0,1}, 0., v->quat);
    else if (!strcmp (p.camera, "back"))
      gl_axis_to_quat ((float[]){0,1,0}, pi, v->quat);
    else if (!strcmp (p.camera, "iso")) {
      gl_axis_to_quat ((float[]){0,1,0}, pi/4., v->quat);
      float q[4];
      gl_axis_to_quat ((float[]){1,0,0}, - pi/4., q);
      gl_add_quats(q, v->quat, v->quat);
    }
    else {
      fprintf (qstderr(), "view(): unknown camera '%s'\n", p.camera);
      exit (1);
    }
  }
  else if (p.theta || p.phi || p.psi) {
    v->gfsview = false;
    float q[4];
    gl_axis_to_quat ((float[]){1,0,0}, - p.phi, q);
    if (p.relative) {
      float q1[4];
      gl_axis_to_quat ((float[]){0,1,0}, p.theta, q1);
      gl_add_quats(q, q1, q1);
      float q2[4];
      gl_axis_to_quat ((float[]){0,0,1}, p.psi, q2);
      gl_add_quats(q1, q2, q2);
      gl_add_quats(q2, v->quat, v->quat);
    }
    else {
      gl_axis_to_quat ((float[]){0,1,0}, p.theta, v->quat);
      gl_add_quats(q, v->quat, v->quat);
      gl_axis_to_quat ((float[]){0,0,1}, p.psi, q);
      gl_add_quats(q, v->quat, v->quat);
    }
  }

  if (p.map)
    v->map = p.map;

  if (p.p1x || p.p1y || p.p2x || p.p2y) {
    float q[4];
    gl_trackball(q, p.p1x, p.p1y, p.p2x, p.p2y);
    gl_add_quats (q, v->quat, v->quat);
  }

  if (p.res)
    v->res = p.res;

  if ((p.width && p.width != v->width) ||
      (p.height && p.height != v->height) ||
      (p.samples && p.samples != v->samples)) {
    v->width = v->width/v->samples;
    v->height = v->height/v->samples;
    if (p.width) v->width = p.width;
    if (p.height) v->height = p.height;
    if (p.samples) v->samples = p.samples;
    v->width *= v->samples;
    v->height *= v->samples;
    framebuffer_destroy (v->fb);
    v->fb = framebuffer_new (v->width, v->height);
    init_gl();
  }

  clear();
}







struct _translate {
  float x, y, z;
};

void begin_translate (struct _translate p)
{
  bview * view = draw();
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  glTranslatef (p.x, p.y, p.z);
  gl_get_frustum (&view->frustum);
}

void end_translate()
{
  bview * view = draw();
  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();
  gl_get_frustum (&view->frustum);
}
#line 224 "/home/jiarongw/basilisk/src/draw.h"
struct _mirror {
  coord n;
  double alpha;
};

void begin_mirror (struct _mirror p)
{
  bview * view = draw();
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  normalize (&p.n);
  GLfloat s[16], t[16];
  s[0] = 1. - 2.*p.n.x*p.n.x;
  s[1] = - 2.*p.n.x*p.n.y; s[2] = - 2.*p.n.x*p.n.z;
  s[3] = 0.;
  s[4] = s[1];
  s[5] = 1. - 2.*p.n.y*p.n.y; s[6] = - 2.*p.n.y*p.n.z;
  s[7] = 0.;
  s[8] = s[2]; s[9] = s[6]; s[10] = 1. - 2.*p.n.z*p.n.z;
  s[11] = 0.;
  s[12] = 0.; s[13] = 0.; s[14] = 0.;
  s[15] = 1.;

  t[0] = -1.; t[1] = 0.; t[2] = 0.; t[3] = 0.;
  t[4] = 0.; t[5] = -1.; t[6] = 0.; t[7] = 0.;
  t[8] = 0.; t[9] = 0.; t[10] = -1.; t[11] = 0.;
  t[12] = - 2.*p.n.x*p.alpha;
  t[13] = - 2.*p.n.y*p.alpha;
  t[14] = - 2.*p.n.z*p.alpha;
  t[15] = 1.;
  matrix_multiply (s, t);
  glMultMatrixf (s);
  gl_get_frustum (&view->frustum);
}

void end_mirror() {
  end_translate();
}







static void mapped_position (bview * view, coord * p, double * r)
{
  double x = p->x, y = p->y, z = p->z, rm = 0.;
  view->map (p);
  for (int i = -1; i <= 1; i += 2)
    for (int j = -1; j <= 1; j += 2)
      for (int k = -1; k <= 1; k += 2) {
 coord q = {x + i**r, y + j**r, z + k**r};
 view->map (&q);
 double pq = sq(p->x - q.x) + sq(p->y - q.y) + sq(p->z - q.z);
 if (pq > rm)
   rm = pq;
      }
  *r = sqrt (rm);
}

#define foreach_visible(view)\
foreach_cell() {\
\
  double _r = Delta*0.71;\
\
\
\
  coord _p = {x, y, z};\
  if ((view)->map)\
    mapped_position (view, &_p, &_r);\
  if (!sphere_in_frustum (_p.x, _p.y, _p.z, _r, &(view)->frustum))\
    continue;\
  if (is_leaf(cell) ||\
      sphere_diameter (_p.x, _p.y, _p.z, _r/L0, &(view)->frustum)\
      < (view)->res) {\
    if (is_active(cell) && is_local(cell)) {\

#line 301

#define end_foreach_visible()\
    }\
    continue;\
  }\
}\
end_foreach_cell();\

#line 308

#line 348 "/home/jiarongw/basilisk/src/draw.h"
static scalar lookup_field (const char * name)
{
  if (name)
    if (all) for (scalar s = *all, *_i98 = all; ((scalar *)&s)->i >= 0; s = *++_i98)
      if (!strcmp (_attribute[s.i].name, name))
 return s;
  return (scalar){-1};
}

static vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    if (all) for (scalar s = *all, *_i99 = all; ((scalar *)&s)->i >= 0; s = *++_i99)
      if (!strcmp (_attribute[s.i].name, component))
 return _attribute[s.i].v;
  }
  return (vector){{-1}};
}

static void draw_lines (bview * view, float color[3], float lw)
{
  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glTranslatef (0., 0., view->lc*view->fov/24.);
  glColor3f (color[0], color[1], color[2]);
  glLineWidth (view->samples*(lw > 0. ? lw : 1.));
}

static inline double interp (Point point, coord p, scalar col) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 379 "/home/jiarongw/basilisk/src/draw.h"

  struct _interpolate _r = { col, x + p.x*Delta, y + p.y*Delta, z + p.z*Delta };
  return interpolate_linear (point, _r);
}
#line 439 "/home/jiarongw/basilisk/src/draw.h"
static void begin_colorized (float fc[3],
        double cmap[127][3], bool use_texture)
{

  if (use_texture) {
    GLfloat texture[3*256];
    for (int i = 0; i < 256; i++) {
      color c = colormap_color (cmap, i/255., 0, 1);
      texture[3*i] = c.r/255.;
      texture[3*i + 1] = c.g/255.;
      texture[3*i + 2] = c.b/255.;
    }
    glTexImage1D (GL_TEXTURE_1D, 0, GL_RGB, 256,0, GL_RGB, GL_FLOAT, texture);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glEnable (GL_TEXTURE_1D);
  }
  glColor3f (fc[0], fc[1], fc[2]);
}

static void end_colorized() {
  glDisable (GL_TEXTURE_1D);
}
#line 493 "/home/jiarongw/basilisk/src/draw.h"
struct _draw_vof {
  char * c;
  char * s;
  bool edges;
  double larger;
  int filled;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3], lw;
};







static bool cfilter (Point point, scalar c, double cmin)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 514 "/home/jiarongw/basilisk/src/draw.h"

  double cmin1 = 4.*cmin;
  if (val(c,0,0,0) <= cmin) {
    {
#line 517

      if (val(c,1,0,0) >= 1. - cmin1 || val(c,-1,0,0) >= 1. - cmin1)
 return true;
#line 517

      if (val(c,0,1,0) >= 1. - cmin1 || val(c,0,-1,0) >= 1. - cmin1)
 return true;}
    return false;
  }
  if (val(c,0,0,0) >= 1. - cmin) {
    {
#line 523

      if (val(c,1,0,0) <= cmin1 || val(c,-1,0,0) <= cmin1)
 return true;
#line 523

      if (val(c,0,1,0) <= cmin1 || val(c,0,-1,0) <= cmin1)
 return true;}
    return false;
  }
  int n = 0;
  double min = HUGE, max = - HUGE;
   { foreach_neighbor(1) {
    if (val(c,0,0,0) > cmin && val(c,0,0,0) < 1. - cmin && ++n >= (1 << 2))
      return true;
    if (val(c,0,0,0) > max) max = val(c,0,0,0);
    if (val(c,0,0,0) < min) min = val(c,0,0,0);
  } end_foreach_neighbor(); }
  return max - min > 0.5;
}


static void glvertex2d (bview * view, double x, double y) {
  if (view->map) {
    coord p = {x, y, 0.};
    view->map (&p);
    glVertex2d (p.x, p.y);
  }
  else
    glVertex2d (x, y);
}
#line 561 "/home/jiarongw/basilisk/src/draw.h"

bool draw_vof (struct _draw_vof p)
{ trace ("draw_vof", "/home/jiarongw/basilisk/src/draw.h", 563);
  scalar c = lookup_field (p.c);
  if (c.i < 0) {
    fprintf (qstderr(), "draw_vof(): no field named '%s'\n", p.c);
    { bool _ret =  false; end_trace("draw_vof", "/home/jiarongw/basilisk/src/draw.h", 567);  return _ret; }
  }
  vector s = lookup_vector (p.s);

  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = lookup_field (p.color); if (col.i < 0) { fprintf (qstderr(), "colorize_args(): no field named '%s'\n", p.color); { bool _ret =  false; end_trace("draw_vof", "/home/jiarongw/basilisk/src/draw.h", 571);  return _ret; } } } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { double spread = (p.spread ? p.spread : 5.)*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((2 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.;;

  double cmin = 1e-3;



  if (_attribute[c.i].prolongation != fraction_refine) {
    _attribute[c.i].prolongation = _attribute[c.i].refine = fraction_refine;
    boundary (((scalar []){c,{-1}}));
  }


  bview * view = draw();

  if (p.filled) {
    glColor3f (p.fc[0], p.fc[1], p.fc[2]);
     { foreach_visible (view){

#line 587 "/home/jiarongw/basilisk/src/draw.h"
 {
      if ((p.filled > 0 && val(c,0,0,0) >= 1.) || (p.filled < 0 && val(c,0,0,0) <= 0.)) {
 glBegin (GL_QUADS);
 glvertex2d (view, x - Delta/2., y - Delta/2.);
 glvertex2d (view, x + Delta/2., y - Delta/2.);
 glvertex2d (view, x + Delta/2., y + Delta/2.);
 glvertex2d (view, x - Delta/2., y + Delta/2.);
 glEnd();
 view->ni++;
      }
      else if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
 coord n = facet_normal (point, c, s), s = {1.,1.};
 if (p.filled < 0)
   {
#line 600

     n.x = - n.x;
#line 600

     n.y = - n.y;}
 double alpha = line_alpha (p.filled < 0. ? 1. - val(c,0,0,0) : val(c,0,0,0), n);
 alpha += (n.x + n.y)/2.;
 {
#line 604

   if (n.x < 0.) alpha -= n.x, n.x = - n.x, s.x = - 1.;
#line 604

   if (n.y < 0.) alpha -= n.y, n.y = - n.y, s.y = - 1.;}
 coord v[5];
 int nv = 0;
 if (alpha >= 0. && alpha <= n.x) {
   v[nv].x = alpha/n.x, v[nv++].y = 0.;
   if (alpha <= n.y)
     v[nv].x = 0., v[nv++].y = alpha/n.y;
   else if (alpha >= n.y && alpha - n.y <= n.x) {
     v[nv].x = (alpha - n.y)/n.x, v[nv++].y = 1.;
     v[nv].x = 0., v[nv++].y = 1.;
   }
   v[nv].x = 0., v[nv++].y = 0.;
 }
 else if (alpha >= n.x && alpha - n.x <= n.y) {
   v[nv].x = 1., v[nv++].y = (alpha - n.x)/n.y;
   if (alpha >= n.y && alpha - n.y <= n.x) {
     v[nv].x = (alpha - n.y)/n.x, v[nv++].y = 1.;
     v[nv].x = 0., v[nv++].y = 1.;
   }
   else if (alpha <= n.y)
     v[nv].x = 0., v[nv++].y = alpha/n.y;
   v[nv].x = 0., v[nv++].y = 0.;
   v[nv].x = 1., v[nv++].y = 0.;
 }
 glBegin (GL_POLYGON);
 if (s.x*s.y < 0.)
   for (int i = nv - 1; i >= 0; i--)
     glvertex2d (view, x + s.x*(v[i].x - 0.5)*Delta,
   y + s.y*(v[i].y - 0.5)*Delta);
 else
   for (int i = 0; i < nv; i++)
     glvertex2d (view, x + s.x*(v[i].x - 0.5)*Delta,
   y + s.y*(v[i].y - 0.5)*Delta);
 glEnd ();
 view->ni++;
      }
    } } end_foreach_visible(); }
  }
  if (!p.filled || p.edges) {
    draw_lines (view, p.lc, p.lw);
    glBegin (GL_LINES);
     { foreach_visible (view){

#line 646 "/home/jiarongw/basilisk/src/draw.h"

      if (cfilter (point, c, cmin)) {
 coord n = facet_normal (point, c, s);
 double alpha = line_alpha (val(c,0,0,0), n);
 coord segment[2];
 if (facets (n, alpha, segment) == 2) {
   glvertex2d (view, x + segment[0].x*Delta, y + segment[0].y*Delta);
   glvertex2d (view, x + segment[1].x*Delta, y + segment[1].y*Delta);
   view->ni++;
 }
      } } end_foreach_visible(); }
    glEnd ();
    glPopMatrix ();
  }
#line 708 "/home/jiarongw/basilisk/src/draw.h"
  { bool _ret =  true; end_trace("draw_vof", "/home/jiarongw/basilisk/src/draw.h", 708);  return _ret; }
 end_trace("draw_vof", "/home/jiarongw/basilisk/src/draw.h", 709); }
#line 722 "/home/jiarongw/basilisk/src/draw.h"
struct _cells {
  coord n;
  double alpha;
  float lc[3], lw;
};


void cells (struct _cells p)
{ trace ("cells", "/home/jiarongw/basilisk/src/draw.h", 730);
  bview * view = draw();
  draw_lines (view, p.lc, p.lw);

   { foreach_visible (view){

#line 734 "/home/jiarongw/basilisk/src/draw.h"
 {
    glBegin (GL_LINE_LOOP);
    glvertex2d (view, x - Delta/2., y - Delta/2.);
    glvertex2d (view, x + Delta/2., y - Delta/2.);
    glvertex2d (view, x + Delta/2., y + Delta/2.);
    glvertex2d (view, x - Delta/2., y + Delta/2.);
    glEnd();
    view->ni++;
  } } end_foreach_visible(); }
#line 756 "/home/jiarongw/basilisk/src/draw.h"
  glPopMatrix ();
 end_trace("cells", "/home/jiarongw/basilisk/src/draw.h", 757); }
#line 774 "/home/jiarongw/basilisk/src/draw.h"
struct _squares {
  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3];

  coord n;
  double alpha;
};


bool squares (struct _squares p)
{ trace ("squares", "/home/jiarongw/basilisk/src/draw.h", 787);
  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = lookup_field (p.color); if (col.i < 0) { fprintf (qstderr(), "colorize_args(): no field named '%s'\n", p.color); { bool _ret =  false; end_trace("squares", "/home/jiarongw/basilisk/src/draw.h", 788);  return _ret; } } } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { double spread = (p.spread ? p.spread : 5.)*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((2 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.;;
  scalar f = col;

  bview * view = draw();
  glShadeModel (GL_SMOOTH);
  if (p.linear) {
    { begin_colorized (p.fc, cmap, !view->vector && p.color && p.linear && col.i >= 0); {

       { foreach_visible (view){

#line 796 "/home/jiarongw/basilisk/src/draw.h"

        if (val(f,0,0,0) != nodata) {
   glBegin (GL_TRIANGLE_FAN);
   if (p.color && p.linear && col.i >= 0) { if (view->vector) { color b = colormap_color (cmap, (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16., p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } }

                                                  ;
   glvertex2d (view, x, y);
   if (p.color && p.linear && col.i >= 0) { if (view->vector) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
   glvertex2d (view, x - Delta/2., y - Delta/2.);
   if (p.color && p.linear && col.i >= 0) { if (view->vector) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
   glvertex2d (view, x + Delta/2., y - Delta/2.);
   if (p.color && p.linear && col.i >= 0) { if (view->vector) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4., p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
   glvertex2d (view, x + Delta/2., y + Delta/2.);
   if (p.color && p.linear && col.i >= 0) { if (view->vector) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4., p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
   glvertex2d (view, x - Delta/2., y + Delta/2.);
   if (p.color && p.linear && col.i >= 0) { if (view->vector) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
   glvertex2d (view, x - Delta/2., y - Delta/2.);
   glEnd();
   view->ni++;
 } } end_foreach_visible(); }
#line 842 "/home/jiarongw/basilisk/src/draw.h"
    } end_colorized(); }
  }
  else {

    glBegin (GL_QUADS);
     { foreach_visible (view){

#line 847 "/home/jiarongw/basilisk/src/draw.h"

      if (val(f,0,0,0) != nodata) {
 if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x - Delta/2., y - Delta/2.);
 glvertex2d (view, x + Delta/2., y - Delta/2.);
 glvertex2d (view, x + Delta/2., y + Delta/2.);
 glvertex2d (view, x - Delta/2., y + Delta/2.);
 view->ni++;
      } } end_foreach_visible(); }
    glEnd();
#line 873 "/home/jiarongw/basilisk/src/draw.h"
  }
  { bool _ret =  true; end_trace("squares", "/home/jiarongw/basilisk/src/draw.h", 874);  return _ret; }
 end_trace("squares", "/home/jiarongw/basilisk/src/draw.h", 875); }
#line 886 "/home/jiarongw/basilisk/src/draw.h"
struct _box {
  bool notics;
  float lc[3], lw;
};


bool box (struct _box p)
{ trace ("box", "/home/jiarongw/basilisk/src/draw.h", 893);
  bview * view = draw();
  draw_lines (view, p.lc, p.lw);

  float height = 0.5*gl_StrokeHeight();
  float width = gl_StrokeWidth ('1'), scale = L0/(60.*width), length;
  float Z1 = 2 == 2 ? 0. : Z0;
  char label[80];

  glMatrixMode (GL_MODELVIEW);

  if (!p.notics) {
    int nt = 8;
    for (int i = 0; i <= nt; i++) {
      glPushMatrix();
      glTranslatef (X0 + i*L0/nt - height/2.*scale, Y0 - width/3.*scale, Z1);
      glRotatef (-90, 0, 0, 1);
      glScalef (scale, scale, scale);
      sprintf (label, "%g", X0 + i*L0/nt);
      gl_StrokeString (label);
      glPopMatrix();

      glPushMatrix();
      sprintf (label, "%g", Y0 + i*L0/nt);
      length = gl_StrokeLength (label);
      glTranslatef (X0 - (length + width/3.)*scale,
      Y0 + i*L0/nt - height/2.*scale, Z1);
      glScalef (scale, scale, scale);
      gl_StrokeString (label);
      glPopMatrix();
#line 935 "/home/jiarongw/basilisk/src/draw.h"
    }

    glPushMatrix();
    sprintf (label, "%g", X0 + L0/2.);
    length = gl_StrokeLength (label);
    glTranslatef (X0 + L0/2 - height*scale, Y0 - (length + 4.*width)*scale, Z1);
    glScalef (2.*scale, 2.*scale, 2.*scale);
    gl_StrokeString ("X");
    glPopMatrix();


    glPushMatrix();
    sprintf (label, "%g", Y0 + L0/2.);
    length = gl_StrokeLength (label);
    glTranslatef (X0 - (length + 4.*width)*scale,
    Y0 + L0/2. - height*scale, Z1);
    glScalef (2.*scale, 2.*scale, 2.*scale);
    gl_StrokeString ("Y");
    glPopMatrix();
#line 966 "/home/jiarongw/basilisk/src/draw.h"
  }


   { foreach_level (0){

#line 969 "/home/jiarongw/basilisk/src/draw.h"
 {
    glBegin (GL_LINE_LOOP);
    glvertex2d (view, x - Delta/2., y - Delta/2.);
    glvertex2d (view, x + Delta/2., y - Delta/2.);
    glvertex2d (view, x + Delta/2., y + Delta/2.);
    glvertex2d (view, x - Delta/2., y + Delta/2.);
    glEnd ();
    view->ni++;
  } } end_foreach_level(); }
#line 999 "/home/jiarongw/basilisk/src/draw.h"
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();
  { bool _ret =  true; end_trace("box", "/home/jiarongw/basilisk/src/draw.h", 1001);  return _ret; }
 end_trace("box", "/home/jiarongw/basilisk/src/draw.h", 1002); }
#line 1014 "/home/jiarongw/basilisk/src/draw.h"
struct _isosurface {
  char * f;
  double v;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3];
};


bool isosurface (struct _isosurface p)
{ trace ("isosurface", "/home/jiarongw/basilisk/src/draw.h", 1027);
#line 1075 "/home/jiarongw/basilisk/src/draw.h"
  { bool _ret =  true; end_trace("isosurface", "/home/jiarongw/basilisk/src/draw.h", 1075);  return _ret; }
 end_trace("isosurface", "/home/jiarongw/basilisk/src/draw.h", 1076); }
#line 1086 "/home/jiarongw/basilisk/src/draw.h"
struct _travelling {
  double start, end;
  float tx, ty, quat[4], fov;
};




void travelling (struct _travelling p)
{
  static float tx, ty, quat[4], fov;
  static double told = -1.;
  if (told < p.start && t >= p.start) {
    bview * view = get_view();
    tx = view->tx, ty = view->ty, fov = view->fov;
    for (int i = 0; i < 4; i++)
      quat[i] = view->quat[i];
  }
  if (t >= p.start && t <= p.end)
    view ((struct _view_set){.tx = (!p.tx ? tx : ((t - p.start)*(p.tx) + (p.end - t)*(tx))/(p.end - p.start)), .ty = (!p.ty ? ty : ((t - p.start)*(p.ty) + (p.end - t)*(ty))/(p.end - p.start)),
   .fov = (!p.fov ? fov : ((t - p.start)*(p.fov) + (p.end - t)*(fov))/(p.end - p.start)),
   .quat = {(!p.quat[0] ? quat[0] : ((t - p.start)*(p.quat[0]) + (p.end - t)*(quat[0]))/(p.end - p.start)), (!p.quat[1] ? quat[1] : ((t - p.start)*(p.quat[1]) + (p.end - t)*(quat[1]))/(p.end - p.start)),
           (!p.quat[2] ? quat[2] : ((t - p.start)*(p.quat[2]) + (p.end - t)*(quat[2]))/(p.end - p.start)), (!p.quat[3] ? quat[3] : ((t - p.start)*(p.quat[3]) + (p.end - t)*(quat[3]))/(p.end - p.start))}});
  if (told < p.end && t >= p.end) {
    bview * view = get_view();
    tx = view->tx, ty = view->ty, fov = view->fov;
    for (int i = 0; i < 4; i++)
      quat[i] = view->quat[i];
  }
  told = t;
}
#line 1133 "/home/jiarongw/basilisk/src/draw.h"
struct _draw_string {
  char * str;
  int pos;
  float size;
  float lc[3], lw;
};


bool draw_string (struct _draw_string p)
{ trace ("draw_string", "/home/jiarongw/basilisk/src/draw.h", 1142);
  bview * view = draw();

  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glColor3f (p.lc[0], p.lc[1], p.lc[2]);
  glLineWidth (view->samples*(p.lw > 0. ? p.lw : 1.));

  float width = gl_StrokeWidth ('1'), height = gl_StrokeHeight();
  if (!p.size)
    p.size = 40;
  float hscale = 2./(p.size*width), vscale = hscale*view->width/view->height;
  float vmargin = width/2.*vscale;
  if (p.pos == 0)
    glTranslatef (-1., -1. + vmargin, 0.);
  else if (p.pos == 1)
    glTranslatef (-1., 1. - height*vscale, 0.);
  else if (p.pos == 2)
    glTranslatef (1. - strlen(p.str)*width*hscale, 1. - height*vscale, 0.);
  else
    glTranslatef (1. - strlen(p.str)*width*hscale, -1. + vmargin, 0.);
  glScalef (hscale, vscale, 1.);
  gl_StrokeString (p.str);

  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();

  { bool _ret =  true; end_trace("draw_string", "/home/jiarongw/basilisk/src/draw.h", 1177);  return _ret; }
 end_trace("draw_string", "/home/jiarongw/basilisk/src/draw.h", 1178); }
#line 340 "/home/jiarongw/basilisk/src/view.h"
#line 360 "/home/jiarongw/basilisk/src/view.h"
struct _load {
  FILE * fp;
  char * file;
  Array * buf;
  Array * history;
};

bool load (struct _load p);
#line 411 "/home/jiarongw/basilisk/src/view.h"
struct _save {
  char * file, * format, * opt;
  FILE * fp;
  float lw;
  int sort, options;

  Array * history;
  bview * view;
};

static void bview_draw (bview * view)
{
  if (!view->active)
    return;
  view->active = false;
  glFinish ();
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}

static void redraw_feedback (struct _save * p)
{
  bview * view = p->view ? p->view : get_view();
  assert (p->history);
  if (p->history->len) {
    float res = view->res;
    view->res = 0.;
    view->vector = true;
    redraw();


    int list = glGenLists (1);
    glNewList (list, GL_COMPILE);
    load ((struct _load){.buf = p->history});
    glEndList();
    glCallList (list);
    glFinish ();
    glDeleteLists (list, 1);
    enable_fpe (FE_DIVBYZERO|FE_INVALID);
    view->active = false;
    view->vector = false;
    view->res = res;
  }
}




bool save (struct _save p)
{ trace ("save", "/home/jiarongw/basilisk/src/view.h", 459);
  char ppm[] = "ppm";
  if (!p.format) {
    p.format = ppm;
    if (p.file) {
      char * s = strchr (p.file, '.'), * dot = s;
      while (s) {
 dot = s;
 s = strchr (s + 1, '.');
      }
      if (dot)
 p.format = dot + 1;
    }
  }

  bview * view = p.view ? p.view : get_view();

  if (!strcmp (p.format, "png") ||
      !strcmp (p.format, "jpg") ||
      (p.file && is_animation (p.file))) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0) {
      FILE * fp = open_image (p.file, p.opt);
      gl_write_image (fp, image, view->width, view->height, view->samples);
      close_image (p.file, fp);
    }
    { bool _ret =  true; end_trace("save", "/home/jiarongw/basilisk/src/view.h", 486);  return _ret; }
  }

  if (p.file && (p.fp = fopen (p.file, "w")) == NULL) {
    perror (p.file);
    { bool _ret =  false; end_trace("save", "/home/jiarongw/basilisk/src/view.h", 491);  return _ret; }
  }
  if (!p.fp)
    p.fp = qstdout();

  if (!strcmp (p.format, "ppm")) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0)
      gl_write_image (p.fp, image, view->width, view->height, view->samples);
  }

  else if (!strcmp (p.format, "bv")) {
    assert (p.history);
    fprintf (p.fp,
      "view (fov = %g, quat = {%g,%g,%g,%g}, "
      "tx = %g, ty = %g, "
      "bg = {%g,%g,%g}, "
      "width = %d, height = %d, samples = %d"
      ");\n",
      view->fov,
      view->quat[0], view->quat[1], view->quat[2], view->quat[3],
      view->tx, view->ty,
      view->bg[0], view->bg[1], view->bg[2],
      view->width/view->samples, view->height/view->samples,
      view->samples);
    fwrite (p.history->p, 1, p.history->len, p.fp);
  }

  else if (!strcmp (p.format, "gnu") ||
    !strcmp (p.format, "obj") ||
    !strcmp (p.format, "kml")) {
    int format = (!strcmp (p.format, "gnu") ? FEEDBACK_GNU :
    !strcmp (p.format, "obj") ? FEEDBACK_OBJ :
    !strcmp (p.format, "kml") ? FEEDBACK_KML :
    -1);
    unsigned buffsize = 1 << 24;
    bool done = false;
    while (!done && buffsize <= (1 << 28)) {
      float * f = gl_feedback_begin (buffsize);
      redraw_feedback (&p);
      done = gl_feedback_end (f, p.fp, format);
      buffsize *= 2;
    }
    if (!done)
      fprintf (ferr, "save(): error: exceeded maximum feedback buffer size\n");
  }

  else if (!strcmp (p.format, "ps") ||
    !strcmp (p.format, "eps") ||
    !strcmp (p.format, "tex") ||
    !strcmp (p.format, "pdf") ||
    !strcmp (p.format, "svg") ||
    !strcmp (p.format, "pgf")) {
    GLint format = (!strcmp (p.format, "ps") ? GL2PS_PS :
      !strcmp (p.format, "eps") ? GL2PS_EPS :
      !strcmp (p.format, "tex") ? GL2PS_TEX :
      !strcmp (p.format, "pdf") ? GL2PS_PDF :
      !strcmp (p.format, "svg") ? GL2PS_SVG :
      !strcmp (p.format, "pgf") ? GL2PS_PGF :
      -1);
    GLint state = GL2PS_OVERFLOW;
    GLint sort = p.sort ? p.sort : GL2PS_SIMPLE_SORT;
    GLint options = p.options ? p.options : (GL2PS_SIMPLE_LINE_OFFSET |
          GL2PS_SILENT |
          GL2PS_BEST_ROOT |
          GL2PS_OCCLUSION_CULL |
          GL2PS_USE_CURRENT_VIEWPORT |
          GL2PS_TIGHT_BOUNDING_BOX);
    unsigned buffsize = 1 << 24;
    while (state == GL2PS_OVERFLOW && buffsize <= (1 << 28)) {
      gl2psBeginPage ("", "bview",
        NULL,
        format, sort, options,
        GL_RGBA, 0, NULL,
        0, 0, 0,
        buffsize, p.fp, "");
      redraw_feedback (&p);
      disable_fpe (FE_DIVBYZERO|FE_INVALID);
      state = gl2psEndPage();
      enable_fpe (FE_DIVBYZERO|FE_INVALID);
      buffsize *= 2;
    }
    if (state == GL2PS_OVERFLOW)
      fprintf (ferr, "save(): error: exceeded maximum feedback buffer size\n");
  }

  else {
    fprintf (ferr, "save(): unknown format '%s'\n", p.format);
    if (p.file) {
      fclose (p.fp);
      remove (p.file);
    }
    { bool _ret =  false; end_trace("save", "/home/jiarongw/basilisk/src/view.h", 584);  return _ret; }
  }

  fflush (p.fp);
  if (p.file)
    fclose (p.fp);

  { bool _ret =  true; end_trace("save", "/home/jiarongw/basilisk/src/view.h", 591);  return _ret; }
 end_trace("save", "/home/jiarongw/basilisk/src/view.h", 592); }







static char * remove_blanks (char * line)
{
  while (strchr (" \t", *line)) line++;
  char * s = line, * cur = line;
  bool instring = false;
  while (*s != '\0' && *s != '#') {
    if (*s == '"')
      instring = !instring;
    if (instring || !strchr (" \t", *s))
      *cur++ = *s;
    s++;
  }
  *cur = '\0';
  return line;
}

static void fields_stats()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  if (all) for (scalar s = *all, *_i100 = all; ((scalar *)&s)->i >= 0; s = *++_i100)
    fprintf (ferr, " %s", _attribute[s.i].name);
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
    "name", "min", "avg", "stddev", "max");
  if (all) for (scalar s = *all, *_i101 = all; ((scalar *)&s)->i >= 0; s = *++_i101) {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
      _attribute[s.i].name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }
}

static void draw_append (char * buf, Array * history, FILE * interactive)
{
  if (interactive) {
    if (history->len)
      load ((struct _load){.buf = history});
    save ((struct _save){.fp = interactive});
  }
  array_append (history, buf, strlen(buf)*sizeof(char));
}






#line 1 "draw_get.h"
#line 1 "/home/jiarongw/basilisk/src/draw_get.h"


bool _view_set_get (struct _view_set * p) {
  Params params[] = {
    {"tx", pfloat, &p->tx},
    {"ty", pfloat, &p->ty},
    {"fov", pfloat, &p->fov},
    {"quat", pfloat, p->quat, 4},
    {"sx", pfloat, &p->sx},
    {"sy", pfloat, &p->sy},
    {"sz", pfloat, &p->sz},
    {"width", punsigned, &p->width},
    {"height", punsigned, &p->height},
    {"samples", punsigned, &p->samples},
    {"bg", pfloat, p->bg, 3},
    {"theta", pfloat, &p->theta},
    {"phi", pfloat, &p->phi},
    {"psi", pfloat, &p->psi},
    {"relative", pbool, &p->relative},
    {"res", pfloat, &p->res},
    {"camera", pstring, &p->camera},
    {"p2y", pfloat, &p->p2y},
    {"p1x", pfloat, &p->p1x},
    {"p1y", pfloat, &p->p1y},
    {"p2x", pfloat, &p->p2x},
    {NULL}
  };
  return parse_params (params);
}

bool _translate_get (struct _translate * p) {
  Params params[] = {
    {"x", pfloat, &p->x},
    {"y", pfloat, &p->y},
    {"z", pfloat, &p->z},
    {NULL}
  };
  return parse_params (params);
}

bool _mirror_get (struct _mirror * p) {
  Params params[] = {
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {NULL}
  };
  return parse_params (params);
}

bool _draw_vof_get (struct _draw_vof * p) {
  Params params[] = {
    {"c", pstring, &p->c},
    {"s", pstring, &p->s},
    {"edges", pbool, &p->edges},
    {"larger", pdouble, &p->larger},
    {"filled", pint, &p->filled},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _cells_get (struct _cells * p) {
  Params params[] = {
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _squares_get (struct _squares * p) {
  Params params[] = {
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {NULL}
  };
  return parse_params (params);
}

bool _box_get (struct _box * p) {
  Params params[] = {
    {"notics", pbool, &p->notics},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _isosurface_get (struct _isosurface * p) {
  Params params[] = {
    {"f", pstring, &p->f},
    {"v", pdouble, &p->v},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {NULL}
  };
  return parse_params (params);
}

bool _travelling_get (struct _travelling * p) {
  Params params[] = {
    {"start", pdouble, &p->start},
    {"end", pdouble, &p->end},
    {"fov", pfloat, &p->fov},
    {"tx", pfloat, &p->tx},
    {"ty", pfloat, &p->ty},
    {"quat", pfloat, p->quat, 4},
    {NULL}
  };
  return parse_params (params);
}

bool _draw_string_get (struct _draw_string * p) {
  Params params[] = {
    {"str", pstring, &p->str},
    {"pos", pint, &p->pos},
    {"size", pfloat, &p->size},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}
#line 647 "/home/jiarongw/basilisk/src/view.h"

static bool process_line (char * line, Array * history, FILE * interactive)
{
  if (line[0] == '\0')
    return true;
  char * buf = pstrdup (line,__func__,__FILE__,__LINE__);
  char * s = strtok (remove_blanks (line), "(");
  if (!s) {
    pfree (buf,__func__,__FILE__,__LINE__);
    return true;
  }

  if (!strcmp (s, "restore")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      if (!restore ((struct Dump){.file = file, .list = all}))
 fprintf (ferr, "could not restore from '%s'\n", file);
      else {
 restriction (all);
 fields_stats();
 clear();

 if (history->len && load ((struct _load){.buf = history}) && interactive)
   save ((struct _save){.fp = interactive});
      }
    }
  }

  else if (!strcmp (s, "dump")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    dump ((struct Dump){.file = file});
  }

  else if (!strcmp (s, "input_gfs")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      input_gfs ((struct OutputGfs){.file = file, .list = all});
      restriction (all);
      fields_stats();
      clear();

      if (history->len && load ((struct _load){.buf = history}) && interactive)
 save ((struct _save){.fp = interactive});
    }
  }

  else if (!strcmp (s, "save")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      save ((struct _save){.file = file, .history = history});
  }

  else if (!strcmp (s, "load")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file && load ((struct _load){.file = file, .history = history}) && interactive) {
      load ((struct _load){.buf = history});
      save ((struct _save){.fp = interactive});
    }
  }

  else if (!strcmp (s, "cells")) {
    struct _cells p = {{0}};
    _cells_get (&p);
    cells (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "draw_vof")) {
    struct _draw_vof p = {0};
    _draw_vof_get (&p);
    if (draw_vof (p))
      draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "squares")) {
    struct _squares p = {0};
    _squares_get (&p);
    squares (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "begin_translate")) {
    struct _translate p = {0};
    _translate_get (&p);
    begin_translate (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "end_translate")) {
    end_translate();
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "begin_mirror")) {
    struct _mirror p = {{0}};
    _mirror_get (&p);
    begin_mirror (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "end_mirror")) {
    end_mirror();
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "squares")) {
    struct _squares p = {0};
    _squares_get (&p);
    squares (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "isosurface")) {
    struct _isosurface p = {0};
    _isosurface_get (&p);
    isosurface (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "draw_string")) {
    struct _draw_string p = {0};
    _draw_string_get (&p);
    draw_string (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "display")) {
    if (interactive && history->len && load ((struct _load){.buf = history}))
      save ((struct _save){.fp = interactive});
  }

  else if (!strcmp (s, "clear")) {
    clear();
    if (interactive)
      save ((struct _save){.fp = interactive});
    history->len = 0;
  }

  else if (!strcmp (s, "show")) {
    if (interactive && history->len)
      save ((struct _save){.fp = ferr, .format = "bv", .history = history});
  }

  else if (!strcmp (s, "box")) {
    struct _box p = {0};
    _box_get (&p);
    box (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "view")) {
    struct _view_set p = {0};
    _view_set_get (&p);
    view (p);
    if (p.width || p.height || p.samples) {

      if (history->len && load ((struct _load){.buf = history}) && p.samples && interactive)
 save ((struct _save){.fp = interactive});
    }
  }

  else if (!strcmp (s, "quit")) {
    pfree (buf,__func__,__FILE__,__LINE__);
    return false;
  }

  else if (s[0] != '\n')
    fprintf (ferr, "load(): syntax error: '%s'\n", s);

  pfree (buf,__func__,__FILE__,__LINE__);
  return true;
}

bool load (struct _load p) {
  if (p.file) {
    p.fp = fopen (p.file, "r");
    if (!p.fp) {
      perror (p.file);
      return false;
    }
  }

  Array * history = array_new();
  if (p.fp) {
    char line[256];
    while (fgets (line, 256, p.fp) && process_line (line, history, NULL));
  }
  else if (p.buf) {
    int i = 0;
    char * s = (char *) p.buf->p;
    while (i < p.buf->len) {
      char * start = s;
      while (i < p.buf->len && *s != '\n')
 s++, i++;
      if (*s == '\n' && ++s > start) {
 char line[s - start + 1];
 strncpy (line, start, s - start);
 line[s - start] = '\0';
 process_line (line, history, NULL);
      }
    }
  }
  if (p.history)
    array_append (p.history, history->p, history->len);
  array_free (history);

  return true;
}
#line 15 "wavewind.c"
#line 1 "tag.h"
#line 1 "/home/jiarongw/basilisk/src/tag.h"
#line 13 "/home/jiarongw/basilisk/src/tag.h"
static void restriction_tag (Point point, scalar t)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 14 "/home/jiarongw/basilisk/src/tag.h"

  double min = HUGE;
   { foreach_child()
    if (val(t,0,0,0) < min)
      min = val(t,0,0,0); end_foreach_child(); }
  val(t,0,0,0) = min;
}







static long lookup_tag (Array * a, double tag)
{
  long len = a->len/sizeof(double);
  double * p = (double *) a->p;
  if (tag == p[0])
    return 0;
  if (tag == p[len - 1])
    return len - 1;

  long s = 0, e = len - 1;
  while (s < e - 1) {
    long m = (s + e)/2;
    if (p[m] <= tag)
      s = m;
    else
      e = m;
  }
  return s;
}


static int compar_double (const void * p1, const void * p2)
{
  const double * a = p1, * b = p2;
  return *a > *b;
}








int tag (scalar t)
{ trace ("tag", "/home/jiarongw/basilisk/src/tag.h", 63);




  _attribute[t.i].restriction = restriction_tag;

  _attribute[t.i].refine = _attribute[t.i].prolongation = refine_injection;
#line 81 "/home/jiarongw/basilisk/src/tag.h"
  scalar index= new_scalar("index");
  z_indexing (index, true);
   { foreach(){

#line 83 "/home/jiarongw/basilisk/src/tag.h"

    val(t,0,0,0) = (val(t,0,0,0) != 0)*(val(index,0,0,0) + 1); } end_foreach(); }
#line 93 "/home/jiarongw/basilisk/src/tag.h"
  boundary (((scalar []){t,{-1}}));





  bool changed;
  do {






    restriction (((scalar []){t,{-1}}));





    changed = false;
    for (int l = 1; l <= grid->maxdepth; l++) {






       { foreach_level(l){

#line 121 "/home/jiarongw/basilisk/src/tag.h"

 if (coarse(t,0,0,0))
   val(t,0,0,0) = coarse(t,0,0,0); } end_foreach_level(); }
      boundary_level (((scalar []){t,{-1}}), l);







       { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _changed = changed; 
#line 132
foreach_level (l){

#line 132 "/home/jiarongw/basilisk/src/tag.h"

        if (val(t,0,0,0)) {
   double min = val(t,0,0,0);
    { foreach_neighbor(1)
     if (val(t,0,0,0) && val(t,0,0,0) < min)
       min = val(t,0,0,0); end_foreach_neighbor(); }






   {
#line 144

     for (int i = -1; i <= 2; i += 3)
       if ((!is_leaf (neighbor((2*i - 1)/3,0,)) && neighbor((2*i - 1)/3,0,).neighbors && neighbor((2*i - 1)/3,0,).pid >= 0))
  for (int j = 0; j <= 1; j++)
    for (int k = 0; k <= 1; k++)
      if (fine(t,i,j,k) && fine(t,i,j,k) < min)
        min = fine(t,i,j,k);
#line 144

     for (int i = -1; i <= 2; i += 3)
       if ((!is_leaf (neighbor(0,(2*i - 1)/3,)) && neighbor(0,(2*i - 1)/3,).neighbors && neighbor(0,(2*i - 1)/3,).pid >= 0))
  for (int j = 0; j <= 1; j++)
    for (int k = 0; k <= 1; k++)
      if (fine(t,j,i,k) && fine(t,j,i,k) < min)
        min = fine(t,j,i,k);}


   if (val(t,0,0,0) != min) {
     _changed = true;
     val(t,0,0,0) = min;
   }
 } } end_foreach_level();OMP(omp critical) if (_changed > changed) changed = _changed;
mpi_all_reduce_double (changed, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 157
 }
      boundary_level (((scalar []){t,{-1}}), l);
    }
  } while (changed);
#line 171 "/home/jiarongw/basilisk/src/tag.h"
  Array * a = array_new();
   { foreach_leaf(){

#line 172 "/home/jiarongw/basilisk/src/tag.h"

    if (val(t,0,0,0) > 0) {







      double tag = val(t,0,0,0), * ap = (double *) a->p;
      long s = -1;
      if (a->len == 0 || tag > ap[a->len/sizeof(double) - 1])
 s = a->len/sizeof(double);
      else if (tag < ap[0])
 s = 0;
      else {






 s = lookup_tag (a, tag) + 1;
 if (tag == ap[s - 1] || tag == ap[s])
   s = -1;
      }
      if (s >= 0) {





 array_append (a, &tag, sizeof(double)), ap = (double *) a->p;
 for (int i = a->len/sizeof(double) - 1; i > s; i--)
   ap[i] = ap[i-1];
 ap[s] = tag;
      }
    } } end_foreach_leaf(); }
#line 222 "/home/jiarongw/basilisk/src/tag.h"
  long lmax = a->len;
  mpi_all_reduce (lmax, MPI_LONG, MPI_MAX);
  a->p = prealloc (a->p, lmax,__func__,__FILE__,__LINE__);
  lmax /= sizeof(double);






  double * q = a->p;
  for (int i = a->len/sizeof(double); i < lmax; i++)
    q[i] = -1;
  double p[lmax*npe()];
  MPI_Allgather (a->p, lmax, MPI_DOUBLE, p, lmax, MPI_DOUBLE, MPI_COMM_WORLD);
  qsort (p, lmax*npe(), sizeof(double), compar_double);







  array_free (a);
  a = array_new();
  double last = -1;
  for (int i = 0; i < lmax*npe(); i++)
    if (p[i] != last) {
      array_append (a, &p[i], sizeof(double));
      last = p[i];
    }






   { foreach(){

#line 259 "/home/jiarongw/basilisk/src/tag.h"

    if (val(t,0,0,0) > 0)
      val(t,0,0,0) = lookup_tag (a, val(t,0,0,0)) + 1; } end_foreach(); }
  boundary (((scalar []){t,{-1}}));




  int n = a->len/sizeof(double);
  array_free (a);
  { int _ret =  n; delete (((scalar []){index,{-1}}));  end_trace("tag", "/home/jiarongw/basilisk/src/tag.h", 269);  return _ret; }
 delete (((scalar []){index,{-1}}));  end_trace("tag", "/home/jiarongw/basilisk/src/tag.h", 270); }
#line 16 "wavewind.c"
#line 28 "wavewind.c"
double ak = 0.05;
double BO = 200.;
double RE = 40000.;




int LEVEL = 2 == 2 ? 10 : 6;




double uemax = 0.001;
double femax = 0.00001;




double RATIO = (1./850.);
double MURATIO = (17.4e-6/8.9e-4);



int DIRAC = 0;



int countvelocity = 0;




double k_ = 2.*pi;
double h_ = 0.5;
double g_ = 1.;



double m = 5.;
double B = 0.;
double Karman = 0.41;
double UstarRATIO = 1;







int main (int argc, char * argv[])
{ _init_solver();
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  if (argc > 2)
    ak = atof(argv[2]);
  if (argc > 3)
    BO = atof(argv[3]);
  if (argc > 4)
    RE = atof(argv[4]);
  if (argc > 5)
    m = atof(argv[5]);
  if (argc > 6)
    B = atof(argv[6]);
  if (argc > 7)
    UstarRATIO = atof(argv[7]);
  if (argc > 8)
    DIRAC = atof(argv[8]);






  origin ((struct _origin){-L0/2, -L0/2, -L0/2});
  tree_periodic(right);
  _attribute[u.x.i].boundary[top] = _boundary4; _attribute[u.x.i].boundary_homogeneous[top] = _boundary4_homogeneous;
  _attribute[u.y.i].boundary[top] = _boundary5; _attribute[u.y.i].boundary_homogeneous[top] = _boundary5_homogeneous;
#line 113 "wavewind.c"
  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0/RE;
  mu2 = 1.0/RE*MURATIO;
  _attribute[f.i].sigma = 1./(BO*sq(k_));
  G.y = -g_;






  N = 32;



  run();
 free_solver(); }
#line 139 "wavewind.c"
double wave (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa -
   3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) +
    3./64.*(8.*cube(alpa)*cube(alpa) +
     (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3 - y;
}

double eta (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa -
   3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) +
    3./64.*(8.*cube(alpa)*cube(alpa) +
     (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3;
}






double gaus (double y, double yc, double T){
  double deltaw = sqrt(2.0/RE)/k_;
  double deltaa = sqrt(2.0/RE*MURATIO/RATIO)/k_;
  double r = y - yc;
  return 2.0/(sqrt(2.0*pi*sq(deltaa)) + sqrt(2.0*pi*sq(deltaw))) *
    (T*exp(-sq(r)/(2.0*sq(deltaw))) + (1.0 - T)*exp(-sq(r)/(2.0*sq(deltaa))));
}





static int init_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init_0 (const int i, const double t, Event * _ev) { trace ("init_0", "wavewind.c", 180); 
{


  double Ustar = sqrt(g_/k_)*UstarRATIO;
  double y1 = m*mu2/rho2/Ustar;
  double Udrift = B*Ustar;
  fprintf(qstderr(), "UstarRATIO=%g B=%g\n m=%g ak=%g", UstarRATIO, B, m, ak);

  if (!restore ((struct Dump){"restart"})) {
    do {
      do { scalar phi= new_vertex_scalar("phi");  { foreach_vertex(){

#line 191 "wavewind.c"
 val(phi,0,0,0) = wave(x,y); } end_foreach_vertex(); } fractions ((struct Fractions){phi, f});  delete (((scalar []){phi,{-1}})); } while(0);




      scalar Phi= new_scalar("Phi");
       { foreach(){

#line 197 "wavewind.c"
 {
       double alpa = 1./tanh(k_*h_);
       double a_ = ak/k_;
       double sgma = sqrt(g_*k_*tanh(k_*h_)*
            (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
                 (sq(alpa) - 1.) + sq(alpa))));
       double A_ = a_*g_/sgma;
       double phi1 = A_*cosh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x);
       double phi2 = 3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
         cosh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_);
       double phi3 = 1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
         (9.*sq(alpa) - 13.)*cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*k_*x);
       val(Phi,0,0,0) = phi1 + ak*phi2 + ak*ak*phi3;
      } } end_foreach(); }
      boundary (((scalar []){Phi,{-1}}));
      if (DIRAC){



 scalar vort2= new_scalar("vort2");
 scalar psi= new_scalar("psi");
  { foreach(){

#line 218 "wavewind.c"
 {
   val(vort2,0,0,0) = -2.0*gaus(y,wave(x,y)+y,val(f,0,0,0))*(val(Phi,1,0,0)-val(Phi,-1,0,0))/(2.*Delta);
   val(psi,0,0,0) = 0.0;
 } } end_foreach(); }

 boundary(((scalar []){vort2,psi,{-1}}));
 _attribute[psi.i].boundary[top] = _boundary6; _attribute[psi.i].boundary_homogeneous[top] = _boundary6_homogeneous;
 _attribute[psi.i].boundary[bottom] = _boundary7; _attribute[psi.i].boundary_homogeneous[bottom] = _boundary7_homogeneous;



 poisson((struct Poisson){psi, vort2});




  { foreach(){

#line 234 "wavewind.c"
{
   val(u.x,0,0,0) = (val(psi,0,1,0) - val(psi,0,-1,0))/(2.*Delta);
   val(u.y,0,0,0) = -(val(psi,1,0,0) - val(psi,-1,0,0))/(2.*Delta);
 } } end_foreach(); }
       delete (((scalar []){psi,vort2,{-1}})); }
      else{



  { foreach(){

#line 243 "wavewind.c"
{
   {
#line 244

     val(u.x,0,0,0) = (val(Phi,1,0,0) - val(Phi,-1,0,0))/(2.0*Delta) * val(f,0,0,0);
#line 244

     val(u.y,0,0,0) = (val(Phi,0,1,0) - val(Phi,0,-1,0))/(2.0*Delta) * val(f,0,0,0);}
 } } end_foreach(); }
  { foreach(){

#line 247 "wavewind.c"
{
   if ((y-eta(x,y))<y1){
     val(u.x,0,0,0) += Udrift + sq(Ustar)/(mu2/rho2)* (y-eta(x,y)) * (1-val(f,0,0,0));
   }
   else{
     double beta = 2*Karman*Ustar/mu2*rho2*((y-eta(x,y))-y1);
     double alpha = log(beta+sqrt(sq(beta)+1));
     double tanhtemp = (exp(alpha/2)-exp(-alpha/2))/(exp(alpha/2)+exp(-alpha/2));
     val(u.x,0,0,0) += (Udrift + m*Ustar + Ustar/Karman*(alpha-tanhtemp)) * (1-val(f,0,0,0));
   }
 } } end_foreach(); }


      }
      boundary ((scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}));
     delete (((scalar []){Phi,{-1}})); }






    while (adapt_wavelet ((struct Adapt){((scalar []){f,u.x,u.y,{-1}}),
     (double[]){0.01,0.05,0.05,0.05}, LEVEL, 5}).nf);



  fprintf(qstderr(), "break3!");
  }
 end_trace("init_0", "wavewind.c", 276); } return 0; } 
#line 289 "wavewind.c"
int dissipation_rate (double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _rateWater = rateWater; double _rateAir = rateAir; 
#line 293

if (!is_constant(cm) && !is_constant(rho)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#line 293
foreach (){

#line 293 "wavewind.c"
 {
    double dudx = (val(u.x,1,0,0) - val(u.x,-1,0,0) )/(2.*Delta);
    double dudy = (val(u.x,0,1,0) - val(u.x,0,-1,0) )/(2.*Delta);
    double dudz = (val(u.x,0,0,1) - val(u.x,0,0,-1))/(2.*Delta);
    double dvdx = (val(u.y,1,0,0) - val(u.y,-1,0,0) )/(2.*Delta);
    double dvdy = (val(u.y,0,1,0) - val(u.y,0,-1,0) )/(2.*Delta);
    double dvdz = (val(u.y,0,0,1) - val(u.y,0,0,-1))/(2.*Delta);
    double dwdx = (_val_higher_dimension(u.z,1,0,0) - _val_higher_dimension(u.z,-1,0,0) )/(2.*Delta);
    double dwdy = (_val_higher_dimension(u.z,0,1,0) - _val_higher_dimension(u.z,0,-1,0) )/(2.*Delta);
    double dwdz = (_val_higher_dimension(u.z,0,0,1) - _val_higher_dimension(u.z,0,0,-1))/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz;
    double sqterm = 2.*(sq(Delta)*val_cm(cm,0,0,0))*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
        sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
        sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
    _rateWater += mu1/val_rho(rho,0,0,0)*val(f,0,0,0)*sqterm;
    _rateAir += mu2/val_rho(rho,0,0,0)*(1. - val(f,0,0,0))*sqterm;
  } } end_foreach(); }
if (is_constant(cm) && !is_constant(rho)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#line 293
foreach (){

#line 293 "wavewind.c"
 {
    double dudx = (val(u.x,1,0,0) - val(u.x,-1,0,0) )/(2.*Delta);
    double dudy = (val(u.x,0,1,0) - val(u.x,0,-1,0) )/(2.*Delta);
    double dudz = (val(u.x,0,0,1) - val(u.x,0,0,-1))/(2.*Delta);
    double dvdx = (val(u.y,1,0,0) - val(u.y,-1,0,0) )/(2.*Delta);
    double dvdy = (val(u.y,0,1,0) - val(u.y,0,-1,0) )/(2.*Delta);
    double dvdz = (val(u.y,0,0,1) - val(u.y,0,0,-1))/(2.*Delta);
    double dwdx = (_val_higher_dimension(u.z,1,0,0) - _val_higher_dimension(u.z,-1,0,0) )/(2.*Delta);
    double dwdy = (_val_higher_dimension(u.z,0,1,0) - _val_higher_dimension(u.z,0,-1,0) )/(2.*Delta);
    double dwdz = (_val_higher_dimension(u.z,0,0,1) - _val_higher_dimension(u.z,0,0,-1))/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz;
    double sqterm = 2.*(sq(Delta)*val_cm(cm,0,0,0))*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
        sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
        sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
    _rateWater += mu1/val_rho(rho,0,0,0)*val(f,0,0,0)*sqterm;
    _rateAir += mu2/val_rho(rho,0,0,0)*(1. - val(f,0,0,0))*sqterm;
  } } end_foreach(); }
if (!is_constant(cm) && is_constant(rho)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 293
foreach (){

#line 293 "wavewind.c"
 {
    double dudx = (val(u.x,1,0,0) - val(u.x,-1,0,0) )/(2.*Delta);
    double dudy = (val(u.x,0,1,0) - val(u.x,0,-1,0) )/(2.*Delta);
    double dudz = (val(u.x,0,0,1) - val(u.x,0,0,-1))/(2.*Delta);
    double dvdx = (val(u.y,1,0,0) - val(u.y,-1,0,0) )/(2.*Delta);
    double dvdy = (val(u.y,0,1,0) - val(u.y,0,-1,0) )/(2.*Delta);
    double dvdz = (val(u.y,0,0,1) - val(u.y,0,0,-1))/(2.*Delta);
    double dwdx = (_val_higher_dimension(u.z,1,0,0) - _val_higher_dimension(u.z,-1,0,0) )/(2.*Delta);
    double dwdy = (_val_higher_dimension(u.z,0,1,0) - _val_higher_dimension(u.z,0,-1,0) )/(2.*Delta);
    double dwdz = (_val_higher_dimension(u.z,0,0,1) - _val_higher_dimension(u.z,0,0,-1))/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz;
    double sqterm = 2.*(sq(Delta)*val_cm(cm,0,0,0))*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
        sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
        sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
    _rateWater += mu1/val_rho(rho,0,0,0)*val(f,0,0,0)*sqterm;
    _rateAir += mu2/val_rho(rho,0,0,0)*(1. - val(f,0,0,0))*sqterm;
  } } end_foreach(); }
if (is_constant(cm) && is_constant(rho)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 293
foreach (){

#line 293 "wavewind.c"
 {
    double dudx = (val(u.x,1,0,0) - val(u.x,-1,0,0) )/(2.*Delta);
    double dudy = (val(u.x,0,1,0) - val(u.x,0,-1,0) )/(2.*Delta);
    double dudz = (val(u.x,0,0,1) - val(u.x,0,0,-1))/(2.*Delta);
    double dvdx = (val(u.y,1,0,0) - val(u.y,-1,0,0) )/(2.*Delta);
    double dvdy = (val(u.y,0,1,0) - val(u.y,0,-1,0) )/(2.*Delta);
    double dvdz = (val(u.y,0,0,1) - val(u.y,0,0,-1))/(2.*Delta);
    double dwdx = (_val_higher_dimension(u.z,1,0,0) - _val_higher_dimension(u.z,-1,0,0) )/(2.*Delta);
    double dwdy = (_val_higher_dimension(u.z,0,1,0) - _val_higher_dimension(u.z,0,-1,0) )/(2.*Delta);
    double dwdz = (_val_higher_dimension(u.z,0,0,1) - _val_higher_dimension(u.z,0,0,-1))/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz;
    double sqterm = 2.*(sq(Delta)*val_cm(cm,0,0,0))*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
        sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
        sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
    _rateWater += mu1/val_rho(rho,0,0,0)*val(f,0,0,0)*sqterm;
    _rateAir += mu2/val_rho(rho,0,0,0)*(1. - val(f,0,0,0))*sqterm;
  } } end_foreach(); }OMP(omp critical) rateWater += _rateWater;
mpi_all_reduce_double (rateWater, MPI_SUM);
OMP(omp critical) rateAir += _rateAir;
mpi_all_reduce_double (rateAir, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 317
 }
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}






static int graphs_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int graphs (const int i, const double t, Event * _ev) { trace ("graphs", "wavewind.c", 328);  {
  static FILE * fpwater =NULL; if (!fpwater || i == 0) fpwater = pid() > 0 ? fopen("/dev/null", "w") :  fopen("budgetWaterwind.dat", "w");
  static FILE * fpair =NULL; if (!fpair || i == 0) fpair = pid() > 0 ? fopen("/dev/null", "w") :  fopen("budgetAirwind.dat", "w");
  double ke = 0., gpe = 0.;
  double keAir = 0., gpeAir = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _ke = ke; double _gpe = gpe; double _keAir = keAir; double _gpeAir = gpeAir; 
#line 333

if (!is_constant(rho) && !is_constant(cm)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 333
foreach(){

#line 334 "wavewind.c"
 {
    double norm2 = 0.;
    {
#line 336

      norm2 += sq(val(u.x,0,0,0));
#line 336

      norm2 += sq(val(u.y,0,0,0));}
    _ke += val_rho(rho,0,0,0)*norm2*val(f,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
    _keAir += val_rho(rho,0,0,0)*norm2*(1.0-val(f,0,0,0))*(sq(Delta)*val_cm(cm,0,0,0));
    _gpe += rho1*g_*y*val(f,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
    _gpeAir += rho2*g_*y*(1.0-val(f,0,0,0))*(sq(Delta)*val_cm(cm,0,0,0));
  } } end_foreach(); }
if (is_constant(rho) && !is_constant(cm)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 333
foreach(){

#line 334 "wavewind.c"
 {
    double norm2 = 0.;
    {
#line 336

      norm2 += sq(val(u.x,0,0,0));
#line 336

      norm2 += sq(val(u.y,0,0,0));}
    _ke += val_rho(rho,0,0,0)*norm2*val(f,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
    _keAir += val_rho(rho,0,0,0)*norm2*(1.0-val(f,0,0,0))*(sq(Delta)*val_cm(cm,0,0,0));
    _gpe += rho1*g_*y*val(f,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
    _gpeAir += rho2*g_*y*(1.0-val(f,0,0,0))*(sq(Delta)*val_cm(cm,0,0,0));
  } } end_foreach(); }
if (!is_constant(rho) && is_constant(cm)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 333
foreach(){

#line 334 "wavewind.c"
 {
    double norm2 = 0.;
    {
#line 336

      norm2 += sq(val(u.x,0,0,0));
#line 336

      norm2 += sq(val(u.y,0,0,0));}
    _ke += val_rho(rho,0,0,0)*norm2*val(f,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
    _keAir += val_rho(rho,0,0,0)*norm2*(1.0-val(f,0,0,0))*(sq(Delta)*val_cm(cm,0,0,0));
    _gpe += rho1*g_*y*val(f,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
    _gpeAir += rho2*g_*y*(1.0-val(f,0,0,0))*(sq(Delta)*val_cm(cm,0,0,0));
  } } end_foreach(); }
if (is_constant(rho) && is_constant(cm)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 333
foreach(){

#line 334 "wavewind.c"
 {
    double norm2 = 0.;
    {
#line 336

      norm2 += sq(val(u.x,0,0,0));
#line 336

      norm2 += sq(val(u.y,0,0,0));}
    _ke += val_rho(rho,0,0,0)*norm2*val(f,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
    _keAir += val_rho(rho,0,0,0)*norm2*(1.0-val(f,0,0,0))*(sq(Delta)*val_cm(cm,0,0,0));
    _gpe += rho1*g_*y*val(f,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
    _gpeAir += rho2*g_*y*(1.0-val(f,0,0,0))*(sq(Delta)*val_cm(cm,0,0,0));
  } } end_foreach(); }OMP(omp critical) ke += _ke;
mpi_all_reduce_double (ke, MPI_SUM);
OMP(omp critical) gpe += _gpe;
mpi_all_reduce_double (gpe, MPI_SUM);
OMP(omp critical) keAir += _keAir;
mpi_all_reduce_double (keAir, MPI_SUM);
OMP(omp critical) gpeAir += _gpeAir;
mpi_all_reduce_double (gpeAir, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 342
 }
  double rates[2];
  dissipation_rate(rates);
  double dissWater = rates[0];
  double dissAir = rates[1];
  if (i == 0) {
    fprintf (fpwater, "t ke gpe dissipation\n");
    fprintf (fpair, "t ke gpe dissipation\n");
  }
  fprintf (fpwater, "%g %g %g %g\n",
    t/(k_/sqrt(g_*k_)), ke/2., gpe + 0.125, dissWater);
  fprintf (fpair, "%g %g %g %g\n",
    t/(k_/sqrt(g_*k_)), keAir/2., gpeAir + 0.125, dissAir);
  fprintf (ferr, "%g %g %g %g\n",
    t/(k_/sqrt(g_*k_)), ke/2., gpe + 0.125, dissWater);
 end_trace("graphs", "wavewind.c", 357); } return 0; } 
#line 372 "wavewind.c"
static int movies_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t += 0.1);   *ip = i; *tp = t;   return ret; } static int movies (const int i, const double t, Event * _ev) { trace ("movies", "wavewind.c", 372);  {





  {
    static FILE * fp =NULL; if (!fp || i == 0) fp = pid() > 0 ? fopen("/dev/null", "w") :  fopen ("f" ".ppm", "w");
    output_ppm ((struct OutputPPM){f, fp, .min = 0, .max = 1, .n = 512});
  }


  {
    scalar l= new_scalar("l");
     { foreach(){

#line 386 "wavewind.c"

      val(l,0,0,0) = level; } end_foreach(); }
    static FILE * fp =NULL; if (!fp || i == 0) fp = pid() > 0 ? fopen("/dev/null", "w") :  fopen ("level" ".ppm", "w");
    output_ppm ((struct OutputPPM){l, fp, .min = 5, .max = LEVEL, .n = 512});
   delete (((scalar []){l,{-1}})); }
#line 404 "wavewind.c"
  scalar omega= new_scalar("omega");
  vorticity (u, omega);

  view ((struct _view_set){.width = 800, .height = 600, .fov = 18.8});
  clear();




  for (double x = -L0; x <= L0; x += L0)
    { begin_translate ((struct _translate){x}); {
      draw_vof ((struct _draw_vof){"f"});
      squares ((struct _squares){"omega", .linear = true});
    } end_translate(); }
#line 462 "wavewind.c"
  {
    static FILE * fp =NULL; if (!fp || i == 0) fp = pid() > 0 ? fopen("/dev/null", "w") :  fopen ("movie" ".ppm", "w");
    save ((struct _save){.fp = fp});
  }
 delete (((scalar []){omega,{-1}}));  end_trace("movies", "wavewind.c", 466); } return 0; } 







static int snapshot_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i += 200);   *ip = i; *tp = t;   return ret; } static int snapshot (const int i, const double t, Event * _ev) { trace ("snapshot", "wavewind.c", 474);  {
  dump ((struct Dump){"dump"});
 end_trace("snapshot", "wavewind.c", 476); } return 0; } 







static int end_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = 2.*k_/sqrt(g_*k_));   *ip = i; *tp = t;   return ret; } static int end (const int i, const double t, Event * _ev) { trace ("end", "wavewind.c", 484);  {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ((struct Dump){"end"});
 end_trace("end", "wavewind.c", 487); } return 0; } 

static int dumpstep_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t += k_/sqrt(g_*k_)/32);   *ip = i; *tp = t;   return ret; } static int dumpstep (const int i, const double t, Event * _ev) { trace ("dumpstep", "wavewind.c", 489);  {
  char dname[100];
  sprintf (dname, "dump%g", t/(k_/sqrt(g_*k_)));
  dump ((struct Dump){dname});
 end_trace("dumpstep", "wavewind.c", 493); } return 0; } 
#line 502 "wavewind.c"
static int adapt_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int adapt_0 (const int i, const double t, Event * _ev) { trace ("adapt_0", "wavewind.c", 502);  {
  adapt_wavelet ((struct Adapt){((scalar []){f,u.x,u.y,{-1}}), (double[]){femax,uemax,uemax,uemax}, LEVEL, 5});
 end_trace("adapt_0", "wavewind.c", 504); } return 0; } 
#line 79 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
static double _boundary0 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 78 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); } return 0.; } static double _boundary0_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 78 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 79
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 79
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 79
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 79
return  neumann_homogeneous(); } return 0.; }
#line 80 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
static double _boundary1 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 79 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); } return 0.; } static double _boundary1_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 79 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 80
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 80
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 80
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 80
return  neumann_homogeneous(); } return 0.; }
#line 88 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
static double _boundary2 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 87 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)); } return 0.; } static double _boundary2_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 87 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 88
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 88
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 88
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 88
return  neumann_homogeneous(); } return 0.; }
#line 89 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"
static double _boundary3 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 88 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); } return 0.; } static double _boundary3_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 88 "/home/jiarongw/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 89
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 89
return  neumann_homogeneous(); } return 0.; }
#line 103 "wavewind.c"
static double _boundary4 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 102 "wavewind.c"
return  dirichlet(0); return 0.; } static double _boundary4_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 102 "wavewind.c"
return  dirichlet_homogeneous(); return 0.; }
#line 104 "wavewind.c"
static double _boundary5 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 103 "wavewind.c"
return  neumann(0); return 0.; } static double _boundary5_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 103 "wavewind.c"
return  neumann_homogeneous(); return 0.; }
#line 224 "wavewind.c"
static double _boundary6 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 223 "wavewind.c"
return  dirichlet(0.); return 0.; } static double _boundary6_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 223 "wavewind.c"
return  dirichlet_homogeneous(); return 0.; }
#line 225 "wavewind.c"
static double _boundary7 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 224 "wavewind.c"
return  dirichlet(0.); return 0.; } static double _boundary7_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 224 "wavewind.c"
return  dirichlet_homogeneous(); return 0.; }
size_t datasize = 12*sizeof (double);
static int defaults (const int i, const double t, Event * _ev);
static int defaults_expr0 (int * ip, double * tp, Event * _ev);
static int init (const int i, const double t, Event * _ev);
static int init_expr0 (int * ip, double * tp, Event * _ev);
static int set_dtmax (const int i, const double t, Event * _ev);
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev);
static int stability (const int i, const double t, Event * _ev);
static int stability_expr0 (int * ip, double * tp, Event * _ev);
static int vof (const int i, const double t, Event * _ev);
static int vof_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_advection (const int i, const double t, Event * _ev);
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_diffusion (const int i, const double t, Event * _ev);
static int tracer_diffusion_expr0 (int * ip, double * tp, Event * _ev);
static int properties (const int i, const double t, Event * _ev);
static int properties_expr0 (int * ip, double * tp, Event * _ev);
static int advection_term (const int i, const double t, Event * _ev);
static int advection_term_expr0 (int * ip, double * tp, Event * _ev);
static int viscous_term (const int i, const double t, Event * _ev);
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration (const int i, const double t, Event * _ev);
static int acceleration_expr0 (int * ip, double * tp, Event * _ev);
static int projection (const int i, const double t, Event * _ev);
static int projection_expr0 (int * ip, double * tp, Event * _ev);
static int adapt (const int i, const double t, Event * _ev);
static int adapt_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_0 (const int i, const double t, Event * _ev);
static int defaults_0_expr0 (int * ip, double * tp, Event * _ev);
static int stability_0 (const int i, const double t, Event * _ev);
static int stability_0_expr0 (int * ip, double * tp, Event * _ev);
static int vof_0 (const int i, const double t, Event * _ev);
static int vof_0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_1 (const int i, const double t, Event * _ev);
static int defaults_1_expr0 (int * ip, double * tp, Event * _ev);
static int properties_0 (const int i, const double t, Event * _ev);
static int properties_0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_2 (const int i, const double t, Event * _ev);
static int defaults_2_expr0 (int * ip, double * tp, Event * _ev);
static int vof_1 (const int i, const double t, Event * _ev);
static int vof_1_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_advection_0 (const int i, const double t, Event * _ev);
static int tracer_advection_0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_3 (const int i, const double t, Event * _ev);
static int defaults_3_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_0 (const int i, const double t, Event * _ev);
static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev);
static int stability_1 (const int i, const double t, Event * _ev);
static int stability_1_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_1 (const int i, const double t, Event * _ev);
static int acceleration_1_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_2 (const int i, const double t, Event * _ev);
static int acceleration_2_expr0 (int * ip, double * tp, Event * _ev);
static int init_0 (const int i, const double t, Event * _ev);
static int init_0_expr0 (int * ip, double * tp, Event * _ev);
static int graphs (const int i, const double t, Event * _ev);
static int graphs_expr0 (int * ip, double * tp, Event * _ev);
static int movies (const int i, const double t, Event * _ev);
static int movies_expr0 (int * ip, double * tp, Event * _ev);
static int snapshot (const int i, const double t, Event * _ev);
static int snapshot_expr0 (int * ip, double * tp, Event * _ev);
static int end (const int i, const double t, Event * _ev);
static int end_expr0 (int * ip, double * tp, Event * _ev);
static int dumpstep (const int i, const double t, Event * _ev);
static int dumpstep_expr0 (int * ip, double * tp, Event * _ev);
static int adapt_0 (const int i, const double t, Event * _ev);
static int adapt_0_expr0 (int * ip, double * tp, Event * _ev);
static void _set_boundary0 (void);
static void _set_boundary1 (void);
static void _set_boundary2 (void);
static void _set_boundary3 (void);
void _init_solver (void) {
  void init_solver();
  init_solver();
  Events = (Event *) pmalloc (sizeof (Event), __func__, __FILE__, __LINE__);
  Events[0].last = 1;
  event_register ((Event){ 0, 1, defaults, {defaults_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 100, "defaults"});
  event_register ((Event){ 0, 1, defaults_0, {defaults_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/vof.h", 51, "defaults"});
  event_register ((Event){ 0, 1, defaults_1, {defaults_1_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/two-phase.h", 25, "defaults"});
  event_register ((Event){ 0, 1, defaults_2, {defaults_2_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/conserving.h", 37, "defaults"});
  event_register ((Event){ 0, 1, defaults_3, {defaults_3_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/iforce.h", 30, "defaults"});
  event_register ((Event){ 0, 1, init, {init_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 139, "init"});
  event_register ((Event){ 0, 1, init_0, {init_0_expr0}, ((int *)0), ((double *)0),
    "wavewind.c", 180, "init"});
  event_register ((Event){ 0, 1, graphs, {graphs_expr0}, ((int *)0), ((double *)0),
    "wavewind.c", 328, "graphs"});
  event_register ((Event){ 0, 1, movies, {movies_expr0}, ((int *)0), ((double *)0),
    "wavewind.c", 372, "movies"});
  event_register ((Event){ 0, 1, snapshot, {snapshot_expr0}, ((int *)0), ((double *)0),
    "wavewind.c", 474, "snapshot"});
  event_register ((Event){ 0, 1, end, {end_expr0}, ((int *)0), ((double *)0),
    "wavewind.c", 484, "end"});
  event_register ((Event){ 0, 1, dumpstep, {dumpstep_expr0}, ((int *)0), ((double *)0),
    "wavewind.c", 489, "dumpstep"});
  event_register ((Event){ 0, 1, set_dtmax, {set_dtmax_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 167, "set_dtmax"});
  event_register ((Event){ 0, 1, stability, {stability_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 169, "stability"});
  event_register ((Event){ 0, 1, stability_0, {stability_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/vof.h", 63, "stability"});
  event_register ((Event){ 0, 1, stability_1, {stability_1_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/tension.h", 36, "stability"});
  event_register ((Event){ 0, 1, vof, {vof_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 179, "vof"});
  event_register ((Event){ 0, 1, vof_0, {vof_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/vof.h", 313, "vof"});
  event_register ((Event){ 0, 1, vof_1, {vof_1_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/conserving.h", 107, "vof"});
  event_register ((Event){ 0, 1, tracer_advection, {tracer_advection_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 180, "tracer_advection"});
  event_register ((Event){ 0, 1, tracer_advection_0, {tracer_advection_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/conserving.h", 178, "tracer_advection"});
  event_register ((Event){ 0, 1, tracer_diffusion, {tracer_diffusion_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 181, "tracer_diffusion"});
  event_register ((Event){ 0, 1, properties, {properties_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 188, "properties"});
  event_register ((Event){ 0, 1, properties_0, {properties_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/two-phase.h", 59, "properties"});
  event_register ((Event){ 0, 1, advection_term, {advection_term_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 249, "advection_term"});
  event_register ((Event){ 0, 1, viscous_term, {viscous_term_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 279, "viscous_term"});
  event_register ((Event){ 0, 1, acceleration, {acceleration_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 314, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_0, {acceleration_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/iforce.h", 44, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_1, {acceleration_1_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/tension.h", 70, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_2, {acceleration_2_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/reduced.h", 36, "acceleration"});
  event_register ((Event){ 0, 1, projection, {projection_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 329, "projection"});
  event_register ((Event){ 0, 1, adapt, {adapt_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/navier-stokes/centered.h", 364, "adapt"});
  event_register ((Event){ 0, 1, adapt_0, {adapt_0_expr0}, ((int *)0), ((double *)0),
    "wavewind.c", 502, "adapt"});
  _attribute = (_Attributes *) pcalloc (datasize/sizeof(double), sizeof (_Attributes), __func__, __FILE__, __LINE__);
  all = (scalar *) pmalloc (sizeof (scalar)*13,__func__, __FILE__, __LINE__);
  for (int i = 0; i < 12; i++)
    all[i].i = i;
  all[12].i = -1;
  set_fpe();
  quadtree_methods();
  init_scalar ((scalar){11}, "rhov");
  init_face_vector ((vector){{9},{10}}, "alphav");
  init_scalar ((scalar){8}, "f");
  init_face_vector ((vector){{6},{7}}, "uf");
  init_scalar ((scalar){5}, "pf");
  init_vector ((vector){{3},{4}}, "g");
  init_vector ((vector){{1},{2}}, "u");
  init_scalar ((scalar){0}, "p");
  init_const_scalar ((scalar){_NVARMAX+5}, "zeroc",  0.);
  init_const_scalar ((scalar){_NVARMAX+4}, "unity",  1.);
  init_const_vector ((vector){{_NVARMAX+2},{_NVARMAX+3}}, "unityf", (double []) {1.,1.,1.});
  init_const_vector ((vector){{_NVARMAX+0},{_NVARMAX+1}}, "zerof", (double []) {0.,0.,0.});
  _set_boundary0();
  _set_boundary1();
  _set_boundary2();
  _set_boundary3();
}
