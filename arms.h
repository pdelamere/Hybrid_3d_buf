/* header file for arms function */

#ifndef ARMS
#define ARMS
typedef struct point {    /* a point in the x,y plane */
  double x,y;             /* x and y coordinates */
  double ey;              /* exp(y-ymax+YCEIL) */
  double cum;             /* integral up to x of rejection envelope */
  int f;                  /* is y an evaluated point of log-density */
  struct point *pl,*pr;   /* envelope points to left and right of x */
} POINT;

/* *********************************************************************** */

typedef struct envelope {  /* attributes of the entire rejection envelope */
  int cpoint;              /* number of POINTs in current envelope */
  int npoint;              /* max number of POINTs allowed in envelope */
  int *neval;              /* number of function evaluations performed */
  double ymax;             /* the maximum y-value in the current envelope */
  POINT *p;                /* start of storage of envelope POINTs */
  double *convex;          /* adjustment for convexity */
} ENVELOPE;

/* *********************************************************************** */

typedef struct funbag { /* everything for evaluating log density          */
  void *mydata;      /* user-defined structure holding data for density */
  double (*myfunc)(double x, void *mydata);
                     /* user-defined function evaluating log density at x */
} FUNBAG;

/* *********************************************************************** */

typedef struct metropolis { /* for metropolis step */
  int on;            /* whether metropolis is to be used */
  double xprev;      /* previous Markov chain iterate */
  double yprev;      /* current log density at xprev */
} METROPOLIS;

/* *********************************************************************** */

/* Struct to contain the state of the arms algorithm */
typedef struct arms_state {
    ENVELOPE *env;
    POINT pwork;
    FUNBAG lpdf;
    METROPOLIS *metrop;
} ARMS_STATE;

/* *********************************************************************** */
int arms_setup (double *xinit, int ninit, double *xl, double *xr,
                double (*myfunc)(double x, void *mydata), void *mydata,
                double *convex, int npoint, int dometrop, double *xprev,
                int *neval, ARMS_STATE **state_ref);

int arms_simple_setup (int ninit, double *xl, double *xr, 
                 double (*myfunc)(double x, void *mydata), void *mydata,
                 int dometrop, double *xprev, ARMS_STATE **state_ref);

int arms_sample(double *xsamp, int nsamp, ARMS_STATE *state);

void arms_cleanup(ARMS_STATE *state);

int arms (double *xinit, int ninit, double *xl, double *xr, 
          double (*myfunc)(double x, void *mydata), void *mydata,
          double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
          int nsamp, double *qcent, double *xcent, int ncent,
          int *neval);


int glob_neval;


#define YCEIL 50.                /* maximum y avoiding overflow in exp(y) */

#endif
