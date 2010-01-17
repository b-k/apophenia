/** \file arms.h  */
/* Copyright (c) 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#ifndef apop_arms_h
#define apop_arms_h

#ifdef	__cplusplus
extern "C" {
#endif

#include "model.h"
#include "settings.h"

    /** \cond doxy_ignore */
typedef struct point {    /* a point in the x,y plane */
  double x,y;             /* x and y coordinates */
  double ey;              /* exp(y-ymax+YCEIL) */
  double cum;             /* integral up to x of rejection envelope */
  int f;                  /* is y an evaluated point of log-density */
  struct point *pl,*pr;   /* envelope points to left and right of x */
} POINT;

/* This includes the envelope info and the metropolis steps. */
typedef struct {  /* attributes of the entire rejection envelope */
  int cpoint;              /* number of POINTs in current envelope */
  int npoint;              /* max number of POINTs allowed in envelope */
  double ymax;             /* the maximum y-value in the current envelope */
  POINT *p;                /* start of storage of envelope POINTs */
  double *convex;          /* adjustment for convexity */
  double metro_xprev;      /* previous Markov chain iterate */
  double metro_yprev;      /* current log density at xprev */
} arms_state;
    /** \endcond */

/** to perform derivative-free adaptive rejection sampling with metropolis step */
typedef struct {
    double *xinit;  /**< Starting values for x in ascending order */
    double  xl;     /**< left bound */
    double  xr;     /**< right bound */
    double convex;  /**< adjustment for convexity */
    int ninit;      /**< number of starting values supplied */
    int npoint;     /**< maximum number of envelope points */
   char do_metro;   /**< whether metropolis step is required */
   double xprev;    /**< previous value from markov chain */
   int neval;       /**< on exit, the number of function evaluations performed */
   arms_state *state;
   apop_model *model;
} apop_arms_settings;

Apop_settings_declarations(apop_arms)


void apop_arms_draw (double *out, gsl_rng *r, apop_model *m);

#define YCEIL 50.                /* maximum y avoiding overflow in exp(y) */

#ifdef	__cplusplus
}
#endif
#endif
