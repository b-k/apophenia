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
    double *xinit;  /**< A <tt>double*</tt> giving starting values for x in ascending order. Default: -1, 0, 1. If this isn't \c NULL, I need at least three items. */
    double  xl;     /**< Left bound. If you don't give me one, I'll use min[min(xinit)/10, min(xinit)*10].*/
    double  xr;     /**< Right bound. If you don't give me one, I'll use max[max(xinit)/10, max(xinit)*10]. */
    double convex;  /**< Adjustment for convexity */
    int ninit;      /**< Number of starting values supplied (i.e. number of elements in \c xinit)*/
    int npoint;     /**< Maximum number of envelope points. I \c malloc space for this many <tt>double</tt>s at the outset. Default = 1e5. */
   char do_metro;   /**< Whether metropolis step is required. (I.e., set to one if you're not sure if the function is log-concave). Set  to <tt>'y'</tt>es or <tt>'n'</tt>o*/
   double xprev;    /**< Previous value from Markov chain */
   int neval;       /**< On exit, the number of function evaluations performed */
   arms_state *state;
   apop_model *model; /**< The model from which I will draw. Mandatory. Must have either a \c log_likelihood or \c p method.*/
} apop_arms_settings;

Apop_settings_declarations(apop_arms)


void apop_arms_draw (double *out, gsl_rng *r, apop_model *m);

#define YCEIL 50.                /* maximum y avoiding overflow in exp(y) */

#ifdef	__cplusplus
}
#endif
#endif
