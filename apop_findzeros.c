#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <apophenia/headers.h>

/** \file apop_findzeros.c
 This just includes the root-finding routine. It is \c \#included in apop_mle.c, because I expect you to call it via that. 

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */


/* This function is cut/pasted/modified from the GSL documentation. It
 calls the various GSL root-finding algorithms to find the zero of the score.
*/
static apop_model * find_roots (infostruct p) {
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  apop_model *dist = p.model;
  int           vsize       =(dist->parameters->vector ? dist->parameters->vector->size :0),
                msize1      =(dist->parameters->matrix ? dist->parameters->matrix->size1 :0),
                msize2      =(dist->parameters->matrix ? dist->parameters->matrix->size2:0);
  int status, betasize      = vsize + msize1* msize2;
  size_t  iter = 0;
  gsl_vector *x;
    apop_mle_settings *mlep   = dist->method_settings;
    dist->status = 1;    //assume failure until we score a success.
    if (!mlep || mlep->starting_pt==NULL){
        x = gsl_vector_alloc(betasize);
        gsl_vector_set_all (x,  2);
    } else
        x   = apop_array_to_vector(mlep->starting_pt, betasize);
  gsl_multiroot_function f = {dnegshell, betasize, &p};
    if (mlep->method == APOP_RF_NEWTON)
        T = gsl_multiroot_fsolver_dnewton;
    else if (mlep->method == APOP_RF_HYBRID_NOSCALE)
        T = gsl_multiroot_fsolver_hybrids;
    else //if (mlep->method == APOP_RF_HYBRID)        --default
        T = gsl_multiroot_fsolver_hybrid;
    s = gsl_multiroot_fsolver_alloc (T, betasize);
    gsl_multiroot_fsolver_set (s, &f, x);
    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);
        if (!mlep || mlep->verbose)
            printf ("iter = %3u x = % .3f f(x) = % .3e\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->f, 0));
        if (status)   /* check if solver is stuck */
            break;
        status = gsl_multiroot_test_residual (s->f, mlep->tolerance);
     } while (status == GSL_CONTINUE && iter < 1000);
     if(GSL_SUCCESS) dist->status = 0;
  printf ("status = %s\n", gsl_strerror (status));
  //dist->parameters   = apop_data_unpack(s->x, vsize, msize1, msize2);
  apop_data_unpack(s->x, dist->parameters);
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  return dist;
}
