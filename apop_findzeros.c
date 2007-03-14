#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <apophenia/headers.h>

/** \file apop_findzero.c
 This just includes the root-finding routine. It is #included in apop_mle, because I expect you to call it via that. */

typedef struct {
  apop_model *model;
  void *model_params;
  apop_data* data;
} score_struct;

int score_shell (const gsl_vector * in, void *params, gsl_vector * out) {
  score_struct *p = params;
    p->model->score(in, p->data, out, p->model_params);
    return GSL_SUCCESS;
}

/** This function is cut/pasted/modified from the GSL documentation. It
 calls the various GSL root-finding algorithms to find the zero of the score.
*/
static apop_estimate * find_roots (apop_data * data, apop_model dist, apop_ep *est_params) {
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status, betasize  = dist.parameter_ct;
  size_t  iter = 0;
  gsl_vector *x;
    if (betasize == -1) {
        dist.parameter_ct   =
        betasize            = data->matrix->size2 - 1;
    }
    apop_estimate *est = apop_estimate_alloc(data, dist, est_params);
    est->status = 1;    //assume failure until we score a success.
    if (!est_params || est_params->starting_pt==NULL){
        x = gsl_vector_alloc(betasize);
        gsl_vector_set_all (x,  2);
    } else
        x   = apop_array_to_vector(est_params->starting_pt, betasize);
  score_struct      p;
    p.data            = data;
    p.model           = &dist;
    p.model_params    = est_params;
  gsl_multiroot_function f = {score_shell, betasize, &p};
    if (est_params->method == 13)
        T = gsl_multiroot_fsolver_dnewton;
    else if (est_params->method == 11)
        T = gsl_multiroot_fsolver_broyden;
    else if (est_params->method == 12)
        T = gsl_multiroot_fsolver_hybrids;
    else //if (est_params->method == 10)        --default
        T = gsl_multiroot_fsolver_hybrid;
    s = gsl_multiroot_fsolver_alloc (T, betasize);
    gsl_multiroot_fsolver_set (s, &f, x);
    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);
        if (!est_params || est_params->verbose)
            printf ("iter = %3u x = % .3f f(x) = % .3e\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->f, 0));
        if (status)   /* check if solver is stuck */
            break;
        status = gsl_multiroot_test_residual (s->f, 1e-7);
     } while (status == GSL_CONTINUE && iter < 1000);
     if(GSL_SUCCESS) est->status = 0;
  printf ("status = %s\n", gsl_strerror (status));
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  return est;
}

