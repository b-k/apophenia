/** \file apop_zipf_rank.c

  The Zipf distribution, rank data.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/


//The default list. Probably don't need them all.
#include "db.h" //apop_opts
#include "types.h"
#include "model.h"
#include "output.h"
#include "conversions.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <assert.h>


static apop_model * zipf_estimate(apop_data * data, apop_model *parameters){
    return apop_maximum_likelihood(data, *parameters);
}

static double beta_greater_than_x_constraint(const apop_data *beta, apop_data *returned_beta, apop_model *v){
    //constraint is 1 < beta_1
  static apop_data *constraint = NULL;
    if (!constraint) constraint = apop_data_calloc(1,1,1);
    apop_data_set(constraint, 0, 0, 1);
    apop_data_set(constraint, 0, -1, 1);
    return apop_linear_constraint(beta->vector, constraint, 1e-3, returned_beta->vector);
}

#include <gsl/gsl_sf_zeta.h>

static double zipf_log_likelihood(const apop_data *beta, apop_data *d, apop_model *params){
  long double     like    = 0, 
                  a       = gsl_vector_get(beta->vector, 0);
  int             j;
    for(j=0; j< d->matrix->size2; j++){
        APOP_COL(d, j, v);
        like   -= apop_sum(v) * log(j+1);
    }
    like    *= a;
    like    -= log(gsl_sf_zeta(a)) * d->matrix->size1 * d->matrix->size2;
    return like;
}    

static double zipf_p(const apop_data *beta, apop_data *d, apop_model *params){
    return exp(zipf_log_likelihood(beta, d, params));
}

static void zipf_dlog_likelihood(const apop_data *beta, apop_data *d, gsl_vector *gradient, apop_model *params){
  long double     a       = gsl_vector_get(beta->vector, 0),
                  dlike   = 0;
  int             j;
    for(j=0; j< d->matrix->size2; j++){
        APOP_COL(d, j, v);
        dlike   -= apop_sum(v) * log(j+1);
    }
    dlike   -= (gsl_sf_zeta(a-1)/(a*gsl_sf_zeta(a))) * d->matrix->size1 * d->matrix->size2;
    gsl_vector_set(gradient,0,dlike);
}    


/** Draw from a Zipf distribution with parameter \f$ a \f$

Returns a ranking: If the population were Zipf distributed, you're most
likely to get the 1st most common item, so this produces a lot of ones,
a great deal of twos, and so on.


For example:
\code
gsl_rng *       r;
gsl_rng_env_setup();
r=gsl_rng_alloc(gsl_rng_taus);    //for example. 
apop_zipf.rng(r, 1.4);
\endcode

Cribbed from <a href="http://cgm.cs.mcgill.ca/~luc/mbookindex.html>Devroye (1986)</a>, Chapter 10, p 551.  */
static void  zipf_rng( double *out, gsl_rng* r, apop_model *a){
  double    p   =a->parameters->vector->data[0];
    if (p <= 1){
//        if (apop_opts.verbose)
            printf("apop_zipf.rng: Zipf needs a parameter >=1. Returning 0.\n"); 
            *out    = 0;
        return;
    }
  int     x;
  double  u, v, t, 
        b       = pow(2, p -1), 
        ainv    = -(1.0/(p -1));
    do {
        u    = gsl_rng_uniform(r);
        v    = gsl_rng_uniform(r);
    //    x    = GSL_MIN(pow(u, ainv), 1e8); //prevent overflows.
        x    = pow(u, ainv);
        t    = pow((1.0 + 1.0/x), (p -1));
    } while (v * x * (t-1.0)/(b-1) > t/b);
    *out = x;
}



/** The Zipf distribution.
Wikipedia has notes on the <a href="http://en.wikipedia.org/wiki/Zipf_distribution">Zipf distribution</a>. 

The data set needs to be in rank-form. The first column is the frequency of the most common item, the second is the frequency of the second most common item, &c.

apop_zipf.estimate() is an MLE, so feed it appropriate \ref apop_params.

\f$Z(a)        = {1\over \zeta(a) * i^a}        \f$

\f$lnZ(a)    = -(\log(\zeta(a)) + a \log(i))    \f$

\f$dlnZ(a)/da    = -{\zeta(a)\over a \log(\zeta(a-1))} -  \log(i)        \f$
\ingroup models
*/
apop_model apop_zipf_rank = {"Zipf, rank data", 1,0,0, 
    .estimate = zipf_estimate, .p = zipf_p, .log_likelihood = zipf_log_likelihood, 
    .score = zipf_dlog_likelihood, .constraint = beta_greater_than_x_constraint, 
    .draw = zipf_rng};
