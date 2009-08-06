/** \file apop_zipf.c

  The Zipf distribution.

\f$Z(a)        = {1\over \zeta(a) * i^a}        \f$<br>

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

//The default list. Probably don't need them all.
#include "db.h" //apop_opts
#include "types.h"
#include "model.h"
#include "output.h"
#include "mapply.h"
#include "settings.h"
#include "conversions.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_zeta.h>

static double zrank_apply(double in, int k) { return in * log(k+1); }
#define apop_data_size(d) ((d) ? (((d)->vector ? (d)->vector->size : 0) + ((d)->matrix ? (d)->matrix->size1 * (d)->matrix->size2 : 0)) : 0)

//There are two versions: the rank version and the non-rank. Here's
//Rank's fns:
static double zipf_log_likelihood_rank(const apop_data *d, apop_model *m){
  long double  a       = gsl_vector_get(m->parameters->vector, 0);
  long double like = -a * apop_map_sum((apop_data *)d, .fn_di=zrank_apply, .part ='c');
    like    -= log(gsl_sf_zeta(a)) * d->matrix->size1 * d->matrix->size2;
    return like;
}    

static void zipf_dlog_likelihood_rank(const apop_data *d, gsl_vector *gradient, apop_model *m){
  long double     a       = gsl_vector_get(m->parameters->vector, 0);
  int         size = apop_data_size(d);
  long double dlike =  dlike = apop_map_sum((apop_data *)d, .fn_di = zrank_apply, .part='c');
    dlike   -= (gsl_sf_zeta(a-1)/(a*gsl_sf_zeta(a))) * size;
    gsl_vector_set(gradient,0,dlike);
}    

static double beta_greater_than_x_constraint(apop_data *returned_beta, apop_model *m){
    //constraint is 1 < beta_1
  static apop_data *constraint = NULL;
    if (!constraint) {
        constraint = apop_data_calloc(1,1,1);
        apop_data_set(constraint, 0, 0, 1);
        apop_data_set(constraint, 0, -1, 1);
    }
    return apop_linear_constraint(m->parameters->vector, constraint, 1e-3);
}

static double zipf_log_likelihood(apop_data *d, apop_model *m){
  apop_assert(m->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
    if (apop_settings_get_group(m, "apop_rank"))
        return zipf_log_likelihood_rank(d, m);
  int size = apop_data_size(d);
  long double   bb      = gsl_vector_get(m->parameters->vector, 0);
  double like = -apop_map_sum(d, log);
    like    *= bb;
    like    -= log(gsl_sf_zeta(bb)) * size;
    return like;
}    

static double zipf_p(apop_data *d, apop_model *v){
    if (apop_settings_get_group(v, "apop_rank"))
        return exp(zipf_log_likelihood_rank(d, v));
    return exp(zipf_log_likelihood(d, v));
}    

static void zipf_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
  apop_assert_void(m->parameters, 0,'s', "You asked me to evaluate an un-parametrized model.");
    if (apop_settings_get_group(m, "apop_rank"))
        return zipf_dlog_likelihood_rank(d, gradient, m);
  int         size      = apop_data_size(d);
  double      bb        = gsl_vector_get(m->parameters->vector, 0);
  long double dlike     =  -apop_map_sum(d, log);
    dlike   -= gsl_sf_zeta(bb-1)/(bb*gsl_sf_zeta(bb))  * size;
    gsl_vector_set(gradient,0,dlike);
}    

/** Draw from a Zipf distribution with parameter \f$ a \f$

Call this fn using \ref apop_zipf.rng().

Returns a ranking: If the population were Zipf distributed, you're most
likely to get the 1st most common item, so this produces a lot of ones,
a great deal of twos, and so on.

In the interest of avoiding overflows, the RNG is capped at 1e8.

For example:
\code
gsl_rng *r  = apop_rng_alloc(r);
double  d   = 1.4
gsl_vector *params = apop_array_to_vector(&d,1);
apop_zipf.draw(r, 1.4, NULL);
\endcode

Cribbed from <a href="http://cgm.cs.mcgill.ca/~luc/mbookindex.html>Devroye (1986)</a>, Chapter 10, p 551.  */
static void zipf_rng(double *out, gsl_rng* r, apop_model *param){
  double a  = gsl_vector_get(param->parameters->vector, 0);
  apop_assert_void(a >= 1, 0, 's', "Zipf needs a parameter >=1. Stopping."); 
  int     x;
  double  u, v, t, 
          b       = pow(2, a-1), 
          ainv    = -(1.0/(a-1));
    do {
        u    = gsl_rng_uniform(r);
        v    = gsl_rng_uniform(r);
        x    = pow(u,ainv); //GSL_MIN(pow(u, ainv), 1e8); //prevent overflows.
        t    = pow((1.0 + 1.0/x), (a-1));
    } while (v * x * (t-1.0)/(b-1) > t/b);
    *out = x;
}


/** The Zipf distribution.
Wikipedia has notes on the <a href="http://en.wikipedia.org/wiki/Zipf_distribution">Zipf distribution</a>. 

Ignores the matrix structure of the input data, so send in a 1 x N, an N x 1, or an N x M.

apop_zipf.estimate() is an MLE, so feed it appropriate \ref apop_mle_settings.
\f$Z(a)        = {1\over \zeta(a) * i^a}        \f$

\f$lnZ(a)    = -(\log(\zeta(a)) + a \log(i))    \f$

\f$dlnZ(a)/da    = -{a \zeta(a)\over\log(\zeta(a-1))} -  \log(i)        \f$

To specify that you have frequency or ranking data, use 
\code
Apop_settings_add_group(your_model, apop_rank, NULL);
\endcode

\ingroup models
*/
apop_model apop_zipf = {"Zipf", 1,0,0, 
    .p = zipf_p, .log_likelihood = zipf_log_likelihood,
    .score = zipf_dlog_likelihood, .constraint = beta_greater_than_x_constraint, 
    .draw = zipf_rng};
