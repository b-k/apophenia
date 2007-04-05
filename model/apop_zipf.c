/** \file apop_zipf.c

  The Zipf distribution.

\f$Z(a)        = {1\over \zeta(a) * i^a}        \f$<br>

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/


//The default list. Probably don't need them all.
#include "db.h" //apop_opts
#include "types.h"
#include "model.h"
#include "output.h"
#include "mapply.h"
#include "conversions.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_zeta.h>


static apop_params * zipf_estimate(apop_data * data, apop_params *parameters){
    return apop_maximum_likelihood(data, apop_zipf, parameters);
}

static double beta_greater_than_x_constraint(const apop_data *beta, apop_data *returned_beta, apop_params *v){
    //constraint is 1 < beta_1
  static apop_data *constraint = NULL;
    if (!constraint) constraint = apop_data_calloc(1,1,1);
    apop_data_set(constraint, 0, 0, 1);
    apop_data_set(constraint, 0, -1, 1);
    return apop_linear_constraint(beta->vector, constraint, 1e-3, returned_beta->vector);
}

static double oneline_log(gsl_vector *v){
  int       j;
  double    like    = 0;
    for(j=0; j< v->size; j++)
        like    -= log(gsl_vector_get(v,j));
    return like;
}

static double zipf_log_likelihood(const apop_data *beta, apop_data *d, apop_params *v){
  gsl_matrix    *data   = d->matrix;
  long double   bb      = gsl_vector_get(beta->vector, 0);
  gsl_vector    *logs   = apop_matrix_map(data, oneline_log);
  long double   like    = apop_vector_sum(logs);
    like    *= bb;
    like    -= log(gsl_sf_zeta(bb)) * data->size1 * data->size2;
    gsl_vector_free(logs);
    return like;
}    

static double zipf_p(const apop_data *beta, apop_data *d, apop_params *v){
    return exp(zipf_log_likelihood(beta, d, v));
}    

static void zipf_dlog_likelihood(const apop_data *beta, apop_data *d, gsl_vector *gradient, apop_params *v){
  double      bb        = gsl_vector_get(beta->vector, 0);
  gsl_matrix  *data     = d->matrix;
  gsl_vector  *logs     = apop_matrix_map(data, oneline_log);
  long double dlike     = apop_vector_sum(logs);
    dlike   -= gsl_sf_zeta(bb-1)/(bb*gsl_sf_zeta(bb))  * data->size1 * data->size2;
    gsl_vector_set(gradient,0,dlike);
    gsl_vector_free(logs);
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
static void zipf_rng( double *out, gsl_rng* r, apop_params *param){
  double a  = gsl_vector_get(param->parameters->vector, 0);
    if (a  <= 1){
//        if (apop_opts.verbose)
            fprintf(stderr, "apop_zipf.rng: Zipf needs a parameter >=1. Returning 0.\n"); 
            *out    = 0;
        return;
    }
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

apop_zipf.estimate() is an MLE, so feed it appropriate \ref apop_params.

\f$Z(a)        = {1\over \zeta(a) * i^a}        \f$

\f$lnZ(a)    = -(\log(\zeta(a)) + a \log(i))    \f$

\f$dlnZ(a)/da    = -{a \zeta(a)\over\log(\zeta(a-1))} -  \log(i)        \f$
\ingroup models
*/
apop_model apop_zipf = {"Zipf", 1,0,0, 
    zipf_estimate, zipf_p, zipf_log_likelihood, zipf_dlog_likelihood, beta_greater_than_x_constraint, zipf_rng};
