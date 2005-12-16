/** \file apop_zipf.c

  The Zipf distribution.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/


//The default list. Probably don't need them all.
#include "types.h"
#include "model.h"
#include "output.h"
#include "conversions.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <assert.h>




static apop_estimate * zipf_estimate(gsl_matrix * data, apop_inventory *uses, void *parameters){
    apop_inventory_filter(uses, apop_zipf.inventory_filter);
    return apop_maximum_likelihood(data, uses, apop_zipf, *(apop_estimation_params *)parameters);
}

static double beta_greater_than_x_constraint(gsl_vector *beta, void * d, gsl_vector *returned_beta){
double  limit       = 1,
        tolerance   = 1e-1;
double  mu          = gsl_vector_get(beta, 0);
    if (mu > limit) 
        return 0;
    //else:
    gsl_vector_memcpy(returned_beta, beta);
    gsl_vector_set(returned_beta, 0, limit + tolerance);
    return limit - mu;    
}

///////////////////////
//The Zipf distribution
///////////////////////
#include <gsl/gsl_sf_zeta.h>

/* The Zipf distribution.

\f$Z(a)        = {1\over \zeta(a) * i^a}        \f$<br>

 \todo link this fn in with the object 
static double apop_zipf_likelihood(double a, int i){
double        z    = 1/(gsl_sf_zeta(a) * pow(i, a));
    return z;
}
*/

static double zipf_log_likelihood(const gsl_vector *beta, void *d){
gsl_matrix  *data   = d;
double      like    = 0, 
            bb      = gsl_vector_get(beta, 0);
int         i, j;
    for(j=0; j< data->size2; j++)
        for(i=0; i< data->size1; i++)
            like    -= log(gsl_matrix_get(data,i,j));
    like    *= bb;
    like    += -log(gsl_sf_zeta(bb)) * data->size1 * data->size2;
    return like;
}    

static void zipf_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
double      a       = gsl_vector_get(beta, 0);
gsl_matrix  *data   = d;
int         i, j;
double      dlike   = 0;
    for(j=0; j< data->size2; j++)
        for(i=0; i< data->size1; i++)
            dlike   += a*gsl_sf_zeta(a+1)/gsl_sf_zeta(a) - log(gsl_matrix_get(data,i,j));
    gsl_vector_set(gradient,0,dlike);
}    


/** Draw from a Zipf distribution with parameter \f$ a \f$

Call this fn using \ref apop_zipf.rng().

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

Cribbed from <a href="http://cgm.cs.mcgill.ca/~luc/mbookindex.html>Devroye (1986)</a>, p 551.  */
static double zipf_rng(gsl_rng* r, double * a){
    if (*a  <= 1){
        if (apop_verbose)
            printf("apop_zipf.rng: Zipf needs a parameter >=1. Returning 0.\n"); 
        return 0;
        }
int     x;
double  u, v, t, 
        b       = pow(2, *a-1), 
        ainv    = -(1.0/(*a-1));
    do {
        u    = gsl_rng_uniform(r);
        v    = gsl_rng_uniform(r);
        x    = pow(u, ainv);
        t    = pow((1.0 + 1.0/x), (*a-1));
    } while (v * x * (t-1.0)/(b-1) > t/b);
    return x;
}




/** The Zipf distribution.
Wikipedia has notes on the <a href="http://en.wikipedia.org/wiki/Zipf_distribution">Zipf distribution</a>. 

Ignores the matrix structure of the input data, so send in a 1 x N, an N x 1, or an N x M.

apop_zipf.estimate() is an MLE, so feed it appropriate \ref apop_estimation_params.

\f$Z(a)        = {1\over \zeta(a) * i^a}        \f$

\f$lnZ(a)    = -(\log(\zeta(a)) + a \log(i))    \f$

\f$dlnZ(a)/da    = -{\zeta(a)\over a \log(\zeta(a-1))} -  \log(i)        \f$
\ingroup models
*/
apop_model apop_zipf = {"Zipf", 1,  {
    1,    //parameters
    1,    //covariance
    1,    //confidence
    0,    //predicted
    0,    //residuals
    1,    //log_likelihood
    0    //names;
},         
    zipf_estimate, zipf_log_likelihood, zipf_dlog_likelihood,  NULL, beta_greater_than_x_constraint, zipf_rng};
