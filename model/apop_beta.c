/** \file apop_binomial.c 
 
Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"

//The default list. You probably don't need them all.
#include "types.h"
#include "conversions.h"
#include "likelihoods.h"
#include "model.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <stdio.h>
#include <assert.h>


/** The Beta distribution is useful for modeling because it is bounded
between zero and one, and can be either unimodal (if the variance is low)
or bimodal (if the variance is high), and can have either a slant toward
the bottom or top of the range (depending on the mean).

The distribution has two parameters, typically named \f$\alpha\f$ and \f$\beta\f$, which
can be difficult to interpret. However, there is a one-to-one mapping
between (alpha, beta) pairs and (mean, variance) pairs. Since we have
good intuition about the meaning of means and variances, this function
takes in a mean and variance, calculates alpha and beta behind the scenes,
and returns a random draw from the appropriate Beta distribution.

\param m
The mean the Beta distribution should have. Notice that m
is in [0,1].

\param v
The variance which the Beta distribution should have. It is in (0, 1/12),
where (1/12) is the variance of a Uniform(0,1) distribution. The closer
to 1/12, the worse off you are.

\param r
An already-declared and already-initialized {{{gsl_rng}}}.

\return
Returns an \c apop_beta model with its parameters appropriately set.

*/
apop_model *apop_beta_from_mean_var(double m, double v){
    double k     = (m * (1- m)/ v) -1 ;
    double alpha = m*k;
    double beta  = k * (1-m);
    return apop_model_set_parameters(apop_beta, alpha, beta);
}


static double beta_log_likelihood(apop_data *d, apop_model *p);

/* Every model should have an estimate function. If you are doing an
 MLE, then it should be the one line function here, and then you will need
 to fill in the beta_log_likelihood function below. If you are not
 doing an MLE, then you won't need the beta_log_likelihood function,
 but will probably instead do some substantial math here in this function.

 a
*/
static apop_model * beta_estimate(apop_data * data,  apop_model *parameters){
  apop_model 	*est= parameters ? parameters : apop_model_copy(apop_beta);
  apop_model_clear(data, est);
  double		mean, var, alpha, beta;
	apop_matrix_mean_and_var(data->matrix, &mean, &var);	
    alpha   = gsl_pow_2(mean) * ((1-mean)/var - 1/mean);
    beta    = alpha * (1-mean)/mean;
	gsl_vector_set(est->parameters->vector, 0, alpha);
	gsl_vector_set(est->parameters->vector, 1, beta);
    est->llikelihood	= beta_log_likelihood(data, parameters);
    //apop_numerical_covariance_matrix(apop_beta, est, data);
	return est;
}


/*
Often, the input data is a gsl_matrix, and you need to go line by line
through the matrix, calculating the log likelihood. That's the sample
code here.  
*/
static double beta_log_likelihood(apop_data *d, apop_model *p){
    if (!p->parameters)
        apop_error(0,'s', "%s: You asked me to evaluate an un-parametrized model.", __func__);
int		    i,j;
double	    x, 
            loglike    	= 0,
            alpha       = apop_data_get(p->parameters,0,-1),
            beta        = apop_data_get(p->parameters,1,-1);
gsl_matrix 	*data 		= d->matrix;
	for(i=0;i< data->size1; i++)
        for(j=0; j< data->size2; j++){
            x       = gsl_matrix_get(data, i, j);
		    loglike += alpha * log(x) + beta *log(1-x);
        }
	return loglike  + gsl_sf_lnbeta (alpha, beta)*data->size1*data->size2;
}

static double beta_p(apop_data *d, apop_model *p){
    return exp(beta_log_likelihood(d, p));
}

/* The derivative of the beta distribution, for use in likelihood
  minimization. 
  The format is often the same as above: go line by line through a gsl_matrix.
  The sample is a three-dimensional parameter vector.

You can delete this function entirely if so inclined. If so, remember
to replace this function with NULL in the model definition below.
static void beta_dlog_likelihood(gsl_vector *beta, apop_data *d, gsl_vector *gradient){
int		    i,j;
double	    dtotal[3];
gsl_matrix 	*data 	= d->matrix;
    dtotal[0]  = 0,
    dtotal[1]  = 0,
    dtotal[2]  = 0;
	for(i=0; i< data->size1; i++){
		dtotal[0]  += 0; //PLACE MATH HERE
		dtotal[1]  += 0; //PLACE MATH HERE
		dtotal[2]  += 0; //PLACE MATH HERE
	}
	for(j=0; j< beta->size; j++){
	    gsl_vector_set(gradient,j,dtotal[j]);
    }
}
 */


/* For constrained optimizations, you will need a constraint function.
 You can use \c apop_linear_constraint if your constraint can be expressed as a set of hyperplanes.

You can delete this function entirely if so inclined. If so, remember
to replace this function with NULL in the model definition below.
 */
static double beta_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_1 and  0 < beta_2
  static apop_data *constraint = NULL;
    if (!constraint){
        constraint= apop_data_alloc(2,2,1);
        apop_data_fill(constraint, 0., 1., 0.,
                                   0., 0., 1.);
    }
    return apop_linear_constraint(v->parameters->vector, constraint, 1e-3);
}


static void beta_rng(double *out, gsl_rng *r, apop_model* eps){
    *out = gsl_ran_beta(r, apop_data_get(eps->parameters,0,-1), apop_data_get(eps->parameters,1,-1));
}

/** The beta model.
You should describe the format of the input data here.

--The second parameter is either the number of parameters your model has, or a negative number indicating a special case; see the manual.
--If you deleted any of the functions above, replace their names with NULL here.

\ingroup models
*/
apop_model apop_beta = {"Beta distribution", 2,0,0,
	.estimate = beta_estimate, .p = beta_p, .log_likelihood = beta_log_likelihood, 
    .constraint = beta_constraint, .draw = beta_rng};
