/** \file apop_probit.c

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"

//The default list. Probably don't need them all.
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


/*
//cut 'n' pasted from apop_regression.c
static void prep_names (apop_model *e){
  int i;
	if (e->data->names->colct > 0) {		
		//apop_name_add(n, n->column[0], 'd');
        apop_name_add(e->expected->names, e->data->names->column[0], 'c');
        apop_name_add(e->expected->names, "predicted", 'c');
        apop_name_add(e->expected->names, "residual", 'c');
        if (e->parameters)
            snprintf(e->parameters->names->title, 100, "Regression of %s", e->data->names->column[0]);
		sprintf(e->data->names->column[0], "1");
        apop_name_add(e->parameters->names, "1", 'r');
        apop_name_add(e->parameters->names, "parameters", 'v');
        for(i=1; i< e->data->names->colct; i++)
            apop_name_add(e->parameters->names, e->data->names->column[i], 'r');
        if (p->want_cov){
            if (e->data->names){
                apop_name_stack(e->covariance->names, e->data->names, 'c');
                apop_name_cross_stack(e->covariance->names, e->data->names, 'c', 'r');
            }
		    sprintf(e->covariance->names->column[0], "1");
		    sprintf(e->covariance->names->row[0], "1");
        }
	}
}
*/


static void probit_prep(apop_data *d, apop_model *m){
    if (!d->vector){
        APOP_COL(d, 0, independent);
        d->vector = apop_vector_copy(independent);
        gsl_vector_set_all(independent, 1);
        if (d->names->colct > 0) {		
            apop_name_add(d->names, d->names->column[0], 'v');
            sprintf(d->names->column[0], "1");
        }
    }
    void *mpt = m->prep; //and use the defaults.
    m->prep = NULL;
    apop_model_prep(d, m);
    m->prep = mpt;
}

static apop_model * probit_estimate(apop_data * data,  apop_model *parameters){
    //modify_in_data(data);
	return apop_maximum_likelihood(data, *parameters);
}

//////////////////
//The probit model
//////////////////

//This section includes some trickery to avoid recalculating beta dot x.
//
static gsl_vector *beta_dot_x 		= NULL;
static gsl_matrix  *last_data_set   = NULL;
static int	beta_dot_x_is_current	= 0;

static void	dot(apop_data *beta, gsl_matrix *data){
	//if (beta_dot_x) printf("comparing %i with %i ",data->size1, beta_dot_x->size); fflush(NULL);
	if ((data != last_data_set) && beta_dot_x){
		gsl_vector_free(beta_dot_x); 
		beta_dot_x = NULL;
    }
	if (!beta_dot_x){
        last_data_set   = data;
		beta_dot_x 	    = gsl_vector_alloc(data->size1);		//global var
    }
    gsl_blas_dgemv (CblasNoTrans, 1.0, data, beta->vector, 0.0, beta_dot_x);	//dot product
}

/*
find (data dot beta'), then find the integral of the \f$\cal{N}(0,1)\f$
up to that point. Multiply likelihood either by that or by 1-that, depending 
on the choice the data made.
*/
static double probit_log_likelihood(apop_data *d, apop_model *p){
  if (!p->parameters)
      apop_error(0,'s', "%s: You asked me to evaluate an un-parametrized model.", __func__);
  int		    i;
  long double	n, total_prob	= 0;
  gsl_matrix    *data           = d->matrix;
	dot(p->parameters, data);
	for(i=0; i< data->size1; i++){
		n	        = gsl_cdf_gaussian_P(gsl_vector_get(beta_dot_x,i),1);
        n = n ? n : 1e-10; //prevent -inf in the next step.
        n = n<1 ? n : 1-1e-10; 
        total_prob += apop_data_get(d, i, -1)==0 ?  log(n): log(1 - n);
	}
	return total_prob;
}

static double probit_p(apop_data *d, apop_model *p){
    return exp(probit_log_likelihood(d, p));
}

/* The derivative of the probit distribution, for use in likelihood
  minimization. You'll probably never need to call this directly.*/
static void probit_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
	//derivative of the above. 
  if (!p->parameters)
      apop_error(0,'s', "%s: You asked me to evaluate an un-parametrized model.", __func__);
  int		i, j;
  long double	one_term, beta_term_sum, cdf;
  gsl_matrix 	*data 		= d->matrix;
	if (!beta_dot_x_is_current) 	
		dot(p->parameters,data); 
	for(j=0; j< p->parameters->vector->size; j++){
		beta_term_sum	= 0;
		for(i=0; i< data->size1; i++){
			one_term	     = gsl_matrix_get(data, i,j)
						        * gsl_ran_gaussian_pdf(gsl_vector_get(beta_dot_x,i),1);
            cdf              = gsl_cdf_gaussian_P(gsl_vector_get(beta_dot_x,i),1);
        cdf = cdf ? cdf : 1e-10; //prevent -inf in the next step.
        cdf = cdf<1 ? cdf : 1-1e-10; 
            one_term        /= apop_data_get(d, i, -1)==0 ? cdf : cdf-1;
			beta_term_sum	+= one_term;
		}
        gsl_vector_set(gradient,j,beta_term_sum);
	}
	gsl_vector_free(beta_dot_x);
	beta_dot_x	= NULL;
}


/** Saves some time in calculating both log likelihood and dlog
likelihood for probit.	
static void probit_fdf(gsl_vector *beta, apop_data *d, double *f, gsl_vector *df){
	*f	= probit_log_likelihood(beta, d, NULL);
	beta_dot_x_is_current	=1;
	probit_dlog_likelihood(beta, d, df, NULL);
	beta_dot_x_is_current	=0;
}*/

/** The Probit model.
 The first column of the data matrix this model expects is ones and zeros;
 the remaining columns are values of the independent variables. Thus,
 the model will return (data columns)-1 parameters.

\ingroup models
*/
apop_model apop_probit = {"Probit", -1,0,0, 
    .estimate = probit_estimate, .p = probit_p, .log_likelihood = probit_log_likelihood, 
    .score = probit_dlog_likelihood, .prep = probit_prep};
