/** \file apop_probit.c

Copyright (c) 2005--2008 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

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
    apop_name_cross_stack(m->parameters->names, d->names, 'r', 'c');
}

static double probit_log_likelihood(apop_data *d, apop_model *p){
  apop_assert(p->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
  int		    i;
  long double	n, total_prob	= 0;
  apop_data *betadotx = apop_dot(d, p->parameters, 0, 'v'); 
	for(i=0; i< d->matrix->size1; i++){
		n	        = gsl_cdf_gaussian_P(gsl_vector_get(betadotx->vector,i),1);
        n = n ? n : 1e-10; //prevent -inf in the next step.
        n = n<1 ? n : 1-1e-10; 
        total_prob += apop_data_get(d, i, -1) ?  log(n): log(1 - n);
	}
    apop_data_free(betadotx);
	return total_prob;
}

static void probit_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
  apop_assert_void(p->parameters, 0,'s', "You asked me to evaluate an un-parametrized model.");
  int		i, j;
  long double	cdf, betax, deriv_base;
  apop_data *betadotx = apop_dot(d, p->parameters, 0, 'v'); 
    gsl_vector_set_all(gradient,0);
    for (i=0; i< d->matrix->size1; i++){
        betax            = gsl_vector_get(betadotx->vector, i);
        cdf              = gsl_cdf_gaussian_P(betax, 1);
        cdf = cdf ? cdf : 1e-10; //prevent -inf in the next step.
        cdf = cdf<1 ? cdf : 1-1e-10; 
        if (apop_data_get(d, i, -1))
            deriv_base      = gsl_ran_gaussian_pdf(betax, 1) /  cdf;
        else
            deriv_base      = -gsl_ran_gaussian_pdf(betax, 1) /(1-cdf);
        for (j=0; j< p->parameters->vector->size; j++)
            apop_vector_increment(gradient, j, apop_data_get(d, i, j) * deriv_base);
	}
	apop_data_free(betadotx);
}



/** The Probit model.
 The first column of the data matrix this model expects is ones and zeros;
 the remaining columns are values of the independent variables. Thus,
 the model will return (data columns)-1 parameters.

\ingroup models
*/
apop_model apop_probit = {"Probit", -1,0,0, .log_likelihood = probit_log_likelihood, 
    .score = probit_dlog_likelihood, .prep = probit_prep};
//estimate via the default MLE.
