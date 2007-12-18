/** \file apop_probit.c

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"

//The default list. Probably don't need them all.
#include "types.h"
#include "conversions.h"
#include "likelihoods.h"
#include "model.h"
#include "mapply.h"
#include "linear_algebra.h"
#include <gsl/gsl_sort.h>

static int multiprobit_count(apop_data *d, apop_model *m){
    int count = 0;
    int i, j;
    double *vals = NULL;
    for(i=0; i < d->vector->size; i ++){
        double current = gsl_vector_get(d->vector, i);
        for(j=0; j < count && current != vals[j]; j ++)
            ;
        if (j==count){//not found
            vals = realloc(vals, sizeof(double) * (count+1));
            vals[count++] = current;
        }
    }
    vals ++;//drop the first entry---a one-double memory leak!
    gsl_sort(vals, 1, --count);
    m->more = vals;
    return count;
}

static void multiprobit_prep(apop_data *d, apop_model *m){
    apop_probit.prep(d, m); 
    apop_data_free(m->parameters);
    int count = multiprobit_count(d, m);
    m->parameters = apop_data_alloc(0, d->matrix->size2, count);
}


static double val;
static double onerow(double in){
    return in >= val;
}

/*
This is just a for loop that runs a probit on each row.
*/
static double multiprobit_log_likelihood(apop_data *d, apop_model *p){
  if (!p->parameters)
      apop_error(0,'s', "%s: You asked me to evaluate an un-parametrized model.", __func__);

  static apop_model *spare_probit = NULL;
    if (!spare_probit){
        spare_probit = apop_model_copy(apop_probit);
        spare_probit->parameters = apop_data_alloc(0,0,0);
        spare_probit->prepared++;
    }
  static apop_data *working_data = NULL;
  if (!working_data)
      working_data = apop_data_alloc(0,0,0);

  working_data->matrix = d->matrix;

  gsl_vector *original_outcome = d->vector;
  int    i;
  double ll    = 0;
  double *vals = p->more;
    for(i=0; i < p->parameters->matrix->size2; i++){
        APOP_COL(p->parameters, i, param);
        val = vals[i];
        working_data->vector = apop_vector_map(original_outcome, onerow);
        spare_probit->parameters->vector = param;
        ll  += apop_log_likelihood(working_data, spare_probit);
        gsl_vector_free(working_data->vector); //yup. It's inefficient.
    }
	return ll;
}


/** The Multinomial Probit model.
 The first column of the data matrix this model expects is ones and zeros;
 the remaining columns are values of the independent variables. Thus,
 the model will return (data columns) parameters.

\ingroup models
*/
apop_model apop_multinomial_probit = {"Multinomial probit", -1,0,0, 
     .log_likelihood = multiprobit_log_likelihood, .prep = multiprobit_prep};
