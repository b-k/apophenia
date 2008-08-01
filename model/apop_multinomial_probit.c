/** \file apop_multinomial_probit.c

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
    //probit uses a numeraire; logit doesn't.
    //if (!strcmp(m->name, "Multinomial probit")){
        gsl_sort(vals, 1, count);
        vals  ++;//drop the first entry---a one-double memory leak!
        count --;
    //} else
    //    gsl_sort(vals, 1, count);
    m->more = vals;
    return count;
}

static void multiprobit_prep(apop_data *d, apop_model *m){
  int i;
  char *s = NULL;
    apop_probit.prep(d, m); 
    apop_data_free(m->parameters);
    int count = multiprobit_count(d, m);
    m->parameters = apop_data_alloc(0, d->matrix->size2, count);
    apop_name_cross_stack(m->parameters->names, d->names, 'r', 'c');
    double *vals = m->more;
    for (i=0; i< count; i++){
        asprintf(&s, "option %g", vals[i]);
        apop_name_add(m->parameters->names, s, 'c');
    }
    gsl_matrix_set_all(m->parameters->matrix, 1);
}


static apop_data *multilogit_expected(apop_data *in, apop_model *m){
  int i, j;
  apop_assert(m->parameters, NULL, 0, 's', "You're asking me to provide expected values of an un-parameterized model. Please run apop_estimate first.");
    apop_model_prep(in, m);
    gsl_matrix *params = m->parameters->matrix;
    apop_data *out = apop_data_alloc(in->matrix->size1, in->matrix->size1, params->size2+1);
    for (i=0; i < in->matrix->size1; i ++){
        Apop_row(in, i, observation);
        Apop_row(out, i, outrow);
        double oneterm;
        int    bestindex  = 0;
        double bestscore  = 0;
        gsl_vector_set(outrow, 0, 1);
        for(j=0; j < params->size2+1; j ++){
            if (j==0)
                oneterm = 1;
            else{
                Apop_matrix_col(params, j-1, p);
                gsl_blas_ddot(observation, p, &oneterm);
            }
            if (oneterm > bestscore){
                bestindex = j;
                bestscore = oneterm;
            }
            gsl_vector_set(outrow, j, exp(oneterm));
        }
        double total = apop_sum(outrow);
        gsl_vector_scale(outrow, 1/total);
        apop_data_set(out, i, -1, bestindex);
    }
    apop_name_stack(out->names, m->parameters->names, 'c');
    return out;
}


static double val;
static double onerow(double in){ return in >= val; }

/*
This is just a for loop that runs a probit on each row.
*/
static double multiprobit_log_likelihood(apop_data *d, apop_model *p){
  apop_assert(p->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
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
 The first column of the data matrix this model expects a number
 indicating the preferred category; the remaining columns are values of
 the independent variables. Thus, the model will return N-1 columns of
 parameters, where N is the number of categories chosen.

\ingroup models
*/
apop_model apop_multinomial_probit = {"Multinomial probit",
     .log_likelihood = multiprobit_log_likelihood, .prep = multiprobit_prep};


static size_t find_index(double in, double *m, size_t max){
  size_t i = 0;
    while (in !=m[i] && i<max) i++;
    return i;
}

/**

  The likelihood of choosing item $j$ is:
  $e^{x\beta_j}/ (\sum_i{e^{x\beta_i}})$

  so the log likelihood is 
  $x\beta_j  - ln(\sum_i{e^{x\beta_i}})$

  A nice trick: let $y_i = x\beta_i$.
$ln(\sum_i{e^{x\beta_i}}) = max(y_i) - ln(\sum_i{e^{y_i - max(y_i)}})$
the elements of the sum are all now exp(something negative), so we don't
have to worry about overflow, and if there's underflow, then that term
must not have been very important. [This trick is attributed to Tom
Minka, who implemented it in his Lightspeed Matlab toolkit.]
*/
static double multilogit_log_likelihood(apop_data *d, apop_model *p){
  apop_assert(p->parameters, 0, 0, 's', "You asked me to evaluate an un-parametrized model.");
  size_t i, index, choicect = p->parameters->matrix->size2;

  //Find X\beta_i for each row of X and each column of \beta.
    apop_data  *xbeta = apop_dot(d, p->parameters, 0, 0);
    long double ll    = 0;

    //get the $x\beta_j$ numerator for the appropriate choice:
    for(i=0; i < d->vector->size; i++){
        index   = find_index(gsl_vector_get(d->vector, i), p->more, choicect);
        if (index< choicect)  //otherwise it's beta_0, which is fixed at zero.
            ll += apop_data_get(xbeta, i, index);
    }

    //Get the denominator, using the subtract-the-max trick above.
    //Don't forget the implicit beta_0, fixed at zero (so exp(beta_0)=1).
    for(i=0; i < xbeta->matrix->size1; i++){
        APOP_ROW(xbeta, i, thisrow);
        double max = gsl_vector_max(thisrow);
        gsl_vector_add_constant(thisrow, -max);
        apop_vector_exp(thisrow);
        ll -= max + log(apop_vector_sum(thisrow) +exp(-max)) ;
    }
    apop_data_free(xbeta);
	return ll;
}

/** The Logit model.
 The first column of the data matrix this model expects a number
 indicating the preferred category; the remaining columns are values of
 the independent variables. Thus, the model will return N-1 columns of
 parameters, where N is the number of categories chosen.

\ingroup models
*/
apop_model apop_logit = {"Logit",
     .log_likelihood = multilogit_log_likelihood, .expected_value=multilogit_expected, .prep = multiprobit_prep};
