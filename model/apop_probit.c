/** \file */ //Tell Doxygen to process this.

/* Copyright (c) 2005--2008, 2010 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "mapply.h"
#include "internal.h"
#include "likelihoods.h"

/////////  Part II: plain old probit

static apop_data *get_category_table(apop_data *d){
    int first_col = d->vector ? -1 : 0;
    apop_data *out = apop_data_get_page(d, "<categories");
    if (!out) 
        out = apop_data_to_factors(d, .intype='d', .incol=first_col, .outcol=first_col);
    return out;
}

static void probit_prep(apop_data *d, apop_model *m){
  int       count;
  apop_data *factor_list;
    apop_ols.prep(d, m);//also runs the default apop_model_clear.
    factor_list = get_category_table(d);
    count = factor_list->textsize[0];
    m->parameters = apop_data_alloc(0, d->matrix->size2, count-1);
    apop_name_stack(m->parameters->names, d->names, 'r', 'c');
    for (int i=1; i< count; i++) 
        apop_name_add(m->parameters->names, factor_list->text[i][0], 'c');
    gsl_matrix_set_all(m->parameters->matrix, 1);
    char *tmp = strdup(m->name);
    snprintf(m->name, 100, "%s with %s as numeraire", tmp, factor_list->text[0][0]);
    free(tmp);
}

double biprobit_ll_row(apop_data *r){
    long double n = gsl_cdf_gaussian_P(-gsl_matrix_get(r->matrix, 0, 0),1);
    n = n ? n : 1e-10; //prevent -inf in the next step.
    n = n<1 ? n : 1-1e-10; 
    return r->vector->data[0] ?  log(1-n): log(n);
}

//The case where outcome is a single zero/one option.
static double biprobit_log_likelihood(apop_data *d, apop_model *p){
  Nullcheck_m(p); Nullcheck_p(p);
    apop_data *betadotx = apop_dot(d, p->parameters); 
    betadotx->vector = d->vector;
    double total_prob = apop_map_sum(betadotx, .fn_r=biprobit_ll_row);
    betadotx->vector = NULL;
    apop_data_free(betadotx);
	return total_prob;
}

static void probit_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
  Nullcheck_m(p); Nullcheck_p(p);
    gsl_vector *val_vector = get_category_table(p->data)->vector;
    if (val_vector->size!=2){
        gsl_vector * numeric_default = apop_numerical_gradient(d, p);
        gsl_vector_memcpy(gradient, numeric_default);
        gsl_vector_free(numeric_default);
    }
  long double	cdf, betax, deriv_base;
  apop_data *betadotx = apop_dot(d, p->parameters); 
    gsl_vector_set_all(gradient,0);
    for (size_t i=0; i< d->matrix->size1; i++){
        betax            = apop_data_get(betadotx, i, 0);
        cdf              = gsl_cdf_gaussian_P(-betax, 1);
        cdf = cdf ? cdf : 1e-10; //prevent -inf in the next step.
        cdf = cdf<1 ? cdf : 1-1e-10; 
        if (apop_data_get(d, i, -1))
            deriv_base      = gsl_ran_gaussian_pdf(-betax, 1) /(1-cdf);
        else
            deriv_base      = -gsl_ran_gaussian_pdf(-betax, 1) / cdf;
        for (size_t j=0; j< d->matrix->size2; j++)
            apop_vector_increment(gradient, j, apop_data_get(d, i, j) * deriv_base);
	}
	apop_data_free(betadotx);
}

/////////  Part III: Multinomial Logit (plain logit is a special case)

static apop_data *multilogit_expected(apop_data *in, apop_model *m){
    Nullcheck_m(m); Nullcheck_p(m);
    apop_prep(in, m);
    gsl_matrix *params = m->parameters->matrix;
    apop_data *out = apop_data_alloc(in->matrix->size1, in->matrix->size1, params->size2+1);
    for (size_t i=0; i < in->matrix->size1; i ++){
        Apop_row(in, i, observation);
        Apop_row(out, i, outrow);
        double oneterm;
        int    bestindex  = 0;
        double bestscore  = 0;
        gsl_vector_set(outrow, 0, 1);
        for (size_t j=0; j < params->size2+1; j ++){
            if (j == 0)
                oneterm = 0;
            else {
                Apop_matrix_col(params, j-1, p);
                gsl_blas_ddot(observation, p, &oneterm);
                gsl_vector_set(outrow, j, exp(oneterm));
            }
            if (oneterm > bestscore){
                bestindex = j;
                bestscore = oneterm;
            }
        }
        double total = apop_sum(outrow);
        gsl_vector_scale(outrow, 1/total);
        apop_data_set(out, i, -1, bestindex);
    }
    apop_data *factor_list = get_category_table(m->data);
    apop_name_add(out->names, factor_list->text[0][0], 'c');
    apop_name_stack(out->names, m->parameters->names, 'c');
    return out;
}


static double val;
static double unordered(double in){ return in == val; }
//static double ordered(double in){ return in >= val; }

// This is just a for loop that runs a probit on each row.
static double multiprobit_log_likelihood(apop_data *d, apop_model *p){
  Nullcheck_m(p); Nullcheck_p(p);
    gsl_vector *val_vector = get_category_table(p->data)->vector;
    if (val_vector->size==2)
        return biprobit_log_likelihood(d, p);
    //else, multinomial loop
    static apop_model *spare_probit = NULL;
    if (!spare_probit){
        spare_probit = apop_model_copy(apop_probit);
        spare_probit->parameters = apop_data_alloc(0,0,0);
    }
    static apop_data *working_data = NULL;
    if (!working_data)
        working_data = apop_data_alloc(0,0,0);

    working_data->matrix = d->matrix;
    gsl_vector *original_outcome = d->vector;
    double ll    = 0;
    double *vals = val_vector->data;
    for(size_t i=0; i < p->parameters->matrix->size2; i++){
        APOP_COL(p->parameters, i, param);
        val = vals[i];
        working_data->vector = apop_vector_map(original_outcome, unordered);
        //gsl_matrix_set_col(spare_probit->parameters->matrix, 0, param);
        spare_probit->parameters->matrix = apop_vector_to_matrix(param);
        ll  += apop_log_likelihood(working_data, spare_probit);
        gsl_vector_free(working_data->vector); //yup. It's inefficient.
        gsl_matrix_free(spare_probit->parameters->matrix);
    }
	return ll;
}

static size_t find_index(double in, double *m, size_t max){
  size_t i = 0;
    while (in !=m[i] && i<max) i++;
    return i;
}

static double multilogit_log_likelihood(apop_data *d, apop_model *p){
  Nullcheck_m(p); Nullcheck_p(p);
  size_t i, j, index, choicect = p->parameters->matrix->size2;
  double* factor_list = get_category_table(p->data)->vector->data;

    //Find X\beta_i for each row of X and each column of \beta.
    apop_data  *xbeta = apop_dot(d, p->parameters, 0, 0);

    //get the $x\beta_j$ numerator for the appropriate choice:
    long double ll    = 0;
    for(i=0; i < d->vector->size; i++){
        index   = find_index(gsl_vector_get(d->vector, i), factor_list, choicect);
        ll += (index==0) ? 0 : apop_data_get(xbeta, i, index-1);
    }

    //Get the denominator, using the subtract-the-max trick mentioned in the documentation.
    //Don't forget the implicit beta_0, fixed at zero (so we need to add exp(0-max)).
    for(j=0; j < xbeta->matrix->size1; j++){
        APOP_ROW(xbeta, j, thisrow);
        double max = gsl_vector_max(thisrow);
        gsl_vector_add_constant(thisrow, -max);
        apop_vector_exp(thisrow);
        ll -= max + log(apop_vector_sum(thisrow) +exp(-max)) ;
    }
    apop_data_free(xbeta);
	return ll;
}

apop_model *logit_estimate(apop_data *d, apop_model *m){
  apop_model *out = apop_maximum_likelihood(d, m);

    //That's the whole estimation. But now we need to add a row of zeros
    //for the numeraire. This is just tedious matrix-shunting.
    apop_data *p = out->parameters;
    apop_data *zeros = apop_data_calloc(0, p->matrix->size1, 1);
    apop_data *newparams = apop_data_stack(zeros, p, 'c');
    apop_name_stack(newparams->names, p->names, 'r');

    apop_data_free(p);
    apop_data_free(zeros);
    out->parameters = newparams;
    return out;
}

apop_model apop_logit = {"Logit", .log_likelihood = multilogit_log_likelihood, .dsize=-1,
     .predict=multilogit_expected, .prep = probit_prep};

apop_model apop_probit = {"Probit", .log_likelihood = multiprobit_log_likelihood, .dsize=-1,
    .score = probit_dlog_likelihood, .prep = probit_prep};
