/** \file apop_binomial.c 
 
  The binomial distribution as an \c apop_model.*/
/*Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "mapply.h"
#include "internal.h"
#include "likelihoods.h"

static double binomial_log_likelihood(apop_data*, apop_model*);

static void get_hits_and_misses(const apop_data *data, char method, double *hitcount, double *misscount){
    if (method == 't'){
        size_t        i, j;
        *hitcount=0, *misscount=0;
        if (data->vector)
            for(i=0; i < data->vector->size; i ++){
                if (gsl_vector_get(data->vector, i))
                    (*hitcount)    ++;
                else
                    (*misscount)   ++;
            }
        if (data->matrix)
            for(i=0; i < data->matrix->size1; i ++)
                for(j=0; j < data->matrix->size2; j ++){
                    if (gsl_matrix_get(data->matrix, i, j))
                        (*hitcount)    ++;
                    else
                        (*misscount)   ++;
                }
    } else {
        APOP_COL(data, 0, misses);
        APOP_COL(data, 1, hits);
        *hitcount = apop_vector_sum(hits);
        *misscount = apop_vector_sum(misses);
    }
}

static void make_covar(apop_model *est){
    int size = est->parameters->vector->size;
    //the trick where we turn the params into a p-vector
    double * pv = est->parameters->vector->data;
    int n = pv[0];
    pv[0] = 1 - (apop_sum(est->parameters->vector)-n);

    apop_data *cov = apop_data_alloc(0, size, size);
    for (int i=0; i < size; i++){
        double p = apop_data_get(est->parameters, i, -1);
        apop_data_set(cov, i, i, n * p *(1-p));
        for (int j=i+1; j < size; j++){
            double pj = apop_data_get(est->parameters, j, -1);
            apop_data_set(cov, i, j, -n*p*pj);
            apop_data_set(cov, j, i, -n*p*pj);
        }
    }
    apop_data_add_page(est->parameters, cov, "Covariance");
    pv[0]=n;
}

static apop_model * binomial_estimate(apop_data * data,  apop_model *est){
  Nullcheck_d(data)
  double hitcount, misscount;
  char method = apop_settings_get_group(est, "apop_rank") ? 'b' : 't';
    get_hits_and_misses(data, method, &hitcount, &misscount);   
    int n = hitcount + misscount;
    apop_name_add(est->parameters->names, "n", 'r');
    apop_name_add(est->parameters->names, "p", 'r');
    apop_data_set(est->parameters, 0, -1, n);
    apop_data_set(est->parameters, 1, -1, hitcount/(hitcount + misscount));
    est->dsize	        = (method == 'b') ? 2 : n;
    apop_data_add_named_elmt(est->info, "log likelihood", binomial_log_likelihood(data, est));
    make_covar(est);
    return est;
}

static double binomial_log_likelihood(apop_data *d, apop_model *params){
  Nullcheck_m(params) Nullcheck_p(params) Nullcheck_d(d)
  double	  n       = apop_data_get(params->parameters, 0, -1),
              p       = apop_data_get(params->parameters, 1, -1);
  double hitcount, misscount, ll = 0;
  char method = apop_settings_get_group(params, "apop_rank") ? 'b' : 't';
    if (method == 't'){
        get_hits_and_misses(d, method, &hitcount, &misscount);
        return log(gsl_ran_binomial_pdf(hitcount, p, n));
    } else {
        for (size_t i=0; i< d->matrix->size1; i++){
            hitcount = gsl_matrix_get(d->matrix, i, 1);
            ll += log(gsl_ran_binomial_pdf(hitcount, p, n));
        }
        return ll;
    }
}

static double binomial_cdf(apop_data *d, apop_model *est){
  Nullcheck_m(est) Nullcheck_p(est) Nullcheck_d(d)
  double hitcount, misscount, psum = 0;
  char method = apop_settings_get_group(est, "apop_rank") ? 'b' : 't';
    get_hits_and_misses(d, method, &hitcount, &misscount);   
    double n = gsl_vector_get(est->parameters->vector, 0);
    double p = gsl_vector_get(est->parameters->vector, 1);
    for (int i=0; i< hitcount; i++)
        psum += gsl_ran_binomial_pdf(hitcount, p, n);
    return psum;
}

static double multinomial_constraint(apop_data *data, apop_model *b){
  //constraint is that 0 < all elmts 
    return apop_linear_constraint(b->parameters->vector, .margin = 1e-3);
}

static void binomial_rng(double *out, gsl_rng *r, apop_model* est){
  double n = gsl_vector_get(est->parameters->vector, 0);
  double p = gsl_vector_get(est->parameters->vector, 1);
  char method = apop_settings_get_group(est, "apop_rank") ? 'b' : 't';
    if (method == 'b'){
        out[1] =  gsl_ran_binomial(r, n ,est->parameters->vector->data[0]); 
        out[0] =  n - out[1];
    } else
        for (int i=0; i < n; i++)
            out[i] = (gsl_rng_uniform(r) <= p) ? 1 : 0; //one Bernoulli draw.
}


static double is_nonzero(double in){return in != 0;}
static double sum_vector_nonzeros(gsl_vector *in){return apop_vector_map_sum(in, is_nonzero); }

static gsl_vector * get_multinomial_hitcount(const apop_data *data, char method){
    size_t        i, j;
    gsl_vector *out;
    if (method == 't'){
        out = gsl_vector_alloc(1+GSL_MAX(data->vector ? gsl_vector_max(data->vector) : 0,
                                       data->matrix ? gsl_matrix_max(data->matrix) : 0));
        if (data->vector)
            for(i=0; i < data->vector->size; i ++)
                (*gsl_vector_ptr(out, apop_data_get(data, i, -1)))++;
        if (data->matrix)
            for(i=0; i < data->matrix->size1; i ++)
                for(j=0; j < data->matrix->size2; j ++)
                    (*gsl_vector_ptr(out, apop_data_get(data, i, j)))++;
    } else {//just count nozeros in each column
        apop_data *outd = apop_map((apop_data *)data, .fn_v=sum_vector_nonzeros, .part='c');
        out = outd->vector;
        outd->vector=NULL;
        apop_data_free(outd);
    }
    return out;
}

static double multinomial_log_likelihood(apop_data *d, apop_model *params){
  apop_assert(params->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
    double *pv = params->parameters->vector->data;
    size_t size = params->parameters->vector->size;
    char method = apop_settings_get_group(params, "apop_rank") ? 'b' : 't';

    //The GSL wants our hit count in an int*.
    gsl_vector *hits = get_multinomial_hitcount(d, method);
    unsigned int *hv = malloc(hits->size * sizeof(unsigned int));
    for(size_t i=0; i < hits->size; i ++)
        hv[i] = hits->data[i];
    gsl_vector_free(hits);

    double n = pv[0]; //making the params a p-vector. Put n back at the end.
    pv[0] = 1 - (apop_sum(params->parameters->vector)-n);
    double out =  gsl_ran_multinomial_lnpdf(size, pv, hv);

    pv[0]=n;
    free(hv);
    return out;
}

static apop_model * multinomial_estimate(apop_data * data,  apop_model *est){
    apop_assert(data,  0, 0,'s', "You asked me to estimate the parameters of a model but sent NULL data.");
    char method = apop_settings_get_group(est, "apop_rank") ? 'b' : 't';
    gsl_vector * count = get_multinomial_hitcount(data, method);
    //int n = apop_sum(count); //potential double-to-int precision issues.
    int n = 0;
    for (int i=0; i< count->size; i++)
        n += gsl_vector_get(count, i);
    apop_vector_normalize(count);
    gsl_vector_set(count, 0, n);
    est->parameters=apop_data_alloc(0,0,0);
    est->parameters->vector = count;
    apop_name_add(est->parameters->names, "n", 'c');
    char name[100];
    for(int i=1; i < count->size; i ++){
        sprintf(name, "p%i", i);
        apop_name_add(est->parameters->names, name, 'c');
    }
    est->dsize	        = (method == 'b') ? count->size : n;
    make_covar(est);
    apop_data_add_named_elmt(est->info, "log likelihood", multinomial_log_likelihood(data, est));
    return est;
}

static void multinomial_rng(double *out, gsl_rng *r, apop_model* est){
    //After the intro, cut/pasted/modded from the GSL. Copyright them.
    Nullcheck_pv(est);
    char method = apop_settings_get_group(est, "apop_rank") ? 'b' : 't';
    double * p = est->parameters->vector->data;
    size_t k = est->parameters->vector->size;
    //the trick where we turn the params into a p-vector
    int N = p[0];
    p[0] = 1 - (apop_sum(est->parameters->vector)-N);
    double sum_p = 0.0;
    int draw, ctr = 0;
    unsigned int sum_n = 0;

  /* p[i] may contain non-negative weights that do not sum to 1.0.
   * Even a probability distribution will not exactly sum to 1.0
   * due to rounding errors. 
   */
    double norm = apop_sum(est->parameters->vector);

    for (int i = 0; i < k; i++) {
        if (p[i] > 0.0)
            draw = gsl_ran_binomial (r, p[i] / (norm - sum_p), N - sum_n);
        else
            draw = 0;
        if (method == 'b')
            out[i] = draw;
        else
            for (int k=0; k < draw; k++)
                out[ctr++] = i;
        sum_p += p[i];
        sum_n += out[i];
    }
    p[0]=N;
}

apop_model apop_binomial = {"Binomial distribution", 2,0,0,
	.estimate = binomial_estimate, .log_likelihood = binomial_log_likelihood, 
   .constraint = multinomial_constraint, .draw = binomial_rng, .cdf =binomial_cdf};

apop_model apop_multinomial = {"Multinomial distribution", -1,0,0,
	.estimate = multinomial_estimate, .log_likelihood = multinomial_log_likelihood, 
   .constraint = multinomial_constraint, .draw = multinomial_rng};
