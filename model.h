/** \file model.h  */
/* Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#ifndef __apop_models_h__
#define __apop_models_h__

#include "types.h"
#include "variadic.h"
#include "stats.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#ifdef	__cplusplus
extern "C" {
#endif

extern apop_model apop_beta;
extern apop_model apop_bernoulli;
extern apop_model apop_binomial;
extern apop_model apop_chi_squared;
extern apop_model apop_exponential;
extern apop_model apop_f_distribution;
extern apop_model apop_gamma;
extern apop_model apop_histogram;
extern apop_model apop_improper_uniform;
extern apop_model apop_iv;
extern apop_model apop_kernel_density;
extern apop_model apop_logit;
extern apop_model apop_lognormal;
extern apop_model apop_multinomial;
extern apop_model apop_multinomial_probit;
extern apop_model apop_multivariate_normal;
extern apop_model apop_normal;
extern apop_model apop_ols;
extern apop_model apop_poisson;
extern apop_model apop_probit;
extern apop_model apop_t_distribution;
extern apop_model apop_uniform;
extern apop_model apop_waring;
extern apop_model apop_wishart;
extern apop_model apop_wls;
extern apop_model apop_yule;
extern apop_model apop_zipf;

#define apop_gaussian apop_normal
#define apop_OLS apop_ols
#define apop_F_distribution apop_f_distribution
#define apop_WLS apop_wls
#define apop_IV apop_iv


/////////Settings

/** Settings for least-squares type models */
typedef struct {
    int destroy_data;
    gsl_vector *weights;
    apop_data *instruments;
    int want_cov;
    int want_expected_value;
} apop_ls_settings;

apop_ls_settings * apop_ls_settings_alloc(apop_data *data);
apop_ls_settings * apop_ls_settings_init(apop_ls_settings);
void * apop_ls_settings_copy(apop_ls_settings *in);
void apop_ls_settings_free(apop_ls_settings *in);

// Find apop_category_settings routines in apop_probit.c
/** for dependent-category models, send in this settings struct to
  specify which column is the dependent variable. 

  If you don't use it, these models will assume that the vector or first
  numeric column is already a coherent set of factors, but by sending
  this in, those functions have a little more information, such as names
  to use in the output.

  See also the \ref apop_category_settings_alloc function.
 */
typedef struct {
    apop_data *factors; 
    char source_type; 
    char source_column; 
    apop_data *source_data;
} apop_category_settings;

apop_category_settings *apop_category_settings_alloc(apop_data *d, int source_column, char source_type);
apop_category_settings *apop_category_settings_copy(apop_category_settings *in);
void apop_category_settings_free(apop_category_settings *in);

/** If this settings group is present, models that can take rank data
  will read the input data as such.  Allocation is thus very simple, e.g.
  \code
  Apop_settings_group_add(your_model, apop_rank, NULL);
  \endcode */
typedef struct {
    char rank_data;
} apop_rank_settings;

//in apop_exponential.c
apop_rank_settings *apop_rank_settings_alloc(void *ignoreme);
apop_rank_settings *apop_rank_settings_init(apop_rank_settings);
void apop_rank_settings_free(apop_rank_settings *in);
void *apop_rank_settings_copy(apop_rank_settings *in);

#include <gsl/gsl_histogram.h>
typedef struct{
    gsl_histogram       *pdf;
    gsl_histogram_pdf   *cdf;
    apop_model          *histobase;
    apop_model          *kernelbase;
} apop_histogram_settings;

#define apop_kernel_density_settings apop_histogram_settings

apop_histogram_settings *apop_histogram_settings_alloc(apop_data *data, int bins);
void  apop_histogram_settings_free(apop_histogram_settings *in);
void * apop_histogram_settings_copy(apop_histogram_settings *in);


#define apop_model_set_parameters(in, ...) apop_model_set_parameters_base((in), (double []) {__VA_ARGS__})
apop_model *apop_model_set_parameters_base(apop_model in, double ap[]);
apop_histogram_settings *apop_kernel_density_settings_alloc(apop_data *data, 
        apop_model *histobase, apop_model *kernelbase, void (*set_params)(double, apop_model*));

#define apop_kernel_density_settings_copy apop_histogram_settings_copy
#define apop_kernel_density_settings_free apop_histogram_settings_free

apop_model * apop_model_copy(apop_model in); //in apop_model.c
apop_model * apop_model_clear(apop_data * data, apop_model *model);

apop_model * apop_estimate(apop_data *d, apop_model m);
void apop_score(apop_data *d, gsl_vector *out, apop_model *m);
double apop_log_likelihood(apop_data *d, apop_model *m);
double apop_p(apop_data *d, apop_model *m);
void apop_draw(double *out, gsl_rng *r, apop_model *m);
void apop_model_prep(apop_data *d, apop_model *m);
apop_data * apop_expected_value(apop_data *d, apop_model *m);

//in apop_beta.c
apop_model *apop_beta_from_mean_var(double m, double v);

#ifdef	__cplusplus
}
#endif
#endif
