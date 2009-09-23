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
extern apop_model apop_dirichlet;
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
extern apop_model apop_pmf;
extern apop_model apop_poisson;
extern apop_model apop_probit;
extern apop_model apop_t_distribution;
extern apop_model apop_uniform;
extern apop_model apop_waring;
extern apop_model apop_wishart;
extern apop_model apop_wls;
extern apop_model apop_yule;
extern apop_model apop_zipf;

/** Alias for the \ref apop_normal distribution, qv.
\hideinitializer */
#define apop_gaussian apop_normal
#define apop_OLS apop_ols
#define apop_PMF apop_pmf
#define apop_F_distribution apop_f_distribution
#define apop_WLS apop_wls
#define apop_IV apop_iv


void apop_model_free (apop_model * free_me);
void apop_model_show (apop_model * print_me);
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

#define apop_model_set_parameters(in, ...) apop_model_set_parameters_base((in), (double []) {__VA_ARGS__})
apop_model *apop_model_set_parameters_base(apop_model in, double ap[]);

#ifdef	__cplusplus
}
#endif
#endif
