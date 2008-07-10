/** \file model.h 
 
 Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
 */
#ifndef __apop_models_h__
#define __apop_models_h__

#include "types.h"
#include "regression.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#undef __BEGIN_DECLS    /* extern "C" stuff cut 'n' pasted from the GSL. */
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

extern apop_model apop_beta;
extern apop_model apop_bernoulli;
extern apop_model apop_binomial; //on hiatus.
extern apop_model apop_exponential;
extern apop_model apop_gamma;
extern apop_model apop_gaussian;//synonym for apop_normal
//extern apop_model apop_GLS;
extern apop_model apop_histogram;
extern apop_model apop_improper_uniform;
extern apop_model apop_iv;
extern apop_model apop_kernel_density;
extern apop_model apop_logit;
extern apop_model apop_lognormal;
extern apop_model apop_multinomial_probit;
extern apop_model apop_multivariate_normal;
extern apop_model apop_normal;
extern apop_model apop_ols;
extern apop_model apop_poisson;
extern apop_model apop_probit;
extern apop_model apop_uniform;
extern apop_model apop_waring;
extern apop_model apop_wls;
extern apop_model apop_yule;
extern apop_model apop_zipf;

#define apop_OLS apop_ols
#define apop_WLS apop_wls
#define apop_IV apop_iv



/////////Settings


apop_ls_settings * apop_ls_settings_alloc(apop_data *data);
void * apop_ls_settings_copy(apop_ls_settings *in);
void apop_ls_settings_free(apop_ls_settings *in);

typedef struct {
    int want_cov;
    void *copy;
    void *free;
} apop_normal_settings;

apop_normal_settings *apop_normal_settings_alloc(int want_cov);
apop_normal_settings *apop_normal_settings_copy(apop_normal_settings *in);
void apop_normal_settings_free(apop_normal_settings *in);

/** This is serious overkill for a single character of data---we could
 simply check for the presence of this struct and be done with it.
 But this allows for future expansion if so desired. */
typedef struct {
    char rank_data;
    void *copy;
    void *free;
}apop_rank_settings;

//in apop_exponential.c
apop_rank_settings *apop_rank_settings_alloc(void *ignoreme);
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



apop_model *apop_model_set_parameters(apop_model in, ...);
void apop_histogram_plot(apop_model *in, char *outfile);
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

__END_DECLS
#endif

/* \defgroup mle  Maximum likelihood estimation

Most of the action with regards to maximum likelihood estimation is in
the function \ref apop_maximum_likelihood and the \ref models "model objects".

The likelihood objects describe anything which one would want to fit
with an MLE. Usually this involves finding the most likely parameters
for a distribution, but this can also involve more elaborate models such
as the \ref apop_probit. 

Because the model is often a probability distribution, the apop_model
object is also Apophenia's means of describing distributions. E.g.,
the PDF of the Waring distribution at the data given the parameters is
exp(-apop_waring.log_likelihood(beta, data)). Where at all possible,
there are also random number generators for the distributions, e.g.
apop_waring.rng(r, beta), where \c r  is an allocated and initialized gsl_rng.


<b>example</b><br>
Here is a simple example; see also \ref mle for other examples.


\code
apop_estimate   * waring_parameters;
double          starting_pt[2] = {3, 0};
double          likelihood;
apop_ep params;
        params.starting_pt	= starting_pt;
        params.method          	= 1;
        params.step_size       	= 1e-2;
        params.tolerance       	= 1e-3;
        params.verbose         	= 0;
	waring_parameters      	= apop_maximum_likelihood(data, apop_waring, params);
printf("Your most likely waring parameters are %g and %g, with likelihood %g",
                        gsl_vector_get(waring_parameter->parameters, 0) gsl_vector_get(waring_parameter->parameters, 1), likelihood);
\endcode

\section vuong Comparing models
The distribution objects make it very easy to test competing models.
Vuong (1989) (<a
href="http://links.jstor.org/sici?sici=0012-9682%28198903%2957%3A2%3C307%3ALRTFMS%3E2.0.CO%3B2-J">Jstor
link</a>) shows that in most cases, the log likelihood ratio is asymptotically normally
distributed, so it is reasonable to apply the following paired t-test:

\code
//A function to produce two ML estimates and compare the output. 
//I had the Waring and Gamma distributions in mind when I wrote this (thus the starting points),
//e.g. call with: compare_two_distributions(data, apop_waring, apop_gamma);
//In the field, you would probably pass in est1 and est2 instead of calculating them here.
void compare_two_distributions(gsl_matrix *data, apop_model d1, apop_model d2){
gsl_vector      *llone, *lltwo;
double          mean, t_stat,
                starting_pt_w[2]= {2.12, .40},
                //starting_pt_w[2]= {2.9795, .01},
                starting_pt_g[2] = {0.12, .40};
apop_estimate   *est1, *est2;

        printf("\n%s estimate:", d1.name);
        est1    = apop_maximum_likelihood(data, NULL, d1, starting_pt_w, .001, 0);
        apop_estimate_print(est1);
        printf("\n%s estimate:", d2.name);
        est2    = apop_maximum_likelihood(data, NULL, d2, starting_pt_g, .001, 0);
        apop_estimate_print(est2);

        //Produce two vectors giving the likelihood of each row in the data set under the two models.
        apop_make_likelihood_vector(data, &lltwo, d1, est1->parameters);
        apop_make_likelihood_vector(data, &llone, d2, est2->parameters);

        gsl_vector_scale(lltwo, -1);
        gsl_vector_add(llone, lltwo);
        mean    = apop_mean(llone);
        t_stat  = apop_paired_t_test(llone,lltwo);
        if (mean > 0)
           printf("The %s is a better fit than the %s with %g%% certainty.\n", d1.name, d2.name, t_stat*100);
        else
           printf("The %s is a better fit than the %s with %g%% certainty.\n", d2.name, d1.name, t_stat*100);
}
\endcode
\todo revise this code. It's about six months behind.
*/
