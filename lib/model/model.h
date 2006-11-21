#ifndef __apop_models_h__
#define __apop_models_h__

#include <apophenia/types.h>
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

extern apop_model apop_exponential;
extern apop_model apop_exponential_rank;
extern apop_model apop_gamma;
extern apop_model apop_gamma_rank;
extern apop_model apop_gaussian;//synonym for apop_normal
extern apop_model apop_GLS;
extern apop_model apop_logit;
extern apop_model apop_normal;
extern apop_model apop_OLS;
extern apop_model apop_WLS;
extern apop_model apop_poisson;
extern apop_model apop_probit;
extern apop_model apop_waring;
extern apop_model apop_waring_rank;
extern apop_model apop_yule;
extern apop_model apop_yule_rank;
extern apop_model apop_zipf;
extern apop_model apop_zipf_rank;


//This is in asst.c.
double apop_generalized_harmonic(int N, double s);

__END_DECLS
#endif

/** \defgroup mle  Maximum likelihood estimation

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
