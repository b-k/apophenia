#ifndef __apop_models_h__
#define __apop_models_h__

#include <apophenia/estimate.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/** This is an object to describe a model whose parameters are to be
estimated. It would primarily be used for maximum likelihood estimation,
but is intended to have anything else you would want a probability
distribution to have too, like a random number generator.  */
typedef struct apop_model{
	/** The model name. You have 100 characters. */
	char	name[101]; 

	/** The number of parameters. If this is -1, it will be
	 dynamically set to the size of the given data set minus one.
	  */
	int	parameter_ct;

	/** the estimator fn, which is all most users will care about.*/ 
	apop_estimate *	(*estimate)(gsl_matrix * data, apop_inventory *uses, void *parameters);

	/** the likelihood fn given data*/ 
	double 	(*log_likelihood)(const gsl_vector *beta, void *d);

	/** the derivative of the likelihood fn */
	void 	(*dlog_likelihood)(const gsl_vector *beta, void *d, gsl_vector *gradient);

	/** Do both of the above at once. Can be NULL if it'd just call them separately. */
	void 	(*fdf)( const gsl_vector *beta, void *d, double *f, gsl_vector *df);

//	/** The constraint count */
//	int 	constraint_ct;
//
//	/** The vector of constraints */
//	double *(* constraint) (const gsl_vector *beta, void *d, gsl_vector *inside_constraint);
	
	/** a random number generator. */
	double (*rng)(gsl_rng* r, double *a);
	//to add: 
	//the 2nd derivative of the likelihood fn
	//the MLE
} apop_model;


extern apop_model apop_exponential;
extern apop_model apop_gamma;
extern apop_model apop_gamma_rank;
extern apop_model apop_gaussian;//synonym for apop_normal
extern apop_model apop_normal;
extern apop_model apop_probit;
extern apop_model apop_waring;
extern apop_model apop_yule;
extern apop_model apop_zipf;



//This is in asst.c.
double apop_generalized_harmonic(int N, double s);

#endif

/** \defgroup likelihood_fns  Likelihood fns 

The \ref apop_model objects are to Apophenia as the 'model' object is to
most other statistics packages: it is a summary of the author's claims
about the real world. Hand an apop_model plus a data set 
to the \ref apop_maximum_likelihood function, and that function returns
the most likely parameters.

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
waring_parameters      = apop_maximum_likelihood(data, apop_waring, starting_pt, 1e-4, 0);
printf("Your most likely waring parameters are %g and %g, with likelihood %g",
                        gsl_vector_get(waring_parameter->parameters, 0) gsl_vector_get(waring_parameter->parameters, 1), likelihood);
\endcode

\section write_likelihoods Writing your own
Writing apop_model objects is easy:

\li Write a likelihood function. Its header will look like this:
\code
double apop_new_log_likelihood(const gsl_vector *beta, void *d)
\endcode 
where \c beta will be the parameters to be maximized, and \c
d is the fixed parameters---the data. In every case currently included
with Apophenia, \c d is a \c gsl_matrix, but you do not have to conform
to that. This function will return the value of the log likelihood function at the given parameters.
\li Write the object. In your header file, include 
\code
apop_model apop_new_likelihood = {"The Me distribution", number_of_parameters, apop_new_log_likelihood, NULL, NULL, NULL};
\endcode
\c number_of_parameters is probably a positive integer like \c 2, but
it is often (the number of columns in your data set) -1, in which case,
set \c number_of_parameters to \c -1.
\li Test. Debug. Retest.
\li (optional) Write a gradient for the log likelihood function. This
typically involves calculating a derivative by hand, which is an easy
problem in high-school calculus. The function's header will look like: 
\code
void apop_new_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient)
\endcode 
where \c beta and \c d are fixed as above, and \c gradient is a \c gsl_vector with dimension matching \c beta. 
At the end of this function, you will have to assign the appropriate derivative to every element of the gradient vector:
\code
gsl_vector_set(gradient,0, d_a);
gsl_vector_set(gradient,1, d_b);
\endcode 
Now add the resulting dlog likelihood function to your object:
\code
apop_model apop_new_likelihood = {"The Me distribution", number_of_parameters, apop_new_log_likelihood, apop_new_dlog_likelihood, NULL, NULL};
\endcode
\li Send the code to the maintainer for inclusion in future versions of Apophenia.



\ingroup mle 
*/

/** \defgroup mle  Maximum likelihood estimation
\ingroup likelihood_fns

Most of the action with regards to maximum likelihood estimation is in
the function \ref apop_maximum_likelihood and the \ref likelihood_fns "distribution objects".

The likelihood objects describe anything which one would want to fit
with an MLE. Usually this involves finding the most likely parameters
for a distribution, but this can also involve more elaborate models such
as the \ref apop_probit. 

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
*/


