/** \file settings.h */ 
/* Copyright (c) 2009-10 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#ifndef __apop_settings_h__
#define __apop_settings_h__

#include "types.h"
#include "asst.h"

#ifdef	__cplusplus
extern "C" {
#endif

    //Part I: macros and fns for getting/setting settings groups and elements

void * apop_settings_get_grp(apop_model *m, char *type, char fail);
void apop_settings_remove_group(apop_model *m, char *delme);
void apop_settings_copy_group(apop_model *outm, apop_model *inm, char *copyme);
void *apop_settings_group_alloc(apop_model *model, char *type, void *free_fn, void *copy_fn, void *the_group);

/** Retrieves a settings group from a model.  See \ref Apop_settings_get
 to just pull a single item from within the settings group.

  If it isn't found, then it returns NULL, so you can easily put it in a conditional like 
  \code 
  if (!apop_settings_get_group(m, "apop_ols")) ...
  \endcode
\hideinitializer \ingroup settings
 */
#define Apop_settings_get_group(m, type) apop_settings_get_grp(m, #type, 'c')

/** Removes a settings group from a model's list. 
\hideinitializer \ingroup settings
 */
#define Apop_settings_rm_group(m, type) apop_settings_remove_group(m, #type)

/** Add a settings group. The first two arguments (the model you are
 attaching to and the settings group name) are mandatory, and then you
 can use the \ref designated syntax to specify default values (if any).
 Returns a pointer to the newly-prepped group.
\hideinitializer \ingroup settings
 */
#define Apop_model_add_group(model, type, ...)  \
    apop_settings_group_alloc(model, #type, type ## _settings_free, type ## _settings_copy, type ##_settings_init ((type ## _settings) {__VA_ARGS__})); 

/** A convenience for your settings group init functions. 
 Gives the output item either the default value if there is one, or the value you specify. 
\hideinitializer \ingroup settings
 */
#define apop_varad_setting(in, out, name, value) (out)->name = (in).name ? (in).name : (value);

/** Retrieves a setting from a model.  See \ref Apop_settings_get_group to pull the entire group.
\hideinitializer \ingroup settings
 */
#define Apop_settings_get(model, type, setting)  \
    (((type ## _settings *) apop_settings_get_grp(model, #type, 'f'))->setting)

/** Modifies a single element of a settings group to the given value. 
\hideinitializer \ingroup settings
 */
#define Apop_settings_set(model, type, setting, data)  \
    do { type ## _settings *apop_tmp_settings = apop_settings_get_grp(model, #type, 'c');  \
    apop_assert_s(apop_tmp_settings, "You're trying to modify a setting in " \
                        #model "'s setting group of type " #type " but that model doesn't have such a group."); \
    apop_tmp_settings->setting = (data);    \
    } while (0);

#define Apop_settings_add Apop_settings_set
#define APOP_SETTINGS_ADD Apop_settings_set
#define apop_settings_set Apop_settings_set
#define APOP_SETTINGS_GET Apop_settings_get
#define apop_settings_get Apop_settings_get
#define APOP_MODEL_ADD_GROUP Apop_model_add_group
#define apop_model_add_group Apop_model_add_group
#define APOP_SETTINGS_GET_GROUP Apop_settings_get_group
#define apop_settings_get_group Apop_settings_get_group
#define APOP_SETTINGS_RM_GROUP Apop_settings_rm_group
#define apop_settings_rm_group Apop_settings_rm_group

#define Apop_settings_declarations(ysg) \
   ysg##_settings * ysg##_settings_init(ysg##_settings); \
   void * ysg##_settings_copy(ysg##_settings *); \
   void ysg##_settings_free(ysg##_settings *);

//see deprecated.h for the apop_settings_add_group

        //Part II: the details of extant settings groups.

typedef enum {
    APOP_SIMPLEX_NM     =0, /**< Nelder-Mead simplex (gradient handling rule is irrelevant) */
    APOP_SIMPLEX_NMJ    =20, /**< Nelder-Mead simplex with occasional jiggering (for when the plain N-M gets stuck in a loop) */
    APOP_CG_FR     =1,      /**<  Conjugate gradient (Fletcher-Reeves) (default) */
    APOP_CG_BFGS   =2,      /**<  Conjugate gradient (BFGS: Broyden-Fletcher-Goldfarb-Shanno) */
    APOP_CG_PR     =3,      /**<  Conjugate gradient (Polak-Ribiere) */
    APOP_SIMAN      =5,         /**<  \ref simanneal "simulated annealing" */
    APOP_RF_NEWTON  =10,        /**<  Find a root of the derivative via Newton's method */
//    APOP_RF_BROYDEN =11,        //  Find a root of the derivative via the Broyden Algorithm
    APOP_RF_HYBRID  =12,        /**<  Find a root of the derivative via the Hybrid method */
    APOP_RF_HYBRID_NOSCALE  =13 /**<  Find a root of the derivative via the Hybrid method; no internal scaling */
} apop_optimization_enum;

/** The settings for maximum likelihood estimation (including simulated annealing).
\ingroup settings */
typedef struct{
    double      *starting_pt;   /**< An array of doubles (i.e., <tt>double*</tt>) suggesting a starting point. 
                                  If NULL, use zero.  Note that if \c v is a \c gsl_vector, then 
                                  \c v->data is of the right form (provided \c v is not a slice of a matrix).*/
    apop_optimization_enum method; /**< See the  \ref apop_optimization_enum documentation for options. */
    double      step_size, /**< the initial step size. */
                tolerance, /**< the precision the minimizer uses. Only vaguely related to the precision of the actual variables. */
delta;
    int         max_iterations; /**< Ignored by simulated annealing. Other methods halt if
                                 they do this many iterations without finding an optimum. */
    int         verbose; /**<	Give status updates as we go.  This is orthogonal to the 
                                <tt>apop_opts.verbose</tt> setting. */
    char        want_cov; /**< Deprecated. Please use \ref apop_parts_wanted_settings. */
    double      dim_cycle_tolerance; /**< If zero (the default), the usual procedure.
                             If \f$>0\f$, cycle across dimensions: fix all but the first dimension at the starting
                             point, optimize only the first dim. Then fix the all but the second dim, and optimize the
                             second dim. Continue through all dims, until the log likelihood at the outset of one cycle
                             through the dimensions is within this amount of the previous cycle's log likelihood. There
                             will be at least two cycles.
                             */
//simulated annealing (also uses step_size);
    int         n_tries, use_score, iters_fixed_T;
    double      k, t_initial, mu_t, t_min ;
    gsl_rng     *rng;
    char        *trace_path; ///< See \ref trace_path
    apop_model  *parent;
} apop_mle_settings;

/** Settings for least-squares type models 
\ingroup settings */
typedef struct {
    int destroy_data; /**< If 'y', then the input data set may be normalized or otherwise mangled */
    apop_data *instruments; /**< Use for the \ref apop_iv regression, qv. */
    char want_cov; /**< Deprecated. Please use \ref apop_parts_wanted_settings. */
    char want_expected_value; /**< Deprecated. Please use \ref apop_parts_wanted_settings. */
    apop_model *input_distribution; /**< The distribution of \f$P(Y|X)\f$ is specified by the model, but the distribution of \f$X\f$ is not.  */
} apop_lm_settings;

/** The default is for the estimation routine to give some auxiliary information,
  such as a covariance matrix, predicted values, and common hypothesis tests.
  Some uses of a model depend on these items, but if they are a waste
  of time for your purposes, this settings group gives a quick way to bypass them all.

  Simply adding this settings group to your model without changing any default values---
  \code
  Apop_model_add_group(your_model, apop_parts_wanted);
  \endcode
  ---will turn off all of the auxiliary calculations covered, because the default value
  for all the switches is <tt>'n'</tt>, indicating that all elements are not wanted.

  From there, you can change some of the default <tt>'n'</tt>s to <tt>'y'</tt>s to retain some but not all auxiliary elements.  If you just want the parameters themselves and the covariance matrix:
  \code
  Apop_model_add_group(your_model, apop_parts_wanted, .covariance='y');
  \endcode

  \li Not all models support this, although the models with especially compute-intensive
  auxiliary info do (e.g., the maximum likelihood estimation system). Check the model's documentation. 

  \li Tests may depend on covariance, so <tt>.covariance='n', .tests='y'</tt> may be 
  treated as <tt>.covariance='y', .tests='y'</tt>.

*/
typedef struct {
    //init/copy/free are in apop_mle.c
    char covariance;    /*< If 'y', calculate the covariance matrix. Default 'n'. */
    char predicted;/*< If 'y', calculate the predicted values. This is typically as many
                     items as rows in your data set. Default 'n'. */
    char tests;/*< If 'y', run any hypothesis tests offered by the model's estimation routine. Default 'n'. */
    char info;/*< If 'y', add an info table with elements such as log likelihood or AIC. Default 'n'. */
} apop_parts_wanted_settings;

/** Some CDFs use random draws; some use closed-form models. 
  \ingroup settings */
typedef struct {
    int draws;  /**< For random draw methods, how many draws? Default: 10,000.*/
    gsl_rng *rng; /**< For random draw methods. See \ref autorng on the default. */
    apop_model *cdf_model; /**< For use by individual models as they see fit. Default=\c NULL. */
    int rng_owner; /**< For internal use. */
} apop_cdf_settings;

/** Settings for getting parameter models (i.e. the distribution of parameter estimates)
  \ingroup settings */
typedef struct {
    apop_model *base;
    int index;
    gsl_rng *rng;
    int draws;
    int own_rng;
} apop_pm_settings;

#include <gsl/gsl_histogram.h>
/** Settings for the histogram and kernel density structures, mostly opaque. 
  On setup, you must set the \c .data and \c .bins_in items. 
  
  \ingroup settings */
typedef struct{
    apop_data           *data;
    gsl_histogram       *pdf; /**< Where the histogram is kept */
    gsl_histogram_pdf   *cdf; /**< If you make random draws, I need a CDF aggregation of the main PDF. I keep it here. */
    apop_model          *histobase;
    apop_model          *kernelbase;
    int                 bins_in; /**< Used as input. May not equal the final number of bins (\c pdf->bins) due to the infinibins.*/
} apop_histogram_settings;

/** Settings for the \ref apop_kernel_density model. 

  \ingroup settings */
typedef struct{
    apop_data *base_data; /**< The data that will be smoothed by the KDE. */
    apop_model *base_pmf; /**< I actually need the data in a \ref apop_pmf. You can give
                            that to me explicitly, or I can wrap the .base_data in a PMF.  */
    apop_model *kernel; /**< The distribution to be centered over each data point. Default, 
                                    \ref apop_normal with std dev 1. */
    void (*set_fn)(apop_data*, apop_model*); /**< The function I will use for each data
                                                  point to center the kernel over each point.*/
    int own_pmf, own_kernel; /**< For internal use only. */
}apop_kernel_density_settings;


/** Method settings for a model to be put through Bayesian updating. 
\ingroup settings 
 */
typedef struct{
    apop_data *data;
    apop_data *starting_pt; /**< The first parameter to check in the MCMC routine */
    long int periods; /**< For how many steps should the MCMC chain run? */
    double burnin; /**< What <em>percentage</em> of the periods should be ignored
                         as initialization. That is, this is a number between zero and one. */
    int histosegments; /**< If outputting a \ref apop_histogram, how many segments should it have? */
    char method;
} apop_update_settings;

apop_update_settings *apop_update_settings_init(apop_update_settings);
#define apop_update_settings_copy NULL
#define apop_update_settings_free NULL

//Loess, including the old FORTRAN-to-C.
struct loess_struct {
	struct {
		long    n, p;
        double  *y, *x;
		double	*weights;
	} in;
	struct {
	        double  span;
	        long    degree;
	        long    normalize;
	        long    parametric[8];
	        long    drop_square[8];
	        char    *family;
	} model;
	struct {
	        char    *surface;
	        char    *statistics;
	        double  cell;
	        char    *trace_hat;
	        long    iterations;
	} control;
	struct {
		long	*parameter, *a;
		double	*xi, *vert, *vval;
	} kd_tree;
	struct {
		double	*fitted_values;
        double  *fitted_residuals;
		double  enp, s;
		double  one_delta, two_delta;
		double	*pseudovalues;
		double	trace_hat;
		double	*diagonal;
		double	*robust;
		double  *divisor;
	} out;
};

struct anova_struct {
	double	dfn;
	double	dfd;
	double  F_value;
	double  Pr_F;
};


/** The code for the loess system is based on FORTRAN code from 1988,
overhauled in 1992, linked in to Apophenia in 2009. The structure that
does all the work, then, is a \c loess_struct that you should
basically take as opaque. 

The useful settings from that struct re-appear in the \ref
apop_loess_settings struct so you can set them directly, and then the
settings init function will copy your preferences into the working struct.

The documentation for the elements is cut/pasted/modified from Cleveland,
Grosse, and Shyu.

Because it's a cut and paste job, doxygen turns this into word salad. 
If you're reading this, my apologies, and it'll stay this way for a few
more days. see the source of \c settings.h for a legible version in the
mean time.

.data: mandatory: your input data set.

	.lo_s.model.span:		smoothing parameter. Default is 0.75.

	.lo_s.model.degree:		overall degree of locally-fitted polynomial. 1 is 
			locally-linear fitting and 2 is locally-quadratic 
			fitting.  Default is 2.

	.lo_s.normalize:	Should numeric predictors 
			be normalized?  If 'y' - the default - the 
			standard normalization is used. If 'n', no
			normalization is carried out.

            TO DO
	parametric:	for two or more numeric predictors, this argument
			specifies those variables that should be 
			conditionally-parametric. The argument should be a 
			logical vector of length p, specified in the order 
			of the predictor group ordered in x.
			Default is a vector of 0's of length p.

            TO DO
	drop_square:	for cases with degree = 2, and with two or more 
			numeric predictors, this argument specifies those 
			numeric predictors whose squares should be dropped 
			from the set of fitting variables. The method of 
			specification is the same as for parametric.
			Default is a vector of 0's of length p.

	.lo_s.model.family:		the assumed distribution of the errors. The values 
			are <tt>"gaussian"</tt> or <tt>"symmetric"</tt>. The first value is 
			the default.  If the second value is specified, 
			a robust fitting procedure is used.

	lo_s.control.surface:	determines whether the fitted surface is computed 
        <tt>"directly"</tt> at all points  or whether an 
        <tt>"interpolation"</tt> method is used. 
			The default, interpolation, is what most users should 
			use unless special circumstances warrant.

    lo_s.control.statistics:	determines whether the statistical quantities are 
        computed <tt>"exactly"</tt> or approximately 
        , where <tt>"approximate"</tt> is the default. The former
        should only be used for testing the approximation in 
        statistical development and is not meant for routine 
        usage because computation time can be horrendous.

        lo_s.control.cell:		if interpolation is used to compute the surface, this
        argument specifies the maximum cell size of the k-d 
        tree.  Suppose k = floor(n*cell*span) where n is the 
        number of observations.  Then a cell is further 
        divided if the number of observations within it
        is greater than or equal to k. default=0.2

	lo_s.control.trace_hat: Options are <tt>"approximate"</tt>, <tt>"exact"</tt>, and <tt>"wait.to.decide"</tt>.	
        When lo_s.control.surface is <tt>"approximate"</tt>, determines
        the computational method used to compute the trace of the hat
        matrix, which is used in the computation of the statistical
        quantities.  If "exact", an exact computation is done; normally
        this goes quite fast on the fastest machines until n, the number
        of observations is 1000 or more, but for very slow machines,
        things can slow down at n = 300.  If "wait.to.decide" is selected,
        then a default is chosen in loess();  the default is "exact" for
        n < 500 and "approximate" otherwise.  If surface is "exact", an
        exact computation is always done for the trace. Set trace_hat to
        "approximate" for large dataset will substantially reduce the
        computation time.

	lo_s.model.iterations:	if family is <tt>"symmetric"</tt>, the number of iterations 
        of the robust fitting method.  Default is 0 for
        lo_s.model.family = gaussian; 4 for family=symmetric.

        That's all you can set. Here are some output parameters:

out	
	fitted_values:	fitted values of the local regression model

	fitted_residuals:	residuals of the local regression fit

        enp:		equivalent number of parameters.

        s:		estimate of the scale of the residuals.

        one_delta:	a statistical parameter used in the computation of standard errors.

        two_delta:	a statistical parameter used in the computation of standard errors.

        pseudovalues:	adjusted values of the response when robust estimation is used.

	trace_hat:	trace of the operator hat matrix.

        diagonal:	diagonal of the operator hat matrix.

        robust:		robustness weights for robust fitting.

        divisor:	normalization divisor for numeric predictors.


struct  anova_struct	*aov;
	
	dfn:		degrees of freedom of the numerator.
	dfd:		degrees of freedom of the denominator.
	F_values:	F statistic.
	Pr_F:		probability F_value is exceeded if null hypothesis is true.

    \ingroup settings
*/
typedef struct {
    apop_data *data;
    struct  loess_struct lo_s;
    int     want_predict_ci; /**< If 'y' (the default), calculate the
                                confidence bands for predicted values */
    double  ci_level; /**< If running a prediction, the level at which
                        to calculate the confidence interval. default:
                        0.95 */
} apop_loess_settings;


/** For the multinomial model. The primary item of interest is the \c bin_ct. For the
  special case of the binomial distribution, you would need:
\code
apop_model *bi = apop_model_copy(apop_multinomial);
apop_model_add_group(bi, apop_multinomial, .bin_ct=2);
apop_estimate(yourdata, bi);
\endcode

The default is to use the maximum value in your data set plus one as the bin count
[the plus one is because there is always a zero bin.]
  */
typedef struct {
    int bin_ct; /**< may be \c 0, meaning autodetect; otherwise, the number of bins. The
          default is \c 0, but you are strongly encouraged to set this to a
          positive integer giving a firm bin count.
          You are so encouraged because your data may simply lack observations
          from the top bin. If your data actually has five bins, but your data
          has draws of <tt> 0 0 2 1 3 2 3</tt>, then I would autodetect a four-bin model,
          because luck gave us no draws of \c 4.
         */
    int compress; /**< default = \c 'y': say there are two bins, zero and one. Say that an
        observation = \c 3. If this is \c 'y', then the observation is placed in the
        one bin. If this is \c 'n', then I emit an error and stop. For a binomial
        model, this means that zeros are put in the zero bin, and everything else
        is put in the one bin, which is probably what you want. */
} apop_multinomial_settings;




/** \defgroup settings Settings*/

//Doxygen is doing funny things right now; having these down here seems to help.
Apop_settings_declarations(apop_histogram)
Apop_settings_declarations(apop_kernel_density)
Apop_settings_declarations(apop_loess)
Apop_settings_declarations(apop_lm)
Apop_settings_declarations(apop_mle)
Apop_settings_declarations(apop_cdf)
Apop_settings_declarations(apop_pm)
Apop_settings_declarations(apop_parts_wanted)
Apop_settings_declarations(apop_multinomial)

#ifdef	__cplusplus
}
#endif
#endif
