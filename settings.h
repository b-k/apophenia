/** \file settings.h */ 
/* Copyright (c) 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#ifndef __apop_settings_h__
#define __apop_settings_h__

#include "types.h"
#include "asst.h"

#ifdef	__cplusplus
extern "C" {
#endif

    //Part I: macros and fns for getting/setting settings groups and elements

void * apop_settings_get_group(apop_model *m, char *type);
void apop_settings_rm_group(apop_model *m, char *delme);
void apop_settings_copy_group(apop_model *outm, apop_model *inm, char *copyme);
void *apop_settings_group_alloc(apop_model *model, char *type, void *free_fn, void *copy_fn, void *the_group);

/** Retrieves a settings group from a model.  See \ref Apop_settings_get
 to just pull a single item from within the settings group.*/
#define Apop_settings_get_group(m, type) apop_settings_get_group(m, #type)

/** Removes a settings group from a model's list. */
#define Apop_settings_rm_group(m, type) apop_settings_rm_group(m, #type)

/** Add a settings group. The first two arguments (the model you are
 attaching to and the settings group name) are mandatory, and then you
 can use the \ref designated syntax to specify default values (if any).
 Returns a pointer to the newly-prepped group.
 */
#define Apop_model_add_group(model, type, ...)  \
    apop_settings_group_alloc(model, #type, type ## _settings_free, type ## _settings_copy, type ##_settings_init ((type ## _settings) {__VA_ARGS__})); 

/* A convenience for your settings group init functions. 
 Gives the output item either the defaut value if there is one, or the value you specify. */
#define apop_varad_setting(in, out, name, value) (out)->name = (in).name ? (in).name : (value);

/** Retrieves a setting from a model.  See \ref Apop_settings_get_group pull the entire group.*/
#define Apop_settings_get(model, type, setting)  \
    (((type ## _settings *) apop_settings_get_group(model, #type))->setting)

/** Modifies a single element of a settings group to the given value. */
#define Apop_settings_set(model, type, setting, data)  \
    do { type ## _settings *apop_tmp_settings = apop_settings_get_group(model, #type);  \
    apop_assert_void(apop_tmp_settings, 0, 's', "You're trying to modify a setting in " \
                        #model "'s setting group of type " #type " but that model doesn't have such a group."); \
    apop_tmp_settings->setting = (data);    \
    } while (0);
/*#define Apop_settings_set(model, type, setting, data)  \
    do {                                                \
    apop_assert_void(apop_settings_get_group(model, #type), 0, 's', "You're trying to modify a setting in " \
                        #model "'s setting group of type " #type " but that model doesn't have such a group."); \
    ((type ## _settings *) apop_settings_get_group(model, #type))->setting = (data);    \
    } while (0);*/

#define Apop_settings_add Apop_settings_set
#define APOP_SETTINGS_ADD Apop_settings_set
#define APOP_SETTINGS_GET Apop_settings_get
#define APOP_MODEL_ADD_GROUP Apop_model_add_group
#define APOP_SETTINGS_GET_GROUP Apop_settings_get_group
#define APOP_SETTINGS_RM_GROUP Apop_settings_rm_group

#define Apop_settings_declarations(ysg) \
   ysg##_settings * ysg##_settings_init(ysg##_settings); \
   void * ysg##_settings_copy(ysg##_settings *); \
   void ysg##_settings_free(ysg##_settings *);

        //Part II: the details of extant settings groups.

typedef enum {
    APOP_SIMPLEX_NM     =0, /**< 0: Nelder-Mead simplex (gradient handling rule is irrelevant) */
    APOP_CG_FR     =1,      /**<  1: conjugate gradient (Fletcher-Reeves) (default) */
    APOP_CG_BFGS   =2,      /**<  2: conjugate gradient (BFGS: Broyden-Fletcher-Goldfarb-Shanno) */
    APOP_CG_PR     =3,      /**<  3: conjugate gradient (Polak-Ribiere) */
    APOP_SIMAN      =5,         /**<  5: \ref simanneal "simulated annealing" */
    APOP_RF_NEWTON  =10,        /**<  10: Find a root of the derivative via Newton's method */
//    APOP_RF_BROYDEN =11,        //  11: Find a root of the derivative via the Broyden Algorithm
    APOP_RF_HYBRID  =12,        /**<  12: Find a root of the derivative via the Hybrid method */
    APOP_RF_HYBRID_NOSCALE  =13 /**<  13: Find a root of the derivative via the Hybrid method; no internal scaling */
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
    char        want_cov; /**< Should I calculate a covariance matrix?  Default: 'y', but this can be the most 
                                time-consuming part of the process. */
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
    char want_cov; /**< The covariance can be computationally expensive, so if this is \c 'n' I won't bother with it. */
    char want_expected_value; /**< If 'y', fill the expected/actual/residual part of the output model. */
} apop_ls_settings;

#if 0
// Find apop_category_settings routines in apop_probit.c
/** For dependent-category models, send in this settings struct to specify which column is the dependent variable. 

 If you don't use it, these models will assume that the vector or first numeric column is already a coherent set of factors, but by sending this in, those functions have a little more information, such as names to use in the output.

See also the \ref apop_category_settings_init function.
\ingroup settings
*/
typedef struct {
    apop_data *factors; 
    char source_type; /**< source_type \c 't' = text; anything else (\c 'd' is a good choice) is  numeric data. */
    int source_column; /**<  The number of the column to convert to factors.  As usual, the vector is -1. */
    apop_data *source_data; /**< The input data set that you're probably about to run a regression on */
} apop_category_settings;
#endif


//in apop_exponential.c
/** If this settings group is present, models that can take rank data
  will read the input data as such.  Allocation is thus very simple, e.g.
  \code
  Apop_model_group_add(your_model, apop_rank);
  \endcode 
\ingroup settings
 */
typedef struct {
    char rank_data;
} apop_rank_settings;


#include <gsl/gsl_histogram.h>
/** Settings for the histogram and kernel density structures, mostly opaque. 
  On setup, you must set the \c .data and \c .bins_in items. 
  
  \ingroup settings */
typedef struct{
    apop_data           *data;
    gsl_histogram       *pdf; /**< Where the histogram is kept */
    gsl_histogram_pdf   *cdf; /**< If you make random draws, I need a CDF aggrgation of the main PDF. I keep it here. */
    apop_model          *histobase;
    apop_model          *kernelbase;
    int                 bins_in; /**< Used as input. May not equal the final number of bins (\c pdf->bins) due to the infinibins.*/
} apop_histogram_settings;

#define apop_kernel_density_settings apop_histogram_settings

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

apop_histogram_settings *apop_kernel_density_settings_alloc(apop_data *data, 
        apop_model *histobase, apop_model *kernelbase, void (*set_params)(double, apop_model*));

#define apop_kernel_density_settings_copy apop_histogram_settings_copy
#define apop_kernel_density_settings_free apop_histogram_settings_free


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
        , whiere <tt>"approximate"</tt> is the default. The former
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


/** \defgroup settings Settings*/

//Doxygen is doing funny things right now; having these down 
//here seems to help.
//Apop_settings_declarations(apop_category)
Apop_settings_declarations(apop_histogram)
Apop_settings_declarations(apop_loess)
Apop_settings_declarations(apop_ls)
Apop_settings_declarations(apop_mle)
Apop_settings_declarations(apop_rank)


    ///All of the below is deprecated.

/** Add a settings group. 
  \deprecated{Use \ref Apop_model_add_group instead.}
  You will need to provide arguments for the specific settings group you
  are dealing with, such as <tt>apop_mle_settings_alloc</tt>, <tt>apop_ls_settings_alloc</tt>, <tt>apop_histogram_settings_alloc</tt>. */
#define Apop_settings_add_group(model, type, ...)  \
    apop_settings_group_alloc(model, #type, type ## _settings_free, type ## _settings_copy, type ##_settings_alloc (__VA_ARGS__)); 

#define APOP_SETTINGS_ADD_GROUP Apop_settings_add_group

apop_rank_settings *apop_rank_settings_alloc(void *ignoreme);
apop_histogram_settings *apop_histogram_settings_alloc(apop_data *data, int bins);
apop_mle_settings *apop_mle_settings_alloc(apop_model *model);

#ifdef	__cplusplus
}
#endif
#endif
