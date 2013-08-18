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
apop_model *apop_settings_group_alloc_wm(apop_model *model, char *type, void *free_fn, void *copy_fn, void *the_group);

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
 
  If the so-named group is not found, do nothing.
\hideinitializer \ingroup settings
 */
#define Apop_settings_rm_group(m, type) apop_settings_remove_group(m, #type)

/** Add a settings group. The first two arguments (the model you are
 attaching to and the settings group name) are mandatory, and then you
 can use the \ref designated syntax to specify default values (if any).
 \return A pointer to the newly-prepped group.
\hideinitializer \ingroup settings
 */
#define Apop_model_add_group(model, type, ...)  \
    apop_settings_group_alloc(model, #type, type ## _settings_free, type ## _settings_copy, type ##_settings_init ((type ## _settings) {__VA_ARGS__}))

/** Copy a model and add a settings group. Useful for models that require a settings group to function. See \ref Apop_model_add_group.

 \return A pointer to the newly-prepped model.
\hideinitializer \ingroup settings
 */
#define apop_model_copy_set(model, type, ...)  \
    apop_settings_group_alloc_wm(apop_model_copy(model), #type, type ## _settings_free, type ## _settings_copy, type ##_settings_init ((type ## _settings) {__VA_ARGS__}))

/** Retrieves a setting from a model.  See \ref Apop_settings_get_group to pull the entire group.
\hideinitializer \ingroup settings
 */
#define Apop_settings_get(model, type, setting)  \
    (((type ## _settings *) apop_settings_get_grp(model, #type, 'f'))->setting)

/** Modifies a single element of a settings group to the given value. 

\li If <tt>model==NULL</tt>, fails silently. 
\li If <tt>model!=NULL</tt> but the given settings group is not found attached to the model, set <tt>model->error='s'</tt>.
\hideinitializer \ingroup settings
 */
#define Apop_settings_set(model, type, setting, data)   \
    do {                                                \
        if (!(model)) continue; /* silent fail. */      \
        type ## _settings *apop_tmp_settings = apop_settings_get_grp(model, #type, 'c');  \
        Apop_stopif(!apop_tmp_settings, (model)->error='s', 0, "You're trying to modify a setting in " \
                        #model "'s setting group of type " #type " but that model doesn't have such a group."); \
    apop_tmp_settings->setting = (data);                \
    } while (0);

/** \cond doxy_ignore */
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
#define Apop_model_copy_set apop_model_copy_set

/** \endcond */ //End of Doxygen ignore.

#define Apop_settings_declarations(ysg) \
   ysg##_settings * ysg##_settings_init(ysg##_settings); \
   void * ysg##_settings_copy(ysg##_settings *); \
   void ysg##_settings_free(ysg##_settings *);

/** A convenience macro for declaring the initialization function for a new settings group.
 See the documentation outline -> models -> model settings -> writing new settings group for details.

  This sets the defaults for every element in the structure, so you will want a line for every element of your structure (except the ones that default to NULL, which have already been set as such).

  \code
  Apop_settings_init (ysg, 
        Apop_varad_set(size1, 99);
        Apop_varad_set(size2, 2.3);
        Apop_varad_set(dataset, apop_data_alloc(out->size1, out->size2));
    )
  \endcode
  If you need them, the input is a structure named \c in, and the output a pointer-to-struct named \c out.
*/
#define Apop_settings_init(name, ...)   \
    name##_settings *name##_settings_init(name##_settings in) {       \
        name##_settings *out = malloc(sizeof(name##_settings));     \
        *out = in; \
        __VA_ARGS__;            \
        return out; \
    }

#define Apop_varad_set(var, value) (out)->var = (in).var ? (in).var : (value);

/** A convenience macro for declaring the copy function for a new settings group.
 See the documentation outline -> models -> model settings -> writing new settings group for details.

  To just do a direct copy, the default works; let your settings group be named \c ysg:
  \code
Apop_settings_copy (ysg, )
  \endcode
  generates a function that allocates space for a new settings group and copies all elements from the input group to the output group.

  The space after the comma indicates that there is no new procedural code. If you want to add some, feel free. E.g.,
  \code
Apop_settings_copy (ysg, 
    if (!in->score)
        out->score = 1;
    out->data_owner = 0;
)
  \endcode
  The names \c in and \c out are built into the macro.
*/
#define Apop_settings_copy(name, ...) \
    void * name##_settings_copy(name##_settings *in) {\
        name##_settings *out = malloc(sizeof(name##_settings)); \
        *out = *in; \
        __VA_ARGS__;    \
        return out;     \
    }

/** A convenience macro for declaring the delete function for a new settings group.
 See the documentation outline -> models -> model settings -> writing new settings group for details.

If you don't have internal structure elements to free, let your settings group be named \c ysg:
  \code
  Apop_settings_free (ysg, )
  \endcode
  generates a function that simply frees the input settings group.

  If your structure is pointing to other structures that need to be freed first, then add them after that comma:
  \code
Apop_settings_copy (ysg, 
    apop_data_free(in->dataset);
)
  \endcode
  The name \c in is built into the macro.
*/
#define Apop_settings_free(name, ...) \
    void name##_settings_free(name##_settings *in) {\
        __VA_ARGS__;    \
        free(in);  \
    }

//see deprecated.h for the apop_settings_add_group


        //Part II: the details of extant settings groups.

typedef enum {
    APOP_SIMPLEX_NM     =0, /**< Nelder-Mead simplex (gradient handling rule is irrelevant) */
    APOP_CG_FR     =1,      /**<  Conjugate gradient (Fletcher-Reeves) (default) */
    APOP_CG_BFGS   =2,      /**<  Conjugate gradient (BFGS: Broyden-Fletcher-Goldfarb-Shanno) */
    APOP_CG_PR     =3,      /**<  Conjugate gradient (Polak-Ribiere) */
    APOP_SIMAN      =5,         /**<  \ref simanneal "simulated annealing" */
    APOP_RF_NEWTON  =10,        /**<  Find a root of the derivative via Newton's method */
//    APOP_RF_BROYDEN =11,        //  Find a root of the derivative via the Broyden Algorithm
    APOP_RF_HYBRID  =12,        /**<  Find a root of the derivative via the Hybrid method */
    APOP_RF_HYBRID_NOSCALE = 13, /**<  Find a root of the derivative via the Hybrid method; no internal scaling */
    APOP_UNKNOWN_ML = -1       /**<  For internal use */
} apop_optimization_enum;

/** The settings for maximum likelihood estimation (including simulated annealing).
\ingroup settings */
typedef struct{
    double      *starting_pt;   /**< An array of doubles (i.e., <tt>double*</tt>) suggesting a starting point. 
                                  If NULL, use an all-ones vector.  Note that if \c v is a \c gsl_vector, then 
                                  \c v->data is of the right form (provided \c v is not a slice of a matrix).*/
    apop_optimization_enum method; /**< See the  \ref apop_optimization_enum documentation for options. */
    double      step_size, /**< the initial step size. */
                tolerance, /**< the precision the minimizer uses. Only vaguely related to the precision of the actual variables. */
delta;
    int         max_iterations; /**< Ignored by simulated annealing. Other methods halt if
                                 they do this many iterations without finding an optimum. */
    int         verbose; /**<	Give status updates as we go.  This is orthogonal to the 
                                <tt>apop_opts.verbose</tt> setting. */
    double      dim_cycle_tolerance; /**< If zero (the default), the usual procedure.
                             If \f$>0\f$, cycle across dimensions: fix all but the first dimension at the starting
                             point, optimize only the first dim. Then fix the all but the second dim, and optimize the
                             second dim. Continue through all dims, until the log likelihood at the outset of one cycle
                             through the dimensions is within this amount of the previous cycle's log likelihood. There
                             will be at least two cycles.
                             */
//simulated annealing (also uses step_size);
    int         n_tries, iters_fixed_T;
    double      k, t_initial, mu_t, t_min ;
    gsl_rng     *rng;
    char        want_path; ///< If 'y', record the points tried by the optimizer in path
    apop_data   *path;      /**< if want_path='y', record each vector tried by the optimizer as one row of this \ref apop_data set.
                              If already allocated, free what is here and reallocate. This data set has no names; add them as desired.*/
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
\ingroup settings
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
    gsl_matrix *draws_made; /**< A store of random draws that I will count up to report the CDF. Need only be generated once, and so stored here. */
    int rng_owner; /**< For internal use. Should I free the RNG when this copy of the settings group is freed? */
    int draws_owner; /**< For internal use.  Should I free \c draws_made when this copy of the settings group is freed?*/
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



/** Settings to accompany the \ref apop_pmf. */
typedef struct {
    gsl_vector *cmf;  /**< A cumulative mass function, for the purposes of making random draws.*/
    char draw_index;  /**< If \c 'y', then draws from the PMF return the integer index of the row drawn. 
                           If \c 'n' (the default), then return the data in the vector/matrix elements of the data set. */
    long double total_weight; /**< Keep the total weight, in case the input weights aren't normalized to sum to one. */
    int cmf_refct;    /**< For internal use, so I can garbage-collect the CMF when needed. */
} apop_pmf_settings;


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
    apop_data *starting_pt; /**< Deprecated and ignored. Starting point is drawn from your prior. */
    long int periods; /**< For how many steps should the MCMC chain run? */
    double burnin; /**< What <em>percentage</em> of the periods should be ignored
                         as initialization. That is, this is a number between zero and one. */
    int histosegments; /**< If outputting a binned PMF, how many segments should it have? */
    char method;
} apop_update_settings;

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

/** The code for the loess system is based on FORTRAN code from 1988,
overhauled in 1992, linked in to Apophenia in 2009. The structure that
does all the work, then, is a \c loess_struct that you should
basically take as opaque. 

The useful settings from that struct re-appear in the \ref
apop_loess_settings struct so you can set them directly, and then the
settings init function will copy your preferences into the working struct.

The documentation for the elements is cut/pasted/modified from Cleveland,
Grosse, and Shyu.

<tt>.data</tt>: Mandatory. Your input data set.

	<tt>.lo_s.model.span</tt>:	smoothing parameter. Default is 0.75.

	<tt>.lo_s.model.degree</tt>: overall degree of locally-fitted polynomial. 1 is
			locally-linear fitting and 2 is locally-quadratic fitting. Default is 2.

	<tt>.lo_s.normalize</tt>:	Should numeric predictors
			be normalized?	If 'y' - the default - the standard normalization
			is used. If 'n', no normalization is carried out.

	\c .lo_s.model.parametric:	for two or more numeric predictors, this argument
			specifies those variables that should be
			conditionally-parametric. The argument should be a logical
			vector of length p, specified in the order of the predictor
			group ordered in x.  Default is a vector of 0's of length p.

	\c .lo_s.model.drop_square:	for cases with degree = 2, and with two or more
			numeric predictors, this argument specifies those numeric
			predictors whose squares should be dropped from the set of
			fitting variables. The method of specification is the same as
			for parametric.  Default is a vector of 0's of length p.

	\c .lo_s.model.family: the assumed distribution of the errors. The values are
	        <tt>"gaussian"</tt> or <tt>"symmetric"</tt>. The first value is the default.
            If the second value is specified, a robust fitting procedure is used.

	\c lo_s.control.surface:	determines whether the fitted surface is computed
            <tt>"directly"</tt> at all points  or whether an <tt>"interpolation"</tt>
            method is used. The default, interpolation, is what most users should use
			unless special circumstances warrant.

    \c lo_s.control.statistics:	determines whether the statistical quantities are 
        computed <tt>"exactly"</tt> or approximately, where <tt>"approximate"</tt>
        is the default. The former should only be used for testing the approximation in
        statistical development and is not meant for routine usage because computation
        time can be horrendous.

        \c lo_s.control.cell: if interpolation is used to compute the surface,
        this argument specifies the maximum cell size of the k-d tree. Suppose k =
        floor(n*cell*span) where n is the number of observations.  Then a cell is
        further divided if the number of observations within it is greater than or
        equal to k. default=0.2

	\c lo_s.control.trace_hat: Options are <tt>"approximate"</tt>, <tt>"exact"</tt>, and <tt>"wait.to.decide"</tt>.	
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

	\c lo_s.model.iterations:	if family is <tt>"symmetric"</tt>, the number of iterations 
        of the robust fitting method.  Default is 0 for
        lo_s.model.family = gaussian; 4 for family=symmetric.

        That's all you can set. Here are some output parameters:

	\c fitted_values:	fitted values of the local regression model

	\c fitted_residuals:	residuals of the local regression fit

       \c  enp:		equivalent number of parameters.

       \c  s:		estimate of the scale of the residuals.

       \c  one_delta:	a statistical parameter used in the computation of standard errors.

       \c  two_delta:	a statistical parameter used in the computation of standard errors.

       \c  pseudovalues:	adjusted values of the response when robust estimation is used.

	\c trace_hat:	trace of the operator hat matrix.

       \c  diagonal:	diagonal of the operator hat matrix.

       \c  robust:		robustness weights for robust fitting.

       \c  divisor:	normalization divisor for numeric predictors.

    \ingroup settings
*/
typedef struct {
    apop_data *data;
    struct  loess_struct lo_s;
    int     want_predict_ci; /**< If 'y' (the default), calculate the
                                confidence bands for predicted values */
    double  ci_level; /**< If running a prediction, the level at which
                        to calculate the confidence interval. default: 0.95 */
} apop_loess_settings;


    /** \cond doxy_ignore */
typedef struct point {    /* a point in the x,y plane */
  double x,y;             /* x and y coordinates */
  double ey;              /* exp(y-ymax+YCEIL) */
  double cum;             /* integral up to x of rejection envelope */
  int f;                  /* is y an evaluated point of log-density */
  struct point *pl,*pr;   /* envelope points to left and right of x */
} POINT;

/* This includes the envelope info and the metropolis steps. */
typedef struct {  /* attributes of the entire rejection envelope */
  int cpoint;              /* number of POINTs in current envelope */
  int npoint;              /* max number of POINTs allowed in envelope */
  double ymax;             /* the maximum y-value in the current envelope */
  POINT *p;                /* start of storage of envelope POINTs */
  double *convex;          /* adjustment for convexity */
  double metro_xprev;      /* previous Markov chain iterate */
  double metro_yprev;      /* current log density at xprev */
} arms_state;
    /** \endcond */

/** to perform derivative-free adaptive rejection sampling with metropolis step */
typedef struct {
    double *xinit;  /**< A <tt>double*</tt> giving starting values for x in ascending order. Default: -1, 0, 1. If this isn't \c NULL, I need at least three items. */
    double  xl;     /**< Left bound. If you don't give me one, I'll use min[min(xinit)/10, min(xinit)*10].*/
    double  xr;     /**< Right bound. If you don't give me one, I'll use max[max(xinit)/10, max(xinit)*10]. */
    double convex;  /**< Adjustment for convexity */
    int ninit;      /**< Number of starting values supplied (i.e. number of elements in \c xinit)*/
    int npoint;     /**< Maximum number of envelope points. I \c malloc space for this many <tt>double</tt>s at the outset. Default = 1e5. */
   char do_metro;   /**< Whether metropolis step is required. (I.e., set to one if you're not sure if the function is log-concave). Set  to <tt>'y'</tt>es or <tt>'n'</tt>o*/
   double xprev;    /**< Previous value from Markov chain */
   int neval;       /**< On exit, the number of function evaluations performed */
   arms_state *state;
   apop_model *model; /**< The model from which I will draw. Mandatory. Must have either a \c log_likelihood or \c p method.*/
} apop_arms_settings;


typedef struct {
    char *splitpage;    /**< The name of the page at which to split the data. If \c NULL, I send the entire data set to both models as needed. */
    apop_model *model1; /**< The first model in the stack.*/
    apop_model *model2; /**< The second model.*/
} apop_stack_settings;

typedef struct {
    apop_data *(*base_to_transformed)(apop_data*);
    apop_data *(*transformed_to_base)(apop_data*);
    double (*jacobian_to_base)(apop_data*);
    apop_model *base_model;
} apop_ct_settings;/**< All of the elements of this struct should be considered private.*/

typedef struct {
    apop_model *base_model; /**< The model, before constraint. */
    double (*constraint)(apop_data *, apop_model *); /**< The constraint. Return 1 if the data is in the constraint; zero if out. */
    double (*scaling)(apop_model *); /**< Optional. Return the percent of the model density inside the constraint. */
    gsl_rng *rng; /**< If you don't provide a \c scaling function, I calculate the in-constraint model density via random draws.
                       If no \c rng is provided, I use a default RNG; see \ref autorng. */
    int draw_ct; /**< How many draws to make for calculating the in-constraint model density via random draws. Current default: 1e4. */
} apop_dconstrain_settings; /**< For use with the \ref apop_dconstrain model. See its documentation for an example. */

typedef struct {
    apop_model *generator_m;
    apop_model *ll_m;
    gsl_rng *rng;
    int draw_ct;
} apop_composition_settings;/**< All of the elements of this struct should be considered private.*/

/** \defgroup settings Settings*/

//Doxygen drops whatever is after these declarations, so I put them last.
Apop_settings_declarations(apop_ct)
Apop_settings_declarations(apop_lm)
Apop_settings_declarations(apop_pm)
Apop_settings_declarations(apop_pmf)
Apop_settings_declarations(apop_mle)
Apop_settings_declarations(apop_cdf)
Apop_settings_declarations(apop_arms)
Apop_settings_declarations(apop_loess)
Apop_settings_declarations(apop_stack)
Apop_settings_declarations(apop_update)
Apop_settings_declarations(apop_dconstrain)
Apop_settings_declarations(apop_composition)
Apop_settings_declarations(apop_parts_wanted)
Apop_settings_declarations(apop_kernel_density)

#ifdef	__cplusplus
}
#endif
#endif
