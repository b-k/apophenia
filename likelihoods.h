/** \file likelihoods.h	 */
 /* Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#ifndef apop_likelihoods_h
#define  apop_likelihoods_h

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include "stats.h"
#include "types.h"
#include "settings.h"
#include "conversions.h"

#define MAX_ITERATIONS 		5000
#define MAX_ITERATIONS_w_d	5000

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

typedef enum {
    APOP_SIMPLEX_NM     =0, // 0: Nelder-Mead simplex (gradient handling rule is irrelevant)
    APOP_CG_FR     =1,      //  1: conjugate gradient (Fletcher-Reeves) (default)
    APOP_CG_BFGS   =2,      //  2: conjugate gradient (BFGS: Broyden-Fletcher-Goldfarb-Shanno)
    APOP_CG_PR     =3,      //  3: conjugate gradient (Polak-Ribiere)
    APOP_SIMAN      =5,         //  5: \ref simanneal "simulated annealing"
    APOP_RF_NEWTON  =10,        //  10: Find a root of the derivative via Newton's method
//    APOP_RF_BROYDEN =11,        //  11: Find a root of the derivative via the Broyden Algorithm
    APOP_RF_HYBRID  =12,        //  12: Find a root of the derivative via the Hybrid method
    APOP_RF_HYBRID_NOSCALE  =13 //  13: Find a root of the derivative via the Hybrid method; no internal scaling
} apop_optimization_enum;

/** The settings for maximum likelihood estimation (including simulated annealing).*/
typedef struct{
//traditional
    double      *starting_pt;
    double      step_size;
    double      tolerance;
    double      delta;
    apop_optimization_enum method;
    int         verbose;
    int         want_cov;
//simulated annealing (also uses step_size);
    int         n_tries;
    int         use_score;
    int         iters_fixed_T;
    double      k, t_initial, mu_t, t_min ;
    gsl_rng     *rng;
    char        *trace_path;
    apop_model  *parent;
} apop_mle_settings;

apop_mle_settings *apop_mle_settings_alloc(apop_model *model);
apop_mle_settings *apop_mle_settings_init(apop_mle_settings in);
void *apop_mle_settings_copy(apop_mle_settings * in);
void apop_mle_settings_free(void * in);


typedef double 	(*apop_fn_with_params) (apop_data *, apop_model *);
gsl_vector * apop_numerical_gradient(apop_data *data, apop_model*);
//gsl_matrix * apop_numerical_second_derivative(apop_model dist, gsl_vector *beta, apop_data * d);
//gsl_matrix * apop_numerical_hessian(apop_model dist, gsl_vector *beta, apop_data * d);

/* Find the var/covar matrix via the hessian. */
//void apop_numerical_covariance_matrix(apop_model dist, apop_model *est, apop_data *data);
//void apop_numerical_var_covar_matrix(apop_model dist, apop_model *est, apop_data *data);

apop_model *	apop_maximum_likelihood(apop_data * data, apop_model dist);

APOP_VAR_DECLARE apop_model * apop_estimate_restart (apop_model *e, apop_model *copy, char * starting_pt, double boundary);

//in apop_linear_constraint.c
APOP_VAR_DECLARE double  apop_linear_constraint(gsl_vector *beta, apop_data * constraint, double margin);

//in apop_model_fix_params.c
apop_model *apop_model_fix_params(apop_data *data, apop_data *paramvals, apop_data *mask, apop_model model_in);
__END_DECLS
#endif
