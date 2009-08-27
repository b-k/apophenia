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
#include "variadic.h"
#include "conversions.h"

#define MAX_ITERATIONS 		5000
#define MAX_ITERATIONS_w_d	5000

#ifdef	__cplusplus
extern "C" {
#endif

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
    apop_optimization_enum method;
    double      step_size, tolerance, delta;
    int         verbose;
    char        want_cov;
//simulated annealing (also uses step_size);
    int         n_tries, use_score, iters_fixed_T;
    double      k, t_initial, mu_t, t_min ;
    gsl_rng     *rng;
    /** See \ref trace_path */
    char        *trace_path;
    apop_model  *parent;
} apop_mle_settings;

apop_mle_settings *apop_mle_settings_alloc(apop_model *model);
#ifndef SWIG
Apop_settings_declarations(apop_mle)
#endif


typedef double 	(*apop_fn_with_params) (apop_data *, apop_model *);
APOP_VAR_DECLARE gsl_vector * apop_numerical_gradient(apop_data * data, apop_model* model, double delta);
//gsl_vector * apop_numerical_gradient(apop_data * data, apop_model* model);
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
apop_model * apop_model_fix_params(apop_model *model_in);

#ifdef	__cplusplus
}
#endif
#endif
