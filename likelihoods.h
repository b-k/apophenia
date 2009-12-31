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

APOP_VAR_DECLARE gsl_vector * apop_numerical_gradient(apop_data * data, apop_model* model, double delta);
APOP_VAR_DECLARE apop_data * apop_model_hessian(apop_data * data, apop_model *model, double delta);
APOP_VAR_DECLARE apop_data * apop_model_numerical_covariance(apop_data * data, apop_model *model, double delta);

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
