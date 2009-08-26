/** \file apop_multivariate_normal.c  The multivariate Normal distribution.*/
/* Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
 
#include "asst.h"
#include "model.h"
#include "stats.h"
#include "conversions.h"

apop_model apop_multivariate_normal;

static double x_prime_sigma_x(gsl_vector *x, gsl_matrix *sigma){
  gsl_vector *  sigma_dot_x = gsl_vector_calloc(x->size);
  double        the_result;
    gsl_blas_dsymv(CblasUpper, 1, sigma, x, 0, sigma_dot_x); //sigma should be symmetric
    gsl_blas_ddot(x, sigma_dot_x, &the_result);
    gsl_vector_free(sigma_dot_x);
    return the_result;
}

static double apop_multinormal_ll(apop_data *data, apop_model * m){
  apop_assert(m->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
  double    determinant = 0;
  gsl_matrix* inverse   = NULL;
  int       i, dimensions  = data->matrix->size2;
  gsl_vector* x_minus_mu= gsl_vector_alloc(data->matrix->size2);
  double        ll      = 0;
    determinant = apop_det_and_inv(m->parameters->matrix, &inverse, 1,1);
    if (determinant == 0) {
        fprintf(stderr, "apop_multivariate_normal.p: the determinant of the given covariance is zero. Returning GSL_NEGINF.  \n"); 
        gsl_vector_free(x_minus_mu);
        return GSL_NEGINF; //tell maximizers to look elsewhere.
    }
    for (i=0; i< data->matrix->size1; i++){
        APOP_ROW(data,i, vv);
        gsl_vector_memcpy(x_minus_mu, vv);
        gsl_vector_sub(x_minus_mu, m->parameters->vector);
       ll  += - x_prime_sigma_x(x_minus_mu, inverse) / 2;
       ll  -= log(2 * M_PI)* dimensions/2. + .5 *  log(determinant);
    }
    gsl_matrix_free(inverse);
    gsl_vector_free(x_minus_mu);
    return ll;
}

#include "mapply.h"
#include "stats.h"
static double a_mean(gsl_vector * in){ return apop_vector_mean(in); }

static apop_model * multivariate_normal_estimate(apop_data * data, apop_model *p){
  apop_model *out = p ? apop_model_copy(*p) : apop_model_copy(apop_multivariate_normal);
    out->parameters         = apop_map(data, .fn_v=a_mean, .part='c');
    out->covariance         =  apop_data_covariance(data);
    out->parameters->matrix =  out->covariance->matrix;
    out->llikelihood = apop_multinormal_ll(data, out);
    return out;
}

/** The nice, easy method from Devroye, p 565 */
static void mvnrng(double *out, gsl_rng *r, apop_model *eps){
  apop_data *params = eps->parameters;
  int i, j;
  gsl_vector *v     = gsl_vector_alloc(params->vector->size);
  gsl_vector *dotted= gsl_vector_calloc(params->vector->size);
    for (i=0; i< params->vector->size; i++)
        gsl_vector_set(v, i, gsl_ran_gaussian(r, 1));
  gsl_matrix *copy  = apop_matrix_copy(params ->matrix);
        gsl_linalg_cholesky_decomp(copy); //returns upper and lower triangle; we want just one.
    for (i=0; i< copy->size1; i++)
        for (j=i+1; j< copy->size2; j++)
            gsl_matrix_set(copy, i, j, 0);
    gsl_blas_dgemv(CblasNoTrans, 1, copy, v, 0, dotted);
    for (i=0; i< params->vector->size; i++)
        out[i]  = gsl_vector_get(dotted, i) + gsl_vector_get(params->vector,i);
    gsl_vector_free(v);
    gsl_vector_free(dotted);
    gsl_matrix_free(copy);
}

/** This is the multivarate generalization of the Normal distribution.
  The probability/log_likelihood methods take in an \c apop_data set whose vector element is the vector of
  means, and whose matrix is the covariances; the estimate method
  returns parameters of that form.

  The RNG fills an input array whose length is based on the input parameters.

\hideinitializer
  \ingroup models
 */
apop_model apop_multivariate_normal= {"Multivariate normal distribution", -1,-1,-1,
     .estimate = multivariate_normal_estimate, .log_likelihood = apop_multinormal_ll, .draw = mvnrng};
