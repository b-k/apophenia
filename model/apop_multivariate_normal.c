/* apop_multivariate_normal.c  The multivariate Normal distribution.
 Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.

\amodel apop_multivariate_normal This is the multivariate generalization of the Normal distribution.

\adoc    Input_format     Each row of the matrix is an observation.
\adoc    Parameter_format  An \c apop_data set whose vector element is the vector of
                            means, and whose matrix is the covariances.

If you had only one dimension, the mean would be a vector of size one, and the covariance matrix a \f$1\times 1\f$ matrix. This differs from the setup for \ref apop_normal, which outputs a single vector with \f$\mu\f$ in element zero and \f$\sigma\f$ in element one.
\adoc    Settings   None.    */
 
#include "apop_internal.h"

static double x_prime_sigma_x(gsl_vector *x, gsl_matrix *sigma){
    gsl_vector * sigma_dot_x = gsl_vector_calloc(x->size);
    double the_result;
    gsl_blas_dsymv(CblasUpper, 1, sigma, x, 0, sigma_dot_x); //sigma should be symmetric
    gsl_blas_ddot(x, sigma_dot_x, &the_result);
    gsl_vector_free(sigma_dot_x);
    return the_result;
}

static long double apop_multinormal_ll(apop_data *data, apop_model * m){
    Nullcheck_mpd(data, m, GSL_NAN);
    double determinant = 0;
    gsl_matrix* inverse = NULL;
    int i, dimensions  = data->matrix->size2;
    gsl_vector* x_minus_mu= gsl_vector_alloc(data->matrix->size2);
    double ll = 0;
    determinant = apop_det_and_inv(m->parameters->matrix, &inverse, 1,1);
    Apop_stopif(isnan(determinant) || determinant == 0,
        gsl_vector_free(x_minus_mu); return GSL_NEGINF, //tell maximizers to look elsewhere.
         1, "the determinant of the given covariance is zero or NaN. Returning GSL_NEGINF."); 
    Apop_stopif(determinant <= 0, return NAN, 0, "The determinant of the covariance matrix you gave me "
            "is negative, but a covariance matrix must always be positive semidefinite "
            "(and so have nonnegative determinant). Maybe run apop_matrix_to_positive_semidefinite?");
    for (i=0; i< data->matrix->size1; i++){
        Apop_row_v(data, i, vv);
        gsl_vector_memcpy(x_minus_mu, vv);
        gsl_vector_sub(x_minus_mu, m->parameters->vector);
        ll += - x_prime_sigma_x(x_minus_mu, inverse) / 2;
    }
    ll -= data->matrix->size1 * (log(2 * M_PI)* dimensions/2. + .5 * log(determinant));
    gsl_matrix_free(inverse);
    gsl_vector_free(x_minus_mu);
    return ll;
}

static double a_mean(gsl_vector * in){ return apop_vector_mean(in); }

/*\adoc  estimated_parameters  Format as above. The <tt>\<Covariance\></tt> page gives
the covariance matrix of the means.

\adoc estimated_info   Reports <tt>log likelihood</tt>.  */ 
static void multivariate_normal_estimate(apop_data * data, apop_model *p){
    p->parameters = apop_map(data, .fn_v=a_mean, .part='c'); 
    apop_data *cov =  apop_data_covariance(data);
    p->parameters->matrix =  cov->matrix;
    cov->matrix = NULL; apop_data_free(cov);
    apop_data_add_named_elmt(p->info, "log likelihood", apop_multinormal_ll(data, p));
}

/* \adoc    RNG  The RNG fills an input array whose length is based on the input parameters.

 The nice, easy method from Devroye, p 565 */
static int mvnrng(double *out, gsl_rng *r, apop_model *eps){
    apop_data *params = eps->parameters;
    gsl_vector *v = gsl_vector_alloc(params->vector->size);
    gsl_vector *dotted = gsl_vector_calloc(params->vector->size);
    for (size_t i=0; i< params->vector->size; i++)
        gsl_vector_set(v, i, gsl_ran_gaussian(r, 1));
    gsl_matrix *copy  = apop_matrix_copy(params->matrix);
        gsl_linalg_cholesky_decomp(copy); //returns upper and lower triangle; we want just one.
    for (size_t i=0; i< copy->size1; i++)
        for (size_t j=i+1; j< copy->size2; j++)
            gsl_matrix_set(copy, i, j, 0);
    gsl_blas_dgemv(CblasNoTrans, 1, copy, v, 0, dotted);
    for (size_t i=0; i< params->vector->size; i++)
        out[i]  = gsl_vector_get(dotted, i) + gsl_vector_get(params->vector,i);
    gsl_vector_free(v);
    gsl_vector_free(dotted);
    gsl_matrix_free(copy);
    return 0;
}

static void mvn_prep(apop_data *d, apop_model *m){
    if (d && d->matrix)    m->dsize = d->matrix->size2; 
    else if (m->vsize > 0) m->dsize = m->vsize;
    apop_model_clear(d, m);
}

static long double mvn_constraint(apop_data *d, apop_model *m){
    return apop_matrix_to_positive_semidefinite(m->parameters->matrix);
}

apop_model *apop_multivariate_normal= &(apop_model){"Multivariate normal distribution", -1,-1,-1, .dsize=-2,
     .estimate = multivariate_normal_estimate, .log_likelihood = apop_multinormal_ll, 
     .draw = mvnrng, .prep=mvn_prep, .constraint = mvn_constraint};
