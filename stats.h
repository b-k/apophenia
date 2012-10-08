/** \file stats.h */
/* Copyright (c) 2005--2010 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2. */
#ifndef APOP_STATS_H
#define APOP_STATS_H
#include <math.h>
#include "types.h"
#include "variadic.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_randist.h>
//#include <gsl/gsl_multimin.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>

#ifdef	__cplusplus
extern "C" {
#endif

    //First, some linear algebra utilities

double apop_det_and_inv(const gsl_matrix *in, gsl_matrix **out, int calc_det, int calc_inv);
APOP_VAR_DECLARE apop_data * apop_dot(const apop_data *d1, const apop_data *d2, char form1, char form2);
APOP_VAR_DECLARE int         apop_vector_bounded(const gsl_vector *in, long double max);
APOP_VAR_DECLARE void apop_vector_increment(gsl_vector * v, int i, double amt);
APOP_VAR_DECLARE void apop_matrix_increment(gsl_matrix * m, int i, int j, double amt);
gsl_matrix * apop_matrix_inverse(const gsl_matrix *in) ;
double      apop_matrix_determinant(const gsl_matrix *in) ;
//apop_data*  apop_sv_decomposition(gsl_matrix *data, int dimensions_we_want);
APOP_VAR_DECLARE apop_data *  apop_matrix_pca(gsl_matrix *data, int const dimensions_we_want);
APOP_VAR_DECLARE gsl_vector * apop_vector_stack(gsl_vector *v1, gsl_vector * v2, char inplace);
APOP_VAR_DECLARE gsl_matrix * apop_matrix_stack(gsl_matrix *m1, gsl_matrix * m2, char posn, char inplace);
gsl_matrix * apop_matrix_rm_columns(gsl_matrix *in, int *drop);

void apop_vector_log(gsl_vector *v);
void apop_vector_log10(gsl_vector *v);
void apop_vector_exp(gsl_vector *v);

#define APOP_SUBMATRIX(m, srow, scol, nrows, ncols, o) gsl_matrix apop_mm_##o = gsl_matrix_submatrix((m), (srow), (scol), (nrows),(ncols)).matrix;\
gsl_matrix * o = &( apop_mm_##o );

#define APOP_MATRIX_ROW(m, row, v) gsl_vector apop_vv_##v = gsl_matrix_row((m), (row)).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_MATRIX_COL(m, col, v) gsl_vector apop_vv_##v = gsl_matrix_column((m), (col)).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_ROW_T(m, row, v) gsl_vector apop_vv_##v = gsl_matrix_row((m)->matrix, apop_name_find((m)->names, row, 'r')).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_COL_T(m, col, v) gsl_vector apop_vv_##v = gsl_matrix_column((m)->matrix, apop_name_find((m)->names, col, 'c')).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_ROW(m, row, v) gsl_vector apop_vv_##v = gsl_matrix_row((m)->matrix, (row)).vector;\
gsl_vector * v = &( apop_vv_##v );

//-1st column is the vector.
#define APOP_COL(m, col, v) gsl_vector apop_vv_##v = (col < 0) ? *((m)->vector) : gsl_matrix_column((m)->matrix, (col)).vector;\
gsl_vector * v = &( apop_vv_##v );

#define Apop_data_rows(d, rownum, len, outd) \
    gsl_vector apop_dd_##outd##_v = ((d)->vector && (d)->vector->size > (rownum)+(len)-1)  \
                                    ? gsl_vector_subvector((d)->vector, (rownum), (len)).vector\
                                    : (gsl_vector) { };\
    gsl_vector apop_dd_##outd##_w = ((d)->weights && (d)->weights->size > (rownum)+(len)-1)  \
                                    ? gsl_vector_subvector((d)->weights, (rownum), (len)).vector\
                                    : (gsl_vector) { };\
    gsl_matrix apop_dd_##outd##_m = ((d)->matrix && (d)->matrix->size1 > (rownum)+(len)-1)  \
                                ? gsl_matrix_submatrix((d)->matrix, rownum, 0, (len), (d)->matrix->size2).matrix\
                                : (gsl_matrix) { };             \
    apop_name apop_dd_##outd##_n = !((d)->names) ? (apop_name) {} :              \
            (apop_name){                                                         \
                .vector = (d)->names->vector,                                    \
                .column = (d)->names->column,                                    \
                .row = ((d)->names->row && (d)->names->rowct > rownum) ? &((d)->names->row[rownum]) : NULL,  \
                .text = (d)->names->text,                                        \
                .colct = (d)->names->colct,                                      \
                .rowct = (d)->names->row ? (GSL_MIN(len, GSL_MAX((d)->names->rowct - rownum, 0)))      \
                                          : 0,                                   \
                .textct = (d)->names->textct };                                  \
    apop_data apop_dd_##outd = (apop_data){                                      \
                .vector= apop_dd_##outd##_v.size ? &apop_dd_##outd##_v : NULL,   \
                .weights=apop_dd_##outd##_w.size ? &apop_dd_##outd##_w : NULL ,  \
                .matrix = apop_dd_##outd##_m.size1 ? &apop_dd_##outd##_m : NULL, \
                .textsize[0]=(d)->textsize[0] ? (len) : 0, .textsize[1]=(d)->textsize[1],   \
                .text = (d)->text ? &((d)->text[rownum]) : NULL,                 \
                .names= (d)->names ? &apop_dd_##outd##_n : NULL };               \
    apop_data *outd =  &apop_dd_##outd;

#define Apop_data_row(d, row, outd) Apop_data_rows(d, row, 1, outd)

#define Apop_col APOP_COL 
#define apop_col APOP_COL 
#define Apop_row APOP_ROW
#define apop_row APOP_ROW
#define Apop_col_t APOP_COL_T
#define Apop_row_t APOP_ROW_T
#define Apop_matrix_col APOP_MATRIX_COL 
#define Apop_matrix_row APOP_MATRIX_ROW
#define Apop_submatrix APOP_SUBMATRIX
#define apop_data_rows Apop_data_rows
#define apop_data_row Apop_data_row

long double apop_vector_sum(const gsl_vector *in);
double apop_var(const gsl_vector *in);
double apop_vector_var_m(const gsl_vector *in, const double mean);
double apop_vector_cov(const gsl_vector *ina, const gsl_vector *inb);
double apop_vector_correlation(const gsl_vector *ina, const gsl_vector *inb);
double apop_vector_kurtosis_pop(const gsl_vector *in);
double apop_vector_kurtosis(const gsl_vector *in);
double apop_vector_skew(const gsl_vector *in);
double apop_vector_skew_pop(const gsl_vector *in);
double apop_vector_weighted_mean(const gsl_vector *, const gsl_vector *);
double apop_vector_weighted_var(const gsl_vector *v, const gsl_vector *w);
double apop_vector_weighted_cov(const gsl_vector *, const gsl_vector *, const gsl_vector *);
double apop_vector_weighted_skew(const gsl_vector *v, const gsl_vector *w);
double apop_vector_weighted_kurt(const gsl_vector *v, const gsl_vector *w);

#define apop_sum(in) apop_vector_sum(in)
#define apop_var(in) apop_vector_var(in) 
#define apop_vector_covar(in) apop_vector_cov(in) 
#define apop_mean(in) apop_vector_mean(in)

/** Find the mean of the input vector.

*/
#define apop_vector_mean(in)  gsl_stats_mean((in)->data,(in)->stride, (in)->size)

#define apop_vector_var(in)  gsl_stats_variance((in)->data,(in)->stride, (in)->size)
#define apop_vector_kurt(in) apop_vector_kurtosis(in)

APOP_VAR_DECLARE double apop_vector_distance(const gsl_vector *ina, const gsl_vector *inb, const char metric, const double norm);

APOP_VAR_DECLARE void apop_vector_normalize(gsl_vector *in, gsl_vector **out, const char normalization_type);
void apop_matrix_normalize(gsl_matrix *data, const char row_or_col, const char normalization);

APOP_VAR_DECLARE gsl_matrix * apop_matrix_covariance(gsl_matrix *in, const char normalize);
APOP_VAR_DECLARE  gsl_matrix * apop_matrix_correlation(gsl_matrix *in, const char normalize);
apop_data * apop_data_covariance(const apop_data *in);
apop_data * apop_data_correlation(const apop_data *in);
long double apop_matrix_sum(const gsl_matrix *m);
double apop_matrix_mean(const gsl_matrix *data);
void apop_matrix_mean_and_var(const gsl_matrix *data, double *mean, double *var);
apop_data * apop_data_summarize(apop_data *data);

apop_data *apop_test_fisher_exact(apop_data *intab); //in apop_fisher.c

//from apop_t_f_chi.c:
APOP_VAR_DECLARE int apop_matrix_is_positive_semidefinite(gsl_matrix *m, char semi);
double apop_matrix_to_positive_semidefinite(gsl_matrix *m);
double apop_multivariate_gamma(double a, double p);
double apop_multivariate_lngamma(double a, double p);

//apop_tests.c
apop_data *	apop_t_test(gsl_vector *a, gsl_vector *b);
apop_data *	apop_paired_t_test(gsl_vector *a, gsl_vector *b);
APOP_VAR_DECLARE apop_data* apop_anova(char *table, char *data, char *grouping1, char *grouping2);
#define apop_ANOVA apop_anova
APOP_VAR_DECLARE apop_data * apop_f_test (apop_model *est, apop_data *contrast);
#define apop_F_test apop_f_test

//from the regression code:
#define apop_estimate_r_squared(in) apop_estimate_coefficient_of_determination(in)

apop_data * apop_text_unique_elements(const apop_data *d, size_t col);
gsl_vector * apop_vector_unique_elements(const gsl_vector *v);
APOP_VAR_DECLARE apop_data * apop_data_to_factors(apop_data *data, char intype, int incol, int outcol);
APOP_VAR_DECLARE apop_data * apop_data_get_factor_names(apop_data *data, int col, char type);

APOP_VAR_DECLARE apop_data * apop_data_to_dummies(apop_data *d, int col, char type, int keep_first, char append, char remove);

APOP_VAR_DECLARE double apop_kl_divergence(apop_model *top, apop_model *bottom, int draw_ct, gsl_rng *rng);

apop_data *apop_estimate_coefficient_of_determination (apop_model *);
void apop_estimate_parameter_tests (apop_model *est);


//Bootstrapping & RNG
apop_data * apop_jackknife_cov(apop_data *data, apop_model model);
APOP_VAR_DECLARE apop_data * apop_bootstrap_cov(apop_data *data, apop_model model, gsl_rng* rng, int iterations, char keep_boots, char ignore_nans);
gsl_rng *apop_rng_alloc(int seed);
double apop_rng_GHgB3(gsl_rng * r, double* a); //in apop_asst.c


void apop_arms_draw (double *out, gsl_rng *r, apop_model *m); //apop_arms.h


    // maximum likelihod estimation related functions

APOP_VAR_DECLARE gsl_vector * apop_numerical_gradient(apop_data * data, apop_model* model, double delta);
APOP_VAR_DECLARE apop_data * apop_model_hessian(apop_data * data, apop_model *model, double delta);
APOP_VAR_DECLARE apop_data * apop_model_numerical_covariance(apop_data * data, apop_model *model, double delta);

apop_model * apop_maximum_likelihood(apop_data * data, apop_model *dist);

APOP_VAR_DECLARE apop_model * apop_estimate_restart (apop_model *e, apop_model *copy, char * starting_pt, double boundary);

//in apop_linear_constraint.c
APOP_VAR_DECLARE double  apop_linear_constraint(gsl_vector *beta, apop_data * constraint, double margin);

//in apop_model_fix_params.c
apop_model * apop_model_fix_params(apop_model *model_in);
apop_model * apop_model_fix_params_get_base(apop_model *model_in);

#ifdef	__cplusplus
}
#endif
#endif

/** \def Apop_submatrix(m, srow, scol, nrow, ncol)
Pull a pointer to a submatrix into a \c gsl_matrix 

 \param m The root matrix
 \param srow the first row (in the root matrix) of the top of the submatrix
 \param scol the first column (in the root matrix) of the left edge of the submatrix
 \param nrow number of rows in the submatrix
 \param ncol number of columns in the submatrix
\hideinitializer */

/** \def Apop_row_t(m, row_name, v)
 After this call, \c v will hold a vector view of the <tt>row</tt>th row of \c m.
 Unlike \ref Apop_row, the second argument is a row name, that I'll look up using \ref apop_name_find.
\hideinitializer */

/** \def Apop_col_t(m, col_name, v)
 After this call, \c v will hold a vector view of the <tt>col</tt>th column of \c m.
 Unlike \ref Apop_col, the second argument is a column name, that I'll look up using \ref apop_name_find.
\hideinitializer */

/** \def Apop_matrix_row(m, row, v)
 After this call, \c v will hold a vector view of the <tt>row</tt>th row of \c m.
\hideinitializer */

/** \def Apop_matrix_col(m, col, v)
 After this call, \c v will hold a vector view of the <tt>col</tt>th column of \c m.
\hideinitializer */

/** \def Apop_row(d, row, v)
 After this call, \c v will hold a vector view of the <tt>row</tt>th row of \ref apop_data set \c d.
\hideinitializer */

/** \def Apop_col(d, col, v)
 After this call, \c v will hold a vector view of the <tt>col</tt>th column of \ref apop_data set \c d.
\hideinitializer */

/** \def Apop_data_rows(d, row, len, outd)
A macro to generate a temporary view of \ref apop_data set \c d, beginning at row \c row
and having length \c len. 
After this call, \c outd will be a pointer to this temporary
view, that you can use as you would any \ref apop_data set. However, 
it expires as soon as the program leaves the current scope (like with the usual statically declared vars). 
\hideinitializer */

/** \def Apop_data_row(d, row, outd)
A macro to generate a temporary one-row view of \ref apop_data set \c d, pulling out only
row \c row. 
After this call, \c outd will be a pointer to this temporary
view, that you can use as you would any \ref apop_data set. This macro expands to
<tt>Apop_data_rows(d, row, 1, outd)</tt>.
\hideinitializer */

/** \def apop_mean(v)
 Returns the mean of the elements of the vector \c v.
\hideinitializer */

/** \def apop_vector_mean(v)
 Returns the mean of the elements of the vector \c v.
\hideinitializer */

/** \def apop_vector_var(v)
 Returns the sample variance of the elements of the vector \c v.
\hideinitializer */
