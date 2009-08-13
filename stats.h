/** \file stats.h
Copyright (c) 2005--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.
*/
#ifndef APOP_STATS_H
#define APOP_STATS_H
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
#include "linear_algebra.h"

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

#define APOP_SUBMATRIX(m, srow, scol, nrows, ncols, o) gsl_matrix apop_mm_##o = gsl_matrix_submatrix(m, (srow), (scol), (nrows),(ncols)).matrix;\
gsl_matrix * o = &( apop_mm_##o );

#define APOP_MATRIX_ROW(m, row, v) gsl_vector apop_vv_##v = gsl_matrix_row(m, (row)).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_MATRIX_COL(m, col, v) gsl_vector apop_vv_##v = gsl_matrix_column(m, (col)).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_ROW_T(m, row, v) gsl_vector apop_vv_##v = gsl_matrix_row((m)->matrix, apop_name_find((m)->names, row, 'r')).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_COL_T(m, col, v) gsl_vector apop_vv_##v = gsl_matrix_column((m)->matrix, apop_name_find((m)->names, col, 'c')).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_ROW(m, row, v) gsl_vector apop_vv_##v = gsl_matrix_row((m)->matrix, (row)).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_COL(m, col, v) gsl_vector apop_vv_##v = gsl_matrix_column((m)->matrix, (col)).vector;\
gsl_vector * v = &( apop_vv_##v );

#define Apop_col APOP_COL 
#define Apop_row APOP_ROW
#define Apop_col_t APOP_COL_T
#define Apop_row_t APOP_ROW_T
#define Apop_matrix_col APOP_MATRIX_COL 
#define Apop_matrix_row APOP_MATRIX_ROW
#define Apop_submatrix APOP_SUBMATRIX
#define apop_vector_kurt(in) apop_vector_kurtosis(in)

	//The following are just convenient hooks to gsl vector functions.
#define __PURE __attribute__((pure))
inline long double apop_vector_sum(const gsl_vector *in) __PURE;
inline double apop_vector_mean(const gsl_vector *in) __PURE;
inline double apop_vector_var(const gsl_vector *in) __PURE;
inline double apop_var(const gsl_vector *in) __PURE;
inline double apop_vector_var_m(const gsl_vector *in, const double mean) __PURE;
inline double apop_vector_cov(const gsl_vector *ina, const gsl_vector *inb) __PURE;
inline double apop_vector_correlation(const gsl_vector *ina, const gsl_vector *inb) __PURE;
inline double apop_vector_kurtosis_pop(const gsl_vector *in) __PURE;
inline double apop_vector_kurtosis(const gsl_vector *in) __PURE;
inline double apop_vector_skew(const gsl_vector *in) __PURE;
inline double apop_vector_skew_pop(const gsl_vector *in) __PURE;
double apop_vector_weighted_mean(const gsl_vector *, const gsl_vector *) __PURE;
double apop_vector_weighted_var(const gsl_vector *v, const gsl_vector *w) __PURE;
double apop_vector_weighted_cov(const gsl_vector *, const gsl_vector *, const gsl_vector *) __PURE;
double apop_vector_weighted_skew(const gsl_vector *v, const gsl_vector *w) __PURE;
double apop_vector_weighted_kurt(const gsl_vector *v, const gsl_vector *w) __PURE;

#define apop_sum(in) apop_vector_sum(in)
#define apop_var(in) apop_vector_var(in) 
#define apop_vector_covar(in) apop_vector_cov(in) 
#define apop_mean(in) apop_vector_mean(in)
#define apop_vector_mean(in)  gsl_stats_mean((in)->data,(in)->stride, (in)->size)
#define apop_vector_var(in)  gsl_stats_variance((in)->data,(in)->stride, (in)->size)

APOP_VAR_DECLARE double apop_vector_distance(const gsl_vector *ina, const gsl_vector *inb, const char metric, const double norm);
double apop_vector_grid_distance(const gsl_vector *ina, const gsl_vector *inb);


APOP_VAR_DECLARE void apop_vector_normalize(gsl_vector *in, gsl_vector **out, const char normalization_type);
void apop_matrix_normalize(gsl_matrix *data, const char row_or_col, const char normalization);

inline double apop_test_chi_squared_var_not_zero(const gsl_vector *in);
	//As described: give it a vector, and it'll tell you the confidence 
	//with which you can say that the vector is not zero.

APOP_VAR_DECLARE double apop_random_double(double min, double max, gsl_rng *r);
APOP_VAR_DECLARE int apop_random_int(double min, double max, const gsl_rng *r);

APOP_VAR_DECLARE gsl_matrix * apop_matrix_covariance(gsl_matrix *in, const char normalize);
APOP_VAR_DECLARE  gsl_matrix * apop_matrix_correlation(gsl_matrix *in, const char normalize);
apop_data * apop_data_covariance(const apop_data *in);
apop_data * apop_data_correlation(const apop_data *in);
long double apop_matrix_sum(const gsl_matrix *m) __PURE;
double apop_matrix_mean(const gsl_matrix *data) __PURE;
double apop_matrix_var_m(const gsl_matrix *data, double mean) __PURE;
void apop_matrix_mean_and_var(const gsl_matrix *data, double *mean, double *var);
double apop_rng_GHgB3(gsl_rng * r, double* a); //in asst.c
apop_data * apop_data_summarize(apop_data *data);

//from apop_fisher.c:
apop_data *apop_test_fisher_exact(apop_data *intab);

//from apop_t_f_chi.c:
APOP_VAR_DECLARE int apop_matrix_is_positive_semidefinite(gsl_matrix *m, char semi);
double apop_matrix_to_positive_semidefinite(gsl_matrix *m);
double apop_multivariate_gamma(double a, double p);
double apop_multivariate_lngamma(double a, double p);

//from apop_sort.c:
APOP_VAR_DECLARE double * apop_vector_percentiles(gsl_vector *data, char rounding); 
APOP_VAR_DECLARE apop_data * apop_data_sort(apop_data *data, int sortby, char asc);

//from the regression code:
#define apop_F_test apop_f_test

apop_data *	apop_t_test(gsl_vector *a, gsl_vector *b);
apop_data *	apop_paired_t_test(gsl_vector *a, gsl_vector *b);


apop_data * apop_text_unique_elements(const apop_data *d, size_t col);
gsl_vector * apop_vector_unique_elements(const gsl_vector *v);
apop_data *apop_text_to_factors(apop_data *d, size_t textcol, int datacol);

APOP_VAR_DECLARE apop_data * apop_data_to_dummies(apop_data *d, int col, char type, int keep_first);
APOP_VAR_DECLARE apop_data * apop_f_test (apop_model *est, apop_data *contrast);

apop_model *apop_estimate_fixed_effects_OLS(apop_data *data, gsl_vector *categories);

apop_data *apop_estimate_coefficient_of_determination (apop_model *in);
apop_data *apop_estimate_r_squared (apop_model *in);
void apop_estimate_parameter_t_tests (apop_model *est);

//apop_testing.c
apop_data* apop_anova(char *table, char *data, char *grouping1, char *grouping2);

#define apop_ANOVA(table, data, grouping1, grouping2) apop_anova(table, data, grouping1, grouping2)

__END_DECLS
#endif
