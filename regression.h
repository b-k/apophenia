//regression.h	Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.

#ifndef apop_regression_h
#define  apop_regression_h

#include <gsl/gsl_matrix.h>
#include <apophenia/types.h>

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

typedef struct {
    int destroy_data;
    gsl_vector *weights;
    apop_model *model;
    apop_data *instruments;
    int want_cov;
    int want_expected_value;
} apop_ls_settings;

apop_ls_settings * apop_ls_settings_alloc(apop_data *data, apop_model model);
apop_model * apop_estimate_OLS(apop_data *set, apop_model *ep);
apop_model * apop_estimate_GLS(apop_data *set, gsl_matrix *sigma);
apop_model *apop_fixed_effects_OLS(apop_data *data, gsl_vector *categories);
//Returns GLS/OLS parameter estimates.
//Destroys the data in the process.

apop_data *apop_F_test (apop_model *est, apop_data *contrast);
apop_data *apop_f_test (apop_model *est, apop_data *contrast);

apop_data *	apop_t_test(gsl_vector *a, gsl_vector *b);
apop_data *	apop_paired_t_test(gsl_vector *a, gsl_vector *b);

apop_data * apop_data_to_dummies(apop_data *d, int col, char type, int keep_first);

double apop_two_tailify(double in);
//My convenience fn to turn the results from a symmetric one-tailed table lookup
//into a two-tailed confidence interval.

apop_model *apop_estimate_fixed_effects_OLS(apop_data *data, gsl_vector *categories);

apop_data *apop_estimate_correlation_coefficient (apop_model *in);
apop_data *apop_estimate_r_squared (apop_model *in);
void apop_estimate_parameter_t_tests (apop_model *est);

apop_model * apop_estimate_IV(apop_data *inset, apop_model *ep);

//apop_testing.c
apop_data* apop_ANOVA(char *table, char *data, char *grouping1, char *grouping2);

__END_DECLS
#endif
