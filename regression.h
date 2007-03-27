//regression.h			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.

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

//apop_estimate * apop_estimate_OLS(apop_data *set, apop_ep *ep);
apop_estimate * apop_estimate_OLS(apop_data *inset, void *epin);
apop_estimate * apop_estimate_GLS(apop_data *set, gsl_matrix *sigma);
apop_estimate *apop_fixed_effects_OLS(apop_data *data, gsl_vector *categories);
//Returns GLS/OLS parameter estimates.
//Destroys the data in the process.

apop_data *apop_F_test(apop_estimate *est, apop_data *contrast);
apop_data *apop_f_test(apop_estimate *est, apop_data *contrast);

apop_data *	apop_t_test(gsl_vector *a, gsl_vector *b);
apop_data *	apop_paired_t_test(gsl_vector *a, gsl_vector *b);

apop_data * apop_produce_dummies(gsl_vector *in, int keep_first);

double apop_two_tailify(double in);
//My convenience fn to turn the results from a symmetric one-tailed table lookup
//into a two-tailed confidence interval.

apop_estimate *apop_estimate_fixed_effects_OLS(apop_data *data, gsl_vector *categories);

apop_data *apop_estimate_correlation_coefficient(apop_estimate *in);
apop_data *apop_estimate_r_squared(apop_estimate *in);
void apop_estimate_parameter_t_tests(apop_estimate *est);

__END_DECLS
