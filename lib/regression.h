//regression.h			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.

#include <gsl/gsl_matrix.h>
#include <apophenia/types.h>

apop_estimate * apop_estimate_OLS(apop_data *set, apop_inventory *uses, void *dummy);
apop_estimate * apop_estimate_GLS(apop_data *set, apop_inventory *uses, gsl_matrix *sigma);
apop_estimate *apop_fixed_effects_OLS(apop_data *data, apop_inventory *uses, gsl_vector *categories);
//Returns GLS/OLS parameter estimates.
//Destroys the data in the process.

double	apop_t_test(gsl_vector *a, gsl_vector *b);
double	apop_paired_t_test(gsl_vector *a, gsl_vector *b);
//A nice, easy t test. With what confidence can we reject the hypothesis
//that the mean of vector A equals the mean of vector B?

apop_data * apop_produce_dummies(gsl_vector *in, int keep_first);

double two_tailify(double in);
//My convenience fn to turn the results from a symmetric one-tailed table lookup
//into a two-tailed confidence interval.

apop_estimate *apop_estimate_fixed_effects_OLS(apop_data *data, apop_inventory *uses, gsl_vector *categories);
