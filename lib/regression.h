//regression.h			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.

#include <gsl/gsl_matrix.h>
#include "estimate.h"

apop_estimate * apop_OLS(gsl_matrix *data, apop_name *n, apop_inventory *uses);
apop_estimate * apop_GLS(gsl_matrix *data, gsl_matrix *sigma, apop_name *n, apop_inventory *uses);
//Returns GLS/OLS parameter estimates.
//Destroys the data in the process.

double	apop_t_test(gsl_vector *a, gsl_vector *b);
double	apop_paired_t_test(gsl_vector *a, gsl_vector *b);
//A nice, easy t test. With what confidence can we reject the hypothesis
//that the mean of vector A equals the mean of vector B?



double two_tailify(double in);
//My convenience fn to turn the results from a symmetric one-tailed table lookup
//into a two-tailed confidence interval.
