//regression.h			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.

#include <gsl/gsl_matrix.h>
#include "model.h"

apop_model * apop_OLS(gsl_matrix *data, int only_want_parameter_estimates);
apop_model * apop_GLS(gsl_matrix *data, gsl_matrix *sigma, int only_want_parameter_estimates);

//Returns GLS/OLS parameter estimates.
//Destroys the data in the process.
