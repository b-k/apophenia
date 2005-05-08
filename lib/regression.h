//regression.h			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.

#include <gsl/gsl_matrix.h>
#include "estimate.h"

apop_estimate * apop_OLS(gsl_matrix *data, apop_name *n, apop_inventory *uses);
apop_estimate * apop_GLS(gsl_matrix *data, gsl_matrix *sigma, apop_name *n, apop_inventory *uses);

//Returns GLS/OLS parameter estimates.
//Destroys the data in the process.
