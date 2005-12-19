/** \file apop_OLS.c

  OLS models. Much of the real work is done in regression.c.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/

#include "model.h"
#include "regression.h"
#include <assert.h>


/** The OLS model

  This is basically a wrapper for the OLS regression function, \ref apop_estimate_OLS.
\ingroup models
*/
apop_model apop_OLS = {"OLS", -1, 
{
	1,	//parameters
	1,	//covariance
	1,	//confidence
	0,	//predicted
	1,	//residuals
	0,	//log_likelihood
	1	//names;
},	 
	apop_estimate_OLS, NULL, NULL, NULL, NULL, NULL};

/** The GLS model

  This is basically a wrapper for the GLS regression function, \ref apop_estimate_GLS.
\ingroup models
*/
apop_model apop_GLS = {"GLS", -1, 
{
	1,	//parameters
	1,	//covariance
	1,	//confidence
	0,	//predicted
	1,	//residuals
	0,	//log_likelihood
	1	//names;
},	 
	apop_estimate_GLS, NULL, NULL, NULL, NULL, NULL};
