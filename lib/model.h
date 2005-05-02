//model.h			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
#ifndef __apop_model__
#define __apop_model__

#include <gsl/gsl_matrix.h>

typedef struct apop_m{
	gsl_vector 	*parameters, *params, *confidence;
	gsl_matrix 	*covariance, *cov;
	int		size;
} apop_model;

apop_model * 	apop_model_alloc(int size);
void 		apop_model_free(apop_model * free_me);
void 		apop_print_model(apop_model * print_me, char *out);
#endif
