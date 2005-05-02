//model.c			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
#include <gsl/gsl_matrix.h>
#include "model.h"
#include "linear_algebra.h"


apop_model * apop_model_alloc(int size){
apop_model * prep_me;
	prep_me	= malloc(sizeof(apop_model));
	prep_me->parameters	= gsl_vector_alloc(size);
	prep_me->confidence	= gsl_vector_alloc(size);
	prep_me->covariance	= gsl_matrix_alloc(size,size);
	prep_me->params		= prep_me->parameters;
	prep_me->cov		= prep_me->covariance;
	prep_me->size		= size;
	return prep_me;
}

void apop_model_free(apop_model * free_me){
	if (free_me->parameters != NULL)
		gsl_vector_free(free_me->parameters);
	if (free_me->confidence != NULL)
		gsl_vector_free(free_me->confidence);
	if (free_me->covariance != NULL)
		gsl_matrix_free(free_me->covariance);
	free(free_me);
}

void apop_print_model(apop_model * print_me, char *out){
	printf("Parameter estimates:\n");
	apop_print_vector(print_me->parameters, "\t", out);
	printf("The variance/covariance matrix:\n");
	apop_print_matrix(print_me->covariance, "\t", out);
	printf("Confidence intervals (H_0: beta == 0):\n");
	apop_print_vector(print_me->confidence, "\t", out);
}
