//estimate.c			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
#include <gsl/gsl_matrix.h>
#include "name.h"
#include "estimate.h"
#include "linear_algebra.h"

void apop_inventory_copy(apop_inventory in, apop_inventory *out){
	out->parameters	= in.parameters;
	out->predicted	= in.predicted;
	out->confidence	= in.confidence;
	out->covariance	= in.covariance;
	out->log_likelihood = in.log_likelihood;
}

void apop_inventory_set(apop_inventory *out, int value){
	out->parameters	= 
	out->predicted	= 
	out->confidence	= 
	out->covariance	= 
	out->log_likelihood = value;
}

apop_estimate * apop_estimate_alloc(int data_size, int param_size, apop_name * n, apop_inventory uses){
apop_estimate * prep_me;
	prep_me	= malloc(sizeof(apop_estimate));
	if (uses.parameters)
		prep_me->parameters	= gsl_vector_alloc(param_size);
	if (uses.confidence)
		prep_me->confidence	= gsl_vector_alloc(param_size);
	if (uses.predicted)
		prep_me->predicted	= gsl_vector_alloc(data_size);
	if (uses.residuals)
		prep_me->residuals	= gsl_vector_alloc(data_size);
	if (uses.covariance)
		prep_me->covariance	= gsl_matrix_alloc(param_size,param_size);
	if (uses.names) {
		if (n != NULL) 	prep_me->names		= n;
		else 		prep_me->names		= apop_name_alloc();
	}
	apop_inventory_copy(uses, &(prep_me->uses));
	return prep_me;
}

void apop_estimate_free(apop_estimate * free_me){
	if (free_me->uses.predicted)
		gsl_vector_free(free_me->predicted);
	if (free_me->uses.predicted)
		gsl_vector_free(free_me->predicted);
	if (free_me->uses.residuals)
		gsl_vector_free(free_me->residuals);
	if (free_me->uses.confidence)
		gsl_vector_free(free_me->confidence);
	if (free_me->uses.covariance)
		gsl_matrix_free(free_me->covariance);
	if (free_me->uses.names)
		apop_name_free(free_me->names);
	free(free_me);
}

void apop_estimate_print(apop_estimate * print_me){
	if (print_me->uses.names){
		apop_name_print(print_me->names);
	}
	if (print_me->uses.parameters){
		printf("Parameter estimates:\t");
		apop_print_vector(print_me->parameters, "\t", NULL);
	}
	if (print_me->uses.covariance){
		printf("The variance/covariance matrix:\n");
		apop_print_matrix(print_me->covariance, "\t", NULL);
	}
	if (print_me->uses.confidence){
		printf("Confidence intervals (H_0: beta == 0):\t");
		apop_print_vector(print_me->confidence, "\t", NULL);
	}
	if (print_me->uses.log_likelihood)
		printf("log likelihood: \t%g\n", print_me->log_likelihood);
}
