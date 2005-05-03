//estimate.c			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
#include <gsl/gsl_matrix.h>
#include "estimate.h"
#include "linear_algebra.h"

void apop_copy_inventory(apop_inventory in, apop_inventory *out){
	out->parameters	= in.parameters;
	out->predicted	= in.predicted;
	out->confidence	= in.confidence;
	out->covariance	= in.covariance;
	out->log_likelihood = in.log_likelihood;
}

void apop_set_inventory(apop_inventory *out, int value){
	out->parameters	= 
	out->predicted	= 
	out->confidence	= 
	out->covariance	= 
	out->log_likelihood = value;
}

apop_estimate * apop_estimate_alloc(int data_size, int param_size, apop_inventory uses){
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
	apop_copy_inventory(uses, &(prep_me->uses));
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
	free(free_me);
}

void apop_print_estimate(apop_estimate * print_me, char *out){
	if (print_me->uses.parameters){
		apop_print_to_file(out, "Parameter estimates:\t");
		apop_print_vector(print_me->parameters, "\t", out);
	}
	if (print_me->uses.covariance){
		apop_print_to_file(out, "The variance/covariance matrix:\n");
		apop_print_matrix(print_me->covariance, "\t", out);
	}
	if (print_me->uses.confidence){
		apop_print_to_file(out, "Confidence intervals (H_0: beta == 0):\t");
		apop_print_vector(print_me->confidence, "\t", out);
	}
	if (print_me->uses.log_likelihood)
		apop_print_to_file(out, "log likelihood: \t%g\n", print_me->log_likelihood);
}
