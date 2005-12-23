/** \file apop_estimate.c	 sets up the estimate structure which outputs from the various regressions and MLEs.


Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL.
*/

#include <gsl/gsl_matrix.h>
#include "apophenia/output.h"
#include "apophenia/linear_algebra.h"

/** \defgroup inv_and_est  Using inventories and estimates 

The \ref apop_estimate structure is the output for every regression or maximum likelihood technique here.
It gathers together the parameter estimates, the variance-covariance matrix of the parameter
estimates, the likelihood of the data if a max. likelihood technique,
and whatever else comes to mind.

\todo Rewrite the inventory and estimate documentation to cohere better. 
 */

/** Allocate an \ref apop_inventory and set all of its elements to some
value. If the inventory is already allocated, use \ref apop_inventory_set.

\param value	probably one or zero.  
\return A ready-to-use inventory
 *
 \ingroup inv_and_est 
 */
apop_inventory * apop_inventory_alloc(int value){
apop_inventory  *setme;
    setme  = malloc(sizeof(apop_inventory));
    apop_inventory_set(setme, value);
    return setme;
}

/** Copy one inventory to another.
 \ingroup inv_and_est 
 */
void apop_inventory_copy(apop_inventory in, apop_inventory *out){
	out->parameters	= in.parameters;
	out->covariance	= in.covariance;
	out->confidence	= in.confidence;
	out->predicted	= in.predicted;
	out->residuals	= in.residuals;
	out->names  	= in.names;
	out->log_likelihood = in.log_likelihood;
}

/** set all the values of the inventory to a certain value. Like 
 \ref apop_inventory_alloc but doesn't call malloc().

\param out	a pointer to the inventory to be set
\param value	probably one or zero.  
\ingroup inv_and_est */
void apop_inventory_set(apop_inventory *out, int value){
	out->parameters	= 
	out->predicted	= 
	out->confidence	= 
	out->covariance	= 
	out->residuals	= 
	out->names  	= 
	out->log_likelihood = value;
}

/** Apply a filter to an input inventory.

Inventories sent in to most estimation functions are just the wish list; it wouldn't make sense (and is often not implemented) that every estimator return every element of the \ref apop_estimate.

\param out	a pointer to the inventory desired. If null, the filter will be copied to this spot.
\param filter	an \ref apop_inventory where the values are one if the value will be calculated, zero if not.
\ingroup inv_and_est */
apop_inventory apop_inventory_filter(apop_inventory *in, apop_inventory filter){
apop_inventory  out;
	if (in==NULL){
		apop_inventory_copy(filter, &out);
        return out;
    }// else:
	apop_inventory_copy(*in, &out);
	out.parameters	    &= filter.parameters;
	out.covariance	    &= filter.covariance;
	out.confidence	    &= filter.confidence;
	out.predicted	    &= filter.predicted;
	out.residuals	    &= filter.residuals;
	out.log_likelihood &= filter.log_likelihood;
	out.names          &= filter.names;
    return out;
}

/** Allocate an \ref apop_estimate.

Are you sure you need to use this? Every regression and MLE function
provided by Apophenia creates this for you. [You will, however, have to
free the estimate yourself.]

\param data_size	If you want predicted values or residuals, we need to know how many slots to declare.
\param	param_size	Notice that in the typical data matrix (not always), data_size is data_matrix->size1 and param_size is data_matrix->size2
\param n		name structure, if you have one
\param uses		Tell me what you want declared.
\ingroup inv_and_est  */
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

/** Free an \ref apop_estimate.

\ingroup inv_and_est */
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

/** Print the results of an estimation

\ingroup output */
void apop_estimate_print(apop_estimate * print_me){
int		i;
	printf("\n");
	if (print_me->uses.names) 	printf("\t");
	if (print_me->uses.parameters) 	printf("value\t\t");
	if (print_me->uses.confidence) 	
            printf("Confidence\n");
    else    printf("\n");
	for (i=0; i<print_me->parameters->size; i++){
		if (print_me->uses.names)	printf("%s\t", print_me->names->colnames[i]);
		if (print_me->uses.parameters)	printf("% 7f\t", gsl_vector_get(print_me->parameters,i));
		if (print_me->uses.confidence)	printf("% 7f\t", gsl_vector_get(print_me->confidence,i));
		printf("\n");
		/*
		if (print_me->uses.parameters){
			printf("Parameter estimates:\t");
			apop_print_vector(print_me->parameters, "\t", NULL);
		}
		if (print_me->uses.confidence){
			printf("Confidence intervals (H_0: beta == 0):\t");
			apop_print_vector(print_me->confidence, "\t", NULL);
		}
		*/
	}
	if (print_me->uses.covariance){
		printf("\nThe variance/covariance matrix:\n");
        apop_data   *covdata    = apop_matrix_to_data(print_me->covariance);
        //We want to show the column names on both axes.
        if (print_me->uses.names && print_me->names !=NULL && print_me->names->colnamect){
            free(covdata->names->colnames);    //prevent a 2-byte memory leak here.
            free(covdata->names->rownames);    
            covdata->names->colnames    = malloc(sizeof(char*) * print_me->names->colnamect);
            covdata->names->rownames    = malloc(sizeof(char*) * print_me->names->colnamect);
            memcpy(covdata->names->colnames, print_me->names->colnames, sizeof(char*) * print_me->names->colnamect);
            memcpy(covdata->names->rownames, print_me->names->colnames, sizeof(char*) * print_me->names->colnamect);
            covdata->names->rownamect   =
            covdata->names->colnamect   =   print_me->names->colnamect;
        }
        apop_data_print(covdata, "\t", NULL);
	}
	if (print_me->uses.log_likelihood)
		printf("\nlog likelihood: \t%g\n", print_me->log_likelihood);
}
