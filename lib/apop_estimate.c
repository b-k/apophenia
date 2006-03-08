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
\ingroup types
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
	out->dependent	= in.dependent;
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
	out->dependent	= 
	out->names  	= 
	out->log_likelihood = value;
}

/** Apply a filter to an input inventory.

Inventories sent in to most estimation functions are just the wish list; it wouldn't make sense (and is often not implemented) that every estimator return every element of the \ref apop_estimate.

\param in   	a pointer to the inventory desired. If null, the filter will be copied to this spot.
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
	out.dependent	    &= filter.dependent;
	out.log_likelihood &= filter.log_likelihood;
	out.names          &= filter.names;
    return out;
}

/** Allocate an \ref apop_estimate.

Are you sure you need to use this? Every regression and MLE function
provided by Apophenia creates this for you. [You will, however, have to
free the estimate yourself.]

It takes in all the inputs to the model: data, inventory, parameters, plus the model itself. Those elements of the inventory that are valid are copied in, and you also get pointers to the data, model, and parameters. Those pointers are never really used by Apophenia, but are there for your reference. They're just pointers, so if you destroy your data eleswhere, these pointers will faithfully point to garbage.

\param data	        A pointer to an input apop_data set
\param model	    A pointer to the model you're estimating
\param  params      An \ref apop_estimation_params structure. May be NULL. Don't forget that this may include an apop_inventory (or that field may be NULL)

\ingroup inv_and_est  */
apop_estimate * apop_estimate_alloc(apop_data * data, apop_model model, apop_estimation_params *params){
apop_estimate * prep_me;
	prep_me	= malloc(sizeof(apop_estimate));
    if (params){
        memcpy(&(prep_me->estimation_params), params, sizeof(apop_estimation_params));
	    apop_inventory_copy(apop_inventory_filter(&(params->uses), model.inventory_filter), &(prep_me->estimation_params.uses));
    } else {
        apop_estimation_params *delme =  apop_estimation_params_alloc();
        memcpy(&(prep_me->estimation_params), delme, sizeof(apop_estimation_params));
	    apop_inventory_filter(&(prep_me->estimation_params.uses), model.inventory_filter);
        free(delme);
    }
	if (prep_me->estimation_params.uses.parameters)
		prep_me->parameters	= gsl_vector_alloc(model.parameter_ct);
	if (prep_me->estimation_params.uses.confidence)
		prep_me->confidence	= gsl_vector_alloc(model.parameter_ct);
	if (prep_me->estimation_params.uses.dependent ||
	                prep_me->estimation_params.uses.predicted){
		prep_me->dependent	= apop_data_alloc(data->data->size1,3);
        apop_name_add(prep_me->dependent->names, "actual", 'c');
        apop_name_add(prep_me->dependent->names, "predicted", 'c');
        apop_name_add(prep_me->dependent->names, "residual", 'c');
        apop_name_stack(prep_me->dependent->names, data->names, 'r');
    }
	if (prep_me->estimation_params.uses.covariance){
		prep_me->covariance	= apop_data_alloc(model.parameter_ct,model.parameter_ct);
        apop_name_stack(prep_me->covariance->names, data->names, 'c');
        apop_name_cross_stack(prep_me->covariance->names, data->names, 'c', 'r');
        }
	if (prep_me->estimation_params.uses.names) {
		if (data->names != NULL) 	apop_name_memcpy(&(prep_me->names), data->names);
		else 		                prep_me->names		= apop_name_alloc();
	}
    prep_me->data               = data;
    prep_me->model              = apop_model_copy(model);
	return prep_me;
}

/** Free an \ref apop_estimate.

\ingroup inv_and_est */
void apop_estimate_free(apop_estimate * free_me){
	if (free_me->estimation_params.uses.predicted
	    || free_me->estimation_params.uses.dependent)
		apop_data_free(free_me->dependent);
	if (free_me->estimation_params.uses.confidence)
		gsl_vector_free(free_me->confidence);
	if (free_me->estimation_params.uses.covariance)
		apop_data_free(free_me->covariance);
	if (free_me->estimation_params.uses.names)
		apop_name_free(free_me->names);
	free(free_me);
}

/** Print the results of an estimation

\ingroup output */
void apop_estimate_print(apop_estimate * print_me){
int		i;
	printf("\n");
	if (print_me->estimation_params.uses.names) 	printf("\t");
	if (print_me->estimation_params.uses.parameters) 	printf("value\t\t");
	if (print_me->estimation_params.uses.confidence) 	
            printf("Confidence\n");
    else    printf("\n");
	for (i=0; i<print_me->parameters->size; i++){
		if (print_me->estimation_params.uses.names)	printf("%s\t", print_me->names->colnames[i]);
		if (print_me->estimation_params.uses.parameters)	printf("% 7f\t", gsl_vector_get(print_me->parameters,i));
		if (print_me->estimation_params.uses.confidence)	printf("% 7f\t", gsl_vector_get(print_me->confidence,i));
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
	if (print_me->estimation_params.uses.covariance){
		printf("\nThe variance/covariance matrix:\n");
        apop_data_show(print_me->covariance);
        /*
        apop_data   *covdata    = apop_matrix_to_data(print_me->covariance);
        int         i;
        //We want to show the column names on both axes.
        if (print_me->estimation_params.uses.names && print_me->names !=NULL && print_me->names->colnamect){
            apop_name_stack(covdata->names, print_me->names, 'c');
            for (i=0; i< print_me->names->colnamect; i++)
                apop_name_add(covdata->names, print_me->names->colnames[i], 'r');
        }
        apop_data_show(covdata);
        free(covdata);  //slightly leaky. I don't care.
        */
	}
	if (print_me->estimation_params.uses.log_likelihood)
		printf("\nlog likelihood: \t%g\n", print_me->log_likelihood);
}

/** Outputs a copy of the \ref apop_model input.
\param in   the model to be copied
\return a pointer to a copy of the original, which you can mangle as you see fit. 
*/
apop_model * apop_model_copy(apop_model in){
apop_model * out = malloc(sizeof(apop_model));
    strcpy(out->name, in.name);
    out->parameter_ct   = in.parameter_ct;
    memcpy(&(out->inventory_filter), &(in.inventory_filter), sizeof(apop_inventory));
	out->estimate           = in.estimate;
	out->log_likelihood     = in.log_likelihood;
	out->dlog_likelihood    = in.dlog_likelihood;
	out->fdf                = in.fdf;
    out->constraint         = in.constraint;
	out->rng                = in.rng;
    return out;
}


apop_estimation_params *apop_estimation_params_alloc(){
apop_estimation_params *setme = calloc(sizeof(apop_estimation_params),1);
    apop_inventory_set(&(setme->uses),1);
    return setme;
}
