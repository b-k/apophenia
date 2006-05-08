/** \file apop_estimate.c	 sets up the estimate structure which outputs from the various regressions and MLEs.


Copyright (c) 2006 by Ben Klemens. Licensed under the GNU GPL v2.
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

It takes in all the inputs to the model: data, 
parameters, plus the model itself. Those elements of the inventory
that are valid are copied in, the parameters are copied, and you also get
pointers to the data and model. Those pointers are  primarily
there for your reference. They are just pointers to the original data,
so if you destroy your data eleswhere, these pointers will faithfully
point to garbage.

Also, the parameters, and anything else specified in the inventory, is allocated and ready for you to use.

\param data	        A pointer to an input apop_data set
\param model	    A pointer to the model you're estimating
\param  params      An \ref apop_estimation_params structure. May be NULL. Don't forget that this may include an apop_inventory (or that field may be NULL)

\ingroup inv_and_est  */
apop_estimate * apop_estimate_alloc(apop_data * data, apop_model model, apop_estimation_params *params){
apop_estimate * prep_me;
	prep_me	= malloc(sizeof(apop_estimate));
    if (params){
        memcpy(&(prep_me->estimation_params), params, sizeof(apop_estimation_params));
	   // apop_inventory_copy(apop_inventory_filter(&(params->uses), model.inventory_filter), &(prep_me->estimation_params.uses));
    } else {
        apop_estimation_params *delme =  apop_estimation_params_alloc();
        memcpy(&(prep_me->estimation_params), delme, sizeof(apop_estimation_params));
        apop_estimation_params_free(delme);
	  //  apop_inventory_filter(&(prep_me->estimation_params.uses), model.inventory_filter);
    }
	if (prep_me->estimation_params.uses.parameters)
		prep_me->parameters	= apop_data_alloc(model.parameter_ct,-1);
	if (prep_me->estimation_params.uses.dependent ||
	                prep_me->estimation_params.uses.predicted){
        if (data && data->matrix)
		    prep_me->dependent	= apop_data_alloc(data->matrix->size1,3);
        else if (data && data->vector)
		    prep_me->dependent	= apop_data_alloc(data->vector->size,3);
        else if (data && data->categories)
		    prep_me->dependent	= apop_data_alloc(data->catsize[0],3);
        else
		    prep_me->dependent	= NULL;
        if (prep_me->dependent){
            apop_name_add(prep_me->dependent->names, "actual", 'c');
            apop_name_add(prep_me->dependent->names, "predicted", 'c');
            apop_name_add(prep_me->dependent->names, "residual", 'c');
        }
        if (data && data->names && data->names->rownamect > 0)
            apop_name_stack(prep_me->dependent->names, data->names, 'r');
    }
	if (prep_me->estimation_params.uses.covariance){
		prep_me->covariance	= apop_data_alloc(model.parameter_ct,model.parameter_ct);
        if (data && data->names){
            apop_name_stack(prep_me->covariance->names, data->names, 'c');
            apop_name_cross_stack(prep_me->covariance->names, data->names, 'c', 'r');
        }
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
	if (free_me->estimation_params.uses.covariance)
		apop_data_free(free_me->covariance);
	if (free_me->estimation_params.uses.parameters)
        apop_data_free(free_me->parameters);
	free(free_me);
}

/** Print the results of an estimation

\ingroup output */
void apop_estimate_print(apop_estimate * print_me){
	if (print_me->estimation_params.uses.parameters)	
        apop_data_show(print_me->parameters);
	if (print_me->estimation_params.uses.covariance){
		printf("\nThe variance/covariance matrix:\n");
        apop_data_show(print_me->covariance);
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
	out->estimate           = in.estimate;
	out->log_likelihood     = in.log_likelihood;
	out->dlog_likelihood    = in.dlog_likelihood;
	out->fdf                = in.fdf;
    out->constraint         = in.constraint;
	out->rng                = in.rng;
    return out;
}


/** Neatly allocate an \ref apop_estimation_params structure. Sets a
few defaults, so you can change just one or two values and everything
else will be predictable.

 */
apop_estimation_params *apop_estimation_params_alloc(){
apop_estimation_params *setme = calloc(sizeof(apop_estimation_params),1);
    setme->starting_pt  = NULL;
    setme->step_size    = 1;
    apop_inventory_set(&(setme->uses),1); 
    return setme;
}

/** Neatly allocate an \ref apop_estimation_params structure. Sets a
few defaults, so you can change just one or two values and everything
else will be predictable.

 */
void apop_estimation_params_free(apop_estimation_params *freeme){
    free(freeme);
}
