/** \file apop_estimate.c	 sets up the estimate structure which outputs from the various regressions and MLEs.


Copyright (c) 2006 by Ben Klemens. Licensed under the GNU GPL v2.
*/

#include <gsl/gsl_matrix.h>
#include "apophenia/output.h"
#include "apophenia/linear_algebra.h"

/** \defgroup inv_and_est  Using estimates (and their parameters)

The \ref apop_estimate structure is the output for every regression or maximum likelihood technique here.
It gathers together the parameter estimates, the variance-covariance matrix of the parameter
estimates, the likelihood of the data if a max. likelihood technique,
and whatever else comes to mind. 

See the page on \ref types for more.

\ingroup types
 */

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
\param  params      An \ref apop_ep structure. May be NULL. 

\ingroup inv_and_est  */
apop_estimate * apop_estimate_alloc(apop_data * data, apop_model model, apop_ep *params){
apop_estimate * prep_me;
	prep_me	= malloc(sizeof(apop_estimate));
    if (params){
        memcpy(&(prep_me->ep), params, sizeof(apop_ep));
    } else {
        apop_ep *delme =  apop_ep_alloc();
        memcpy(&(prep_me->ep), delme, sizeof(apop_ep));
        apop_ep_free(delme);
    }
	if (prep_me->ep.uses.parameters)
		prep_me->parameters	= apop_data_alloc(model.vsize, model.msize1, model.msize2);
	if (prep_me->ep.uses.dependent ||
	                prep_me->ep.uses.predicted){
        if (data && data->matrix)
		    prep_me->dependent	= apop_data_alloc(0, data->matrix->size1,3);
        else if (data && data->vector)
		    prep_me->dependent	= apop_data_alloc(0, data->vector->size,3);
        else if (data && data->categories)
		    prep_me->dependent	= apop_data_alloc(0, data->catsize[0],3);
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
	if (prep_me->ep.uses.covariance){
		prep_me->covariance	= apop_data_alloc(0, model.vsize,model.vsize);
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
	if (free_me->ep.uses.predicted
	    || free_me->ep.uses.dependent)
        if (free_me->dependent)
		    apop_data_free(free_me->dependent);
	if (free_me->ep.uses.covariance)
		apop_data_free(free_me->covariance);
	if (free_me->ep.uses.parameters)
        apop_data_free(free_me->parameters);
	free(free_me);
}

/** Print the results of an estimation

\ingroup output */
void apop_estimate_show(apop_estimate * print_me){
	if (print_me->ep.uses.parameters)
        apop_data_show(print_me->parameters);
	if (print_me->ep.uses.covariance){
		printf("\nThe variance/covariance matrix:\n");
        apop_data_show(print_me->covariance);
	}
	if (print_me->ep.uses.log_likelihood)
		printf("\nlog likelihood: \t%g\n", print_me->log_likelihood);
}

/** Currently an alias for \ref apop_estimate_show, but when I get
  around to it, it will conform better with the other apop_..._print
  fns.*/
void apop_estimate_print(apop_estimate * print_me){
    apop_estimate_show(print_me);
}

/** Outputs a copy of the \ref apop_model input.
\param in   the model to be copied
\return a pointer to a copy of the original, which you can mangle as you see fit. 
*/
apop_model * apop_model_copy(apop_model in){
  apop_model * out = malloc(sizeof(apop_model));
    strcpy(out->name, in.name);
    out->vsize              = in.vsize;
    out->msize1             = in.msize1;
    out->msize2             = in.msize2;
	out->estimate           = in.estimate;
	out->p                  = in.p;
	out->log_likelihood     = in.log_likelihood;
	out->score              = in.score;
    out->constraint         = in.constraint;
	out->draw               = in.draw;
    return out;
}


/** Neatly allocate an \ref apop_ep structure. Sets a
few defaults, so you can change just one or two values and everything
else will be predictable.

 */
apop_ep *apop_ep_alloc(){
  apop_ep *setme = calloc(sizeof(apop_ep),1);
    setme->starting_pt          = NULL;
    setme->weights              = NULL;
    setme->more                 = NULL;
    setme->step_size            = 1;
    memset(&(setme->uses),1, sizeof(setme->uses)); 
    return setme;
}

/** Neatly allocate an \ref apop_ep structure. Sets a
few defaults, so you can change just one or two values and everything
else will be predictable.

 */
void apop_ep_free(apop_ep *freeme){
    free(freeme);
}
