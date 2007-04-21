/** \file apop_params.c	 sets up the estimate structure which outputs from the various regressions and MLEs.


Copyright (c) 2006 by Ben Klemens. Licensed under the GNU GPL v2.
*/

#include <gsl/gsl_matrix.h>
#include "apophenia/output.h"
#include "apophenia/linear_algebra.h"

/** \defgroup inv_and_est  Using \ci apop_paramss 

The data to accompany an \c apop_model, including the input settings and the output parameters, expected values, et cetera.

See the page on \ref types for more.

\ingroup types
 */

/** Allocate an \ref apop_params.

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

\ingroup inv_and_est  */
apop_params * apop_params_alloc(apop_data * data, apop_model *model, void *method_params, void *model_params){
apop_params * prep_me = malloc(sizeof(apop_params));
    /* This has some cruft from when there were separate apop_ep and
       apop_params structs. Harmless, so I didn't mess with it.*/
    int vsize  = model->vbase == -1 ? data->matrix->size2 : model->vbase;
    int msize1 = model->m1base == -1 ? data->matrix->size2 : model->m1base ;
    int msize2 = model->m2base == -1 ? data->matrix->size2 : model->m2base ;
    prep_me->parameters	= apop_data_alloc(vsize, msize1, msize2);
    memset(&(prep_me->uses), 1, sizeof(prep_me->uses));
	//if (prep_me->uses.expected ||
	 //               prep_me->uses.predicted){
        if (data && data->matrix)
		    prep_me->expected	= apop_data_alloc(0, data->matrix->size1,3);
        else if (data && data->vector)
		    prep_me->expected	= apop_data_alloc(0, data->vector->size,3);
        else if (data && data->categories)
		    prep_me->expected	= apop_data_alloc(0, data->catsize[0],3);
        else
		    prep_me->expected	= NULL;
        if (prep_me->expected){
            apop_name_add(prep_me->expected->names, "actual", 'c');
            apop_name_add(prep_me->expected->names, "predicted", 'c');
            apop_name_add(prep_me->expected->names, "residual", 'c');
            if (data && data->names && data->names->rownamect > 0)
                apop_name_stack(prep_me->expected->names, data->names, 'r');
        }
    //}
	if (prep_me->uses.covariance){
		prep_me->covariance	= apop_data_alloc(0, vsize,vsize);
        if (data && data->names){
            apop_name_stack(prep_me->covariance->names, data->names, 'c');
            apop_name_cross_stack(prep_me->covariance->names, data->names, 'c', 'r');
        }
    } else
		prep_me->covariance	= NULL;
    prep_me->data               = data;
    prep_me->model              = model;
    prep_me->method_name[0]     = '\0';
    prep_me->method_params      = method_params;
    prep_me->model_params       = model_params;
    prep_me->status             = 0;
	return prep_me;
}

/** Sometimes, you have parameters in mind, and just want to turn them
  into parameters as quickly as possible. 

  \param params An \c apop_data struct with the parameters. It will be copied in to the \c apop_params structure.
  \return An \c apop\_params structure that is basically empty, except for the \c parameters element. This should already be enough for most \c log_likelihood or \c draw methods.
 
 */
apop_params *apop_params_alloc_p(apop_data *params){
  apop_params * prep_me = malloc(sizeof(apop_params));
    prep_me->parameters = apop_data_copy(params);
    prep_me->method_name[0] = '\0';
    prep_me->data           =
    prep_me->expected       =
    prep_me->covariance     = NULL;
    prep_me->model_params   =
    prep_me->method_params  = NULL;
    return prep_me;
}


/** Copy an \c apop_params structure.


The structure includes several pointers, so this
function just returns a new bundle of pointers pointing at the same data.
Changes to the copy's elements affect the original until you change individual pointers to point to new locations.

See also \c apop_params_complete_copy, which also copies most of the internal pointers.
  
  \param in The apop_param to copy
  \return  A pointer to a new copy.
\ingroup inv_and_est  */
apop_params *apop_params_copy(apop_params *in){
    if (!in) return NULL;
  apop_params *out  = malloc(sizeof(apop_params));
    memcpy(out, &in, sizeof(apop_params));
    return out;
}

/** Make a copy an \c apop_params structure, and most of its constituent elements.

This copies all the static elements, does a \c memcpy of the
method_params, model_params, more, parameters, expected, and covariance
elements. It does not copy the model, under the assumption that the model
is static (use \c apop_model_copy if you need), and it does not copy the
data, under the assumption that it is huge.

\param *in  The \c apop_params structure to be copied.
\param method_size  The size of the \c method_params structure. E.g. \c sizeof(apop_mle_params). Zero if there is no such structure.
\param model_size The size of the \c model_params structure. Zero if there is no such structure.
\param more_size The size of the \c more structure. Zero if there is no such structure.
\ingroup inv_and_est  */
apop_params *apop_params_clone(apop_params *in, size_t method_size, size_t model_size, size_t more_size){
    if (!in) return NULL;
  apop_params *out  = malloc(sizeof(apop_params));
    memcpy(out, in, sizeof(apop_params));
    if (method_size && in->method_params){
        out->method_params  = malloc(method_size);
        memcpy(out->method_params, in->method_params, method_size);
    }
    if (model_size && in->model_params){
        out->model_params  = malloc(model_size);
        memcpy(out->model_params, in->model_params, model_size);
    }
    if (more_size && in->more){
        out->more  = malloc(more_size);
        memcpy(out->more, in->more, more_size);
    }
    out->parameters = apop_data_copy(in->parameters);
    out->expected = apop_data_copy(in->expected);
    out->covariance = apop_data_copy(in->covariance);
    return out;
}


/** Free an \ref apop_params structure.

\ingroup inv_and_est */
void apop_params_free (apop_params * free_me){
	if (free_me->uses.predicted
	    || free_me->uses.expected)
        if (free_me->expected)
		    apop_data_free(free_me->expected);
	if (free_me->uses.covariance)
		apop_data_free(free_me->covariance);
	if (free_me->uses.parameters)
        apop_data_free(free_me->parameters);
	free(free_me);
}

/** Print the results of an estimation

\ingroup output */
void apop_params_show (apop_params * print_me){
    if (strlen(print_me->model->name))
        printf (print_me->model->name);
    if (strlen(print_me->model->name) && strlen(print_me->method_name))
        printf(" estimated via ");
    if (strlen(print_me->method_name))
        printf (print_me->method_name);
    printf("\n\n");
	if (print_me->uses.parameters)
        apop_data_show(print_me->parameters);
	if (print_me->uses.covariance){
		printf("\nThe covariance matrix:\n");
        apop_data_show(print_me->covariance);
	}
	if (print_me->uses.log_likelihood)
		printf("\nlog likelihood: \t%g\n", print_me->log_likelihood);
}

/** Currently an alias for \ref apop_params_show, but when I get
  around to it, it will conform better with the other apop_..._print
  fns.*/
void apop_params_print (apop_params * print_me){
    apop_params_show(print_me);
}

/** Outputs a copy of the \ref apop_model input.
\param in   The model to be copied
\return A pointer to a copy of the original, which you can mangle as you see fit. Includes a copy of the pointer to the original model's \c apop_params if any.
*/
apop_model * apop_model_copy(apop_model in){
  apop_model * out = malloc(sizeof(apop_model));
    memcpy(out, &in, sizeof(apop_model));
    return out;
}
