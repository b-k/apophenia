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

This sets up the output elements of the \c apop_model: the parameters,
covarinace, and expected data sets. 

At close, the input model has parameters of the correct size, the
covariance and expected elements are \c NULL, and the \c status element
is zero, indicating no estimation has been done yet.

The input model is modified, so you probably want to call this after you call \c apop_model_copy.

\param data If your params vary with the size of the data set, then the function needs a data set to calibrate against. Otherwise, it's OK to set this to \c NULL
\param model    The model whose output elements will be modified.
\return A pointer to the same model, should you need it.

\ingroup inv_and_est  */
apop_model * apop_model_clear(apop_data * data, apop_model *model){
  int vsize  = model->vbase == -1 ? data->matrix->size2 : model->vbase;
  int msize1 = model->m1base == -1 ? data->matrix->size2 : model->m1base ;
  int msize2 = model->m2base == -1 ? data->matrix->size2 : model->m2base ;
    apop_data_free(model->parameters);
    apop_data_free(model->covariance);
    apop_data_free(model->expected);
    model->parameters	    = apop_data_alloc(vsize, msize1, msize2);
    model->data             = data;
    model->status           = 0;
    model->covariance	    = NULL;
    model->expected	        = NULL;
	return model;
}

/** Free an \ref apop_model structure.

  If \c free_me is \c NULL, this does nothing.
\ingroup inv_and_est */
void apop_params_free (apop_model * free_me){
    if (!free_me) return;
	if (free_me->parameters)
        apop_data_free(free_me->parameters);
	if (free_me->covariance)
		apop_data_free(free_me->covariance);
    if (free_me->expected)
        apop_data_free(free_me->expected);
	free(free_me);
}

/** Print the results of an estimation

\ingroup output */
void apop_params_show (apop_model * print_me){
    if (strlen(print_me->name))
        printf (print_me->name);
    if (strlen(print_me->name) && strlen(print_me->method_name))
        printf(" estimated via ");
    if (strlen(print_me->method_name))
        printf (print_me->method_name);
    printf("\n\n");
	if (print_me->parameters)
        apop_data_show(print_me->parameters);
	if (print_me->covariance){
		printf("\nThe covariance matrix:\n");
        apop_data_show(print_me->covariance);
	}
//under the false presumption that if it is calculated it is never quite ==0.
	if (print_me->log_likelihood)
		printf("\nlog likelihood: \t%g\n", print_me->llikelihood);
}

/** Currently an alias for \ref apop_params_show, but when I get
  around to it, it will conform better with the other apop_..._print
  fns.*/
void apop_params_print (apop_model * print_me){
    apop_params_show(print_me);
}

/** Outputs a copy of the \ref apop_model input.
\param in   The model to be copied
\return A pointer to a copy of the original, which you can mangle as you see fit. Includes a copy of the pointer to the original model's \c apop_model if any.
*/
apop_model * apop_model_copy(apop_model in){
  apop_model * out = malloc(sizeof(apop_model));
    memcpy(out, &in, sizeof(apop_model));
    if (in.method_params_size){
        out->method_params  = malloc(in.method_params_size);
        //out->method_params_size  = in.method_params_size;
        memcpy(out->method_params, in.method_params, in.method_params_size);
    }
    if (in.model_params_size){
        out->model_params  = malloc(in.model_params_size);
        //out->model_params_size  = in.model_params_size;
        memcpy(out->model_params, in.model_params, in.model_params_size);
    }
    if (in.more_size){
        out->more  = malloc(in.more_size);
        //out->more_size  = in.more_size;
        memcpy(out->more, in.more, in.more_size);
    }
    out->parameters = apop_data_copy(in.parameters);
    out->expected   = apop_data_copy(in.expected);
    out->covariance = apop_data_copy(in.covariance);
    return out;
}

/* Frees the memory taken by an \c apop_model. 
   Notably the \c method_params, \c model_params, \c more, \c parameters,
   \c expected, and \c covariance elements are freed.  These are all
   the things that are completely copied, by \c apop_model_copy, so the
   parent model is still safe after this is called.

   The system has no idea what the \c method_params, \c model_params,
   and \c more elements contain, so if they point to other things,
   they need to be freed before calling this function.

   \param in A pointer to the model to be freed.

 */
void apop_model_free(apop_model *in){
    free (in->method_params);
    free (in->model_params);
    free (in->more);
    apop_data_free(in->parameters);
    apop_data_free(in->expected);
    apop_data_free(in->covariance);
    free(in);
}

apop_model *apop_estimate(apop_data *d, apop_model m){
    return m.estimate(d, &m); }
