/** \file apop_model.c	 sets up the estimate structure which outputs from the various regressions and MLEs.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include <gsl/gsl_matrix.h>
#include "apophenia/output.h"
#include "apophenia/linear_algebra.h"

/** \defgroup inv_and_est  Using \ci apop_paramss 

The data to accompany an \c apop_model, including the input settings and the output parameters, expected values, et cetera.

See the page on \ref types for more.

\ingroup types
 */

/** Allocate an \ref apop_model.

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
    model->prepared         = 0;
    model->status           = 0;
    model->covariance	    = NULL;
    model->expected	        = NULL;
	return model;
}

/** Free an \ref apop_model structure.

   The  \c parameters, \c expected, and \c covariance elements are
   freed.  These are all the things that are completely copied, by
   \c apop_model_copy, so the parent model is still safe after this is
   called. \c data is not freed, because the odds are you still need it.

   The system has no idea what the \c method_settings, \c model_settings,
   and \c more elements contain, so if they point to other things,
   they need to be freed before calling this function.

  If \c free_me is \c NULL, this does nothing.

   \param free_me A pointer to the model to be freed.

\ingroup inv_and_est */
void apop_model_free (apop_model * free_me){
    if (!free_me) return;
    apop_data_free(free_me->parameters);
    apop_data_free(free_me->covariance);
    apop_data_free(free_me->expected);
	free(free_me);
}

/** Print the results of an estimation. If your model has a \c print
 method, then I'll use that, else I'll use a default setup.

 Your \c print method can use both by masking itself for a second:
 \code
void print_method(apop_model *in){
  void *temp = in->print;
  in->print = NULL;
  apop_model_show(in);
  in->print = temp;

  printf("Additional info:\n");
  ...
}
 \endcode

\ingroup output models */
void apop_model_show (apop_model * print_me){
    if (print_me->print){
        print_me->print(print_me);
        return;
    }
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

/** Currently an alias for \ref apop_model_show, but when I get
  around to it, it will conform better with the other apop_..._print
  fns.*/
void apop_model_print (apop_model * print_me){
    apop_model_show(print_me);
}

/** Outputs a copy of the \ref apop_model input.
\param in   The model to be copied
\return A pointer to a copy of the original, which you can mangle as you see fit. 

\ingroup models
*/
apop_model * apop_model_copy(apop_model in){
  apop_model * out = malloc(sizeof(apop_model));
    memcpy(out, &in, sizeof(apop_model));
    if (in.method_settings_size){
        out->method_settings  = malloc(in.method_settings_size);
        //out->method_settings_size  = in.method_settings_size;
        memcpy(out->method_settings, in.method_settings, in.method_settings_size);
    }
    if (in.model_settings_size){
        out->model_settings  = malloc(in.model_settings_size);
        //out->model_settings_size  = in.model_settings_size;
        memcpy(out->model_settings, in.model_settings, in.model_settings_size);
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

/** Take in an unparametrized \c apop_model and return a
  new \c apop_model with the given parameters. This would have been
  called apop_model_parametrize, but the OED lists four acceptable
  spellings for parameterise, so it's not a great candidate for a function name.

For example, if you need a N(0,1) quickly: 
\code
apop_model *std_normal = apop_model_set_parameters(apop_normal, 0.0, 1.0);
\endcode

Warning: Your parameters need to be <tt>double</tt>s, not <tt>int</tt>s. 
If you were to use <tt> apop_model_set_parameters(apop_normal, 0, 1);</tt>, you'd wind up with a N(0,0). This is an unfortunate feature of C's variadic function handling.

This doesn't take in data, so it won't work with models that take the number of parameters from the data, and it will only set the vector of the model's parameter apop_data set. This is most standard models, so that's not a real problem either.
If you have a situation where these options are out, you'll have to do something like
<tt>apop_model *new = apop_model_copy(in); apop_model_clear(your_data, in);</tt> and then set in->parameters using your data.

  \param in An unparametrized model, like \c apop_normal or \c apop_poisson.
  \param ... The list of parameters.

  \ingroup models
  */
apop_model *apop_model_set_parameters(apop_model in, ...){
  va_list  ap;
    if (in.vbase == -1 || in.m1base == -1 || in.m2base == -1)
        apop_error(0, 's', "%s only works with models whose number of params does not depend on data size. You'll have to use apop_model *new = apop_model_copy(in); apop_model_clear(your_data, in); and then set in->parameters using your data.\n");
    apop_model *out = apop_model_copy(in);
    apop_model_clear(NULL, out);
    va_start(ap, in);
    apop_vector_vfill(out->parameters->vector, ap);
    va_end(ap);
    return out; 
}

/* estimate the parameters of a model given data.

   This is a brief convenience function, which expands to \c m.estimate(d,&m). If your model has no \c estimate method, then I assume \c apop_maximum_likelihood(d, m), with the default MLE params.


\param d    The data
\param m    The model
\return     A pointer to an output model, which typically matches the input model but has its \c parameters element filled in.

\ingroup models
*/
apop_model *apop_estimate(apop_data *d, apop_model m){
    if (m.estimate)
        return m.estimate(d, &m); 
    return apop_maximum_likelihood(d, m);
}

/* Find the probability of a data/parametrized model pair.

\param d    The data
\param m    The parametrized model, which must have either a \c log_likelihood or a \c p method.

\ingroup models
*/
double apop_p(apop_data *d, apop_model *m){
    if (!m->parameters){
        apop_error(0, 's', "%s: You gave me a function that has no parameters.\n", __func__);
        return 0;
    }
    if (m->prep && !m->prepared){
        m->prep(d, m);
        m->prepared++;
    }
    if (m->p)
        return m->p(d, m);
    else if (m->log_likelihood)
        return exp(m->log_likelihood(d, m));
    apop_error(0, 's', "%s: You asked for the log likelihood of a model that has neither p nor log_likelihood methods.\n", __func__);
    return 0;
}

/* Find the log likelihood of a data/parametrized model pair.

\param d    The data
\param m    The parametrized model, which must have either a \c log_likelihood or a \c p method.

\ingroup models
*/
double apop_log_likelihood(apop_data *d, apop_model *m){
    if (!m->parameters){
        apop_error(0, 's', "%s: You gave me a function that has no parameters.\n", __func__);
        return 0;
    }
    if (m->prep && !m->prepared){
        m->prep(d, m);
        m->prepared++;
    }
    if (m->log_likelihood)
        return m->log_likelihood(d, m);
    else if (m->p)
        return log(m->p(d, m));
    apop_error(0, 's', "%s: You asked for the log likelihood of a model that has neither p nor log_likelihood methods.\n", __func__);
    return 0;
}

/** Find the vector of derivatives of the log likelihood of a data/parametrized model pair.

\param d    The data
\param m    The parametrized model, which must have either a \c log_likelihood or a \c p method.

\ingroup models
*/
void apop_score(apop_data *d, gsl_vector *out, apop_model *m){
    if (!m->parameters){
        apop_error(0, 's', "%s: You gave me a function that has no parameters.\n", __func__);
        return;
    }
    if (m->prep && !m->prepared){
        m->prep(d, m);
        m->prepared++;
    }
    if (m->score){
        m->score(d, out, m);
        return;
    }
    gsl_vector * numeric_default = apop_numerical_gradient(d, m);
    gsl_vector_memcpy(out, numeric_default);
    gsl_vector_free(numeric_default);
}


/** draw from a model. If the model has its own RNG, then you're good to
 go; if not, then do a simulated annealing run to generate a million draws from the data. 

 That second half is actually forthcoming...

 \ingroup models
*/
void apop_draw(double *out, gsl_rng *r, apop_model *m){
    if (m->draw){
        m->draw(out,r, m); 
        return;
    } 
}


/** Some models have a \c model_settings element that consists of a
 string. In that case, you can copy off a new model and set that setting 
 at the same time with this function. 
\param m The base model to be copied.
\param param The string to be placed in the <tt>model_settings</tt> slot.
\return A copy of \c m, with the appropriately set <tt>model_settings</tt> element.
 

\ingroup models
 */
apop_model *apop_model_copy_set_string(apop_model m, char* param){
  apop_model *out = apop_model_copy(m);
    out->model_settings = malloc(strlen(param)+1);
        strcpy((char *) out->model_settings ,param);
    out->model_settings_size = strlen(param)+1;
    return out;
}

/** The default prep is to simply call \c apop_model_clear. If the
 function has a prep method, then that gets called instead.

\ingroup models
 */
void apop_model_prep(apop_data *d, apop_model *m){
    if (m->prepared)
        return;
    if (m->prep)
        m->prep(d, m);
    else
        apop_model_clear(d, m);
    m->prepared++;
}
