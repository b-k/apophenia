/** \file apop_model.c	 sets up the estimate structure which outputs from the various regressions and MLEs.*/
/* Copyright (c) 2006--2011 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#define Declare_type_checking_fns
#include "apop_internal.h"

/** Allocate an \ref apop_model.

This sets up the output elements of the \c apop_model: the parameters and info. 

At close, the input model has parameters of the correct size.

\li This is the default action for \ref apop_prep. If your model has its own \c prep method, then that gets used instead, but most don't (or call \ref apop_model_clear at the end of their prep routine).

\li \ref apop_estimate calls \ref apop_prep internally. 

The above two points mean that you probably don't need to call this function directly.

\param data If your params vary with the size of the data set, then the function needs a data set to calibrate against. Otherwise, it's OK to set this to \c NULL.
\param model    The model whose output elements will be modified.
\return A pointer to the same model, should you need it.
\exception outmodel->error=='d' dimension error.
\ingroup models  */
apop_model * apop_model_clear(apop_data * data, apop_model *model){
    Get_vmsizes(data)
    int width = msize2 ? msize2 : -firstcol;//use the vector only if there's no matrix.
    Apop_stopif(model->dsize==-1 && !width, model->error='d', 0, "The model's dsize==-1, meaning size=data width, but the input data has NULL vector and matrix.");
    Apop_stopif(model->vsize==-1 && !width, model->error='d', 0, "The model's vsize==-1, meaning size=data width, but the input data has NULL vector and matrix.");
    Apop_stopif(model->msize1==-1 && !width, model->error='d', 0, "The model's msize1==-1, meaning size=data width, but the input data has NULL vector and matrix.");
    Apop_stopif(model->msize2==-1 && !width, model->error='d', 0, "The model's msize2==-1, meaning size=data width, but the input data has NULL vector and matrix.");

    model->dsize  = (model->dsize == -1 ? width : model->dsize);
    vsize  = model->vsize  == -1 ? width : model->vsize;
    msize1 = model->msize1 == -1 ? width : model->msize1 ;
    msize2 = model->msize2 == -1 ? width : model->msize2 ;
    if (!model->parameters) model->parameters = apop_data_alloc(vsize, msize1, msize2);
    if (!model->info) model->info = apop_data_alloc();
    free(model->info->names->title);
    Asprintf(&model->info->names->title, "<Info>");
    model->data = data;
	return model;
}

/** Free an \ref apop_model structure.

The  \c parameters element is freed.  These are all the things that are completely copied, by \c apop_model_copy, so the parent model is still safe after this is called. \c data is not freed, because the odds are you still need it.

If <tt>free_me->more_size</tt> is positive, the function runs
<tt>free(free_me->more)</tt>. But it has no idea what the \c more element contains;
if it points to other structures (like an \ref apop_data set), you need to free them
before calling this function.

\li If \c free_me is \c NULL, this does nothing.

\param free_me A pointer to the model to be freed.

\ingroup models */
void apop_model_free (apop_model * free_me){
    if (!free_me) return;
    apop_data_free(free_me->parameters);
    if (free_me->settings){
        int   i=0;
        while (free_me->settings[i].name[0]){
            if (free_me->settings[i].free)
                ((void (*)(void*))(free_me->settings[i].free))(free_me->settings[i].setting_group);
            i++;
        }
        free(free_me->settings);
    }
    if (free_me->more_size)
        free(free_me->more);
    if (free_me->info)
        apop_data_free(free_me->info);
	free(free_me);
}

/** Print the results of an estimation for a human to look over.

\param print_me The model whose information should be displayed
\param ap  The output stream. If \c NULL, use \c stdout. If you'd like something else, use \c fopen. E.g.:
\code
FILE *out =fopen("outfile.txt", "w"); //or "a" to append.
apop_model_print(the_model, out);
fclose(out);  //optional in most cases.
\endcode

\li The default prints the name, parameters, info, &c. but I check a vtable for
alternate methods you define; see \ref vtables for details. The typedef new functions
must conform to and the hash used for lookups are:

\code
typedef void (*apop_model_print_type)(apop_model *params, FILE *out);
#define apop_model_print_hash(m1) ((m1)->log_likelihood ? (size_t)(m1)->log_likelihood : \
            (m1)->p ? (size_t)(m1)->p*33 : \
            (m1)->estimate ? (size_t)(m1)->estimate*33*33 : \
            (m1)->draw ? (size_t)(m1)->draw*33*27  : \
            (m1)->cdf ? (size_t)(m1)->cdf*27*27 : 27)
\endcode

\li All output should \c fprintf to the input \c FILE* handle. 
  Apophenia's output routines also accept a file handle; e.g., if the file handle is
  named \c out, then if the \c thismodel print method uses \c apop_data_print to
  print the parameters, it must do so via a form like <tt>apop_data_print(thismodel->parameters,
  .output_pipe=ap)</tt>.

\li Your \c print method can use both by masking itself for a second:
 \code
void print_method(apop_model *in, FILE* ap){
  void *temp = in->estimate;
  in->estimate = NULL;
  apop_model_print(in, ap);
  in->estimate = temp;

  printf("Additional info:\n");
  ...
}
 \endcode

\li Print methods are intended for human consumption and are subject to change.

\ingroup output */
void apop_model_print (apop_model * print_me, FILE *ap){
    if (!ap) ap = stdout;
    apop_model_print_type mpf = apop_model_print_vtable_get(print_me);
    if (mpf){
        mpf(print_me, ap);
        return;
    }
    if (strlen(print_me->name)) fprintf (ap, "%s", print_me->name);
    fprintf(ap, "\n\n");
	if (print_me->parameters) apop_data_print(print_me->parameters, .output_pipe=ap);
    Get_vmsizes(print_me->info); //maxsize
    if (print_me->info && maxsize) apop_data_print(print_me->info, .output_pipe=ap);
}

/* Alias for \ref apop_model_print. Use that one. */
void apop_model_show (apop_model * print_me){
    apop_model_print(print_me, NULL);
}

/** Outputs a copy of the \ref apop_model input.
\param in The model to be copied
\return A pointer to a copy of the original, which you can mangle as you see fit. Includes copies of all settings groups, and the \c parameters (if not \c NULL, copied via \ref apop_data_copy).

\li If <tt>in.more_size > 0</tt> I <tt>memcpy</tt> the \c more pointer from the original data set.

\exception out->error=='a' Allocation error. In extreme cases, where there aren't even a few hundred bytes available, I will return \c NULL.
\exception out->error=='s' Error copying settings groups.
\exception out->error=='p' Error copying parameters or info page; the given \ref apop_data struct may be \c NULL or may have its own <tt>->error</tt> element.
\ingroup models
*/
apop_model * apop_model_copy(apop_model *in){
    Apop_stopif(!in, return NULL, 1, "Copying a NULL input; returning NULL.");
    apop_model * out = malloc(sizeof(apop_model));
    Apop_stopif(!out, return NULL, 0, "Serious allocation error; returning NULL.");
    memcpy(out, in, sizeof(apop_model));
    if (in->more_size){
        out->more  = malloc(in->more_size);
        Apop_stopif(!out->more, out->error='a'; return out, 0, "Allocation error setting up the ->more pointer.");
        memcpy(out->more, in->more, in->more_size);
    }
    int i=0; 
    out->settings = NULL;
    if (in->settings)
        do 
            apop_settings_copy_group(out, in, in->settings[i].name);
        while (strlen(in->settings[i++].name));
    out->parameters = apop_data_copy(in->parameters);
    Apop_stopif(in->parameters && (!out->parameters || out->parameters->error), 
                    out->error='p'; return out, 0, "Error copying the model parameters.");
    out->info = apop_data_copy(in->info);
    Apop_stopif(in->info && (!out->info || out->info->error), 
                    out->error='p'; return out, 0, "Error copying the info segment.");
    return out;
}

/** \def apop_model_set_parameters(in, ...)
Take in an unparameterized \c apop_model and return a new \c apop_model with the given parameters.  
For example, if you need a N(0,1) quickly:
\code
apop_model *std_normal = apop_model_set_parameters(apop_normal, 0, 1);
\endcode

This doesn't take in data, so it won't work with models that take the number of parameters from the data, and it will only set the vector of the model's parameter \ref apop_data set. This is most standard models, so that's not a real problem either.
If you have a situation where these options are out, you'll have to do something like
<tt>apop_model *new = apop_model_copy(in); apop_model_clear(your_data, in);</tt> and then set \c in->parameters using your data.

\param in An unparameterized model, like \ref apop_normal or \ref apop_poisson.
\param ... The list of parameters.
\return A copy of the input model, with parameters set.
\exception out->error=='d' dimension error: you gave me a model with an indeterminate number of parameters. Set .vsize or .msize1 and .msize2 first, then call this fn, or use apop_model *new = apop_model_copy(in); apop_model_clear(your_data, in); and then call this (because apop_model_clear sets the dimension based on your data size).
\hideinitializer   
\ingroup models
\li This would have been called apop_model_parametrize, but the OED lists four acceptable spellings for parameterise, so it's not a great candidate for a function name.
*/
apop_model *apop_model_set_parameters_base(apop_model *in, double ap[]){
    apop_model *out = apop_model_copy(in);
    apop_prep(NULL, out);
    Apop_stopif((in->vsize == -1) || (in->msize1 == -1) || (in->msize2 == -1), out->error='d', 
            0, "This function only works with models whose number of params does not "
            "depend on data size. You'll have to use apop_model *new = apop_model_copy(in); "
           " apop_model_clear(your_data, in); and then set in->parameters using your data.");
    apop_data_fill_base(out->parameters, ap);
    return out; 
}

/** estimate the parameters of a model given data.

This function copies the input model, preps it, and calls \c
m.estimate(d, m). If your model has no \c estimate method, then I
assume \c apop_maximum_likelihood(d, m), with the default MLE params.

I assume that you are using this function rather than directly calling the
model's the \c estimate method. For example, the \c estimate
method may assume that \c apop_prep has already been called.

\param d    The data
\param m    The model
\return     A pointer to an output model, which typically matches the input model but has its \c parameters element filled in.

\ingroup models
*/
apop_model *apop_estimate(apop_data *d, apop_model *m){
    apop_model *out = apop_model_copy(m);
    apop_prep(d, out);
    if (out->estimate) out->estimate(d, out); 
    else               apop_maximum_likelihood(d, out);
    return out;
}

/** Find the probability of a data/parametrized model pair.

\param d The data
\param m The parametrized model, which must have either a \c log_likelihood or a \c p method.

\ingroup models
*/
double apop_p(apop_data *d, apop_model *m){
    Nullcheck_m(m, GSL_NAN);
    if (m->p)
        return m->p(d, m);
    else if (m->log_likelihood)
        return exp(m->log_likelihood(d, m));
    Apop_stopif(0, , 0, "You asked for the probability of a model that has neither p nor log_likelihood methods.");
    return GSL_NAN;
}

/** Find the log likelihood of a data/parametrized model pair.

\param d    The data
\param m    The parametrized model, which must have either a \c log_likelihood or a \c p method.

\ingroup models
*/
double apop_log_likelihood(apop_data *d, apop_model *m){
    Nullcheck_m(m, GSL_NAN); //Nullcheck_p(m); //Too many models don't use the params.
    if (m->log_likelihood)
        return m->log_likelihood(d, m);
    else if (m->p)
        return log(m->p(d, m));
    Apop_stopif(0, , 0, "You asked for the log likelihood of a model that has neither p nor log_likelihood methods.");
    return GSL_NAN;
}

/** Find the vector of first derivatives (aka the gradient) of the log likelihood of a data/parametrized model pair.

\param d    The data
\param out  The score to be returned. I expect you to have allocated this already.
\param m    The parametrized model, which must have either a \c log_likelihood or a \c p method.

\li The default is to use \ref apop_numerical_gradient, but special-case calculations
for certain models are held in a vtable; see \ref vtables for details. The typedef
new functions must conform to and the hash used for lookups are:

\code
typedef void (*apop_score_type)(apop_data *d, gsl_vector *gradient, apop_model *m);
#define apop_score_hash(m1) ((size_t)((m1).log_likelihood ? (m1).log_likelihood : (m1).p))
\endcode

As input to your function, you can expect that the model \c m is sufficiently prepped
that the log likelihood can be evaluated; see \ref psubsection for details.
On output, the a \c gsl_vector input to the function must be filled with the gradients
(or <tt>NaN</tt>s on errors). If the model parameters have a more complex shape
than a simple vector, then the vector must be in \c apop_data_pack order; use \c
apop_data_unpack to reformat to the preferred shape.

\ingroup models
*/
void apop_score(apop_data *d, gsl_vector *out, apop_model *m){
    Nullcheck_m(m, );
    apop_score_type ms = apop_score_vtable_get(m);
    if (ms){
        ms(d, out, m);
        return;
    }
    gsl_vector * numeric_default = apop_numerical_gradient(d, m);
    gsl_vector_memcpy(out, numeric_default);
    gsl_vector_free(numeric_default);
}

Apop_settings_init(apop_pm,
    //defaults include base=NULL, index=0, own_rng=0
    Apop_varad_set(rng, NULL);
    Apop_varad_set(draws, 1e4);
)

Apop_settings_copy(apop_pm,)

Apop_settings_free(apop_pm, )

void distract_doxygen(){/*Doxygen gets thrown by the settings macros. This decoy function is a workaround. */}

/** Get a model describing the distribution of the given parameter estimates.

For many models, the parameter estimates are well-known, such as the
\f$t\f$-distribution of the parameters for OLS.

For models where the distribution of \f$\hat{}p\f$ is not known, if you give me data, I
will return a \ref apop_normal or \ref apop_multivariate_normal model, using the parameter estimates as mean and \ref apop_bootstrap_cov for the variances.

If you don't give me data, then I will assume that this is a stochastic model where 
re-running the model will produce different parameter estimates each time. In this case, I will
run the model 1e4 times and return a \ref apop_pmf model with the resulting parameter
distributions.

Before calling this, I expect that you have already run \ref apop_estimate to produce \f$\hat{}p\f$.

The \ref apop_pm_settings structure dictates details of how the model is generated.
For example, if you want only the distribution of the third parameter, and you know the
distribution will be a PMF generated via random draws, then set settings and call the
model via:
\code
  apop_model_group_add(your_model, apop_pm, .index =3, .draws=3e5);
  apop_model *dist = apop_parameter_model(your_data, your_model);
\endcode

\li \c index gives the position of the parameter (in \ref apop_data_pack order)
in which you are interested. Thus, if this is zero or more, then you will get a
univariate output distribution describing a single parameter. If <tt>index == -1</tt>,
then I will give you the multivariate distribution across all parameters.  The default
is zero (i.e. the univariate distribution of the zeroth parameter).
\li \c rng If the method requires random draws (as the default bootstrap will), then use this. If you provide \c NULL and one is needed, I provide one for you via <tt>apop_rng_alloc(apop_opts.rng_seed++)</tt>.
\li \c draws If there is no closed-form solution and bootstrap is inappropriate, then
the last resort is a large numbr of random draws of the model, summarized into a PMF. Default: 1,000 draws.
\li The default is via resampling as above, but special-case calculations for certain models are held in a vtable; see \ref vtables for details. The typedef new functions must conform to and the hash used for lookups are:

\code
typedef apop_model* (*apop_parameter_model_type)(apop_data *, apop_model *);
#define apop_parameter_model_hash(m1) ((size_t)((m1).log_likelihood ? (m1).log_likelihood : (m1).p)*33 + (m1).estimate ? (size_t)(m1).estimate: 27)
\endcode
\ingroup models
*/ 
apop_model *apop_parameter_model(apop_data *d, apop_model *m){
    apop_pm_settings *settings = apop_settings_get_group(m, apop_pm);
    if (!settings)
        settings = Apop_settings_add_group(m, apop_pm, .base= m);
    apop_parameter_model_type pm = apop_parameter_model_vtable_get(m);
    if (pm) return pm(d, m);
    else if (d){
        Get_vmsizes(m->parameters);//vsize, msize1, msize2
        apop_model *out = apop_model_copy(apop_multivariate_normal);
        out->msize1 = out->vsize = out->msize2 = out->dsize = vsize+msize1+msize2;
        out->parameters = apop_bootstrap_cov(d, m, settings->rng, settings->draws);
        out->parameters->vector = apop_data_pack(m->parameters);
        if (settings->index == -1)
            return out;
        else {
            apop_model *out2 = apop_model_set_parameters(apop_normal, 
                    apop_data_get(out->parameters, settings->index, -1), //mean
                    apop_data_get(out->parameters, settings->index, settings->index)//var
                    );
            apop_model_free(out);
            return out2;
        }
    } //else
    Get_vmsizes(m->parameters);//vsize, msize1, msize2
    apop_data *param_draws = apop_data_alloc(0, settings->draws, vsize+msize1+msize2);
    for (int i=0; i < settings->draws; i++){
        apop_model *mm = apop_estimate (NULL, m);//If you're here, d==NULL.
        Apop_row_v(param_draws, i, onerow);
        apop_data_pack(mm->parameters, onerow);
        apop_model_free(mm);
    }
    if (settings->index == -1)
        return apop_estimate(param_draws, apop_pmf);
    else {
        apop_data *param_draws1 = apop_data_alloc(settings->draws, 0,0);
        Apop_col_v(param_draws, settings->index, the_draws);
        gsl_vector_memcpy(param_draws1->vector, the_draws);
        apop_data_free(param_draws);
        return apop_estimate(param_draws1, apop_pmf);
    }
}

extern apop_model *apop_swap_model; //apop_missing_data.c
int apop_model_metropolis_draw(double *out, gsl_rng* rng, apop_model *params);//apop_update.c

/** Draw from a model. 

\param out An already-allocated array of <tt>double</tt>s to be filled by the draw method. This probably has size <tt>your_model->dsize</tt>.
\param r   A \c gsl_rng, probably allocated via \ref apop_rng_alloc. Optional; if \c NULL, then I will call \ref apop_rng_get_thread for an RNG.
\param m   The model from which to make draws.

\li If the model has its own \c draw method, then use that.
\li If the model is univariate, use \ref apop_arms_draw to generate random draws.
\li If the model is multivariate, use \ref apop_model_metropolis to generate random draws.

\li This makes a single draw of the given size. See \ref apop_model_draws to fill a matrix with draws.

\return Zero on success; nozero on failure. <tt>out[0]</tt> is probably \c NAN on failure.
\ingroup models
*/
int apop_draw(double *out, gsl_rng *r, apop_model *m){
    if (!r) r = apop_rng_get_thread();
    if (m->draw)
        return m->draw(out, r, m); 
    else if (m->dsize == 1)
        return apop_arms_draw(out, r, m);
    //Else, MCMC, possibly setting it up first.
    //generate a model with data/params reversed
    //estimate mcmc. Swapped model will be stored as settings->base_model.
    OMP_critical (apop_draw)
    if (!Apop_settings_get_group(m, apop_mcmc)){
        apop_model *swapped = apop_model_copy(apop_swap_model);
        swapped->more = m;
        swapped->msize1 = 1;
        swapped->msize2 = m->dsize;
        swapped->data = m->parameters;
        Apop_settings_add_group(swapped, apop_mcmc, .burnin=0.999, .periods=1000);
        apop_model *est = apop_model_metropolis(m->parameters, swapped, r); //leak.
        m->draw = apop_model_metropolis_draw;
        apop_settings_copy_group(m, est, "apop_mcmc");
    }
    return apop_draw(out, r, m);
}

/** The default prep is to simply call \ref apop_model_clear. If the
 function has a prep method, then that gets called instead.

\ingroup models
 */
void apop_prep(apop_data *d, apop_model *m){
    if (m->prep) m->prep(d, m);
    else         apop_model_clear(d, m);
}

static double disnan(double in) {return gsl_isnan(in);}

/** A prediction supplies E(a missing value | original data, already-estimated parameters, and other supplied data elements ).

For a regression, one would first estimate the parameters of the model, then supply a row of predictors <b>X</b>. The value of the dependent variable \f$y\f$ is unknown, so the system would predict that value. [In some models, this may not be the expected value, but is a best value for the missing item using some other meaning of `best'.]

For a univariate model (i.e. a model in one-dimensional data space), there is only one variable to omit and fill in, so the prediction problem reduces to the expected value: E(a missing value | original data, already-estimated parameters).

In other cases, prediction is the missing data problem: for three-dimensional data, you may supply the input (34, \c NaN, 12), and the parameterized model provides the most likely value of the middle
parameter.

\li If you give me a \c NULL data set, I will assume you want all values filled in---the expected value.

\li If you give me data with \c NaNs, I will take those as the points to
be predicted given the provided data.

If the model has no \c predict method, the default is to use the \ref apop_ml_impute function to do the work.

\return If you gave me a non-\c NULL data set, I will return that, with the zeroth column or the \c NaNs filled in.  If \c NULL input, I will allocate an \ref apop_data set and fill it with the expected values.

There may be a second page (i.e., a \ref apop_data set attached to the <tt>->more</tt> pointer of the main) listing confidence and standard error information. See your specific model documentation for details.

This segment of the framework is in beta---subject to revision of the details.

\li The default is to use \ref apop_ml_impute, but special-case calculations for certain models are held in a vtable; see \ref vtables for details. The typedef new functions must conform to and the hash used for lookups are:

\code
typedef apop_data * (*apop_predict_type)(apop_data *d, apop_model *params);
#define apop_predict_hash(m1) ((size_t)((m1).log_likelihood ? (m1).log_likelihood : (m1).p)*33 + (m1).estimate ? (size_t)(m1).estimate: 27)
\endcode
\ingroup models
*/
apop_data *apop_predict(apop_data *d, apop_model *m){
    apop_data *prediction = NULL;
    apop_data *out = d ? d : apop_data_alloc(0, 1, m->dsize);
    if (!d) gsl_matrix_set_all(out->matrix, GSL_NAN);
    apop_predict_type mp = apop_predict_vtable_get(m);
    if (mp) prediction = mp(out, m);
    if (prediction) return prediction;
    if (!apop_map_sum(out, disnan)) return out;
    //default:
    apop_model *f = apop_ml_imputation(out, m);
    apop_model_free(f);
    return out;
}

/* Are all the elements of v less than or equal to the corresponding elements of the reference vector? */
static int lte(gsl_vector *v, gsl_vector *ref){
    for (int i=0; i< v->size; i++) 
        if(v->data[i] > gsl_vector_get(ref, i))
            return 0;
    return 1;
}

/** Input a data point in canonical form and a model; returns the area of the model's PDF beneath the given point.

  By default, I just make random draws from the PDF and return the percentage of those
  draws beneath or equal to the given point. Many models have closed-form solutions that
  make no use of random draws. 

See also \ref apop_cdf_settings, which is the structure I use to store draws already made (which means the second, third, ... calls to this function will take much less time than the first), the \c gsl_rng, and the number of draws to be made. These are handled without your involvement, but if you would like to change the number of draws from the default, add this group before calling \ref apop_cdf :

\code
Apop_model_add_group(your_model, apop_cdf, .draws=1e5, .rng=my_rng);
double cdf_value = apop_cdf(your_data_point, your_model);
\endcode

Here are many examples using common, mostly symmetric distributions.

\include some_cdfs.c

\ingroup models
  */
double apop_cdf(apop_data *d, apop_model *m){
    if (m->cdf) return m->cdf(d, m);
    apop_cdf_settings *cs = Apop_settings_get_group(m, apop_cdf);
    if (!cs) cs = Apop_model_add_group(m, apop_cdf);
    long int tally = 0; 
    Apop_row(d, 0, row);
    gsl_vector *ref = apop_data_pack(row);
    if (!cs->draws_made){
        if (m->dsize == -1) apop_prep(d, m);
        Apop_stopif(m->dsize==0, return GSL_NAN, 0, "I need to make random draws from your model, but it has dsize==0. Returning NaN");
        cs->draws_made = gsl_matrix_alloc(cs->draws, m->dsize);
        for (int i=0; i< cs->draws; i++){
            Apop_matrix_row(cs->draws_made, i, onerow);
            apop_draw(onerow->data, cs->rng, m);
        }
    }
    for (int i=0; i< cs->draws_made->size1; i++){
        Apop_matrix_row(cs->draws_made, i, onerow);
        tally += lte(onerow, ref);
    }
    gsl_vector_free(ref);
    return tally/(double)cs->draws_made->size1;
}

Apop_settings_init(apop_cdf,
    Apop_varad_set(draws, 1e4);
    Apop_varad_set(rng, NULL);
    out->draws_refcount = malloc(sizeof(int));
    *out->draws_refcount = 1;
)

Apop_settings_free(apop_cdf,
    if (in->draws_made && !--*in->draws_refcount)
        gsl_matrix_free(in->draws_made);
    apop_model_free(in->cdf_model);
)

Apop_settings_copy(apop_cdf,
    ++*out->draws_refcount;
)
