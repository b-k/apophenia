#include "asst.h"
#include "likelihoods.h"
#include "model.h"
/** \file apop_model_fix_params.c 
 Set some of the parameters of a model to fixed values.

 There's only one function here. Its header is in likelihoods.h
 
Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */


typedef struct {
    apop_data *paramvals, *mask, *gradient_for_base;
    apop_model *base_model;
} apop_model_fixed_params_settings;

apop_model_fixed_params_settings *apop_model_fixed_params_settings_alloc (apop_model base_model, 
            apop_data *paramvals, apop_data *mask, int size){
    apop_model_fixed_params_settings *out = malloc(sizeof(apop_model_fixed_params_settings));
    out->base_model  = apop_model_copy(base_model);
    out->paramvals  = paramvals;
    out->mask       = mask;
    out->base_model->parameters   = NULL;
    out->gradient_for_base       = apop_data_alloc(size,0,0); //use full size.
    return out;
}

void *apop_model_fixed_params_settings_copy (apop_model_fixed_params_settings *in ){ 
    apop_model_fixed_params_settings *out = malloc(sizeof(apop_model_fixed_params_settings));
    out->base_model  = in->base_model;
    out->paramvals  = in->paramvals;
    out->mask       = in->mask;
    out->gradient_for_base       = in->gradient_for_base;
    return out;
 } 

void apop_model_fixed_params_settings_free (apop_model_fixed_params_settings *in ){ free(in); }

static apop_model fixed_param_model;

static void  fixed_params_pack(const apop_data *in, apop_data  *out, apop_data  *mask, int is_gradient){
  int ctr   = 0;
    if (mask->vector)
        for(size_t i=0; i< mask->vector->size; i++)
            if (!apop_data_get(mask, i, -1))
                apop_data_set(out, ctr++, -1, apop_data_get(in, i, -1));
    if (mask->matrix)
        for(size_t i=0; i< mask->matrix->size1; i++)
            for(size_t j=0; j< mask->matrix->size2; j++)
                if (!apop_data_get(mask, i, j))
                    apop_data_set(out, ctr++, -1, apop_data_get(in, i, is_gradient ? -1 : j));
}

static void  fixed_param_unpack(const gsl_vector *in, apop_model_fixed_params_settings *p, int is_gradient){
  int ctr   = 0;
    if (!is_gradient && !p->base_model->parameters)
        p->base_model->parameters  = apop_data_copy(p->paramvals);
    apop_data *target = (is_gradient ? p->gradient_for_base: p->base_model->parameters);
    if (p->mask->vector)
        for(size_t i=0; i< p->mask->vector->size; i++)
            if (!apop_data_get(p->mask, i, -1))
                apop_data_set(target, i, -1, gsl_vector_get(in, ctr++));
    if (p->mask->matrix)
        for(size_t i=0; i< p->mask->matrix->size1; i++)
            for(size_t j=0; j< p->mask->matrix->size2; j++)
                if (!apop_data_get(p->mask, i, j))
                    apop_data_set(target, i, is_gradient ? -1 : j, gsl_vector_get(in, ctr++));
}

static double i_ll(apop_data *d, apop_model *fixed_model){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(fixed_model, "apop_model_fixed_params");
    fixed_param_unpack(fixed_model->parameters->vector, p, 0);
    return apop_log_likelihood(d, p->base_model);
}

static double i_p(apop_data *d, apop_model *params){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(params, "apop_model_fixed_params");
    fixed_param_unpack(params->parameters->vector, p, 0);
    return apop_p(d, p->base_model);
}

static void i_score(apop_data *d, gsl_vector *gradient, apop_model *params){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(params, "apop_model_fixed_params");
  static apop_data *gdummy   = NULL;
    fixed_param_unpack(params->parameters->vector, p, 0);
    fixed_param_unpack(gradient, p, 1);
    p->base_model->score(d, p->gradient_for_base->vector, p->base_model);

    if (!gdummy) gdummy =apop_data_alloc(0,0,0);
    gdummy->vector = gradient;
    gsl_vector_set_all(gradient, 0);
    fixed_params_pack(p->gradient_for_base, gdummy, p->mask, 1);
}

static double  i_constraint(apop_data *data, apop_model *params){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(params, "apop_model_fixed_params");
    fixed_param_unpack(params->parameters->vector, p, 0);
  double out = p->base_model->constraint(data, p->base_model);
    if (out) 
        fixed_params_pack(p->base_model->parameters, params->parameters, p->mask, 0);
    return out;
}

static void i_draw(double *out, gsl_rng* r, apop_model *eps){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(eps, "apop_model_fixed_params");
  apop_data             *tmp    = p->base_model->parameters;
    fixed_param_unpack(eps->parameters->vector, p, 0);
    p->base_model->draw(out, r, p->base_model);
    p->base_model->parameters   = tmp;
}

static apop_model *fixed_est(apop_data * data, apop_model *params){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(params, "apop_model_fixed_params");
    if (!data)
        data    = params->data;
    apop_model *e = apop_maximum_likelihood(data,  *params);
    fixed_param_unpack(e->parameters->vector, p, 0);
    apop_data_free(e->parameters);
    e->parameters   = apop_data_copy(p->base_model->parameters);
    return e;
}

static apop_model fixed_param_model = {"Fill me", .estimate=fixed_est, .p = i_p, .log_likelihood=i_ll, 
                                    .score=i_score, .constraint= i_constraint, .draw=i_draw};


/** Produce a model based on another model, but with some of the 
  parameters fixed at a given value. You input one \c apop_data set with the values at which you want 
  the values fixed, and one apop_data set which is zero at the location of the free parameters and 
  one at the location of the fixed parameters.

  You also need to input the base model (which should have any settings groups attached before calling this function), and two sets of parameters:
  the parameters for your model, which will be blindly passed on, and MLE-style parameters for the 
  \c estimate method. The \c estimate method always uses an MLE, and it never calls the base model's \c estimate method.

  The output is an \c apop_model that can be estimated, Bayesian updated, et cetera.

  Here is a sample program. It produces a few thousand draws from a multivariate normal distribution,
  and then tries to recover the means given a var/covar matrix fixed at the correct variance.

  \code


#include <apop.h>

int main(){
  apop_data *pv   = apop_data_alloc(2,2,2);
  apop_data *mask = apop_data_calloc(2,2,2);
  size_t    i, ct = 1000;
  apop_data *d  = apop_data_alloc(0,ct,2);
  gsl_rng *r    = apop_rng_alloc(10);
  double    draw[2];
  apop_data_fill(pv, 8.,  1., 0.5,
                     2., 0.5, 1.);
    apop_model *pvm = apop_model_copy(apop_multivariate_normal);
    pvm->parameters = pv;
    for(i=0; i< ct; i++){
        apop_draw(draw, r, pvm);
        apop_data_set(d, i, 0, draw[0]);
        apop_data_set(d, i, 1, draw[1]);
    }

    gsl_matrix_set_all(mask->matrix, 1);
    apop_model *mep1   = apop_model_fix_params(d, pv, mask, apop_multivariate_normal);
    apop_model *e1  = apop_estimate(d, *mep1);
    printf("original params: ");
    apop_vector_show(pv->vector);
    printf("estimated params: ");
    apop_vector_show(e1->parameters->vector);
}
  \endcode
  
  \param data       The model data
  \param paramvals  An \c apop_data set with the values of the variables to be fixed.
  \param mask       Set to zero for free values, one for values to be read from \c paramvals
  \param model_in   The base model
 \return an \c apop_mle_settings structure for you to fill as desired. If this is named \c m, then \c m->ep->estimate(NULL, m->ep) will run the estimation.
  */

apop_model *apop_model_fix_params(apop_data *data, apop_data *paramvals, apop_data *mask, apop_model model_in){
    int size  = (mask->vector? mask->vector->size : 0 )+
            (mask->matrix? mask->matrix->size1 * mask->matrix->size2 : 0);
    
  apop_model                *model_out  = apop_model_copy(fixed_param_model);
  Apop_settings_add_group(model_out, apop_model_fixed_params, model_in, paramvals, mask, size);

    if (mask->vector)
        for(size_t i=0; i< mask->vector->size; i++)
            if (apop_data_get(mask, i, -1))
                size    --;
    if (mask->matrix)
        for(size_t i=0; i< mask->matrix->size1; i++)
            for(size_t j=0; j< mask->matrix->size2; j++)
                if (apop_data_get(mask, i, j))
                    size    --;
    if (!size)
        apop_error(0,'s',"You're asking me to estimate a model where every single parameter is fixed.\n");
    model_out->vbase            = size;
  
    apop_data *starting_pt = apop_data_alloc(size, 0, 0);
    fixed_params_pack(paramvals, starting_pt, mask, 0);
    double *sp           = malloc(sizeof(double) * size);
    memcpy(sp, starting_pt->vector->data, sizeof(double) * size);
    Apop_model_add_group(model_out, apop_mle, .starting_pt=sp,.parent= model_out);
    apop_data_free(starting_pt);

/*    mle_out->starting_pt        = malloc(sizeof(double) * size);
    memcpy(mle_out->starting_pt, startingpt->vector->data, sizeof(double) * size);
    apop_data_free(startingpt);*/

    if (!model_in.p) model_out->p = NULL;
    if (!model_in.log_likelihood) model_out->log_likelihood = NULL;
    if (!model_in.score) model_out->score = NULL;
    if (!model_in.constraint) model_out->constraint = NULL;
    if (!model_in.draw) model_out->draw = NULL;
    snprintf(model_out->name, 100, "%s, with some params fixed", model_in.name);
    return model_out;
}
