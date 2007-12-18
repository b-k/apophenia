#include <apop.h>
/** \file apop_model_fix_params.c 
 Set some of the parameters of a model to fixed values.

 There's only one function here. Its header is in asst.h
 
Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */


typedef struct apop_model_fixed_params{
    apop_data *paramvals, *mask, *gradient_for_base;
    apop_model * base_model;
} apop_model_fixed_params;

static apop_model fixed_param_model;

static void  fixed_params_pack(const apop_data *in, apop_data  *out, apop_data  *mask, int is_gradient){
  int ctr   = 0;
  int i, j;
    if (mask->vector)
        for(i=0; i< mask->vector->size; i++)
            if (!apop_data_get(mask, i, -1))
                apop_data_set(out, ctr++, -1, apop_data_get(in, i, -1));
    if (mask->matrix)
        for(i=0; i< mask->matrix->size1; i++)
            for(j=0; j< mask->matrix->size2; j++)
                if (!apop_data_get(mask, i, j))
                    apop_data_set(out, ctr++, -1, apop_data_get(in, i, is_gradient ? -1 : j));
}

static void  fixed_param_unpack(const gsl_vector *in, apop_model_fixed_params *p, int is_gradient){
  int ctr   = 0;
  int i, j;
    if (!is_gradient && !p->base_model->parameters)
        p->base_model->parameters  = apop_data_copy(p->paramvals);
    apop_data *target = (is_gradient ? p->gradient_for_base: p->base_model->parameters);
    if (p->mask->vector)
        for(i=0; i< p->mask->vector->size; i++)
            if (!apop_data_get(p->mask, i, -1))
                apop_data_set(target, i, -1, gsl_vector_get(in, ctr++));
    if (p->mask->matrix)
        for(i=0; i< p->mask->matrix->size1; i++)
            for(j=0; j< p->mask->matrix->size2; j++)
                if (!apop_data_get(p->mask, i, j))
                    apop_data_set(target, i, is_gradient ? -1 : j, gsl_vector_get(in, ctr++));
}

static double i_ll(apop_data *d, apop_model *fixed_model){
  apop_model_fixed_params *p    = fixed_model->model_settings;
    fixed_param_unpack(fixed_model->parameters->vector, p, 0);
    return p->base_model->log_likelihood(d, p->base_model);
}

static double i_p(apop_data *d, apop_model *params){
  apop_model_fixed_params *p    = params->model_settings;
    fixed_param_unpack(params->parameters->vector, p, 0);
    return p->base_model->p(d, p->base_model);
}

static void i_score(apop_data *d, gsl_vector *gradient, apop_model *params){
  apop_model_fixed_params *p    = params->model_settings;
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
  apop_model_fixed_params *p    = params->model_settings;
    fixed_param_unpack(params->parameters->vector, p, 0);
  double out = p->base_model->constraint(data, p->base_model);
    if (out) 
        fixed_params_pack(p->base_model->parameters, params->parameters, p->mask, 0);
    return out;
}

static void i_draw(double *out, gsl_rng* r, apop_model *eps){
  apop_model_fixed_params *p    = eps->model_settings;
  apop_data             *tmp    = p->base_model->parameters;
    fixed_param_unpack(eps->parameters->vector, p, 0);
    p->base_model->draw(out, r, p->base_model);
    p->base_model->parameters   = tmp;
}

static apop_model *fixed_est(apop_data * data, apop_model *params){
    if (!params)
        apop_error(0,'c',"Fixed parameter model, closure error: you need to send the parameters as the second argument of the estimate function.\n");
  apop_model_fixed_params *p    = params->model_settings;
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

  You also need to input the base model, and two sets of parameters:
  the parameters for your model, which will be blindly passed on, and MLE-style parameters for the 
  \c estimate method. The \c estimate method always uses an MLE, and it never calls the base model's \c estimate method. So if that method does prep work before doing anything, that prep work won't get done.

  The output is an \c apop_mle_settings struct.

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
    apop_data_set(pv, 0, -1, 8);
    apop_data_set(pv, 1, -1, 2);
    gsl_matrix_set(pv->matrix, 0,0,1);
    gsl_matrix_set(pv->matrix, 1,1,1);
    gsl_matrix_set(pv->matrix, 1,0,0.5);
    gsl_matrix_set(pv->matrix, 0,1,0.5);
    apop_model *pvm = apop_model_copy(apop_multivariate_normal);
    pvm->parameters = pv;
    for(i=0; i< ct; i++){
        apop_draw(draw, r, pvm);
        apop_data_set(d, i, 0, draw[0]);
        apop_data_set(d, i, 1, draw[1]);
    }

    gsl_matrix_set_all(mask->matrix, 1);
    apop_mle_settings *mep1   = apop_model_fix_params(d, pv, mask, apop_multivariate_normal);
    apop_model   *e1  = apop_estimate(d, *mep1->model);
    gsl_vector_sub(e1->parameters->vector, pv->vector);
    assert(apop_vector_sum(e1->parameters->vector) < 1e-2);
}
  \endcode
  
  \param data       The model data
  \param paramvals  An \c apop_data set with the values of the variables to be fixed.
  \param mask       Set to zero for free values, one for values to be read from \c paramvals
  \param model_in   The base model
 \return an \c apop_mle_settings structure for you to fill as desired. If this is named \c m, then \c m->ep->estimate(NULL, m->ep) will run the estimation.
  */

apop_mle_settings *apop_model_fix_params(apop_data *data, apop_data *paramvals, apop_data *mask, apop_model model_in){
/*
--copy the model. Change the params to an appropriately-sized vector.
--create MLE_params mle_out. mle_out->ep are the base apop_model (B).
--B->model_settings is an apop_model_fixed_params. So set b->model_settings = p; set model->ep = B.
--p->base_model_settings are as input by the user.

So given mle_out:
base params = mle_out->ep
fixed model params = mle_out->ep->model_settings
original model = mle_out->ep->model_settings->m
original model params = mle_out->ep->model_settings->base_model_settings
*/
  apop_model_fixed_params   *p          = malloc(sizeof(*p));
  apop_model                *model_out  = apop_model_copy(fixed_param_model);
  int        i, j;
    if (!model_in.p) model_out->p = NULL;
    if (!model_in.log_likelihood) model_out->log_likelihood = NULL;
    if (!model_in.score) model_out->score = NULL;
    if (!model_in.constraint) model_out->constraint = NULL;
    if (!model_in.draw) model_out->draw = NULL;
    snprintf(model_out->name, 100, "%s, with some params fixed", model_in.name);
    int size  = (mask->vector? mask->vector->size : 0 )+
            (mask->matrix? mask->matrix->size1 * mask->matrix->size2 : 0);
    p->gradient_for_base       = apop_data_alloc(size,0,0); //use full size.
    if (mask->vector)
        for(i=0; i< mask->vector->size; i++)
            if (apop_data_get(mask, i, -1))
                size    --;
    if (mask->matrix)
        for(i=0; i< mask->matrix->size1; i++)
            for(j=0; j< mask->matrix->size2; j++)
                if (apop_data_get(mask, i, j))
                    size    --;
    if (!size)
        apop_error(0,'s',"You're asking me to estimate a model where every single parameter is fixed.\n");
    model_out->vbase            = size;
    model_out->m1base           = 
    model_out->m2base           = 0;
  apop_mle_settings   *mle_out    = apop_mle_settings_alloc(data, *model_out);
    mle_out->model->model_settings   = p;
    apop_data *startingpt = apop_data_alloc(size, 0, 0);
    fixed_params_pack(paramvals, startingpt, mask, 0);
    mle_out->starting_pt        = malloc(sizeof(double) * size);
    memcpy(mle_out->starting_pt, startingpt->vector->data, sizeof(double) * size);
    apop_data_free(startingpt);
    p->mask                     = mask;
    p->paramvals                = paramvals;
    p->base_model               = apop_model_copy(model_in);
    p->base_model->parameters   = NULL;
    return mle_out;
}
