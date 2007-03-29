#include <apop.h>
/** \file apop_model_fix_params.c 
 Fix the parameters of a model.

 There's only one function here. Its header is in asst.h
 
 (c) 2007 Ben Klemens. Licensed under the GNU GPL v2. */




typedef struct apop_model_fixed_params{
    apop_data *paramvals;
    apop_data *mask;
    apop_data *filled_beta;
    apop_model *m;
    void * base_model_params;
    apop_ep *selfep;
    struct apop_model_fixed_params *selfparams;
} apop_model_fixed_params;

static apop_model fixed_param_model;

void  fixed_param_unpack(gsl_vector *in, apop_model_fixed_params *p){
  int ctr   = 0;
  int i, j;
    if (!p->filled_beta)
        p->filled_beta  = apop_data_copy(p->paramvals);
    if (p->mask->vector)
        for(i=0; i< p->mask->vector->size; i++)
            if (!apop_data_get(p->mask, i, -1))
                apop_data_set(p->filled_beta, i, -1, gsl_vector_get(in, ctr++));
    if (p->mask->matrix)
        for(i=0; i< p->mask->matrix->size1; i++)
            for(j=0; j< p->mask->matrix->size2; j++)
                if (!apop_data_get(p->mask, i, j))
                    apop_data_set(p->filled_beta, i, j, gsl_vector_get(in, ctr++));
}

static double i_ll(const apop_data *beta, apop_data *d, void *params){
  apop_model_fixed_params *p    = ((apop_ep*)params)->more;
    fixed_param_unpack(beta->vector, p);
    return p->m->log_likelihood(p->filled_beta, d, p->base_model_params);
}

static double i_p(const apop_data *beta, apop_data *d, void *params){
  apop_model_fixed_params *p    = ((apop_ep*)params)->more;
    fixed_param_unpack(beta->vector, p);
    return p->m->p(p->filled_beta, d, p->base_model_params);
}

static void i_score(const apop_data *beta, apop_data *d, gsl_vector *gradient, void *params){
  apop_model_fixed_params *p    = ((apop_ep*)params)->more;
    fixed_param_unpack(beta->vector, p);
    p->m->score(p->filled_beta, d, gradient, p->base_model_params);
}

static double  i_constraint(const apop_data *beta, void * d, apop_data *returned_beta, void *params){
  apop_model_fixed_params *p    = ((apop_ep*)params)->more;
    fixed_param_unpack(beta->vector, p);
    return p->m->constraint(p->filled_beta, d, returned_beta, p->base_model_params);
}

static void i_draw(double *out, apop_data *beta, apop_ep *eps, gsl_rng* r){
  apop_model_fixed_params *p    = eps->more;
    fixed_param_unpack(beta->vector, p);
    p->m->draw(out, p->filled_beta, p->base_model_params, r);
}

static apop_estimate *fixed_est(apop_data * data, void *params){
  apop_model_fixed_params *p    = ((apop_ep*)params)->more;
    apop_estimate *e = apop_maximum_likelihood(data,  *(p->selfep->model), p->selfep);
    fixed_param_unpack(e->parameters->vector, p);
    apop_data_free(e->parameters);
    e->parameters   = apop_data_copy(p->filled_beta);
    return e;
}

static apop_model fixed_param_model = {"Fill me", 0, 0, 0, fixed_est, i_p, i_ll, i_score, i_constraint, i_draw};


/* Produce a model based on another model, but with some of the 
  parameters fixed at a given value. You input one \c apop_data set with the values at which you want 
  the values fixed, and one apop_data set which is zero at the location of the free parameters and 
  one at the location of the fixed parameters.

  You also need to input the base model, and two sets of parameters:
  the parameters for your model, which will be blindly passed on, and MLE-style parameters for the 
  \c estimate method. The \c estimate method always uses an MLE.

  The output is an \c apop_ep. That structure includes a \c apop_model element, which you can call directly.
  Here is a sample program. It produces a few thousand draws from a multivariate normal distribution,
  and then tries to recover the means given a var/covar matrix fixed at the correct variance.

  \code
#include <apop.h>
int main(){
  apop_data *pv   = apop_data_alloc(2,2,2);
  apop_data *mask = apop_data_calloc(2,2,2);
  size_t    i, ct = 1000;
  apop_data *d  = apop_data_alloc(0,ct,2);
  apop_ep *ep   = apop_ep_alloc();
  gsl_rng *r    = apop_rng_alloc(10);
  double    draw[2];
    apop_data_set(pv, 0, -1, 8);
    apop_data_set(pv, 1, -1, 2);
    gsl_matrix_set(pv->matrix, 0,0,1);
    gsl_matrix_set(pv->matrix, 1,1,1);
    gsl_matrix_set(pv->matrix, 1,0,0.5);
    gsl_matrix_set(pv->matrix, 0,1,0.5);
    ep->method      = 0;
    for(i=0; i< ct; i++){
        apop_multivariate_normal.draw(draw, pv, NULL, r);
        apop_data_set(d, i, 0, draw[0]);
        apop_data_set(d, i, 1, draw[1]);
    }

    gsl_matrix_set_all(mask->matrix, 1);
    apop_ep *mep1   = apop_model_fix_params(pv, mask, apop_multivariate_normal, NULL, ep);
    apop_estimate   *e1  = mep1->model->estimate(d, mep1);
    gsl_vector_sub(e1->parameters->vector, pv->vector);
    assert(apop_vector_sum(e1->parameters->vector) < 1e-2);
}

  \endcode
  
  \param paramvals  An \c apop_data set with the values of the variables to be fixed.
  \param mask       Set to zero for free values, one for values to be read from \c paramvals
  \param model_in   The base model
  \param params_for_model   the parameters to be passed to the model (if any)
  \param mle_params The parameters to send to the MLE that will be used for the \c estimate method.
  */

apop_ep *apop_model_fix_params(apop_data *paramvals, apop_data *mask, apop_model model_in, apop_ep *params_for_model, apop_ep *mle_params){
  apop_model_fixed_params *p    = malloc(sizeof(*p));
  apop_model *model_out = apop_model_copy(fixed_param_model);
  int        i, j;
    if (!model_in.estimate) model_out->estimate = NULL;
    if (!model_in.p) model_out->p = NULL;
    if (!model_in.log_likelihood) model_out->log_likelihood = NULL;
    if (!model_in.score) model_out->score = NULL;
    if (!model_in.constraint) model_out->constraint = NULL;
    if (!model_in.draw) model_out->draw = NULL;
    snprintf(model_out->name, 100, "%s, with some params fixed", model_in.name);
    model_out->more     = p;
    model_out->vsize    = (mask->vector? mask->vector->size: 0)+
                           (mask->matrix? mask->matrix->size1+ mask->matrix->size2: 0);
    model_out->msize1   = 
    model_out->msize2   = 0;
    if (mask->vector)
        for(i=0; i< mask->vector->size; i++)
            if (apop_data_get(mask, i, -1))
                model_out->vsize    --;
    if (mask->matrix)
        for(i=0; i< mask->matrix->size1; i++)
            for(j=0; j< mask->matrix->size2; j++)
                if (apop_data_get(mask, i, j))
                    model_out->vsize    --;
    p->mask             = mask;
    p->paramvals        = paramvals;
    p->base_model_params= params_for_model;
    p->selfep           = mle_params ? mle_params : apop_ep_alloc();
    p->selfep->more     = p;
    p->m                =  apop_model_copy(model_in);
    p->filled_beta      =  NULL;
    p->selfep->model    =  apop_model_copy(*model_out);
    return p->selfep;
}
