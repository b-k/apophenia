/* Stacking distributions.
 Copyright (c) 2013 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  

\amodel apop_stack A stack of models.

For the case when you need to bundle two uncorrelated models into one larger model. For example, the prior for a multivariate normal (whose parameters are a vector of means and a covariance matrix) is a Multivariate Normal-Wishart pair.


\adoc    Input_format     There are two means of handling the input format. If the settings group attached to the data set has a non-\c NULL \c splitpage element, then 
Append the second data set as an additional page to the first data set, and name the second set with the name you listed in \c splitpage; see the example.  

If \c splitpage is \c NULL, then I will send the same data set to both models.


\adoc    Settings   \ref apop_stack_settings

\adoc    Parameter_format  
currently \c NULL; check the sub-models for their parameters.

\adoc    Example 
\include stack_models.c
*/

#include "apop_internal.h"

Apop_settings_init(apop_stack, )
Apop_settings_copy(apop_stack, )
Apop_settings_free(apop_stack, )


/* For (almost) all methods, the workings are:
  --get the settings group; fail if missing
  --find the first and second data sets. This may require changing a 
    ->more pointer in the input data set to split it into two parts
  --call the submodels using our two data sets 
  --restore that ->more pointer, if needed.
   */

typedef struct {
    apop_data *d1, *d2, *dangly_bit;
} twop_s;

static twop_s get_second(apop_data *d, char *splitpage){
    twop_s out = {.d1=d, .d2=d};
    if (splitpage) {
        apop_data *ctr = d;
        if (!ctr ||(ctr->names && !strcasecmp(ctr->names->title, splitpage))){
            out.d1 = NULL;
            out.d2 = d;
            return out;
        }
        for ( ; ctr->more && (!ctr->more->names || strcasecmp(ctr->more->names->title, splitpage)); ) 
            ctr = ctr->more; 
        out.d2 = ctr->more;
        out.dangly_bit = ctr;
        ctr->more = NULL; //the only change to the original data set.
    }
    return out;
}

static void repaste(twop_s dd){
    if (dd.dangly_bit) dd.dangly_bit->more = dd.d2;
}

#define check_settings(ret) Apop_stopif(!s, m->error='s'; return ret, 0, "This model wasn't set up right. Maybe use apop_model_stack to set up your pair of models.");

#define Preliminaries(ret)          \
    apop_stack_settings *s = Apop_settings_get_group(m, apop_stack);    \
    check_settings(ret);        \
    twop_s datas = get_second(d, s->splitpage);

static apop_model *stack_est(apop_data *d, apop_model *m){
    Preliminaries(m);

    s->model1 = apop_estimate(datas.d1, *s->model1);
    s->model2 = apop_estimate(datas.d2, *s->model2);

    repaste(datas);
    return m;
}

static double stack_ll(apop_data *d, apop_model *m){
    Preliminaries(GSL_NAN);

    double out =  apop_log_likelihood(datas.d1, s->model1)
                 +apop_log_likelihood(datas.d2, s->model2);
    repaste(datas);
    return out;
}

static double stack_p(apop_data *d, apop_model *m){
    Preliminaries(GSL_NAN)

    double out =  apop_p(datas.d1, s->model1) *apop_p(datas.d2, s->model2);
    repaste(datas);
    return out;
}

void stack_draw(double *d, gsl_rng *r, apop_model *m){
    apop_stack_settings *s = Apop_settings_get_group(m, apop_stack);
    check_settings(); 
    apop_draw(d, r, s->model1);
    double *d2 = d+ s->model1->dsize;
    apop_draw(d2, r, s->model2);
}

apop_model apop_stack = {"Stack of models", .p=stack_p, .log_likelihood=stack_ll, 
    .estimate=stack_est, .draw=stack_draw
};


/** Generate a model consisting of two models bound together. The output \ref apop_model
 is a copy of \ref apop_stack; see that model's documentation for details.

 Do I even need this function?

*/
apop_model *apop_model_stack(apop_model *m1, apop_model *m2){
    apop_model *out = apop_model_copy(apop_stack);
    Apop_model_add_group(out, apop_stack, .model1=m1, .model2=m2);
    if (m1->dsize >=0 && m2->dsize >=0) out->dsize = m1->dsize + m2->dsize;
    return out;
}
