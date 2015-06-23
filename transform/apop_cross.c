/* Cross product of distributions.
 Copyright (c) 2013 by Ben Klemens.  Licensed under the GPLv2; see COPYING.  

\amodel apop_cross A cross product of models.  Generate via \ref apop_model_cross .

For the case when you need to bundle two uncorrelated models into one larger model. For example, the prior for a multivariate normal (whose parameters are a vector of means and a covariance matrix) is a Multivariate Normal-Wishart pair.

\adoc    Input_format     There are two means of handling the input format. If the settings group attached to the data set has a non-\c NULL \c splitpage element, then 
append the second data set as an additional page to the first data set, and name the second set with the name you listed in \c splitpage; see the example.  

If \c splitpage is \c NULL, then I will send the same data set to both models.

\adoc    Settings   \ref apop_cross_settings

\adoc    Parameter_format  
currently \c NULL; check the sub-models for their parameters.

\adoc    For an example, see \ref apop_model_cross .
*/

#include "apop_internal.h"

Apop_settings_init(apop_cross, )
Apop_settings_copy(apop_cross, )
Apop_settings_free(apop_cross, )


/* For (almost) all methods, the workings are:
  --get the settings group; fail if missing
  --find the first and second data sets. This may require changing a 
    ->more pointer in the input data set to split it into two parts
  --call the submodels using our two data sets 
  --restore that ->more pointer, if needed.
   */

typedef struct {
    apop_data *d1, *d2, *dangly_bit;
    _Bool need_to_free;
} twop_s;


/* A model must accept data in the form of each observation being a single row of
   data---i.e., each row is what you'd get from apop_draw.
    If the length of the data matrix is longer than the first model's stated msize2,
    then we have to unpack the draw into two new apop_data sets.
*/
twop_s unpack_a_draw(apop_data *d, apop_cross_settings *s){
    twop_s out = (twop_s){
        .d1 = (s->model1 ? apop_data_alloc(s->model1->vsize, s->model1->msize1, s->model1->msize2) : NULL),
        .d2 = (s->model2 ? apop_data_alloc(s->model2->vsize, s->model2->msize1, s->model2->msize2) : NULL),
        .need_to_free=1};
    int len1 = s->model1->vsize + s->model1->msize1 * s->model1->msize2;
    for (int i=0; i< d->matrix->size1; i++){
        apop_data_unpack(Apop_subvector(Apop_rv(d, i), 0, len1), Apop_r(out.d1, i));
        apop_data_unpack(Apop_subvector(Apop_rv(d, i), len1, d->matrix->size1), Apop_r(out.d2, i));
      }
    return out;
}

static twop_s get_second(apop_data *d, char *splitpage, apop_cross_settings *s){
    twop_s out = {.d1=d, .d2=d};
    if (splitpage) {
        if (d->matrix && (d->matrix->size2 > s->model1->msize2))
            return unpack_a_draw(d, s);
        apop_data *ctr = d;
        if (!ctr ||(ctr->names && ctr->names->title && !strcasecmp(ctr->names->title, splitpage))){
            out.d1 = NULL;
            out.d2 = d;
            return out;
        }
        for ( ; ctr->more && (!ctr->more->names || !ctr->more->names->title || strcasecmp(ctr->more->names->title, splitpage)); ) 
            ctr = ctr->more; 
        out.d2 = ctr->more;
        out.dangly_bit = ctr;
        ctr->more = NULL; //the only change to the original data set.
    }
    return out;
}

//post-work cleanup: reattach 2nd data set, or if we built a newly-shaped data set, free it.
static void repaste(twop_s dd){
    if (dd.dangly_bit) dd.dangly_bit->more = dd.d2;
    if (dd.need_to_free) {apop_data_free(dd.d1); apop_data_free(dd.d2);}
}

#define check_settings(ret) Apop_stopif(!s, m->error='s'; return ret, 0, "This model wasn't set up right. Maybe use apop_model_cross to set up your pair of models.");

void cross_print(apop_model *m, FILE *out){
    apop_cross_settings *s = Apop_settings_get_group(m, apop_cross);
    check_settings( );
    apop_model *m1 = s->model1;
    apop_model *m2 = s->model2;
    fprintf(out, "Cross product of %s and %s models:\n", m1->name, m2->name);
    apop_model_print(m1, out);
    apop_model_print(m2, out);
}
    
#define Preliminaries(ret)          \
    apop_cross_settings *s = Apop_settings_get_group(m, apop_cross);    \
    check_settings(ret);        \
    twop_s datas = get_second(d, s->splitpage, s);

static void cross_est(apop_data *d, apop_model *m){
    Preliminaries();

    s->model1 = apop_estimate(datas.d1, s->model1);
    s->model2 = apop_estimate(datas.d2, s->model2);

    repaste(datas);
}

static long double cross_ll(apop_data *d, apop_model *m){
    Preliminaries(GSL_NAN);

    double out =  apop_log_likelihood(datas.d1, s->model1)
                 +apop_log_likelihood(datas.d2, s->model2);
    repaste(datas);
    return out;
}

static long double cross_p(apop_data *d, apop_model *m){
    Preliminaries(GSL_NAN)

    double out =  apop_p(datas.d1, s->model1) *apop_p(datas.d2, s->model2);
    repaste(datas);
    return out;
}

static int cross_draw(double *d, gsl_rng *r, apop_model *m){
    apop_cross_settings *s = Apop_settings_get_group(m, apop_cross);
    check_settings(1); 
    Apop_stopif(apop_draw(d, r, s->model1), return 1, 0, "draw from first model failed.");
    double *d2 = d+ s->model1->dsize;
    Apop_stopif(apop_draw(d2, r, s->model2), return 1, 0, "draw from second model failed.");
    return 0;
}

apop_model *apop_cross = &(apop_model){"Cross product of models", .p=cross_p, .log_likelihood=cross_ll, 
    .estimate=cross_est, .draw=cross_draw
};

apop_model *apop_model_cross_base(apop_model *mlist[]){
    apop_model_print_vtable_add(cross_print, apop_cross);
    Apop_stopif(!mlist[0], apop_model *oute = apop_model_copy(apop_cross); oute->error='i', 
                            0, "No inputs. Returning blank model with outmodel->error=='n'.");
    Apop_stopif(!mlist[1], return apop_model_copy(mlist[1]), 
                            1, "Only one model input; returning a copy of that model.");
    apop_model *m2 = mlist[2] ? apop_model_cross_base(mlist+1): mlist[1];
    apop_model *out = apop_model_copy(apop_cross);
    Apop_model_add_group(out, apop_cross, .model1=mlist[0], .model2=m2);
    if (mlist[0]->dsize >=0 && m2->dsize >=0) out->dsize = mlist[0]->dsize + m2->dsize;
    return out;
}

/** \def apop_model_cross
Generate a model consisting of the cross product of several independent models. The output \ref apop_model
is a copy of \ref apop_cross; see that model's documentation for details.

\li If you input only one model, return a copy of that model; print a warning iff <tt>apop_opts.verbose >= 1</tt>.
\exception error=='n' First model input is \c NULL.

Examples:

\include cross_models.c
*/
