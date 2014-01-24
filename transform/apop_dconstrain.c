#include "apop_internal.h"

/* \amodel apop_dconstrain A model that constrains the base model to within some
data constraint. E.g., truncate \f$P(d)\f$ to zero for all \f$d\f$ outside of a given
constraint. Generate using \ref apop_model_dconstrain .

\li This is still in beta. Expect the interface to change.

The log likelihood works by using the \c base_model log likelihood, and then scaling
it based on the part of the base model's density that is within the constraint. If you
have an easy means of specifying what that density is, please do, as in the example. If
you do not, the log likelihood will calculate it by making \c draw_ct random draws from
the base model and checking whether they are in or out of the constraint. Because this
default method is stochastic, there is some loss of precision, and conjugate gradient
methods may get confused.

Here is an example that makes a few draws and estimations from data-constrained models.

\adoc Examples
\include dconstrain.c

\adoc Input_format That of the base model.
\adoc Parameter_format That of the base model. In fact, the \c parameters element is a pointer to the base model \c parameters, so both are modified simultaneously.
\adoc Settings   \ref apop_dconstrain_settings
*/

#define Get_set(inmodel, outval) \
    apop_dconstrain_settings *cs = Apop_settings_get_group(inmodel, apop_dconstrain); \
    Apop_stopif(!cs, return outval, 0, "At this point, I expect your model to" \
            "have an apop_dconstrain_settings group.");

//what percent of the model density is inside the constraint?
static double get_scaling(apop_model *m){
    Get_set(m, GSL_NAN)
    apop_data *d = apop_data_alloc(1, cs->base_model->dsize);
    int tally = 0;
    for (int i=0; i< cs->draw_ct; i++){
        apop_draw(d->matrix->data, cs->rng, cs->base_model);
        tally += !!cs->constraint(d, cs->base_model);
    }
    apop_data_free(d);
    return (tally+0.0)/cs->draw_ct;
}

Apop_settings_init(apop_dconstrain,
    Apop_stopif(!in.base_model, , 0, "I need a .base_model.");
    Apop_stopif(!in.constraint, , 0, "I need a .constraint.");
    if (!in.draw_ct) out->draw_ct = 1e4;
    if (!in.rng && !in.scaling) out->rng = apop_rng_alloc(apop_opts.rng_seed++);
    if (!in.scaling) out->scaling = get_scaling;
)

Apop_settings_copy(apop_dconstrain,)
Apop_settings_free(apop_dconstrain,)

static void dc_prep(apop_data *d, apop_model *m){
    apop_dconstrain_settings *cs = Apop_settings_get_group(m, apop_dconstrain); 
    Apop_stopif(!cs, m->error='s', 0, "missing apop_dconstrainct_settings group. "
            "Maybe initialize this with apop_model_dcompose?");
    apop_prep(d, cs->base_model);
    m->parameters=cs->base_model->parameters;
    m->constraint=cs->base_model->constraint;
    m->vsize = cs->base_model->vsize;
    m->msize1 = cs->base_model->msize1;
    m->msize2 = cs->base_model->msize2;
    m->dsize=cs->base_model->dsize;
}

/* \adoc RNG Draw from the base model; if the draw is outside the constraint, throw it out and try again. */
static int dc_rng(double *out, gsl_rng *r, apop_model *m){
    Get_set(m, 1)
    gsl_matrix_view mv;
    do {
        apop_draw(out, r, cs->base_model);
        mv = gsl_matrix_view_array(out, 1, cs->base_model->dsize);
    } while (!cs->constraint(&(apop_data){.matrix=&(mv.matrix)}, cs->base_model));
    return 0;
}

static double constr(apop_data *d, void *csin){
    apop_dconstrain_settings* cs = csin;
    return !cs->constraint(d, cs->base_model);
}

static long double dc_ll(apop_data *indata, apop_model* m){
    Get_set(m, GSL_NAN)
    Apop_stopif(!cs->base_model, return GSL_NAN, 0, "No base model.");
    double any_outside = apop_map_sum(indata, .fn_rp=constr, .param=cs);
    if (any_outside) return -INFINITY;

    double scaling = cs->scaling((cs->scaling == get_scaling) ? m : cs->base_model);
    Get_vmsizes(indata); //maxsize
    return apop_log_likelihood(indata, cs->base_model) - log(scaling)*maxsize;
}

apop_model *apop_dconstrain = &(apop_model){"Data-constrained model", .log_likelihood=dc_ll, .draw=dc_rng, .prep=dc_prep};

/** \def apop_model_dconstrain
Build an \c apop_dconstrain model, q.v., which applies a data constraint to the data set. For example, this is how one would truncate a model to have data above zero.

\return An \ref apop_model that is a copy of \ref apop_dconstrain and is appropriately set up.

\li Uses the \ref apop_dconstrain_settings group. This macro takes elements of that struct as inputs.

\li This function uses the \ref designated syntax for inputs.
\ingroup model_transformations
*/
