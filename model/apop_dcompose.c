#include <apop.h>

/* \amodel apop_dcomposition Use random draws or parameter estimates output from a first model as input data for a second model.

\li This is still in beta. Expect the interface to change.

\adoc Input_format The input data is sent to the first model, so use the input format for that model.
\adoc Settings   \ref apop_composition_settings

*/

typedef struct {
    apop_model *generator_m;
    apop_model *ll_m;
    gsl_rng *rng;
    int draw_ct;
} apop_composition_settings;/**< All of the elements of this struct should be considered private.*/

Apop_settings_copy(apop_composition,)

Apop_settings_free(apop_composition,)

Apop_settings_init(apop_composition,)

#define Get_cs(inmodel, outval) \
    apop_composition_settings *cs = Apop_settings_get_group(inmodel, apop_composition); \
    Apop_stopif(!cs, return outval, 0, "At this point, I expect your model to" \
            "have an apop_compose settings group. Perhaps set up the " \
            "model using apop_model_dcompose");

static void compose_prep(apop_data *d, apop_model *m){
    Get_cs(m, )
    apop_prep(d, cs->ll_m);
    m->parameters=cs->ll_m->parameters;
    m->parameters = cs->ll_m->parameters;
    m->vbase = cs->ll_m->vbase;
    m->m1base = cs->ll_m->m1base;
    m->m2base = cs->ll_m->m2base;
    m->dsize = cs->ll_m->dsize;
}

static double compose_ll(apop_data *in, apop_model*composition){
    Get_cs(composition, GSL_NAN)
    apop_data *draws = apop_data_alloc(cs->draw_ct, cs->generator_m->dsize);
    for (int i=0; i< cs->draw_ct; i++){
        Apop_row(draws, i, onerow);
        apop_draw(onerow->data, cs->rng, cs->generator_m);
    }
    //apop_model *est = apop_estimate(draws, *cs->ll_m);
    double ll = apop_log_likelihood(draws, cs->ll_m);
    //apop_model_free(est);
    apop_data_free(draws);
    return ll;
}

apop_model apop_composition = {"Data-composed model", .prep=compose_prep, .log_likelihood=compose_ll};


/** <em>Data composition</em> is using either random draws or parameter estimates from
the output of one model as the input data for another model. 

\li The \ref apop_dcomposition model relies on the \ref apop_composition_settings struct, qv.

\return An \ref apop_model that is a copy of \ref apop_composition.

\ingroup model_transformations
*/
apop_model *apop_model_dcompose(apop_model *datamod, apop_model *post, int draw_ct, gsl_rng *rng){
    apop_model *out = apop_model_copy(apop_composition);
    Apop_model_add_group(out, apop_composition, 
                .generator_m = datamod,
                .ll_m = post,
                .draw_ct = draw_ct,
                .rng = rng ? rng : apop_rng_alloc(apop_opts.rng_seed++));
    return out;
}
