#include "apop_internal.h"

/* \amodel apop_coordinate_transform Apply a coordinate transformation of the data to
produce a distribution over the transformed data space. This is sometimes called a Jacobian transformation.

\li This is still in beta. Expect the interface to change.

Here is an example that replicates the Lognormal distribution.

\include jacobian.c

\adoc Input_format The input data is sent to the first model, so use the input format for that model.
\adoc Settings   \ref apop_ct_settings

*/

Apop_settings_init(apop_ct,
    Apop_stopif(!in.base_model, , 0, "I need a .base_model.");
)
Apop_settings_copy(apop_ct,)
Apop_settings_free(apop_ct,)

#define Get_cs(inmodel, outval) \
    apop_ct_settings *cs = Apop_settings_get_group(inmodel, apop_ct); \
    Apop_stopif(!cs, return outval, 0, "At this point, I expect your model to" \
            "have an apop_ct_settings group.");

static void jacobian_prep(apop_data *d, apop_model *m){
    apop_ct_settings *cs = Apop_settings_get_group(m, apop_ct); 
    Apop_stopif(!cs, m->error='s', 0, "missing apop_ct_settings group. "
            "Maybe initialize this with apop_model_dcompose?");
    apop_prep(d, cs->base_model);
    m->parameters=cs->base_model->parameters;
    m->dsize=cs->base_model->dsize;
}

static long double ct_ll(apop_data *indata, apop_model* mj){
    Get_cs(mj, GSL_NAN)
    Apop_stopif(!cs->base_model, return GSL_NAN, 0, "No base model to transform back to.");
    Apop_stopif(!cs->transformed_to_base, return GSL_NAN, 0, "No reverse transformation function.");
    Apop_stopif(!cs->jacobian_to_base, return GSL_NAN, 0, "No Jacobian for the reverse transformation function, "
                                                          "and I haven't yet implemented the numeric derivative.");

    apop_data *rev = cs->transformed_to_base(indata);
    double ll = apop_log_likelihood(rev, cs->base_model);
    apop_data_free(rev);
    ll += log(cs->jacobian_to_base(indata));
    return ll;
}

apop_model *apop_coordinate_transform = &(apop_model){"Jacobian-transformed model", .log_likelihood=ct_ll, .prep=jacobian_prep};

typedef apop_data *(*d_to_d)(apop_data*);
typedef double (*d_to_f)(apop_data*);

/** \def apop_model_coordinate_transform
Build an \ref apop_coordinate_transform model, qv.

\return An \ref apop_model that is a copy of \ref apop_coordinate_transform and is appropriately set up.

\li Uses the \ref apop_ct_settings group. This macro takes elements of that struct as inputs.

\li This function uses the \ref designated syntax for inputs.
\ingroup model_transformations
*/
