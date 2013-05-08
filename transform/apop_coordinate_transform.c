#include "apop_internal.h"

/* \amodel apop_coordinate_transform Apply a coordinate transformation of the data to
produce a distribution over the transformed data space. This is sometimes called a Jacobian transformation.

\li This is still in beta. Expect the interface to change.

Here is an example that replicates the Lognormal distribution.

\include jacobian.c

\adoc Input_format The input data is sent to the first model, so use the input format for that model.
\adoc Settings   \ref apop_ct_settings

*/

typedef struct {
    apop_data *(*base_to_transformed)(apop_data*);
    apop_data *(*transformed_to_base)(apop_data*);
    double (*jacobian_to_base)(apop_data*);
    apop_model *base_model;
} apop_ct_settings;/**< All of the elements of this struct should be considered private.*/

Apop_settings_declarations(apop_ct)
Apop_settings_init(apop_ct,)//use defaults.
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

static double ct_ll(apop_data *indata, apop_model* mj){
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

apop_model apop_coordinate_transform = {"Jacobian-transformed model", .log_likelihood=ct_ll, .prep=jacobian_prep};

typedef apop_data *(*d_to_d)(apop_data*);
typedef double (*d_to_f)(apop_data*);

/** Build an \ref apop_coordinate_transform model, qv.

\return An \ref apop_model that is a copy of \ref apop_coordinate_transform and is appropriately set up.

\ingroup model_transformations
*/
APOP_VAR_HEAD apop_model* apop_model_coordinate_transform(apop_model *base_model, d_to_d transformed_to_base, d_to_d base_to_transformed, d_to_f jacobian_to_base){ 
    apop_model * apop_varad_var(base_model, NULL);
    Apop_stopif(!base_model, apop_model *out=apop_model_copy(apop_coordinate_transform); out->error='b'; return out,
                    0, "No base model specified. Returning a composition model with error='b'.");
    //type-checked on the way in:
    void * apop_varad_var(transformed_to_base, NULL);
    void * apop_varad_var(base_to_transformed, NULL);
    void * apop_varad_var(jacobian_to_base, NULL);
APOP_VAR_ENDHEAD
    apop_model *out = apop_model_copy(apop_coordinate_transform);
    Apop_model_add_group(out, apop_ct,
        .base_to_transformed = base_to_transformed,
        .transformed_to_base= transformed_to_base,
        .jacobian_to_base = jacobian_to_base,
        .base_model= base_model
    );
    return out;
}
