/* A Lognormal distribution is a transform of the Normal distribution, where 
 the data space of the Normal is exponentiated. Thus, to get back to the original data space, take the log of the current data.
 */
#include <apop.h>
#define Diff(a, b) assert(fabs((a)-(b)) < 1e-2);

apop_data *draw_exponentiated_normal(double mu, double sigma, double draws){
    apop_model *n01 = apop_model_set_parameters(apop_normal, mu, sigma);
    apop_data *d = apop_data_alloc(draws);
    gsl_rng *r = apop_rng_alloc(13);
    for (int i=0; i< draws; i++) apop_draw(gsl_vector_ptr(d->vector,i), r, n01);
    apop_vector_exp(d->vector);
    return d;
}

// The transformed->base function and its derivative for the Jacobian:
apop_data *rev(apop_data *in){ return apop_map(in, .fn_d=log, .part='a'); }

/*The derivative of the transformed_to_base function. */
double inv(double in){return 1./in;} 
double rev_j(apop_data *in){ return fabs(apop_map_sum(in, .fn_d=inv, .part='a')); }

int main(){
    apop_model *ct = apop_model_coordinate_transform(
            .transformed_to_base= rev, .jacobian_to_base=rev_j,
            .base_model=apop_normal);
    Apop_model_add_group(ct, apop_parts_wanted);//Speed up the MLE.

    //make fake data
    double mu=2, sigma=1;
    apop_data *d = draw_exponentiated_normal(mu, sigma, 2e5);

    //If we correctly replicated a Lognormal, mu and sigma will be right:
    apop_model *est = apop_estimate(d, ct);
    apop_model_free(ct);
    Diff(apop_data_get(est->parameters, 0, -1), mu);
    Diff(apop_data_get(est->parameters, 1, -1), sigma);

    /*The K-L divergence between our Lognormal and the stock Lognormal
      should be small. I try it with both the original params and the estimated ones. */
    apop_model *ln = apop_model_set_parameters(apop_lognormal, mu, sigma);
    apop_model *ln2 = apop_model_copy(apop_lognormal);
    ln2->parameters = est->parameters;
    Diff(apop_kl_divergence(ln, ln2,.draw_ct=1000), 0);
    Diff(apop_kl_divergence(ln, est,.draw_ct=1000), 0);
}
