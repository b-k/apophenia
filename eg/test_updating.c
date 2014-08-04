#include <apop.h>

//For the test suite.
void distances(gsl_vector *v1, gsl_vector *v2, double tol){
    double error = apop_vector_distance(v1, v2, .metric='m');
    double updated_size = apop_vector_sum(v1);
    Apop_stopif(error/updated_size > tol, exit(1), 0, "The error is %g, which is too big.", error/updated_size);
}

int main(){
    gsl_rng *r = apop_rng_alloc(2468);
    double binom_start = 0.6;
    double beta_start_a = 0.3;
    double beta_start_b = 0.5;
    double n = 4000;
    //First, the easy estimation using the conjugate distribution table.
    apop_model *bin = apop_model_set_parameters(apop_binomial, n, binom_start);
    apop_model *beta = apop_model_set_parameters(apop_beta, beta_start_a, beta_start_b);
    apop_model *updated = apop_update(.prior= beta, .likelihood=bin,.rng=r);

    //Now estimate via MCMC. 
    //Requires a one-parameter binomial, with n fixed,
    //and a data set of n data points with the right p.
    apop_model *bcopy = apop_model_set_parameters(apop_binomial, n, GSL_NAN);
    apop_data *bin_draws = apop_data_falloc((1,2), n*(1-binom_start), n*binom_start);
    bin = apop_model_fix_params(bcopy);
    Apop_settings_add_group(beta, apop_mcmc, .burnin=.2, .periods=1e5);

    apop_model *out_h = apop_update(bin_draws, beta, bin, NULL);
    apop_model *out_beta = apop_estimate(out_h->data, apop_beta);

    //Finally, we can compare the conjugate and Gibbs results:
    distances(updated->parameters->vector, out_beta->parameters->vector, 0.01);

    //The apop_update function used apop_model_metropolis to generate
    //a batch of draws. Let's use apop_model_metropolis_draw to get draws.
    int i, draws = 1.3e5;
    apop_data *d = apop_data_alloc(draws, 1);
    for(i=0; i < draws; i ++)
        apop_draw(apop_data_ptr(d, i, 0), r, out_h);
    apop_model *drawn = apop_estimate(d, apop_beta);
    distances(updated->parameters->vector, drawn->parameters->vector, 0.02);
}
