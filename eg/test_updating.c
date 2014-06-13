#include <apop.h>
int main(){
    gsl_rng *r = apop_rng_alloc(2468);
    double binom_start = 0.6;
    double beta_start_a = 0.3;
    double beta_start_b = 0.5;
    int i, draws = 1500;
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
    Apop_settings_add_group(beta, apop_mcmc, .burnin=.1, .periods=1e4);
    apop_model *out_h = apop_update(bin_draws, beta, bin, NULL);

    //We now have a histogram of values for p. What's the closest beta
    //distribution?
    apop_data *d = apop_data_alloc(draws, 1);
    for(i=0; i < draws; i ++)
        apop_draw(apop_data_ptr(d, i, 0), r, out_h);
    apop_model *out_beta = apop_estimate(d, apop_beta);
    //Finally, we can compare the conjugate and Gibbs results:
    apop_vector_normalize(updated->parameters->vector);
    apop_vector_normalize(out_beta->parameters->vector);
    double error = apop_vector_distance(updated->parameters->vector, out_beta->parameters->vector, .metric='m');
    double updated_size = apop_vector_sum(updated->parameters->vector);
    Apop_assert(error/updated_size < 0.01, "The error is %g, which is too big.", error/updated_size);
}
