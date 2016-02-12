#include <apop.h>

long double fake_p (apop_data *d, apop_model *m){
    return apop_pmf->p(d, m);
}

int main(){
    gsl_rng *r = apop_rng_alloc(2312311);
    int empirical_size = 5e3;
    apop_model *expo = apop_model_set_parameters(apop_exponential, 1.7);
    //divergence from self should be zero.
    assert (fabs(apop_kl_divergence(expo, expo)) < 1e-4);

    apop_data *empirical = apop_model_draws(expo, .count=empirical_size);
    //Double the odds of half the data, so likelihoods aren't uniform.
    int half =empirical_size/2;
    apop_data *start = Apop_rs(empirical, 0, half);
    empirical = apop_data_stack(empirical, start);
    apop_data_pmf_compress(empirical);

    //Compare the PMF calculator to the everything else calculator
    apop_model *pmf = apop_estimate(empirical, apop_pmf);
    double div= apop_kl_divergence(pmf,expo);
    pmf->p = fake_p;
    double div2= apop_kl_divergence(pmf,expo);
    printf("%g %g\n", div, div2);
    assert(fabs(div-div2)<5e-3);
    apop_data_free(empirical);
}
