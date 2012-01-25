void test_kl_divergence(gsl_rng *r){
    int i, empirical_size = 5e3;
    apop_model *expo = apop_model_set_parameters(apop_exponential, 1.7);
    assert (apop_kl_divergence(expo, expo) < 1e-4);
    apop_data *empirical = apop_data_alloc(0,  empirical_size, 1);
    for (i=0; i<empirical_size; i++)
        apop_draw(apop_data_ptr(empirical, i, 0), r, expo);
    apop_model *pmf = apop_estimate(empirical, apop_pmf);
    assert(apop_kl_divergence(expo,pmf) < 1e-4);
    apop_data_free(empirical);
}
