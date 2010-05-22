void test_updating(gsl_rng *r){
  double binom_start = 0.6;
  double beta_start_a = 0.3;
  double beta_start_b = 0.5;
  int i, draws = 15000;
  double n = 800;
  //First, the easy estimation using the conjugate distribution table.
  apop_model *bin = apop_model_set_parameters(apop_binomial, n, binom_start);
  apop_model *beta = apop_model_set_parameters(apop_beta, beta_start_a, beta_start_b);
    apop_model *updated = apop_update(.prior= beta, .likelihood=bin,.rng=r);
    apop_model *bcopy = apop_model_set_parameters(apop_binomial, n, GSL_NAN);

    //Now estimate via Gibbs sampling. 
    //Requires a one-parameter binomial, with n fixed,
    //and a data set of n data points with the right p.
    apop_data *bin_draws = apop_data_calloc(0, n, 1);
    for(i=0; i < n*binom_start; i ++)
        apop_matrix_increment(bin_draws->matrix, i, 0);
    bin = apop_model_fix_params(bcopy);
    apop_model *out_h = apop_update(bin_draws, beta, bin, NULL);

    //We now have a histogram of values for p. What's the closest beta
    //distribution?
    apop_data *d = apop_data_alloc(0, draws, 1);
    for(i=0; i < draws; i ++)
        apop_draw(apop_data_ptr(d, i, 0), r, out_h);
    apop_model *out_beta = apop_estimate(d, apop_beta);

    //Finally, we can compare the conjugate and Gibbs results:
    double updated_size = apop_vector_sum(updated->parameters->vector);
    double error = apop_vector_grid_distance(updated->parameters->vector, out_beta->parameters->vector);
    Apop_assert_s(error/updated_size < 0.1, "The error is %g, which is too big.", error/updated_size);
}
