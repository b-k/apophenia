gsl_histogram * apop_vector_to_histogram(gsl_vector *data, int bins);
gsl_histogram_pdf * apop_vector_to_cmf(gsl_vector *data, int bins);
gsl_histogram ** apop_vectors_to_histograms(gsl_vector *v1, gsl_vector *v2, int bins);
gsl_histogram * apop_model_to_histogram(apop_model *m, gsl_histogram *h, int draws, double *params, gsl_rng *r);
apop_data *apop_model_test_goodness_of_fit(gsl_vector *v1, apop_model *m,
int bins, long int draws, double *params, gsl_rng *r);
apop_data *apop_vectors_test_goodness_of_fit(gsl_vector *v0, gsl_vector *v1);
apop_data *apop_histograms_test_goodness_of_fit(gsl_histogram *h0, gsl_histogram *h1, int bins);
