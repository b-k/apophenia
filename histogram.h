#undef __BEGIN_DECLS    /* extern "C" stuff cut 'n' pasted from the GSL. */
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS
apop_model *apop_histogram_refill_with_vector(apop_model *template, gsl_vector *indata);
apop_model *apop_histogram_refill_with_model(apop_model *template, apop_model *m, long int draws, gsl_rng *r);
apop_data *apop_histograms_test_goodness_of_fit(apop_model *h0, apop_model *h1);
apop_data *apop_test_kolmogorov(apop_model *m1, apop_model *m2);
void apop_histogram_normalize(apop_model *m);
__END_DECLS
