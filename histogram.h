#ifndef APOP_HISTOGRAM_H
#define APOP_HISTOGRAM_H
#undef __BEGIN_DECLS    /* extern "C" stuff cut 'n' pasted from the GSL. */
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
#include "variadic.h"

__BEGIN_DECLS
apop_model * apop_histogram_vector_reset(apop_model *template, gsl_vector *indata);
APOP_VAR_DECLARE apop_model * apop_histogram_model_reset(apop_model *template, apop_model *m, long int draws, gsl_rng *rng);
apop_data * apop_histograms_test_goodness_of_fit(apop_model *h0, apop_model *h1);
apop_data * apop_test_kolmogorov(apop_model *m1, apop_model *m2);
void apop_histogram_normalize(apop_model *m);
__END_DECLS

#endif
