#ifndef APOP_BOOT_H
#define APOP_BOOT_H

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

apop_data * apop_jackknife_cov(apop_data *data, apop_model model);
APOP_VAR_DECLARE apop_data * apop_bootstrap_cov(apop_data *data, apop_model model, gsl_rng* rng, int iterations);
gsl_rng *apop_rng_alloc(int seed);

__END_DECLS
#endif
