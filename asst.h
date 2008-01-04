/** \file asst.h		Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#ifndef __apop_asst__
#define __apop_asst__

#include <assert.h>
#include <apophenia/types.h>
#include <apophenia/likelihoods.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

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

double apop_generalized_harmonic(int N, double s);
void apop_error(int level, char stop, char *message, ...);

apop_model * apop_update(apop_data *data, apop_model prior, apop_model likelihood, 
                        apop_data *starting_pt, gsl_rng *r, int periods, double burnin, int histosegments);

apop_data * apop_test_ANOVA_independence(apop_data *d);

int apop_system(const char *fmt, ...) __attribute__ ((format (printf,1,2)));

gsl_vector * apop_vector_moving_average(gsl_vector *, size_t);

/** A convenient front for \ref apop_error, that tests the first element and
 basically runs \ref apop_error if it is false. See also \ref apop_assert_void and \ref apop_error*/
#define apop_assert(test, returnval, level, stop, ...) do \
    if (!(test)) {  \
        if (apop_opts.verbose >= level) { fprintf(stderr, "%s: ", __func__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");}   \
        if (stop == 's' || stop == 'h') assert(test);   \
        return returnval;  \
} while (0);

/** Like \ref apop_assert, but useful in void functions. */
#define apop_assert_void(test,  level, stop, ...) do \
    if (!(test)) {  \
        if (apop_opts.verbose >= level) { fprintf(stderr, "%s: ", __func__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");}   \
        if (stop == 's' || stop == 'h') assert(test);   \
} while (0);


__END_DECLS
#endif
