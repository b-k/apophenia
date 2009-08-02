/** \file asst.h		Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#ifndef __apop_asst__
#define __apop_asst__

#include <assert.h>
#include "types.h"
#include "variadic.h"
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

apop_data * apop_test_anova_independence(apop_data *d);
#define apop_test_ANOVA_independence(d) apop_test_anova_independence(d)

int apop_system(const char *fmt, ...) __attribute__ ((format (printf,1,2)));

gsl_vector * apop_vector_moving_average(gsl_vector *, size_t);
apop_model *apop_histogram_moving_average(apop_model *m, size_t bandwidth);

/** A convenient front for \ref apop_error, that tests the first element and
 basically runs \ref apop_error if it is false. See also \ref Apop_assert_void and \ref apop_error.
 
 Following the tradition regarding assert functions, this is a macro but is not in all caps.

 \param test The expression that you are asserting is nonzero.
 \param returnval If the assertion fails, return this. If you want to halt on error, this is irrelevant, but still has to match your function's return type.
 \param level Print the warning message only if \ref apop_opts_type "apop_opts.verbose" is greater than or equal to this. Zero usually works, but for minor infractions use one.
 \param stop If 's', halt the program (using the standard C \c assert); if 'c', continue by returning the return value and printing an error message if appropriate.
 \param ... The error message in printf form, plus any arguments to be inserted into the printf string. I'll provide the function name and a carriage return.
 */
#define Apop_assert(test, returnval, level, stop, ...) do \
    if (!(test)) {  \
        if (apop_opts.verbose >= level) { fprintf(stderr, "%s: ", __func__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");}   \
        if (stop == 's' || stop == 'h') assert(test);   \
        return returnval;  \
} while (0);

#define apop_assert(test, returnval, level, stop, ...) Apop_assert(test, returnval, level, stop, __VA_ARGS__)
#define APOP_ASSERT(test, returnval, level, stop, ...) Apop_assert(test, returnval, level, stop, __VA_ARGS__)

/** Like \ref Apop_assert, but no return step. It is thus useful in void functions.

 Following the tradition regarding assert functions, this is a macro but is not in all caps.

 \param test The expression that you are asserting is nonzero.
 \param level Print the warning message only if \ref apop_opts_type "apop_opts.verbose" is greater than or equal to this. Zero usually works, but for minor infractions use one.
 \param stop If 's', halt the program (using the standard C \c assert); if 'c', continue by returning the return value and printing an error message if appropriate.
 \param ... The error message in printf form, plus any arguments to be inserted into the printf string. I'll provide the function name and a carriage return.
 */
#define Apop_assert_void(test,  level, stop, ...) do \
    if (!(test)) {  \
        if (apop_opts.verbose >= level) { fprintf(stderr, "%s: ", __func__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");}   \
        if (stop == 's' || stop == 'h') assert(test);   \
} while (0);

#define apop_assert_void(test, level, stop, ...) Apop_assert_void(test, level, stop, __VA_ARGS__)
#define APOP_ASSERT_VOID(test, level, stop, ...) Apop_assert_void(test, level, stop, __VA_ARGS__)

/** An unfortunate side-effect of using standardized naming is that you
 often get unpleasantly verbose declarations like
 \code
 apop_mle_settings *ms = apop_mle_settings_alloc(data, model);
 \endcode
 which is clear and simple but a darn lot of typing. I asked Roget's
 thesaurus for something shorter than ``settings'' that would get the
 idea across, but the closest it could come was ``aura''.

 However, because the names are so standard, it's easy to write a macro
 \code
 Apop_settings_alloc(mle, ms, data, model);
 \endcode
 that expands to the above. As of this writing, options for the first argument include \ref apop_mle_settings_alloc "mle", \ref apop_histogram_settings_alloc "histogram", \ref apop_ls_settings_alloc "ls", and \ref apop_update_settings_alloc "update". See the respective documentations for the arguments to be sent to the respective allocation functions.

ps: we just capitalize the first letter to remind you that it's a macro, but so that it doesn't look like we're yelling. I mean, as long as this is a syntactic-sugar macro for aesthetic purposes, it might as well look good.

 */
#define Apop_settings_alloc(type, out, ...) apop_ ##type ##_settings *out = apop_ ##type ##_settings_alloc(__VA_ARGS__);

#define APOP_SETTINGS_ALLOC(type, out, ...) Apop_settings_alloc(type, out, __VA_ARGS__)

////////Bootstrapping & RNG
apop_data * apop_jackknife_cov(apop_data *data, apop_model model);
APOP_VAR_DECLARE apop_data * apop_bootstrap_cov(apop_data *data, apop_model model, gsl_rng* rng, int iterations);
gsl_rng *apop_rng_alloc(int seed);

////////Missing data
apop_data * apop_data_listwise_delete(apop_data *d);
apop_model * apop_ml_imputation(apop_data *d, apop_model* meanvar);

/** Method settings for a model to be put through Bayesian updating. 
\li starting_pt      The first parameter to check in the MCMC routine
\li periods How many steps should the MCMC chain run?
\li burnin  What <em>percentage</em> of the periods should be ignored as initialization. That is, this is a number between zero and one.
\li histosegments If outputting a \ref apop_histogram, how many segments should it have?
 
 */
typedef struct{
    apop_data *data;
    apop_data *starting_pt;
    long int periods;
    double burnin;
    int histosegments;
    char method;
} apop_update_settings;

apop_update_settings *apop_update_settings_alloc(apop_data *d);
apop_update_settings *apop_update_settings_init(apop_update_settings);
#define apop_update_settings_copy NULL
#define  apop_update_settings_free NULL

APOP_VAR_DECLARE apop_model * apop_update(apop_data *data, apop_model *prior, apop_model *likelihood, gsl_rng *rng);

APOP_VAR_DECLARE double apop_test(double statistic, char *distribution, double p1, double p2, char tail);

__END_DECLS
#endif
