/** \file asst.h  */
/* Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#ifndef __apop_asst__
#define __apop_asst__

#include <assert.h>
#include "types.h"
#include "variadic.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

#ifdef	__cplusplus
extern "C" {
#endif

double apop_generalized_harmonic(int N, double s);
void apop_error(int level, char stop, char *message, ...);

apop_data * apop_test_anova_independence(apop_data *d);
#define apop_test_ANOVA_independence(d) apop_test_anova_independence(d)

int apop_strcmp(char *, char*);
APOP_VAR_DECLARE int  apop_regex(const char *string, const char* regex, apop_data **substrings, const char use_case);

gsl_vector * apop_vector_moving_average(gsl_vector *, size_t);
apop_model *apop_histogram_moving_average(apop_model *m, size_t bandwidth);

apop_model * apop_histogram_vector_reset(apop_model *, gsl_vector *);
APOP_VAR_DECLARE apop_model * apop_histogram_model_reset(apop_model *base, apop_model *m, long int draws, gsl_rng *rng);
apop_data * apop_histograms_test_goodness_of_fit(apop_model *h0, apop_model *h1);
apop_data * apop_test_kolmogorov(apop_model *m1, apop_model *m2);
void apop_histogram_normalize(apop_model *m);

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


/** A simplified version of \ref Apop_assert that always stops. Therefore, you don't need
  the \c returnval or the stop/continue options. Because a stop is invariably a true
  error, the verbosity level is also omitted, so just give the test and the error.*/
#define Apop_assert_s(test, ...) do \
    if (!(test)) {  \
        fprintf(stderr, "%s: ", __func__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");   \
        assert(test);   \
} while (0);

#define apop_assert(test, returnval, level, stop, ...) Apop_assert(test, returnval, level, stop, __VA_ARGS__)
#define apop_assert_s Apop_assert_s
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

/** \deprecated Use \ref Apop_model_add_group.
 
  For what it's worth, this is a convenience macro. Expands:
 \code
 Apop_settings_alloc(mle, ms, data, model);
 \endcode
to:
 \code
 apop_mle_settings *ms = apop_mle_settings_alloc(data, model);
 \endcode
 As of this writing, options for the first argument include \ref apop_mle_settings_init "mle", \ref apop_histogram_settings_init "histogram", and \ref apop_update_settings_init "update". See the respective documentations for the arguments to be sent to the respective allocation functions. Because this is an obsolete function, that list may shrink.

 */
#define Apop_settings_alloc(type, out, ...) apop_ ##type ##_settings *out = apop_ ##type ##_settings_alloc(__VA_ARGS__);

#define APOP_SETTINGS_ALLOC(type, out, ...) Apop_settings_alloc(type, out, __VA_ARGS__)

//Bootstrapping & RNG
apop_data * apop_jackknife_cov(apop_data *data, apop_model model);
APOP_VAR_DECLARE apop_data * apop_bootstrap_cov(apop_data *data, apop_model model, gsl_rng* rng, int iterations);
gsl_rng *apop_rng_alloc(int seed);

//Missing data
APOP_VAR_DECLARE apop_data * apop_data_listwise_delete(apop_data *d, char inplace);
apop_model * apop_ml_impute(apop_data *d, apop_model* meanvar);
#define apop_ml_imputation(d, m) apop_ml_impute(d, m)
apop_data * apop_multiple_imputation_variance(apop_data *(*stat)(apop_data *), apop_data *base_data, apop_data *fill_ins);


APOP_VAR_DECLARE apop_model * apop_update(apop_data *data, apop_model *prior, apop_model *likelihood, gsl_rng *rng);

APOP_VAR_DECLARE double apop_test(double statistic, char *distribution, double p1, double p2, char tail);


//PMF (model/apop_pmf.c)
apop_model *apop_crosstab_to_pmf(apop_data *d);

//Sorting (apop_asst.c)
APOP_VAR_DECLARE double * apop_vector_percentiles(gsl_vector *data, char rounding); 
APOP_VAR_DECLARE apop_data * apop_data_sort(apop_data *data, int sortby, char asc);

//asprintf, vararg, &c
#include <stdarg.h>
extern int asprintf (char **res, const char *format, ...)
       __attribute__ ((__format__ (__printf__, 2, 3)));
extern int vasprintf (char **res, const char *format, va_list args)
       __attribute__ ((__format__ (__printf__, 2, 0)));

#ifdef	__cplusplus
}
#endif

#endif
