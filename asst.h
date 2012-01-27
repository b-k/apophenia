/** \file asst.h  */
/* Copyright (c) 2006--2007, 2010 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
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

double apop_generalized_harmonic(int N, double s) __attribute__ ((__pure__));

apop_data * apop_test_anova_independence(apop_data *d);
#define apop_test_ANOVA_independence(d) apop_test_anova_independence(d)

int apop_strcmp(char const *, char const *);
APOP_VAR_DECLARE int  apop_regex(const char *string, const char* regex, apop_data **substrings, const char use_case);

int apop_system(const char *fmt, ...) __attribute__ ((format (printf,1,2)));

//Histograms and PMFs
apop_model *apop_crosstab_to_pmf(apop_data *d);
gsl_vector * apop_vector_moving_average(gsl_vector *, size_t);
apop_model *apop_histogram_moving_average(apop_model *m, size_t bandwidth);
apop_model * apop_histogram_vector_reset(apop_model *, gsl_vector *);
APOP_VAR_DECLARE apop_model * apop_histogram_model_reset(apop_model *base, apop_model *m, long int draws, gsl_rng *rng);
apop_data * apop_histograms_test_goodness_of_fit(apop_model *h0, apop_model *h1);
apop_data * apop_test_kolmogorov(apop_model *m1, apop_model *m2);
void apop_histogram_normalize(apop_model *m);
apop_data *apop_data_pmf_compress(apop_data *in);
APOP_VAR_DECLARE apop_data * apop_data_to_bins(apop_data *indata, apop_data *binspec, int bin_count, char close_top_bin);
APOP_VAR_DECLARE apop_model * apop_model_to_pmf(apop_model *model, apop_data *binspec, long int draws, int bin_count, gsl_rng *rng);

//text conveniences
char * apop_strip_dots(char const *in, char strip_type);
APOP_VAR_DECLARE char* apop_text_paste(apop_data *strings, char *between, char *before, char *after, char *between_cols, int (*prune)(apop_data* ! int ! int ! void*), void* prune_parameter);
/** Notify the user of errors, warning, or debug info. 

 \param verbosity   At what verbosity level should the user be warned? E.g., if level==2, then print iff apop_opts.verbosity >= 2.
 \param ... The message to write to STDERR (presuming the verbosity level is high enough). This can be a printf-style format with following arguments. You can produce much more informative error messages this way, e.g., \c apop_notify(0, "Beta is %g but should be greater than zero.", beta);.
*/
#define Apop_notify(verbosity, ...) {\
    if (apop_opts.verbose >= verbosity) {  \
        fprintf(stderr, "%s: ", __func__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");   \
        fflush(stderr); \
} }

/** Tests whether a condition is true, and if it is not, prints an error to \c stderr and 
  exits the function. 
 
\param test The expression that you are asserting is nonzero.
\param returnval If the assertion fails, return this. If your assertion is inside a function that returns nothing, then leave this blank, as in <tt>apop_assert_c(a==b, , 0, "a should equal b");</tt>
\param level Print the warning message only if \ref apop_opts_type "apop_opts.verbose" is greater than or equal to this. Zero usually works, but for minor infractions use one.
\param ... The error message in printf form, plus any arguments to be inserted into the printf string. I'll provide the function name and a carriage return.

\li If \ref apop_opts.stop_on_warning is nonzero and not <tt>'v'</tt>, then a failed test halts via \c abort(), even if the <tt>apop_opts.verbose</tt> level is set so that the warning message doesn't print to screen. Use this when running via debugger.
\li If \ref apop_opts.stop_on_warning is <tt>'v'</tt>, then a failed test halts via \c abort() iff the verbosity level is high enough to print the error.

\see \ref Apop_assert, which always halts on error.
 */
#define Apop_assert_c(test, returnval, level, ...) {\
     if (!(test)) {  \
        Apop_notify(level,  __VA_ARGS__);   \
        if ((apop_opts.verbose >= level && apop_opts.stop_on_warning == 'v') \
            || (apop_opts.stop_on_warning && apop_opts.stop_on_warning !='v') ) \
                abort(); \
        return returnval;  \
    } }

/** This is just a slightly more user-friendly version of the C-standard \c assert(). The
  program halts if the first argument evaluates to false, and the remaining arguments are
  a printf-style message to display to \c stderr in such an event.

  This is how Apophenia does (almost) all of its assertions, and is made public as a
  convenience to you. You can see that it isn't hard to re-implement.

\param test An expression that, if false, halts the program
\param ... A printf-style message to display on \c stderr on halt. I'll provide the function name and a carriage return.

\see \ref Apop_assert_c, which continues with a message rather than shutting down.

  */
#define Apop_assert(test, ...)  {\
    if (!(test)) {  \
        fprintf(stderr, "%s: ", __func__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");   \
        abort();   \
} }

#define apop_assert_s Apop_assert
#define apop_assert Apop_assert
#define Apop_assert_s Apop_assert
#define apop_assert_c Apop_assert_c

//Missing data
APOP_VAR_DECLARE apop_data * apop_data_listwise_delete(apop_data *d, char inplace);
apop_model * apop_ml_impute(apop_data *d, apop_model* meanvar);
#define apop_ml_imputation(d, m) apop_ml_impute(d, m)
APOP_VAR_DECLARE apop_data * apop_multiple_imputation_variance(apop_data *base_data, apop_data *fill_ins, apop_data *(*stat)(apop_data *), char * row_name, char * col_name, char * value_name, char * imputation_name);


APOP_VAR_DECLARE apop_model * apop_update(apop_data *data, apop_model *prior, apop_model *likelihood, gsl_rng *rng);

APOP_VAR_DECLARE double apop_test(double statistic, char *distribution, double p1, double p2, char tail);

//Sorting (apop_asst.c)
APOP_VAR_DECLARE double * apop_vector_percentiles(gsl_vector *data, char rounding); 
APOP_VAR_DECLARE apop_data * apop_data_sort(apop_data *data, int sortby, char asc);

//raking
APOP_VAR_DECLARE apop_data * apop_rake(char *table_name, char *all_vars, char **contrasts, int contrast_ct, char *structural_zeros, int max_iterations, double tolerance, char *count_col, int run_number, char *init_table, char *init_count_col, double nudge);

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
