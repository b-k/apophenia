/* These are functions used here and there to write Apophenia. They're
 not incredibly useful, or even very good form, so they're not public. Cut
 & paste `em into your own code if you'd like.
 */

/* Many Apop functions try to treat the vector and matrix equally, which
 requires knowing which exists and what the sizes are. */
#define Get_vmsizes(d) \
    int firstcol = d && (d)->vector ? -1 : 0; \
    int vsize = d && (d)->vector ? (d)->vector->size : 0; \
    int wsize = d && (d)->weights ? (d)->weights->size : 0; \
    int msize1 = d && (d)->matrix ? (d)->matrix->size1 : 0; \
    int msize2 = d && (d)->matrix ? (d)->matrix->size2 : 0; \
    int tsize = vsize + msize1*msize2; \
    int maxsize = GSL_MAX(vsize, GSL_MAX(msize1, d?d->textsize[0]:0));\
    (void)(tsize||wsize||firstcol||maxsize) /*prevent unused variable complaints */;

// Define a static variable, and initialize on first use.
#define Staticdef(type, name, def) static type (name) = NULL; if (!(name)) (name) = (def);

// Check for NULL and complain if so.
#define Nullcheck(in, errval) Apop_assert_c(in, errval, apop_errorlevel, "%s is NULL.", #in);
#define Nullcheck_m(in, errval) Apop_assert_c(in, errval, apop_errorlevel, "%s is a NULL model.", #in);
#define Nullcheck_mp(in, errval) Nullcheck_m(in, errval); Apop_assert_c((in)->parameters, errval, apop_errorlevel, "%s is a model with NULL parameters. Please set the parameters and try again.", #in);
#define Nullcheck_d(in, errval) Apop_assert_c(in, errval, apop_errorlevel, "%s is a NULL data set.", #in);
//And because I do them all so often:
#define Nullcheck_mpd(data, model, errval) Nullcheck_m(model, errval); Nullcheck_p(model, errval); Nullcheck_d(data, errval);
//deprecated:
#define Nullcheck_p(in, errval) Nullcheck_mp(in, errval);

//in apop_conversions.c Extend a string.
void xprintf(char **q, char *format, ...);
#define XN(in) ((in) ? (in) : "")

//For a pedantic compiler. Continues on error, because there's not much else to do: the computer is clearly broken.
#define Asprintf(...) Apop_stopif(asprintf(__VA_ARGS__)==-1, , 0, "Error printing to a string.")

#include <sqlite3.h>
#include <stddef.h>
int apop_use_sqlite_prepared_statements(size_t col_ct);
int apop_prepare_prepared_statements(char const *tabname, size_t col_ct, sqlite3_stmt **statement);
char *prep_string_for_sqlite(int prepped_statements, char const *astring);//apop_conversions.c
void apop_gsl_error(char const *reason, char const *file, int line, int gsl_errno); //apop_linear_algebra.c

//For when we're forced to use a global variable.
#undef threadlocal
#if __STDC_VERSION__ > 201100L
    #define threadlocal _Thread_local
#elif defined(__APPLE__) 
    #define threadlocal
#elif defined(__GNUC__) && !defined(threadlocal)
    #define threadlocal __thread
#else
    #define threadlocal
#endif

#ifdef _OPENMP
#define PRAGMA(x) _Pragma(#x)
#define OMP_critical(tag) PRAGMA(omp critical ( tag ))
#define OMP_for(...) _Pragma("omp parallel for") for(__VA_ARGS__)
#else
#define OMP_critical(tag)
#define OMP_for(...) for(__VA_ARGS__)
#endif

#include "config.h"
#ifndef HAVE___ATTRIBUTE__
#define __attribute__(...)
#endif

#ifndef HAVE_ASPRINTF
#include <stdarg.h>

//asprintf, vararg, &c
extern int asprintf (char **res, const char *format, ...)
       __attribute__ ((__format__ (__printf__, 2, 3)));
extern int vasprintf (char **res, const char *format, va_list args)
       __attribute__ ((__format__ (__printf__, 2, 0)));
#endif

#include "apop.h"
void add_info_criteria(apop_data *d, apop_model *m, apop_model *est, double ll); //In apop_mle.c
