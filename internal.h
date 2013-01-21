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
#define Nullcheck(in, errval) Apop_assert_c(in, errval, apop_errorlevel, "%s is NULL.", #in)
#define Nullcheck_m(in, errval) Apop_assert_c(in, errval, apop_errorlevel, "%s is a NULL model.", #in)
#define Nullcheck_mp(in, errval) Nullcheck_m(in, errval); Apop_assert_c((in)->parameters, errval, apop_errorlevel, "%s is a model with NULL parameters. Please set the parameters and try again.", #in)
#define Nullcheck_d(in, errval) Apop_assert_c(in, errval, apop_errorlevel, "%s is a NULL data set.", #in)
//And because I do them all so often:
#define Nullcheck_mpd(data, model, errval) Nullcheck_m(model, errval); Nullcheck_p(model, errval); Nullcheck_d(data, errval);
//deprecated:
#define Nullcheck_p(in, errval) Nullcheck_mp(in, errval) 

//in apop_conversions.c Extend a string.
void xprintf(char **q, char *format, ...);
#define XN(in) ((in) ? (in) : "")

#include <sqlite3.h>
#include <stddef.h>
int apop_use_sqlite_prepared_statements(size_t col_ct);
int apop_prepare_prepared_statements(char const *tabname, size_t col_ct, sqlite3_stmt **statement);
char *prep_string_for_sqlite(int prepped_statements, char const *astring);//apop_conversions.c
void apop_gsl_error(char const *reason, char const *file, int line, int gsl_errno); //apop_linear_algebra.c

//For when we're forced to use a global variable.
#undef threadlocal
#ifdef _ISOC11_SOURCE 
    #define threadlocal _Thread_local
#elif defined(__APPLE__) 
    #define threadlocal
#elif defined(__GNUC__) && !defined(threadlocal)
    #define threadlocal __thread
#else
    #define threadlocal
#endif
