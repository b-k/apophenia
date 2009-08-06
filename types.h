//types.h   Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.
#ifndef __apop_estimate__
#define __apop_estimate__

#include <assert.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include "variadic.h"

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

/**\defgroup types Types defined by Apophenia. 

The basic story for a statistical analysis is that the researcher
assembles a data set into an \ref apop_data structure, then sends it to
an \ref apop_model so that the model's parameters can be estimated,
and that is returned in an \ref apop_model structure.

Supporting these main structures are a few more structures you'd only have to worry about 
for fine tuning.
The \ref apop_name
structure is an organized list of row and column names; functions that
take an \ref apop_data set try to automatically handle the names for you.
The more elaborate models, such as MLEs, require some parameters to run,
in which case you will need to fill out an \ref apop_model form and hand it in to the model.


\li Data 
The \ref apop_data structure adds a touch of metadata on
top of the basic \c gsl_matrix and \c gsl_vector. It includes an
\ref apop_name structure (see below), and a table for non-numeric
variables. See \ref data_struct.

\li Models 
The \ref _apop_model "apop_model" structure encapsulates a description of the world
in which the data and the parameters produce observed outcomes. The
\ref apop_estimate function takes in data and an un-parametrizes model and outputs a parametrized model.
See \ref models, or the full declaration of the structure on the \ref _apop_model "apop_model" page.

\li Names 
The \ref apop_name structure has three components: a list of column
names, a list of row names, and a list of dependent variable names. It
is intended to accompany the <tt>gsl_matrix</tt> structure, which holds
all the other information about a data aaray such as the number of rows
and columns.  See \ref names.

*/


/** A data set is assumed to be a matrix where each row is a single
observation and each column is a variable. Usually there is only one
dependent variable (the value to be predicted), which is the first column;
the independent variables (the predictors) follow thereafter.

This structure holds the names of these variables. It is pretty much always associated with (and generated along with) an \c apop_data set.
\ingroup names
*/
typedef struct{
	char * vector;
	char ** column;
	char ** row;
	char ** text;
	int colct, rowct, textct;
    char title[101];
} apop_name;

/**
Gathers together a <tt>gsl_vector</tt>, a <tt>gsl_matrix</tt>, an \ref apop_name structure, and a space for a table of non-numeric data.
Allocate using \c apop_data_alloc, free via \c apop_data_free, or more generally, see the \c apop_data_... section of the index (in the header links) for the many other functions that operate on this struct.
\ingroup data_struct
*/
typedef struct {
    gsl_vector  *vector;
    gsl_matrix  *matrix;
    apop_name   *names;
    char        ***text;
    int         textsize[2];
    gsl_vector  *weights;
} apop_data;

/** A description of a parametrized statistical model, including the
input settings and the output parameters, expected values, et cetera.
The full declaration is given in the \c _apop_model page, see the longer discussion on the \ref models page, or see 
the \ref apop_ols page for a sample program that uses an \ref apop_model.
\ingroup types
*/
typedef struct _apop_model apop_model;


typedef struct {
    char name[101];
    void *setting_group;
    void *copy;
    void *free;
} apop_settings_type;


/**
\param parameters 	The vector of coefficients or parameters estimated by the regression/MLE. Usually has as many dimensions as your data set has columns.
\param expected	An \ref apop_data structure with
three columns. If this is a model with a single dependent and lots of
independent vars, then the first column is the actual data. Let our model be \f$ Y = \beta X + \epsilon\f$. Then the second column is the predicted values: \f$\beta X\f$, and the third column is the residuals: \f$\epsilon\f$. The third column is therefore always the first minus the second, and this is probably how that column was calculated internally. There is thus currently no way to get just the predicted but not the residuals or vice versa.
\param covariance 	The variance-covariance matrix.
\param status		The return status from the estimate that had populated this apop_model, if any.
\ingroup inv_and_est
*/
struct _apop_model{
    char        name[101]; 
    int         vbase, m1base, m2base;
    apop_settings_type *settings;
    apop_data   *parameters, *expected, *covariance;
    double      llikelihood;
    int         prepared, status;
    apop_data   *data;
    apop_model * (*estimate)(apop_data * data, apop_model *params);
    double  (*p)(apop_data *d, apop_model *params);
    double  (*log_likelihood)(apop_data *d, apop_model *params);
    void    (*score)(apop_data *d, gsl_vector *gradient, apop_model *params);
    double  (*constraint)(apop_data *data, apop_model *params);
    apop_data*  (*expected_value)(apop_data *d, apop_model *params);
    void (*draw)(double *out, gsl_rng* r, apop_model *params);
    void (*prep)(apop_data *data, apop_model *params);
    void (*print)(apop_model *params);
    void    *more;
    size_t  more_size;
} ;

/** The global options.
  \ingroup global_vars */
typedef struct{
            /** Set this to zero for silent mode, one for errors and warnings. default = 0. */
    int verbose;
            /** 's'   = to screen
                'f'   = to file
                'd'   = to db. 
                'p'   = to pipe (specifically, apop_opts.output_pipe). 
             If 1 or 2, then you'll need to set output_name in the apop_..._print fn. default = 0. */
    char output_type;
            /** If printing to a pipe or FILE, set it here. */
    FILE *output_pipe;
            /** The separator between elements of output tables. The
             default is "\t", but for LaTeX, use "&\t", or use "|" to
             get pipe-delimited output. */
    char output_delimiter[100];
            /** Append to output files(1), or overwrite(0)? default = 0 */
    int output_append;
            /** What other people have put between your columns. Default = "|,\t" */
    char input_delimiters[100];
            /** If set, the name of the column in your tables that holds row names. */
    char db_name_column[300];
            /** The string that the database takes to indicate NaN. May be a regex. */
    char db_nan[100];
            /** If this is 'm', use mySQL, else use SQLite. */
    char db_engine;
            /** Username for database login. Max 100 chars. */
    char db_user[101];
            /** Password for database login. Max 100 chars. */
    char db_pass[101];
            /** Threads to use internally. See \ref apop_matrix_apply and family. */
    int  thread_count;
    int  rng_seed;
    float version;
} apop_opts_type;

extern apop_opts_type apop_opts;

apop_name * apop_name_alloc(void);
int apop_name_add(apop_name * n, char *add_me, char type);
void  apop_name_free(apop_name * free_me);
void  apop_name_print(apop_name * n);
APOP_VAR_DECLARE void  apop_name_stack(apop_name * n1, apop_name *nadd, char type1, char typeadd);
void  apop_name_cross_stack(apop_name * n1, apop_name *n2, char type1, char type2);
apop_name * apop_name_copy(apop_name *in);
int  apop_name_find(apop_name *n, char *findme, char type);

void 		apop_model_free (apop_model * free_me);
void 		apop_model_show (apop_model * print_me);

void        apop_data_free(apop_data *freeme);
apop_data * apop_matrix_to_data(gsl_matrix *m);
apop_data * apop_vector_to_data(gsl_vector *v);
apop_data * apop_data_alloc(const size_t, const size_t, const int);
apop_data * apop_data_calloc(const size_t, const size_t, const int);
APOP_VAR_DECLARE apop_data * apop_data_stack(apop_data *m1, apop_data * m2, char posn, char inplace);
apop_data ** apop_data_split(apop_data *in, int splitpoint, char r_or_c);
apop_data * apop_data_copy(const apop_data *in);
void        apop_data_rm_columns(apop_data *d, int *drop);
void apop_data_memcpy(apop_data *out, const apop_data *in);
double * apop_data_ptr(const apop_data *data, const int i, const int j);
double * apop_data_ptr_it(const apop_data *in, size_t row, char* col);
double * apop_data_ptr_ti(const apop_data *in, char* row, int col);
double * apop_data_ptr_tt(const apop_data *in, char *row, char* col);
double apop_data_get(const apop_data *in, size_t row, int  col);
double apop_data_get_it(const apop_data *in, size_t row, char* col);
double apop_data_get_ti(const apop_data *in, char* row, int col);
double apop_data_get_tt(const apop_data *in, char *row, char* col);
void apop_data_set(apop_data *in, size_t row, int col, double data);
void apop_data_set_ti(apop_data *in, char* row, int col, double data);
void apop_data_set_it(apop_data *in, size_t row, char* col, double data);
void apop_data_set_tt(apop_data *in, char *row, char* col, double data);
void apop_data_add_named_elmt(apop_data *d, char *name, double val);
void apop_text_add(apop_data *in, const size_t row, const size_t col, const char *fmt, ...);
apop_data * apop_text_alloc(apop_data *in, const size_t row, const size_t col);
apop_data *apop_data_transpose(apop_data *in);
gsl_matrix * apop_matrix_realloc(gsl_matrix *m, size_t newheight, size_t newwidth);
gsl_vector * apop_vector_realloc(gsl_vector *v, size_t newheight);

void apop_text_free(char ***freeme, int rows, int cols); //in apop_data.c

void apop_opts_memcpy(apop_opts_type *out, apop_opts_type *in); //in apop_output.c

__END_DECLS
#endif
