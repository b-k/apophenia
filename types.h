//types.h			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
#ifndef __apop_estimate__
#define __apop_estimate__

#include <assert.h>
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

/**\defgroup types Types defined by Apophenia. 

The basic story for a statistical analysis is that the researcher
assembles a data set into an \ref apop_data structure, then sends it to
an \ref apop_model so that the model's parameters can be estimated,
and that is returned in an \ref apop_params structure.

Supporting these main structures are a few more structures you'd only have to worry about 
for fine tuning.
The \ref apop_name
structure is an organized list of row and column names; functions that
take an \ref apop_data set try to automatically handle the names for you.
The more elaborate models, such as MLEs, require some parameters to run,
in which case you will need to fill out an \ref apop_params form and hand it in to the model.


\li Data 
The \ref apop_data structure adds a touch of metadata on
top of the basic \c gsl_matrix and \c gsl_vector. It includes an
\ref apop_name structure (see below), and a table for non-numeric
variables. See \ref data_struct.

\li Models 
The \ref apop_model structure encapsulates a description of the world
in which the data and the parameters produce observed outcomes. The
apop_model.estimate() method takes in data and produces an \ref
apop_estimate. See \ref models.

\li Estimates
The \ref apop_params structure complements the \c apop_model by
providing input params (like the method and tolerance for an ML
estimation) and the output parameters.

\li Names 
The \ref apop_name structure has three components: a list of column
names, a list of row names, and a list of dependent variable names. It
is intended to accompany the <tt>gsl_matrix</tt> structure, which holds
all the other information about a data aaray such as the number of rows
and columns.  See \ref names.

*/

/** The structure has two uses. The first is to tell the regression/MLE
functions what you would like to receive in return. Alternatively, you
can just send in a <tt>NULL</tt> pointer, and the functions will return
everything apropos.

The second is for the internal workings of the \ref apop_estimate
structure, giving a list of the elements of the structure which are
actually in use. For example, the regressions won't return a log
likelihood, and the ML estimates won't return an  R^2.


\b the elements 
\verbatim
int     parameters, covariance, confidence, dependent, predicted, log_likelihood;
\endverbatim

There is one element for each element of the \ref apop_params structure.

The \ref apop_params structure
has an element named <tt>uses</tt> embedded
within it. Those elements for which <tt>uses.elmt</tt> are zero are
unallocated pointers (so be careful: precede all dereferences with an
<tt>if(est->uses.element)</tt> clause).


It may sometimes be useful to manipulate the ["apop_estimate"] structure's
<tt>uses</tt> element to your own benefit. For
example, if you set <tt>est->ep.uses.predicted =
0</tt> before calling <tt>apop_print_estimate(est, NULL)</tt>, then
the predicted values won't get printed. But be careful: if you then call
<tt>apop_estimate_free(est)</tt>, then the predicted values won't get freed,
either.  
*/

#include <string.h>


/** A data set is assumed to be a matrix where each row is a single
observation and each column is a variable. Usually there is only one
dependent variable (the value to be predicted), which is the first column;
the independent variables (the predictors) follow thereafter.

This structure holds the names of these variables. You can fill it quickly
with \ref apop_db_get_names after running a query, or add names manually
with \ref apop_name_add .

Typically, the row names are not used, but they are there for your convenience.  
\ingroup names
*/
typedef struct{
	char * vecname;
	char ** colnames;
	char ** rownames;
	char ** catnames;
	char ** textnames;
	int colnamect, rownamect, catnamect, textnamect;
} apop_name;

/**
Gathers together a <tt>gsl_vector</tt>, a <tt>gsl_matrix</tt>, an \ref apop_name structure, and a space for a table of non-numeric data.
\ingroup data_struct
*/
typedef struct {
    gsl_vector  *vector;
    gsl_matrix  *matrix;
    apop_name   *names;
    char        ***categories;
    char        ***text;
    int         catsize[2];
    int         textsize[2];
    gsl_vector  *weights;
} apop_data;

/** The data to accompany an \c apop_model, including the input settings and the output parameters, expected values, et cetera.

<b>An example</b><br>

The \ref apop_OLS page has a sample program which uses an <tt>apop_estimate</tt> structure.

\param parameters 	The vector of coefficients or parameters estimated by the regression/MLE. Usually has as many dimensions as your data set has columns.
\param dependent	An \ref apop_data structure with
three columns. If this is a model with a single dependent and lots of
independent vars, then the first column is the actual data. Let our model be \f$ Y = \beta X + \epsilon\f$. Then the second column is the predicted values: \f$\beta X\f$, and the third column is the residuals: \f$\epsilon\f$. The third column is therefore always the first minus the second, and this is probably how that column was calculated internally. There is thus currently no way to get just the predicted but not the residuals or vice versa.
\param covariance 	The variance-covariance matrix.
\param confidence 	The two-tailed test of the hypothesis that the variable is zero. One element for each parameter.
\param status		The return status from the estimate that had populated this apop_estimate, if any.
\ingroup inv_and_est
*/
typedef struct{
    char        method_name[101];
    void        *method_params;
    void        *model_params;
    struct {
        char    parameters, covariance, confidence, expected, predicted, log_likelihood;
    } uses;
    apop_data   *parameters, *expected, *covariance;
    double      log_likelihood;
    int         status;
    apop_data   *data;
    struct apop_model  *model;
    void        *more;
} apop_params;

/** This is an object to describe a model whose parameters are to be
estimated. It would primarily be used for maximum likelihood estimation,
but is intended to have anything else you would want a probability
distribution to have too, like a random number generator.  

\param name	The model name. You have 100 characters. 
\param parameter_ct	The number of parameters. If this is 0, it will be dynamically set to the number of columns in the given data set's matrix; if -1 it will be set to columns minus one.
\param estimate		the estimator fn, which is all most users will care about.
\param log_likelihood	the likelihood fn given data 
\param 	dlog_likelihood	the derivative of the likelihood fn
\param 	fdf	Do both of the above at once. Can be NULL if it'd just call them separately. 
\param 	constraint	The constraints to the parameters, if any. Really only necessary for MLEs.
\param rng 	a random number generator. 

\ingroup models
 */
typedef struct apop_model{
    char    name[101]; 
    int     vbase, m1base, m2base;
    apop_params * (*estimate)(apop_data * data, apop_params *params);
    double  (*p)(const apop_data *beta, apop_data *d, apop_params *params);
    double  (*log_likelihood)(const apop_data *beta, apop_data *d, apop_params *params);
    void    (*score)(const apop_data *beta, apop_data *d, gsl_vector *gradient, apop_params *params);
    double  (*constraint)(const apop_data *beta, apop_data *returned_beta, apop_params *params);
    void (*draw)(double *out, apop_data *beta, gsl_rng* r, apop_params *params);
    void    *more;
    apop_params *ep;
} apop_model;

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
            /** If writing to a file, its name. Limit: 1000 chars. */
    char output_delimiter[100];
            /** Append to output files(1), or overwrite(0)? default = 0 */
    int output_append;
            /** What other people have put between your columns. Default = "|,\t" */
    char input_delimiters[100];
            /** If set, the name of the column in your tables that holds row names. */
    char db_name_column[300];
            /** If set, plot the path of the max. likelihood search. */
    char db_nan[100];
            /** If this is 'm', use mySQL, else use SQLite. */
    char db_engine;
            /** If set, plot the path of the max. likelihood search. */
    char  mle_trace_path[1000];
            /** Threads to use internally. See \ref apop_apply. */
    int  thread_count;
} apop_opts_type;

extern apop_opts_type apop_opts;

apop_name * apop_name_alloc(void);
int apop_name_add(apop_name * n, char *add_me, char type);
void  apop_name_free(apop_name * free_me);
void  apop_name_print(apop_name * n);
void  apop_name_stack(apop_name * n1, apop_name *n2, char type);
void  apop_name_cross_stack(apop_name * n1, apop_name *n2, char type1, char type2);
void apop_name_rm_columns(apop_name *n, int *drop);
void apop_name_memcpy(apop_name **out, apop_name *in);
apop_name * apop_name_copy(apop_name *in);
size_t  apop_name_find(apop_name *n, char *findme, char type);

apop_params * apop_params_alloc(apop_data * data, apop_model *model, void *method_params, void *model_params);
apop_params *apop_params_copy(apop_params *in);
apop_params *apop_params_clone(apop_params *in, size_t method_size, size_t model_size, size_t more_size);
void 		apop_params_free (apop_params * free_me);
void 		apop_params_print (apop_params * print_me);
void 		apop_params_show (apop_params * print_me);

void        apop_data_free(apop_data *freeme);
apop_data * apop_matrix_to_data(gsl_matrix *m);
apop_data * apop_data_from_matrix(gsl_matrix *m);
apop_data * apop_vector_to_data(gsl_vector *v);
apop_data * apop_data_from_vector(gsl_vector *v);
apop_data * apop_data_alloc(const size_t, const size_t, const int);
apop_data * apop_data_calloc(const size_t, const size_t, const int);
apop_data * apop_data_stack(apop_data *m1, apop_data * m2, char posn);
apop_data ** apop_data_split(apop_data *in, int splitpoint, char r_or_c);
apop_data * apop_data_copy(const apop_data *in);
void        apop_data_rm_columns(apop_data *d, int *drop);
void apop_data_memcpy(apop_data *out, const apop_data *in);
double * apop_data_ptr(const apop_data *data, const size_t i, const size_t j);
double apop_data_get(const apop_data *in, size_t row, int  col);
double apop_data_get_nt(const apop_data *in, size_t row, char* col);
double apop_data_get_tn(const apop_data *in, char* row, int col);
double apop_data_get_tt(const apop_data *in, char *row, char* col);
void apop_data_set(apop_data *in, size_t row, int col, double data);
void apop_data_set_tn(apop_data *in, char* row, int col, double data);
void apop_data_set_nt(apop_data *in, size_t row, char* col, double data);
void apop_data_set_tt(apop_data *in, char *row, char* col, double data);
void apop_data_add_named_elmt(apop_data *d, char *name, double val);

void apop_text_free(char ***freeme, int rows, int cols); //in apop_data.c

apop_model * apop_model_copy(apop_model in); //in apop_estimate.c.

void apop_opts_memcpy(apop_opts_type *out, apop_opts_type *in); //in apop_output.c

__END_DECLS
#endif
