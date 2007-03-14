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
and that is returned in an \ref apop_estimate structure.

Supporting these main structures are a few more structures you'd only have to worry about 
for fine tuning.
The \ref apop_name
structure is an organized list of row and column names; functions that
take an \ref apop_data set try to automatically handle the names for you.
The \ref apop_inventory structure lists the elements that an \ref
apop_estimate has generated (no model produces everything). Model
estimates accept an inventory, but all work fine if you send in
NULL. Finally, the more elaborate models, such as MLEs, require some parameters to run,
in which case you will need to fill out an \ref apop_ep form and hand it in to the model.


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
The \ref apop_estimate structure returns all the data one would want
from a regression or ML estimation, including the parameters estimated,
the variance/covariance matrix, the residuals, et cetera. The structure
includes instances of both of the strucutres below. See \ref inv_and_est.

\li Names 
The \ref apop_name structure has three components: a list of column
names, a list of row names, and a list of dependent variable names. It
is intended to accompany the <tt>gsl_matrix</tt> structure, which holds
all the other information about a data aaray such as the number of rows
and columns.  See \ref names.

\li Inventory
The \ref apop_inventory structure serves two purposes. It is an input
to a regression or ML estimation, tells the function what output you
would like the <tt>apop_estimate</tt> output to include. It is also an
output from these functions, since the returned <tt>apop_estimate</tt>
will include its own <tt>apop_inventory</tt>,  which can be used later on
to test whether any given element is in use. See \ref inv_and_est.

\li Estimate parameters
The \ref apop_ep are the details for how an \ref
apop_estimate should do its work; currently it is just the specifications
for tolerances, step sizes, starting points, et cetera, for \ref apop_maximum_likelihood.
 See \ref inv_and_est.

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

There is one element for each element of the \ref apop_estimate structure.

If the <tt>apop_inventory</tt> will be sent in to a regression/MLE
function, set the appropriate element to either zero or one if you would
like the function to return the designated \ref apop_estimate element.

The \ref apop_ep of the \ref apop_estimate structure
has an <tt>apop_inventory</tt> element named <tt>uses</tt> embedded
within it. Those elements for which <tt>uses.elmt</tt> are zero are
unallocated pointers (so be careful: precede all dereferences with an
<tt>if(est->uses.element)</tt> clause).

<b>functions</b><br>
\code
void apop_inventory_copy(apop_inventory in, apop_inventory *out);
\endcode
Copy the input inventory list to a new output list. Notice that the input list is an actual inventory, while the output is a pointer to an inventory (since it will be modified).

\code
void apop_inventory_set(apop_inventory *out, int value);
\endcode
Set all of the elements of the inventory to the value given, e.g.,
<tt>apop_set_inventory(&want_all, 1)</tt>. Clearly, <tt>value</tt>
should either be zero or one.

<b>notes </b><br>

Unlike almost everything else in the GSL and Apophenia, it is
generally assumed that <tt>apop_inventory</tt>s are not pointers, but
are automatically allocated. Notably, this is true of the <tt>uses</tt>
element of the \ref apop_ep structure; therefore, to check
whether the variance-covariance matrix of an <tt>apop_estimate*</tt>
is present, for example, you would look at 
<tt>est->ep.uses.covariance</tt>.


It may sometimes be useful to manipulate the ["apop_estimate"] structure's
internal <tt>apop_inventory</tt> element to your own benefit. For
example, if you set <tt>est->ep.uses.predicted =
0</tt> before calling <tt>apop_print_estimate(est, NULL)</tt>, then
the predicted values won't get printed. But be careful: if you then call
<tt>apop_estimate_free(est)</tt>, then the predicted values won't get freed,
either.  
*/
typedef struct {
	int	parameters, covariance, confidence, dependent, predicted, log_likelihood, names;
} apop_inventory;

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

/** Parameters for running estimations. No estimation uses all of them.
  E.g., the MLE functions don't look at preserve_data but OLS and GLS do, while OLS and GLS ignore all the other params.
 \ingroup inv_and_est
 */
typedef struct{
	double      *starting_pt; 
	double 	    step_size; 
	double 	    tolerance; 
	double 	    resolution; 
	int 	    method;
	int 	    verbose;
	int 	    destroy_data;
    int         params_per_column;
	apop_inventory	uses;
    void        *parameters;
    gsl_vector  *weights;
    struct apop_model  *model;
    void        *more;
} apop_ep;

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

/** Regression and MLE functions return this structure, which includes
the various elements that one would want from a model estimate.

If you need control of the types of information these functions return,
see the \ref apop_inventory page. [If you don't, just send <tt>NULL</tt>
every time a function asks for an <tt>apop_inventory*</tt> structure.]

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
	apop_data 	*parameters, *dependent, *covariance;
	double		log_likelihood;
	int		    status;
    apop_data   *data;
    struct apop_model  *model;
    apop_ep     ep;
} apop_estimate;

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
	char	name[101]; 
	int	parameter_ct;
	apop_estimate *	(*estimate)(apop_data * data, void *params);
	double 	(*p)(const gsl_vector *beta, apop_data *d, void *params);
	double 	(*log_likelihood)(const gsl_vector *beta, apop_data *d, void *params);
	void 	(*score)(const gsl_vector *beta, apop_data *d, gsl_vector *gradient, void *params);
    double  (*constraint)(gsl_vector *beta, void * d, gsl_vector *returned_beta, void *params);
	double (*draw)(gsl_rng* r, gsl_vector *a, void *params);
    void    *more;
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

apop_estimate * apop_estimate_alloc(apop_data * data, apop_model model, apop_ep *params);
void 		apop_estimate_free(apop_estimate * free_me);
void 		apop_estimate_print(apop_estimate * print_me);
void 		apop_estimate_show(apop_estimate * print_me);

apop_ep *apop_ep_alloc();
void apop_ep_free(apop_ep *freeme);

apop_inventory * apop_inventory_alloc(int value);
void 		apop_inventory_copy(apop_inventory in, apop_inventory *out);
void 		apop_inventory_set(apop_inventory *out, int value);
apop_inventory apop_inventory_filter(apop_inventory *in, apop_inventory filter);

void        apop_data_free(apop_data *freeme);
apop_data * apop_matrix_to_data(gsl_matrix *m);
apop_data * apop_data_from_matrix(gsl_matrix *m);
apop_data * apop_vector_to_data(gsl_vector *v);
apop_data * apop_data_from_vector(gsl_vector *v);
apop_data * apop_data_alloc(int size1, int size2);
apop_data * apop_data_calloc(int size1, int size2);
apop_data * apop_data_stack(apop_data *m1, apop_data * m2, char posn);
apop_data ** apop_data_split(apop_data *in, int splitpoint, char r_or_c);
apop_data * apop_data_copy(apop_data *in);
void        apop_data_rm_columns(apop_data *d, int *drop);
void apop_data_memcpy(apop_data *out, apop_data *in);
double apop_data_get(apop_data *in, size_t row, int  col);
double apop_data_get_nt(apop_data *in, size_t row, char* col);
double apop_data_get_tn(apop_data *in, char* row, int col);
double apop_data_get_tt(apop_data *in, char *row, char* col);
void apop_data_set(apop_data *in, size_t row, int col, double data);
void apop_data_set_tn(apop_data *in, char* row, int col, double data);
void apop_data_set_nt(apop_data *in, size_t row, char* col, double data);
void apop_data_set_tt(apop_data *in, char *row, char* col, double data);
void apop_data_add_named_elmt(apop_data *d, char *name, double val);

void apop_text_free(char ***freeme, int rows, int cols); //in apop_data.c

apop_model * apop_model_copy(apop_model in); //this is in apop_estimate.c.

void apop_opts_memcpy(apop_opts_type *out, apop_opts_type *in); //in apop_output.c

__END_DECLS
#endif
