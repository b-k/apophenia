/** \file types.h */
/* Copyright (c) 2005--2007, 2010 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2. */
#ifndef __apop_estimate__
#define __apop_estimate__

#include <assert.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include "variadic.h"

#ifdef	__cplusplus
extern "C" {
#endif

/** This structure holds the names of the components of the \ref apop_data set. You may never have to worry about it directly, because most operations on \ref apop_data sets will take care of the names for you.
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

/** The \ref apop_data structure represents a data set. It primarily joins together a gsl_vector, a gsl_matrix, and a table of strings, then gives them all row and column names. It tries to be minimally intrusive, so you can use it everywhere you would use a \c gsl_matrix or a \c gsl_vector.

If you are viewing the HTML documentation, here is a diagram showing a sample data set with all of the elements in place. Together, they represet a data set where each row is an observation, which includes both numeric and text values, and where each row/column is named.

\htmlinclude apop_data_fig.html

Allocate using \c apop_data_alloc, free via \c apop_data_free, or more generally, see the \c apop_data_... section of the index (in the header links) for the many other functions that operate on this struct.

See also the Data Sets section of the outline page (also in the header links) for further notes on getting and manipulating the elements of an \ref apop_data set.
*/
typedef struct apop_data{
    gsl_vector  *vector;
    gsl_matrix  *matrix;
    apop_name   *names;
    char        ***text;
    size_t      textsize[2];
    gsl_vector  *weights;
    struct apop_data   *more;
    char        error;
} apop_data;

/* Settings groups. For internal use only; see apop_settings.c and 
   settings.h for related machinery. */
typedef struct {
    char name[101];
    unsigned long name_hash;
    void *setting_group;
    void *copy;
    void *free;
} apop_settings_type;

/** A statistical model. */
typedef struct apop_model apop_model;

/** The elements of the \ref apop_model type, representing a statistical model. */
struct apop_model{
    char        name[101]; 
    int         vbase, m1base, m2base, dsize; /**< The size of the parameter set.
                     If a dimension is -1, then use yourdata->matrix->size2. For
                    anything more complex, allocate the parameter set in the prep
                    method. \c dsize is for the canonical form, and is
                    the size of the data the RNG will return. */
    apop_settings_type *settings;
    apop_data   *parameters; /**< The coefficients or parameters estimated by the model. */
    apop_data   *data; /**< The input data. Typically a link to what you sent to \ref apop_estimate */
    apop_data   *info; /**< Several pages of assorted info, perhaps including the log likelihood, AIC, BIC,
                        covariance matrix, confidence intervals, expected score. See your
                        specific model's documentation for what it puts here.
                        */
    apop_model * (*estimate)(apop_data * data, apop_model *params); 
                /**< The estimation routine. Call via \ref apop_estimate */
    double  (*p)(apop_data *d, apop_model *params);
                /**< Probability of the given data and parameterized model. Call via \ref apop_p */
    double  (*log_likelihood)(apop_data *d, apop_model *params);
                /**< Log likelihood of the given data and parameterized model. Call via \ref apop_log_likelihood */
    void    (*score)(apop_data *d, gsl_vector *gradient, apop_model *params);
                /**< Derivative of the log likelihood. Call via \ref apop_score */
    apop_data*  (*predict)(apop_data *d, apop_model *params);
    apop_model * (*parameter_model)(apop_data *, apop_model *);
    double  (*cdf)(apop_data *d, apop_model *params); /**< Cumulative distribution function: 
                            the integral up to the single data point you provide.  Call via \ref apop_cdf */
    double  (*constraint)(apop_data *data, apop_model *params);
    void (*draw)(double *out, gsl_rng* r, apop_model *params);
                /**< Random draw from a parametrized model. Call via \ref apop_draw */
    void (*prep)(apop_data *data, apop_model *params);
    void (*print)(apop_model *params);
    void    *more; /**< This element is copied and freed as necessary by Apophenia's
                     model-handling functions, but is otherwise untouched. Put whatever
                     information you want here. */
    size_t  more_size; /**< If setting \c more, set this to \c sizeof(your_more_type) so
                         \ref apop_model_copy can do the \c memcpy as necessary. */
    char        error;
};

/** The global options.
  \ingroup global_vars */
typedef struct{
    int verbose; /**< Set this to zero for silent mode, one for errors and warnings. default = 0. */
    char stop_on_warning; /**< See outline page on error handling. */
    char output_type;
           /**< 's'   = to screen
                'f'   = to file
                'd'   = to db. 
                'p'   = to pipe (specifically, apop_opts.output_pipe). 
             If 'f' or 'd', then you'll need to set output_name in the apop_..._print fn. default = 's'. */
    FILE *output_pipe; /**< If printing to a pipe or FILE, set it here.  */
    char output_delimiter[100]; /**< The separator between elements of output tables. The default is "\t", but 
                                for LaTeX, use "&\t", or use "|" to get pipe-delimited output. */
    int output_append; /**< Append to output files(1), or overwrite(0)?  default = 0 */
    char input_delimiters[100]; /**< Deprecated. Please use per-function inputs to \ref apop_text_to_db and \ref apop_text_to_data. Default = "|,\t" */
    char db_name_column[300]; /**< If set, the name of the column in your tables that holds row names. */
    char db_nan[100]; /**< The string that the database takes to indicate NaN. */
    char db_engine; /**< If this is 'm', use mySQL, else use SQLite. */
    char db_user[101]; /**< Username for database login. Max 100 chars.  */
    char db_pass[101]; /**< Password for database login. Max 100 chars.  */
    FILE *log_file;  /**< The file handle for the log. Defaults to \c stderr, but change it with, e.g.,
                           <tt>apop_opts.log_file = fopen("outlog", "w");</tt> */
    int  thread_count; /**< Threads to use internally. See \ref apop_map and family.  */
    int  rng_seed;
    float version;
} apop_opts_type;

apop_opts_type apop_opts;

apop_name * apop_name_alloc(void);
int apop_name_add(apop_name * n, char const *add_me, char type);
void  apop_name_free(apop_name * free_me);
void  apop_name_print(apop_name * n);
APOP_VAR_DECLARE void  apop_name_stack(apop_name * n1, apop_name *nadd, char type1, char typeadd);
apop_name * apop_name_copy(apop_name *in);
int  apop_name_find(const apop_name *n, const char *findme, const char type);

void apop_data_add_names_base(apop_data *d, const char type, char const ** names);

/** Add a list of names to a data set.

\li Use this with a list of names that you type in yourself, like
\code
apop_data_add_names(mydata, 'c', "age", "sex", "height");
\endcode
Notice the lack of curly braces around the list.

\li You may have an array of names, probably autogenerated, that you would like to
add. In this case, make certain that the last element of the array is \c NULL, and
call the base function:
\code
char **[] colnames = {"age", "sex", "height", NULL};
apop_data_add_names_base(mydata, 'c', colnames);
\endcode
If you forget the \c NULL marker, this has good odds of segfaulting. You may prefer to use a \c for loop that inserts each name in turn using \ref apop_name_add.

\see \ref apop_name_add, although \ref apop_data_add_names will be more useful in most cases. 

 \ingroup data_struct 
*/
#define apop_data_add_names(dataset, type, ...) apop_data_add_names_base((dataset), (type), (char const*[]) {__VA_ARGS__, NULL}) 


/** Free an \ref apop_data structure.
 
As with \c free(), it is safe to send in a \c NULL pointer (in which case the function does nothing).

If the \c more pointer is not \c NULL, I will free the pointed-to data set first.
If you don't want to free data sets down the chain, set <tt>more=NULL</tt> before calling this.

\li This is actually a macro (that calls \ref apop_data_free_base to do the real work). It
sets \c freeme to \c NULL when it's done, because there's nothing safe you can do with the
freed location, and you can later safely test conditions like <tt>if (data) ...</tt>.

 \ingroup data_struct
  */
#define apop_data_free(freeme) (apop_data_free_base(freeme) ? 0 : ((freeme)= NULL))

char        apop_data_free_base(apop_data *freeme);
apop_data * apop_matrix_to_data(gsl_matrix *m);
apop_data * apop_vector_to_data(gsl_vector *v);
APOP_VAR_DECLARE apop_data * apop_data_alloc(const size_t size1, const size_t size2, const int size3);
APOP_VAR_DECLARE apop_data * apop_data_calloc(const size_t size1, const size_t size2, const int size3);
APOP_VAR_DECLARE apop_data * apop_data_stack(apop_data *m1, apop_data * m2, char posn, char inplace);
apop_data ** apop_data_split(apop_data *in, int splitpoint, char r_or_c);
apop_data * apop_data_copy(const apop_data *in);
void        apop_data_rm_columns(apop_data *d, int *drop);
void apop_data_memcpy(apop_data *out, const apop_data *in);
APOP_VAR_DECLARE double * apop_data_ptr(apop_data *data, const int row, const int col, const char *rowname, const char *colname, const char *page);
APOP_VAR_DECLARE double apop_data_get(const apop_data *data, const size_t row, const int  col, const char *rowname, const char *colname, const char *page);
APOP_VAR_DECLARE int apop_data_set(apop_data *data, const size_t row, const int col, const double val, const char *rowname, const char * colname, const char *page);
void apop_data_add_named_elmt(apop_data *d, char *name, double val);
int apop_text_add(apop_data *in, const size_t row, const size_t col, const char *fmt, ...);
apop_data * apop_text_alloc(apop_data *in, const size_t row, const size_t col);
void apop_text_free(char ***freeme, int rows, int cols);
APOP_VAR_DECLARE apop_data * apop_data_transpose(apop_data const *in, char transpose_text);
gsl_matrix * apop_matrix_realloc(gsl_matrix *m, size_t newheight, size_t newwidth);
gsl_vector * apop_vector_realloc(gsl_vector *v, size_t newheight);

#define apop_data_prune_columns(in, ...) apop_data_prune_columns_base((in), (char *[]) {__VA_ARGS__, NULL})
void apop_data_prune_columns_base(apop_data *d, char **colnames);

APOP_VAR_DECLARE apop_data * apop_data_get_page(const apop_data * data, const char * title, const char match);
apop_data * apop_data_add_page(apop_data * dataset, apop_data *newpage,const char *title);
APOP_VAR_DECLARE apop_data* apop_data_rm_page(apop_data * data, const char *title, const char free_p);
APOP_VAR_DECLARE void apop_data_rm_rows(apop_data *in, int *drop, int (*do_drop)(apop_data* ! void*), void* drop_parameter);

//in apop_asst.c:
APOP_VAR_DECLARE apop_data * apop_model_draws(apop_model *model, int count, gsl_rng *rng, apop_data *draws);


/* Convenience functions to convert among vectors (gsl_vector), matrices (gsl_matrix), 
  arrays (double **), and database tables */

//From vector
gsl_vector *apop_vector_copy(const gsl_vector *in);
APOP_VAR_DECLARE gsl_matrix * apop_vector_to_matrix(const gsl_vector *in, char row_col);

//From matrix
gsl_matrix *apop_matrix_copy(const gsl_matrix *in);
apop_data  *apop_db_to_crosstab(char *tabname, char *r1, char *r2, char *datacol);

//From array
APOP_VAR_DECLARE gsl_vector * apop_array_to_vector(double *in, int size);
#define apop_line_to_vector apop_array_to_vector

//From line
gsl_matrix * apop_line_to_matrix(double *line, int rows, int cols);
apop_data * apop_line_to_data(double *in, int vsize, int rows, int cols);

//From text
APOP_VAR_DECLARE apop_data * apop_text_to_data(char const *text_file, int has_row_names, int has_col_names, int const *field_ends, char const *delimiters);
APOP_VAR_DECLARE int apop_text_to_db(char const *text_file, char *tabname, int has_row_names, int has_col_names, char **field_names, int const *field_ends, apop_data *field_params, char *table_params, char const *delimiters);

//rank data
apop_data *apop_data_rank_expand (apop_data *in);
apop_data *apop_data_rank_compress (apop_data *in);

//From crosstabs
void apop_crosstab_to_db(apop_data *in, char *tabname, char *row_col_name, 
						char *col_col_name, char *data_col_name);

//packing data into a vector
APOP_VAR_DECLARE gsl_vector * apop_data_pack(const apop_data *in, gsl_vector *out, char all_pages, char use_info_pages);
APOP_VAR_DECLARE void apop_data_unpack(const gsl_vector *in, apop_data *d, char use_info_pages);

#define apop_vector_fill(avfin, ...) apop_vector_fill_base((avfin), (double []) {__VA_ARGS__})
#define apop_data_fill(adfin, ...) apop_data_fill_base((adfin), (double []) {__VA_ARGS__})
#define apop_text_fill(dataset, ...)   apop_text_fill_base((dataset), (char* []) {__VA_ARGS__, NULL})
    
apop_data *apop_data_fill_base(apop_data *in, double []);
gsl_vector *apop_vector_fill_base(gsl_vector *in, double []);
apop_data *apop_text_fill_base(apop_data *data, char* text[]);

int apop_data_set_row(apop_data * row, apop_data *d, int row_number);

    // Models and model support functions

apop_model apop_beta;
apop_model apop_bernoulli;
apop_model apop_binomial;
apop_model apop_chi_squared;
apop_model apop_dirichlet;
apop_model apop_exponential;
apop_model apop_f_distribution;
apop_model apop_gamma;
apop_model apop_improper_uniform;
apop_model apop_iv;
apop_model apop_kernel_density;
apop_model apop_loess;
apop_model apop_logit;
apop_model apop_lognormal;
apop_model apop_multinomial;
apop_model apop_multivariate_normal;
apop_model apop_normal;
apop_model apop_ols;
apop_model apop_pmf;
apop_model apop_poisson;
apop_model apop_probit;
apop_model apop_t_distribution;
apop_model apop_uniform;
apop_model apop_waring;
apop_model apop_wishart;
apop_model apop_wls;
apop_model apop_yule;
apop_model apop_zipf;

//model transformations
apop_model apop_coordinate_transform;
apop_model apop_mixture;
apop_model apop_stack;

/** Alias for the \ref apop_normal distribution, qv.
\hideinitializer */
#define apop_gaussian apop_normal
#define apop_OLS apop_ols
#define apop_PMF apop_pmf
#define apop_F_distribution apop_f_distribution
#define apop_WLS apop_wls
#define apop_IV apop_iv


void apop_model_free (apop_model * free_me);
void apop_model_print (apop_model * print_me);
void apop_model_show (apop_model * print_me); //deprecated
apop_model * apop_model_copy(apop_model in); //in apop_model.c
apop_model * apop_model_clear(apop_data * data, apop_model *model);

apop_model * apop_estimate(apop_data *d, apop_model m);
void apop_score(apop_data *d, gsl_vector *out, apop_model *m);
double apop_log_likelihood(apop_data *d, apop_model *m);
double apop_p(apop_data *d, apop_model *m);
double apop_cdf(apop_data *d, apop_model *m);
void apop_draw(double *out, gsl_rng *r, apop_model *m);
void apop_prep(apop_data *d, apop_model *m);
apop_model *apop_parameter_model(apop_data *d, apop_model *m);
apop_data * apop_predict(apop_data *d, apop_model *m);

apop_model *apop_beta_from_mean_var(double m, double v); //in apop_beta.c

#define apop_model_set_parameters(in, ...) apop_model_set_parameters_base((in), (double []) {__VA_ARGS__})
apop_model *apop_model_set_parameters_base(apop_model in, double ap[]);

//apop_mixture.c
/** Produce a model as a linear combination of other models. See the documentation for the \ref apop_mixture model.  */
#define apop_model_mixture(...) apop_model_mixture_base((apop_model *[]){__VA_ARGS__, NULL})
apop_model *apop_model_mixture_base(apop_model **inlist);

//apop_model_stack.c Subject to change.
apop_model *apop_model_stack(apop_model *m1, apop_model *m2);

        // Map and apply
#include <pthread.h>

    //The variadic versions, with lots of options to input extra parameters to the
    //function being mapped/applied
APOP_VAR_DECLARE apop_data * apop_map(apop_data *in, double (*fn_d)(double), double (*fn_v)(gsl_vector*), double (*fn_r)(apop_data *), double (*fn_dp)(double! void *), double (*fn_vp)(gsl_vector*! void *), double (*fn_rp)(apop_data *! void *), double (*fn_dpi)(double! void *! int), double (*fn_vpi)(gsl_vector*! void *! int), double (*fn_rpi)(apop_data*! void *! int), double (*fn_di)(double! int), double (*fn_vi)(gsl_vector*! int), double (*fn_ri)(apop_data*! int), void *param, int inplace, char part, int all_pages);
APOP_VAR_DECLARE double apop_map_sum(apop_data *in, double (*fn_d)(double), double (*fn_v)(gsl_vector*), double (*fn_r)(apop_data *), double (*fn_dp)(double! void *), double (*fn_vp)(gsl_vector*! void *), double (*fn_rp)(apop_data *! void *), double (*fn_dpi)(double! void *! int), double (*fn_vpi)(gsl_vector*! void *! int), double (*fn_rpi)(apop_data*! void *! int), double (*fn_di)(double! int), double (*fn_vi)(gsl_vector*! int), double (*fn_ri)(apop_data*! int), void *param, char part, int all_pages);

    //the specific-to-a-type versions, quicker and easier when appropriate.
gsl_vector *apop_matrix_map(const gsl_matrix *m, double (*fn)(gsl_vector*));
gsl_vector *apop_vector_map(const gsl_vector *v, double (*fn)(double));
void apop_matrix_apply(gsl_matrix *m, void (*fn)(gsl_vector*));
void apop_vector_apply(gsl_vector *v, void (*fn)(double*));
gsl_matrix * apop_matrix_map_all(const gsl_matrix *in, double (*fn)(double));
void apop_matrix_apply_all(gsl_matrix *in, void (*fn)(double *));

double apop_vector_map_sum(const gsl_vector *in, double(*fn)(double));
double apop_matrix_map_sum(const gsl_matrix *in, double (*fn)(gsl_vector*));
double apop_matrix_map_all_sum(const gsl_matrix *in, double (*fn)(double));


        // Some output routines

APOP_VAR_DECLARE void apop_plot_line_and_scatter(apop_data *data, apop_model *est, char const * output_file, FILE *output_pipe, char output_type, char output_append);
APOP_VAR_DECLARE void apop_plot_histogram(gsl_vector *data, size_t bin_count, char *with, char const *output_file, FILE *output_pipe, char output_type, char output_append);
APOP_VAR_DECLARE  void apop_plot_lattice(const apop_data *d, char const *output_file, FILE *output_pipe, char output_type, char output_append);
APOP_VAR_DECLARE void apop_plot_qq(gsl_vector *v, apop_model *m, char const *output_file, FILE *output_pipe, char output_type, char output_append, size_t bins, gsl_rng *r);
APOP_VAR_DECLARE void apop_plot_triangle(apop_data *in, char const *output_file, FILE *output_pipe, char output_type, char output_append);

APOP_VAR_DECLARE void apop_matrix_print(const gsl_matrix *data, char const *output_file, FILE *output_pipe, char output_type, char output_append);
APOP_VAR_DECLARE void apop_vector_print(gsl_vector *data, char const *output_file, FILE *output_pipe, char output_type, char output_append);
APOP_VAR_DECLARE void apop_data_print(const apop_data *data, char const *output_file, FILE *output_pipe, char output_type, char output_append);

void apop_matrix_show(const gsl_matrix *data);
void apop_vector_show(const gsl_vector *data);
void apop_data_show(const apop_data *data);

#define apop_model_coordinate_transform(...) Apop_model_copy_set(apop_coordinate_transform, apop_ct, __VA_ARGS__);
#define apop_model_dcompose(...) Apop_model_copy_set(apop_composition, apop_composition, __VA_ARGS__);

#ifdef	__cplusplus
}
#endif
#endif
