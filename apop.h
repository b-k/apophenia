
/** \file  */
/* Copyright (c) 2005--2014 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2. */

/* Here are the headers for all of apophenia's functions, typedefs, static variables and
macros. All of these begin with the apop_ (or Apop_ or APOP_) prefix.

There used to be a series of sub-headers, but they never provided any serious
benefit. Please use your text editor's word-search feature to find any elements you
may be looking for. About a third of the file is comments and doxygen documentation,
so syntax highlighting that distinguishes code from comments will also help to make
this more navigable.

*/
#pragma once
#ifdef	__cplusplus
extern "C" {
#endif

#ifndef _GNU_SOURCE
#define  _GNU_SOURCE //for asprintf
#endif

#include <assert.h>
#include <signal.h> //raise(SIGTRAP)
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>


            //////Optional arguments

/* A means of providing more script-like means of sending arguments to a function.

These macros are intended as internal.  See the documentation if you're interested in
using this mechanism in out-of-Apophenia work.
*/

#define apop_varad_head(type, name) type variadic_##name(variadic_type_##name varad_in)

#define apop_varad_declare(type, name, ...) \
    typedef struct {                        \
                __VA_ARGS__ ;               \
            } variadic_type_##name;         \
    apop_varad_head(type, name);

#define apop_varad_var(name, value) name = varad_in.name ? varad_in.name : (value);
#define apop_varad_link(name,...) variadic_##name((variadic_type_##name) {__VA_ARGS__})



            //////The types


/** This structure holds the names of the components of the \ref apop_data set. You may never have to worry about it directly, because most operations on \ref apop_data sets will take care of the names for you.
\ingroup names
*/
typedef struct{
    char *title;
	char * vector;
	char ** col;
	char ** row;
	char ** text;
	int colct, rowct, textct;
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
    char name[101]; 
    int vsize, msize1, msize2, dsize; /**< The size of the parameter set.
                     If a dimension is -1, then use yourdata->matrix->size2. For
                    anything more complex, allocate the parameter set in the prep
                    method. \c dsize is for the canonical form, and is
                    the size of the data the RNG will return. */
    apop_data *data; /**< The input data. Typically a link to what you sent to \ref apop_estimate */
    apop_data *parameters; /**< The coefficients or parameters estimated by the model. */
    apop_data *info; /**< Several pages of assorted info, perhaps including the log likelihood, AIC, BIC,
                        covariance matrix, confidence intervals, expected score. See your
                        specific model's documentation for what it puts here.
                        */
    void (*estimate)(apop_data * data, apop_model *params); 
                /**< The estimation routine. Call via \ref apop_estimate */
    long double (*p)(apop_data *d, apop_model *params);
                /**< Probability of the given data and parameterized model. Call via \ref apop_p */
    long double (*log_likelihood)(apop_data *d, apop_model *params);
                /**< Log likelihood of the given data and parameterized model. Call via \ref apop_log_likelihood */
    long double (*cdf)(apop_data *d, apop_model *params); /**< Cumulative distribution function: 
                            the integral up to the single data point you provide.  Call via \ref apop_cdf */
    long double (*constraint)(apop_data *data, apop_model *params);
    int (*draw)(double *out, gsl_rng* r, apop_model *params);
                /**< Random draw from a parametrized model. Call via \ref apop_draw */
    void (*prep)(apop_data *data, apop_model *params);
    apop_settings_type *settings;
    void *more; /**< This element is copied and freed as necessary by Apophenia's
                     model-handling functions, but is otherwise untouched. Put whatever
                     information you want here. */
    size_t more_size; /**< If setting \c more, set this to \c sizeof(your_more_type) so
                         \ref apop_model_copy can do the \c memcpy as necessary. */
    char error;
};

/** The global options.
  \ingroup global_vars */
typedef struct{
    int verbose; /**< Set this to zero for silent mode, one for errors and warnings. default = 0. */
    char stop_on_warning; /**< See outline page on error handling. */
    char output_delimiter[100]; /**< The separator between elements of output tables. The default is "\t", but 
                                for LaTeX, use "&\t", or use "|" to get pipe-delimited output. */
    char input_delimiters[100]; /**< Deprecated. Please use per-function inputs to \ref apop_text_to_db and \ref apop_text_to_data. Default = "|,\t" */
    char db_name_column[300]; /**< If set, the name of the column in your tables that holds row names. */
    char *nan_string; /**< The string used to indicate NaN. */
    char db_engine; /**< If this is 'm', use mySQL, else use SQLite. */
    char db_user[101]; /**< Username for database login. Max 100 chars.  */
    char db_pass[101]; /**< Password for database login. Max 100 chars.  */
    FILE *log_file;  /**< The file handle for the log. Defaults to \c stderr, but change it with, e.g.,
                           <tt>apop_opts.log_file = fopen("outlog", "w");</tt> */
    int  thread_count; /**< Deprecated. Use \c omp_set_num_threads(n).  */
    #if __STDC_VERSION__ > 201100L && !defined(__STDC_NO_ATOMICS__)
        _Atomic(int) rng_seed;
    #else
        int rng_seed;
    #endif
    float version;
} apop_opts_type;

extern apop_opts_type apop_opts;

apop_name * apop_name_alloc(void);
int apop_name_add(apop_name * n, char const *add_me, char type);
void  apop_name_free(apop_name * free_me);
void  apop_name_print(apop_name * n);
#ifdef APOP_NO_VARIADIC
 void  apop_name_stack(apop_name * n1, apop_name *nadd, char type1, char typeadd) ;
#else
 void apop_name_stack_base(apop_name * n1, apop_name *nadd, char type1, char typeadd) ;
 apop_varad_declare(void, apop_name_stack, apop_name * n1; apop_name *nadd; char type1; char typeadd);
#define apop_name_stack(...) apop_varad_link(apop_name_stack, __VA_ARGS__)
#endif

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
#ifdef APOP_NO_VARIADIC
 apop_data * apop_data_alloc(const size_t size1, const size_t size2, const int size3) ;
#else
 apop_data * apop_data_alloc_base(const size_t size1, const size_t size2, const int size3) ;
 apop_varad_declare(apop_data *, apop_data_alloc, const size_t size1; const size_t size2; const int size3);
#define apop_data_alloc(...) apop_varad_link(apop_data_alloc, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 apop_data * apop_data_calloc(const size_t size1, const size_t size2, const int size3) ;
#else
 apop_data * apop_data_calloc_base(const size_t size1, const size_t size2, const int size3) ;
 apop_varad_declare(apop_data *, apop_data_calloc, const size_t size1; const size_t size2; const int size3);
#define apop_data_calloc(...) apop_varad_link(apop_data_calloc, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 apop_data * apop_data_stack(apop_data *m1, apop_data * m2, char posn, char inplace) ;
#else
 apop_data * apop_data_stack_base(apop_data *m1, apop_data * m2, char posn, char inplace) ;
 apop_varad_declare(apop_data *, apop_data_stack, apop_data *m1; apop_data * m2; char posn; char inplace);
#define apop_data_stack(...) apop_varad_link(apop_data_stack, __VA_ARGS__)
#endif

apop_data ** apop_data_split(apop_data *in, int splitpoint, char r_or_c);
apop_data * apop_data_copy(const apop_data *in);
void        apop_data_rm_columns(apop_data *d, int *drop);
void apop_data_memcpy(apop_data *out, const apop_data *in);
#ifdef APOP_NO_VARIADIC
 double * apop_data_ptr(apop_data *data, int row, int col, const char *rowname, const char *colname, const char *page) ;
#else
 double * apop_data_ptr_base(apop_data *data, int row, int col, const char *rowname, const char *colname, const char *page) ;
 apop_varad_declare(double *, apop_data_ptr, apop_data *data; int row; int col; const char *rowname; const char *colname; const char *page);
#define apop_data_ptr(...) apop_varad_link(apop_data_ptr, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 double apop_data_get(const apop_data *data, size_t row, int  col, const char *rowname, const char *colname, const char *page) ;
#else
 double apop_data_get_base(const apop_data *data, size_t row, int  col, const char *rowname, const char *colname, const char *page) ;
 apop_varad_declare(double, apop_data_get, const apop_data *data; size_t row; int  col; const char *rowname; const char *colname; const char *page);
#define apop_data_get(...) apop_varad_link(apop_data_get, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 int apop_data_set(apop_data *data, size_t row, int col, const double val, const char *rowname, const char * colname, const char *page) ;
#else
 int apop_data_set_base(apop_data *data, size_t row, int col, const double val, const char *rowname, const char * colname, const char *page) ;
 apop_varad_declare(int, apop_data_set, apop_data *data; size_t row; int col; const double val; const char *rowname; const char * colname; const char *page);
#define apop_data_set(...) apop_varad_link(apop_data_set, __VA_ARGS__)
#endif

void apop_data_add_named_elmt(apop_data *d, char *name, double val);
int apop_text_add(apop_data *in, const size_t row, const size_t col, const char *fmt, ...);
apop_data * apop_text_alloc(apop_data *in, const size_t row, const size_t col);
void apop_text_free(char ***freeme, int rows, int cols);
#ifdef APOP_NO_VARIADIC
 apop_data * apop_data_transpose(apop_data *in, char transpose_text, char inplace) ;
#else
 apop_data * apop_data_transpose_base(apop_data *in, char transpose_text, char inplace) ;
 apop_varad_declare(apop_data *, apop_data_transpose, apop_data *in; char transpose_text; char inplace);
#define apop_data_transpose(...) apop_varad_link(apop_data_transpose, __VA_ARGS__)
#endif

gsl_matrix * apop_matrix_realloc(gsl_matrix *m, size_t newheight, size_t newwidth);
gsl_vector * apop_vector_realloc(gsl_vector *v, size_t newheight);

#define apop_data_prune_columns(in, ...) apop_data_prune_columns_base((in), (char *[]) {__VA_ARGS__, NULL})
apop_data* apop_data_prune_columns_base(apop_data *d, char **colnames);

#ifdef APOP_NO_VARIADIC
 apop_data * apop_data_get_page(const apop_data * data, const char * title, const char match) ;
#else
 apop_data * apop_data_get_page_base(const apop_data * data, const char * title, const char match) ;
 apop_varad_declare(apop_data *, apop_data_get_page, const apop_data * data; const char * title; const char match);
#define apop_data_get_page(...) apop_varad_link(apop_data_get_page, __VA_ARGS__)
#endif

apop_data * apop_data_add_page(apop_data * dataset, apop_data *newpage,const char *title);
#ifdef APOP_NO_VARIADIC
 apop_data* apop_data_rm_page(apop_data * data, const char *title, const char free_p) ;
#else
 apop_data* apop_data_rm_page_base(apop_data * data, const char *title, const char free_p) ;
 apop_varad_declare(apop_data*, apop_data_rm_page, apop_data * data; const char *title; const char free_p);
#define apop_data_rm_page(...) apop_varad_link(apop_data_rm_page, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 apop_data * apop_data_rm_rows(apop_data *in, int *drop, int (*do_drop)(apop_data* , void*), void* drop_parameter) ;
#else
 apop_data * apop_data_rm_rows_base(apop_data *in, int *drop, int (*do_drop)(apop_data* , void*), void* drop_parameter) ;
 apop_varad_declare(apop_data *, apop_data_rm_rows, apop_data *in; int *drop; int (*do_drop)(apop_data* , void*); void* drop_parameter);
#define apop_data_rm_rows(...) apop_varad_link(apop_data_rm_rows, __VA_ARGS__)
#endif


//in apop_asst.c:
#ifdef APOP_NO_VARIADIC
 apop_data * apop_model_draws(apop_model *model, int count, apop_data *draws) ;
#else
 apop_data * apop_model_draws_base(apop_model *model, int count, apop_data *draws) ;
 apop_varad_declare(apop_data *, apop_model_draws, apop_model *model; int count; apop_data *draws);
#define apop_model_draws(...) apop_varad_link(apop_model_draws, __VA_ARGS__)
#endif



/* Convenience functions to convert among vectors (gsl_vector), matrices (gsl_matrix), 
  arrays (double **), and database tables */

//From vector
gsl_vector *apop_vector_copy(const gsl_vector *in);
#ifdef APOP_NO_VARIADIC
 gsl_matrix * apop_vector_to_matrix(const gsl_vector *in, char row_col) ;
#else
 gsl_matrix * apop_vector_to_matrix_base(const gsl_vector *in, char row_col) ;
 apop_varad_declare(gsl_matrix *, apop_vector_to_matrix, const gsl_vector *in; char row_col);
#define apop_vector_to_matrix(...) apop_varad_link(apop_vector_to_matrix, __VA_ARGS__)
#endif


//From matrix
gsl_matrix *apop_matrix_copy(const gsl_matrix *in);
apop_data  *apop_db_to_crosstab(char *tabname, char *r1, char *r2, char *datacol);

//From array
#ifdef APOP_NO_VARIADIC
 gsl_vector * apop_array_to_vector(double *in, int size) ;
#else
 gsl_vector * apop_array_to_vector_base(double *in, int size) ;
 apop_varad_declare(gsl_vector *, apop_array_to_vector, double *in; int size);
#define apop_array_to_vector(...) apop_varad_link(apop_array_to_vector, __VA_ARGS__)
#endif

#define apop_line_to_vector apop_array_to_vector

//From text
#ifdef APOP_NO_VARIADIC
 apop_data * apop_text_to_data(char const *text_file, int has_row_names, int has_col_names, int const *field_ends, char const *delimiters) ;
#else
 apop_data * apop_text_to_data_base(char const *text_file, int has_row_names, int has_col_names, int const *field_ends, char const *delimiters) ;
 apop_varad_declare(apop_data *, apop_text_to_data, char const *text_file; int has_row_names; int has_col_names; int const *field_ends; char const *delimiters);
#define apop_text_to_data(...) apop_varad_link(apop_text_to_data, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 int apop_text_to_db(char const *text_file, char *tabname, int has_row_names, int has_col_names, char **field_names, int const *field_ends, apop_data *field_params, char *table_params, char const *delimiters, char if_table_exists) ;
#else
 int apop_text_to_db_base(char const *text_file, char *tabname, int has_row_names, int has_col_names, char **field_names, int const *field_ends, apop_data *field_params, char *table_params, char const *delimiters, char if_table_exists) ;
 apop_varad_declare(int, apop_text_to_db, char const *text_file; char *tabname; int has_row_names; int has_col_names; char **field_names; int const *field_ends; apop_data *field_params; char *table_params; char const *delimiters; char if_table_exists);
#define apop_text_to_db(...) apop_varad_link(apop_text_to_db, __VA_ARGS__)
#endif


//rank data
apop_data *apop_data_rank_expand (apop_data *in);
apop_data *apop_data_rank_compress (apop_data *in);

//From crosstabs
void apop_crosstab_to_db(apop_data *in, char *tabname, char *row_col_name, 
						char *col_col_name, char *data_col_name);

//packing data into a vector
#ifdef APOP_NO_VARIADIC
 gsl_vector * apop_data_pack(const apop_data *in, gsl_vector *out, char all_pages, char use_info_pages) ;
#else
 gsl_vector * apop_data_pack_base(const apop_data *in, gsl_vector *out, char all_pages, char use_info_pages) ;
 apop_varad_declare(gsl_vector *, apop_data_pack, const apop_data *in; gsl_vector *out; char all_pages; char use_info_pages);
#define apop_data_pack(...) apop_varad_link(apop_data_pack, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 void apop_data_unpack(const gsl_vector *in, apop_data *d, char use_info_pages) ;
#else
 void apop_data_unpack_base(const gsl_vector *in, apop_data *d, char use_info_pages) ;
 apop_varad_declare(void, apop_data_unpack, const gsl_vector *in; apop_data *d; char use_info_pages);
#define apop_data_unpack(...) apop_varad_link(apop_data_unpack, __VA_ARGS__)
#endif


#define apop_vector_fill(avfin, ...) apop_vector_fill_base((avfin), (double []) {__VA_ARGS__})
#define apop_data_fill(adfin, ...) apop_data_fill_base((adfin), (double []) {__VA_ARGS__})
#define apop_text_fill(dataset, ...)   apop_text_fill_base((dataset), (char* []) {__VA_ARGS__, NULL})

#define apop_data_falloc(sizes, ...) apop_data_fill(apop_data_alloc sizes, __VA_ARGS__)
    
apop_data *apop_data_fill_base(apop_data *in, double []);
gsl_vector *apop_vector_fill_base(gsl_vector *in, double []);
apop_data *apop_text_fill_base(apop_data *data, char* text[]);

int apop_data_set_row(apop_data * row, apop_data *d, int row_number);

// Models and model support functions

extern apop_model *apop_beta;
extern apop_model *apop_bernoulli;
extern apop_model *apop_binomial;
extern apop_model *apop_chi_squared;
extern apop_model *apop_dirichlet;
extern apop_model *apop_exponential;
extern apop_model *apop_f_distribution;
extern apop_model *apop_gamma;
extern apop_model *apop_improper_uniform;
extern apop_model *apop_iv;
extern apop_model *apop_kernel_density;
extern apop_model *apop_loess;
extern apop_model *apop_logit;
extern apop_model *apop_lognormal;
extern apop_model *apop_multinomial;
extern apop_model *apop_multivariate_normal;
extern apop_model *apop_normal;
extern apop_model *apop_ols;
extern apop_model *apop_pmf;
extern apop_model *apop_poisson;
extern apop_model *apop_probit;
extern apop_model *apop_t_distribution;
extern apop_model *apop_uniform;
//extern apop_model *apop_wishart;
extern apop_model *apop_wls;
extern apop_model *apop_yule;
extern apop_model *apop_zipf;

//model transformations
extern apop_model *apop_coordinate_transform;
extern apop_model *apop_composition;
extern apop_model *apop_dconstrain;
extern apop_model *apop_mixture;
extern apop_model *apop_stack;

/** Alias for the \ref apop_normal distribution, qv.
\hideinitializer */
#define apop_gaussian apop_normal
#define apop_OLS apop_ols
#define apop_PMF apop_pmf
#define apop_F_distribution apop_f_distribution
#define apop_WLS apop_wls
#define apop_IV apop_iv


void apop_model_free (apop_model * free_me);
void apop_model_print (apop_model * print_me, FILE *out);
void apop_model_show (apop_model * print_me); //deprecated
apop_model * apop_model_copy(apop_model *in); //in apop_model.c
apop_model * apop_model_clear(apop_data * data, apop_model *model);

apop_model * apop_estimate(apop_data *d, apop_model *m);
void apop_score(apop_data *d, gsl_vector *out, apop_model *m);
double apop_log_likelihood(apop_data *d, apop_model *m);
double apop_p(apop_data *d, apop_model *m);
double apop_cdf(apop_data *d, apop_model *m);
int apop_draw(double *out, gsl_rng *r, apop_model *m);
void apop_prep(apop_data *d, apop_model *m);
apop_model *apop_parameter_model(apop_data *d, apop_model *m);
apop_data * apop_predict(apop_data *d, apop_model *m);

apop_model *apop_beta_from_mean_var(double m, double v); //in apop_beta.c

#define apop_model_set_parameters(in, ...) apop_model_set_parameters_base((in), (double []) {__VA_ARGS__})
apop_model *apop_model_set_parameters_base(apop_model *in, double ap[]);

//apop_mixture.c
/** Produce a model as a linear combination of other models. See the documentation for the \ref apop_mixture model. 
\param ... A list of models, either all parameterized or all unparameterized. See examples in the \ref apop_mixture documentation.
 */
#define apop_model_mixture(...) apop_model_mixture_base((apop_model *[]){__VA_ARGS__, NULL})
apop_model *apop_model_mixture_base(apop_model **inlist);

//transform/apop_model_stack.c.
apop_model *apop_model_stack_base(apop_model *mlist[]);
#define apop_model_stack(...) apop_model_stack_base((apop_model *[]){__VA_ARGS__, NULL})

    //The variadic versions, with lots of options to input extra parameters to the
    //function being mapped/applied
#ifdef APOP_NO_VARIADIC
 apop_data * apop_map(apop_data *in, double (*fn_d)(double), double (*fn_v)(gsl_vector*),
                double (*fn_r)(apop_data *), double (*fn_dp)(double, void *), double (*fn_vp)(gsl_vector*, void *),
                double (*fn_rp)(apop_data *, void *), double (*fn_dpi)(double, void *, int),
                double (*fn_vpi)(gsl_vector*, void *, int), double (*fn_rpi)(apop_data*, void *, int),
                double (*fn_di)(double, int), double (*fn_vi)(gsl_vector*, int), double (*fn_ri)(apop_data*, int),
                void *param, int inplace, char part, int all_pages) ;
#else
 apop_data * apop_map_base(apop_data *in, double (*fn_d)(double), double (*fn_v)(gsl_vector*),
                double (*fn_r)(apop_data *), double (*fn_dp)(double, void *), double (*fn_vp)(gsl_vector*, void *),
                double (*fn_rp)(apop_data *, void *), double (*fn_dpi)(double, void *, int),
                double (*fn_vpi)(gsl_vector*, void *, int), double (*fn_rpi)(apop_data*, void *, int),
                double (*fn_di)(double, int), double (*fn_vi)(gsl_vector*, int), double (*fn_ri)(apop_data*, int),
                void *param, int inplace, char part, int all_pages) ;
 apop_varad_declare(apop_data *, apop_map, apop_data *in; double (*fn_d)(double); double (*fn_v)(gsl_vector*);
                double (*fn_r)(apop_data *); double (*fn_dp)(double, void *); double (*fn_vp)(gsl_vector*, void *);
                double (*fn_rp)(apop_data *, void *); double (*fn_dpi)(double, void *, int);
                double (*fn_vpi)(gsl_vector*, void *, int); double (*fn_rpi)(apop_data*, void *, int);
                double (*fn_di)(double, int); double (*fn_vi)(gsl_vector*, int); double (*fn_ri)(apop_data*, int);
                void *param; int inplace; char part; int all_pages);
#define apop_map(...) apop_varad_link(apop_map, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 double apop_map_sum(apop_data *in, double (*fn_d)(double), double (*fn_v)(gsl_vector*),
                double (*fn_r)(apop_data *), double (*fn_dp)(double, void *), double (*fn_vp)(gsl_vector*, void *),
                double (*fn_rp)(apop_data *, void *), double (*fn_dpi)(double, void *, int),
                double (*fn_vpi)(gsl_vector*, void *, int), double (*fn_rpi)(apop_data*, void *, int),
                double (*fn_di)(double, int), double (*fn_vi)(gsl_vector*, int), double (*fn_ri)(apop_data*, int),
                void *param, char part, int all_pages) ;
#else
 double apop_map_sum_base(apop_data *in, double (*fn_d)(double), double (*fn_v)(gsl_vector*),
                double (*fn_r)(apop_data *), double (*fn_dp)(double, void *), double (*fn_vp)(gsl_vector*, void *),
                double (*fn_rp)(apop_data *, void *), double (*fn_dpi)(double, void *, int),
                double (*fn_vpi)(gsl_vector*, void *, int), double (*fn_rpi)(apop_data*, void *, int),
                double (*fn_di)(double, int), double (*fn_vi)(gsl_vector*, int), double (*fn_ri)(apop_data*, int),
                void *param, char part, int all_pages) ;
 apop_varad_declare(double, apop_map_sum, apop_data *in; double (*fn_d)(double); double (*fn_v)(gsl_vector*);
                double (*fn_r)(apop_data *); double (*fn_dp)(double, void *); double (*fn_vp)(gsl_vector*, void *);
                double (*fn_rp)(apop_data *, void *); double (*fn_dpi)(double, void *, int);
                double (*fn_vpi)(gsl_vector*, void *, int); double (*fn_rpi)(apop_data*, void *, int);
                double (*fn_di)(double, int); double (*fn_vi)(gsl_vector*, int); double (*fn_ri)(apop_data*, int);
                void *param; char part; int all_pages);
#define apop_map_sum(...) apop_varad_link(apop_map_sum, __VA_ARGS__)
#endif


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

#ifdef APOP_NO_VARIADIC
 void apop_plot_line_and_scatter(apop_data *data, apop_model *est, char const * output_name, FILE *output_pipe, char output_type, char output_append) ;
#else
 void apop_plot_line_and_scatter_base(apop_data *data, apop_model *est, char const * output_name, FILE *output_pipe, char output_type, char output_append) ;
 apop_varad_declare(void, apop_plot_line_and_scatter, apop_data *data; apop_model *est; char const * output_name; FILE *output_pipe; char output_type; char output_append);
#define apop_plot_line_and_scatter(...) apop_varad_link(apop_plot_line_and_scatter, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 void apop_plot_histogram(gsl_vector *data, size_t bin_count, char *with, char const *output_name, FILE *output_pipe, char output_type, char output_append) ;
#else
 void apop_plot_histogram_base(gsl_vector *data, size_t bin_count, char *with, char const *output_name, FILE *output_pipe, char output_type, char output_append) ;
 apop_varad_declare(void, apop_plot_histogram, gsl_vector *data; size_t bin_count; char *with; char const *output_name; FILE *output_pipe; char output_type; char output_append);
#define apop_plot_histogram(...) apop_varad_link(apop_plot_histogram, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 void apop_plot_lattice(const apop_data *d, char const *output_name, FILE *output_pipe, char output_type, char output_append) ;
#else
 void apop_plot_lattice_base(const apop_data *d, char const *output_name, FILE *output_pipe, char output_type, char output_append) ;
 apop_varad_declare(void, apop_plot_lattice, const apop_data *d; char const *output_name; FILE *output_pipe; char output_type; char output_append);
#define apop_plot_lattice(...) apop_varad_link(apop_plot_lattice, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 void apop_plot_qq(gsl_vector *v, apop_model *m, char const *output_name, FILE *output_pipe, char output_type, char output_append, size_t bins, gsl_rng *r) ;
#else
 void apop_plot_qq_base(gsl_vector *v, apop_model *m, char const *output_name, FILE *output_pipe, char output_type, char output_append, size_t bins, gsl_rng *r) ;
 apop_varad_declare(void, apop_plot_qq, gsl_vector *v; apop_model *m; char const *output_name; FILE *output_pipe; char output_type; char output_append; size_t bins; gsl_rng *r);
#define apop_plot_qq(...) apop_varad_link(apop_plot_qq, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 void apop_plot_triangle(apop_data *in, char const *output_name, FILE *output_pipe, char output_type, char output_append) ;
#else
 void apop_plot_triangle_base(apop_data *in, char const *output_name, FILE *output_pipe, char output_type, char output_append) ;
 apop_varad_declare(void, apop_plot_triangle, apop_data *in; char const *output_name; FILE *output_pipe; char output_type; char output_append);
#define apop_plot_triangle(...) apop_varad_link(apop_plot_triangle, __VA_ARGS__)
#endif


#ifdef APOP_NO_VARIADIC
 void apop_matrix_print(const gsl_matrix *data, char const *output_name, FILE *output_pipe, char output_type, char output_append) ;
#else
 void apop_matrix_print_base(const gsl_matrix *data, char const *output_name, FILE *output_pipe, char output_type, char output_append) ;
 apop_varad_declare(void, apop_matrix_print, const gsl_matrix *data; char const *output_name; FILE *output_pipe; char output_type; char output_append);
#define apop_matrix_print(...) apop_varad_link(apop_matrix_print, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 void apop_vector_print(gsl_vector *data, char const *output_name, FILE *output_pipe, char output_type, char output_append) ;
#else
 void apop_vector_print_base(gsl_vector *data, char const *output_name, FILE *output_pipe, char output_type, char output_append) ;
 apop_varad_declare(void, apop_vector_print, gsl_vector *data; char const *output_name; FILE *output_pipe; char output_type; char output_append);
#define apop_vector_print(...) apop_varad_link(apop_vector_print, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 void apop_data_print(const apop_data *data, char const *output_name, FILE *output_pipe, char output_type, char output_append) ;
#else
 void apop_data_print_base(const apop_data *data, char const *output_name, FILE *output_pipe, char output_type, char output_append) ;
 apop_varad_declare(void, apop_data_print, const apop_data *data; char const *output_name; FILE *output_pipe; char output_type; char output_append);
#define apop_data_print(...) apop_varad_link(apop_data_print, __VA_ARGS__)
#endif


void apop_matrix_show(const gsl_matrix *data);
void apop_vector_show(const gsl_vector *data);
void apop_data_show(const apop_data *data);


        //statistics
#ifdef APOP_NO_VARIADIC
 double apop_vector_mean(gsl_vector const *v, gsl_vector const *weights);
#else
 double apop_vector_mean_base(gsl_vector const *v, gsl_vector const *weights);
 apop_varad_declare(double, apop_vector_mean, gsl_vector const *v; gsl_vector const *weights);
#define apop_vector_mean(...) apop_varad_link(apop_vector_mean, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 double apop_vector_var(gsl_vector const *v, gsl_vector const *weights);
#else
 double apop_vector_var_base(gsl_vector const *v, gsl_vector const *weights);
 apop_varad_declare(double, apop_vector_var, gsl_vector const *v; gsl_vector const *weights);
#define apop_vector_var(...) apop_varad_link(apop_vector_var, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 double apop_vector_skew_pop(gsl_vector const *v, gsl_vector const *weights);
#else
 double apop_vector_skew_pop_base(gsl_vector const *v, gsl_vector const *weights);
 apop_varad_declare(double, apop_vector_skew_pop, gsl_vector const *v; gsl_vector const *weights);
#define apop_vector_skew_pop(...) apop_varad_link(apop_vector_skew_pop, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 double apop_vector_kurtosis_pop(gsl_vector const *v, gsl_vector const *weights);
#else
 double apop_vector_kurtosis_pop_base(gsl_vector const *v, gsl_vector const *weights);
 apop_varad_declare(double, apop_vector_kurtosis_pop, gsl_vector const *v; gsl_vector const *weights);
#define apop_vector_kurtosis_pop(...) apop_varad_link(apop_vector_kurtosis_pop, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 double apop_vector_cov(gsl_vector const *v1, gsl_vector const *v2,
                                         gsl_vector const *weights);
#else
 double apop_vector_cov_base(gsl_vector const *v1, gsl_vector const *v2,
                                         gsl_vector const *weights);
 apop_varad_declare(double, apop_vector_cov, gsl_vector const *v1; gsl_vector const *v2;
                                         gsl_vector const *weights);
#define apop_vector_cov(...) apop_varad_link(apop_vector_cov, __VA_ARGS__)
#endif


#ifdef APOP_NO_VARIADIC
 double apop_vector_distance(const gsl_vector *ina, const gsl_vector *inb, const char metric, const double norm) ;
#else
 double apop_vector_distance_base(const gsl_vector *ina, const gsl_vector *inb, const char metric, const double norm) ;
 apop_varad_declare(double, apop_vector_distance, const gsl_vector *ina; const gsl_vector *inb; const char metric; const double norm);
#define apop_vector_distance(...) apop_varad_link(apop_vector_distance, __VA_ARGS__)
#endif


#ifdef APOP_NO_VARIADIC
 void apop_vector_normalize(gsl_vector *in, gsl_vector **out, const char normalization_type) ;
#else
 void apop_vector_normalize_base(gsl_vector *in, gsl_vector **out, const char normalization_type) ;
 apop_varad_declare(void, apop_vector_normalize, gsl_vector *in; gsl_vector **out; const char normalization_type);
#define apop_vector_normalize(...) apop_varad_link(apop_vector_normalize, __VA_ARGS__)
#endif

void apop_matrix_normalize(gsl_matrix *data, const char row_or_col, const char normalization);

apop_data * apop_data_covariance(const apop_data *in);
apop_data * apop_data_correlation(const apop_data *in);
long double apop_matrix_sum(const gsl_matrix *m);
double apop_matrix_mean(const gsl_matrix *data);
void apop_matrix_mean_and_var(const gsl_matrix *data, double *mean, double *var);
apop_data * apop_data_summarize(apop_data *data);

apop_data *apop_test_fisher_exact(apop_data *intab); //in apop_fisher.c

//from apop_t_f_chi.c:
#ifdef APOP_NO_VARIADIC
 int apop_matrix_is_positive_semidefinite(gsl_matrix *m, char semi) ;
#else
 int apop_matrix_is_positive_semidefinite_base(gsl_matrix *m, char semi) ;
 apop_varad_declare(int, apop_matrix_is_positive_semidefinite, gsl_matrix *m; char semi);
#define apop_matrix_is_positive_semidefinite(...) apop_varad_link(apop_matrix_is_positive_semidefinite, __VA_ARGS__)
#endif

double apop_matrix_to_positive_semidefinite(gsl_matrix *m);
long double apop_multivariate_gamma(double a, int p);
long double apop_multivariate_lngamma(double a, int p);

//apop_tests.c
apop_data *	apop_t_test(gsl_vector *a, gsl_vector *b);
apop_data *	apop_paired_t_test(gsl_vector *a, gsl_vector *b);
#ifdef APOP_NO_VARIADIC
 apop_data* apop_anova(char *table, char *data, char *grouping1, char *grouping2) ;
#else
 apop_data* apop_anova_base(char *table, char *data, char *grouping1, char *grouping2) ;
 apop_varad_declare(apop_data*, apop_anova, char *table; char *data; char *grouping1; char *grouping2);
#define apop_anova(...) apop_varad_link(apop_anova, __VA_ARGS__)
#endif

#define apop_ANOVA apop_anova
#ifdef APOP_NO_VARIADIC
 apop_data * apop_f_test (apop_model *est, apop_data *contrast) ;
#else
 apop_data * apop_f_test_base(apop_model *est, apop_data *contrast) ;
 apop_varad_declare(apop_data *, apop_f_test, apop_model *est; apop_data *contrast);
#define apop_f_test(...) apop_varad_link(apop_f_test, __VA_ARGS__)
#endif

#define apop_F_test apop_f_test

//from the regression code:
#define apop_estimate_r_squared(in) apop_estimate_coefficient_of_determination(in)

apop_data * apop_text_unique_elements(const apop_data *d, size_t col);
gsl_vector * apop_vector_unique_elements(const gsl_vector *v);
#ifdef APOP_NO_VARIADIC
 apop_data * apop_data_to_factors(apop_data *data, char intype, int incol, int outcol) ;
#else
 apop_data * apop_data_to_factors_base(apop_data *data, char intype, int incol, int outcol) ;
 apop_varad_declare(apop_data *, apop_data_to_factors, apop_data *data; char intype; int incol; int outcol);
#define apop_data_to_factors(...) apop_varad_link(apop_data_to_factors, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 apop_data * apop_data_get_factor_names(apop_data *data, int col, char type) ;
#else
 apop_data * apop_data_get_factor_names_base(apop_data *data, int col, char type) ;
 apop_varad_declare(apop_data *, apop_data_get_factor_names, apop_data *data; int col; char type);
#define apop_data_get_factor_names(...) apop_varad_link(apop_data_get_factor_names, __VA_ARGS__)
#endif


#ifdef APOP_NO_VARIADIC
 apop_data * apop_data_to_dummies(apop_data *d, int col, char type, int keep_first, char append, char remove) ;
#else
 apop_data * apop_data_to_dummies_base(apop_data *d, int col, char type, int keep_first, char append, char remove) ;
 apop_varad_declare(apop_data *, apop_data_to_dummies, apop_data *d; int col; char type; int keep_first; char append; char remove);
#define apop_data_to_dummies(...) apop_varad_link(apop_data_to_dummies, __VA_ARGS__)
#endif


#ifdef APOP_NO_VARIADIC
 double apop_kl_divergence(apop_model *from, apop_model *to, int draw_ct, gsl_rng *rng, apop_model *top, apop_model *bottom) ;
#else
 double apop_kl_divergence_base(apop_model *from, apop_model *to, int draw_ct, gsl_rng *rng, apop_model *top, apop_model *bottom) ;
 apop_varad_declare(double, apop_kl_divergence, apop_model *from; apop_model *to; int draw_ct; gsl_rng *rng; apop_model *top; apop_model *bottom);
#define apop_kl_divergence(...) apop_varad_link(apop_kl_divergence, __VA_ARGS__)
#endif


apop_data *apop_estimate_coefficient_of_determination (apop_model *);
void apop_estimate_parameter_tests (apop_model *est);


//Bootstrapping & RNG
apop_data * apop_jackknife_cov(apop_data *data, apop_model *model);
#ifdef APOP_NO_VARIADIC
 apop_data * apop_bootstrap_cov(apop_data *data, apop_model *model, gsl_rng* rng, int iterations, char keep_boots, char ignore_nans) ;
#else
 apop_data * apop_bootstrap_cov_base(apop_data *data, apop_model *model, gsl_rng* rng, int iterations, char keep_boots, char ignore_nans) ;
 apop_varad_declare(apop_data *, apop_bootstrap_cov, apop_data *data; apop_model *model; gsl_rng* rng; int iterations; char keep_boots; char ignore_nans);
#define apop_bootstrap_cov(...) apop_varad_link(apop_bootstrap_cov, __VA_ARGS__)
#endif

gsl_rng *apop_rng_alloc(int seed);
double apop_rng_GHgB3(gsl_rng * r, double* a); //in apop_asst.c

#define apop_rng_get_thread(thread_in) apop_rng_get_thread_base(#thread_in[0]=='\0' ? -1: (thread_in+0))
gsl_rng *apop_rng_get_thread_base(int thread);

int apop_arms_draw (double *out, gsl_rng *r, apop_model *m); //apop_arms.h


    // maximum likelihod estimation related functions

#ifdef APOP_NO_VARIADIC
 gsl_vector * apop_numerical_gradient(apop_data * data, apop_model* model, double delta) ;
#else
 gsl_vector * apop_numerical_gradient_base(apop_data * data, apop_model* model, double delta) ;
 apop_varad_declare(gsl_vector *, apop_numerical_gradient, apop_data * data; apop_model* model; double delta);
#define apop_numerical_gradient(...) apop_varad_link(apop_numerical_gradient, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 apop_data * apop_model_hessian(apop_data * data, apop_model *model, double delta) ;
#else
 apop_data * apop_model_hessian_base(apop_data * data, apop_model *model, double delta) ;
 apop_varad_declare(apop_data *, apop_model_hessian, apop_data * data; apop_model *model; double delta);
#define apop_model_hessian(...) apop_varad_link(apop_model_hessian, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 apop_data * apop_model_numerical_covariance(apop_data * data, apop_model *model, double delta) ;
#else
 apop_data * apop_model_numerical_covariance_base(apop_data * data, apop_model *model, double delta) ;
 apop_varad_declare(apop_data *, apop_model_numerical_covariance, apop_data * data; apop_model *model; double delta);
#define apop_model_numerical_covariance(...) apop_varad_link(apop_model_numerical_covariance, __VA_ARGS__)
#endif


void apop_maximum_likelihood(apop_data * data, apop_model *dist);

#ifdef APOP_NO_VARIADIC
 apop_model * apop_estimate_restart (apop_model *e, apop_model *copy, char * starting_pt, double boundary) ;
#else
 apop_model * apop_estimate_restart_base(apop_model *e, apop_model *copy, char * starting_pt, double boundary) ;
 apop_varad_declare(apop_model *, apop_estimate_restart, apop_model *e; apop_model *copy; char * starting_pt; double boundary);
#define apop_estimate_restart(...) apop_varad_link(apop_estimate_restart, __VA_ARGS__)
#endif


//in apop_linear_constraint.c
#ifdef APOP_NO_VARIADIC
 long double  apop_linear_constraint(gsl_vector *beta, apop_data * constraint, double margin) ;
#else
 long double apop_linear_constraint_base(gsl_vector *beta, apop_data * constraint, double margin) ;
 apop_varad_declare(long double, apop_linear_constraint, gsl_vector *beta; apop_data * constraint; double margin);
#define apop_linear_constraint(...) apop_varad_link(apop_linear_constraint, __VA_ARGS__)
#endif


//in apop_model_fix_params.c
apop_model * apop_model_fix_params(apop_model *model_in);
apop_model * apop_model_fix_params_get_base(apop_model *model_in);






            //////vtables

/*
This declares the vtable macros for each procedure that uses the mechanism.

--We want to have type-checking on the functions put into the vtables. Type checking
happens only with functions, not macros, so we need a type_check function for every
vtable.

--Only once in your codebase, you'll need to #define Declare_type_checking_fns to
actually define the type checking function. Everywhere else, the function is merely
declared.

--All other uses point to having a macro, such as using __VA_ARGS__ to allow any sort
of inputs to the hash.

--We want to have such a macro for every vtable. That means that we need a macro
to write macros. We can't do that with C macros.  Thus, this file uses m4 macros to
generate C macros.

--After the m4 definition of make_vtab_fns, each new vtable requires a typedef, a hash
definition, and a call to make_vtab_fns to do the rest.
*/


int apop_vtable_add(char const *tabname, void *fn_in, unsigned long hash);
void *apop_vtable_get(char const *tabname, unsigned long hash);
int apop_vtable_drop(char const *tabname, unsigned long hash);

typedef apop_model *(*apop_update_type)(apop_data *, apop_model* , apop_model*);
#define apop_update_hash(m1, m2) (          \
           ((m1)->log_likelihood ? (size_t)(m1)->log_likelihood : \
            (m1)->p              ? (size_t)(m1)->p*33 : \
            (m1)->draw           ? (size_t)(m1)->draw*33*27 \
                                 : 33*27*19) \
          +((m2)->log_likelihood ? (size_t)(m2)->log_likelihood : \
            (m2)->p              ? (size_t)(m2)->p*33 : \
            (m2)->draw           ? (size_t)(m2)->draw*33*27 \
                                 : 33*27*19 \
           ) * 37)
#ifdef Declare_type_checking_fns
void apop_update_type_check(apop_update_type in){ };
#else
void apop_update_type_check(apop_update_type in);
#endif
#define apop_update_vtable_add(fn, ...) apop_update_type_check(fn), apop_vtable_add("apop_update", fn, apop_update_hash(__VA_ARGS__))
#define apop_update_vtable_get(...) apop_vtable_get("apop_update", apop_update_hash(__VA_ARGS__))
#define apop_update_vtable_drop(...) apop_vtable_drop("apop_update", apop_update_hash(__VA_ARGS__))

typedef void (*apop_score_type)(apop_data *d, gsl_vector *gradient, apop_model *params);
#define apop_score_hash(m1) ((size_t)((m1)->log_likelihood ? (m1)->log_likelihood : (m1)->p))
#ifdef Declare_type_checking_fns
void apop_score_type_check(apop_score_type in){ };
#else
void apop_score_type_check(apop_score_type in);
#endif
#define apop_score_vtable_add(fn, ...) apop_score_type_check(fn), apop_vtable_add("apop_score", fn, apop_score_hash(__VA_ARGS__))
#define apop_score_vtable_get(...) apop_vtable_get("apop_score", apop_score_hash(__VA_ARGS__))
#define apop_score_vtable_drop(...) apop_vtable_drop("apop_score", apop_score_hash(__VA_ARGS__))

typedef apop_model* (*apop_parameter_model_type)(apop_data *, apop_model *);
#define apop_parameter_model_hash(m1) ((size_t)((m1)->log_likelihood ? (m1)->log_likelihood : (m1)->p)*33 + (m1)->estimate ? (size_t)(m1)->estimate: 27)
#ifdef Declare_type_checking_fns
void apop_parameter_model_type_check(apop_parameter_model_type in){ };
#else
void apop_parameter_model_type_check(apop_parameter_model_type in);
#endif
#define apop_parameter_model_vtable_add(fn, ...) apop_parameter_model_type_check(fn), apop_vtable_add("apop_parameter_model", fn, apop_parameter_model_hash(__VA_ARGS__))
#define apop_parameter_model_vtable_get(...) apop_vtable_get("apop_parameter_model", apop_parameter_model_hash(__VA_ARGS__))
#define apop_parameter_model_vtable_drop(...) apop_vtable_drop("apop_parameter_model", apop_parameter_model_hash(__VA_ARGS__))

typedef apop_data * (*apop_predict_type)(apop_data *d, apop_model *params);
#define apop_predict_hash(m1) ((size_t)((m1)->log_likelihood ? (m1)->log_likelihood : (m1)->p)*33 + (m1)->estimate ? (size_t)(m1)->estimate: 27)
#ifdef Declare_type_checking_fns
void apop_predict_type_check(apop_predict_type in){ };
#else
void apop_predict_type_check(apop_predict_type in);
#endif
#define apop_predict_vtable_add(fn, ...) apop_predict_type_check(fn), apop_vtable_add("apop_predict", fn, apop_predict_hash(__VA_ARGS__))
#define apop_predict_vtable_get(...) apop_vtable_get("apop_predict", apop_predict_hash(__VA_ARGS__))
#define apop_predict_vtable_drop(...) apop_vtable_drop("apop_predict", apop_predict_hash(__VA_ARGS__))

typedef void (*apop_model_print_type)(apop_model *params, FILE *out);
#define apop_model_print_hash(m1) ((m1)->log_likelihood ? (size_t)(m1)->log_likelihood : \
            (m1)->p ? (size_t)(m1)->p*33 : \
            (m1)->estimate ? (size_t)(m1)->estimate*33*33 : \
            (m1)->draw ? (size_t)(m1)->draw*33*27  : \
            (m1)->cdf ? (size_t)(m1)->cdf*27*27  \
            : 27)
#ifdef Declare_type_checking_fns
void apop_model_print_type_check(apop_model_print_type in){ };
#else
void apop_model_print_type_check(apop_model_print_type in);
#endif
#define apop_model_print_vtable_add(fn, ...) apop_model_print_type_check(fn), apop_vtable_add("apop_model_print", fn, apop_model_print_hash(__VA_ARGS__))
#define apop_model_print_vtable_get(...) apop_vtable_get("apop_model_print", apop_model_print_hash(__VA_ARGS__))
#define apop_model_print_vtable_drop(...) apop_vtable_drop("apop_model_print", apop_model_print_hash(__VA_ARGS__))





        //////Asst


double apop_generalized_harmonic(int N, double s) __attribute__ ((__pure__));

apop_data * apop_test_anova_independence(apop_data *d);
#define apop_test_ANOVA_independence(d) apop_test_anova_independence(d)

#ifdef APOP_NO_VARIADIC
 int apop_regex(const char *string, const char* regex, apop_data **substrings, const char use_case) ;
#else
 int apop_regex_base(const char *string, const char* regex, apop_data **substrings, const char use_case) ;
 apop_varad_declare(int, apop_regex, const char *string; const char* regex; apop_data **substrings; const char use_case);
#define apop_regex(...) apop_varad_link(apop_regex, __VA_ARGS__)
#endif


int apop_system(const char *fmt, ...) __attribute__ ((format (printf,1,2)));

//Histograms and PMFs
gsl_vector * apop_vector_moving_average(gsl_vector *, size_t);
apop_data * apop_histograms_test_goodness_of_fit(apop_model *h0, apop_model *h1);
apop_data * apop_test_kolmogorov(apop_model *m1, apop_model *m2);
apop_data *apop_data_pmf_compress(apop_data *in);
#ifdef APOP_NO_VARIADIC
 apop_data * apop_data_to_bins(apop_data *indata, apop_data *binspec, int bin_count, char close_top_bin) ;
#else
 apop_data * apop_data_to_bins_base(apop_data *indata, apop_data *binspec, int bin_count, char close_top_bin) ;
 apop_varad_declare(apop_data *, apop_data_to_bins, apop_data *indata; apop_data *binspec; int bin_count; char close_top_bin);
#define apop_data_to_bins(...) apop_varad_link(apop_data_to_bins, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 apop_model * apop_model_to_pmf(apop_model *model, apop_data *binspec, long int draws, int bin_count, gsl_rng *rng) ;
#else
 apop_model * apop_model_to_pmf_base(apop_model *model, apop_data *binspec, long int draws, int bin_count, gsl_rng *rng) ;
 apop_varad_declare(apop_model *, apop_model_to_pmf, apop_model *model; apop_data *binspec; long int draws; int bin_count; gsl_rng *rng);
#define apop_model_to_pmf(...) apop_varad_link(apop_model_to_pmf, __VA_ARGS__)
#endif


//text conveniences
#ifdef APOP_NO_VARIADIC
 char* apop_text_paste(apop_data const*strings, char *between, char *before, char *after, char *between_cols, int (*prune)(apop_data* , int , int , void*), void* prune_parameter) ;
#else
 char* apop_text_paste_base(apop_data const*strings, char *between, char *before, char *after, char *between_cols, int (*prune)(apop_data* , int , int , void*), void* prune_parameter) ;
 apop_varad_declare(char*, apop_text_paste, apop_data const*strings; char *between; char *before; char *after; char *between_cols; int (*prune)(apop_data* , int , int , void*); void* prune_parameter);
#define apop_text_paste(...) apop_varad_link(apop_text_paste, __VA_ARGS__)
#endif

/** Notify the user of errors, warning, or debug info. 

 \param verbosity   At what verbosity level should the user be warned? E.g., if level==2, then print iff apop_opts.verbosity >= 2.
 \param ... The message to write to STDERR (presuming the verbosity level is high enough). This can be a printf-style format with following arguments. You can produce much more informative error messages this way, e.g., \c apop_notify(0, "Beta is %g but should be greater than zero.", beta);.
*/
#define Apop_notify(verbosity, ...) {\
    if (apop_opts.verbose != -1 && apop_opts.verbose >= verbosity) {  \
        if (!apop_opts.log_file) apop_opts.log_file = stderr; \
        fprintf(apop_opts.log_file, "%s: ", __func__); fprintf(apop_opts.log_file, __VA_ARGS__); fprintf(apop_opts.log_file, "\n");   \
        fflush(apop_opts.log_file); \
} }

#define Apop_maybe_abort(level) \
            {if ((level == -5 && apop_opts.stop_on_warning!='n')                \
            || (apop_opts.verbose >= level && apop_opts.stop_on_warning == 'v') \
            || (apop_opts.stop_on_warning=='w') ) \
                raise(SIGTRAP);}

/** Execute an action and print a message to \c stderr (or the current \c FILE handle held by <tt>apop_opts.log_file</tt>).
 Intended for leaving a function on failure.
 
\param test The expression that, if true, triggers all the action.
\param onfail If the assertion fails, do this. E.g., <tt>out->error='x'; return GSL_NAN</tt>. Notice that it is OK to include several lines of semicolon-separated code here, but if you have a lot to do, the most readable option may be <tt>goto outro</tt>, plus an appropriately-labeled section at the end of your function.
\param level Print the warning message only if \ref apop_opts_type "apop_opts.verbose" is greater than or equal to this. Zero usually works, but for minor infractions use one.
\param ... The error message in printf form, plus any arguments to be inserted into the printf string. I'll provide the function name and a carriage return.

\li If \c apop_opts.stop_on_warning is nonzero and not <tt>'v'</tt>, then a failed test halts via \c abort(), even if the <tt>apop_opts.verbose</tt> level is set so that the warning message doesn't print to screen. Use this when running via debugger.
\li If \c apop_opts.stop_on_warning is <tt>'v'</tt>, then a failed test halts via \c abort() iff the verbosity level is high enough to print the error.
*/
#define Apop_stopif(test, onfail, level, ...) do {\
     if (test) {  \
        Apop_notify(level,  __VA_ARGS__);   \
        Apop_maybe_abort(level)  \
        onfail;  \
    } } while(0)

#define apop_errorlevel -5

//For use in stopif, to return a blank apop_data set with an error attached.
#define apop_return_data_error(E) {apop_data *out=apop_data_alloc(); out->error='E'; return out;}

/* The Apop_stopif macro is currently favored, but there's a long history of prior
   error-handling setups. Consider all of the Assert... macros below to be deprecated.
*/
#define Apop_assert_c(test, returnval, level, ...) \
    Apop_stopif(!(test), return returnval, level, __VA_ARGS__)

#define Apop_assert(test, ...) Apop_assert_c((test), 0, apop_errorlevel, __VA_ARGS__)

//For things that return void. Transitional and deprecated at birth.
#define Apop_assert_n(test, ...) Apop_assert_c((test),  , apop_errorlevel, __VA_ARGS__)
#define Apop_assert_nan(test, ...) Apop_assert_c((test), GSL_NAN, apop_errorlevel, __VA_ARGS__)
#define Apop_assert_negone(test, ...) Apop_assert_c((test), -1, apop_errorlevel, __VA_ARGS__)

//Missing data
#ifdef APOP_NO_VARIADIC
 apop_data * apop_data_listwise_delete(apop_data *d, char inplace) ;
#else
 apop_data * apop_data_listwise_delete_base(apop_data *d, char inplace) ;
 apop_varad_declare(apop_data *, apop_data_listwise_delete, apop_data *d; char inplace);
#define apop_data_listwise_delete(...) apop_varad_link(apop_data_listwise_delete, __VA_ARGS__)
#endif

apop_model * apop_ml_impute(apop_data *d, apop_model* meanvar);
#define apop_ml_imputation(d, m) apop_ml_impute(d, m)

#ifdef APOP_NO_VARIADIC
 apop_model *apop_model_metropolis(apop_data *d, apop_model *m, gsl_rng* rng);
#else
 apop_model * apop_model_metropolis_base(apop_data *d, apop_model *m, gsl_rng* rng);
 apop_varad_declare(apop_model *, apop_model_metropolis, apop_data *d; apop_model *m; gsl_rng* rng);
#define apop_model_metropolis(...) apop_varad_link(apop_model_metropolis, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 apop_model * apop_update(apop_data *data, apop_model *prior, apop_model *likelihood, gsl_rng *rng) ;
#else
 apop_model * apop_update_base(apop_data *data, apop_model *prior, apop_model *likelihood, gsl_rng *rng) ;
 apop_varad_declare(apop_model *, apop_update, apop_data *data; apop_model *prior; apop_model *likelihood; gsl_rng *rng);
#define apop_update(...) apop_varad_link(apop_update, __VA_ARGS__)
#endif


#ifdef APOP_NO_VARIADIC
 double apop_test(double statistic, char *distribution, double p1, double p2, char tail) ;
#else
 double apop_test_base(double statistic, char *distribution, double p1, double p2, char tail) ;
 apop_varad_declare(double, apop_test, double statistic; char *distribution; double p1; double p2; char tail);
#define apop_test(...) apop_varad_link(apop_test, __VA_ARGS__)
#endif


//Sorting (apop_asst.c)
#ifdef APOP_NO_VARIADIC
 double * apop_vector_percentiles(gsl_vector *data, char rounding)  ;
#else
 double * apop_vector_percentiles_base(gsl_vector *data, char rounding)  ;
 apop_varad_declare(double *, apop_vector_percentiles, gsl_vector *data; char rounding);
#define apop_vector_percentiles(...) apop_varad_link(apop_vector_percentiles, __VA_ARGS__)
#endif


//apop_sort.c
#ifdef APOP_NO_VARIADIC
 apop_data *apop_data_sort(apop_data *data, apop_data *sort_order, char asc, char inplace, double *col_order);
#else
 apop_data * apop_data_sort_base(apop_data *data, apop_data *sort_order, char asc, char inplace, double *col_order);
 apop_varad_declare(apop_data *, apop_data_sort, apop_data *data; apop_data *sort_order; char asc; char inplace; double *col_order);
#define apop_data_sort(...) apop_varad_link(apop_data_sort, __VA_ARGS__)
#endif


//raking
#ifdef APOP_NO_VARIADIC
 apop_data * apop_rake(char const *margin_table, char * const*var_list, 
                    int var_ct, char const *all_vars, char * const *contrasts, int contrast_ct, 
                    char const *structural_zeros, int max_iterations, double tolerance, 
                    char const *count_col, int run_number, char const *init_table, 
                    char const *init_count_col, double nudge, char const* table_name) ;
#else
 apop_data * apop_rake_base(char const *margin_table, char * const*var_list, 
                    int var_ct, char const *all_vars, char * const *contrasts, int contrast_ct, 
                    char const *structural_zeros, int max_iterations, double tolerance, 
                    char const *count_col, int run_number, char const *init_table, 
                    char const *init_count_col, double nudge, char const* table_name) ;
 apop_varad_declare(apop_data *, apop_rake, char const *margin_table; char * const*var_list; 
                    int var_ct; char const *all_vars; char * const *contrasts; int contrast_ct; 
                    char const *structural_zeros; int max_iterations; double tolerance; 
                    char const *count_col; int run_number; char const *init_table; 
                    char const *init_count_col; double nudge; char const* table_name);
#define apop_rake(...) apop_varad_link(apop_rake, __VA_ARGS__)
#endif



#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>


    //First, some linear algebra utilities

double apop_det_and_inv(const gsl_matrix *in, gsl_matrix **out, int calc_det, int calc_inv);
#ifdef APOP_NO_VARIADIC
 apop_data * apop_dot(const apop_data *d1, const apop_data *d2, char form1, char form2) ;
#else
 apop_data * apop_dot_base(const apop_data *d1, const apop_data *d2, char form1, char form2) ;
 apop_varad_declare(apop_data *, apop_dot, const apop_data *d1; const apop_data *d2; char form1; char form2);
#define apop_dot(...) apop_varad_link(apop_dot, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 int         apop_vector_bounded(const gsl_vector *in, long double max) ;
#else
 int apop_vector_bounded_base(const gsl_vector *in, long double max) ;
 apop_varad_declare(int, apop_vector_bounded, const gsl_vector *in; long double max);
#define apop_vector_bounded(...) apop_varad_link(apop_vector_bounded, __VA_ARGS__)
#endif

gsl_matrix * apop_matrix_inverse(const gsl_matrix *in) ;
double      apop_matrix_determinant(const gsl_matrix *in) ;
//apop_data*  apop_sv_decomposition(gsl_matrix *data, int dimensions_we_want);
#ifdef APOP_NO_VARIADIC
 apop_data *  apop_matrix_pca(gsl_matrix *data, int const dimensions_we_want) ;
#else
 apop_data * apop_matrix_pca_base(gsl_matrix *data, int const dimensions_we_want) ;
 apop_varad_declare(apop_data *, apop_matrix_pca, gsl_matrix *data; int const dimensions_we_want);
#define apop_matrix_pca(...) apop_varad_link(apop_matrix_pca, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 gsl_vector * apop_vector_stack(gsl_vector *v1, gsl_vector * v2, char inplace) ;
#else
 gsl_vector * apop_vector_stack_base(gsl_vector *v1, gsl_vector * v2, char inplace) ;
 apop_varad_declare(gsl_vector *, apop_vector_stack, gsl_vector *v1; gsl_vector * v2; char inplace);
#define apop_vector_stack(...) apop_varad_link(apop_vector_stack, __VA_ARGS__)
#endif

#ifdef APOP_NO_VARIADIC
 gsl_matrix * apop_matrix_stack(gsl_matrix *m1, gsl_matrix * m2, char posn, char inplace) ;
#else
 gsl_matrix * apop_matrix_stack_base(gsl_matrix *m1, gsl_matrix * m2, char posn, char inplace) ;
 apop_varad_declare(gsl_matrix *, apop_matrix_stack, gsl_matrix *m1; gsl_matrix * m2; char posn; char inplace);
#define apop_matrix_stack(...) apop_varad_link(apop_matrix_stack, __VA_ARGS__)
#endif

gsl_matrix * apop_matrix_rm_columns(gsl_matrix *in, int *drop);

void apop_vector_log(gsl_vector *v);
void apop_vector_log10(gsl_vector *v);
void apop_vector_exp(gsl_vector *v);

#define APOP_SUBMATRIX(m, srow, scol, nrows, ncols, o) gsl_matrix apop_mm_##o = gsl_matrix_submatrix((m), (srow), (scol), (nrows),(ncols)).matrix;\
gsl_matrix * o = &( apop_mm_##o );

#define Apop_row_v(m, row, v) Apop_matrix_row((m)->matrix, row, v)

#define Apop_col_v(m, col, v) gsl_vector apop_vv_##v = ((col) == -1) ? (gsl_vector){} : gsl_matrix_column((m)->matrix, (col)).vector;\
gsl_vector * v = ((col)==-1) ? (m)->vector : &( apop_vv_##v );

#define Apop_row_tv(m, row, v) gsl_vector apop_vv_##v = gsl_matrix_row((m)->matrix, apop_name_find((m)->names, row, 'r')).vector;\
gsl_vector * v = &( apop_vv_##v );

#define Apop_col_tv(m, col, v) gsl_vector apop_vv_##v = gsl_matrix_column((m)->matrix, apop_name_find((m)->names, col, 'c')).vector;\
gsl_vector * v = &( apop_vv_##v );

/* //-1st column is the vector.
#define Apop_matrix_col(m, col, v) gsl_vector apop_vv_##v = (col < 0) ? *((m)->vector) : gsl_matrix_column((m)->matrix, (col)).vector;\
gsl_vector * v = &( apop_vv_##v );*/

#define Apop_rows(d, rownum, len, outd) \
    gsl_vector apop_dd_##outd##_v = ((d)->vector && (d)->vector->size > (rownum)+(len)-1)  \
                                    ? gsl_vector_subvector((d)->vector, (rownum), (len)).vector\
                                    : (gsl_vector) { };\
    gsl_vector apop_dd_##outd##_w = ((d)->weights && (d)->weights->size > (rownum)+(len)-1)  \
                                    ? gsl_vector_subvector((d)->weights, (rownum), (len)).vector\
                                    : (gsl_vector) { };\
    gsl_matrix apop_dd_##outd##_m = ((d)->matrix && (d)->matrix->size1 > (rownum)+(len)-1)  \
                                ? gsl_matrix_submatrix((d)->matrix, rownum, 0, (len), (d)->matrix->size2).matrix\
                                : (gsl_matrix) { };             \
    apop_name apop_dd_##outd##_n = !((d)->names) ? (apop_name) {} :              \
            (apop_name){                                                         \
                .title = (d)->names->title,                                      \
                .vector = (d)->names->vector,                                    \
                .col = (d)->names->col,                                          \
                .row = ((d)->names->row && (d)->names->rowct > rownum) ? &((d)->names->row[rownum]) : NULL,  \
                .text = (d)->names->text,                                        \
                .colct = (d)->names->colct,                                      \
                .rowct = (d)->names->row ? (GSL_MIN(len, GSL_MAX((d)->names->rowct - rownum, 0)))      \
                                          : 0,                                   \
                .textct = (d)->names->textct };                                  \
    apop_data apop_dd_##outd = (apop_data){                                      \
                .vector= apop_dd_##outd##_v.size ? &apop_dd_##outd##_v : NULL,   \
                .weights=apop_dd_##outd##_w.size ? &apop_dd_##outd##_w : NULL ,  \
                .matrix = apop_dd_##outd##_m.size1 ? &apop_dd_##outd##_m : NULL, \
                .textsize[0]=(d)->textsize[0] ? (len) : 0,                       \
                .textsize[1]=(d)->textsize[1],                                   \
                .text = (d)->text ? &((d)->text[rownum]) : NULL,                 \
                .names= (d)->names ? &apop_dd_##outd##_n : NULL };               \
    apop_data *outd =  &apop_dd_##outd;

#define Apop_row(d, row, outd) Apop_rows(d, row, 1, outd)

#define Apop_row_t(d, rowname, outd) int apop_row_##outd = apop_name_find((d)->names, rowname, 'r'); Apop_rows(d, apop_row##outd, 1, outd)

#define Apop_col_t(d, colname, outd) int apop_col_##outd = apop_name_find((d)->names, colname, 'c'); Apop_cols(d, apop_col##outd, 1, outd)

#define Apop_cols(d, colnum, len, outd) \
    gsl_matrix apop_dd_##outd##_m = ((d)->matrix && (d)->matrix->size2 > (colnum)+(len)-1)  \
                                ? gsl_matrix_submatrix((d)->matrix, 0, colnum, (d)->matrix->size1, (len)).matrix\
                                : (gsl_matrix) { };             \
    apop_name apop_dd_##outd##_n = !((d)->names) ? (apop_name) {} :              \
            (apop_name){                                                         \
                .title = (d)->names->title,                                      \
                .vector = NULL,                                    \
                .row = (d)->names->row,                                          \
                .col = ((d)->names->col && (d)->names->colct > colnum) ? &((d)->names->col[colnum]) : NULL,  \
                .text = NULL,                                                    \
                .rowct = (d)->names->rowct,                                      \
                .colct = (d)->names->col ? (GSL_MIN(len, GSL_MAX((d)->names->colct - colnum, 0)))      \
                                          : 0,                                   \
                .textct = (d)->names->textct };                                  \
    apop_data apop_dd_##outd = (apop_data){                                      \
                .vector= NULL,                                                   \
                .weights= (d)->weights,                                          \
                .matrix = apop_dd_##outd##_m.size1 ? &apop_dd_##outd##_m : NULL, \
                .textsize[0] = 0,                                                \
                .textsize[1] = 0,                                                \
                .text = NULL,                                                    \
                .names= (d)->names ? &apop_dd_##outd##_n : NULL };               \
    apop_data *outd =  &apop_dd_##outd;

#define Apop_col(d, col, outd) Apop_cols(d, col, 1, outd)

#define APOP_COL Apop_col
#define apop_col Apop_col
#define APOP_COLS Apop_cols
#define apop_cols Apop_cols
#define APOP_COL_T Apop_col_t
#define apop_col_t Apop_col_t
#define APOP_COL_TV Apop_col_tv
#define apop_col_tv Apop_col_tv
#define APOP_COL_V Apop_col_v
#define apop_col_v Apop_col_v

#define APOP_ROW Apop_row
#define apop_row Apop_row
#define Apop_data_row Apop_row   #deprecated
#define APOP_ROWS Apop_rows
#define apop_rows Apop_rows
#define APOP_ROW_T Apop_row_t
#define apop_row_t Apop_row_t
#define APOP_ROW_TV Apop_row_tv
#define apop_row_tv Apop_row_tv
#define APOP_ROW_V Apop_row_v
#define apop_row_v Apop_row_v

/** View a single row of a \c gsl_matrix as a \c gsl_vector. This 
 is a convenience macro wrapping \c gsl_matrix_row. 
 
\param m The \c gsl_matrix
\param row The number of the desired row. 
\param v The name of the vector view that will be created.

An: example
\code 
gsl_matrix *m = [fill matrix here];
Apop_matrix_row(m, 2, rowtwo);
Apop_matrix_row(m, 3, rowthree);
printf("The correlation coefficient between rows two "
       "and three is %g.\n", apop_vector_correlation(rowtwo, rowthree));
\endcode 
\see Apop_rows, Apop_row, Apop_row_v, Apop_row_tv, Apop_row_t
*/
#define Apop_matrix_row(m, row, v) gsl_vector apop_vv_##v = gsl_matrix_row((m), (row)).vector;\
gsl_vector * v = &( apop_vv_##v );

/** View a single column of a \c gsl_matrix as a \c gsl_vector. This 
 is a convenience macro wrapping \c gsl_matrix_column. 
 
\param m The \c gsl_matrix
\param col The number of the desired column.
\param v The name of the vector view that will be created.

An: example
\code 
gsl_matrix *m = [fill matrix here];
Apop_matrix_col(m, 2, coltwo);
Apop_matrix_col(m, 3, colthree);
printf("The correlation coefficient between columns two "
       "and three is %g.\n", apop_vector_correlation(coltwo, colthree));
\endcode 
\see Apop_cols, Apop_col, Apop_col_v, Apop_col_tv, Apop_col_t
*/
#define Apop_matrix_col(m, col, v) gsl_vector apop_vv_##v = gsl_matrix_column((m), (col)).vector;\
gsl_vector * v = &( apop_vv_##v );

#define Apop_submatrix APOP_SUBMATRIX
#define APOP_MATRIX_ROW Apop_matrix_row 
#define apop_matrix_row Apop_matrix_row 
#define APOP_MATRIX_COL Apop_matrix_col 
#define apop_matrix_col Apop_matrix_col 


long double apop_vector_sum(const gsl_vector *in);
double apop_vector_var_m(const gsl_vector *in, const double mean);
double apop_vector_correlation(const gsl_vector *ina, const gsl_vector *inb);
double apop_vector_kurtosis(const gsl_vector *in);
double apop_vector_skew(const gsl_vector *in);

#define apop_sum(in) apop_vector_sum(in)
#define apop_var(in) apop_vector_var(in) 
#define apop_mean(in) apop_vector_mean(in)

/** Find the mean of the input vector.

*/

/** \def Apop_submatrix(m, srow, scol, nrow, ncol)
Generate a subview of a submatrix within a \c gsl_matrix. Like \ref Apop_matrix_row, \ref Apop_row, et al., the view is an automatically-allocated variable that is lost once the program flow leaves the scope in which it is declared.

 \param m The root matrix
 \param srow the first row (in the root matrix) of the top of the submatrix
 \param scol the first column (in the root matrix) of the left edge of the submatrix
 \param nrow number of rows in the submatrix
 \param ncol number of columns in the submatrix
\hideinitializer */

/** \def Apop_row_t(m, row_name, v)
 After this call, \c v will hold a view of an \ref apop_data set consisting only of the <tt>row</tt>th row of the \ref apop_data set \c m.
 Unlike \ref Apop_row, the second argument is a row name, that I'll look up using \ref apop_name_find.
\see Apop_rows, Apop_row, Apop_row_v, Apop_row_tv, Apop_matrix_row
\hideinitializer */

/** \def Apop_col_t(m, col_name, v)
 After this call, \c v will hold a view of an \ref apop_data set consisting only of  vector view of the <tt>col</tt>th column of the \ref apop_data set \c m.
 Unlike \ref Apop_col, the second argument is a column name, that I'll look up using \ref apop_name_find.
\see Apop_cols, Apop_col, Apop_col_v, Apop_col_tv, Apop_matrix_col
\hideinitializer */

/** \def Apop_row_tv(m, row_name, v)
 After this call, \c v will hold a vector view of the <tt>row</tt>th row of the \ref apop_data set \c m.
 Unlike \ref Apop_row_v, the second argument is a row name, that I'll look up using \ref apop_name_find.
\see Apop_rows, Apop_row, Apop_row_v, Apop_row_t, Apop_matrix_row
\hideinitializer */

/** \def Apop_col_tv(m, col_name, v)
 After this call, \c v will hold a vector view of the <tt>col</tt>th column of the \ref apop_data set \c m.
 Unlike \ref Apop_col_v, the second argument is a column name, that I'll look up using \ref apop_name_find.
\see Apop_cols, Apop_col, Apop_col_v, Apop_col_t, Apop_matrix_col
\hideinitializer */

/** \def Apop_row_v(m, row, v)
 After this call, \c v will hold a vector view of the <tt>row</tt>th row of the \ref apop_data set \c m.
 This is like \ref Apop_row, but the output is a \c gsl_vector, not a full \ref apop_data set. It is like \ref Apop_matrix_row, but the input is a \ref apop_data set, not a \c gsl_matrix.
\see Apop_rows, Apop_row, Apop_row_tv, Apop_row_t, Apop_matrix_row
\hideinitializer */

/** \def Apop_col_v(m, col, v)
 After this call, \c v will hold a vector view of the <tt>col</tt>th column of the \ref apop_data set \c m.
 This is like \ref Apop_col, but the output is a \c gsl_vector, not a full \ref apop_data set. It is like \ref Apop_matrix_col, but the input is a \ref apop_data set, not a \c gsl_matrix.
\see Apop_cols, Apop_col, Apop_col_tv, Apop_col_t, Apop_matrix_col
\hideinitializer */

/** \def Apop_cols(d, col, len, outd)
A macro to generate a temporary view of \ref apop_data set \c d, beginning at column \c col and having length \c len. 
After this call, \c outd will be a pointer to this temporary
view, that you can use as you would any \ref apop_data set. However, 
it expires as soon as the program leaves the current scope (like with the usual automatically declared vars). 
\see Apop_col, Apop_col_v, Apop_col_tv, Apop_col_t, Apop_matrix_col
\hideinitializer */

/** \def Apop_col(d, col, v)
A macro to generate a temporary one-column view of \ref apop_data set \c d, pulling out only
column \c col. 
After this call, \c outd will be a pointer to this temporary
view, that you can use as you would any \ref apop_data set. This macro expands to
<tt>Apop_cols(d, col, 1, outd)</tt>.
\see Apop_cols, Apop_col_v, Apop_col_tv, Apop_col_t, Apop_matrix_col
\hideinitializer */

/** \def Apop_rows(d, row, len, outd)
A macro to generate a temporary view of \ref apop_data set \c d, beginning at row \c row
and having length \c len. 
After this call, \c outd will be a pointer to this temporary
view, that you can use as you would any \ref apop_data set. However, 
it expires as soon as the program leaves the current scope (like with the usual automatically declared vars). 
\see Apop_row, Apop_row_v, Apop_row_tv, Apop_row_t, Apop_matrix_row
\hideinitializer */

/** \def Apop_row(d, row, outd)
A macro to generate a temporary one-row view of \ref apop_data set \c d, pulling out only
row \c row. 
After this call, \c outd will be a pointer to this temporary
view, that you can use as you would any \ref apop_data set. This macro expands to
<tt>Apop_rows(d, row, 1, outd)</tt>.
\see Apop_rows, Apop_row_v, Apop_row_tv, Apop_row_t, Apop_matrix_row
\hideinitializer */

/** \def apop_mean(v)
 Returns the mean of the elements of the vector \c v.
\hideinitializer */



        //////database utilities

#ifdef APOP_NO_VARIADIC
 int apop_table_exists(char const *name, char remove) ;
#else
 int apop_table_exists_base(char const *name, char remove) ;
 apop_varad_declare(int, apop_table_exists, char const *name; char remove);
#define apop_table_exists(...) apop_varad_link(apop_table_exists, __VA_ARGS__)
#endif


int apop_db_open(char const *filename);
#ifdef APOP_NO_VARIADIC
 int apop_db_close(char vacuum) ;
#else
 int apop_db_close_base(char vacuum) ;
 apop_varad_declare(int, apop_db_close, char vacuum);
#define apop_db_close(...) apop_varad_link(apop_db_close, __VA_ARGS__)
#endif


int apop_query(const char *q, ...) __attribute__ ((format (printf,1,2)));
gsl_matrix * apop_query_to_matrix(const char * fmt, ...) __attribute__ ((format (printf,1,2)));
apop_data * apop_query_to_text(const char * fmt, ...) __attribute__ ((format (printf,1,2)));
apop_data * apop_query_to_data(const char * fmt, ...) __attribute__ ((format (printf,1,2)));
apop_data * apop_query_to_mixed_data(const char *typelist, const char * fmt, ...) __attribute__ ((format (printf,2,3)));
gsl_vector * apop_query_to_vector(const char * fmt, ...) __attribute__ ((format (printf,1,2)));
double apop_query_to_float(const char * fmt, ...) __attribute__ ((format (printf,1,2)));

int apop_data_to_db(const apop_data *set, const char *tabname, char);

double apop_db_t_test(char * tab1, char *col1, char *tab2, char *col2);
double apop_db_paired_t_test(char * tab1, char *col1, char *col2);


        //////Settings groups

    //Part I: macros and fns for getting/setting settings groups and elements

void * apop_settings_get_grp(apop_model *m, char *type, char fail);
void apop_settings_remove_group(apop_model *m, char *delme);
void apop_settings_copy_group(apop_model *outm, apop_model *inm, char *copyme);
void *apop_settings_group_alloc(apop_model *model, char *type, void *free_fn, void *copy_fn, void *the_group);
apop_model *apop_settings_group_alloc_wm(apop_model *model, char *type, void *free_fn, void *copy_fn, void *the_group);

/** Retrieves a settings group from a model.  See \ref Apop_settings_get
 to just pull a single item from within the settings group.

  If it isn't found, then it returns NULL, so you can easily put it in a conditional like 
  \code 
  if (!apop_settings_get_group(m, "apop_ols")) ...
  \endcode
\hideinitializer \ingroup settings
 */
#define Apop_settings_get_group(m, type) apop_settings_get_grp(m, #type, 'c')

/** Removes a settings group from a model's list. 
 
  If the so-named group is not found, do nothing.
\hideinitializer \ingroup settings
 */
#define Apop_settings_rm_group(m, type) apop_settings_remove_group(m, #type)

/** Add a settings group. The first two arguments (the model you are
 attaching to and the settings group name) are mandatory, and then you
 can use the \ref designated syntax to specify default values (if any).
 \return A pointer to the newly-prepped group.

\li If a settings group of the given type is already attached to the model, 
the previous version is removed. Use \ref Apop_settings_get to check whether a group
of the given type is already attached to a model, and \ref Apop_settings_set to modify
an existing group.

\hideinitializer \ingroup settings
 */
#define Apop_settings_add_group(model, type, ...)  \
    apop_settings_group_alloc(model, #type, type ## _settings_free, type ## _settings_copy, type ##_settings_init ((type ## _settings) {__VA_ARGS__}))

/** Copy a model and add a settings group. Useful for models that require a settings group to function. See \ref Apop_settings_add_group.

 \return A pointer to the newly-prepped model.
\hideinitializer \ingroup settings
 */
#define apop_model_copy_set(model, type, ...)  \
    apop_settings_group_alloc_wm(apop_model_copy(model), #type, type ## _settings_free, type ## _settings_copy, type ##_settings_init ((type ## _settings) {__VA_ARGS__}))

/** Retrieves a setting from a model.  See \ref Apop_settings_get_group to pull the entire group.
\hideinitializer \ingroup settings
 */
#define Apop_settings_get(model, type, setting)  \
    (((type ## _settings *) apop_settings_get_grp(model, #type, 'f'))->setting)

/** Modifies a single element of a settings group to the given value. 

\li If <tt>model==NULL</tt>, fails silently. 
\li If <tt>model!=NULL</tt> but the given settings group is not found attached to the model, set <tt>model->error='s'</tt>.
\hideinitializer \ingroup settings
 */
#define Apop_settings_set(model, type, setting, data)   \
    do {                                                \
        if (!(model)) continue; /* silent fail. */      \
        type ## _settings *apop_tmp_settings = apop_settings_get_grp(model, #type, 'c');  \
        Apop_stopif(!apop_tmp_settings, (model)->error='s', 0, "You're trying to modify a setting in " \
                        #model "'s setting group of type " #type " but that model doesn't have such a group."); \
    apop_tmp_settings->setting = (data);                \
    } while (0);

/** \cond doxy_ignore */
#define Apop_settings_add Apop_settings_set
#define APOP_SETTINGS_ADD Apop_settings_set
#define apop_settings_set Apop_settings_set
#define APOP_SETTINGS_GET Apop_settings_get
#define apop_settings_get Apop_settings_get
#define APOP_SETTINGS_ADD_GROUP Apop_settings_add_group
#define apop_settings_add_group Apop_settings_add_group
#define APOP_SETTINGS_GET_GROUP Apop_settings_get_group
#define apop_settings_get_group Apop_settings_get_group
#define APOP_SETTINGS_RM_GROUP Apop_settings_rm_group
#define apop_settings_rm_group Apop_settings_rm_group
#define Apop_model_copy_set apop_model_copy_set

//deprecated:
#define Apop_model_add_group Apop_settings_add_group

/** \endcond */ //End of Doxygen ignore.

#define Apop_settings_declarations(ysg) \
   ysg##_settings * ysg##_settings_init(ysg##_settings); \
   void * ysg##_settings_copy(ysg##_settings *); \
   void ysg##_settings_free(ysg##_settings *);

/** A convenience macro for declaring the initialization function for a new settings group.
 See the documentation outline -> models -> model settings -> writing new settings group for details.

  This sets the defaults for every element in the structure, so you will want a line for every element of your structure (except the ones that default to NULL, which have already been set as such).

  \code
  Apop_settings_init (ysg, 
        Apop_varad_set(size1, 99);
        Apop_varad_set(size2, 2.3);
        Apop_varad_set(dataset, apop_data_alloc(out->size1, out->size2));
    )
  \endcode
  If you need them, the input is a structure named \c in, and the output a pointer-to-struct named \c out.
*/
#define Apop_settings_init(name, ...)   \
    name##_settings *name##_settings_init(name##_settings in) {       \
        name##_settings *out = malloc(sizeof(name##_settings));     \
        *out = in; \
        __VA_ARGS__;            \
        return out; \
    }

#define Apop_varad_set(var, value) (out)->var = (in).var ? (in).var : (value);

/** A convenience macro for declaring the copy function for a new settings group.
 See the documentation outline -> models -> model settings -> writing new settings group for details.

  To just do a direct copy, the default works; let your settings group be named \c ysg:
  \code
Apop_settings_copy (ysg, )
  \endcode
  generates a function that allocates space for a new settings group and copies all elements from the input group to the output group.

  The space after the comma indicates that there is no new procedural code. If you want to add some, feel free. E.g.,
  \code
Apop_settings_copy (ysg, 
    if (!in->score)
        out->score = 1;
    out->data_owner = 0;
)
  \endcode
  The names \c in and \c out are built into the macro.
*/
#define Apop_settings_copy(name, ...) \
    void * name##_settings_copy(name##_settings *in) {\
        name##_settings *out = malloc(sizeof(name##_settings)); \
        *out = *in; \
        __VA_ARGS__;    \
        return out;     \
    }

/** A convenience macro for declaring the delete function for a new settings group.
 See the documentation outline -> models -> model settings -> writing new settings group for details.

If you don't have internal structure elements to free, let your settings group be named \c ysg:
  \code
  Apop_settings_free (ysg, )
  \endcode
  generates a function that simply frees the input settings group.

  If your structure is pointing to other structures that need to be freed first, then add them after that comma:
  \code
Apop_settings_copy (ysg, 
    apop_data_free(in->dataset);
)
  \endcode
  The name \c in is built into the macro.
*/
#define Apop_settings_free(name, ...) \
    void name##_settings_free(name##_settings *in) {\
        __VA_ARGS__;    \
        free(in);  \
    }

        //Part II: the details of extant settings groups.


/** The settings for maximum likelihood estimation (including simulated annealing).
\ingroup settings */
typedef struct{
    double      *starting_pt;   /**< An array of doubles (i.e., <tt>double*</tt>) suggesting a starting point. 
                                  If NULL, use an all-ones vector.  Note that if \c v is a \c gsl_vector, then 
                                  \c v->data is of the right form (provided \c v is not a slice of a matrix).*/
    char *method; /**< The method to be used for the optimization. All strings are case-insensitive.

        <table>
<tr>
<td> String <td></td> Name  <td></td>  Notes
</td> </tr>
                                     
<tr><td> "NM simplex" </td><td> Nelder-Mead simplex </td><td> Does not use gradients at all. Can sometimes get stuck.</td></tr>

<tr><td> "FR cg"  </td><td> Conjugate gradient (Fletcher-Reeves) (default) </td><td> CG methods use derivatives. The converge to the optimum of a quadratic function in one step; performance degrades as the objective digresses from quadratic.</td></tr>

<tr><td> "BFGS cg" </td><td> Broyden-Fletcher-Goldfarb-Shanno conjugate gradient        </td><td>  </td></tr>

<tr><td> "PR cg"  </td><td> Polak-Ribiere conjugate gradient  </td><td>  </td></tr>

<tr><td> "Annealing"  </td><td> \ref simanneal "simulated annealing"         </td><td> Slow but works for objectives of arbitrary complexity, including stochastic objectives.</td></tr>

<tr><td> "Newton"</td><td> Newton's method  </td><td> Search by finding a root of the derivative. Expects that gradient is reasonably well-behaved. </td></tr>

<tr><td> "Newton hybrid"</td><td> Newton's method/gradient descent hybrid        </td><td>  Find a root of the derivative via the Hybrid method </td> If Newton proposes stepping outside of a certain interval, use an alternate method. See <a href="https://www.gnu.org/software/gsl/manual/gsl-ref_35.html#SEC494">the GSL manual</a> for discussion.</tr>

<tr><td> "Newton hybrid no scale"</td><td>  Newton's method/gradient descent hybrid with spherical scale</td><td>  As above, but use a simplified trust region. </td></tr>
</table> */
    double      step_size, /**< the initial step size. */
                tolerance, /**< the precision the minimizer uses. Only vaguely related to the precision of the actual variables. */
delta;
    int         max_iterations; /**< Ignored by simulated annealing. Other methods halt if
                                 they do this many iterations without finding an optimum. */
    int         verbose; /**<	Give status updates as we go.  This is orthogonal to the 
                                <tt>apop_opts.verbose</tt> setting. */
    double      dim_cycle_tolerance; /**< If zero (the default), the usual procedure.
                             If \f$>0\f$, cycle across dimensions: fix all but the first dimension at the starting
                             point, optimize only the first dim. Then fix the all but the second dim, and optimize the
                             second dim. Continue through all dims, until the log likelihood at the outset of one cycle
                             through the dimensions is within this amount of the previous cycle's log likelihood. There
                             will be at least two cycles.
                             */
//simulated annealing (also uses step_size);
    int         n_tries, iters_fixed_T;
    double      k, t_initial, mu_t, t_min ;
    gsl_rng     *rng;
    apop_data   **path;    /**< If not \c NULL, record each vector tried by the optimizer as one row of this \ref apop_data set.
                              Each row of the \c matrix element holds the vector tried; the corresponding element in the \c vector is the evaluated value at that vector (after out-of-constraints penalties have been subtracted).
                              A new \ref apop_data set is allocated at the pointer you send in. This data set has no names; add them as desired. Sample use:
\code                              
apop_data *mypath;
Apop_model_add_group(mymodel, apop_mle, .path=&mypath);
apop_model *out = apop_estimate(mydata, mymodel);
apop_data_print(mypath, .output_name="search");
apop_data_free(mypath);
\endcode                              
                              
*/
} apop_mle_settings;

/** Settings for least-squares type models 
\ingroup settings */
typedef struct {
    int destroy_data; /**< If 'y', then the input data set may be normalized or otherwise mangled */
    apop_data *instruments; /**< Use for the \ref apop_iv regression, qv. */
    char want_cov; /**< Deprecated. Please use \ref apop_parts_wanted_settings. */
    char want_expected_value; /**< Deprecated. Please use \ref apop_parts_wanted_settings. */
    apop_model *input_distribution; /**< The distribution of \f$P(Y|X)\f$ is specified by the model, but the distribution of \f$X\f$ is not.  */
} apop_lm_settings;

/** The default is for the estimation routine to give some auxiliary information,
  such as a covariance matrix, predicted values, and common hypothesis tests.
  Some uses of a model depend on these items, but if they are a waste
  of time for your purposes, this settings group gives a quick way to bypass them all.

  Simply adding this settings group to your model without changing any default values---
  \code
  Apop_model_add_group(your_model, apop_parts_wanted);
  \endcode
  ---will turn off all of the auxiliary calculations covered, because the default value
  for all the switches is <tt>'n'</tt>, indicating that all elements are not wanted.

  From there, you can change some of the default <tt>'n'</tt>s to <tt>'y'</tt>s to retain some but not all auxiliary elements.  If you just want the parameters themselves and the covariance matrix:
  \code
  Apop_model_add_group(your_model, apop_parts_wanted, .covariance='y');
  \endcode

  \li Not all models support this, although the models with especially compute-intensive
  auxiliary info do (e.g., the maximum likelihood estimation system). Check the model's documentation. 

  \li Tests may depend on covariance, so <tt>.covariance='n', .tests='y'</tt> may be 
  treated as <tt>.covariance='y', .tests='y'</tt>.
\ingroup settings
*/
typedef struct {
    //init/copy/free are in apop_mle.c
    char covariance;    /*< If 'y', calculate the covariance matrix. Default 'n'. */
    char predicted;/*< If 'y', calculate the predicted values. This is typically as many
                     items as rows in your data set. Default 'n'. */
    char tests;/*< If 'y', run any hypothesis tests offered by the model's estimation routine. Default 'n'. */
    char info;/*< If 'y', add an info table with elements such as log likelihood or AIC. Default 'n'. */
} apop_parts_wanted_settings;

/** Some CDFs use random draws; some use closed-form models. 
  \ingroup settings */
typedef struct {
    int draws;  /**< For random draw methods, how many draws? Default: 10,000.*/
    gsl_rng *rng; /**< For random draw methods. See \ref apop_rng_get_thread on the default. */
    apop_model *cdf_model; /**< For use by individual models as they see fit. Default=\c NULL. */
    gsl_matrix *draws_made; /**< A store of random draws that I will count up to report the CDF. Need only be generated once, and so stored here. */
    int *draws_refcount; /**< For internal use.*/
} apop_cdf_settings;


/** Settings for getting parameter models (i.e. the distribution of parameter estimates)
  \ingroup settings */
typedef struct {
    apop_model *base;
    int index;
    gsl_rng *rng;
    int draws;
} apop_pm_settings;


/** Settings to accompany the \ref apop_pmf. */
typedef struct {
    gsl_vector *cmf;  /**< A cumulative mass function, for the purposes of making random draws.*/
    char draw_index;  /**< If \c 'y', then draws from the PMF return the integer index of the row drawn. 
                           If \c 'n' (the default), then return the data in the vector/matrix elements of the data set. */
    long double total_weight; /**< Keep the total weight, in case the input weights aren't normalized to sum to one. */
    int *cmf_refct;    /**< For internal use, so I can garbage-collect the CMF when needed. */
} apop_pmf_settings;


/** Settings for the \ref apop_kernel_density model. 

  \ingroup settings */
typedef struct{
    apop_data *base_data; /**< The data that will be smoothed by the KDE. */
    apop_model *base_pmf; /**< I actually need the data in a \ref apop_pmf. You can give
                            that to me explicitly, or I can wrap the .base_data in a PMF.  */
    apop_model *kernel; /**< The distribution to be centered over each data point. Default, 
                                    \ref apop_normal with std dev 1. */
    void (*set_fn)(apop_data*, apop_model*); /**< The function I will use for each data
                                                  point to center the kernel over each point.*/
    int own_pmf, own_kernel; /**< For internal use only. */
}apop_kernel_density_settings;


/** Method settings for a model to be put through Bayesian updating. 
\ingroup settings 
 */
typedef struct{
    apop_data *data;
    long int periods; /**< For how many steps should the MCMC chain run? */
    double burnin; /**< What <em>percentage</em> of the periods should be ignored
                         as initialization. That is, this is a number between zero and one. */
    int histosegments; /**< If outputting a binned PMF, how many segments should it have? */
    char method;
    apop_model *proposal; /**< The distribution from which test parameters will be drawn. After getting the draw using the \c draw method of the proposal, the base model's \c parameters element is filled using \ref apop_data_fill.
If \c NULL, \ref apop_model_metropolis will use a Multivariate Normal with the appropriate dimension, mean zero, and covariance matrix I. If not \c NULL, be sure to parameterize your model with an initial position. */
    void (*step_fn)(apop_data *, apop_model*); /**< Modifies the parameters of the
        proposal distribution given a successful draw. Thus, this function writes the
        drawn data point to the parameter set. If the draw is a scalar, the default
        function sets the 0th element of the model's \c parameter set with the draw
        (works for the \ref apop_normal and other models). If the draw has multiple
        dimensions, they are all copied to the parameter set, which must have the same
        size. */
    apop_data *starting_pt; /**< Deprecated and ignored. Starting point is drawn from your proposal distribution. */
    double last_ll; /**< If you have already run mcmc, the last log likelihood in the chain.*/
    apop_model *pmf; /**< If you have already run mcmc, I keep a pointer to the model so far here. Use \ref apop_model_metropolis_draw to get one more draw.*/
    apop_model *base_model; /**< The model you provided with a \c log_likelihood or \c p element (which need not sum to one). You do not have to set this: if it is \c NULL on input to \ref apop_model_metropolis, I will fill it in.*/
    int proposal_is_cp; /**< For internal use. */
} apop_mcmc_settings;

//Loess, including the old FORTRAN-to-C.
struct loess_struct {
	struct {
		long    n, p;
        double  *y, *x;
		double	*weights;
	} in;
	struct {
	        double  span;
	        long    degree;
	        long    normalize;
	        long    parametric[8];
	        long    drop_square[8];
	        char    *family;
	} model;
	struct {
	        char    *surface;
	        char    *statistics;
	        double  cell;
	        char    *trace_hat;
	        long    iterations;
	} control;
	struct {
		long	*parameter, *a;
		double	*xi, *vert, *vval;
	} kd_tree;
	struct {
		double	*fitted_values;
        double  *fitted_residuals;
		double  enp, s;
		double  one_delta, two_delta;
		double	*pseudovalues;
		double	trace_hat;
		double	*diagonal;
		double	*robust;
		double  *divisor;
	} out;
};

/** The code for the loess system is based on FORTRAN code from 1988,
overhauled in 1992, linked in to Apophenia in 2009. The structure that
does all the work, then, is a \c loess_struct that you should
basically take as opaque. 

The useful settings from that struct re-appear in the \ref
apop_loess_settings struct so you can set them directly, and then the
settings init function will copy your preferences into the working struct.

The documentation for the elements is cut/pasted/modified from Cleveland,
Grosse, and Shyu.

<tt>.data</tt>: Mandatory. Your input data set.

	<tt>.lo_s.model.span</tt>:	smoothing parameter. Default is 0.75.

	<tt>.lo_s.model.degree</tt>: overall degree of locally-fitted polynomial. 1 is
			locally-linear fitting and 2 is locally-quadratic fitting. Default is 2.

	<tt>.lo_s.normalize</tt>:	Should numeric predictors
			be normalized?	If 'y' - the default - the standard normalization
			is used. If 'n', no normalization is carried out.

	\c .lo_s.model.parametric:	for two or more numeric predictors, this argument
			specifies those variables that should be
			conditionally-parametric. The argument should be a logical
			vector of length p, specified in the order of the predictor
			group ordered in x.  Default is a vector of 0's of length p.

	\c .lo_s.model.drop_square:	for cases with degree = 2, and with two or more
			numeric predictors, this argument specifies those numeric
			predictors whose squares should be dropped from the set of
			fitting variables. The method of specification is the same as
			for parametric.  Default is a vector of 0's of length p.

	\c .lo_s.model.family: the assumed distribution of the errors. The values are
	        <tt>"gaussian"</tt> or <tt>"symmetric"</tt>. The first value is the default.
            If the second value is specified, a robust fitting procedure is used.

	\c lo_s.control.surface:	determines whether the fitted surface is computed
            <tt>"directly"</tt> at all points  or whether an <tt>"interpolation"</tt>
            method is used. The default, interpolation, is what most users should use
			unless special circumstances warrant.

    \c lo_s.control.statistics:	determines whether the statistical quantities are 
        computed <tt>"exactly"</tt> or approximately, where <tt>"approximate"</tt>
        is the default. The former should only be used for testing the approximation in
        statistical development and is not meant for routine usage because computation
        time can be horrendous.

        \c lo_s.control.cell: if interpolation is used to compute the surface,
        this argument specifies the maximum cell size of the k-d tree. Suppose k =
        floor(n*cell*span) where n is the number of observations.  Then a cell is
        further divided if the number of observations within it is greater than or
        equal to k. default=0.2

	\c lo_s.control.trace_hat: Options are <tt>"approximate"</tt>, <tt>"exact"</tt>, and <tt>"wait.to.decide"</tt>.	
        When lo_s.control.surface is <tt>"approximate"</tt>, determines
        the computational method used to compute the trace of the hat
        matrix, which is used in the computation of the statistical
        quantities.  If "exact", an exact computation is done; normally
        this goes quite fast on the fastest machines until n, the number
        of observations is 1000 or more, but for very slow machines,
        things can slow down at n = 300.  If "wait.to.decide" is selected,
        then a default is chosen in loess();  the default is "exact" for
        n < 500 and "approximate" otherwise.  If surface is "exact", an
        exact computation is always done for the trace. Set trace_hat to
        "approximate" for large dataset will substantially reduce the
        computation time.

	\c lo_s.model.iterations:	if family is <tt>"symmetric"</tt>, the number of iterations 
        of the robust fitting method.  Default is 0 for
        lo_s.model.family = gaussian; 4 for family=symmetric.

        That's all you can set. Here are some output parameters:

	\c fitted_values:	fitted values of the local regression model

	\c fitted_residuals:	residuals of the local regression fit

       \c  enp:		equivalent number of parameters.

       \c  s:		estimate of the scale of the residuals.

       \c  one_delta:	a statistical parameter used in the computation of standard errors.

       \c  two_delta:	a statistical parameter used in the computation of standard errors.

       \c  pseudovalues:	adjusted values of the response when robust estimation is used.

	\c trace_hat:	trace of the operator hat matrix.

       \c  diagonal:	diagonal of the operator hat matrix.

       \c  robust:		robustness weights for robust fitting.

       \c  divisor:	normalization divisor for numeric predictors.

    \ingroup settings
*/
typedef struct {
    apop_data *data;
    struct  loess_struct lo_s;
    int     want_predict_ci; /**< If 'y' (the default), calculate the
                                confidence bands for predicted values */
    double  ci_level; /**< If running a prediction, the level at which
                        to calculate the confidence interval. default: 0.95 */
} apop_loess_settings;


    /** \cond doxy_ignore */
typedef struct point {    /* a point in the x,y plane */
  double x,y;             /* x and y coordinates */
  double ey;              /* exp(y-ymax+YCEIL) */
  double cum;             /* integral up to x of rejection envelope */
  int f;                  /* is y an evaluated point of log-density */
  struct point *pl,*pr;   /* envelope points to left and right of x */
} POINT;

/* This includes the envelope info and the metropolis steps. */
typedef struct {  /* attributes of the entire rejection envelope */
  int cpoint;              /* number of POINTs in current envelope */
  int npoint;              /* max number of POINTs allowed in envelope */
  double ymax;             /* the maximum y-value in the current envelope */
  POINT *p;                /* start of storage of envelope POINTs */
  double *convex;          /* adjustment for convexity */
  double metro_xprev;      /* previous Markov chain iterate */
  double metro_yprev;      /* current log density at xprev */
} arms_state;
    /** \endcond */

/** to perform derivative-free adaptive rejection sampling with metropolis step */
typedef struct {
    double *xinit;  /**< A <tt>double*</tt> giving starting values for x in ascending order. Default: -1, 0, 1. If this isn't \c NULL, I need at least three items. */
    double  xl;     /**< Left bound. If you don't give me one, I'll use min[min(xinit)/10, min(xinit)*10].*/
    double  xr;     /**< Right bound. If you don't give me one, I'll use max[max(xinit)/10, max(xinit)*10]. */
    double convex;  /**< Adjustment for convexity */
    int ninit;      /**< Number of starting values supplied (i.e. number of elements in \c xinit)*/
    int npoint;     /**< Maximum number of envelope points. I \c malloc space for this many <tt>double</tt>s at the outset. Default = 1e5. */
   char do_metro;   /**< Whether metropolis step is required. (I.e., set to one if you're not sure if the function is log-concave). Set  to <tt>'y'</tt>es or <tt>'n'</tt>o*/
   double xprev;    /**< Previous value from Markov chain */
   int neval;       /**< On exit, the number of function evaluations performed */
   arms_state *state;
   apop_model *model; /**< The model from which I will draw. Mandatory. Must have either a \c log_likelihood or \c p method.*/
} apop_arms_settings;


typedef struct {
    char *splitpage;    /**< The name of the page at which to split the data. If \c NULL, I send the entire data set to both models as needed. */
    apop_model *model1; /**< The first model in the stack.*/
    apop_model *model2; /**< The second model.*/
} apop_stack_settings;

typedef struct {
    apop_data *(*base_to_transformed)(apop_data*);
    apop_data *(*transformed_to_base)(apop_data*);
    double (*jacobian_to_base)(apop_data*);
    apop_model *base_model;
} apop_ct_settings;/**< All of the elements of this struct should be considered private.*/

/** For use with the \ref apop_dconstrain model. See its documentation for an example. 
\hideinitializer \ingroup settings */
typedef struct {
    apop_model *base_model; /**< The model, before constraint. */
    double (*constraint)(apop_data *, apop_model *); /**< The constraint. Return 1 if the data is in the constraint; zero if out. */
    double (*scaling)(apop_model *); /**< Optional. Return the percent of the model density inside the constraint. */
    gsl_rng *rng; /**< If you don't provide a \c scaling function, I calculate the in-constraint model density via random draws.
                       If no \c rng is provided, I use a default RNG; see \ref apop_rng_get_thread. */
    int draw_ct; /**< How many draws to make for calculating the in-constraint model density via random draws. Current default: 1e4. */
} apop_dconstrain_settings;

typedef struct {
    apop_model *generator_m;
    apop_model *ll_m;
    int draw_ct;
} apop_composition_settings;/**< All of the elements of this struct should be considered private.*/

/** For mixture distributions, typically set up using \ref apop_model_mixture. See
\ref apop_mixture for discussion. Please consider all elements but \c model_list and \c
weights as private and subject to change. See the examples for use of these elements.  
\hideinitializer \ingroup settings */
typedef struct {
    gsl_vector *weights;     /**< The likelihood of a draw from each component. */
    apop_model **model_list; /**< A \c NULL-terminated list of component models. */
    int model_count;
    int *param_sizes;  /**< The number of parameters for each model. Useful for unpacking the params. */
    apop_model *cmf;   /**< For internal use by the draw method. */
    int *cmf_refct;    /**< For internal use, so I can garbage-collect the CMF when needed. */
} apop_mixture_settings;

    //Models built via call to apop_model_copy_set.

#define apop_model_coordinate_transform(...) Apop_model_copy_set(apop_coordinate_transform, apop_ct, __VA_ARGS__)
#define apop_model_dcompose(...) Apop_model_copy_set(apop_composition, apop_composition, __VA_ARGS__)
#define apop_model_dconstrain(...) Apop_model_copy_set(apop_dconstrain, apop_dconstrain, __VA_ARGS__)

/** \defgroup settings Settings*/

//Doxygen drops whatever is after these declarations, so I put them last.
Apop_settings_declarations(apop_ct)
Apop_settings_declarations(apop_lm)
Apop_settings_declarations(apop_pm)
Apop_settings_declarations(apop_pmf)
Apop_settings_declarations(apop_mle)
Apop_settings_declarations(apop_cdf)
Apop_settings_declarations(apop_arms)
Apop_settings_declarations(apop_mcmc)
Apop_settings_declarations(apop_loess)
Apop_settings_declarations(apop_stack)
Apop_settings_declarations(apop_mixture)
Apop_settings_declarations(apop_dconstrain)
Apop_settings_declarations(apop_composition)
Apop_settings_declarations(apop_parts_wanted)
Apop_settings_declarations(apop_kernel_density)

#ifdef	__cplusplus
}
#endif

//Part of the intent of a convenience header like this is that you
//don't have to remember what else you're including. So here are 
//some other common GSL headers:
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>
