/** \file  */
/* Copyright (c) 2005--2014 by Ben Klemens.  Licensed under the GPLv2; see COPYING. */

/* Here are the headers for all of apophenia's functions, typedefs, static variables and
macros. All of these begin with the apop_ (or Apop_ or APOP_) prefix.

There used to be a series of sub-headers, but they never provided any serious
benefit. Please use your text editor's word-search feature to find any elements you
may be looking for. About a third of the file is comments and doxygen documentation,
so syntax highlighting that distinguishes code from comments will also help to make
this more navigable.*/

/** \defgroup all_public Public functions, structs, and types
\addtogroup all_public
@{

*/
#pragma once
#ifdef	__cplusplus
extern "C" {
#endif

/** \cond doxy_ignore */
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

These macros are intended as internal. Grep docs/documentation.h for optionaldetails
to find notes on how these are used (Doxygen doesn't use that page), if you are
interested in using this mechanism in out-of-Apophenia work.

*/

#define apop_varad_head(type, name) type variadic_##name(variadic_type_##name varad_in)

#define apop_varad_declare(type, name, ...) \
    typedef struct {                        \
                __VA_ARGS__ ;               \
            } variadic_type_##name;         \
    apop_varad_head(type, name);

#define apop_varad_var(name, value) name = varad_in.name ? varad_in.name : (value);
#define apop_varad_link(name,...) variadic_##name((variadic_type_##name) {__VA_ARGS__})

/** \endcond */ //End of Doxygen ignore.


            //////The types


/** This structure holds the names of the components of the \ref apop_data set. You may never have to worry about it directly, because most operations on \ref apop_data sets will take care of the names for you.
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

Here is a diagram showing a sample data set with all of the elements in place. Together, they represet a data set where each row is an observation, which includes both numeric and text values, and where each row/column is named.

\htmlinclude apop_data_fig.html
\latexinclude apop_data_fig.tex

Allocate using \c apop_data_alloc, free via \c apop_data_free, or more generally, see the \c apop_data_... section of the index (in the header links) for the many other functions that operate on this struct.

See also \ref dataoverview for further notes on getting and manipulating the elements of an \ref apop_data set.
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

/** A statistical model. See \ref modelsec for details. */
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

/** The global options. */
typedef struct{
    int verbose; /**< Set this to zero for silent mode, one for errors and warnings. default = 0. */
    char stop_on_warning; /**< See \ref debugging . */
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
Apop_var_declare( void  apop_name_stack(apop_name * n1, apop_name *nadd, char type1, char typeadd) )
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
But if you forget the \c NULL marker, this has good odds of segfaulting. You may prefer to use a \c for loop that inserts each name in turn using \ref apop_name_add.

\see \ref apop_name_add, although \ref apop_data_add_names will be more useful in most cases. 
*/
#define apop_data_add_names(dataset, type, ...) apop_data_add_names_base((dataset), (type), (char const*[]) {__VA_ARGS__, NULL}) 


/** Free an \ref apop_data structure.
 
\li As with \c free(), it is safe to send in a \c NULL pointer (in which case the function does nothing).
\li If the \c more pointer is not \c NULL, I will free the pointed-to data set first.
If you don't want to free data sets down the chain, set <tt>more=NULL</tt> before calling this.
\li This is actually a macro (that calls \ref apop_data_free_base to do the real work). It
sets \c freeme to \c NULL when it's done, because there's nothing safe you can do with the
freed location, and you can later safely test conditions like <tt>if (data) ...</tt>.
*/
#define apop_data_free(freeme) (apop_data_free_base(freeme) ? 0 : ((freeme)= NULL))

char        apop_data_free_base(apop_data *freeme);
Apop_var_declare( apop_data * apop_data_alloc(const size_t size1, const size_t size2, const int size3) )
Apop_var_declare( apop_data * apop_data_calloc(const size_t size1, const size_t size2, const int size3) )
Apop_var_declare( apop_data * apop_data_stack(apop_data *m1, apop_data * m2, char posn, char inplace) )
apop_data ** apop_data_split(apop_data *in, int splitpoint, char r_or_c);
apop_data * apop_data_copy(const apop_data *in);
void        apop_data_rm_columns(apop_data *d, int *drop);
void apop_data_memcpy(apop_data *out, const apop_data *in);
Apop_var_declare( double * apop_data_ptr(apop_data *data, int row, int col, const char *rowname, const char *colname, const char *page) )
Apop_var_declare( double apop_data_get(const apop_data *data, size_t row, int  col, const char *rowname, const char *colname, const char *page) )
Apop_var_declare( int apop_data_set(apop_data *data, size_t row, int col, const double val, const char *rowname, const char * colname, const char *page) )
void apop_data_add_named_elmt(apop_data *d, char *name, double val);
int apop_text_add(apop_data *in, const size_t row, const size_t col, const char *fmt, ...);
apop_data * apop_text_alloc(apop_data *in, const size_t row, const size_t col);
void apop_text_free(char ***freeme, int rows, int cols);
Apop_var_declare( apop_data * apop_data_transpose(apop_data *in, char transpose_text, char inplace) )
gsl_matrix * apop_matrix_realloc(gsl_matrix *m, size_t newheight, size_t newwidth);
gsl_vector * apop_vector_realloc(gsl_vector *v, size_t newheight);

#define apop_data_prune_columns(in, ...) apop_data_prune_columns_base((in), (char *[]) {__VA_ARGS__, NULL})
apop_data* apop_data_prune_columns_base(apop_data *d, char **colnames);

Apop_var_declare( apop_data * apop_data_get_page(const apop_data * data, const char * title, const char match) )
apop_data * apop_data_add_page(apop_data * dataset, apop_data *newpage,const char *title);
Apop_var_declare( apop_data* apop_data_rm_page(apop_data * data, const char *title, const char free_p) )
Apop_var_declare( apop_data * apop_data_rm_rows(apop_data *in, int *drop, int (*do_drop)(apop_data* ! void*), void* drop_parameter) )

//in apop_asst.c:
Apop_var_declare( apop_data * apop_model_draws(apop_model *model, int count, apop_data *draws) )


/* Convenience functions to convert among vectors (gsl_vector), matrices (gsl_matrix), 
  arrays (double **), and database tables */

//From vector
gsl_vector *apop_vector_copy(const gsl_vector *in);
Apop_var_declare( gsl_matrix * apop_vector_to_matrix(const gsl_vector *in, char row_col) )

//From matrix
gsl_matrix *apop_matrix_copy(const gsl_matrix *in);
apop_data  *apop_db_to_crosstab(char *tabname, char *r1, char *r2, char *datacol);

//From array
Apop_var_declare( gsl_vector * apop_array_to_vector(double *in, int size) )
#define apop_line_to_vector apop_array_to_vector

//From text
Apop_var_declare( apop_data * apop_text_to_data(char const *text_file, int has_row_names, int has_col_names, int const *field_ends, char const *delimiters) )
Apop_var_declare( int apop_text_to_db(char const *text_file, char *tabname, int has_row_names, int has_col_names, char **field_names, int const *field_ends, apop_data *field_params, char *table_params, char const *delimiters, char if_table_exists) )

//rank data
apop_data *apop_data_rank_expand (apop_data *in);
apop_data *apop_data_rank_compress (apop_data *in);

//From crosstabs
void apop_crosstab_to_db(apop_data *in, char *tabname, char *row_col_name, 
						char *col_col_name, char *data_col_name);

//packing data into a vector
Apop_var_declare( gsl_vector * apop_data_pack(const apop_data *in, gsl_vector *out, char more_pages, char use_info_pages) )
Apop_var_declare( void apop_data_unpack(const gsl_vector *in, apop_data *d, char use_info_pages) )

#define apop_vector_fill(avfin, ...) apop_vector_fill_base((avfin), (double []) {__VA_ARGS__})
#define apop_data_fill(adfin, ...) apop_data_fill_base((adfin), (double []) {__VA_ARGS__})
#define apop_text_fill(dataset, ...)   apop_text_fill_base((dataset), (char* []) {__VA_ARGS__, NULL})

#define apop_data_falloc(sizes, ...) apop_data_fill(apop_data_alloc sizes, __VA_ARGS__)
    
apop_data *apop_data_fill_base(apop_data *in, double []);
gsl_vector *apop_vector_fill_base(gsl_vector *in, double []);
apop_data *apop_text_fill_base(apop_data *data, char* text[]);

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
extern apop_model *apop_cross;

/** Alias for the \ref apop_normal distribution, qv. */
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

//transform/apop_cross.c.
apop_model *apop_model_cross_base(apop_model *mlist[]);
#define apop_model_cross(...) apop_model_cross_base((apop_model *[]){__VA_ARGS__, NULL})

    //The variadic versions, with lots of options to input extra parameters to the
    //function being mapped/applied
Apop_var_declare( apop_data * apop_map(apop_data *in, double (*fn_d)(double), double (*fn_v)(gsl_vector*),
                double (*fn_r)(apop_data *), double (*fn_dp)(double! void *), double (*fn_vp)(gsl_vector*! void *),
                double (*fn_rp)(apop_data *! void *), double (*fn_dpi)(double! void *! int),
                double (*fn_vpi)(gsl_vector*! void *! int), double (*fn_rpi)(apop_data*! void *! int),
                double (*fn_di)(double! int), double (*fn_vi)(gsl_vector*! int), double (*fn_ri)(apop_data*! int),
                void *param, int inplace, char part, int all_pages) )
Apop_var_declare( double apop_map_sum(apop_data *in, double (*fn_d)(double), double (*fn_v)(gsl_vector*),
                double (*fn_r)(apop_data *), double (*fn_dp)(double! void *), double (*fn_vp)(gsl_vector*! void *),
                double (*fn_rp)(apop_data *! void *), double (*fn_dpi)(double! void *! int),
                double (*fn_vpi)(gsl_vector*! void *! int), double (*fn_rpi)(apop_data*! void *! int),
                double (*fn_di)(double! int), double (*fn_vi)(gsl_vector*! int), double (*fn_ri)(apop_data*! int),
                void *param, char part, int all_pages) )

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

Apop_var_declare( void apop_plot_histogram(gsl_vector *data, size_t bin_count, char *with, char const *output_name, FILE *output_pipe, char output_type, char output_append) )

Apop_var_declare( void apop_matrix_print(const gsl_matrix *data, char const *output_name, FILE *output_pipe, char output_type, char output_append) )
Apop_var_declare( void apop_vector_print(gsl_vector *data, char const *output_name, FILE *output_pipe, char output_type, char output_append) )
Apop_var_declare( void apop_data_print(const apop_data *data, char const *output_name, FILE *output_pipe, char output_type, char output_append) )

void apop_matrix_show(const gsl_matrix *data);
void apop_vector_show(const gsl_vector *data);
void apop_data_show(const apop_data *data);


        //statistics
Apop_var_declare( double apop_vector_mean(gsl_vector const *v, gsl_vector const *weights))
Apop_var_declare( double apop_vector_var(gsl_vector const *v, gsl_vector const *weights))
Apop_var_declare( double apop_vector_skew_pop(gsl_vector const *v, gsl_vector const *weights))
Apop_var_declare( double apop_vector_kurtosis_pop(gsl_vector const *v, gsl_vector const *weights))
Apop_var_declare( double apop_vector_cov(gsl_vector const *v1, gsl_vector const *v2,
                                         gsl_vector const *weights))

Apop_var_declare( double apop_vector_distance(const gsl_vector *ina, const gsl_vector *inb, const char metric, const double norm) )

Apop_var_declare( void apop_vector_normalize(gsl_vector *in, gsl_vector **out, const char normalization_type) )
void apop_matrix_normalize(gsl_matrix *data, const char row_or_col, const char normalization);

apop_data * apop_data_covariance(const apop_data *in);
apop_data * apop_data_correlation(const apop_data *in);
long double apop_vector_entropy(gsl_vector *in);
long double apop_matrix_sum(const gsl_matrix *m);
double apop_matrix_mean(const gsl_matrix *data);
void apop_matrix_mean_and_var(const gsl_matrix *data, double *mean, double *var);
apop_data * apop_data_summarize(apop_data *data);

apop_data *apop_test_fisher_exact(apop_data *intab); //in apop_fisher.c

//from apop_t_f_chi.c:
Apop_var_declare( int apop_matrix_is_positive_semidefinite(gsl_matrix *m, char semi) )
double apop_matrix_to_positive_semidefinite(gsl_matrix *m);
long double apop_multivariate_gamma(double a, int p);
long double apop_multivariate_lngamma(double a, int p);

//apop_tests.c
apop_data *	apop_t_test(gsl_vector *a, gsl_vector *b);
apop_data *	apop_paired_t_test(gsl_vector *a, gsl_vector *b);
Apop_var_declare( apop_data* apop_anova(char *table, char *data, char *grouping1, char *grouping2) )
#define apop_ANOVA apop_anova
Apop_var_declare( apop_data * apop_f_test (apop_model *est, apop_data *contrast) )
#define apop_F_test apop_f_test

//from the regression code:
#define apop_estimate_r_squared(in) apop_estimate_coefficient_of_determination(in)

apop_data * apop_text_unique_elements(const apop_data *d, size_t col);
gsl_vector * apop_vector_unique_elements(const gsl_vector *v);
Apop_var_declare( apop_data * apop_data_to_factors(apop_data *data, char intype, int incol, int outcol) )
Apop_var_declare( apop_data * apop_data_get_factor_names(apop_data *data, int col, char type) )

Apop_var_declare( apop_data * apop_data_to_dummies(apop_data *d, int col, char type, int keep_first, char append, char remove) )

Apop_var_declare( long double apop_model_entropy(apop_model *in, int draws) )
Apop_var_declare( double apop_kl_divergence(apop_model *from, apop_model *to, int draw_ct, gsl_rng *rng) )

apop_data *apop_estimate_coefficient_of_determination (apop_model *);
void apop_estimate_parameter_tests (apop_model *est);


//Bootstrapping & RNG
apop_data * apop_jackknife_cov(apop_data *data, apop_model *model);
Apop_var_declare( apop_data * apop_bootstrap_cov(apop_data *data, apop_model *model, gsl_rng* rng, int iterations, char keep_boots, char ignore_nans) )
gsl_rng *apop_rng_alloc(int seed);
double apop_rng_GHgB3(gsl_rng * r, double* a); //in apop_asst.c

#define apop_rng_get_thread(thread_in) apop_rng_get_thread_base(#thread_in[0]=='\0' ? -1: (thread_in+0))
gsl_rng *apop_rng_get_thread_base(int thread);

int apop_arms_draw (double *out, gsl_rng *r, apop_model *m);


    // maximum likelihod estimation related functions

Apop_var_declare( gsl_vector * apop_numerical_gradient(apop_data * data, apop_model* model, double delta) )
Apop_var_declare( apop_data * apop_model_hessian(apop_data * data, apop_model *model, double delta) )
Apop_var_declare( apop_data * apop_model_numerical_covariance(apop_data * data, apop_model *model, double delta) )

void apop_maximum_likelihood(apop_data * data, apop_model *dist);

Apop_var_declare( apop_model * apop_estimate_restart (apop_model *e, apop_model *copy, char * starting_pt, double boundary) )

//in apop_linear_constraint.c
Apop_var_declare( long double  apop_linear_constraint(gsl_vector *beta, apop_data * constraint, double margin) )

//in apop_model_fix_params.c
apop_model * apop_model_fix_params(apop_model *model_in);
apop_model * apop_model_fix_params_get_base(apop_model *model_in);






            //////vtables
/** \cond doxy_ignore */

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
m4_define(make_vtab_fns, <|m4_dnl
#ifdef Declare_type_checking_fns
void $1_type_check($1_type in){ };
#else
void $1_type_check($1_type in);
#endif
#define $1_vtable_add(fn, ...) $1_type_check(fn), apop_vtable_add("$1", fn, $1_hash(__VA_ARGS__))
#define $1_vtable_get(...) apop_vtable_get("$1", $1_hash(__VA_ARGS__))
#define $1_vtable_drop(...) apop_vtable_drop("$1", $1_hash(__VA_ARGS__))m4_dnl
|>)

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
make_vtab_fns(apop_update)

typedef long double (*apop_entropy_type)(apop_model *model);
#define apop_entropy_hash(m1) ((size_t)(m1)->log_likelihood + 33 * (size_t)((m1)->p) + 27*(size_t)((m1)->draw))
make_vtab_fns(apop_entropy)

typedef void (*apop_score_type)(apop_data *d, gsl_vector *gradient, apop_model *params);
#define apop_score_hash(m1) ((size_t)((m1)->log_likelihood ? (m1)->log_likelihood : (m1)->p))
make_vtab_fns(apop_score)

typedef apop_model* (*apop_parameter_model_type)(apop_data *, apop_model *);
#define apop_parameter_model_hash(m1) ((size_t)((m1)->log_likelihood ? (m1)->log_likelihood : (m1)->p)*33 + (m1)->estimate ? (size_t)(m1)->estimate: 27)
make_vtab_fns(apop_parameter_model)

typedef apop_data * (*apop_predict_type)(apop_data *d, apop_model *params);
#define apop_predict_hash(m1) ((size_t)((m1)->log_likelihood ? (m1)->log_likelihood : (m1)->p)*33 + (m1)->estimate ? (size_t)(m1)->estimate: 27)
make_vtab_fns(apop_predict)

typedef void (*apop_model_print_type)(apop_model *params, FILE *out);
#define apop_model_print_hash(m1) ((m1)->log_likelihood ? (size_t)(m1)->log_likelihood : \
            (m1)->p ? (size_t)(m1)->p*33 : \
            (m1)->estimate ? (size_t)(m1)->estimate*33*33 : \
            (m1)->draw ? (size_t)(m1)->draw*33*27  : \
            (m1)->cdf ? (size_t)(m1)->cdf*27*27  \
            : 27)
make_vtab_fns(apop_model_print)

/** \endcond */ //End of Doxygen ignore.




        //////Asst


double apop_generalized_harmonic(int N, double s) __attribute__ ((__pure__));

apop_data * apop_test_anova_independence(apop_data *d);
#define apop_test_ANOVA_independence(d) apop_test_anova_independence(d)

Apop_var_declare( int apop_regex(const char *string, const char* regex, apop_data **substrings, const char use_case) )

int apop_system(const char *fmt, ...) __attribute__ ((format (printf,1,2)));

//Histograms and PMFs
gsl_vector * apop_vector_moving_average(gsl_vector *, size_t);
apop_data * apop_histograms_test_goodness_of_fit(apop_model *h0, apop_model *h1);
apop_data * apop_test_kolmogorov(apop_model *m1, apop_model *m2);
apop_data *apop_data_pmf_compress(apop_data *in);
Apop_var_declare( apop_data * apop_data_to_bins(apop_data *indata, apop_data *binspec, int bin_count, char close_top_bin) )
Apop_var_declare( apop_model * apop_model_to_pmf(apop_model *model, apop_data *binspec, long int draws, int bin_count, gsl_rng *rng) )

//text conveniences
Apop_var_declare( char* apop_text_paste(apop_data const*strings, char *between, char *before, char *after, char *between_cols, int (*prune)(apop_data* ! int ! int ! void*), void* prune_parameter) )
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

/** Execute an action and print a message to the current \c FILE handle held by <tt>apop_opts.log_file</tt> (default: \c stderr).
 
\param test The expression that, if true, triggers the action.
\param onfail If the assertion fails, do this. E.g., <tt>out->error='x'; return GSL_NAN</tt>. Notice that it is OK to include several lines of semicolon-separated code here, but if you have a lot to do, the most readable option may be <tt>goto outro</tt>, plus an appropriately-labeled section at the end of your function.
\param level Print the warning message only if \ref apop_opts_type "apop_opts.verbose" is greater than or equal to this. Zero usually works, but for minor infractions use one, or for more verbose debugging output use 2.
\param ... The error message in printf form, plus any arguments to be inserted into the printf string. I'll provide the function name and a carriage return.

Some examples:

\code
//the typical case, stopping function execution:
Apop_stopif(isnan(x), return NAN, 0, "x is NAN; failing");

//Mark a flag, go to a cleanup step
Apop_stopif(x < 0, needs_cleanup=1; goto cleanup, 0, "x is %g; cleaning up and exiting.", x);

//Print a diagnostic iff <tt>apop_opts.verbose>=1</tt> and continue
Apop_stopif(x < 0,  , 1, "warning: x is %g.", x);
\endcode

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
/** \cond doxy_ignore */
#define Apop_assert_c(test, returnval, level, ...) \
    Apop_stopif(!(test), return returnval, level, __VA_ARGS__)

#define Apop_assert(test, ...) Apop_assert_c((test), 0, apop_errorlevel, __VA_ARGS__)

//For things that return void. Transitional and deprecated at birth.
#define Apop_assert_n(test, ...) Apop_assert_c((test),  , apop_errorlevel, __VA_ARGS__)
#define Apop_assert_negone(test, ...) Apop_assert_c((test), -1, apop_errorlevel, __VA_ARGS__)
/** \endcond */ //End of Doxygen ignore.

//Missing data
Apop_var_declare( apop_data * apop_data_listwise_delete(apop_data *d, char inplace) )
apop_model * apop_ml_impute(apop_data *d, apop_model* meanvar);
#define apop_ml_imputation(d, m) apop_ml_impute(d, m)

Apop_var_declare(apop_model *apop_model_metropolis(apop_data *d, gsl_rng* rng, apop_model *m))
Apop_var_declare( apop_model * apop_update(apop_data *data, apop_model *prior, apop_model *likelihood, gsl_rng *rng) )

Apop_var_declare( double apop_test(double statistic, char *distribution, double p1, double p2, char tail) )

//Sorting (apop_asst.c)
Apop_var_declare( double * apop_vector_percentiles(gsl_vector *data, char rounding)  )

//apop_sort.c
Apop_var_declare( apop_data *apop_data_sort(apop_data *data, apop_data *sort_order, char asc, char inplace, double *col_order))

//raking
Apop_var_declare( apop_data * apop_rake(char const *margin_table, char * const*var_list, 
                    int var_ct, char const *all_vars, char * const *contrasts, int contrast_ct, 
                    char const *structural_zeros, int max_iterations, double tolerance, 
                    char const *count_col, int run_number, char const *init_table, 
                    char const *init_count_col, double nudge, char const* table_name) )


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
Apop_var_declare( apop_data * apop_dot(const apop_data *d1, const apop_data *d2, char form1, char form2) )
Apop_var_declare( int         apop_vector_bounded(const gsl_vector *in, long double max) )
gsl_matrix * apop_matrix_inverse(const gsl_matrix *in) ;
double      apop_matrix_determinant(const gsl_matrix *in) ;
//apop_data*  apop_sv_decomposition(gsl_matrix *data, int dimensions_we_want);
Apop_var_declare( apop_data *  apop_matrix_pca(gsl_matrix *data, int const dimensions_we_want) )
Apop_var_declare( gsl_vector * apop_vector_stack(gsl_vector *v1, gsl_vector * v2, char inplace) )
Apop_var_declare( gsl_matrix * apop_matrix_stack(gsl_matrix *m1, gsl_matrix * m2, char posn, char inplace) )
gsl_matrix * apop_matrix_rm_columns(gsl_matrix *in, int *drop);

void apop_vector_log(gsl_vector *v);
void apop_vector_log10(gsl_vector *v);
void apop_vector_exp(gsl_vector *v);

/** \cond doxy_ignore */
/** Deprecated. Use \ref Apop_subm. */
#define APOP_SUBMATRIX(m, srow, scol, nrows, ncols, o) gsl_matrix apop_mm_##o = gsl_matrix_submatrix((m), (srow), (scol), (nrows),(ncols)).matrix;\
gsl_matrix * o = &( apop_mm_##o );
#define Apop_submatrix APOP_SUBMATRIX

/** \def Deprecated. Use \ref Apop_rv. */
#define Apop_row_v(m, row, v) Apop_matrix_row((m)->matrix, row, v)

/** Deprecated. Use \ref Apop_cv. */
#define Apop_col_v(m, col, v) gsl_vector apop_vv_##v = ((col) == -1) ? (gsl_vector){} : gsl_matrix_column((m)->matrix, (col)).vector;\
gsl_vector * v = ((col)==-1) ? (m)->vector : &( apop_vv_##v );

/** Deprecated. Use \ref Apop_rs.  */ 
#define Apop_rows(d, rownum, len, outd) apop_data *outd = Apop_rs(d, rownum, len)

/** Deprecated. Use \ref Apop_r.  */ 
#define Apop_row(d, row, outd) Apop_rows(d, row, 1, outd)

/** Deprecated. Use \ref Apop_cs.  */ 
#define Apop_cols(d, colnum, len, outd) apop_data *outd =  Apop_cs(d, colnum, len);
/** \endcond */ //End of Doxygen ignore.

#define Apop_row_tv(m, row, v) gsl_vector apop_vv_##v = gsl_matrix_row((m)->matrix, apop_name_find((m)->names, row, 'r')).vector;\
gsl_vector * v = &( apop_vv_##v );

#define Apop_col_tv(m, col, v) gsl_vector apop_vv_##v = gsl_matrix_column((m)->matrix, apop_name_find((m)->names, col, 'c')).vector;\
gsl_vector * v = &( apop_vv_##v );

#define Apop_row_t(d, rowname, outd) int apop_row_##outd = apop_name_find((d)->names, rowname, 'r'); Apop_rows(d, apop_row_##outd, 1, outd)

#define Apop_col_t(d, colname, outd) int apop_col_##outd = apop_name_find((d)->names, colname, 'c'); Apop_cols(d, apop_col_##outd, 1, outd)

// The above versions relied on gsl_views, which stick to C as of 1989 CE.
// Better to just create the views via designated initializers.

#define Apop_subm(data_to_view, srow, scol, nrows, ncols)(                  \
        (!(data_to_view)                                                   \
            || (data_to_view)->size1 < (srow)+(nrows) || (srow) < 0        \
            || (data_to_view)->size2 < (scol)+(ncols) || (scol) < 0) ? NULL \
        : &(gsl_matrix){.size1=(nrows), .size2=(ncols),                         \
             .tda=(data_to_view)->tda,                                  \
             .data=gsl_matrix_ptr((data_to_view), (srow), (scol))}      \
        )

#define Apop_rv(data_to_view, row) (                                            \
        ((data_to_view) == NULL || (data_to_view)->matrix == NULL               \
            || (data_to_view)->matrix->size1 <= (row) || (row) < 0) ? NULL        \
        : &(gsl_vector){.size=(data_to_view)->matrix->size2,                    \
             .stride=1, .data=gsl_matrix_ptr((data_to_view)->matrix, (row), 0)} \
        )

#define Apop_cv(data_to_view, col) (                                           \
          !(data_to_view) ? NULL                                               \
        : (col)==-1       ? (data_to_view)->vector                             \
        : (!(data_to_view)->matrix                                             \
            || (data_to_view)->matrix->size2 <= (col) || ((int)(col)) < -1) ? NULL    \
        : &(gsl_vector){.size=(data_to_view)->matrix->size1,                   \
             .stride=(data_to_view)->matrix->tda, .data=gsl_matrix_ptr((data_to_view)->matrix, 0, (col))} \
        )

/** \cond doxy_ignore */
/* Not (yet) for public use. */
#define apop_subvector(v, start, len) (                                          \
        ((v) == NULL || (v)->size < ((start)+(len)) || (start) < 0) ? NULL      \
        : &(gsl_vector){.size=(len), .stride=(v)->stride, .data=(v)->data+(start*(v)->stride)})

/* Not (yet) for public use. */
#define apop_mrow(m, row) (                                       \
        ((m) == NULL || (m)->size1 <= (row) || (row) < 0) ? NULL    \
        : &(gsl_matrix){.size1=1, .size2=(m)->size2, \
             .tda=(m)->tda, .data=gsl_matrix_ptr((m), (row), 0)} \
        )
/** \endcond */

#define Apop_rs(d, rownum, len)(                                                      \
        (!(d) || (rownum) < 0) ? NULL                                       \
        : &(apop_data){                                                          \
         .names= ( !((d)->names) ? NULL :                                        \
            &(apop_name){                                                        \
                .title = (d)->names->title,                                      \
                .vector = (d)->names->vector,                                    \
                .col = (d)->names->col,                                          \
                .row = ((d)->names->row && (d)->names->rowct > (rownum)) ? &((d)->names->row[rownum]) : NULL,  \
                .text = (d)->names->text,                                        \
                .colct = (d)->names->colct,                                      \
                .rowct = (d)->names->row ? (GSL_MIN(1, GSL_MAX((d)->names->rowct - (int)(rownum), 0)))      \
                                          : 0,                                   \
                .textct = (d)->names->textct }),                                 \
        .vector= apop_subvector((d->vector), (rownum), (len)),                   \
        .matrix = Apop_subm(((d)->matrix), (rownum), 0,  (len), (d)->matrix?(d)->matrix->size2:0),    \
        .weights =  apop_subvector(((d)->weights), (rownum), (len)),             \
        .textsize[0]=(d)->textsize[0]> (rownum)+(len)-1 ? (len) : 0,                                   \
        .textsize[1]=(d)->textsize[1],                                           \
        .text = (d)->text ? &((d)->text[rownum]) : NULL,                         \
        })

#define Apop_cs(d, colnum, len) ( \
            (!(d)||!(d)->matrix || (d)->matrix->size2 <= (colnum)+(len)-1        \
             ? NULL                                                              \
             : &(apop_data){                                                     \
                .vector= NULL,                                                   \
                .weights= (d)->weights,                                          \
                .matrix = Apop_subm((d)->matrix, 0, colnum, (d)->matrix->size1, (len)),\
                .textsize[0] = 0,                                                \
                .textsize[1] = 0,                                                \
                .text = NULL,                                                    \
                .names= (d)->names ? &(apop_name){                                                         \
                    .title = (d)->names->title,                                      \
                    .vector = NULL,                                                  \
                    .row = (d)->names->row,                                          \
                    .col = ((d)->names->col && (d)->names->colct > colnum) ? &((d)->names->col[colnum]) : NULL,  \
                    .text = NULL,                                                    \
                    .rowct = (d)->names->rowct,                                      \
                    .colct = (d)->names->col ? (GSL_MIN(len, GSL_MAX((d)->names->colct - colnum, 0)))      \
                                              : 0,                                   \
                    .textct = (d)->names->textct } : NULL \
            })

#define Apop_r(d, rownum) Apop_rs(d, rownum, 1)
#define Apop_c(d, col) Apop_cs(d, col, 1)

/** \cond doxy_ignore */
#define APOP_COL Apop_col
#define apop_col Apop_col
#define APOP_COL_T Apop_col_t
#define apop_col_t Apop_col_t
#define APOP_COL_TV Apop_col_tv
#define apop_col_tv Apop_col_tv

#define APOP_ROW Apop_row
#define apop_row Apop_row
#define APOP_COLS Apop_cols
#define apop_cols Apop_cols
#define APOP_COL_V Apop_col_v
#define apop_col_v Apop_col_v
#define APOP_ROW_V Apop_row_v
#define apop_row_v Apop_row_v
#define APOP_ROWS Apop_rows
#define apop_rows Apop_rows
#define Apop_data_row Apop_row   #deprecated
#define APOP_ROW_T Apop_row_t
#define apop_row_t Apop_row_t
#define APOP_ROW_TV Apop_row_tv
#define apop_row_tv Apop_row_tv
/** \endcond */

/** View a single row of a \c gsl_matrix as a \c gsl_vector. This 
 is a convenience macro wrapping \c gsl_matrix_row. 
 
\param m The \c gsl_matrix
\param row The number of the desired row. 
\param v The name of the vector view that will be created.

See \ref apop_vector_correlation for an example of use.

\see Apop_rs, Apop_r, Apop_row_v, Apop_row_tv, Apop_row_t
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
\see Apop_cs, Apop_c, Apop_cv, Apop_col_tv, Apop_col_t
*/
#define Apop_matrix_col(m, col, v) gsl_vector apop_vv_##v = gsl_matrix_column((m), (col)).vector;\
gsl_vector * v = &( apop_vv_##v );

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

/** \def Apop_subm(data_to_view, srow, scol, nrows, ncols)
Generate a view of a submatrix within a \c gsl_matrix. Like \ref Apop_r, et al., the view is an automatically-allocated variable that is lost once the program flow leaves the scope in which it is declared.

 \param data_to_view The root matrix
 \param srow the first row (in the root matrix) of the top of the submatrix
 \param scol the first column (in the root matrix) of the left edge of the submatrix
 \param nrows number of rows in the submatrix
 \param ncols number of columns in the submatrix
 \return An automatically-allocated view of type c \gsl_matrix.
*/

/** \def Apop_row_t(m, row_name, v)
 After this call, \c v will hold an \ref apop_data view of an \ref apop_data set \c m. The view will consist only of the row with name \c row_name.
 Unlike \ref Apop_r, the second argument is a row name, that I'll look up using \ref apop_name_find, and the third is the name of the view to be generated.
\see Apop_rs, Apop_r, Apop_rv, Apop_row_tv, Apop_matrix_row
*/

/** \def Apop_col_t(m, col_name, v)
 After this call, \c v will hold a view of the \ref apop_data set \c m. The view will consist only of a \c gsl_vector view of the column of the \ref apop_data set \c m with name \c col_name.
 Unlike \ref Apop_c, the second argument is a column name, that I'll look up using \ref apop_name_find, and the third is the name of the view to be generated.
\see Apop_cs, Apop_c, Apop_cv, Apop_col_tv, Apop_matrix_col
*/

/** \def Apop_row_tv(m, row_name, v)
 After this call, \c v will hold a \c gsl_vector view of an \ref apop_data set \c m. The view will consist only of the row with name \c row_name.
 Unlike \ref Apop_rv, the second argument is a row name, that I'll look up using \ref apop_name_find, and the third is the name of the view to be generated.
\see Apop_rs, Apop_r, Apop_rv, Apop_row_t, Apop_matrix_row
*/

/** \def Apop_col_tv(m, col_name, v)
After this call, \c v will hold a \c gsl_vector view of the \ref apop_data set \c m.
The view will consist only of the column with name \c col_name.
Unlike \ref Apop_cv, the second argument is a column name, that I'll look up using \ref apop_name_find, and the third is the name of the view to be generated.
\see Apop_cs, Apop_c, Apop_cv, Apop_col_t, Apop_matrix_col
*/

/** \def Apop_cs(d, col, len)
A macro to generate a temporary view of \ref apop_data set \c d, beginning at column \c col and having length \c len. 
It expires as soon as the program leaves the current scope (like with the usual automatically declared vars). 
\see Apop_c, Apop_cv, Apop_col_tv, Apop_col_t, Apop_matrix_col
*/

/** \def Apop_c(d, col)
A macro to generate a temporary one-column view of \ref apop_data set \c d, pulling out only
column \c col. 
After this call, \c outd will be a pointer to this temporary
view, that you can use as you would any \ref apop_data set.
\see Apop_cs, Apop_cv, Apop_col_tv, Apop_col_t, Apop_matrix_col
*/

/** \def Apop_rs(d, row, len)
A macro to generate a temporary view of \ref apop_data set \c d, beginning at row \c row
and having length \c len. 
The view expires as soon as the program leaves the current scope (like with the usual automatically declared vars). 
\see Apop_r, Apop_rv, Apop_row_tv, Apop_row_t, Apop_matrix_row
*/

/** \def Apop_rv(d, row)
A macro to generate a temporary one-row view of the matrix in an \ref apop_data set \c d, pulling out only
row \c row. The view is a \c gsl_vector set.

\code
gsl_vector *v = Apop_rv(your_data, i);

for (int i=0; i< your_data->matrix->size1; i++)
    printf("Σ_%i = %g\n", i, apop_vector_sum(Apop_r(your_data, i)));
\endcode

The view is automatically allocated, and disappears as soon as the program leaves the scope in which it is declared.
\see Apop_rows, Apop_row_v, Apop_row_tv, Apop_row_t, Apop_matrix_row
*/

/** \def Apop_cv(d, col)
A macro to generate a temporary one-column view of the matrix in an \ref apop_data
set \c d, pulling out only column \c col. The view is a \c gsl_vector set.

As usual, column -1 is the vector element of the \ref apop_data set.

\code
gsl_vector *v = Apop_cv(your_data, i);

for (int i=0; i< your_data->matrix->size2; i++)
    printf("Σ_%i = %g\n", i, apop_vector_sum(Apop_c(your_data, i)));
\endcode

The view is automatically allocated, and disappears as soon as the program leaves the
scope in which it is declared.

\see Apop_cs, Apop_c, Apop_col_tv, Apop_col_t, Apop_matrix_col
*/

/** \def Apop_r(d, row)
A macro to generate a temporary one-row view of \ref apop_data set \c d, pulling out only
row \c row. The view is also an \ref apop_data set, with names and other decorations.
\code
apop_data *v = Apop_r(your_data, i);

for (int i=0; i< your_data->matrix->size1; i++)
    apop_data_print(Apop_r(your_data, i));
\endcode

The view is automatically allocated, and disappears as soon as the program leaves the scope in which it is declared.
\see Apop_rs, Apop_row_v, Apop_row_tv, Apop_row_t, Apop_matrix_row
*/

/** \def apop_mean(v)
 Returns the mean of the elements of the vector \c v.
\param v A \ref gsl_vector.
*/



        //////database utilities

Apop_var_declare( int apop_table_exists(char const *name, char remove) )

int apop_db_open(char const *filename);
Apop_var_declare( int apop_db_close(char vacuum) )

int apop_query(const char *q, ...) __attribute__ ((format (printf,1,2)));
gsl_matrix * apop_query_to_matrix(const char * fmt, ...) __attribute__ ((format (printf,1,2)));
apop_data * apop_query_to_text(const char * fmt, ...) __attribute__ ((format (printf,1,2)));
apop_data * apop_query_to_data(const char * fmt, ...) __attribute__ ((format (printf,1,2)));
apop_data * apop_query_to_mixed_data(const char *typelist, const char * fmt, ...) __attribute__ ((format (printf,2,3)));
gsl_vector * apop_query_to_vector(const char * fmt, ...) __attribute__ ((format (printf,1,2)));
double apop_query_to_float(const char * fmt, ...) __attribute__ ((format (printf,1,2)));

int apop_data_to_db(const apop_data *set, const char *tabname, char);


        //////Settings groups

    //Part I: macros and fns for getting/setting settings groups and elements

void * apop_settings_get_grp(apop_model *m, char *type, char fail);
void apop_settings_remove_group(apop_model *m, char *delme);
void apop_settings_copy_group(apop_model *outm, apop_model *inm, char *copyme);
void *apop_settings_group_alloc(apop_model *model, char *type, void *free_fn, void *copy_fn, void *the_group);
apop_model *apop_settings_group_alloc_wm(apop_model *model, char *type, void *free_fn, void *copy_fn, void *the_group);

/** Retrieves a settings group from a model.  See \ref Apop_settings_get
to just pull a single item from within the settings group.

This macro returns NULL if a group of type \c type_settings isn't found attached
to model \c m, so you can easily put it in a conditional like
  \code 
  if (!apop_settings_get_group(m, "apop_ols")) ...
  \endcode

\param m An \ref apop_model
\param type A string giving the type of the settings group you are retrieving. E.g., for an \ref apop_mle_settings group, use \c apop_mle.
\return A void pointer to the desired struct (or \c NULL if not found).
*/
#define Apop_settings_get_group(m, type) apop_settings_get_grp(m, #type, 'c')

/** Removes a settings group from a model's list. 
 
\li  If the so-named group is not found, do nothing.
*/
#define Apop_settings_rm_group(m, type) apop_settings_remove_group(m, #type)

/** Add a settings group. The first two arguments (the model you are
attaching to and the settings group name) are mandatory, and then you
can use the \ref designated syntax to specify default values (if any).
\return A pointer to the newly-prepped group.

See \ref modelsettings or \ref maxipage for examples.

\li If a settings group of the given type is already attached to the model, 
the previous version is removed. Use \ref Apop_settings_get to check whether a group
of the given type is already attached to a model, and \ref Apop_settings_set to modify
an existing group.
*/
#define Apop_settings_add_group(model, type, ...)  \
    apop_settings_group_alloc(model, #type, type ## _settings_free, type ## _settings_copy, type ##_settings_init ((type ## _settings) {__VA_ARGS__}))

/** Copy a model and add a settings group. Useful for models that require a settings group to function. See \ref Apop_settings_add_group.

\return A pointer to the newly-prepped model.
*/
#define apop_model_copy_set(model, type, ...)  \
    apop_settings_group_alloc_wm(apop_model_copy(model), #type, type ## _settings_free, type ## _settings_copy, type ##_settings_init ((type ## _settings) {__VA_ARGS__}))

/** Retrieves a setting from a model.  See \ref Apop_settings_get_group to pull the entire group.

\param model An \ref apop_model.
\param type A string giving the type of the settings group you are retrieving, without the \c _settings ending. E.g., for an \ref apop_mle_settings group, use \c apop_mle.
\param setting The struct element you want to retrieve.
*/
#define Apop_settings_get(model, type, setting)  \
    (((type ## _settings *) apop_settings_get_grp(model, #type, 'f'))->setting)

/** Modifies a single element of a settings group to the given value. 

\li If <tt>model==NULL</tt>, fails silently. 
\li If <tt>model!=NULL</tt> but the given settings group is not found attached to the model, set <tt>model->error='s'</tt>.
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
See \ref settingswriting for details and an example.
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
See \ref settingswriting for details and an example.
*/
#define Apop_settings_copy(name, ...) \
    void * name##_settings_copy(name##_settings *in) {\
        name##_settings *out = malloc(sizeof(name##_settings)); \
        *out = *in; \
        __VA_ARGS__;    \
        return out;     \
    }

/** A convenience macro for declaring the delete function for a new settings group.
See \ref settingswriting for details and an example.
*/
#define Apop_settings_free(name, ...) \
    void name##_settings_free(name##_settings *in) {\
        __VA_ARGS__;    \
        free(in);  \
    }

        //Part II: the details of extant settings groups.


/** The settings for maximum likelihood estimation (including simulated annealing). */
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

/** Settings for least-squares type models */
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
*/
typedef struct {
    //init/copy/free are in apop_mle.c
    char covariance;    /*< If 'y', calculate the covariance matrix. Default 'n'. */
    char predicted;/*< If 'y', calculate the predicted values. This is typically as many
                     items as rows in your data set. Default 'n'. */
    char tests;/*< If 'y', run any hypothesis tests offered by the model's estimation routine. Default 'n'. */
    char info;/*< If 'y', add an info table with elements such as log likelihood or AIC. Default 'n'. */
} apop_parts_wanted_settings;

/** Some CDFs use random draws; some use closed-form models.  */
typedef struct {
    int draws;  /**< For random draw methods, how many draws? Default: 10,000.*/
    gsl_rng *rng; /**< For random draw methods. See \ref apop_rng_get_thread on the default. */
    apop_model *cdf_model; /**< For use by individual models as they see fit. Default=\c NULL. */
    gsl_matrix *draws_made; /**< A store of random draws that I will count up to report the CDF. Need only be generated once, and so stored here. */
    int *draws_refcount; /**< For internal use.*/
} apop_cdf_settings;


/** Settings for getting parameter models (i.e. the distribution of parameter estimates) */
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


/** Settings for the \ref apop_kernel_density model. */
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

struct apop_mcmc_settings;

/** A proposal distribution for \ref apop_mcmc_settings and its accompanying functions and
information.  By default, these will be \ref apop_multivariate_normal models. The \c
step_fn and \c adapt_fn have to be written around the model and your preferences.
For the defaults, the step function recenters the mean of the distribution around the
last accepted proposal, and the adapt function widens the Σ for the Normal if the
accept rate is too low; narrows it if the accept rate is too large.

You may provide an array of proposals. The length of the list of proposals
must match the number of chunks, as per the \c gibbs_chunks setting in the \ref
apop_mcmc_settings group that the array of proposals is a part of. Each proposal must
be initialized to include all elements, and the step and adapt functions probably have
to be written anew for each type of model.

This segment of the interface is in beta. A future revision may make it easier to design new proposals.
*/
typedef struct apop_mcmc_proposal_s {
    apop_model *proposal; /**< The distribution from which test parameters will be
        drawn. After getting the draw using the \c draw method of the proposal, the base
        model's \c parameters element is filled using \ref apop_data_fill.
        If \c NULL, \ref apop_model_metropolis will use a Multivariate Normal with the
        appropriate dimension, mean zero, and covariance matrix I. If not \c NULL, be sure to
        parameterize your model with an initial position. */

    void (*step_fn)(double const *, struct apop_mcmc_proposal_s*, struct apop_mcmc_settings *); /**< Modifies the parameters of the
        proposal distribution given a successful draw. Typically, this function writes the
        drawn data point to the parameter set. If the draw is a scalar, the default
        function sets the 0th element of the model's \c parameter set with the draw
        (works for the \ref apop_normal and other models). If the draw has multiple
        dimensions, they are all copied to the parameter set, which must have the same
        size. */

    int (*adapt_fn)(struct apop_mcmc_proposal_s *ps, struct apop_mcmc_settings *ms); /**< Called
        every step, to adapt the proposal distribution using information to this point in
        the chain. */

    int accept_count, reject_count;  /**< These are about this chunk. The \ref apop_mcmc_settings group
                                       has a total for the aggregate across all chunks. */
} apop_mcmc_proposal_s;

/** Method settings for a model to be put through Bayesian updating. */
typedef struct apop_mcmc_settings {
    apop_data *data;
    long int periods; /**< For how many steps should the MCMC chain run? */
    double burnin; /**< What <em>percentage</em> of the periods should be ignored
                         as initialization. That is, this is a number between zero and one. */
    int histosegments; /**< If outputting a binned PMF, how many segments should it have? */
    double last_ll; /**< If you have already run mcmc, the last log likelihood in the chain.*/
    apop_model *pmf; /**< If you have already run mcmc, I keep a pointer to the model
            so far here. Use \ref apop_model_metropolis_draw to get one more draw.*/
    apop_model *base_model; /**< The model you provided with a \c log_likelihood or
            \c p element (which need not sum to one). You do not have to set this: if it is
            \c NULL on input to \ref apop_model_metropolis, I will fill it in.*/
    apop_mcmc_proposal_s *proposals; /**< The list of proposals. You can probably use
            the default of adaptive multivariate normals. See the \ref apop_mcmc_proposal_s
            struct for details. */
    int proposal_count; /**< The number of proposal sets; see \c gibbs_chunks below. */
    double target_accept_rate; /**< The desired acceptance rate, for use by adaptive proposals. Default: .35 */
    int accept_count;   /**< After calling apop_mcmc, this will have the number of accepted proposals.*/
    int reject_count;   /**< After calling apop_mcmc, this will have the number of rejected proposals.*/
    char gibbs_chunks;  /**< 'a': One step draws and accepts/rejects all parameters as a unit<br>

                             'b': draw in blocks: the vector is a block, the matrix
                                is a separate block, the weights are a separate
                                block, and so on through every page of the model
                                parameters. Each block of parameters is drawn and
                                accepted/rejected as a unit. <br>

                             '1': draw each parameter and accept/reject separately. One
                                MCMC step consists of a set of draws for every
                                parameter.<br> */
    size_t *block_starts; /**< For internal use */
    int block_count, proposal_is_cp; /**< For internal use. */

    char start_at; /**< If \c '1' (the default), start with a first proposal of all
        1s. Even when this is a far-from-useful starting point, MCMC typically does a good
        job of crawling to better spots early in the chain.<br>
    If \c 'p', start at the \c parameters of the \ref apop_model sent in to \ref
    apop_model_metropolis.*/
    void (*base_step_fn)(double const *, struct apop_mcmc_proposal_s*, struct apop_mcmc_settings *); /**< If a \ref apop_mcmc_proposal_s has \c NULL \c step_fn, use this. If you don't want a step function, set this to a do-nothing function. */
    int (*base_adapt_fn)(struct apop_mcmc_proposal_s *ps, struct apop_mcmc_settings *ms); /**< If a \ref apop_mcmc_proposal_s has \c NULL \c adapt_fn, use this.  If you don't want an adapt function, set this to a do-nothing function.*/

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
*/
typedef struct {
    apop_data *data;
    struct  loess_struct lo_s; /**< 

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
*/

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

/** For use with \ref apop_arms_draw, to perform derivative-free adaptive rejection sampling with metropolis step. 

That function generates default values for this if you do not attach one to the
model beforehand, via a form like <tt>apop_model_add_group(your_model, apop_arms,
.model=your_model, .xl=8, .xr =14);</tt>.  The \c model element is mandatory; you'll
get a run-time complaint if you forget it.
*/
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
} apop_cross_settings;

typedef struct {
    apop_data *(*base_to_transformed)(apop_data*);
    apop_data *(*transformed_to_base)(apop_data*);
    double (*jacobian_to_base)(apop_data*);
    apop_model *base_model;
} apop_ct_settings;/**< All of the elements of this struct should be considered private.*/

/** For use with the \ref apop_dconstrain model. See its documentation for an example. 
*/
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
*/
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
Apop_settings_declarations(apop_cross)
Apop_settings_declarations(apop_mixture)
Apop_settings_declarations(apop_dconstrain)
Apop_settings_declarations(apop_composition)
Apop_settings_declarations(apop_parts_wanted)
Apop_settings_declarations(apop_kernel_density)

#ifdef	__cplusplus
}
#endif

/** @} */ //End doxygen's all_public grouping

//Part of the intent of a convenience header like this is that you
//don't have to remember what else you're including. So here are 
//some other common GSL headers:
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>
