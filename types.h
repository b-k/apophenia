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

/** This structure holds the names of the components of the \ref apop_data set. You may never have to worry about it directly, because most operation on \ref apop_data sets will take care of the names for you.
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

/** The \ref apop_data structure adds a touch of metadata on top of the basic \c gsl_matrix and \c gsl_vector. It includes an \ref apop_name structure, and a table for non-numeric variables.  Allocate using \c apop_data_alloc, free via \c apop_data_free, or more generally, see the \c apop_data_... section of the index (in the header links) for the many other functions that operate on this struct.
\ingroup data_struct
*/
typedef struct _apop_data apop_data;

struct _apop_data{
    gsl_vector  *vector;
    gsl_matrix  *matrix;
    apop_name   *names;
    char        ***text;
    int         textsize[2];
    gsl_vector  *weights;
    apop_data   *more;
};

/** A description of a parametrized statistical model, including the input settings and the output parameters, predicted/expected values, et cetera.  The full declaration is given in the \c _apop_model page, see the longer discussion on the \ref models page, or see the \ref apop_ols page for a sample program that uses an \ref apop_model.
*/
typedef struct _apop_model apop_model;

typedef struct {
    char name[101];
    void *setting_group;
    void *copy;
    void *free;
} apop_settings_type;

/** The elements of the \ref apop_model type. */
struct _apop_model{
    char        name[101]; 
    int         vbase, m1base, m2base, dsize; /**< The size of the parameter set.
                     If a dimension is -1, then use yourdata->matrix->size2. For
                    anything more complex, allocate the parameter set in the prep
                    method. \c dsize is for the canonical form, and is
                    the size of the data the RNG will return. */
    apop_settings_type *settings;
    apop_data   *parameters; /**< The vector of coefficients or parameters estimated by the model. */
    apop_model **param_dist;
    int         prepared;
    apop_data   *data;
    apop_model * (*estimate)(apop_data * data, apop_model *params);
    double  (*p)(apop_data *d, apop_model *params);
    double  (*log_likelihood)(apop_data *d, apop_model *params);
    void    (*score)(apop_data *d, gsl_vector *gradient, apop_model *params);
    double  (*constraint)(apop_data *data, apop_model *params);
    apop_data*  (*predict)(apop_data *d, apop_model *params);
    double  (*cdf)(apop_data *d, apop_model *params);
    void (*draw)(double *out, gsl_rng* r, apop_model *params);
    void (*prep)(apop_data *data, apop_model *params);
    void (*print)(apop_model *params);
    void    *more;
    size_t  more_size;
} ;

/** The global options.
  \ingroup global_vars */
typedef struct{
    int verbose; /**< Set this to zero for silent mode, one for errors and warnings. default = 0. */
    char output_type;
           /**< 's'   = to screen
                'f'   = to file
                'd'   = to db. 
                'p'   = to pipe (specifically, apop_opts.output_pipe). 
             If 1 or 2, then you'll need to set output_name in the apop_..._print fn. default = 0. */
    FILE *output_pipe; /**< If printing to a pipe or FILE, set it here.  */
    char output_delimiter[100]; /**< The separator between elements of output tables. The default is "\t", but 
                                for LaTeX, use "&\t", or use "|" to get pipe-delimited output. */
    int output_append; /**< Append to output files(1), or overwrite(0)?  default = 0 */
    char input_delimiters[100]; /**< What other people have put between your columns. Default = "|,\t" */
    char db_name_column[300]; /**< If set, the name of the column in your tables that holds row names. */
    char db_nan[100]; /**< The string that the database takes to indicate NaN. May be a regex. */
    char db_engine; /**< If this is 'm', use mySQL, else use SQLite. */
    char db_user[101]; /**< Username for database login. Max 100 chars.  */
    char db_pass[101]; /**< Password for database login. Max 100 chars.  */
    int  thread_count; /**< Threads to use internally. See \ref apop_matrix_apply and family.  */
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
int  apop_name_find(const apop_name *n, const char *findme, const char type);

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
APOP_VAR_DECLARE double * apop_data_ptr(apop_data *data, const int row, const int col, const char *rowname, const char *colname, const char *page);
APOP_VAR_DECLARE double apop_data_get(const apop_data *data, const size_t row, const int  col, const char *rowname, const char *colname, const char *page);
APOP_VAR_DECLARE void apop_data_set(apop_data *data, const size_t row, const int col, const double val, const char *rowname, const char * colname, const char *page);
void apop_data_add_named_elmt(apop_data *d, char *name, double val);
void apop_text_add(apop_data *in, const size_t row, const size_t col, const char *fmt, ...);
apop_data * apop_text_alloc(apop_data *in, const size_t row, const size_t col);
void apop_text_free(char ***freeme, int rows, int cols);
apop_data *apop_data_transpose(apop_data *in);
gsl_matrix * apop_matrix_realloc(gsl_matrix *m, size_t newheight, size_t newwidth);
gsl_vector * apop_vector_realloc(gsl_vector *v, size_t newheight);

APOP_VAR_DECLARE apop_data * apop_data_get_page(const apop_data * data, const char * title);
apop_data * apop_data_add_page(apop_data * dataset, apop_data *newpage,const char *title);
APOP_VAR_DECLARE void apop_data_rm_page(apop_data * data, const char *title, const char free_p);

//Deprecated.
double * apop_data_ptr_it(apop_data *in, const size_t row, const char* col);
double * apop_data_ptr_ti(apop_data *in,const  char* row,const  int col);
double * apop_data_ptr_tt(apop_data *in,const  char *row,const  char* col);
double apop_data_get_it(const apop_data *in,const  size_t row,const  char* col);
double apop_data_get_ti(const apop_data *in,const  char* row,const  int col);
double apop_data_get_tt(const apop_data *in,const  char *row,const  char* col);
void apop_data_set_ti(apop_data *in,const  char* row,const  int col,const  double data);
void apop_data_set_it(apop_data *in,const  size_t row,const  char* col,const  double data);
void apop_data_set_tt(apop_data *in,const  char *row,const  char* col,const  double data);
#ifdef	__cplusplus
}
#endif
#endif
