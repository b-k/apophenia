//apop.i   Copyright (c) 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.

/* This is the interface file for SWIG, and unfortunately, this file is
an absolute mess. We can't just #include the other header files, because
we need to clarify what is from Apophenia and what is just included from
elsewhere. Also, swig doesn't yet like __attribs__. Thus, this is redundant
with the other header files.

You'll also find a series of %extend elements that link appropriate
functions to objects, so you can use the familiar forms like
yourdata.copy() instead of apop_data_copy(yourdata).
*/

//These guys are all now called via the object.verb form.
%ignore apop_data_add_named_elmt;
%ignore apop_data_copy;
%ignore apop_data_correlation;
%ignore apop_data_covariance;
%ignore apop_data_get;
%ignore apop_data_listwise_delete;
%ignore apop_data_memcpy;
%ignore apop_data_ptr;
%ignore apop_data_rm_columns;
%ignore apop_data_set;
%ignore apop_data_split;
%ignore apop_data_stack;
%ignore apop_data_summarize;
%ignore apop_data_transpose;
%ignore apop_model_clear;
%ignore apop_name_add;
%ignore apop_name_cross_stack;
%ignore apop_name_find;
%ignore apop_name_print;
%ignore gsl_matrix_add;
%ignore gsl_matrix_add_constant;
%ignore gsl_matrix_add_diagonal;
%ignore gsl_matrix_column;
%ignore gsl_matrix_div_elements;
%ignore gsl_matrix_get;
%ignore gsl_matrix_isnull;
%ignore gsl_matrix_max;
%ignore gsl_matrix_min;
%ignore gsl_matrix_mul_elements;
%ignore gsl_matrix_ptr;
%ignore gsl_matrix_row;
%ignore gsl_matrix_scale;
%ignore gsl_matrix_set;
%ignore gsl_matrix_sub;
%ignore gsl_matrix_swap_columns;
%ignore gsl_matrix_swap_rowcol;
%ignore gsl_matrix_swap_rows;
%ignore gsl_matrix_transpose;
%ignore gsl_matrix_transpose_memcpy;
%ignore gsl_vector_add;
%ignore gsl_vector_add_constant;
%ignore gsl_vector_div;
%ignore gsl_vector_get;
%ignore gsl_vector_isnull;
%ignore gsl_vector_mul;
%ignore gsl_vector_ptr;
%ignore gsl_vector_scale;
%ignore gsl_vector_set;
%ignore gsl_vector_set_all;
%ignore gsl_vector_set_basis;
%ignore gsl_vector_set_zero;
%ignore gsl_vector_sub;


%module apop
%{
#include <apop.h>
%}

%include "carrays.i"
%array_class(double, apop_double)
%array_class(int, apop_int)

%pythoncode %{
from fpconst import *

def apop_matrix_col(data, colno):
    vive = data.column(colno)
    vive.thisown = 0
    return vive.vector

def apop_matrix_row(data, colno):
    vive = data.row(colno)
    vive.thisown = 0
    return vive.vector

def apop_col(data, colno):
    vive = data.matrix.column(colno)
    vive.thisown = 0
    return vive.vector

def apop_row(data, colno):
    vive = data.matrix.row(colno)
    vive.thisown = 0
    return vive.vector

%}

%rename(add) __add;
%rename(add_constant) __add_constant;
%rename(add_diagonal) __add_diagonal;
%rename(add_named_elmt) __add_named_elmt;
%rename(clear) __clear;
%rename(column) __column;
%rename(copy) __copy;
%rename(correlation) __correlation;
%rename(covariance) __covariance;
%rename(cross_stack) __cross_stack;
%rename(div) __div;
%rename(div_elements) __div_elements;
%rename(find) __find;
%rename(get) __get;
%rename(isnull) __isnull;
%rename(listwise_delete) __listwise_delete;
%rename(max) __max;
%rename(memcpy) __memcpy;
%rename(min) __min;
%rename(mul_elements) __mul_elements;
%rename(mul) __mul;
%rename(print) __print;
%rename(ptr) __ptr;
%rename(rm_columns) __rm_columns;
%rename(row) __row;
%rename(scale) __scale;
%rename(set_all) __set_all;
%rename(set_basis) __set_basis;
%rename(set) __set;
%rename(set_zero) __set_zero;
%rename(split) __split;
%rename(stack) __stack;
%rename(sub) __sub;
%rename(summarize) __summarize;
%rename(swap_columns) __swap_columns;
%rename(swap_rowcol) __swap_rowcol;
%rename(swap_rows) __swap_rows;
%rename(transpose_memcpy) __transpose_memcpy;
%rename(transpose) __transpose;

/* Variadics: */
%rename(apop_text_to_data) apop_text_to_data_base;
%rename(apop_text_to_db) apop_text_to_db_base;
%rename(apop_data_to_dummies) apop_data_to_dummies_base;
%rename(apop_det_and_inv) apop_det_and_inv_base;
%rename(apop_db_close) apop_db_close_base;
%rename(apop_dot) apop_dot_base;

%extend apop_data {
    apop_data(const size_t v, const size_t m1, const int m2){ return  apop_data_alloc(v, m1, m2); }
    ~apop_data()    {apop_data_free($self);}
    char* __str__() {apop_data_show($self); return " ";}
    void __show()   {apop_data_show($self);}
    double __get(size_t row, int  col)              { return apop_data_get($self, row, col); }
    double __get(size_t row, char*  col)            { return apop_data_get_it($self, row, col); }
    double __get(char* row, int col)                { return apop_data_get_ti($self, row, col); }
    double __get(char* row, char* col)              { return apop_data_get_tt($self, row,  col); }
    double *__ptr(size_t row, int  col)             { return apop_data_ptr($self, row, col); }
    double *__ptr(size_t row, char*  col)           { return apop_data_ptr_it($self, row, col); }
    double *__ptr(char* row, int col)               { return apop_data_ptr_ti($self, row, col); }
    double *__ptr(char* row, char* col)             { return apop_data_ptr_tt($self, row,  col); }
    void __set(size_t row, int  col, double data)   { apop_data_set($self, row, col, data); }
    void __set(size_t row, char*  col, double data) { apop_data_set_it($self, row, col, data); }
    void __set(char* row, int col, double data)     { apop_data_set_ti($self, row, col, data); }
    void __set(char* row, char* col, double data)   { apop_data_set_tt($self, row,  col, data); }
    apop_data * __copy()                            { return apop_data_copy($self); }
    apop_data * __stack(apop_data * m2, char posn)  { return apop_data_stack($self, m2, posn); }
    apop_data** __split(int splitpoint, char r_or_c){ return apop_data_split($self, splitpoint, r_or_c); }
    apop_data*  __rm_columns(int *drop)             { apop_data_rm_columns($self, drop); return $self; }
    void __memcpy(apop_data *out)                   { return apop_data_memcpy(out, $self);}
    void __add_named_elmt(char *name, double val)   { apop_data_add_named_elmt($self, name, val); }
    apop_data* __listwise_delete()                  { return apop_data_listwise_delete($self); }
    apop_data* __transpose()                        { return apop_data_transpose($self);}
    apop_data* __covariance()                       { return apop_data_covariance($self);}
    apop_data* __correlation()                      { return apop_data_correlation($self);}
    apop_data* __summarize()                        { return apop_data_summarize($self);}
};

%extend _apop_model{
    ~_apop_model (){ apop_model_free ($self);}

    void show (){ apop_model_show ($self);}
    char* __str__() {apop_model_show($self); return " ";}

    apop_model * __clear(apop_data * data){
        return apop_model_clear(data, $self);}

    /*void score(apop_data *d, gsl_vector *out, apop_model *m){
        apop_score(d, out, $self);}

    double log_likelihood(apop_data *d){
        double apop_log_likelihood(d, $self);}

    double p(apop_data *d, apop_model *m){ return apop_p(d, $self);}

    void apop_draw(double *out, gsl_rng *r){ apop_draw(out, r, $self);}

    void prep(apop_data *d){ apop_model_prep(d, $self);}

    apop_data * expected_value(apop_data *d){ return apop_expected_value(d, $self);}*/
}


%extend apop_name {
            apop_name()     { return apop_name_alloc(); }
            ~apop_name()    { apop_name_free($self); }
    void    __print()       { apop_name_print($self); }
    apop_name* __copy()     { return apop_name_copy($self);}
    int     __add(char *add_me, char type)      { return apop_name_add($self, add_me, type); }
    void    __stack(apop_name *n2, char type)   { apop_name_stack($self, n2, type); }
    void    __cross_stack(apop_name *n2, char type1, char type2){ apop_name_cross_stack($self, n2, type1, type2); }
    int     __find(char *findme, char type)     { return apop_name_find($self, findme, type); }
};

%extend gsl_matrix {
    gsl_matrix(size_t i, size_t j) {return gsl_matrix_alloc(i, j);}
    ~gsl_matrix() {gsl_matrix_free($self);}
    char*   __str__()                   {apop_matrix_show($self); return " ";}
    double  __get(size_t i, size_t j)   { gsl_matrix_get($self, i, j); }
    void    __set(size_t i, size_t j, double in)  { gsl_matrix_set($self, i, j, in); }
    double* __ptr(size_t i, size_t j)   {return gsl_matrix_ptr($self, i, j);}
    apop_data   *to_data()              { return apop_matrix_to_data($self); }
    gsl_matrix  *copy()                 { return apop_matrix_copy($self);}
    gsl_matrix  *inverse()              { return apop_matrix_inverse($self) ;}
    double      determinant()           { return apop_matrix_determinant($self) ;}
    apop_data   *summarize()            { return apop_matrix_summarize($self);}
    long double sum()                   { return apop_matrix_sum($self) ;}
    double      mean()                  { return apop_matrix_mean($self) ;}

    void normalize(const char row_or_col, const char normalization){
        return apop_matrix_normalize($self, row_or_col, normalization); }

    apop_data*  pca( int dimensions_we_want){
    return  apop_matrix_pca($self, dimensions_we_want);}

    inline void increment( int i, int j, double amt){
        return apop_matrix_increment($self, i, j, amt);}

    gsl_matrix *stack( gsl_matrix * m2, char posn){
        return apop_matrix_stack($self,  m2, posn);}

    gsl_matrix *rm_columns( int *drop){
    return apop_matrix_rm_columns($self, drop);}

    gsl_matrix *covariance(const char normalize){
        return apop_matrix_covariance($self, normalize);}

    gsl_matrix * correlation(const char normalize){
        return apop_matrix_correlation($self, normalize);}

    double  var_m(double mean)                          { return apop_matrix_var_m($self, mean) ;}
    void    mean_and_var (double *mean, double *var)    { apop_matrix_mean_and_var($self, mean, var);}
    _gsl_vector_view __row (const size_t i)             {return gsl_matrix_row ($self, i);}
    _gsl_vector_view __column (const size_t j)          {return gsl_matrix_column ($self, j);}

    double __max () {return gsl_matrix_max ($self);}
    double __min () {return gsl_matrix_min ($self);}
    int __isnull () {return gsl_matrix_isnull ($self);}
    int __swap_rows(const size_t i, const size_t j)     {return gsl_matrix_swap_rows($self, i, j);}
    int __swap_columns(const size_t i, const size_t j)  {return gsl_matrix_swap_columns($self, i, j);}
    int __swap_rowcol(const size_t i, const size_t j)   {return gsl_matrix_swap_rowcol($self, i, j);}
    int __transpose ()                                  {return gsl_matrix_transpose ($self);}
    int __transpose_memcpy (gsl_matrix * dest)          {return gsl_matrix_transpose_memcpy (dest, $self);}
    int __add (const gsl_matrix * b)                    {return gsl_matrix_add ($self,  b);}
    int __sub (const gsl_matrix * b)                    {return gsl_matrix_sub ($self,  b);}
    int __mul_elements (const gsl_matrix * b)           {return gsl_matrix_mul_elements ($self,  b);}
    int __div_elements (const gsl_matrix * b)           {return gsl_matrix_div_elements ($self,  b);}
    int __scale (const double x)                        {return gsl_matrix_scale ($self, x);}
    int __add_constant (const double x)                 {return gsl_matrix_add_constant ($self, x);}
    int __add_diagonal (const double x)                 {return gsl_matrix_add_diagonal ($self, x);}
};

%extend gsl_vector {
    gsl_vector(size_t i) {return gsl_vector_alloc(i);}
    ~gsl_vector() {gsl_vector_free($self);}
    char*   __str__()                   {apop_vector_show($self); return " ";}
    void    __set(size_t i, double in)  { gsl_vector_set($self, i, in); }
    double  __get(size_t i)             { gsl_vector_get($self, i); }
    double* __ptr(size_t i)             { gsl_vector_ptr($self, i);}
    long double sum()                   { return apop_sum($self) ; }
    void        log()                   { apop_vector_log($self);}
    void        log10()                 { apop_vector_log10($self);}
    void        exp()                   { apop_vector_exp($self);}
    double      mean()                  { return apop_vector_mean($self) ; }
    double      var()                   { return apop_vector_var($self) ; }
    double      var_m(const double mean){ return apop_vector_var_m($self, mean) ; }
    apop_data * to_data()               { return apop_vector_to_data($self); }
    double      kurtosis_pop()          { return apop_vector_kurtosis_pop($self) ; }
    double      kurtosis()              { return apop_vector_kurtosis($self) ; }
    double      skew()                  { return apop_vector_skew($self) ; }
    double      skew_pop()              { return apop_vector_skew_pop($self) ; }
    double      kurt()                  { return apop_vector_kurt($self) ; }
    double      cov(const gsl_vector *inb)   { return apop_vector_cov($self, inb) ; }

    double var_m( const double mean) {
        return apop_vector_var_m($self, mean) ; }

    gsl_vector * moving_average(size_t window){
        return  apop_vector_moving_average($self, window);}

    double cov(const gsl_vector *inb) { return apop_vector_cov($self, inb) ; }

    double correlation(const gsl_vector *inb) { return apop_vector_correlation($self, inb) ; }
    gsl_vector * moving_average(size_t window){ return  apop_vector_moving_average($self, window);}

    //double * percentiles(char rounding) {return apop_vector_percentiles($self, rounding);}
    apop_double * percentiles(char rounding) {return apop_vector_percentiles($self, rounding);}
    double weighted_mean( const gsl_vector *weight) {
        return apop_vector_weighted_mean($self, weight);
    }

    double weighted_var( const gsl_vector *weight) {
        return apop_vector_weighted_var($self, weight) ;
    }

    double weighted_cov(const gsl_vector *inb, const gsl_vector *weights) {
        return apop_vector_weighted_cov($self, inb, weights);
    }

    double weighted_skew( const gsl_vector *w) {
        return apop_vector_weighted_skew($self, w);
    }

    double weighted_kurt( const gsl_vector *w) {
        return apop_vector_weighted_kurt($self, w);
    }

    double distance(const gsl_vector *inb){
        return apop_vector_distance($self, inb); }

    double grid_distance(const gsl_vector *inb){
        return apop_vector_grid_distance($self, inb);}

    void normalize(gsl_vector **out, const char normalization_type){
        apop_vector_normalize($self, out, normalization_type); }

    void increment(gsl_vector * v, int i, double amt){
        return apop_vector_increment($self, i, amt);}

    int  bounded(long double max){ return apop_vector_bounded($self, max);}
    gsl_vector *stack(gsl_vector * v2){ return apop_vector_stack($self, v2);}

    void __set_zero () {gsl_vector_set_zero ($self);}
    void __set_all (double x) {return gsl_vector_set_all ($self, x);}
    int __set_basis (size_t i) {return gsl_vector_set_basis ($self, i);}
    int __add (const gsl_vector * b) {return gsl_vector_add ($self,  b);}
    int __sub (const gsl_vector * b) {return gsl_vector_sub ($self,  b);}
    int __mul (const gsl_vector * b) {return gsl_vector_mul ($self,  b);}
    int __div (const gsl_vector * b) {return gsl_vector_div ($self,  b);}
    int __scale (const double x) {return gsl_vector_scale ($self, x);}
    int __add_constant (const double x) {return gsl_vector_add_constant ($self, x);}
    int __isnull () {return gsl_vector_isnull ($self);}
};

/* Now declare everything that will be included */

typedef struct _apop_model apop_model;

typedef struct {
    char name[101];
    void *setting_group;
    void *copy;
    void *free;
} apop_settings_type;

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

typedef struct{
	char * vector;
	char ** column;
	char ** row;
	char ** text;
	int colct, rowct, textct;
    char title[101];
} apop_name;

typedef struct {
    gsl_vector  *vector;
    gsl_matrix  *matrix;
    apop_name   *names;
    char        ***text;
    int         textsize[2];
    gsl_vector  *weights;
} apop_data;

typedef struct{
    int verbose;
    char output_type;
    FILE *output_pipe;
    char output_delimiter[100];
    int output_append;
    char input_delimiters[100];
    char db_name_column[300];
    char db_nan[100];
    char db_engine;
    int  thread_count;
} apop_opts_type;

extern apop_opts_type apop_opts;

void 		apop_model_free (apop_model * free_me);
void 		apop_model_show (apop_model * print_me);

void        apop_data_free(apop_data *freeme);
apop_data * apop_matrix_to_data(gsl_matrix *m);
apop_data * apop_vector_to_data(gsl_vector *v);
apop_data * apop_data_alloc(const size_t, const size_t, const int);
apop_data * apop_data_calloc(const size_t, const size_t, const int);
apop_data * apop_data_stack(apop_data *m1, apop_data * m2, char posn);
apop_data ** apop_data_split(apop_data *in, int splitpoint, char r_or_c);
apop_data * apop_data_copy(const apop_data *in);
void        apop_data_rm_columns(apop_data *d, int *drop);
void apop_data_memcpy(apop_data *out, const apop_data *in);
double * apop_data_ptr(const apop_data *data, const int i, const int j);
double * apop_data_ptr_it(const apop_data *in, size_t row, char* col);
double * apop_data_ptr_ti(const apop_data *in, char* row, int col);
double * apop_data_ptr_tt(const apop_data *in, char *row, char* col);
void apop_text_add(apop_data *in, const size_t row, const size_t col, const char *fmt, ...);
apop_data * apop_text_alloc(apop_data *in, const size_t row, const size_t col);
apop_data *apop_data_transpose(apop_data *in);
gsl_matrix * apop_matrix_realloc(gsl_matrix *m, size_t newheight, size_t newwidth);
gsl_vector * apop_vector_realloc(gsl_vector *v, size_t newheight);
apop_name * apop_name_alloc(void);
//int apop_name_add(apop_name * n, char *add_me, char type);
//void  apop_name_free(apop_name * free_me);
//void  apop_name_print(apop_name * n);
//void  apop_name_stack(apop_name * n1, apop_name *n2, char type);
//void  apop_name_cross_stack(apop_name * n1, apop_name *n2, char type1, char type2);
//apop_name * apop_name_copy(apop_name *in);
//int  apop_name_find(apop_name *n, char *findme, char type);

void apop_text_free(char ***freeme, int rows, int cols); //in apop_data.c

void apop_opts_memcpy(apop_opts_type *out, apop_opts_type *in); //in apop_output.c
double apop_generalized_harmonic(int N, double s);
void apop_error(int level, char stop, char *message, ...);

apop_data * apop_test_anova_independence(apop_data *d);
#define apop_test_ANOVA_independence(d) apop_test_anova_independence(d)

int apop_system(const char *fmt, ...);

gsl_vector * apop_vector_moving_average(gsl_vector *, size_t);
apop_model *apop_histogram_moving_average(apop_model *m, size_t bandwidth);

#define Apop_assert(test, returnval, level, stop, ...) do \
    if (!(test)) {  \
        if (apop_opts.verbose >= level) { fprintf(stderr, "%s: ", __func__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");}   \
        if (stop == 's' || stop == 'h') assert(test);   \
        return returnval;  \
} while (0);

#define apop_assert(test, returnval, level, stop, ...) Apop_assert(test, returnval, level, stop, __VA_ARGS__)
#define APOP_ASSERT(test, returnval, level, stop, ...) Apop_assert(test, returnval, level, stop, __VA_ARGS__)

#define Apop_assert_void(test,  level, stop, ...) do \
    if (!(test)) {  \
        if (apop_opts.verbose >= level) { fprintf(stderr, "%s: ", __func__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");}   \
        if (stop == 's' || stop == 'h') assert(test);   \
} while (0);

#define apop_assert_void(test, level, stop, ...) Apop_assert_void(test, level, stop, __VA_ARGS__)
#define APOP_ASSERT_VOID(test, level, stop, ...) Apop_assert_void(test, level, stop, __VA_ARGS__)

#define Apop_settings_alloc(type, out, ...) apop_ ##type ##_settings *out = apop_ ##type ##_settings_alloc(__VA_ARGS__);

#define APOP_SETTINGS_ALLOC(type, out, ...) Apop_settings_alloc(type, out, __VA_ARGS__)

typedef struct{
    apop_data *data;
    apop_data *starting_pt;
    long int periods;
    double burnin;
    int histosegments;
    char method;
} apop_update_settings;

apop_update_settings *apop_update_settings_alloc(apop_data *d);
#define apop_update_settings_copy NULL
#define  apop_update_settings_free NULL

apop_model * apop_update(apop_data *data, apop_model *prior, apop_model *likelihood, gsl_rng *r);

apop_data * apop_jackknife_cov(apop_data *data, apop_model model);
apop_data * apop_bootstrap_cov(apop_data *data, apop_model model, gsl_rng*, int);
gsl_rng *apop_rng_alloc(int seed);

gsl_vector *apop_vector_copy(const gsl_vector *in);
double * apop_vector_to_array(const gsl_vector *in);
gsl_matrix * apop_vector_to_matrix(const gsl_vector *in);
gsl_matrix *apop_matrix_copy(const gsl_matrix *in);
apop_data  *apop_db_to_crosstab(char *tabname, char *r1, char *r2, char *datacol);
//takes a three-column table (dim1, dim2, data) and creates a 2D crosstab.
//Returns the crosstab, and the dimension names (if d1!=NULL and d2!=NULL).
gsl_vector * apop_array_to_vector(const double *in, const int size);
gsl_matrix * apop_array_to_matrix(const double **in, const int rows, const int cols);
apop_data * apop_array_to_data(const double **in, const int rows, const int cols);
gsl_matrix * apop_line_to_matrix(double *line, int rows, int cols);
apop_data * apop_line_to_data(double *in, int vsize, int rows, int cols);
apop_data * apop_text_to_data_base(char *text_file, int has_row_names, int has_col_names);
int apop_text_to_db_base(char *text_file, char *tabname, int has_row_names, int has_col_names, char **field_names);
int apop_crosstab_to_db(apop_data *in, char *tabname, char *row_col_name, 
						char *col_col_name, char *data_col_name);

gsl_vector * apop_data_pack(const apop_data *in);
void apop_data_unpack(const gsl_vector *in, apop_data *d);

char * apop_strip_dots(char *in, char strip_type);

gsl_vector *apop_vector_vfill(gsl_vector *in, va_list ap);
apop_data *apop_data_fill(apop_data *in, ...);
gsl_vector *apop_vector_fill(gsl_vector *in, ...);
gsl_matrix *apop_matrix_fill(gsl_matrix *in, ...);

extern int asprintf (char **result, const char *format, ...);

int apop_table_exists(char *q, char whattodo);

void apop_db_rng_init(int seed);

int apop_count_cols(const char *name);
int apop_db_open(char *filename);

int apop_db_close_base(char vacuum);

int apop_query(const char *q, ...) ;

gsl_matrix * apop_query_to_matrix(const char * fmt, ...) ;
apop_data * apop_query_to_text(const char * fmt, ...) ;
apop_data * apop_query_to_data(const char * fmt, ...) ;
apop_data * apop_query_to_mixed_data(const char *typelist, const char * fmt, ...) ;
gsl_vector * apop_query_to_vector(const char * fmt, ...) ;
double apop_query_to_float(const char * fmt, ...) ;

int apop_matrix_to_db(gsl_matrix *data,char *tabname, char **headers);
int apop_data_to_db(apop_data *set, char *tabname);

void apop_db_merge(char *infile);
void apop_db_merge_table(char *infile, char *tabname);
double apop_db_t_test(char * tab1, char *col1, char *tab2, char *col2);
double apop_db_paired_t_test(char * tab1, char *col1, char *col2);

apop_model *apop_histogram_refill_with_vector(apop_model *template, gsl_vector *indata);
apop_model *apop_histogram_refill_with_model(apop_model *template, apop_model *m, long int draws, gsl_rng *r);
apop_model *apop_histogram_vector_reset(apop_model *template, gsl_vector *indata);
apop_model *apop_histogram_model_reset(apop_model *template, apop_model *m, long int draws, gsl_rng *r);
apop_data *apop_histograms_test_goodness_of_fit(apop_model *h0, apop_model *h1);
apop_data *apop_test_kolmogorov(apop_model *m1, apop_model *m2);
void apop_histogram_normalize(apop_model *m);

typedef enum {
    APOP_SIMPLEX_NM     =0, 
    APOP_CG_FR     =1,     
    APOP_CG_BFGS   =2,    
    APOP_CG_PR     =3,   
    APOP_SIMAN      =5,     
    APOP_RF_NEWTON  =10,    
    APOP_RF_HYBRID  =12,    
    APOP_RF_HYBRID_NOSCALE  =13 
} apop_optimization_enum;

/** The settings for maximum likelihood estimation (including simulated annealing).*/
typedef struct{
//traditional
    double      *starting_pt;
    double      step_size;
    double      tolerance;
    double      delta;
    apop_optimization_enum method;
    int         verbose;
    int         want_cov;
//simulated annealing (also uses step_size);
    int         n_tries;
    int         use_score;
    int         iters_fixed_T;
    double      k, t_initial, mu_t, t_min ;
    gsl_rng     *rng;
    char        *trace_path;
} apop_mle_settings;

apop_mle_settings *apop_mle_settings_alloc(apop_model *model);
void *apop_mle_settings_copy(apop_mle_settings * in);
void apop_mle_settings_free(void * in);


typedef double 	(*apop_fn_with_params) (apop_data *, apop_model *);
gsl_vector * apop_numerical_gradient(apop_data *data, apop_model*);

//void apop_numerical_covariance_matrix(apop_model dist, apop_model *est, apop_data *data);
//void apop_numerical_var_covar_matrix(apop_model dist, apop_model *est, apop_data *data);


apop_model *	apop_maximum_likelihood(apop_data * data, apop_model dist);

apop_model * apop_estimate_restart (apop_model *, apop_model *);

double  apop_linear_constraint(gsl_vector *beta, apop_data * constraint, double margin);

apop_model *apop_model_fix_params(apop_data *data, apop_data *paramvals, apop_data *mask, apop_model model_in);

double      apop_det_and_inv_base(const gsl_matrix *in, gsl_matrix **out, int calc_det, int calc_inv);
gsl_matrix *apop_matrix_inverse(const gsl_matrix *in) ;
double      apop_matrix_determinant(const gsl_matrix *in) ;
apop_data*  apop_matrix_pca(gsl_matrix *data, int dimensions_we_want);
inline void apop_vector_increment(gsl_vector * v, int i, double amt);
inline void apop_matrix_increment(gsl_matrix * m, int i, int j, double amt);
gsl_vector *apop_vector_stack(gsl_vector *v1, gsl_vector * v2);
gsl_matrix *apop_matrix_stack(gsl_matrix *m1, gsl_matrix * m2, char posn);
gsl_matrix *apop_matrix_rm_columns(gsl_matrix *in, int *drop);
int         apop_vector_bounded(gsl_vector *in, long double max);
apop_data * apop_dot_base(const apop_data *d1, const apop_data *d2, char form1, char form2);
void        apop_vector_log(gsl_vector *v);
void        apop_vector_log10(gsl_vector *v);
void        apop_vector_exp(gsl_vector *v);
gsl_vector *apop_matrix_map(const gsl_matrix *m, double (*fn)(gsl_vector*));
gsl_vector *apop_vector_map(const gsl_vector *v, double (*fn)(double));
void apop_matrix_apply(gsl_matrix *m, void (*fn)(gsl_vector*));
void apop_vector_apply(gsl_vector *v, void (*fn)(double*));
gsl_matrix * apop_matrix_map_all(const gsl_matrix *in, double (*fn)(double));
void apop_matrix_apply_all(gsl_matrix *in, void (*fn)(double *));

double apop_vector_map_sum(const gsl_vector *in, double(*fn)(double));
double apop_matrix_map_sum(const gsl_matrix *in, double (*fn)(gsl_vector*));
double apop_matrix_map_all_sum(const gsl_matrix *in, double (*fn)(double));

apop_data * apop_data_listwise_delete(apop_data *d);
apop_model * apop_ml_imputation(apop_data *d, apop_model* meanvar);

void apop_plot_line_and_scatter(apop_data *data, apop_model *est, char *);
void apop_histogram_plot(apop_model *hist, char *outfile);
void apop_plot_histogram(gsl_vector *data, size_t bin_ct, char *outfile);
void apop_histogram_print(apop_model *h, char *outfile);
void apop_plot_lattice(apop_data *d, char filename[]);
void apop_plot_qq(gsl_vector *v, apop_model m, char *outfile);

void apop_matrix_print(gsl_matrix *data, char *file);
void apop_vector_print(gsl_vector *data, char *file);
void apop_data_print(apop_data *data, char *file);

void apop_matrix_show(const gsl_matrix *data);
void apop_vector_show(const gsl_vector *data);
void apop_data_show(const apop_data *data);

/** Settings for least-squares type models */
typedef struct {
    int destroy_data;
    gsl_vector *weights;
    apop_data *instruments;
    int want_cov;
    int want_expected_value;
    void *copy;
    void *free;
} apop_ls_settings;

apop_data *apop_F_test (apop_model *est, apop_data *contrast);
apop_data *apop_f_test (apop_model *est, apop_data *contrast);

apop_data *	apop_t_test(gsl_vector *a, gsl_vector *b);
apop_data *	apop_paired_t_test(gsl_vector *a, gsl_vector *b);

apop_data * apop_text_unique_elements(const apop_data *d, size_t col);
gsl_vector * apop_vector_unique_elements(const gsl_vector *v);
apop_data *apop_text_to_factors(apop_data *d, size_t textcol, int datacol);
apop_data * apop_data_to_dummies_base(apop_data *d, int col, char type, int keep_first);

double apop_two_tailify(double in);

apop_data *apop_estimate_coefficient_of_determination (apop_model *in);
apop_data *apop_estimate_r_squared (apop_model *in);
void apop_estimate_parameter_t_tests (apop_model *est);

apop_data* apop_anova(char *table, char *data, char *grouping1, char *grouping2);

#define apop_ANOVA(table, data, grouping1, grouping2) apop_anova(table, data, grouping1, grouping2)

void * apop_settings_get_group(apop_model *m, char *type);
void apop_settings_rm_group(apop_model *m, char *delme);
void apop_settings_copy_group(apop_model *outm, apop_model *inm, char *copyme);
void apop_settings_group_alloc(apop_model *model, char *type, void *free_fn, void *copy_fn, void *the_group);

#define Apop_settings_get_group(m, type) apop_settings_get_group(m, #type)
#define Apop_settings_rm_group(m, type) apop_settings_rm_group(m, #type)

#define Apop_settings_add_group(model, type, ...)  \
    apop_settings_group_alloc(model, #type, type ## _settings_free, type ## _settings_copy, type ##_settings_alloc (__VA_ARGS__)); 

#define Apop_settings_alloc_add(model, type, setting, data, ...)  \
    do {                                                \
        Apop_settings_add_group(model, type, __VA_ARGS__)           \
        Apop_settings_add(model, type, setting, data)       \
    } while (0);

#define Apop_settings_get(model, type, setting)  \
    (((type ## _settings *) apop_settings_get_group(model, #type))->setting)

#define Apop_settings_add(model, type, setting, data)  \
    do {                                                \
    apop_assert_void(apop_settings_get_group(model, #type), 0, 's', "You're trying to modify a setting in " \
                        #model "'s setting group of type " #type " but that model doesn't have such a group."); \
    ((type ## _settings *) apop_settings_get_group(model, #type))->setting = (data);    \
    } while (0);

#define APOP_SETTINGS_ADD Apop_settings_add
#define APOP_SETTINGS_ALLOC_ADD Apop_settings_alloc_add
#define APOP_SETTINGS_GET Apop_settings_get
#define APOP_SETTINGS_ADD_GROUP Apop_settings_add_group
#define APOP_SETTINGS_GET_GROUP Apop_settings_get_group
#define APOP_SETTINGS_RM_GROUP Apop_settings_rm_group

#define APOP_SUBMATRIX(m, srow, scol, nrows, ncols, o) gsl_matrix apop_mm_##o = gsl_matrix_submatrix(m, (srow), (scol), (nrows),(ncols)).matrix;\
gsl_matrix * o = &( apop_mm_##o );

#define APOP_MATRIX_ROW(m, row, v) gsl_vector apop_vv_##v = gsl_matrix_row(m, (row)).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_MATRIX_COL(m, col, v) gsl_vector apop_vv_##v = gsl_matrix_column(m, (col)).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_ROW_T(m, row, v) gsl_vector apop_vv_##v = gsl_matrix_row((m)->matrix, apop_name_find((m)->names, row, 'r')).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_COL_T(m, col, v) gsl_vector apop_vv_##v = gsl_matrix_column((m)->matrix, apop_name_find((m)->names, col, 'c')).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_ROW(m, row, v) gsl_vector apop_vv_##v = gsl_matrix_row((m)->matrix, (row)).vector;\
gsl_vector * v = &( apop_vv_##v );

#define APOP_COL(m, col, v) gsl_vector apop_vv_##v = gsl_matrix_column((m)->matrix, (col)).vector;\
gsl_vector * v = &( apop_vv_##v );

#define Apop_col APOP_COL 
#define Apop_row APOP_ROW
#define Apop_col_t APOP_COL_T
#define Apop_row_t APOP_ROW_T
#define Apop_matrix_col APOP_MATRIX_COL 
#define Apop_matrix_row APOP_MATRIX_ROW
#define Apop_submatrix APOP_SUBMATRIX

inline long double apop_vector_sum(const gsl_vector *in) ;
inline long double apop_sum(const gsl_vector *in) ;
inline double apop_vector_mean(const gsl_vector *in) ;
inline double apop_mean(const gsl_vector *in) ;
inline double apop_vector_var(const gsl_vector *in) ;
inline double apop_var(const gsl_vector *in) ;
inline double apop_vector_var_m(const gsl_vector *in, const double mean) ;
inline double apop_vector_cov(const gsl_vector *ina, const gsl_vector *inb) ;
inline double apop_vector_correlation(const gsl_vector *ina, const gsl_vector *inb) ;
inline double apop_vector_kurtosis_pop(const gsl_vector *in) ;
inline double apop_vector_kurtosis(const gsl_vector *in) ;
inline double apop_vector_skew(const gsl_vector *in) ;
inline double apop_vector_skew_pop(const gsl_vector *in) ;
inline double apop_vector_kurt(const gsl_vector *in) ;
double apop_vector_weighted_mean(const gsl_vector *, const gsl_vector *) ;
double apop_vector_weighted_var(const gsl_vector *v, const gsl_vector *w) ;
double apop_vector_weighted_cov(const gsl_vector *, const gsl_vector *, const gsl_vector *) ;
double apop_vector_weighted_skew(const gsl_vector *v, const gsl_vector *w) ;
double apop_vector_weighted_kurt(const gsl_vector *v, const gsl_vector *w) ;

//Distances, Euclidian and Manhattan:
double apop_vector_distance(const gsl_vector *ina, const gsl_vector *inb);
double apop_vector_grid_distance(const gsl_vector *ina, const gsl_vector *inb);


void apop_vector_normalize(gsl_vector *in, gsl_vector **out, const char normalization_type);
void apop_matrix_normalize(gsl_matrix *data, const char row_or_col, const char normalization);

inline double apop_test_chi_squared_var_not_zero(const gsl_vector *in);

double apop_random_double(double min, double max, gsl_rng *r);
int apop_random_int(const double min, const double max, const gsl_rng *r);

gsl_matrix *apop_matrix_covariance(gsl_matrix *in, const char normalize);
gsl_matrix *apop_matrix_correlation(gsl_matrix *in, const char normalize);
//apop_data * apop_data_covariance(const apop_data *in);
//apop_data * apop_data_correlation(const apop_data *in);
long double apop_matrix_sum(const gsl_matrix *m) ;
double apop_matrix_mean(const gsl_matrix *data) ;
double apop_matrix_var_m(const gsl_matrix *data, double mean) ;
void apop_matrix_mean_and_var(const gsl_matrix *data, double *mean, double *var);
double apop_rng_GHgB3(gsl_rng * r, double* a); //in asst.c
//apop_data * apop_data_summarize(apop_data *data);
apop_data * apop_matrix_summarize(gsl_matrix *data);

apop_data *apop_test_fisher_exact(apop_data *intab);

double * apop_vector_percentiles(gsl_vector *data, char rounding);
apop_data * apop_data_sort(apop_data *data, int sortby, char asc);

extern apop_model apop_beta;
extern apop_model apop_bernoulli;
extern apop_model apop_binomial; //on hiatus.
extern apop_model apop_exponential;
extern apop_model apop_gamma;
extern apop_model apop_gaussian;//synonym for apop_normal
extern apop_model apop_histogram;
extern apop_model apop_improper_uniform;
extern apop_model apop_iv;
extern apop_model apop_kernel_density;
extern apop_model apop_logit;
extern apop_model apop_lognormal;
extern apop_model apop_multinomial_probit;
extern apop_model apop_multivariate_normal;
extern apop_model apop_normal;
extern apop_model apop_ols;
extern apop_model apop_poisson;
extern apop_model apop_probit;
extern apop_model apop_uniform;
extern apop_model apop_waring;
extern apop_model apop_wls;
extern apop_model apop_yule;
extern apop_model apop_zipf;

#define apop_OLS apop_ols
#define apop_WLS apop_wls
#define apop_IV apop_iv

apop_ls_settings * apop_ls_settings_alloc(apop_data *data);
void * apop_ls_settings_copy(apop_ls_settings *in);
void apop_ls_settings_free(apop_ls_settings *in);

typedef struct {
    int want_cov;
    void *copy;
    void *free;
} apop_normal_settings;

apop_normal_settings *apop_normal_settings_alloc(int want_cov);
apop_normal_settings *apop_normal_settings_copy(apop_normal_settings *in);
void apop_normal_settings_free(apop_normal_settings *in);

typedef struct {
    apop_data *factors; char source_type; char source_column; apop_data
    *source_data;
} apop_category_settings;

apop_category_settings *apop_category_settings_alloc(apop_data *d, int source_column, char source_type);
apop_category_settings *apop_category_settings_copy(apop_category_settings *in);
void apop_category_settings_free(apop_category_settings *in);

typedef struct {
    char rank_data;
    void *copy;
    void *free;
} apop_rank_settings;

apop_rank_settings *apop_rank_settings_alloc(void *ignoreme);
void apop_rank_settings_free(apop_rank_settings *in);
void *apop_rank_settings_copy(apop_rank_settings *in);

typedef struct{
    gsl_histogram       *pdf;
    gsl_histogram_pdf   *cdf;
    apop_model          *histobase;
    apop_model          *kernelbase;
} apop_histogram_settings;

#define apop_kernel_density_settings apop_histogram_settings

apop_histogram_settings *apop_histogram_settings_alloc(apop_data *data, int bins);
void  apop_histogram_settings_free(apop_histogram_settings *in);
void * apop_histogram_settings_copy(apop_histogram_settings *in);

apop_model *apop_model_set_parameters(apop_model in, ...);
apop_histogram_settings *apop_kernel_density_settings_alloc(apop_data *data, 
        apop_model *histobase, apop_model *kernelbase, void (*set_params)(double, apop_model*));

#define apop_kernel_density_settings_copy apop_histogram_settings_copy
#define apop_kernel_density_settings_free apop_histogram_settings_free

apop_model * apop_model_copy(apop_model in); //in apop_model.c
apop_model * apop_model_clear(apop_data * data, apop_model *model);

apop_model * apop_estimate(apop_data *d, apop_model m);
void apop_score(apop_data *d, gsl_vector *out, apop_model *m);
double apop_log_likelihood(apop_data *d, apop_model *m);
double apop_p(apop_data *d, apop_model *m);
void apop_draw(double *out, gsl_rng *r, apop_model *m);
void apop_model_prep(apop_data *d, apop_model *m);
apop_data * apop_expected_value(apop_data *d, apop_model *m);

apop_model *apop_beta_from_mean_var(double m, double v);


/* matrix/gsl_matrix_double.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

typedef struct 
{
  size_t size1;
  size_t size2;
  size_t tda;
  double * data;
  gsl_block * block;
  int owner;
} gsl_matrix;

typedef struct
{
  gsl_matrix matrix;
} _gsl_matrix_view;

typedef _gsl_matrix_view gsl_matrix_view;

typedef struct
{
  gsl_matrix matrix;
} _gsl_matrix_const_view;

typedef const _gsl_matrix_const_view gsl_matrix_const_view;

/* Allocation */

gsl_matrix * gsl_matrix_alloc (const size_t n1, const size_t n2);
gsl_matrix * gsl_matrix_calloc (const size_t n1, const size_t n2);
gsl_matrix * gsl_matrix_alloc_from_block (gsl_block * b, 
                                   const size_t offset, 
                                   const size_t n1, 
                                   const size_t n2, 
                                   const size_t d2);

gsl_matrix * gsl_matrix_alloc_from_matrix (gsl_matrix * m,
                                    const size_t k1, 
                                    const size_t k2,
                                    const size_t n1, 
                                    const size_t n2);

gsl_vector * gsl_vector_alloc_row_from_matrix (gsl_matrix * m,
                                        const size_t i);

gsl_vector * gsl_vector_alloc_col_from_matrix (gsl_matrix * m,
                                        const size_t j);

void gsl_matrix_free (gsl_matrix * m);

/* Views */

_gsl_matrix_view gsl_matrix_submatrix (gsl_matrix * m, 
                            const size_t i, const size_t j, 
                            const size_t n1, const size_t n2);

_gsl_vector_view gsl_matrix_row (gsl_matrix * m, const size_t i);
_gsl_vector_view gsl_matrix_column (gsl_matrix * m, const size_t j);
_gsl_vector_view gsl_matrix_diagonal (gsl_matrix * m);
_gsl_vector_view gsl_matrix_subdiagonal (gsl_matrix * m, const size_t k);
_gsl_vector_view gsl_matrix_superdiagonal (gsl_matrix * m, const size_t k);
_gsl_matrix_view gsl_matrix_view_array (double * base,
                             const size_t n1, 
                             const size_t n2);

_gsl_matrix_view
gsl_matrix_view_array_with_tda (double * base, 
                                      const size_t n1, 
                                      const size_t n2,
                                      const size_t tda);


_gsl_matrix_view
gsl_matrix_view_vector (gsl_vector * v,
                              const size_t n1, 
                              const size_t n2);

_gsl_matrix_view
gsl_matrix_view_vector_with_tda (gsl_vector * v,
                                       const size_t n1, 
                                       const size_t n2,
                                       const size_t tda);


_gsl_matrix_const_view 
gsl_matrix_const_submatrix (const gsl_matrix * m, 
                                  const size_t i, const size_t j, 
                                  const size_t n1, const size_t n2);

_gsl_vector_const_view 
gsl_matrix_const_row (const gsl_matrix * m, 
                            const size_t i);

_gsl_vector_const_view 
gsl_matrix_const_column (const gsl_matrix * m, 
                               const size_t j);

_gsl_vector_const_view
gsl_matrix_const_diagonal (const gsl_matrix * m);

_gsl_vector_const_view 
gsl_matrix_const_subdiagonal (const gsl_matrix * m, 
                                    const size_t k);

_gsl_vector_const_view 
gsl_matrix_const_superdiagonal (const gsl_matrix * m, 
                                      const size_t k);

_gsl_matrix_const_view
gsl_matrix_const_view_array (const double * base,
                                   const size_t n1, 
                                   const size_t n2);

_gsl_matrix_const_view
gsl_matrix_const_view_array_with_tda (const double * base, 
                                            const size_t n1, 
                                            const size_t n2,
                                            const size_t tda);

_gsl_matrix_const_view
gsl_matrix_const_view_vector (const gsl_vector * v,
                                    const size_t n1, 
                                    const size_t n2);

_gsl_matrix_const_view
gsl_matrix_const_view_vector_with_tda (const gsl_vector * v,
                                             const size_t n1, 
                                             const size_t n2,
                                             const size_t tda);

/* Operations */

double   gsl_matrix_get(const gsl_matrix * m, const size_t i, const size_t j);
void    gsl_matrix_set(gsl_matrix * m, const size_t i, const size_t j, const double x);

double * gsl_matrix_ptr(gsl_matrix * m, const size_t i, const size_t j);
const double * gsl_matrix_const_ptr(const gsl_matrix * m, const size_t i, const size_t j);

void gsl_matrix_set_zero (gsl_matrix * m);
void gsl_matrix_set_identity (gsl_matrix * m);
void gsl_matrix_set_all (gsl_matrix * m, double x);

int gsl_matrix_fread (FILE * stream, gsl_matrix * m) ;
int gsl_matrix_fwrite (FILE * stream, const gsl_matrix * m) ;
int gsl_matrix_fscanf (FILE * stream, gsl_matrix * m);
int gsl_matrix_fprintf (FILE * stream, const gsl_matrix * m, const char * format);
 
int gsl_matrix_memcpy(gsl_matrix * dest, const gsl_matrix * src);
int gsl_matrix_swap(gsl_matrix * m1, gsl_matrix * m2);

int gsl_matrix_swap_rows(gsl_matrix * m, const size_t i, const size_t j);
int gsl_matrix_swap_columns(gsl_matrix * m, const size_t i, const size_t j);
int gsl_matrix_swap_rowcol(gsl_matrix * m, const size_t i, const size_t j);
int gsl_matrix_transpose (gsl_matrix * m);
int gsl_matrix_transpose_memcpy (gsl_matrix * dest, const gsl_matrix * src);

double gsl_matrix_max (const gsl_matrix * m);
double gsl_matrix_min (const gsl_matrix * m);
void gsl_matrix_minmax (const gsl_matrix * m, double * min_out, double * max_out);

void gsl_matrix_max_index (const gsl_matrix * m, size_t * imax, size_t *jmax);
void gsl_matrix_min_index (const gsl_matrix * m, size_t * imin, size_t *jmin);
void gsl_matrix_minmax_index (const gsl_matrix * m, size_t * imin, size_t * jmin, size_t * imax, size_t * jmax);
int gsl_matrix_isnull (const gsl_matrix * m);
int gsl_matrix_add (gsl_matrix * a, const gsl_matrix * b);
int gsl_matrix_sub (gsl_matrix * a, const gsl_matrix * b);
int gsl_matrix_mul_elements (gsl_matrix * a, const gsl_matrix * b);
int gsl_matrix_div_elements (gsl_matrix * a, const gsl_matrix * b);
int gsl_matrix_scale (gsl_matrix * a, const double x);
int gsl_matrix_add_constant (gsl_matrix * a, const double x);
int gsl_matrix_add_diagonal (gsl_matrix * a, const double x);
double gsl_matrix_get(const gsl_matrix * m, const size_t i, const size_t j);
void gsl_matrix_set(gsl_matrix * m, const size_t i, const size_t j, const double x):
double * gsl_matrix_ptr(gsl_matrix * m, const size_t i, const size_t j);
const double * gsl_matrix_const_ptr(const gsl_matrix * m, const size_t i, const size_t j);



/* vector/gsl_vector_double.h
 * copyright as for gsl_matrix above.
 */

typedef struct 
{
  size_t size;
  size_t stride;
  double *data;
  gsl_block *block;
  int owner;
} 
gsl_vector;

typedef struct
{
  gsl_vector vector;
} _gsl_vector_view;

typedef _gsl_vector_view gsl_vector_view;

typedef struct
{
  gsl_vector vector;
} _gsl_vector_const_view;

typedef const _gsl_vector_const_view gsl_vector_const_view;


/* Allocation */

gsl_vector *gsl_vector_alloc (const size_t n);
gsl_vector *gsl_vector_calloc (const size_t n);

gsl_vector *gsl_vector_alloc_from_block (gsl_block * b,
                                                     const size_t offset, 
                                                     const size_t n, 
                                                     const size_t stride);

gsl_vector *gsl_vector_alloc_from_vector (gsl_vector * v,
                                                      const size_t offset, 
                                                      const size_t n, 
                                                      const size_t stride);

void gsl_vector_free (gsl_vector * v);

/* Views */

_gsl_vector_view 
gsl_vector_view_array (double *v, size_t n);

_gsl_vector_view 
gsl_vector_view_array_with_stride (double *base,
                                         size_t stride,
                                         size_t n);

_gsl_vector_const_view 
gsl_vector_const_view_array (const double *v, size_t n);

_gsl_vector_const_view 
gsl_vector_const_view_array_with_stride (const double *base,
                                               size_t stride,
                                               size_t n);

_gsl_vector_view 
gsl_vector_subvector (gsl_vector *v, 
                            size_t i, 
                            size_t n);

_gsl_vector_view 
gsl_vector_subvector_with_stride (gsl_vector *v, 
                                        size_t i,
                                        size_t stride,
                                        size_t n);

_gsl_vector_const_view 
gsl_vector_const_subvector (const gsl_vector *v, 
                                  size_t i, 
                                  size_t n);

_gsl_vector_const_view 
gsl_vector_const_subvector_with_stride (const gsl_vector *v, 
                                              size_t i, 
                                              size_t stride,
                                              size_t n);

/* Operations */

double gsl_vector_get (const gsl_vector * v, const size_t i);
void gsl_vector_set (gsl_vector * v, const size_t i, double x);

double *gsl_vector_ptr (gsl_vector * v, const size_t i);
const double *gsl_vector_const_ptr (const gsl_vector * v, const size_t i);

void gsl_vector_set_zero (gsl_vector * v);
void gsl_vector_set_all (gsl_vector * v, double x);
int gsl_vector_set_basis (gsl_vector * v, size_t i);

int gsl_vector_fread (FILE * stream, gsl_vector * v);
int gsl_vector_fwrite (FILE * stream, const gsl_vector * v);
int gsl_vector_fscanf (FILE * stream, gsl_vector * v);
int gsl_vector_fprintf (FILE * stream, const gsl_vector * v,
                              const char *format);
int gsl_vector_memcpy (gsl_vector * dest, const gsl_vector * src);
int gsl_vector_reverse (gsl_vector * v);
int gsl_vector_swap (gsl_vector * v, gsl_vector * w);
int gsl_vector_swap_elements (gsl_vector * v, const size_t i, const size_t j);
double gsl_vector_max (const gsl_vector * v);
double gsl_vector_min (const gsl_vector * v);
void gsl_vector_minmax (const gsl_vector * v, double * min_out, double * max_out);
size_t gsl_vector_max_index (const gsl_vector * v);
size_t gsl_vector_min_index (const gsl_vector * v);
void gsl_vector_minmax_index (const gsl_vector * v, size_t * imin, size_t * imax);
int gsl_vector_add (gsl_vector * a, const gsl_vector * b);
int gsl_vector_sub (gsl_vector * a, const gsl_vector * b);
int gsl_vector_mul (gsl_vector * a, const gsl_vector * b);
int gsl_vector_div (gsl_vector * a, const gsl_vector * b);
int gsl_vector_scale (gsl_vector * a, const double x);
int gsl_vector_add_constant (gsl_vector * a, const double x);
int gsl_vector_isnull (const gsl_vector * v);
double gsl_vector_get (const gsl_vector * v, const size_t i);
void gsl_vector_set (gsl_vector * v, const size_t i, double x);
double * gsl_vector_ptr (gsl_vector * v, const size_t i);
const double * gsl_vector_const_ptr (const gsl_vector * v, const size_t i);
double gsl_rng_uniform (const gsl_rng * r);
