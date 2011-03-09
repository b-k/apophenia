//apop.i   Copyright (c) 2009-10 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.

/* This is the interface file for SWIG, and unfortunately, it is
a mess. There are two systematic difficulties, whose solutions are
something of a slog:

(A) Including apop_data_copy, for example, in the apop_data object will
produce a name clash. So I need to:
    --%ignore the function on the Python/R side.
    --add a %rename directive to change a dummy like __copy to copy
    --In the apop_data object, declare a method named __copy, which the
    above rename directive and Swig's rules will change to your_data.copy
    on the Python/R side.

(B) Swig's system doesn't mesh with the C preprocessor-based method of
implementing C-side named arguments. So I have to start over,
re-declaring the function using Swig's form for default values.
This means redundancy: I have to have the default values set in the code
itself and here; be careful about that.

[Also, swig doesn't yet like C99 __attribs__.]

That said, here is the current table of contents:

--%ignores for issue (A) above
--Apophenia's headers, to be pasted as C code into wrapper.c 
--type mappings, including swig's arrays and some Python conversions like list->apop_data 
--%renames for (A)
--declarations with Swig/c++-style defaults for (B)
--%extend elements that link appropriate functions to objects, so you can use the familiar forms like yourdata.copy() instead of apop_data_copy(yourdata).
--Apophenia's headers, to be parsed by Swig to produce the interfaces
--A subset of the GSL's vector and matrix headers, so Swig will produce wrappers for various GSL functions as well (primarily those that have methods declared in the vector and matrix objects)

*/

//These are all now called via the object.verb form.
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
%ignore apop_data_sort;
%ignore apop_data_split;
%ignore apop_data_stack;
%ignore apop_data_rank_expand;
%ignore apop_data_rank_compress;
%ignore apop_data_summarize;
%ignore apop_data_transpose;
%ignore apop_model_clear;
%ignore apop_model_show;
%ignore apop_name_add;
%ignore apop_name_copy;
%ignore apop_name_stack;
%ignore apop_name_find;
%ignore apop_name_print;

#define APOP_NO_VARIADIC
%module apop

/*Put all the headers in the C code: */
%{
#include "types.h"
#include "db.h"
#include "asst.h"
#include "stats.h"
#include "variadic.h"
#include "settings.h"
#include "deprecated.h"

//Part of the intent of a convenience header like this is that you
//don't have to remember what else you're including. So here are 
//some other common GSL headers:
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>

//And common headers for other uses (such as seeding an RNG):
#include <time.h>
#include <unistd.h>


%}

%include "carrays.i"
%array_class(double, apop_double)
%array_class(int, apop_int)

#ifdef SWIGPYTHON
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

def apop_pylist_to_data(inlist):
    colsize =len(inlist)
    rowsize =len(inlist[1])
    dlist = apop_double(colsize*rowsize)
    if rowsize == 1:
        for i in xrange(colsize):
            dlist[i] = inlist[i]
        out = apop_line_to_data(dlist, rowsize, 0, 0)
    else:
        for i in xrange(colsize):
            for j in xrange(rowsize):
                dlist[i + j*rowsize] = inlist[i][j]
        out = apop_line_to_data(dlist, 0, rowsize, colsize)
    return out
%}
#endif

%rename(add) __add;
%rename(add_constant) __add_constant;
%rename(add_diagonal) __add_diagonal;
%rename(add_named_elmt) __add_named_elmt;
%rename(clear) __clear;
%rename(column) __column;
%rename(copy) __copy;
%rename(correlation) __correlation;
%rename(covariance) __covariance;
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
%rename(sort) __sort;
%rename(stack) __stack;
%rename(sub) __sub;
%rename(summarize) __summarize;
%rename(swap_columns) __swap_columns;
%rename(swap_rowcol) __swap_rowcol;
%rename(swap_rows) __swap_rows;
%rename(transpose_memcpy) __transpose_memcpy;
%rename(transpose) __transpose;
%rename(rank_expand) __rank_expand;
%rename(rank_compress) __rank_compress;

%rename(apop_data) _apop_data;
%rename(apop_model) _apop_model;

#  define __attribute__(Spec) /* empty */

/* Now declare everything that will be included */
%include "types.h"
%include "db.h"
%include "asst.h"
%include "stats.h"
%include "variadic.h"
%include "settings.h"
%include "deprecated.h"

/* Variadics: */
int apop_db_close(char vacuum='q');
void apop_db_merge(char *db_file, char inout='i');
void apop_db_merge_table(char *db_file, char *tabname, char inout='i');
apop_data * apop_text_to_data(char *text_file="-", int has_row_names=0, int has_col_names=1);
int apop_text_to_db(char *text_file="-", char *tabname="t", int has_row_names =0, int has_col_names=1, char **field_names=NULL);
int apop_matrix_is_positive_semidefinite(gsl_matrix *m, char semi='s');
apop_data * apop_f_test (apop_model *est, apop_data *contrast=NULL, int normalize=0);
int apop_table_exists(char *name, char remove='n');
double apop_vector_distance(const gsl_vector *ina, const gsl_vector *inb=NULL, const char metric='e', const double norm=2);
apop_data * apop_dot(const apop_data *d1, const apop_data *d2, char form1=0, char form2=0);
apop_data * apop_data_to_dummies(apop_data *d, int col=0, char type='t', int keep_first=0);
void apop_plot_lattice(const apop_data *d, char *outfile=NULL);


/*
//will require a wrapper function
%rename(test) apop_test;
%rename(plot_qq) apop_plot_qq;
%rename(apop_bootstrap_cov) apop_bootstrap_cov_base;
%rename(apop_plot_histogram) apop_plot_histogram_base;
%rename(apop_linear_constraint) apop_linear_constraint_base;
%rename(apop_model_set_parameters) apop_model_set_parameters_base;
%rename(apop_histogram_model_reset) apop_histogram_model_reset_base;
*/

%extend _apop_data {
    apop_data* _apop_data(const size_t v, const size_t m1, const int m2){ return  apop_data_alloc(v, m1, m2); }
    apop_data* _apop_data(const size_t s1, const int s2){ return  apop_data_alloc(s1, s2); }//just the matrix
    apop_data* _apop_data(const size_t s1){ return  apop_data_alloc(s1); } //just the vector
    ~_apop_data()    {apop_data_free($self);}
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
    apop_data** __split(int splitpoint, char r_or_c){ return apop_data_split($self, splitpoint, r_or_c); }
    apop_data*  __rm_columns(int *drop)             { apop_data_rm_columns($self, drop); return $self; }
    void __memcpy(apop_data *out)                   { return apop_data_memcpy(out, $self);}
    void __add_named_elmt(char *name, double val)   { apop_data_add_named_elmt($self, name, val); }
    apop_data* __listwise_delete()                  { return apop_data_listwise_delete_base($self, 'y'); }
    apop_data* __transpose()                        { return apop_data_transpose($self);}
    apop_data* __covariance()                       { return apop_data_covariance($self);}
    apop_data* __correlation()                      { return apop_data_correlation($self);}
    apop_data* __summarize()                        { return apop_data_summarize($self);}
    apop_data* __rank_expand()                      { return apop_data_rank_expand($self);}
    apop_data* __rank_compress()                    { return apop_data_rank_compress($self);}
    apop_data * __stack(apop_data * m2=NULL, char posn='r', char inplace ='n'){
                                     return apop_data_stack_base($self, m2, posn,inplace); }
    apop_data * __sort(int sortby=0, char asc='a'){
                                     return apop_data_sort($self, sortby, asc); }
};

%extend _apop_model{
    ~_apop_model (){ apop_model_free ($self);}

    void show (){ apop_model_show ($self);}
    char* __str__() {apop_model_show($self); return " ";}

    apop_model * __clear(apop_data * data){
        return apop_model_clear(data, $self);}
}

%extend apop_name {
            apop_name()     { return apop_name_alloc(); }
            ~apop_name()    { apop_name_free($self); }
    void    __print()       { apop_name_print($self); }
    apop_name * __copy()     { return apop_name_copy($self); }
    int     __add(char *add_me, char type)      { return apop_name_add($self, add_me, type); }
    void    __stack(apop_name *n2, char type, char type2)   { apop_name_stack($self, n2, type, type2); }
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
    long double sum()                   { return apop_matrix_sum($self) ;}
    double      mean()                  { return apop_matrix_mean($self) ;}

    void normalize(const char row_or_col, const char normalization){
        return apop_matrix_normalize($self, row_or_col, normalization); }

    int is_positive_semidefinite(char semi){
    return apop_matrix_is_positive_semidefinite_base($self, semi);
    }

    double to_positive_semidefinite( ){
        return apop_matrix_to_positive_semidefinite($self); }

    apop_data*  pca( int dimensions_we_want){
    return  apop_matrix_pca($self, dimensions_we_want);}

    inline void increment( int i, int j, double amt){
        return apop_matrix_increment($self, i, j, amt);}

    gsl_matrix *stack( gsl_matrix * m2=NULL, char posn='r', char inplace=0){
        return apop_matrix_stack($self,  m2, posn, inplace);}

    gsl_matrix *rm_columns( int *drop){
    return apop_matrix_rm_columns($self, drop);}

    gsl_matrix *covariance(const char normalize=0){
        return apop_matrix_covariance_base($self, normalize);}

    gsl_matrix * correlation(const char normalize=0){
        return apop_matrix_correlation_base($self, normalize);}

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
    gsl_matrix * to_matrix(char row_col='c') { return apop_vector_to_matrix($self, row_col); }
    double      kurtosis_pop()          { return apop_vector_kurtosis_pop($self) ; }
    double      kurtosis()              { return apop_vector_kurtosis($self) ; }
    double      skew()                  { return apop_vector_skew($self) ; }
    double      skew_pop()              { return apop_vector_skew_pop($self) ; }
    double      kurt()                  { return apop_vector_kurt($self) ; }
    double      cov(const gsl_vector *inb)   { return apop_vector_cov($self, inb) ; }

    double correlation(const gsl_vector *inb) { return apop_vector_correlation($self, inb) ; }
    gsl_vector * moving_average(size_t window){ return  apop_vector_moving_average($self, window);}

    //double * percentiles(char rounding) {return apop_vector_percentiles_base($self, rounding);}
    apop_double * percentiles(char rounding='d') {return apop_vector_percentiles_base($self, rounding);}
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
        return apop_vector_distance_base($self, inb, 'e', 0); }

    double grid_distance(const gsl_vector *inb){
        return apop_vector_grid_distance($self, inb);}

    void normalize(gsl_vector **out=NULL, const char normalization_type = 'p'){
        apop_vector_normalize($self, out, normalization_type); }

    void increment(gsl_vector * v, int i, double amt){
        return apop_vector_increment($self, i, amt);}

    int bounded(double max=GSL_POSINF){ return apop_vector_bounded($self, max);}
    gsl_vector *apop_vector_stack(gsl_vector *v2=NULL, char inplace=0){ return apop_vector_stack($self, v2, 'n');}

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
//Tightened up a bit by BK because you are not going to read this anyway.

typedef struct 
{ size_t size1;
  size_t size2;
  size_t tda;
  double * data;
  gsl_block * block;
  int owner;
} gsl_matrix;

typedef struct
{ gsl_matrix matrix;
} _gsl_matrix_view;

typedef _gsl_matrix_view gsl_matrix_view;

typedef struct
{ gsl_matrix matrix;
} _gsl_matrix_const_view;

typedef const _gsl_matrix_const_view gsl_matrix_const_view;

/* Allocation */
gsl_matrix * gsl_matrix_calloc (const size_t n1, const size_t n2);
gsl_matrix * gsl_matrix_alloc_from_block (gsl_block * b, const size_t offset, const size_t n1, const size_t n2, const size_t d2);
gsl_matrix * gsl_matrix_alloc_from_matrix (gsl_matrix * m, const size_t k1, const size_t k2, const size_t n1, const size_t n2);
gsl_vector * gsl_vector_alloc_row_from_matrix (gsl_matrix * m, const size_t i); 
gsl_vector * gsl_vector_alloc_col_from_matrix (gsl_matrix * m, const size_t j);

/* Views */
_gsl_matrix_view gsl_matrix_submatrix (gsl_matrix * m, const size_t i, const size_t j, const size_t n1, const size_t n2);
_gsl_vector_view gsl_matrix_diagonal (gsl_matrix * m);
_gsl_vector_view gsl_matrix_subdiagonal (gsl_matrix * m, const size_t k);
_gsl_vector_view gsl_matrix_superdiagonal (gsl_matrix * m, const size_t k);
_gsl_matrix_view gsl_matrix_view_array (double * base, const size_t n1, const size_t n2); 
_gsl_matrix_view gsl_matrix_view_array_with_tda (double * base, const size_t n1, const size_t n2, const size_t tda); 
_gsl_matrix_view gsl_matrix_view_vector (gsl_vector * v, const size_t n1, const size_t n2); 
_gsl_matrix_view gsl_matrix_view_vector_with_tda (gsl_vector * v, const size_t n1, const size_t n2, const size_t tda); 
_gsl_matrix_const_view gsl_matrix_const_submatrix (const gsl_matrix * m, const size_t i, const size_t j, const size_t n1, const size_t n2); 
_gsl_vector_const_view gsl_matrix_const_row (const gsl_matrix * m, const size_t i);
_gsl_vector_const_view gsl_matrix_const_column (const gsl_matrix * m, const size_t j);
_gsl_vector_const_view gsl_matrix_const_diagonal (const gsl_matrix * m);
_gsl_vector_const_view gsl_matrix_const_subdiagonal (const gsl_matrix * m, const size_t k);
_gsl_vector_const_view gsl_matrix_const_superdiagonal (const gsl_matrix * m, const size_t k);
_gsl_matrix_const_view gsl_matrix_const_view_array (const double * base, const size_t n1, const size_t n2); 
_gsl_matrix_const_view gsl_matrix_const_view_array_with_tda (const double * base, const size_t n1, const size_t n2, const size_t tda); 
_gsl_matrix_const_view gsl_matrix_const_view_vector (const gsl_vector * v, const size_t n1, const size_t n2); 
_gsl_matrix_const_view gsl_matrix_const_view_vector_with_tda (const gsl_vector * v, const size_t n1, const size_t n2, const size_t tda);

/* Operations */
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

void gsl_matrix_minmax (const gsl_matrix * m, double * min_out, double * max_out);

void gsl_matrix_max_index (const gsl_matrix * m, size_t * imax, size_t *jmax);
void gsl_matrix_min_index (const gsl_matrix * m, size_t * imin, size_t *jmin);
void gsl_matrix_minmax_index (const gsl_matrix * m, size_t * imin, size_t * jmin, size_t * imax, size_t * jmax);
const double * gsl_matrix_const_ptr(const gsl_matrix * m, const size_t i, const size_t j);

/* vector/gsl_vector_double.h
 * copyright as for gsl_matrix above.
 */

typedef struct 
{ size_t size;
  size_t stride;
  double *data;
  gsl_block *block;
  int owner;
} 
gsl_vector;

typedef struct
{ gsl_vector vector;
} _gsl_vector_view;

typedef _gsl_vector_view gsl_vector_view;

typedef struct
{ gsl_vector vector;
} _gsl_vector_const_view;

typedef const _gsl_vector_const_view gsl_vector_const_view;

/* Allocation */
gsl_vector *gsl_vector_calloc (const size_t n); 
gsl_vector *gsl_vector_alloc_from_block (gsl_block * b, const size_t offset, const size_t n, const size_t stride); 
gsl_vector *gsl_vector_alloc_from_vector (gsl_vector * v, const size_t offset, const size_t n, const size_t stride); 

/* Views */
_gsl_vector_view gsl_vector_view_array (double *v, size_t n); 
_gsl_vector_view gsl_vector_view_array_with_stride (double *base, size_t stride, size_t n); 
_gsl_vector_const_view gsl_vector_const_view_array (const double *v, size_t n); 
_gsl_vector_const_view gsl_vector_const_view_array_with_stride (const double *base, size_t stride, size_t n); 
_gsl_vector_view gsl_vector_subvector (gsl_vector *v, size_t i, size_t n); 
_gsl_vector_view gsl_vector_subvector_with_stride (gsl_vector *v, size_t i, size_t stride, size_t n); 
_gsl_vector_const_view gsl_vector_const_subvector (const gsl_vector *v, size_t i, size_t n); 
_gsl_vector_const_view gsl_vector_const_subvector_with_stride (const gsl_vector *v, size_t i, size_t stride, size_t n);

/* Operations */
const double *gsl_vector_const_ptr (const gsl_vector * v, const size_t i);

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
const double * gsl_vector_const_ptr (const gsl_vector * v, const size_t i);
double gsl_rng_uniform (const gsl_rng * r);
