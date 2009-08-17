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
%ignore apop_data_summarize;
%ignore apop_data_transpose;
%ignore apop_model_clear;
%ignore apop_name_add;
%ignore apop_name_copy;
%ignore apop_name_stack;
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
%ignore apop_data_stack;


%module apop
%{
#include "asst.h"
#include "stats.h"
#include "apop.h"
#include "model.h"
#include "conversions.h"
#include "db.h"
#include "documentation.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include "mapply.h"
#include "output.h"
#include "settings.h"
#include "types.h"
#include "variadic.h"
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
%rename(apop_db_close) apop_db_close_base;
%rename(apop_dot) apop_dot_base;
%rename(apop_f_test) apop_f_test_base;
%rename(apop_histogram_plot) apop_histogram_plot_base;
%rename(apop_plot_histogram) apop_plot_histogram_base;
%rename(apop_histogram_print) apop_histogram_print_base;
%rename(apop_plot_lattice) apop_plot_lattice_base;
%rename(apop_plot_qq) apop_plot_qq_base;
%rename(apop_table_exists) apop_table_exists_base;
%rename(apop_vector_bounded) apop_vector_bounded_base;
%rename(apop_vector_normalize) apop_vector_normalize_base;
%rename(apop_histogram_model_reset) apop_histogram_model_reset_base;
%rename(apop_bootstrap_cov) apop_bootstrap_cov_base;
%rename(apop_linear_constraint) apop_linear_constraint_base;
%rename(apop_data_sort) apop_data_sort_base;
%rename(apop_vector_percentiles) apop_vector_percentiles_base;
%rename(apop_vector_distance) apop_vector_distance_base;
%rename(apop_matrix_covariance) apop_matrix_covariance_base;
%rename(apop_matrix_correlation) apop_matrix_correlation_base;
%rename(apop_test) apop_test_base;
%rename(apop_matrix_stack) apop_matrix_stack_base;
%rename(apop_vector_stack) apop_vector_stack_base;
%rename(apop_vector_to_matrix) apop_vector_to_matrix_base;
%rename(apop_data_stack) apop_data_stack_base;
%rename(apop_db_merge) apop_db_merge_base;
%rename(apop_db_merge_table) apop_db_merge_table_base;
%rename(apop_random_int) apop_random_int_base;
%rename(apop_model_set_parameters) apop_model_set_parameters_base;
%rename(apop_matrix_is_positive_semidefinite) apop_matrix_is_positive_semidefinite_base;

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
#    apop_data * __stack(apop_data * m2, char posn,char inplace)  { return apop_data_stack_base($self, m2, posn,inplace); }
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

    gsl_matrix *stack( gsl_matrix * m2, char posn){
        return apop_matrix_stack($self,  m2, posn);}

    gsl_matrix *rm_columns( int *drop){
    return apop_matrix_rm_columns($self, drop);}

    gsl_matrix *covariance(const char normalize){
        return apop_matrix_covariance_base($self, normalize);}

    gsl_matrix * correlation(const char normalize){
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
    gsl_matrix * to_matrix(char row_col) { return apop_vector_to_matrix($self, row_col); }
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

    //double * percentiles(char rounding) {return apop_vector_percentiles_base($self, rounding);}
    apop_double * percentiles(char rounding) {return apop_vector_percentiles_base($self, rounding);}
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

    void normalize(gsl_vector **out, const char normalization_type){
        apop_vector_normalize($self, out, normalization_type); }

    void increment(gsl_vector * v, int i, double amt){
        return apop_vector_increment($self, i, amt);}

    int  bounded(long double max){ return apop_vector_bounded($self, max);}
    gsl_vector *stack(gsl_vector * v2){ return apop_vector_stack($self, v2, 'n');}

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
#define APOP_NO_VARIADIC
#  define __attribute__(Spec) /* empty */
%include "asst.h"
%include "stats.h"
%include "apop.h"
%include "model.h"
%include "conversions.h"
%include "db.h"
%include "documentation.h"
%include "likelihoods.h"
%include "linear_algebra.h"
%include "mapply.h"
%include "output.h"
%include "settings.h"
%include "types.h"
%include "variadic.h"


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
//Tightened up a bit by BK because you're not going to read this anyway.

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

gsl_matrix * gsl_matrix_alloc (const size_t n1, const size_t n2);
gsl_matrix * gsl_matrix_calloc (const size_t n1, const size_t n2);
gsl_matrix * gsl_matrix_alloc_from_block (gsl_block * b, const size_t offset, const size_t n1, const size_t n2, const size_t d2);
gsl_matrix * gsl_matrix_alloc_from_matrix (gsl_matrix * m, const size_t k1, const size_t k2, const size_t n1, const size_t n2);
gsl_vector * gsl_vector_alloc_row_from_matrix (gsl_matrix * m, const size_t i); 
gsl_vector * gsl_vector_alloc_col_from_matrix (gsl_matrix * m, const size_t j);
void gsl_matrix_free (gsl_matrix * m);

/* Views */

_gsl_matrix_view gsl_matrix_submatrix (gsl_matrix * m, const size_t i, const size_t j, const size_t n1, const size_t n2);
_gsl_vector_view gsl_matrix_row (gsl_matrix * m, const size_t i);
_gsl_vector_view gsl_matrix_column (gsl_matrix * m, const size_t j);
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

gsl_vector *gsl_vector_alloc (const size_t n);
gsl_vector *gsl_vector_calloc (const size_t n); 
gsl_vector *gsl_vector_alloc_from_block (gsl_block * b, const size_t offset, const size_t n, const size_t stride); 
gsl_vector *gsl_vector_alloc_from_vector (gsl_vector * v, const size_t offset, const size_t n, const size_t stride); 
void gsl_vector_free (gsl_vector * v);

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
