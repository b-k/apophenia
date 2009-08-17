/** \file conversions.h 

  Convenience functions to convert among vectors (gsl_vector), matrices (gsl_matrix), 
  arrays (double **), and database tables
*/
/*	Copyright (c) 2005--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#ifndef APOP_CONVERSIONS_H
#define APOP_CONVERSIONS_H

#include "db.h"
#include "variadic.h"
#include <regex.h>
#include <string.h>
#include <stdarg.h>
#include <gsl/gsl_matrix.h>

#ifdef	__cplusplus
extern "C" {
#endif

//From vector
gsl_vector *apop_vector_copy(const gsl_vector *in);
double * apop_vector_to_array(const gsl_vector *in);
APOP_VAR_DECLARE gsl_matrix * apop_vector_to_matrix(const gsl_vector *in, char row_col);

//From matrix
gsl_matrix *apop_matrix_copy(const gsl_matrix *in);
apop_data  *apop_db_to_crosstab(char *tabname, char *r1, char *r2, char *datacol);

//From array
APOP_VAR_DECLARE gsl_vector * apop_array_to_vector(double *in, int size);
gsl_matrix * apop_array_to_matrix(const double **in, const int rows, const int cols);
apop_data * apop_array_to_data(const double **in, const int rows, const int cols);

//From line
gsl_matrix * apop_line_to_matrix(double *line, int rows, int cols);
apop_data * apop_line_to_data(double *in, int vsize, int rows, int cols);

//From text
APOP_VAR_DECLARE apop_data * apop_text_to_data(char *text_file, int has_row_names, int has_col_names);
APOP_VAR_DECLARE int apop_text_to_db(char *text_file, char *tabname, int has_row_names, int has_col_names, char **field_names);

//From crosstabs
void apop_crosstab_to_db(apop_data *in, char *tabname, char *row_col_name, 
						char *col_col_name, char *data_col_name);

//packing data into a vector
gsl_vector * apop_data_pack(const apop_data *in);
void apop_data_unpack(const gsl_vector *in, apop_data *d);

char * apop_strip_dots(char *in, char strip_type);

#define apop_vector_fill(in, ...) apop_vector_fill_base((in), (double []) {__VA_ARGS__})
#define apop_data_fill(in, ...) apop_data_fill_base((in), (double []) {__VA_ARGS__})
#define apop_matrix_fill(in, ...) apop_matrix_fill_base((in), (double []) {__VA_ARGS__})
apop_data *apop_data_fill_base(apop_data *in, double []);
gsl_vector *apop_vector_fill_base(gsl_vector *in, double []);
gsl_matrix *apop_matrix_fill_base(gsl_matrix *in, double []);

#ifdef	__cplusplus
}
#endif
#endif
