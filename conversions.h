/*conversions.h 	Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2. 

  Convenience functions to convert among vectors (gsl_vector), matrices (gsl_matrix), 
  arrays (double **), and database tables

  To the extent that it makes sense, all functions are of the form convert_in_to_out(in, out),
  and allocate the output object for you.
*/

#ifndef APOP_CONVERSIONS_H
#define APOP_CONVERSIONS_H

#include "db.h"
#include "variadic.h"
#include <regex.h>
#include <string.h>
#include <stdarg.h>
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

//From vector
gsl_vector *apop_vector_copy(const gsl_vector *in);
double * apop_vector_to_array(const gsl_vector *in);
gsl_matrix * apop_vector_to_matrix(const gsl_vector *in);

//From matrix
gsl_matrix *apop_matrix_copy(const gsl_matrix *in);
apop_data  *apop_db_to_crosstab(char *tabname, char *r1, char *r2, char *datacol);
//takes a three-column table (dim1, dim2, data) and creates a 2D crosstab.
//Returns the crosstab, and the dimension names (if d1!=NULL and d2!=NULL).

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
int apop_crosstab_to_db(apop_data *in, char *tabname, char *row_col_name, 
						char *col_col_name, char *data_col_name);

//packing data into a vector
gsl_vector * apop_data_pack(const apop_data *in);
void apop_data_unpack(const gsl_vector *in, apop_data *d);

char * apop_strip_dots(char *in, char strip_type);

gsl_vector *apop_vector_vfill(gsl_vector *in, va_list ap);
apop_data *apop_data_fill(apop_data *in, ...);
gsl_vector *apop_vector_fill(gsl_vector *in, ...);
gsl_matrix *apop_matrix_fill(gsl_matrix *in, ...);

__END_DECLS
#endif
