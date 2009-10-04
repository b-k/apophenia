#ifndef APOP_OUTPUT_H
#define APOP_OUTPUT_H
#include <gsl/gsl_matrix.h>
#include "db.h"
#include "types.h"
#include "stats.h"
#include "variadic.h"
#include "linear_algebra.h"

#ifdef	__cplusplus
extern "C" {
#endif

APOP_VAR_DECLARE void apop_plot_line_and_scatter(apop_data *data, apop_model *est, char * output_file, FILE *output_pipe, char output_type, char output_append);
APOP_VAR_DECLARE void apop_histogram_plot(apop_model *hist,char *output_file, FILE *output_pipe, char output_type, char output_append);
APOP_VAR_DECLARE void apop_plot_histogram(gsl_vector *data, size_t bins, char *output_file, FILE *output_pipe, char output_type, char output_append);
APOP_VAR_DECLARE void apop_histogram_print(apop_model *h, char *output_file, FILE *output_pipe, char output_type, char output_append);
APOP_VAR_DECLARE  void apop_plot_lattice(const apop_data *d, char *output_file, FILE *output_pipe, char output_type, char output_append);
APOP_VAR_DECLARE void apop_plot_qq(gsl_vector *v, apop_model *m, char *output_file, FILE *output_pipe, char output_type, char output_append, size_t bins, gsl_rng *r);
APOP_VAR_DECLARE void apop_plot_triangle(apop_data *in, char *output_file, FILE *output_pipe, char output_type, char output_append);

APOP_VAR_DECLARE void apop_matrix_print(const gsl_matrix *data, char *output_file, FILE *output_pipe, char output_type, char output_append);
APOP_VAR_DECLARE void apop_vector_print(gsl_vector *data, char *output_file, FILE *output_pipe, char output_type, char output_append);
APOP_VAR_DECLARE void apop_data_print(const apop_data *data, char *output_file, FILE *output_pipe, char output_type, char output_append);

void apop_matrix_show(const gsl_matrix *data);
void apop_vector_show(const gsl_vector *data);
void apop_data_show(const apop_data *data);

#ifdef	__cplusplus
}
#endif
#endif
