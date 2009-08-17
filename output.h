#ifndef APOP_OUTPUT_H
#define APOP_OUTPUT_H
#include <gsl/gsl_matrix.h>
#include "db.h"
#include "variadic.h"
#include "types.h"
#include "stats.h"
#include "linear_algebra.h"

#ifdef	__cplusplus
extern "C" {
#endif

void apop_plot_line_and_scatter(apop_data *data, apop_model *est, char *);
APOP_VAR_DECLARE void apop_histogram_plot(apop_model *hist, char *outfile);
APOP_VAR_DECLARE void apop_plot_histogram(gsl_vector *data, size_t bins, char *outfile);
APOP_VAR_DECLARE void apop_histogram_print(apop_model *h, char *outfile);
APOP_VAR_DECLARE void apop_plot_lattice(const apop_data *d, char *outfile);
APOP_VAR_DECLARE void apop_plot_qq(gsl_vector *v, apop_model *m, char *outfile, size_t bins, gsl_rng *r);

void apop_matrix_print(gsl_matrix *data, char *file);
void apop_vector_print(gsl_vector *data, char *file);
void apop_data_print(apop_data *data, char *file);

void apop_matrix_show(const gsl_matrix *data);
void apop_vector_show(const gsl_vector *data);
void apop_data_show(const apop_data *data);

#ifdef	__cplusplus
}
#endif
#endif
