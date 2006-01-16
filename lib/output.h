#include <gsl/gsl_matrix.h>
#include <apophenia/db.h>
#include <apophenia/types.h>
#include <apophenia/stats.h>
#include <apophenia/linear_algebra.h>


void apop_plot_line_and_scatter(apop_data *data, apop_estimate *est);
void apop_plot(gsl_matrix *data, char plot_type, int delay);
void apop_plot_histogram(gsl_vector *data, size_t bin_ct, char *outfile);

void apop_print_to_file(char *filename, const char *fmt, ...);
void apop_print_matrix(gsl_matrix *data);
void apop_print_matrix_int(gsl_matrix *data);
void apop_print_vector(gsl_vector *data);
void apop_print_vector_int(gsl_vector *data);

void apop_matrix_print(gsl_matrix *data);
void apop_matrix_print_int(gsl_matrix *data);
void apop_vector_print(gsl_vector *data);
void apop_vector_print_int(gsl_vector *data);

void apop_data_print(apop_data *data);
void apop_data_print_int(apop_data *data);
