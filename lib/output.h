#include <gsl/gsl_matrix.h>
#include <apophenia/db.h>
#include <apophenia/types.h>
#include <apophenia/stats.h>
#include <apophenia/linear_algebra.h>


void apop_plot_line_and_scatter(apop_data *data, apop_estimate *est, char *);
void apop_plot(gsl_matrix *data, char plot_type, int delay);
void apop_plot_histogram(gsl_vector *data, size_t bin_ct, char *outfile);

void apop_matrix_print(gsl_matrix *data, char *file);
void apop_matrix_print_int(gsl_matrix *data, char *file);
void apop_vector_print(gsl_vector *data, char *file);
void apop_vector_print_int(gsl_vector *data, char *file);
void apop_data_print(apop_data *data, char *file);
void apop_data_print_int(apop_data *data, char *file);

void apop_matrix_show(gsl_matrix *data);
void apop_matrix_show_int(gsl_matrix *data);
void apop_vector_show(gsl_vector *data);
void apop_vector_show_int(gsl_vector *data);
void apop_data_show(apop_data *data);
void apop_data_show_int(apop_data *data);
