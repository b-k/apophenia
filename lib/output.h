#include <gsl/gsl_matrix.h>
#include <apophenia/db.h>
#include <apophenia/types.h>
#include <apophenia/stats.h>
#include <apophenia/linear_algebra.h>


void apop_plot_line_and_scatter(gsl_matrix *data, apop_estimate *est, apop_name *n, char *outfile);
void apop_plot(gsl_matrix *data, char plot_type, int delay);
void apop_plot_histogram(gsl_vector *data, size_t bin_ct, char *outfile);

void apop_print_to_file(char *filename, const char *fmt, ...);
void apop_print_matrix(gsl_matrix *data, char *separator, char* file);
void apop_print_matrix_int(gsl_matrix *data, char *separator, char* file);
void apop_print_vector(gsl_vector *data, char *separator, char* file);
void apop_print_vector_int(gsl_vector *data, char *separator, char* file);

void apop_matrix_print(gsl_matrix *data, char *separator, char* file);
void apop_matrix_print_int(gsl_matrix *data, char *separator, char* file);
void apop_vector_print(gsl_vector *data, char *separator, char* file);
void apop_vector_print_int(gsl_vector *data, char *separator, char* file);

gsl_matrix * apop_matrix_summarize(gsl_matrix *data, apop_name *names_in, apop_name **names_out); //void apop_matrix_summarize(gsl_matrix *data, apop_name *names);
