#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

void invert_matrix(gsl_matrix *in, gsl_matrix *out);
double det_and_inv(gsl_matrix *in, gsl_matrix *out, int calc_det, int calc_inv);
double x_prime_sigma_x(gsl_vector *x, gsl_matrix *sigma);
void sv_decomposition(gsl_matrix *data, int dimensions_we_want, gsl_matrix ** pc_space, gsl_vector **eigenvalues);
void print_matrix(gsl_matrix *data);
void print_matrix_int(gsl_matrix *data);
void print_vector(gsl_vector *data);
void print_vector_int(gsl_vector *data);
