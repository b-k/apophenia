#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

void apop_invert_matrix(gsl_matrix *in, gsl_matrix *out);
double apop_det_and_inv(gsl_matrix *in, gsl_matrix *out, int calc_det, int calc_inv);
double apop_x_prime_sigma_x(gsl_vector *x, gsl_matrix *sigma);
void apop_sv_decomposition(gsl_matrix *data, int dimensions_we_want, gsl_matrix ** pc_space, gsl_vector **eigenvalues);
void apop_print_matrix(gsl_matrix *data);
void apop_print_matrix_int(gsl_matrix *data);
void apop_print_vector(gsl_vector *data);
void apop_print_vector_int(gsl_vector *data);
