#include <gsl/gsl_matrix.h>

void apop_GLS(gsl_matrix *data, gsl_matrix *sigma, gsl_vector **beta);
//Returns GLS parameter estimates in beta.
//Destroys the data in the process.

void apop_OLS(gsl_matrix *data, gsl_vector **beta);
