#include "regression.h"
#include <gsl/gsl_blas.h>
#include "stats.h"
#include "linear_algebra.h"

void apop_GLS(gsl_matrix *data, gsl_matrix *sigma, gsl_vector **beta){
//Returns GLS parameter estimates in beta.
//Destroys the data in the process.
gsl_vector 	*y_data		= gsl_vector_alloc(data->size1);
gsl_matrix 	*temp		= gsl_matrix_calloc(data->size2, data->size1);
gsl_vector 	*xsy 		= gsl_vector_calloc(data->size2);
gsl_matrix 	*xsx 		= gsl_matrix_calloc(data->size2, data->size2);
gsl_matrix 	*sigma_inverse	= gsl_matrix_alloc(data->size1, data->size1);
gsl_vector_view	v 		= gsl_matrix_column(data, 0);
	*beta	= gsl_vector_calloc(data->size2);
	gsl_matrix_get_col(y_data, data, 0);
	apop_normalize_matrix(data);		//every column should have mean zero.
	gsl_vector_set_all(&(v.vector), 1);	//affine: first column is ones.
	apop_det_and_inv(sigma, sigma_inverse, 0, 1);					//find sigma^{-1}
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, data, sigma_inverse, 0, temp); 	//temp = X' \sigma^{-1}.
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, temp, data, 0, xsx);    		//(X' \sigma^{-1} X)
	gsl_blas_dgemv(CblasNoTrans, 1, temp, y_data, 0, xsy);     			//(X' \sigma^{-1} y)
	gsl_linalg_HH_solve (xsx, xsy, *beta);

	gsl_matrix_free(sigma_inverse); gsl_matrix_free(xsx); gsl_matrix_free(temp);
	gsl_vector_free(y_data); gsl_vector_free(xsy);
}

void apop_OLS(gsl_matrix *data, gsl_vector **beta){
//Returns GLS parameter estimates in beta.
//Destroys the data in the process.
gsl_vector 	*y_data		= gsl_vector_alloc(data->size1);
gsl_vector 	*xpy 		= gsl_vector_calloc(data->size2);
gsl_matrix 	*xpx 		= gsl_matrix_calloc(data->size2, data->size2);
gsl_vector_view	v 		= gsl_matrix_column(data, 0);
	*beta	= gsl_vector_calloc(data->size2);
	gsl_matrix_get_col(y_data, data, 0);
	apop_normalize_matrix(data);		//every column should have mean zero.
	gsl_vector_set_all(&(v.vector), 1);	//affine: first column is ones.
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, data, data, 0, xpx); //(X'X)
	gsl_blas_dgemv(CblasTrans, 1, data, y_data, 0, xpy);     	//(X'y)
	gsl_linalg_HH_solve (xpx, xpy, *beta);

	gsl_matrix_free(xpx); gsl_vector_free(y_data); gsl_vector_free(xpy);
}
