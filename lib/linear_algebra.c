#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "math.h" //pow!
void invert_matrix(gsl_matrix *in, gsl_matrix *out)
{
    gsl_matrix *invert_me = gsl_matrix_alloc(in->size1, in->size1);
    gsl_permutation * perm = gsl_permutation_alloc(in->size1);
    int dummy;
    invert_me = gsl_matrix_alloc(in->size1, in->size1);
    gsl_matrix_memcpy (invert_me, in);
    gsl_linalg_LU_decomp(invert_me, perm, &dummy);
    gsl_linalg_LU_invert(invert_me, perm, out);
    gsl_matrix_free(invert_me);
    gsl_permutation_free(perm);
}

double det_and_inv(gsl_matrix *in, gsl_matrix *out, int calc_det, int calc_inv)
{
    gsl_matrix *invert_me = gsl_matrix_alloc(in->size1, in->size1);
    gsl_permutation * perm = gsl_permutation_alloc(in->size1);
    int sign;
    double the_determinant = 0;
    invert_me = gsl_matrix_alloc(in->size1, in->size1);
    gsl_matrix_memcpy (invert_me, in);
    gsl_linalg_LU_decomp(invert_me, perm, &sign);
    if (calc_inv)
    	gsl_linalg_LU_invert(invert_me, perm, out);
    if (calc_det)
    	the_determinant	= gsl_linalg_LU_det(invert_me, sign);
    gsl_matrix_free(invert_me);
    gsl_permutation_free(perm);
    return(the_determinant);
}

double x_prime_sigma_x(gsl_vector *x, gsl_matrix *sigma){
//This comes up often enough that it deserves its own convenience function.
gsl_vector * sigma_dot_x	= gsl_vector_calloc(x->size);
double	the_result;
	//gsl_blas_dgemv(CblasNoTrans, 1, sigma, x, 0, sigma_dot_x);
	gsl_blas_dsymv(CblasUpper, 1, sigma, x, 0, sigma_dot_x); //sigma should be symmetric
	gsl_blas_ddot(x, sigma_dot_x, &the_result);
	gsl_vector_free(sigma_dot_x);
	return(the_result);
}

void normalize_for_svd(gsl_matrix *in){
//Greene (2nd ed, p 271) recommends pre- and post-multiplying by sqrt(diag(X'X)) so that X'X = I.
gsl_vector_view 	v;
gsl_vector		*diagonal = gsl_vector_alloc(in->size1);
int 			i;
	//Get the diagonal, take the square root
	v	= gsl_matrix_diagonal(in);
	gsl_vector_memcpy(diagonal, &(v.vector));
	for (i=0; i<diagonal->size; i++)
		gsl_vector_set(diagonal, i, pow(gsl_vector_get(diagonal,i), .5));
	//mulitply each row and column by the diagonal vector.
	for (i=0; i<diagonal->size; i++){
		v	= gsl_matrix_column(in, i);
		gsl_vector_mul(&(v.vector), diagonal);
		v	= gsl_matrix_row(in, i);
		gsl_vector_mul(&(v.vector), diagonal);
	}
	gsl_vector_free(diagonal);
}

void sv_decomposition(gsl_matrix *data, int dimensions_we_want, gsl_matrix ** pc_space, gsl_vector **total_explained) {
//Get X'X
gsl_matrix * 	eigenvectors 	= gsl_matrix_alloc(data->size2, data->size2);
gsl_vector * 	dummy_v 	= gsl_vector_alloc(data->size2);
gsl_vector * 	all_evalues 	= gsl_vector_alloc(data->size2);
gsl_matrix * 	square  	= gsl_matrix_calloc(data->size2, data->size2);
gsl_vector_view v;
int 		i;
double		eigentotals	= 0;
	*pc_space	= gsl_matrix_alloc(data->size2, dimensions_we_want);
	*total_explained= gsl_vector_alloc(dimensions_we_want);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, data, data, 0, square);
	normalize_for_svd(square);	
	gsl_linalg_SV_decomp(square, eigenvectors, all_evalues, dummy_v);
	for (i=0; i< all_evalues->size; i++)
		eigentotals	+= gsl_vector_get(all_evalues, i);
	for (i=0; i<dimensions_we_want; i++){
		v	= gsl_matrix_column(eigenvectors, i);
		gsl_matrix_set_col(*pc_space, i, &(v.vector));
		gsl_vector_set(*total_explained, i, gsl_vector_get(all_evalues, i)/eigentotals);
	}
	gsl_vector_free(dummy_v); 	gsl_vector_free(all_evalues);
	gsl_matrix_free(square); 	gsl_matrix_free(eigenvectors);
}

void print_matrix(gsl_matrix *data){
int i,j;
	for (i=0; i<data->size1; i++){
		for (j=0; j<data->size2; j++)
			printf("% 5f\t", gsl_matrix_get(data, i, j));
		printf("\n");
	}
}

void print_vector(gsl_vector *data){
int i;
	for (i=0; i<data->size; i++)
			printf("% 5f\t", gsl_vector_get(data, i));
	printf("\n");
}

void print_matrix_int(gsl_matrix *data){
int i,j;
	for (i=0; i<data->size1; i++){
		for (j=0; j<data->size2; j++)
			printf("% i\t", (int) gsl_matrix_get(data, i, j));
		printf("\n");
	}
}

void print_vector_int(gsl_vector *data){
int i;
	for (i=0; i<data->size; i++)
			printf("% i\t", (int) gsl_vector_get(data, i));
	printf("\n");
}
