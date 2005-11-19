#include <gsl/gsl_linalg.h>

void sv_decomposition(gsl_matrix *data, int dimensions_we_want, 
		gsl_matrix * pc_space, gsl_vector *eigenvalues)
{
    gsl_matrix * eigenvectors = gsl_matrix_alloc(data->size2, data->size2);
    gsl_vector * dummy_v = gsl_vector_alloc(data->size2);
    gsl_matrix * dummy_m  = gsl_matrix_calloc(data->size1, data->size2);
    gsl_matrix_memcpy(dummy_m, data);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, data, data, 0, dummy_m);
    gsl_linalg_SV_decomp(dummy_m, eigenvectors, eigenvalues, dummy_v);
    int i;
    for (i=0; i<dimensions_we_want; i++){
        gsl_matrix_get_col(dummy_v, eigenvectors, i);
        gsl_matrix_set_col(pc_space, i, dummy_v);
    }
    gsl_vector_free(dummy_v);
    gsl_matrix_free(dummy_m);
}
