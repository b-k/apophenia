#include <gsl/gsl_linalg.h>

void function invert_matrix(gsl_matrix *in, gsl_matrix *out)
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
