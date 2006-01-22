#include <gsl/gsl_linalg.h>

double apop_det_and_inv(gsl_matrix *in, gsl_matrix **out, 
                                    int calc_det, int calc_inv) {
int         sign;
double      the_determinant = 0;
    gsl_matrix *invert_me = gsl_matrix_alloc(in->size1, in->size1);
    gsl_permutation * perm = gsl_permutation_alloc(in->size1);
    invert_me = gsl_matrix_alloc(in->size1, in->size1);
    gsl_matrix_memcpy (invert_me, in);
    gsl_linalg_LU_decomp(invert_me, perm, &sign);
    if (calc_inv){
        *out    = gsl_matrix_alloc(in->size1, in->size1); //square.
        gsl_linalg_LU_invert(invert_me, perm, *out);
        }
    if (calc_det)
        the_determinant = gsl_linalg_LU_det(invert_me, sign);
    gsl_matrix_free(invert_me);
    gsl_permutation_free(perm);
    return(the_determinant);
}
