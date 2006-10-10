#include <apophenia/headers.h>
int main(){
    apop_data *t = apop_text_to_data("data-markov", 0, 0);
    gsl_matrix *out= gsl_matrix_calloc(t->matrix->size1,t->matrix->size1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, t->matrix, t->matrix,0, out);
    apop_matrix_show(out);
    return 0;
}
