#ifndef __apop_mapply__
#define __apop_mapply__

#include "types.h"
#include <assert.h>
#include <pthread.h>
#include "variadic.h"
#include <gsl/gsl_matrix.h>

#undef __BEGIN_DECLS    /* extern "C" stuff cut 'n' pasted from the GSL. */
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
__BEGIN_DECLS

APOP_VAR_DECLARE apop_data * apop_map(apop_data *in, double (*fn_v)(gsl_vector*), double (*fn_d)(double), double (*fn_vp)(gsl_vector*! void *), double (*fn_dp)(double! void *), double (*fn_vpi)(gsl_vector*! void *! size_t), double (*fn_dpi)(double! void *! size_t), double (*fn_vi)(gsl_vector*! size_t), double (*fn_di)(double! size_t), void *param, int inplace, char part);
APOP_VAR_DECLARE double apop_map_sum(apop_data *in, double (*fn_v)(gsl_vector*), double (*fn_d)(double), double (*fn_vp)(gsl_vector*! void *), double (*fn_dp)(double! void *), double (*fn_vpi)(gsl_vector*! void *! size_t), double (*fn_dpi)(double! void *! size_t), double (*fn_vi)(gsl_vector*! size_t), double (*fn_di)(double! size_t), void *param, char part);

gsl_vector *apop_matrix_map(const gsl_matrix *m, double (*fn)(gsl_vector*));
gsl_vector *apop_vector_map(const gsl_vector *v, double (*fn)(double));
void apop_matrix_apply(gsl_matrix *m, void (*fn)(gsl_vector*));
void apop_vector_apply(gsl_vector *v, void (*fn)(double*));
gsl_matrix * apop_matrix_map_all(const gsl_matrix *in, double (*fn)(double));
void apop_matrix_apply_all(gsl_matrix *in, void (*fn)(double *));

double apop_vector_map_sum(const gsl_vector *in, double(*fn)(double));
double apop_matrix_map_sum(const gsl_matrix *in, double (*fn)(gsl_vector*));
double apop_matrix_map_all_sum(const gsl_matrix *in, double (*fn)(double));

__END_DECLS
#endif
