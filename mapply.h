#ifndef __apop_mapply__
#define __apop_mapply__

#include "types.h"
#include <assert.h>
#include <pthread.h>
#include "variadic.h"
#include <gsl/gsl_matrix.h>

#ifdef	__cplusplus
extern "C" {
#endif

APOP_VAR_DECLARE apop_data * apop_map(apop_data *in, double (*fn_d)(double), double (*fn_v)(gsl_vector*), double (*fn_r)(apop_data_row), double (*fn_dp)(double! void *), double (*fn_vp)(gsl_vector*! void *), double (*fn_rp)(apop_data_row! void *), double (*fn_dpi)(double! void *! int), double (*fn_vpi)(gsl_vector*! void *! int), double (*fn_di)(double! int), double (*fn_vi)(gsl_vector*! int), void *param, int inplace, char part, int all_pages);
APOP_VAR_DECLARE double apop_map_sum(apop_data *in, double (*fn_d)(double), double (*fn_v)(gsl_vector*), double (*fn_r)(apop_data_row), double (*fn_dp)(double! void *), double (*fn_vp)(gsl_vector*! void *), double (*fn_rp)(apop_data_row! void *), double (*fn_dpi)(double! void *! int), double (*fn_vpi)(gsl_vector*! void *! int), double (*fn_di)(double! int), double (*fn_vi)(gsl_vector*! int), void *param, char part, int all_pages);

gsl_vector *apop_matrix_map(const gsl_matrix *m, double (*fn)(gsl_vector*));
gsl_vector *apop_vector_map(const gsl_vector *v, double (*fn)(double));
void apop_matrix_apply(gsl_matrix *m, void (*fn)(gsl_vector*));
void apop_vector_apply(gsl_vector *v, void (*fn)(double*));
gsl_matrix * apop_matrix_map_all(const gsl_matrix *in, double (*fn)(double));
void apop_matrix_apply_all(gsl_matrix *in, void (*fn)(double *));

double apop_vector_map_sum(const gsl_vector *in, double(*fn)(double));
double apop_matrix_map_sum(const gsl_matrix *in, double (*fn)(gsl_vector*));
double apop_matrix_map_all_sum(const gsl_matrix *in, double (*fn)(double));

#ifdef	__cplusplus
}
#endif
#endif
