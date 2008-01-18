#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "types.h"

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

double      apop_det_and_inv(const gsl_matrix *in, gsl_matrix **out, int calc_det, int calc_inv);
gsl_matrix *apop_matrix_inverse(const gsl_matrix *in) ;
double      apop_matrix_determinant(const gsl_matrix *in) ;
void apop_normalize_for_svd(gsl_matrix *in);
apop_data*  apop_sv_decomposition(gsl_matrix *data, int dimensions_we_want);
inline void apop_vector_increment(gsl_vector * v, int i, double amt);
inline void apop_matrix_increment(gsl_matrix * m, int i, int j, double amt);
gsl_vector *apop_vector_stack(gsl_vector *v1, gsl_vector * v2);
gsl_matrix *apop_matrix_stack(gsl_matrix *m1, gsl_matrix * m2, char posn);
gsl_matrix *apop_matrix_rm_columns(gsl_matrix *in, int *drop);
int         apop_vector_bounded(gsl_vector *in, long double max);
apop_data * apop_dot(const apop_data *d1, const apop_data *d2, ...);
void        apop_vector_log(gsl_vector *v);
void        apop_vector_log10(gsl_vector *v);
void        apop_vector_exp(gsl_vector *v);
__END_DECLS
