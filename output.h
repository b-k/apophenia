#include <gsl/gsl_matrix.h>
#include "db.h"
#include "types.h"
#include "stats.h"
#include "linear_algebra.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

void apop_plot_line_and_scatter(apop_data *data, apop_model *est, char *);
void apop_plot(gsl_matrix *data, char plot_type, int delay);
//void apop_plot_histogram(gsl_vector *data, size_t bin_ct, char *outfile, char *preface, ...) __attribute__ ((format (printf,4,5)));
void apop_plot_histogram(gsl_vector *data, size_t bin_ct, char *outfile);
void apop_histogram_print(apop_model *h, char *outfile);
void apop_plot_lattice(apop_data *d, char filename[]);
void apop_plot_qq(gsl_vector *v, apop_model m, char *outfile);

void apop_matrix_print(gsl_matrix *data, char *file);
void apop_vector_print(gsl_vector *data, char *file);
void apop_data_print(apop_data *data, char *file);

void apop_matrix_show(const gsl_matrix *data);
void apop_vector_show(const gsl_vector *data);
void apop_data_show(const apop_data *data);
__END_DECLS
