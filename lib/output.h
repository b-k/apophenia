#include <gsl/gsl_matrix.h>
#include <apophenia/db.h>
#include <apophenia/types.h>
#include <apophenia/stats.h>
#include <apophenia/linear_algebra.h>

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

void apop_plot_line_and_scatter(apop_data *data, apop_estimate *est, char *);
void apop_plot(gsl_matrix *data, char plot_type, int delay);
void apop_plot_histogram(gsl_vector *data, size_t bin_ct, char *outfile);
void apop_plot_lattice(apop_data *d, char filename[]);

void apop_matrix_print(gsl_matrix *data, char *file);
void apop_matrix_print_int(gsl_matrix *data, char *file);
void apop_vector_print(gsl_vector *data, char *file);
void apop_vector_print_int(gsl_vector *data, char *file);
void apop_data_print(apop_data *data, char *file);
void apop_data_print_int(apop_data *data, char *file);

void apop_matrix_show(gsl_matrix *data);
void apop_matrix_show_int(gsl_matrix *data);
void apop_vector_show(gsl_vector *data);
void apop_vector_show_int(gsl_vector *data);
void apop_data_show(apop_data *data);
void apop_data_show_int(apop_data *data);
__END_DECLS
