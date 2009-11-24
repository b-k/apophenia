#include <apop.h>

int main(){
  gsl_rng *r    = apop_rng_alloc(10);
  size_t    i, ct = 1e4;
  apop_data *d  = apop_data_alloc(0,ct,2);
  double    draws[2];
  apop_data *params = apop_data_alloc(2,2,2);
  apop_data_fill(params, 8,  1, 0.5,
                         2,  0.5, 1);
  apop_model *pvm = apop_model_copy(apop_multivariate_normal);
  pvm->parameters = apop_data_copy(params);
    for(i=0; i< ct; i++){
        apop_draw(draws, r, pvm);
        apop_data_set(d, i, 0, draws[0]);
        apop_data_set(d, i, 1, draws[1]);
    }
    gsl_vector_set_all(pvm->parameters->vector, GSL_NAN);
    apop_model *mep1   = apop_model_fix_params(pvm);
    apop_model *e1  = apop_estimate(d, *mep1);
    printf("original params: ");
    apop_vector_show(params->vector);
    printf("estimated params: ");
    apop_vector_show(e1->parameters->vector);
    return 0;
}
