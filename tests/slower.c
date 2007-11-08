#include <apop.h>

void test_score(){
  int         i, j, len = 800;
  gsl_rng     *r        = apop_rng_alloc(123);
  apop_data   *data     = apop_data_alloc(0, len,1);
  apop_model  *source   = apop_model_copy(apop_normal);
    apop_model_clear(NULL, source);
    for (i=0; i<10; i++){
        gsl_vector_set(source->parameters->vector, 0, gsl_ran_flat(r, -5, 5));
        gsl_vector_set(source->parameters->vector, 1, gsl_ran_flat(r, .01, 5));
        for (j=0; j< len; j++)
            apop_draw(gsl_matrix_ptr(data->matrix, j, 0), r, source);
        apop_mle_settings *estme = apop_mle_settings_alloc(data, apop_normal);
        estme->method =5;
        apop_model *out = apop_maximum_likelihood(data, *estme->model);
        double sigsqn = gsl_pow_2(out->parameters->vector->data[1])/len;
        assert(fabs(apop_data_get(out->covariance, 0,0)-sigsqn) < 1e-3);
        assert(fabs(apop_data_get(out->covariance, 1,1)-sigsqn/2) < 1e-3);
        assert(apop_data_get(out->covariance, 0,1) + apop_data_get(out->covariance, 0,1) < 1e-3);
        apop_model_free(out);
        printf(".");
    }
    apop_data_free(data);
    apop_model_free(source); 
}
