#include <apop.h>

#ifdef Testing
#define Show_results(m)
#else
#define Show_results(m) apop_model_print(m, NULL);
#endif

//The constraint function.
double over_zero(apop_data *in, apop_model *m){
    return apop_data_get(in) > 0;
}

//The optional scaling function.
double in_bounds(apop_model *m){
    double z = 0;
    gsl_vector_view vv = gsl_vector_view_array(&z, 1);
    return 1- apop_cdf(&((apop_data){.vector=&vv.vector}), m);
}

int main(){
    /*Set up a Normal distribution, with data truncated to be nonnegative.
      This version doesn't use the in_bounds function above, and so the
      default scaling function is used.*/
    gsl_rng *r = apop_rng_alloc(213);
    apop_model *norm = apop_model_set_parameters(apop_normal, 1.2, 0.8);
    apop_model *trunc = apop_model_dconstrain(.base_model=apop_model_copy(norm), 
                            .constraint=over_zero, .draw_ct=5e4, .rng=r);

    //make draws. Currently, you need to prep the model first.
    apop_prep(NULL, trunc);
    apop_data *d = apop_model_draws(trunc, 1e5);

    //Estimate the parameters given the just-produced data:
    apop_model *est = apop_estimate(d, trunc);
    Show_results(est);
    assert(apop_vector_distance(est->parameters->vector, norm->parameters->vector)<1e-1);

    //Generate a data set that is truncated at zero using alternate means
    apop_data *normald = apop_model_draws(apop_model_set_parameters(apop_normal, 0, 1), 5e4);
    for (int i=0; i< normald->matrix->size1; i++){
        double *d = apop_data_ptr(normald, i);
        if (*d < 0) *d *= -1;
    }

    //this time, use an unparameterized model, and the in_bounds fn
    apop_model *re_trunc = apop_model_dconstrain(.base_model=apop_normal, 
                            .constraint=over_zero, .scaling=in_bounds);

    apop_model *re_est = apop_estimate(normald, re_trunc);
    Show_results(re_est)
    assert(apop_vector_distance(re_est->parameters->vector, apop_vector_fill(gsl_vector_alloc(2), 0, 1))<1e-1);
    apop_model_free(trunc);
}
