#include <apop.h>

/* A beta example whose syntax and method is under development, under the presumption that
  a sketch is better than no examples at all.
 
 Uses \ref apop_model_mixture to generate a hump-filled distribution, with the tallest hump at (3,3).

 Part of Apophenia's test suite.
 */

//Produce a 2-D multivariate normal model with unit covariance and given mean 
apop_model *produce_fixed_mvn(double x, double y){
    apop_model *out = apop_model_copy(apop_multivariate_normal);
    out->parameters = apop_data_falloc((2, 2, 2),
                        x, 1, 0,
                        y, 0, 1);
    out->dsize = 2;
    return out;
}

int main(){
    //here's a mean/covariance matrix for a standard multivariate normal.
    apop_model *many_humps = apop_model_mixture(
                        produce_fixed_mvn(5, 6),
                        produce_fixed_mvn(-5, -4),
                        produce_fixed_mvn(0, 1));
    apop_prep(NULL, many_humps);


    gsl_rng *r = apop_rng_alloc(21);
    int len = 100000;
    apop_data *d = apop_model_draws(many_humps, len, r);
    Apop_matrix_col(d->matrix, 0, first);
    #ifndef Testing
    printf("mu=%g\n", apop_mean(first));
    #endif
    assert(fabs(apop_mean(first)- 0) < 5e-2);
    Apop_matrix_col(d->matrix, 1, second);
    #ifndef Testing
    printf("mu=%g\n", apop_mean(second));
    #endif
    assert(fabs(apop_mean(second)- 1) < 5e-2);

 /*  Abuse the ML imputation routine to search for the input value with the highest
 log likelihood.
  Do the search via simulated annealing. */

    apop_data *x = apop_data_alloc(1,2);
    gsl_matrix_set_all(x->matrix, NAN);

    apop_opts.stop_on_warning='v';
    apop_model_add_group(many_humps, apop_mle, .n_tries=20, .iters_fixed_T=10, .k=3,  .method=APOP_SIMAN);
    apop_ml_impute(x, many_humps);
    #ifndef Testing
    printf("Optimum found at:\n");
    apop_data_show(x);
    #endif
    assert(fabs(apop_data_get(x, .col=0)- 0) + fabs(apop_data_get(x, .col=1) - 1) < 1e-2);
}
