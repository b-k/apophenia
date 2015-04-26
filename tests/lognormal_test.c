#include <apop.h>
#define Diff(L, R, eps) {double left=(L), right=(R); Apop_stopif(isnan(left-right) || fabs((left)-(right))>(eps), abort(), 0, "%g is too different from %g (abitrary limit=%g).", (double)(left), (double)(right), eps);}

void test_lognormal(gsl_rng *r){
    apop_model *source = apop_model_copy(apop_normal);
    apop_model_clear(NULL, source);
    double mu = gsl_ran_flat(r, -1, 1);
    double sigma = gsl_ran_flat(r, .01, 1);
    int n = gsl_ran_flat(r,1,8e5);
    apop_data *data = apop_data_alloc(0,1,n);
    gsl_vector_set(source->parameters->vector, 0, mu);
    gsl_vector_set(source->parameters->vector, 1, sigma);
    for (int j=0; j< n; j++){
        double *k   = gsl_matrix_ptr(data->matrix, 0, j);
        apop_draw(k, r, source);
        *k = exp(*k);
    }
    apop_model_free(source);
    apop_model *out = apop_estimate(data, apop_lognormal);
    double muhat = apop_data_get(out->parameters, 0,-1);
    double sigmahat = apop_data_get(out->parameters, 1,-1);
    //if (verbose) printf("mu: %g, muhat: %g, var: %g, varhat: %g\n", mu, muhat,  sigma,sigmahat);
    Diff(mu, muhat, 1e-2);
    Diff(sigma, sigmahat, 1e-2);
    apop_model_free(out);

    apop_model *for_mle= apop_model_copy(apop_lognormal);
    for_mle->estimate=NULL;
    apop_model *out2 = apop_estimate(data, for_mle);
    apop_model_free(for_mle);
    muhat = apop_data_get(out2->parameters, 0,-1);
    sigmahat = apop_data_get(out2->parameters, 1,-1);
    Diff(mu, muhat, 1e-2);
    Diff(sigma, sigmahat, 1e-2);
    apop_model_free(out2);
    apop_data_free(data);
}

int main(){ test_lognormal(apop_rng_alloc(24)); }
