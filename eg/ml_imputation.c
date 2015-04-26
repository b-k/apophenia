#include <apop.h>

static void compare_mvn_estimates(apop_model *L, apop_model *R, double tolerance){
    gsl_vector_sub(L->parameters->vector, R->parameters->vector);
    gsl_matrix_sub(L->parameters->matrix, R->parameters->matrix);
    assert(fabs(apop_sum(L->parameters->vector)) + fabs (apop_matrix_sum(L->parameters->matrix)) < tolerance);
}

void test_ml_imputation(gsl_rng *r){
    size_t len = 4e4;
    int i,j;
    apop_data *fillme = apop_data_alloc(len, 3);
    apop_model *mvn = apop_model_copy(apop_multivariate_normal);
    mvn->parameters = apop_data_alloc(3, 3, 3);
    for(i=0; i < 3; i ++)
        for(j=-1; j < 3; j ++)
            apop_data_set(mvn->parameters, i, j, gsl_rng_uniform(r));
    //now make your random garbage symmetric
    for(i=0; i < 3; i ++)
        for(j=i+1; j < 3; j ++)
            apop_data_set(mvn->parameters, j, i, apop_data_get(mvn->parameters, i, j));
    apop_matrix_to_positive_semidefinite(mvn->parameters->matrix);
    apop_model_draws(mvn, .draws=fillme);
    //apop_data_show(mvn->parameters);
    apop_model *est = apop_estimate(fillme, apop_multivariate_normal);
    //apop_data_show(est->parameters);
    compare_mvn_estimates(est, mvn, 1e-1);

    double pct_to_delete = 0.01;
    int max_to_delete = 7, ctr = 0;
    for(i=0; i < len && ctr < max_to_delete; i ++)
        for(j=0; j < 3; j ++)
            if (gsl_rng_uniform(r) < pct_to_delete){
                apop_data_set(fillme, i, j, GSL_NAN);
                ctr++;
            }
    apop_ml_imputation(fillme, mvn); 
    apop_model *est2 = apop_estimate(fillme, apop_multivariate_normal);
    //apop_data_show(est2->parameters);
    compare_mvn_estimates(est2, mvn, 1e-1);
    apop_data_free(fillme);
}

int main(){
    test_ml_imputation(apop_rng_alloc(42));
}
