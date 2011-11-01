#include <apop.h>

int main(){
    gsl_rng *r    = apop_rng_alloc(10);
    size_t i, ct = 5e4;

    //set up the model & params
    apop_data *d  = apop_data_alloc(ct,2);
    apop_data *params = apop_data_alloc(2,2,2);
    apop_data_fill(params, 8,  1, 0.5,
                           2,  0.5, 1);
    apop_model *pvm = apop_model_copy(apop_multivariate_normal);
    pvm->parameters = apop_data_copy(params);

    //make random draws from the multivar. normal
    //this `pull a row view, fill its data element' form works for rows but not cols.
    for(i=0; i< ct; i++){
        Apop_row(d, i, onerow);
        apop_draw(onerow->data, r, pvm);
    }

    //set up and estimate a model with fixed covariance matrix but free means
    gsl_vector_set_all(pvm->parameters->vector, GSL_NAN);
    apop_model *mep1   = apop_model_fix_params(pvm);
    apop_model *e1  = apop_estimate(d, *mep1);
    
    //compare results
    printf("original params: ");
    apop_vector_show(params->vector);
    printf("estimated params: ");
    apop_vector_show(e1->parameters->vector);
    return 0;
}
