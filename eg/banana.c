#include <apop.h>

typedef struct {
    double scaling;
} coeff_struct;

long double banana (double *params, coeff_struct *in){
    return (gsl_pow_2(1-params[0]) + in->scaling*gsl_pow_2(params[1]-gsl_pow_2(params[0])));
}

long double ll (apop_data *d, apop_model *in){
    return - banana(in->parameters->vector->data, in->more);
}

int main(){
    coeff_struct co = {.scaling=100};
    apop_model *b = &(apop_model) {"¡Bananas!", .log_likelihood= ll, .vsize=2, 
                                .more = &co, .more_size=sizeof(coeff_struct)};
    Apop_model_add_group(b, apop_mle, .verbose='y', .method="NM simplex");
    Apop_model_add_group(b, apop_parts_wanted);
    apop_model *e1 = apop_estimate(NULL, b);
    apop_model_print(e1, NULL);

    Apop_settings_set(b, apop_mle, method, "BFGS cg");
    apop_model *e2 = apop_estimate(NULL, b);
    apop_model_print(e2, NULL);

    gsl_vector *one = apop_vector_fill(gsl_vector_alloc(2), 1, 1);
    assert(apop_vector_distance(e1->parameters->vector, one) < 1e-2);
    assert(apop_vector_distance(e2->parameters->vector, one) < 1e-2);
}
