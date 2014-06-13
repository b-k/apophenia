#include <apop.h>

long double pack_p (apop_data *d, apop_model *m){
    double loss = 0;
    gsl_vector *v = apop_data_pack(m->parameters, NULL, 'y');
    int i;
    for (i=0; i< v->size; i++)
        loss += fabs(v->data[i] - i);
    gsl_vector_free(v);
    return 1/(1+loss);
}

void pack_prep(apop_data *d, apop_model *m){
    m->parameters = apop_data_alloc(0, 2, 2);
    apop_data_add_page(m->parameters, apop_data_alloc(0, 2, 2), "page two");
    if (!Apop_settings_get_group(m, apop_mle))
        Apop_model_add_group(m, apop_mle, .tolerance=1e-6, .step_size=3);
}

long double pack_constraint(apop_data *d, apop_model *m){
    return apop_linear_constraint(apop_data_pack(m->parameters, .all_pages='y'))*1e-5;
    //penalty size must be smaller than p().
}

apop_model *pack_counter = &(apop_model){"Optimum is that each element equals its pack order", .p = pack_p, 
                .prep=pack_prep, .constraint = pack_constraint };

int main(){
    apop_model *list = apop_estimate(NULL, pack_counter);
    #ifndef Testing
    apop_data_show(list->parameters);
    printf("%g", fabs( 1- 1/apop_p(NULL, list)));
    #endif
    assert(fabs( 1- 1/apop_p(NULL, list))< 4e-2); //lousy.
}
