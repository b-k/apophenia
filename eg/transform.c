#include <apop.h>

// For defining the bounds the data-constraining function
// needs to enforce.
double greater_than_zero(apop_data *d, apop_model *m){
    return apop_data_get(d) > 0;
}

int main(){
    apop_model_print (
        apop_estimate(
             apop_update(
                apop_model_draws(
                    apop_model_mixture(
                        apop_model_set_parameters(apop_poisson, 2.8),
                        apop_model_set_parameters(apop_poisson, 2.0),
                        apop_model_set_parameters(apop_poisson, 1.3)
                    ), 
                    1e4
                ),
                apop_model_dconstrain(
                    .base_model=apop_model_set_parameters(apop_normal, 2, 1), 
                    .constraint=greater_than_zero
                ),
                apop_poisson
            )->data,
            apop_normal
        )
    );
}
