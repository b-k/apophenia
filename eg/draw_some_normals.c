#include <apop.h>
#include <time.h>

int main(){
    apop_opts.rng_seed = time(NULL);
    apop_data_print(
            apop_model_draws(
                apop_model_set_parameters(apop_normal, 0, 1), 
                .count=10, 
            )
    );
}
