#include <apop.h>

// Find the μ/σ  of a set of 10 draws from a Uniform(-1, 1)
void sim_step(apop_data *none, apop_model *m){
    int sub_draws = 20;
    static apop_model *unif;
    if (!unif) unif = apop_model_set_parameters(apop_uniform, -1, 1);
    apop_data *draws= apop_model_draws(unif, sub_draws);

    apop_data_set(m->parameters, 0, .val=apop_mean(Apop_cv(draws, 0)));
    apop_data_set(m->parameters, 1, .val=sqrt(apop_var(Apop_cv(draws, 0))));
    apop_data_add_names(m->parameters, 'r', "μ", "σ");
    apop_data_free(draws);
}

apop_model *clt_sim = &(apop_model){.name="CLT simulation", .vsize=2, .estimate=sim_step};

int main(){
    apop_data *boots;
    apop_data * boot_cov = apop_bootstrap_cov(NULL, clt_sim, .iterations=1000, .boot_store=&boots);
    apop_data_print(boot_cov);
    apop_data *means = Apop_c(boots, 0);

    printf("\nStats via Normal model:\n");
    apop_data *np = apop_estimate(means, apop_normal)->parameters;
    np->more = NULL; //rm covariance of statistics.
    apop_data_print(np);

    //σ from the Normal should == sqrt(cov(μ_boot))
    assert(fabs(sqrt(apop_data_get(boot_cov,0,0)) - apop_data_get(np, 1)) < 1e-4);
}
