#include <apop.h>


int main(){
    //Set up an apop_data set with only one number.
    //Most of these functions will only look at the first data point encountered.
    apop_data *onept = apop_data_alloc(1, 1);
    apop_data_set(onept, .val=23);

    apop_model *norm = apop_model_set_parameters(apop_normal, 23, 138.8);
    double val = apop_cdf(onept, norm);
    assert(fabs(val - 0.5) < 1e-4);

    double tolerance = 1e-8;
    //Macroizing the sample routine above:
    #define model_val_cdf(model, value, cdf_result) {   \
        apop_data_set(onept, .val=(value));             \
        assert(fabs((apop_cdf(onept, model))-(cdf_result))< tolerance);   \
    }

    apop_model *uni = apop_model_set_parameters(apop_uniform, 20, 26);
    model_val_cdf(uni, 0, 0);
    model_val_cdf(uni, 20, 0);
    model_val_cdf(uni, 21, 1./6);
    model_val_cdf(uni, 23, 0.5);
    model_val_cdf(uni, 25, 5./6);
    model_val_cdf(uni, 26, 1);
    model_val_cdf(uni, 260, 1);

    //Improper uniform always returns 1/2.
    model_val_cdf(apop_improper_uniform, 0, 0.5);
    model_val_cdf(apop_improper_uniform, 228, 0.5);
    model_val_cdf(apop_improper_uniform, INFINITY, 0.5);

    apop_model *binom = apop_model_set_parameters(apop_binomial, 2001, 0.5);
    model_val_cdf(binom, 0, 0);
    model_val_cdf(binom, 1000, .5);
    model_val_cdf(binom, 2000, 1);

    apop_model *bernie = apop_model_set_parameters(apop_bernoulli, 0.75);
    //p(0)=.25; p(1)=.75; that determines the CDF.
    //Notice that the CDF's integral is over a closed interval.
    model_val_cdf(bernie, -1, 0);
    model_val_cdf(bernie, 0, 0.25);
    model_val_cdf(bernie, 0.1, 0.25);
    model_val_cdf(bernie, .99, 0.25);
    model_val_cdf(bernie, 1, 1);
    model_val_cdf(bernie, INFINITY, 1);

    //alpha=beta -> symmetry
    apop_model *beta = apop_model_set_parameters(apop_beta, 2, 2);
    model_val_cdf(beta, -INFINITY, 0);
    model_val_cdf(beta, 0.5, 0.5);
    model_val_cdf(beta, INFINITY, 1);

    //This beta distribution -> uniform
    apop_model *beta_uni = apop_model_set_parameters(apop_beta, 1, 1);
    model_val_cdf(beta_uni, 0, 0);
    model_val_cdf(beta_uni, 1./6, 1./6);
    model_val_cdf(beta_uni, 0.5, 0.5);
    model_val_cdf(beta_uni, 1, 1);


    beta_uni->cdf = NULL; //With no closed-form CDF; make random draws to estimate the CDF.
    Apop_model_add_group(beta_uni, apop_cdf, .draws=1e6); //extra draws to improve accuracy, but we have to lower our tolerance anyway.
    tolerance=1e-3;
    model_val_cdf(beta_uni, 0, 0);
    model_val_cdf(beta_uni, 1./6, 1./6);
    model_val_cdf(beta_uni, 0.5, 0.5);
    model_val_cdf(beta_uni, 1, 1);


    //sum of three symmetric distributions: still symmetric.
    apop_model *sum_of_three = apop_model_mixture(beta, apop_improper_uniform, beta_uni);
    model_val_cdf(sum_of_three, 0.5, 0.5);


    apop_data *threepts = apop_data_falloc((3,1), -1, 0, 1);
    apop_model *kernels = apop_estimate(threepts, apop_kernel_density);
    model_val_cdf(kernels, -5, 0);
    model_val_cdf(kernels, 0, 0.5);
    model_val_cdf(kernels, 10, 1);
}
