#include <apop.h>

#define print_draws(mm) apop_data_print(apop_model_draws(mm, 20),\
                                        .output_name= "draws-" #mm);

int main(){
    apop_model *uniform_20 = apop_model_set_parameters(apop_uniform, 0, 20);
    apop_data *d = apop_model_draws(uniform_20, 10);

    //Estimate a Normal distribution from the data:
    apop_model *N = apop_estimate(d, apop_normal);
    print_draws(N);

    //estimate a one-dimensional multivariate Normal from the data:
    apop_model *mvN = apop_estimate(d, apop_multivariate_normal);
    print_draws(mvN);


    //fixed parameter list:
    apop_model *std_normal = apop_model_set_parameters(apop_normal, 0, 1);
    print_draws(std_normal);

    //variable-size parameter list:
    apop_model *std_multinormal = apop_model_copy(apop_multivariate_normal);
    std_multinormal->msize1 =
    std_multinormal->msize2 =
    std_multinormal->vsize =
    std_multinormal->dsize = 3;
    std_multinormal->parameters = apop_data_falloc((3, 3, 3),
                                1,  1, 0, 0, 
                                1,  0, 1, 0,
                                1,  0, 0, 1);
    print_draws(std_multinormal);


    //estimate a KDE using the defaults:
    apop_model *k = apop_estimate(d, apop_kernel_density);
    print_draws(k);

    /*the documentation tells us that a KDE estimation consists of filling 
      an apop_kernel_density_settings group, so we can set it to use a 
      Normal(μ, 2) kernel via: */

    apop_model *k2 = apop_model_copy_set(apop_kernel_density, apop_kernel_density, 
                         .base_data=d,
                         .kernel = apop_model_set_parameters(apop_normal, 0, 2));
    print_draws(k2);
}
