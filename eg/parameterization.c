/*Parameterizing or initializing a model

Apophenia ships with several models, which have the requisite procedures for estimation,
making draws, and so on, but have <tt>params==NULL</tt> and <tt>settings==NULL</tt>. The
model is thus, for many purposes, incomplete, and you will need to take some action to
complete the model. There are several possibilities:



\li Estimate it! Almost all models can be sent with a data set as an argument to the
<tt>apop_estimate</tt> function. The input model is unchanged, but the output model has
parameters and settings in place.

\li If your model has a fixed number of numeric parameters, then you can set them with
\ref apop_model_set_parameters.

\li If your model has a variable number of parameters, you can directly set the \c
parameters element via \c apop_data_falloc.  For most purposes, you will also need to
set the \c msize1, \c msize2, \c vsize, and \c dsize elements to the size you want. See
the example below.

\li Some models have disparate, non-numeric settings rather than a simple matrix of
parameters. For example, an kernel density estimate needs a model as a kernel and a base data set, which can be set via \ref apop_model_copy_set.

Here is an example that shows the options for parameterizing a model. After each parameterization, 20 draws are made and written to a file named draws-[modelname].
*/




#include <apop.h>

#define print_draws(mm) apop_data_print(apop_model_draws(mm, 20), .output_name= "draws-" #mm);

apop_data *draw_some_data(){
    apop_model *uniform_0_20 = apop_model_set_parameters(apop_uniform, 0, 20);
    return apop_model_draws(uniform_0_20, 10);
}

int main(){
    apop_data *d = draw_some_data();    

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
