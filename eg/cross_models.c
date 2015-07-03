#include <apop.h>

/* In this initial example, build a cross product of two Normal(2,.1) distributions.
Make 10,000 draws from it.
 
Then, build a cross product of two unparameterized Normals and estimate the parameters
of the combined model; check that they match the (2, .1) we started with.
*/
void cross_normals(){
    double mu = 2;
    double sigma = .1;
    apop_model *n1 = apop_model_set_parameters(apop_normal, mu, sigma);
    apop_model *n2 = apop_model_copy(n1);
    apop_model *two_independent_normals = apop_model_cross(n1, n2);
    //
    //We don't use it, but the cross product of three is just as easy:
    apop_model *n3 = apop_model_copy(n1);
    apop_model *three_independent_normals = apop_model_cross(n1, n2, n3);

    apop_data *draws = apop_model_draws(two_independent_normals, .count=10000);

    //The unparameterized cross product:
    apop_model *two_n = apop_model_cross(
                    apop_model_copy(apop_normal),
                    apop_model_copy(apop_normal)
                    );
    apop_model *estimated_norms = apop_estimate(draws, two_n);

    apop_model_print(estimated_norms);
    apop_data *estp1 = Apop_settings_get(estimated_norms, apop_cross, model1)->parameters;
    apop_data *estp2 = Apop_settings_get(estimated_norms, apop_cross, model2)->parameters;
    assert(fabs(apop_data_get(estp1, 0) - mu)    < 2e-3);
    assert(fabs(apop_data_get(estp2, 0) - mu)    < 2e-3);
    assert(fabs(apop_data_get(estp1, 1) - sigma) < 2e-3);
    assert(fabs(apop_data_get(estp2, 1) - sigma) < 2e-3);
}

//bind together a Poisson and a Normal
void norm_cross_poisson(){
    apop_model *m1 = apop_model_set_parameters(apop_poisson, 3);
    apop_model *m2 = apop_model_set_parameters(apop_normal, -5, 1);
    apop_model *mm = apop_model_cross(m1, m2);
    int len = 1e5;
    apop_data *draws = apop_model_draws(mm, len);
    for (int i=0; i< len; i++){
        Apop_row_v(draws, i, onev);
        assert((int)onev->data[0] == onev->data[0]);
        assert(onev->data[1]<0);
    }

    /*The rest of the test script recovers the parameters.
    Input data to an apop_cross model can take two formats. In cross_normals, the
    draws are in a single matrix. Here, the data for the Poisson (col 0 of the draws)
    will be put in an apop_data set, and the data for the Normal (col 1 of the draws)
    on a second page appended to the first. Then, set the .splitpage element of the
    apop_cross settings group to the name of the second page.
    */
    apop_data *comeback = apop_data_alloc();
    comeback->vector = apop_vector_copy(Apop_cv(draws, 0));
    apop_data_add_page(comeback, apop_data_alloc(), "p2");
    comeback->more->vector = apop_vector_copy(Apop_cv(draws, 1));

    //set up the un-parameterized crossed model, including
    //the name at which to split the data set
    apop_model *estme = apop_model_cross(apop_model_copy(apop_poisson), apop_model_copy(apop_normal));
    Apop_settings_add(estme, apop_cross, splitpage, "p2");
    apop_model *ested = apop_estimate(comeback, estme);

    //test that the parameters are as promised.
    apop_model *m1back = apop_settings_get(ested, apop_cross, model1);
    apop_model *m2back = apop_settings_get(ested, apop_cross, model2);
    assert(fabs(apop_data_get(m1back->parameters, .col=-1) - 3) < 5e-1);
    assert(fabs(apop_data_get(m2back->parameters, .col=-1) - -5) < 5e-1);
    assert(fabs(apop_data_get(m2back->parameters, .col=-1, .row=1) - 1) < 5e-1);

    //You can cross as many models as you'd like.
    apop_model *m3 = apop_model_set_parameters(apop_poisson, 8);
    apop_model *mmm = apop_model_cross(m1, m2, m3);
    apop_data *sum = apop_data_summarize(apop_model_draws(mmm, 1e5));
    assert(fabs(apop_data_get(sum, .row=0, .colname="mean") - 3) < 2e-2);
    assert(fabs(apop_data_get(sum, .row=1, .colname="mean") - -5) < 2e-2);
    assert(fabs(apop_data_get(sum, .row=2, .colname="mean") - 8) < 4e-2);
    assert(apop_data_get(sum, .row=0, .colname="median") == 3);
    assert(apop_data_get(sum, .row=2, .colname="median") == 8);
}

int main(){
    cross_normals();
    norm_cross_poisson();
}
