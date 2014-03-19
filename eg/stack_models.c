#include <apop.h>

int main(){
    //bind together a Poisson and a Normal;
    //make a draw producing a 2-element vector
    apop_model *m1 = apop_model_set_parameters(apop_poisson, 3);
    apop_model *m2 = apop_model_set_parameters(apop_normal, -5, 1);
    apop_model *mm = apop_model_stack(m1, m2);
    int len = 1e5;
    apop_data *draws = apop_model_draws(mm, len);
    for (int i=0; i< len; i++){
        Apop_row_v(draws, i, onev);
        assert((int)onev->data[0] == onev->data[0]);
        assert(onev->data[1]<0);
    }

    //The rest of the test script recovers the parameters.
    //First, set up a two-page data set: poisson data on p1, Normal on p2:
    apop_data *comeback = apop_data_alloc();
    Apop_col_v(draws, 0,fishdraws)
    comeback->vector = apop_vector_copy(fishdraws);
    apop_data_add_page(comeback, apop_data_alloc(), "p2");
    Apop_col_v(draws, 1, meandraws)
    comeback->more->vector = apop_vector_copy(meandraws);

    //set up the un-parameterized stacked model, including
    //the name at which to split the data set
    apop_model *estme = apop_model_stack(apop_model_copy(apop_poisson), apop_model_copy(apop_normal));
    Apop_settings_add(estme, apop_stack, splitpage, "p2");
    apop_model *ested = apop_estimate(comeback, estme);

    //test that the parameters are as promised.
    apop_model *m1back = apop_settings_get(ested, apop_stack, model1);
    apop_model *m2back = apop_settings_get(ested, apop_stack, model2);
    assert(fabs(apop_data_get(m1back->parameters, .col=-1) - 3) < 5e-1);
    assert(fabs(apop_data_get(m2back->parameters, .col=-1) - -5) < 5e-1);
    assert(fabs(apop_data_get(m2back->parameters, .col=-1, .row=1) - 1) < 5e-1);

    //You can stack as many models as you'd like.
    apop_model *m3 = apop_model_set_parameters(apop_poisson, 8);
    apop_model *mmm = apop_model_stack(m1, m2, m3);
    apop_data *sum = apop_data_summarize(apop_model_draws(mmm, 1e5));
    assert(fabs(apop_data_get(sum, .row=0, .colname="mean") - 3) < 2e-2);
    assert(fabs(apop_data_get(sum, .row=1, .colname="mean") - -5) < 2e-2);
    assert(fabs(apop_data_get(sum, .row=2, .colname="mean") - 8) < 4e-2);
    assert(apop_data_get(sum, .row=0, .colname="median") == 3);
    assert(apop_data_get(sum, .row=2, .colname="median") == 8);
}
