#include <apop.h>

/* This replacement for apop_model_show(in) demonstrates retrieval of the useful
settings: the weights (λ) and list of estimated models. */
void show_mix(apop_model *in){
    apop_mixture_settings *ms = Apop_settings_get_group(in, apop_mixture);
    printf("The weights:\n");
    apop_vector_print(ms->weights);
    printf("\nThe models:\n");
    for (apop_model **m = ms->model_list; *m; m++) //model_list is a NULL-terminated list.
        apop_model_print(*m, stdout);
}


int main(){
    apop_text_to_db("faith.data", "ff");
    apop_data *dd = apop_query_to_data("select waiting from ff");
    apop_model *mf = apop_model_mixture(apop_model_copy(apop_normal), apop_model_copy(apop_normal));

    /* The process is famously sensitive to starting points. Try many random points, or
       eyeball the distribution's plot and guess at the starting values. */
    Apop_model_add_group(mf, apop_mle, .starting_pt=(double[]){50, 5, 80, 5}, 
                                       .step_size=3, .tolerance=1e-6);
    apop_model *mfe = apop_estimate(dd, mf);
    apop_model_print(mfe, stdout);
    printf("LL=%g\n", apop_log_likelihood(dd, mfe));


    printf("\n\nValues calculated in the source paper, for comparison.\n");
    apop_model *r_ed = apop_model_mixture(
                         apop_model_set_parameters(apop_normal, 54.61364, 5.869089),
                         apop_model_set_parameters(apop_normal, 80.09031, 5.869089));
    apop_data *wts = apop_data_falloc((2), 0.3608498, 0.6391502);
    Apop_settings_add(r_ed, apop_mixture, weights, wts->vector);
    show_mix(r_ed);
    printf("LL=%g\n", apop_log_likelihood(dd, r_ed));
}
