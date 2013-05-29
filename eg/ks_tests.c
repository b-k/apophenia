#include <apop.h>
//This program finds the p-value of a K-S test between
//500 draws from a N(0, 1) and a N(x, 1), where x grows from 0 to 1.

//Fill a data set with draws from a model.
//This will be added to Apophenia soon.
apop_data *apop_model_draws(apop_model *m, int size, gsl_rng *r){
    apop_data *out = apop_data_alloc(size, m->dsize);
    for (int i=0; i< size; i++){
        Apop_data_row(out, i, onerow);
        apop_draw(onerow->matrix->data, r, m);
    }
    return out;
}

//Produce two models with synced PMFs.
//To do: rewrite the K-S test so that this is unnecessary
void models_to_pmfs(apop_model *m1, apop_model *m2, int size, gsl_rng *r,
        apop_model **out1, apop_model **out2){
    apop_data *outd1 = apop_model_draws(m1, size, r);
    apop_data_sort(outd1);
    apop_data *outd2 = apop_data_copy(outd1);
    outd2->weights = gsl_vector_alloc(size);
    for(int i=0; i< size; i++){
        Apop_data_row(outd1, i, onerow);
        gsl_vector_set(outd2->weights, i, apop_p(onerow, m2));
    }
    *out1= apop_estimate(outd1, apop_pmf);
    *out2= apop_estimate(outd2, apop_pmf);
}

#ifndef Testing
#define cprintf(...) printf(__VA_ARGS__)
#else
#define cprintf(...)
#endif 

int main(){
    apop_model *n1 = apop_model_set_parameters(apop_normal, 0, 1);
    apop_model *pmf1, *pmf2;
    gsl_rng *r = apop_rng_alloc(123);
    apop_data *ktest;

    //as the mean m drifts, the pval for a comparison
    //between a N(0, 1) and N(m, 1) gets smaller.
    cprintf("mean\tpval\n");
    double prior_pval = 18;
    for(double i=0; i<= 1; i+=0.2){
        apop_model *n11 = apop_model_set_parameters(apop_normal, i, 1);
        models_to_pmfs(n1, n11, 5e2, r, &pmf1, &pmf2);
        ktest = apop_test_kolmogorov(pmf1, pmf2);
        double pval = apop_data_get(ktest, .rowname="p value, 2 tail");
        assert(pval < prior_pval);
        cprintf("%g\t%g\n", i, pval);
        prior_pval = pval;
    }
}
