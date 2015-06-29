#include <apop.h>
//This program finds the p-value of a K-S test between
//500 draws from a N(0, 1) and a N(x, 1), where x grows from 0 to 1.

apop_model * model_to_pmfs(apop_model *m1, int size){
    apop_data *outd1 = apop_model_draws(m1, size);
    return apop_estimate(apop_data_sort(outd1), apop_pmf);
}

int main(){
    apop_model *n1 = apop_model_set_parameters(apop_normal, 0, 1);
    apop_model *pmf1 = model_to_pmfs(n1, 5e2);
    apop_data *ktest;

    //first, there should be zero divergence between a PMF and itself:
    apop_model *pmf2 = apop_model_copy(pmf1);
    ktest = apop_test_kolmogorov(pmf1, pmf2);
    double pval = apop_data_get(ktest, .rowname="p value, 2 tail");
    assert(pval > .999);

    //as the mean m drifts, the pval for a comparison
    //between a N(0, 1) and N(m, 1) gets smaller.
    printf("mean\tpval\n");
    double prior_pval = 18;
    for(double i=0; i<= .6; i+=0.2){
        apop_model *n11 = apop_model_set_parameters(apop_normal, i, 1);
        ktest = apop_test_kolmogorov(pmf1, n11);
        apop_data_print(ktest, NULL);
        double pval = apop_data_get(ktest, .rowname="p value, 2 tail");
        assert(pval < prior_pval);
        printf("%g\t%g\n", i, pval);
        prior_pval = pval;
    }
    apop_model_free(pmf1);
}
