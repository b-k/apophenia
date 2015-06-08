#include <apop.h>

#define Diff(L, R, eps) {double left=(L), right=(R); Apop_stopif(isnan(left-right) || fabs((left)-(right))>(eps), abort(), 0, "%g is too different from %g (abitrary limit=%g).", (double)(left), (double)(right), eps);}

/** I claim that the F test calculated via apop_F_test(est, NULL, NULL)
 equals a transformation of R^2 (after a normalization step).
*/
void test_f(apop_model *est){
    apop_data *rsq  = apop_estimate_coefficient_of_determination(est);
    apop_data *constr= apop_data_calloc(est->parameters->vector->size-1, est->parameters->vector->size);
    int i;
    for (i=1; i< est->parameters->vector->size; i++)
        apop_data_set(constr, i-1, i, 1);
    apop_data *ftab = apop_F_test(est, constr);
    apop_data *ftab2 = apop_F_test(est, NULL);
    //apop_data_show(ftab);
    //apop_data_show(ftab2);
    double n = est->data->matrix->size1;
    double K = est->parameters->vector->size-1;
    double r = apop_data_get(rsq, .rowname="R squared");
    double f = apop_data_get(ftab, .rowname="F statistic");
    double f2 = apop_data_get(ftab2, .rowname="F statistic");
    Diff (f , r*(n-K)/((1-r)*K) , 1e-3);
    Diff (f2 , r*(n-K)/((1-r)*K) , 1e-3);
}

int main(){
    apop_data *d = apop_text_to_data("test_data2");
    apop_model *an_ols_model = apop_model_copy(apop_ols);
    Apop_model_add_group(an_ols_model, apop_lm, .want_expected_value= 1);
    apop_model *e  = apop_estimate(d, an_ols_model);
    test_f(e);
}
