#include <apop.h>

int main(void){
    apop_text_to_db(.text_file="data", .tabname="d");
    apop_data *data = apop_query_to_data("select * from d");
    apop_model *est = apop_estimate(data, apop_ols);
    apop_model_show(est);

    apop_model *first_param_distribution = apop_parameter_model(data, est);
//    apop_data *param = apop_data_alloc(0, 1, 1);
    //apop_data_set(param, 0, 0, apop_data_get(est->parameters, 0, -1));
    Apop_data_row(est->parameters, 1, param);
    double area_under_p = apop_cdf(param, first_param_distribution);
    printf("reject the null with %g percent confidence.\n", 2*(area_under_p-0.5));
}
