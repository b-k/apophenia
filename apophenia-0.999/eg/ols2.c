#include <apop.h>

int main(void){
    apop_text_to_db(.text_file="data", .tabname="d");
    apop_data *data = apop_query_to_data("select * from d");
    apop_model *est = apop_estimate(data, apop_ols);
    apop_model_show(est);

    Apop_settings_add_group(est, apop_pm, .index =1);  
    apop_model *first_param_distribution = apop_parameter_model(data, est);
    Apop_row(est->parameters, 1, param);
    double area_under_p = apop_cdf(param, first_param_distribution);
    apop_data_set(param, 0, -1, 0);
    double area_under_zero = apop_cdf(param, first_param_distribution);
    printf("reject the null for x_1 with %g percent confidence.\n",
                                 2*fabs(area_under_p-area_under_zero));
}
