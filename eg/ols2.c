#include <apop.h>
#include <unistd.h>

int main(void){
    char *datafile = (access("ss08pdc.csv", R_OK)!=-1) ? "ss08pdc.csv" : "data";
    apop_text_to_db(.text_file=datafile, .tabname="dc");
    apop_data *data = apop_query_to_data("select log(pincp+10), agep, sex "
                                    "from dc where agep+ pincp+sex is not null and pincp>=0");
    apop_model *est = apop_estimate(data, apop_ols);
    apop_model_show(est);

    Apop_settings_add_group(est, apop_pm, .index =1);  
    apop_model *first_param_distribution = apop_parameter_model(data, est);

    Apop_row(est->parameters, 1, param);
    double area_under_p = apop_cdf(param, first_param_distribution);

    apop_data_set(param, 0, -1, .val=0);
    double area_under_zero = apop_cdf(param, first_param_distribution);
    printf("reject the null for agep with %g percent confidence.\n",
                                 100*(2*fabs(area_under_p-area_under_zero)));
}
