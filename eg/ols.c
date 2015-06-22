#include <apop.h>

int main(){
    apop_text_to_db(.text_file="data", .tabname="d");
    apop_data *data = apop_query_to_data("select * from d");
    apop_model *est = apop_estimate(data, apop_ols);
    apop_model_print(est, NULL);
}
