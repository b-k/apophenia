#include <apop.h>

int main(void){
apop_data       *data;
apop_model   *est;
    apop_text_to_db("data","d");
    data = apop_query_to_data("select * from d");
    est  = apop_estimate(data, apop_ols);
    apop_model_show(est);
    return 0;
}
