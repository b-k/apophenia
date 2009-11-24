#include <apop.h>

int main(void){ 
  apop_data       *data; 
  apop_model   *est;
    apop_db_open(NULL);
    apop_text_to_db("data", "d",.has_row_names= 0,.has_col_names=1);
    data       = apop_query_to_data("select * from d");
    est        = apop_estimate(data, apop_ols);
    printf("The OLS coefficients:\n");
    apop_model_print(est);
} 
