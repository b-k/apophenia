#include <apop.h>
int main(){ apop_model_print(apop_estimate(apop_text_to_data("data"), apop_ols), NULL); }
