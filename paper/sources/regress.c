//A complete C program to read in a text file,
//run an OLS estimation, and show the result.
#include <apophenia/headers.h>

int main(){
    apop_data *dataset = apop_text_to_data("test_data", 0, 1);
    apop_estimate *est = apop_OLS.estimate(dataset, NULL);
    apop_data_show(est);
    return 0;
}
