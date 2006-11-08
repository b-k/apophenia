#include <apophenia/headers.h>
int main(){
    int         i, test_ct  = 10000;
    double      data[]      = { 30, 86,
                       24, 38 };
    apop_data   *testdata   = apop_line_to_data(data, 2, 2);
    for (i = 0; i< test_ct; i++)
        apop_test_fisher_exact(testdata);
    return 0;
}
