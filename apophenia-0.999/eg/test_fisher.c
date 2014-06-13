#include <apop.h>

int main() {
    /* This test is thanks to Nick Eriksson, who sent it to me in the form of a bug report. */
    apop_data * testdata = apop_data_falloc((2, 3),
                              30, 50, 45, 
                              34, 12, 17 );
    apop_data * t2 = apop_test_fisher_exact(testdata);
    assert(fabs(apop_data_get(t2,1) - 0.0001761) < 1e-6);
}
