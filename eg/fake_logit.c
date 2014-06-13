#include <apop.h>
#include <unistd.h>

char *testfile = "logit_test_data";

//generate a fake data set.
//Notice how the first column is the outcome, just as with standard regression.
void write_data(){
    FILE *f = fopen(testfile, "w");
    fprintf(f, "\
        outcome,A, B \n\
        0, 0, 0     \n\
        1, 1, 1     \n\
        1, .7, .5   \n\
        1, .7, .3   \n\
        1, .3, .7   \n\
        \n\
        1, .5, .5   \n\
        0, .4, .4   \n\
        0, .3, .4   \n\
        1, .1, .3   \n\
        1, .3, .1   ");
    fclose(f);
}

int main(){
    write_data();
    apop_data *d = apop_text_to_data(testfile);
    Apop_model_add_group(apop_logit, apop_mle, .tolerance=1e-5);
    apop_model *est = apop_estimate(d, apop_logit);
    unlink(testfile);

    /* Apophenia's test suite checks that this code produces 
       values close to canned values. As a human, you probably 
       just want to print the results to the screen. */
    #ifndef Testing
        apop_model_show(est);
    #else
        assert(fabs(apop_data_get(est->parameters, .rowname="1")- -1.155026) < 1e-6);
        assert(fabs(apop_data_get(est->parameters, .rowname="A")- 4.039903) < 1e-6);
        assert(fabs(apop_data_get(est->parameters, .rowname="B")- 1.494694) < 1e-6);
    #endif
}
