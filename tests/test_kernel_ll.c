/* The kernel log likelihood is somewhat convoluted, to retain numeric precision. This
test compares the result to the direct approach, in the case of a few data sets
that are small enough that we don't have to worry about underflow.
*/
#include <apop.h>

long double sum_of_parts(apop_data *d1, apop_data *target){
    apop_data *row, *trow;
    long double p=1;
    for(int i=0; (trow=Apop_r(target, i)) && trow->vector; i++){
        long double tprob=0;
        for(int j=0; (row=Apop_r(d1, j)) && row->vector; j++)
            tprob+= apop_p(trow, apop_model_set_parameters(apop_normal, *row->vector->data, 1));
        p *=tprob/d1->vector->size;
    }
    return p;
}

void go(apop_data *d1, apop_data *d2){
    apop_model *k = apop_model_copy_set(apop_kernel_density, apop_kernel_density, .base_data=d1);
    assert(fabs(apop_p(d2, k)- sum_of_parts(d1, d2)) < 1e-5);

    apop_model *test_copying = apop_model_copy(k);
    assert(fabs(apop_p(d2, test_copying)- sum_of_parts(d1, d2)) < 1e-8);
    apop_model_free(k);
    apop_model_free(test_copying);
}

int main(){
    apop_data *d1= apop_data_falloc((4), 2,4,6,8);
    go(d1, apop_data_falloc((4), 1,3,5,7));
    go(d1, apop_data_falloc((4), 1,1,1,1));

    apop_data *d2= apop_data_falloc((4, 4, 1), 2,1.1, 4,2.2, 6,3.1, 8,0);
    go(d2, apop_data_falloc((4, 4, 1), 1, 0, 3,0, 5, 0, 7, 0));
}
