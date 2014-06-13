#include <apop.h>
/* This sample code sets the elements of a data set's vector to one if the index is even.
   Then, via the weights vector, it adds up the even indices.

   There is really no need to use the weights vector; this code snippet is an
   element of Apophenia's test suite, and goes the long way to test that the weights are
   correctly handled. */

double set_vector_to_even(apop_data * r, int index){
    apop_data_set(r, 0, -1, 1 - (index %2));
    return 0;
}

double set_weight_to_index(apop_data * r, int index){ 
    gsl_vector_set(r->weights, 0, index); 
    return 0;
}

double weight_given_even(apop_data *r){ 
    return gsl_vector_get(r->vector, 0) ? gsl_vector_get(r->weights, 0) : 0; 
}

int main(){
    apop_data *d = apop_data_alloc(100);
    d->weights = gsl_vector_alloc(100);
    apop_map(d, .fn_ri=set_vector_to_even, .inplace='v');
    apop_map(d, .fn_ri=set_weight_to_index, .inplace='v');
    double sum = apop_map_sum(d, .fn_r = weight_given_even);
    assert(sum == 49*25*2);
}
