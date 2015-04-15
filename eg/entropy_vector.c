#include <apop.h>

#define Diff(left, right, eps) Apop_stopif(fabs((left)-(right))>(eps), abort(), 0, "%g is too different from %g (abitrary limit=%g).", (double)(left), (double)(right), eps)

long double entropy_base_2(gsl_vector *x) {
    return apop_vector_entropy(x)/log(2);
}

apop_data *flip_a_coin(int how_many){
    return apop_model_draws(apop_model_set_parameters(apop_bernoulli, .5), how_many);
}

int main(){
    //zero data => entropy zero
    gsl_vector *v = gsl_vector_calloc(1);
    assert(apop_vector_entropy(v) == 0);

    //negative data => NaN
    gsl_vector_set(v, 0, -1);
    int v1 = apop_opts.verbose;
    apop_opts.verbose = -1;
    assert(isnan(apop_vector_entropy(v)));
    apop_opts.verbose = v1;

    //N equiprobable bins => entropy = log(N)
    v = apop_vector_realloc(v, 100);
    gsl_vector_set_all(v, 1./100);
    Diff(log(100), apop_vector_entropy(v), 1e-5);

    //flip two coins.
    apop_data *coin_flips = flip_a_coin(10000);
    apop_data *c2         = flip_a_coin(10000);
    apop_data_stack(c2, coin_flips, 'c', .inplace='y');

    //entropy of one coin flip in base2 == 1
    apop_data_pmf_compress(coin_flips);
    Diff(entropy_base_2(coin_flips->weights), 1, 1e-3);

    //entropy of two coin flips in base2 == 2
    apop_data_pmf_compress(c2);
    Diff(entropy_base_2(c2->weights), 2, 1e-3);

    apop_data_free(coin_flips);
    apop_data_free(c2);
    gsl_vector_free(v);
}

/*
export LDLIBS="`pkg-config --libs apophenia`"
export CFLAGS="-g -Wall `pkg-config --cflags apophenia` -O3 -std=gnu11"
*/
