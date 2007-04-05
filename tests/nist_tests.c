#include <apop.h>

#define TOL 1e-15
#define TOL2 1e-5
#define TOL3 1e-9
#define TOL4 1e-4

void pontius(){
    apop_text_to_db("pontius.dat","d", 0,1, NULL);
apop_data *d        = apop_query_to_data("select y, x, pow(x,2) as p from d");
apop_params *est  =  apop_OLS.estimate(d, NULL);

    assert(fabs(apop_data_get(est->parameters, 0, -1) - 0.673565789473684E-03) < TOL);
    assert(fabs(apop_data_get(est->parameters, 1, -1) - 0.732059160401003E-06) < TOL);
    assert(fabs(apop_data_get(est->parameters, 2, -1) - -0.316081871345029E-14)    < TOL);
    assert(fabs(apop_data_get(est->covariance, 0, 0) - pow(0.107938612033077E-03,2))    < TOL2);
    assert(fabs(apop_data_get(est->covariance, 1, 1) - pow(0.157817399981659E-09,2))    < TOL2);
    assert(fabs(apop_data_get(est->covariance, 2, 2) - pow(0.486652849992036E-16,2))    < TOL2);
apop_data *cc   = apop_estimate_correlation_coefficient(est);
    assert(fabs(apop_data_get_tn(cc, "R.sq.*", -1) - 0.999999900178537)    < TOL);
    assert(fabs(apop_data_get_tn(cc, "SSR", -1) - 15.6040343244198)    < TOL3);
}

void wampler1(){
    apop_text_to_db("wampler1.dat","w1", 0,1, NULL);
int             i;
apop_data       *d    = apop_query_to_data("select y, x, pow(x,2) as p2, \
                                pow(x,3) as p3, pow(x,4) as p4, pow(x,5) as p5 from w1");
apop_params   *est  =  apop_OLS.estimate(d, NULL);
    for (i=0; i<6; i++)
        assert(fabs(apop_data_get(est->parameters, i, -1) - 1) < TOL4);
    for (i=0; i<6; i++)
        assert(fabs(apop_data_get(est->covariance, i, i)) < TOL2);
apop_data *cc   = apop_estimate_correlation_coefficient(est);
    assert(fabs(apop_data_get_tn(cc, "R.sq.*", -1) - 1)    < TOL);
}

void numacc4(){
apop_data   *d  = apop_text_to_data("numacc4.dat", 0, 0);
gsl_vector  v   = gsl_matrix_column(d->matrix, 0).vector;
    assert(apop_vector_mean(&v) == 10000000.2);
    assert(apop_vector_var(&v) - 0.01 < TOL3);
    //I don't do this yet:
    //Sample Autocorrelation Coefficient (lag 1) r(1):   -0.999     (exact)
}

int nist_tests(){
    pontius();
    wampler1();
    numacc4();
    return 0;
}
