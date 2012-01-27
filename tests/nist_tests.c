/* These are stats tests from NIST. See http://www.itl.nist.gov/div898/strd/
Notice that I use various levels of tolerance, so this gives you an idea
of the relative accuracies of various operations. */


#include <apop.h>

#define TOL 1e-15
#define TOL2 1e-5
#define TOL3 1e-9
#define TOL4 1e-4

void pontius(){
    apop_text_to_db("pontius.dat","d", .delimiters=" ");
apop_data *d        = apop_query_to_data("select y, x, pow(x,2) as p from d");
apop_model *est  =  apop_estimate(d, apop_ols);

    assert(fabs(apop_data_get(est->parameters, 0, -1) - 0.673565789473684E-03) < TOL3);
    assert(fabs(apop_data_get(est->parameters, 1, -1) - 0.732059160401003E-06) < TOL);
    assert(fabs(apop_data_get(est->parameters, 2, -1) - -0.316081871345029E-14)    < TOL);
    apop_data *cov = apop_data_get_page(est->parameters, "cov");
    assert(fabs(apop_data_get(cov, 0, 0) - pow(0.107938612033077E-03,2))    < TOL2);
    assert(fabs(apop_data_get(cov, 1, 1) - pow(0.157817399981659E-09,2))    < TOL2);
    assert(fabs(apop_data_get(cov, 2, 2) - pow(0.486652849992036E-16,2))    < TOL2);
    assert(fabs(apop_data_get(est->info, .rowname="R.sq.*") - 0.999999900178537)    < TOL);
    assert(fabs(apop_data_get(est->info, .rowname="SSR") - 15.6040343244198)    < TOL3);
}

void wampler1(){
    apop_text_to_db("wampler1.dat","w1", .delimiters=" ");
int             i;
apop_data       *d    = apop_query_to_data("select y, x, pow(x,2) as p2, \
                                pow(x,3) as p3, pow(x,4) as p4, pow(x,5) as p5 from w1");
apop_model   *est  =  apop_estimate(d, apop_ols);
    for (i=0; i<6; i++)
        assert(fabs(apop_data_get(est->parameters, i, -1) - 1) < TOL4*10);
    apop_data *cov = apop_data_get_page(est->parameters, "cov");
    for (i=0; i<6; i++)
        assert(fabs(apop_data_get(cov, i, i)) < TOL2);
    assert(fabs(apop_data_get(est->info, .rowname="R.sq.*") - 1)    < TOL);
}

void numacc4(){
apop_data   *d  = apop_text_to_data("numacc4.dat", 0, 0);
    APOP_COL(d, 0, v)
    assert(apop_vector_mean(v) == 10000000.2);
    assert(apop_vector_var(v)*(v->size -1)/v->size - 0.01 < TOL3);
    //I don't do this yet:
    //Sample Autocorrelation Coefficient (lag 1) r(1):   -0.999     (exact)
}

void nist_tests(){
    pontius();
    wampler1();
    numacc4();
}
