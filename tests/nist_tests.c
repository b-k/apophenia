/* These are stats tests from NIST. See http://www.itl.nist.gov/div898/strd/
Notice that I use various levels of tolerance, so this gives you an idea
of the relative accuracies of various operations. */

#include <apop.h>
#include <unistd.h>

#define TOL 1e-15
#define TOL2 1e-5
#define TOL3 1e-9

#define Diff(L, R, eps) Apop_stopif(isnan(L-R) || fabs((L)-(R))>(eps), abort(), 0, "%g is too different from %g (abitrary limit=%g).", (double)(L), (double)(R), eps);

void pontius(){
    apop_text_to_db("pontius.dat","pont", .delimiters=" ");
    apop_data *d = apop_query_to_data("select y, x, pow(x,2) as p from pont");
    apop_model *est =  apop_estimate(d, apop_ols);

    Diff(apop_data_get(est->parameters, 0, -1), 0.673565789473684E-03, TOL3);
    Diff(apop_data_get(est->parameters, 1, -1), 0.732059160401003E-06, TOL);
    Diff(apop_data_get(est->parameters, 2, -1), -0.316081871345029E-14, TOL);
    apop_data *cov = apop_data_get_page(est->parameters, "<covariance>");
    Diff(apop_data_get(cov, 0, 0), pow(0.107938612033077E-03,2), TOL2);
    Diff(apop_data_get(cov, 1, 1), pow(0.157817399981659E-09,2), TOL2);
    Diff(apop_data_get(cov, 2, 2), pow(0.486652849992036E-16,2), TOL2);
    Diff(apop_data_get(est->info, .rowname="R squared"), 0.999999900178537, TOL);
    Diff(apop_data_get(est->info, .rowname="SSR"), 15.6040343244198, TOL3);
}

void wampler1(){
    apop_text_to_db("wampler1.dat","w1", .delimiters=" ");
    apop_data *d = apop_query_to_data("select y, x, pow(x,2) as p2, \
                                pow(x,3) as p3, pow(x,4) as p4, pow(x,5) as p5 from w1");
    apop_model *est = apop_estimate(d, apop_ols);
    for (int i=0; i<6; i++)
        Diff(apop_data_get(est->parameters, i, -1) ,1 , 1e-3);
    apop_data *cov = apop_data_get_page(est->parameters, "<covariance>");
    for (int i=0; i<6; i++)
        Diff(apop_data_get(cov, i, i), 0, TOL2);
    Diff(apop_data_get(est->info, .rowname="R squared"), 1, TOL);
}

void numacc4(){
    apop_data *d  = apop_text_to_data("numacc4.dat");
    Apop_col_v(d, 0, v)
    Diff(apop_vector_mean(v), 10000000.2, 1e-5);
    Diff(apop_vector_var(v)*(v->size -1)/v->size, 0.01, TOL3);
    //I don't do this yet:
    //Sample Autocorrelation Coefficient (lag 1) r(1):   -0.999     (exact)
}

int main(){
    chdir(Datadir); //Datadir is defined via autoconf.

    pontius();
    wampler1();
    numacc4();
}
