#include <apophenia/headers.h>

double  pu, pd, pl, pr;
int     n;

void find_marginals(apop_data *d){
  gsl_vector  v;
    v   = gsl_matrix_row(d->matrix, 0).vector;
    pu  = apop_vector_sum(&v)/n;
    v   = gsl_matrix_column(d->matrix, 0).vector;
    pl  = apop_vector_sum(&v)/n;
    pd  = 1 - pu;
    pr  = 1 - pl;
}

double one_chi_sq(double o, double e){
    return gsl_pow_2(o - e)/e;
}

double calc_chi_squared(apop_data *d){
  double chi_ul = one_chi_sq(apop_data_get(d,0,0), n*pu*pl);
  double chi_ur = one_chi_sq(apop_data_get(d,0,1), n*pu*pr);
  double chi_dl = one_chi_sq(apop_data_get(d,1,0), n*pd*pl);
  double chi_dr = one_chi_sq(apop_data_get(d,1,1), n*pd*pr);
    return chi_ul + chi_ur + chi_dl + chi_dr;
}

int main(){
  double      data[]      = { 30,86,
                    24,38 };
  apop_data   *testdata   = apop_line_to_data(data,2,2);
  double      stat, chisq;
    n   = data[0]+data[1]+data[2]+data[3];
    find_marginals(testdata);
    stat   = calc_chi_squared(testdata);
    chisq   = gsl_cdf_chisq_Q(stat, 1);
    printf("chi squared statistic: %g; p, Chi-squared: %g\n", stat,chisq);
    apop_data_show(apop_test_fisher_exact(testdata));
    return 0;
}
