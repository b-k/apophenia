#include <apop.h>

int main(){
    apop_opts.db_engine='s'; //SQLite only.
    apop_query("create table atab (a numeric)");
    for (int i=0; i< 1e5; i++)
        apop_query("insert into atab values(ran())");
    apop_query("create table powa as "
            "select a, pow(a, 2) as sq, pow(a, 0.5) as sqrt "
            "from atab");

    //compare the std dev of a uniform as reported by the 
    //database routine, the matrix routine, and math.
    double db_pop_stddev = apop_query_to_float("select stddev_pop(a) from powa");
    apop_data *d = apop_query_to_data("select * from powa");
    apop_data *cov = apop_data_covariance(d);
    double matrix_pop_stddev = sqrt(apop_data_get(cov)*(d->matrix->size1/(d->matrix->size1-1.)));
    assert(fabs(db_pop_stddev - matrix_pop_stddev) < 1e-4);
    double actual_stddev = sqrt(2*gsl_pow_3(.5)/3);
    assert(fabs(db_pop_stddev - actual_stddev) < 1e-3);

    float sq_mean = apop_query_to_float("select avg(sq) from powa");
    float actual_sq_mean = 1./3;
    assert(fabs(sq_mean - actual_sq_mean) < 1e-3);

    float sqrt_mean = apop_query_to_float("select avg(sqrt) from powa");
    float actual_sqrt_mean = 2./3;
    assert(fabs(sqrt_mean - actual_sqrt_mean) < 1e-3);
}
