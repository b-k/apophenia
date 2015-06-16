#include <apop.h>

#define Diff(L, R) assert(fabs((L)-(R)<1e-4));
#define Diff2(L, R) assert(fabs((L)-(R)<1e-3));
#define getrow(rowname) apop_data_get(row, .colname=#rowname)

double test_all(apop_data *row){
    Diff(gsl_pow_2(getrow(root)), getrow(rr))
    Diff(getrow(ln), getrow(L10)*log(10))
    Diff(getrow(rr), getrow(rragain))
    Diff(getrow(one), 1)
    return 0;
}

int main(){
    apop_opts.db_engine='s'; //SQLite only.

    //create a table with two rows.
    //We didn't explicitly open a db with apop_db_open,
    //so this will be an in-memory SQLite db.
    apop_query("create table a(b); "
               "insert into a values(1); "
               "insert into a values(1); "

                "create table randoms as "
                "select ran() as rr "
                /* join to create 2^13=8192 rows*/
                "from a,a,a,a,a,a,a,a,a,a,a,a,a;");
    apop_data *d = apop_query_to_data(
            "select rr, sqrt(rr) as root, "
            "log(rr) as ln, log10(rr) as L10, "
            "exp(log(rr)) as rragain, "
            "pow(sin(rr),2)+pow(cos(rr),2) as one "
            "from randoms");
    apop_map(d, .fn_r=test_all);

    //the pop variance of a Uniform[0,1]=1/12; kurtosis=1/80.
    Apop_col_tv(d, "rr", rrow);
    Diff(apop_var(rrow)*8191./8192., 1/12. );
    Diff(apop_vector_kurtosis(rrow)*8191./8192., 1/80.);//approx.

    Diff(apop_query_to_float("select stddev(rr) from randoms"), 
                sqrt(1/12.)*8192./8191);


    //compare the std dev of a uniform as reported by the 
    //database routine, the matrix routine, and math.
    apop_query("create table atab (a numeric)");
    for (int i=0; i< 2e5; i++)
        apop_query("insert into atab values(ran())");
    apop_query("create table powa as "
            "select a, pow(a, 2) as sq, pow(a, 0.5) as sqrt "
            "from atab");

    double db_pop_stddev = apop_query_to_float("select stddev_pop(a) from powa");
    d = apop_query_to_data("select * from powa");
    //get the full covariance matrix, but just use the (0,0)th elmt.
    apop_data *cov = apop_data_covariance(d);
    double matrix_pop_stddev = sqrt(apop_data_get(cov)*(d->matrix->size1/(d->matrix->size1-1.)));
    Diff(db_pop_stddev, matrix_pop_stddev);
    double actual_stddev = sqrt(2*gsl_pow_3(.5)/3);
    Diff2(db_pop_stddev, actual_stddev);

    float sq_mean = apop_query_to_float("select avg(sq) from powa");
    float actual_sq_mean = 1./3;
    Diff2(sq_mean, actual_sq_mean);

    float sqrt_mean = apop_query_to_float("select avg(sqrt) from powa");
    float actual_sqrt_mean = 2./3;
    Diff2(sqrt_mean, actual_sqrt_mean);
}
