#include <apop.h>

#define Diff(L, R) assert(fabs((L)-(R)<1e-4));
#define getrow(rowname) apop_data_get(row, .colname=#rowname)

double test_all(apop_data *row){
    Diff(gsl_pow_2(getrow(root)), getrow(rr))
    Diff(getrow(ln), getrow(l10)*log(10))
    Diff(getrow(rr), getrow(rragain))
    Diff(getrow(one), 1)
    return 0;
}

int main(){
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
            "log(rr) as ln, log10(rr) as l10, "
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
}
