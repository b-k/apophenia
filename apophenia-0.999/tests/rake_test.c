#include <apop.h>

#define Diff(L, R, eps) Apop_assert_n(fabs((L)-(R))<(eps), "%g is too different from %g (abitrary limit=%g).", (double)(L), (double)(R), eps);

void rake_check(apop_model *base, apop_model *fitted){
    Diff(apop_query_to_float("select sum(weights) from raked where first=1"), 12, 1e-4);
    Diff(apop_query_to_float("select sum(weights) from raked where first=2"), 20, 1e-4);
    Diff(apop_query_to_float("select sum(weights) from raked where second=1"), 25, 1e-4);
    Diff(apop_query_to_float("select sum(weights) from raked where second=2"), 7, 1e-4);
    /* Raking minimizes KL divergence given the margin constraints. So nudging the table 
       in a manner that fits the constraints should raise KLdiv. */
    double kl1= apop_kl_divergence(base, fitted);
    *gsl_vector_ptr(fitted->data->weights, 0) += 0.05;
    *gsl_vector_ptr(fitted->data->weights, 1) += -0.05;
    *gsl_vector_ptr(fitted->data->weights, 2) += -0.05;
    *gsl_vector_ptr(fitted->data->weights, 3) += 0.05;
    assert(kl1 < apop_kl_divergence(base, fitted));
}

//these work by checking that K-L divergence shrunk, and that individual margins are correct.
void test_raking_further(){
    apop_table_exists("rake_test", 'd');
    apop_query("create table rake_test (first, second, weights);"
            "insert into rake_test values(1, 1, 10);"
            "insert into rake_test values(1, 2, 2);"
            "insert into rake_test values(2, 1, 15);"
            "insert into rake_test values(2, 2, 5);"
            );

    //Synthetic data, starting at all ones.
    apop_data_print(
            apop_rake(.margin_table="rake_test", .count_col="weights", 
                .contrasts=(char*[]){"first", "second"}, .contrast_ct=2),
        .output_name="raked", .output_type='d');
    apop_model *fitted= apop_estimate(apop_query_to_mixed_data("mmw", "select * from raked"), apop_pmf);
    rake_check(apop_estimate(apop_query_to_mixed_data("mmw", "select first, second, 1 from rake_test"), apop_pmf), fitted);

        //With an alternate init table
    apop_table_exists("raked", 'd');
    apop_query("create table rakeinit (first, second, weights);"
            "insert into rakeinit values(1, 1, 32);"
            "insert into rakeinit values(1, 2, 289);"
            "insert into rakeinit values(2, 1, 19);"
            "insert into rakeinit values(2, 2, 5447);"
            );
    apop_data_print(
            apop_rake(.margin_table="rake_test", .count_col="weights", 
                .contrasts=(char*[]){"first", "second"}, .contrast_ct=2, .init_table="rakeinit", .init_count_col="weights"),
        .output_name="raked", .output_type='d');
    //apop_data_show(apop_query_to_data("select * from raked"));

    apop_model *base= apop_estimate(apop_query_to_mixed_data("mmw", "select * from rakeinit"), apop_pmf);
    fitted= apop_estimate(apop_query_to_mixed_data("mmw", "select * from raked"), apop_pmf);
    rake_check(base, fitted);
}


/* Some OK tests on the raking procedure. We assert that a regression on the raw data 
and a set of dummies is equivalent to a regression on the base data. */

double weights_are_one(apop_data *in){assert(gsl_vector_get(in->weights, 0)==1); return 0;}
double equal_or_absent(apop_data *in){
    assert(apop_data_get(in, .colname="a") != 3);
    assert(apop_data_get(in, .colname="b") != 7);
    assert(gsl_vector_get(in->weights, 0)==100./81);
    return 0;
}

double compare_results(apop_data *in, void *other, int index){
    double other_val= apop_data_get(((apop_data *)other), .row=index, .col=-1);
    assert(fabs(gsl_vector_get(in->vector, 0)-other_val)< 1e-3); 
    return 0;
}
    
int main(){
    //trivial case: if all margins are equal, MLE is to give equal weights.
    int a, b, c;
    apop_query("create table equals (a,b,c)");
    for (a=0; a < 10; a++)
        for (b=0; b < 10; b++)
            for (c=0; c < 10; c++)
                apop_query("insert into equals values (%i, %i, %i)", a,b,c);
    apop_data *equal_weights = apop_rake(.margin_table="equals", .init_table="equals");
    apop_map(equal_weights, .fn_r=weights_are_one);

    //structural zeros should be missing from the output table.
    apop_data *with_zeros = apop_rake("equals", .structural_zeros="a+0.0==3 or b+0.0==7");
    apop_map(with_zeros, .fn_r=equal_or_absent);

    //Not-trivial case.
    //Regression parameters on the inequal weights should match regression parameters on the raw data.
    apop_query("create table inequals (a,b,c, weights)");
    gsl_rng *r = apop_rng_alloc(25);
    for (a=0; a < 3; a++)
        for (b=0; b < 4; b++)
            for (c=0; c < 10; c++)
                apop_query("insert into inequals values (%i, %i, %i, %g)", a,b,c, a/2.+gsl_rng_uniform(r));
    char *contrasts[] ={"a|b"};
    apop_data *inequal_weights = apop_rake("inequals", .var_list=(char*[]){"a", "b", "c"}, .var_ct=3, .contrasts=contrasts, .count_col="weights", .contrast_ct=1, .tolerance=1e-8);


    //regress using the estimates from the raking
    apop_data_rm_columns(inequal_weights, (int[]){0, 0, 1});
    apop_data_pmf_compress(inequal_weights);
    apop_data_to_dummies(inequal_weights, .type='d', .append='y', .remove='y'); //first A column
    apop_data_to_dummies(inequal_weights, .col=apop_name_find(inequal_weights->names, "b", 'c'), 
                                          .type='d', .append='y', .remove='y');
    inequal_weights->vector = inequal_weights->weights;
    inequal_weights->weights = NULL;
    apop_model *raked_ols = apop_estimate(inequal_weights, apop_ols);  

   
    //regress using the original data
    apop_data *t = apop_query_to_mixed_data("vmm", "select sum(weights), a, b from inequals group by a, b");//force linear, not affine.
    apop_data_to_dummies(t, .type='d', .append='y', .remove='y');
    apop_data_to_dummies(t, .col=apop_name_find(t->names, "b", 'c'), .type='d', .append='y', .remove='y');
    apop_model *ols_out = apop_estimate(t, apop_ols);

    apop_map(ols_out->parameters, .fn_rpi=compare_results, .param= raked_ols->parameters);

    test_raking_further();
}
