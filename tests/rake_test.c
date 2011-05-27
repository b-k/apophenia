#include <apop.h>

/* Some OK tests on the raking procedure. We assert that a regression on the raw data 
and a set of dummies is equivalent to a regression on the base data. */

double weights_are_one(apop_data *in){assert(gsl_vector_get(in->weights, 0)==1); return 0;}
double equal_or_absent(apop_data *in){
    assert(apop_data_get(in, .colname="a") != 3);
    assert(apop_data_get(in, .colname="b") != 7);
    assert(gsl_vector_get(in->weights, 0)==1);
    return 0;
}

double compare_results(apop_data *in, void *other, int index){
    double other_val= apop_data_get(((apop_data *)other), .row=index, .col=-1);
    assert(fabs(gsl_vector_get(in->vector, 0)-other_val)< 1e-3); 
    return 0;
}
    
void test_raking(){
    //trivial case: if all margins are equal, MLE is to give equal weights.
    int a, b, c;
    apop_query("create table equals (a,b,c)");
    for (a=0; a < 10; a++)
        for (b=0; b < 10; b++)
            for (c=0; c < 10; c++)
                apop_query("insert into equals values (%i, %i, %i)", a,b,c);
    apop_data *equal_weights = apop_rake("equals");
    apop_map(equal_weights, .fn_r=weights_are_one);

    //structural zeros should be missing from the output table.
    apop_data *with_zeros = apop_rake("equals", .structural_zeros="a+0.0==3 or b+0.0==7");
    //apop_data_show(with_zeros);
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
    apop_data *inequal_weights = apop_rake("inequals", .all_vars="a|b|c",.contrasts=contrasts, .weights_col="weights");


    //regress using the estimates from the raking
    apop_data_rm_columns(inequal_weights, (int[]){0, 0, 1});
    apop_data_pmf_compress(inequal_weights);
    apop_data_to_dummies(inequal_weights, .type='d', .append='y', .remove='y'); //first A column
    apop_data_to_dummies(inequal_weights, .col=apop_name_find(inequal_weights->names, "b", 'c'), 
                                          .type='d', .append='y', .remove='y');
    inequal_weights->vector = inequal_weights->weights;
    inequal_weights->weights = NULL;
    //apop_data_show(inequal_weights);
    apop_model *raked_ols = apop_estimate(inequal_weights, apop_ols);  

   
    //regress using the original data
    apop_data *t = apop_query_to_mixed_data("vmm", "select sum(weights), a, b from inequals group by a, b");//force linear, not affine.
    apop_data_to_dummies(t, .type='d', .append='y', .remove='y');
    apop_data_to_dummies(t, .col=apop_name_find(t->names, "b", 'c'), .type='d', .append='y', .remove='y');
//    apop_data_to_dummies(t, .col=1, .type='d', .append='y', .remove='y'); //affine version.
    apop_model *ols_out = apop_estimate(t, apop_ols);  
    //apop_data_show(t);
//    apop_data_show(ols_out->parameters);

    apop_map(ols_out->parameters, .fn_rpi=compare_results, .param= raked_ols->parameters);

}
