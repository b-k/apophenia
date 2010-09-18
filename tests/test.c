#include <apop.h>
#include "nist_tests.c"

#define Diff(L, R, eps) Apop_assert(fabs((L)-(R)<(eps)), "%g is too different from %g (abitrary limit=%g).", (double)(L), (double)(R), eps);

//I'm using the test script an experiment to see if 
//these macros add any value.
#define APOP_VECTOR_ALLOC(name, r) gsl_vector *name = gsl_vector_alloc(r)
#define APOP_DATA_ALLOC(name, r, c) apop_data *name = apop_data_alloc(0, (r),(c))
#define APOP_RNG_ALLOC(name, seed) gsl_rng *name = apop_rng_alloc(seed)

//These test functions are also displayed in the documentation as examples.
#include "../eg/test_kl_divergence.c"
#include "../eg/test_strip_dots.c"
#include "../eg/test_distances.c"
#include "../eg/test_harmonic.c"
#include "../eg/test_updating.c" 
#include "../eg/test_pruning.c"     // test_prune_cols()
#include "../eg/apop_map_row.c" 
#include "../eg/test_strcmp.c" 
#include "../eg/test_fisher.c" 
#include "../eg/test_ranks.c" 
#include "../eg/test_regex.c" 
#include "../eg/pmf_test.c" 


/*
Some of these tests are mechanical tests that data gets shunted to the right place and
that nothing segfaults. Those are not very conceptually difficult.

Many of these tests are much more computationally-intensive than the norm, both in
terms of compute time, and in terms of the expectations of the algorithm. For example,
let us say that we wish to verify the results of a regression. Some systems have a
canned screenshot of the 'correct' regression results that ships with the test suite,
and compare a screenshot of the run to the canned version. I don't get much confidence
from this---what if the canned screenshot is wrong? Better would be to know something
about the regression results (like the relation between the common F statistic, SSR,
and SSE) and check that the fact always holds.

Some other examples: given a set of parameters for a distribution, make a million
draws from the distribution given those parameters, then estimate the parameters of
the distribution; the "before" and "after" parameters should match. Or, if there are
multiple methods of doing Bayesian updating, the output distributions should match.

Those claims are true as N goes to infinity; for finite N the routines have to strike
a balance. How many draws should I make, and how much user time should I waste, before
measuring the error, and what error tolerance should I set? This is a difficult balance,
and is to some extent the key problem behind all of numeric computing.

There are two types of error bounds here. One is tighter, and therefore more prone
to false alarms, but really forces us to write better numeric code. The other is
much more permissive, and just tells us whether the computation failed to go in the
right direction.  Users who run 'make test' will be running the second type of test,
because I (BK) just got sick of people sending me bug reports that a test failed
because it reported an error of 1e-5 when it should have been 1e-8. There is always
room for better numeric precision; we all know this with or without reminders from
the post-install tests.
*/

//One-liners for mapply:
gsl_rng *r_global;
//void random_draw(double *in) { *in = gsl_rng_uniform(r_global);}
double nan_map(double in){return gsl_isnan(in);}

#ifdef FULL_TOLERANCE
double tol6 = 1e-6;
double tol5 = 1e-5;
double tol3 = 1e-3;
double tol2 = 1e-2;
double tol1 = 1e-1;
#else
double tol6 = 1e-1;
double tol5 = 1e-1;
double tol3 = 1e-1;
double tol2 = 1e-1;
double tol1 = 1e-1;
#endif

int     len                 = 8000;
int     verbose             = 1;

void test_nan_data();
void db_to_text();

static void test_printing(){
    //This compares printed output to the printed output in the attached file. 
    char outfile[] = "print_test.out";

    if (!apop_table_exists("nandata"))
        test_nan_data();
    apop_opts.output_type ='s';
    gsl_matrix *m  = apop_query_to_matrix("select * from nandata");
    apop_matrix_print(m, .output_file=outfile, .output_append='w');

apop_system("cp %s xxx", outfile);

    if (!apop_table_exists("d"))
        db_to_text();
    apop_data *d = apop_query_to_mixed_data ("tvttmmwt", "select * from d");
    FILE *f = fopen(outfile, "a");
    fprintf(f, "\nand a full vector+matrix+text+weights data set, formatted for computer reading:\n");
    strcpy(apop_opts.output_delimiter, "\t| ");

    apop_data_print(d, .output_pipe =f);
    fprintf(f, "\nand just the weights vector:\n");
    strcpy(apop_opts.output_delimiter, "\t");
    apop_opts.output_type = 'p';
    apop_opts.output_pipe = f;
    apop_vector_print(d->weights);
    apop_opts.output_type = 's';
    fclose(f);
    int has_diffs = apop_system("diff -b printing_sample %s", outfile);
    assert(!has_diffs);
    //apop_system("rm %s", outfile);
}

void v_pow10(double *in){ *in = pow(10,*in);}

static void log_and_exp(gsl_rng *r){
    int i;
    apop_data *d = apop_data_alloc(100,2);
    apop_name_add(d->names, "10", 'c');
    apop_name_add(d->names, "e", 'c');
    for (i=0; i< 100; i++){
        apop_data_set(d, i, .colname="10", .val=gsl_rng_uniform(r)*10);
        apop_data_set(d, i, .colname="e", .val=gsl_rng_uniform(r)*10);
    }
    apop_data *d2 = apop_data_copy(d);
    Apop_col_t(d, "10", tencol);
    apop_vector_log10(tencol);
    apop_vector_apply(tencol, v_pow10);
    Apop_col(d2, 0, o_tencol);
    assert(apop_vector_distance(tencol, o_tencol) < 1e-3);
    
    Apop_col_t(d, "e", ecol);
    apop_vector_log(ecol);
    apop_vector_exp(ecol);
    Apop_col(d2, 1, o_ecol);
    assert(apop_vector_distance(ecol, o_ecol) < 1e-3);

}

static void compare_mvn_estimates(apop_model *L, apop_model *R, double tolerance){
    gsl_vector_sub(L->parameters->vector, R->parameters->vector);
    gsl_matrix_sub(L->parameters->matrix, R->parameters->matrix);
    assert(fabs(apop_sum(L->parameters->vector)) + fabs (apop_matrix_sum(L->parameters->matrix)) < tolerance);
}

void test_ml_imputation(gsl_rng *r){
    size_t len = 4e4;
    int i,j;
    apop_data *fillme = apop_data_alloc(len, 3);
    apop_model *mvn = apop_model_copy(apop_multivariate_normal);
    mvn->parameters = apop_data_alloc(3, 3, 3);
    for(i=0; i < 3; i ++)
        for(j=-1; j < 3; j ++)
            apop_data_set(mvn->parameters, i, j, gsl_rng_uniform(r));
    //now make your random garbage symmetric
    for(i=0; i < 3; i ++)
        for(j=i+1; j < 3; j ++)
            apop_data_set(mvn->parameters, j, i, apop_data_get(mvn->parameters, i, j));
    apop_matrix_to_positive_semidefinite(mvn->parameters->matrix);
    for(i=0; i < len; i ++){
        Apop_row(fillme, i, row);  //relying on the rows having stride=1.
        apop_draw(row->data, r, mvn);
    }
    //apop_data_show(mvn->parameters);
    apop_model *est = apop_estimate(fillme, apop_multivariate_normal);
    //apop_data_show(est->parameters);
    compare_mvn_estimates(est, mvn, 1e-1);

    double pct_to_delete = 0.01;
    int max_to_delete = 7, ctr = 0;
    for(i=0; i < len && ctr < max_to_delete; i ++)
        for(j=0; j < 3; j ++)
            if (gsl_rng_uniform(r) < pct_to_delete){
                apop_data_set(fillme, i, j, GSL_NAN);
                ctr++;
            }
    apop_ml_imputation(fillme, mvn); 
    apop_model *est2 = apop_estimate(fillme, apop_multivariate_normal);
    //apop_data_show(est2->parameters);
    compare_mvn_estimates(est2, mvn, 1e-1);
}

void test_percentiles(){
  gsl_vector    *v      = gsl_vector_alloc(307);
  int           i;
    for (i=0; i< 307; i++)
        gsl_vector_set(v, i, i);
  double        *pcts_up    = apop_vector_percentiles(v, 'u');
  double        *pcts_down  = apop_vector_percentiles(v, 'd');
  double        *pcts_avg   = apop_vector_percentiles(v, 'a');
    for (i=0; i< 101; i++){
        assert(pcts_up[i] >= pcts_down[i]);
        assert(pcts_up[i] >= pcts_avg[i]);
        assert(pcts_avg[i] >= pcts_down[i]);
    }
    assert(pcts_up[100] == pcts_down[100] && pcts_avg[100] == pcts_down[100]);
    assert(pcts_up[0] == pcts_down[0] && pcts_avg[0] == pcts_down[0]);
    assert(pcts_avg[50] == (pcts_down[50] + pcts_up[50])/2);
}

void test_score(){
  int         i, j, len = 1e5;
  gsl_rng     *r        = apop_rng_alloc(123);
  apop_data   *data     = apop_data_alloc(len,1);
    for (i=0; i<10; i++){
        apop_model  *source   = apop_model_set_parameters(apop_normal, 
                                    gsl_ran_flat(r, -5, 5), gsl_ran_flat(r, .01, 5));
        for (j=0; j< len; j++)
            apop_draw(gsl_matrix_ptr(data->matrix, j, 0), r, source);
        apop_model *estme = apop_model_copy(apop_normal);
        Apop_model_add_group(estme, apop_mle, .method= APOP_SIMAN,.parent= estme);
        apop_model *out = apop_maximum_likelihood(data, estme);

        apop_model *straight_est = apop_estimate(data, apop_normal);
        Diff (straight_est->parameters->vector->data[0], source->parameters->vector->data[0], tol1);
        Diff (straight_est->parameters->vector->data[1], source->parameters->vector->data[1], tol1);

        double sigsqn = gsl_pow_2(out->parameters->vector->data[1])/len;
        apop_data *cov = apop_data_get_page(out->parameters, "cov");
        Diff (apop_data_get(cov, 0,0),sigsqn , tol3);
        Diff (apop_data_get(cov, 1,1),sigsqn/2 , tol3);
        assert(apop_data_get(cov, 0,1) + apop_data_get(cov, 0,1) < tol3);
        apop_model_free(out);
        printf(".");
        apop_model_free(source); 
    }
    apop_data_free(data);
}

//This tests the database-side functions.
void test_skew_and_kurt(){
  gsl_rng *r  = apop_rng_alloc(time(0));
  int     i;
    apop_table_exists(.remove=1, .name="t");
    apop_query("create table t(vals)");
    for(i=0;i<1e4; i++){
        apop_query("insert into t values(%g)", gsl_rng_uniform(r));
    }
  gsl_vector  *v    = apop_query_to_vector("select * from t");
    Diff (apop_var(v) ,apop_query_to_float("select var(vals) from t"),tol6);
    Diff (apop_vector_skew(v) ,apop_query_to_float("select skew(vals) from t"),tol6);
    Diff (apop_vector_kurt(v) ,apop_query_to_float("select kurt(vals) from t"),tol5);
    apop_table_exists("t",1);
}

void test_listwise_delete(){
  int i, j;
  apop_data *t1 = apop_data_calloc(10,10);
  apop_text_alloc(t1, 10, 10);
  for (i=0; i< 10; i++)
      for (j=0; j< 10; j++)
          apop_text_add(t1, i, j, "%i", i*j);
  apop_data *t1c = apop_data_listwise_delete(t1);
  assert(t1c->matrix->size1==10);
  assert(t1c->matrix->size2==10);
  assert(atoi(t1c->text[3][4])==12);

  apop_data *t2 = apop_data_calloc(10); //check in on this form.
  t1->vector    = t2->vector;
  apop_data_set(t1, 4,-1, GSL_NAN);
  apop_data *t2c = apop_data_listwise_delete(t1);
  assert(t2c->matrix->size1==9);
  apop_data_set(t1, 4,-1, GSL_NAN);
  apop_data_set(t1, 7,-1, GSL_NAN);
  apop_data *t3c = apop_data_listwise_delete(t1);
  assert(atoi(t3c->text[3][4])==12);
  assert(atoi(t3c->text[4][4])==20);
  assert(atoi(t3c->text[7][4])==36);
  assert(t3c->matrix->size1==8);
  APOP_COL(t1, 7, v)
  gsl_vector_set_all(v, GSL_NAN);
  assert(!apop_data_listwise_delete(t1));
}

void test_nan_data(){
    apop_text_to_db("test_data_nans", "nandata");
    strcpy(apop_opts.db_name_column, "head");
    strcpy(apop_opts.db_nan, "(nan|\\.)");
  apop_data *d  = apop_query_to_data("select * from nandata");
    apop_opts.output_type ='d';//check that rownames come in OK, and NaNs written right.
    apop_data_print(d, "nantest");
    apop_data_free(d);
  apop_data *d2  = apop_query_to_data("select * from nantest");
    assert(gsl_isnan(apop_data_get_tt(d2,"second", "c")));
    assert(gsl_isnan(apop_data_get_tt(d2,"third", "b")));
    assert(!apop_data_get_tt(d2,"fourth", "b"));
    apop_data_free(d2);
    strcpy(apop_opts.db_nan, "NaN");
}

static void wmt(gsl_vector *v, gsl_vector *v2, gsl_vector *w, gsl_vector *av, gsl_vector *av2, double mean){
    assert(apop_vector_mean(v) == apop_vector_weighted_mean(v,NULL));
    assert(apop_vector_mean(av) == apop_vector_weighted_mean(v,w));
    assert(apop_vector_weighted_mean(v,w) == mean);
    Diff (apop_vector_var(v), apop_vector_weighted_var(v,NULL), tol5);
    Diff (apop_vector_cov(v,v2), apop_vector_weighted_cov(v,v2,NULL), tol5);
    Diff (apop_vector_var(av), apop_vector_weighted_var(v,w), tol5);
    Diff (apop_vector_cov(av,av2), apop_vector_weighted_cov(v,v2,w), tol5);
    Diff (apop_vector_skew_pop(av), apop_vector_weighted_skew(v,w), tol5);
    Diff (apop_vector_kurtosis_pop(av), apop_vector_weighted_kurt(v,w), tol5);
}

void test_weigted_moments(){
  //double        data[]      = {1,2,3};//checking vector_fill with this and w2
  double        alldata[]   = {1,2,3};
  double        data3[]     = {3,2,1};
  double        alldata3[]  = {3,2,1};
  double        weights[]   = {1,1,1};
  gsl_vector    *v          = gsl_vector_alloc(3);
  apop_vector_fill(v, 1, 2, 3);
  gsl_vector    *v2         = apop_array_to_vector(data3, 3);
  gsl_vector    *w          = apop_array_to_vector(weights, 3);
  gsl_vector    *av         = apop_array_to_vector(alldata, 3);
  gsl_vector    *av2        = apop_array_to_vector(alldata3, 3);
    wmt(v,v2,w,av,av2,2);
  double data2[]       = {0,1,2,3,4};
  double alldata2[]    = {0,0,0,0,1,1,1,2,2,3};
  double data4[]       = {0,1,3,2,4};
  double alldata4[]    = {0,0,0,0,1,1,1,3,3,2};
  //double weights2[]    = {4,3,2,1,0};
    v             = apop_array_to_vector(data2, 5);
    v2            = apop_array_to_vector(data4, 5);
    av            = apop_array_to_vector(alldata2, 10);
    av2           = apop_array_to_vector(alldata4, 10);
    gsl_vector    *w2          = gsl_vector_alloc(5);
    apop_vector_fill(w2, 4, 3, 2, 1, 0);
    wmt(v,v2,w2,av,av2,1);
}

void test_split_and_stack(){
APOP_RNG_ALLOC(r, 19);
APOP_VECTOR_ALLOC(dv, 10);
  apop_data *d1 = apop_data_alloc(0,10,10);
int     i,j, tr, tc;
apop_data   **splits, *dv2;
    d1->vector  = dv;
    for(i=-1; i< 10; i++)
        for(j=0; j< 10; j++)
            apop_data_set(d1, j, i, gsl_rng_uniform(r));
    for(i=-1; i< 13; i++){
        splits  = apop_data_split(d1, i, 'r');
        if (i>0 && i< 10)
            assert(splits[0]->matrix->size1 == i);
        else if (i < 0)
            assert(!splits[0]);
        else if (i >= 10)
            assert(splits[0]->matrix->size1 == 10);
        dv2 = apop_data_stack(splits[0], splits[1], 'r');
        for(j=0; j< 50; j++){
            tr  = (int) gsl_rng_uniform(r)*10;
            tc  = (int) gsl_rng_uniform(r)*11-1;
            assert(apop_data_get(dv2, tr, tc) == apop_data_get(d1,tr,tc));
        }
    }
    for(i=-1; i< 13; i++){
        splits  = apop_data_split(d1, i, 'c');
        if (i>0 && i< 10){
            assert(splits[0]->matrix->size2 == i);
            assert(splits[0]->vector->size == 10);
            assert(!splits[1]->vector);
        }
        else if (i < 0){
            assert(!splits[0]);
            assert(!splits[0]);
            assert(splits[1]->vector->size == 10);
        }
        else if (i >= 10){
            assert(splits[0]->matrix->size1 == 10);
            assert(splits[0]->vector->size == 10);
            assert(!splits[1]);
        }
        dv2 = apop_data_stack(splits[0], splits[1], 'c');
        for(j=0; j< 50; j++){
            tr  = (int) gsl_rng_uniform(r)*10;
            tc  = (int) gsl_rng_uniform(r)*11-1;
            assert(apop_data_get(dv2, tr, tc) == apop_data_get(d1,tr,tc));
        }
        apop_data_free(splits[0]);
        apop_data_free(splits[1]);
        apop_data_free(dv2);
    }
}

/** I claim that the mean residual is near zero, and that the predicted
  value is \f$X'\beta\f$.
  */
void test_predicted_and_residual(apop_model *est){
gsl_vector  v,
            *prediction = gsl_vector_alloc(est->data->matrix->size1);
gsl_matrix  *m          = gsl_matrix_alloc(est->data->matrix->size1,est->data->matrix->size2);
    //prep an affine data matrix.
    gsl_matrix_memcpy(m, est->data->matrix);
    v   = gsl_matrix_column(m, 0).vector;
    gsl_vector_set_all(&v, 1);

    apop_data *predict_tab = apop_data_get_page(est->info, "predict");
    v   = gsl_matrix_column(predict_tab->matrix, apop_name_find(predict_tab->names, "residual", 'c')).vector;
    assert(fabs(apop_mean(&v)) < tol5);

    Apop_col_t(predict_tab, "pred", vv);
    gsl_blas_dgemv(CblasNoTrans, 1, m, est->parameters->vector, 0, prediction);
    gsl_vector_sub(prediction, vv);
    assert(fabs(apop_vector_sum(prediction)) < tol5);
}

/** I claim that the F test calculated via apop_F_test(est, NULL, NULL)
 equals a transformation of R^2 (after a normalization step).
*/
void test_f(apop_model *est){
apop_matrix_normalize(est->data->matrix, 'c', 'm');
apop_data *rsq  = apop_estimate_coefficient_of_determination(est);
apop_data *ftab = apop_F_test(est, NULL);
double    n     = est->data->matrix->size1;
double    K     = est->parameters->vector->size;
double    r     = apop_data_get_ti(rsq, "R.squared", 0);
double    f     = apop_data_get(ftab, .rowname="F.stat");
    Diff (f , r*(n-K)/((1-r)*K) , tol5);
}

void test_OLS(gsl_rng *r){
  apop_data *set = apop_data_alloc(0, len, 2);
  int i;
    for(i=0; i< len; i++){
        apop_data_set(set, i, 1, 100*(gsl_rng_uniform(r)-0.5));
        apop_data_set(set, i, 0, -1.4 + apop_data_get(set,i,1)*2.3);
    }
    apop_data *bkup = apop_data_copy(set);
    apop_model *out = apop_estimate(set, apop_ols);
    Diff (apop_data_get(out->parameters, 0,-1) , -1.4 , tol5);
    Diff (apop_data_get(out->parameters, 1,-1) , 2.3 , tol5);

    gsl_vector *w = gsl_vector_alloc(set->matrix->size1);
    gsl_vector_set_all(w, 14);
    bkup->weights  = w;
    out = apop_estimate(bkup, apop_ols);
    Diff (apop_data_get(out->parameters, 0,-1) , -1.4 , tol5);
    Diff (apop_data_get(out->parameters, 1,-1) , 2.3 , tol5);
}

#define INVERTSIZE 100
void test_inversion(gsl_rng *r){
  gsl_matrix  *invme   = gsl_matrix_alloc(INVERTSIZE, INVERTSIZE);
  gsl_matrix  *inved;
  gsl_matrix  *inved_back;
  int         i,j;
  double      error   = 0;
  apop_data *four     = apop_data_alloc(1);
    apop_data_set(four, 0, -1, 4);
  apop_model *fourp    = apop_model_copy(apop_zipf);
    fourp->parameters  = four;
    for(i=0; i<INVERTSIZE; i++)
        for(j=0; j<INVERTSIZE; j++)
            apop_zipf.draw(gsl_matrix_ptr(invme, i,j),r,  fourp);
    apop_det_and_inv(invme, &inved, 0, 1);
    apop_det_and_inv(inved, &inved_back, 0, 1);
    for(i=0; i<INVERTSIZE; i++)
        for(j=0; j<INVERTSIZE; j++)
            error    += gsl_matrix_get(invme, i,j) - gsl_matrix_get(inved_back, i,j);
    assert (error < 1e-5);
}

void test_summarize(){
gsl_matrix      *m;
apop_data       *s;
double          t, v;
    apop_text_to_db("test_data", .has_row_names= 0,1, .tabname = "td");
    m    = apop_query_to_matrix("select * from td");
    s    = apop_data_summarize(apop_matrix_to_data(m));
    //apop_matrix_print(s,"\t", NULL);
    t    = gsl_matrix_get(s->matrix, 1,0);
    assert (t ==3);
    t    = gsl_matrix_get(s->matrix, 2, 1);
    v    = sqrt((2*2 +3*3 +3*3 +4.*4.)/3.);
    assert (t == v) ;
}

void test_dot(){
apop_data *d1   = apop_text_to_data(.text_file="test_data2",0,1); // 100 x 2
apop_data *d2   = apop_text_to_data("test_data2"); // 100 x 2
apop_data *d3   = apop_dot(d1, d2, 0, 1);
//apop_data *d4   = apop_dot(d1, d2, 'p', 0);
gsl_vector  v1 = gsl_matrix_row(d1->matrix, 0).vector;
apop_data *d5   = apop_vector_to_data(&v1); // 2 x 1
apop_data *d7   = apop_dot(d5, d5, 0, 0);
    assert(apop_data_get(d7, 0, -1) == apop_data_get(d3, 0,0));
apop_data *d8   = apop_dot(d1, d5, 0, 0);
apop_data *d9   = apop_dot(d5, d1, .form2=1);
    gsl_vector_sub(d8->vector, d9->vector);
    assert(!apop_vector_sum(d8->vector));
}
 
static void fill_p(apop_data *d, gsl_rng *r){
    int j, k;
    if (d->vector)
        for (j=0; j< d->vector->size; j++)
            gsl_vector_set(d->vector, j, gsl_rng_uniform(r));
    if (d->matrix)
        for (j=0; j< d->matrix->size1; j++)
            for (k=0; k< d->matrix->size2; k++)
                gsl_matrix_set(d->matrix, j, k, gsl_rng_uniform(r));
    if (d->weights)
        for (j=0; j< d->weights->size; j++)
            gsl_vector_set(d->weights, j, gsl_rng_uniform(r));
}

static void check_p(apop_data *d, apop_data *dout){
    int j, k;
    if (d->vector)
        for (j=0; j< d->vector->size; j++)
            assert(dout->vector->data[j] == d->vector->data[j]);
    if (d->weights)
        for (j=0; j< d->weights->size; j++)
            assert(dout->weights->data[j] == d->weights->data[j]);
    if (d->matrix)
        for (j=0; j< d->matrix->size1; j++)
            for (k=0; k< d->matrix->size2; k++)
                assert(gsl_matrix_get(d->matrix, j, k) == gsl_matrix_get(dout->matrix, j, k));
}

void apop_pack_test(gsl_rng *r){
  int i, v, m1,m2, w;
  apop_data *d, *dout, *p2, *outp2;
  gsl_vector *mid;
    for (i=0; i< 10; i++){
        v   = gsl_rng_uniform(r) > 0.5 ? gsl_rng_uniform(r)*100 : 0;
        m1  = gsl_rng_uniform(r)*100;
        m2  = gsl_rng_uniform(r) > 0.5 ? gsl_rng_uniform(r)*100 : 0;
        w   = gsl_rng_uniform(r) > 0.5 ? gsl_rng_uniform(r)*100 : 0;
        if (!v && !w && (!m1 || !m2))
            continue; //I actually get this unlucky draw.
        d   = apop_data_alloc(v, m1, m2);
        dout    = apop_data_alloc(v, m1, m2);
        if (w) {d->weights = gsl_vector_alloc(w);
                dout->weights = gsl_vector_alloc(w);}
        fill_p(d, r);
        int second_p = i %2;
        if (second_p){
            p2 = apop_data_add_page(d, apop_data_alloc(v, m1, m2), "second p");
            outp2 = apop_data_add_page(dout, apop_data_alloc(v, m1, m2), "second p");
            if (w) {p2->weights = gsl_vector_alloc(w);
                    outp2->weights = gsl_vector_alloc(w);}
            fill_p(p2, r);
        }
        mid     = apop_data_pack(d, .all_pages= second_p ? 'y' : 'n');
        apop_data_unpack(mid, dout);
        //apop_data_unpack(mid, dout, .all_pages= second_p ? 'y' : 'n');
        check_p(d, dout);
        if (second_p)
            check_p(d->more, dout->more);
        if (mid) gsl_vector_free(mid); 
        apop_data_free(d); apop_data_free(dout);
        }
}

void test_model_fix_parameters(gsl_rng *r){
  size_t    i, ct = 1000;
  apop_data *d  = apop_data_alloc(0,ct,2);
  double    draw[2];
  apop_multivariate_normal.vbase =
  apop_multivariate_normal.m1base =
  apop_multivariate_normal.m2base = 2;
  apop_model *pp = apop_model_set_parameters(apop_multivariate_normal,
                                        8, 1, 0.5,
                                        2, 0.5, 1);
    for(i=0; i< ct; i++){
        apop_multivariate_normal.draw(draw, r, pp);
        apop_data_set(d, i, 0, draw[0]);
        apop_data_set(d, i, 1, draw[1]);
    }

    apop_data *pcopy = apop_data_copy(pp->parameters);
    gsl_matrix_set_all(pp->parameters->matrix, GSL_NAN);
    apop_model *mep1   = apop_model_fix_params(pp);
    Apop_settings_add(mep1, apop_mle, starting_pt, ((double[]){1.5, .25, .25, 1.5}));
    apop_model   *e1  = apop_estimate(d, *mep1);
    gsl_vector_sub(e1->parameters->vector, pcopy->vector);
    assert(apop_vector_sum(e1->parameters->vector) < 1e-1);

  double    start2[] = {7,3};
    pp->parameters = apop_data_copy(pcopy);
    gsl_vector_set_all(pp->parameters->vector, GSL_NAN);
    apop_model *mep2   = apop_model_fix_params(pp);
    Apop_settings_add(mep2, apop_mle, starting_pt, start2);
    Apop_settings_add(mep2, apop_mle, method, APOP_CG_PR);
    apop_model   *e2  = apop_estimate(d, *mep2);
    gsl_matrix_sub(e2->parameters->matrix, pcopy->matrix);
    assert(apop_matrix_sum(e2->parameters->matrix) < 1e-2);
}

void test_linear_constraint(){
  gsl_vector *beta      = gsl_vector_alloc(2);
    gsl_vector_set(beta, 0, 7);
    gsl_vector_set(beta, 1, 7);
  apop_data *contrasts  = apop_data_calloc(1,1,2);
    apop_data_set(contrasts, 0, 0, -1);
    apop_data_set(contrasts, 0, 1, -1);
    Diff (apop_linear_constraint(beta, contrasts, 0) , sqrt(2*49) , tol5);
    assert(!apop_vector_sum(beta));
    gsl_vector_set(beta, 0, 0);
    gsl_vector_set(beta, 1, 7);
    Diff (apop_linear_constraint(beta, contrasts, 0) , sqrt(49/2.) , tol5);
    assert(!apop_vector_sum(beta));
    assert(gsl_vector_get(beta,0)==-7/2.);
    //inside corner: find the corner
  gsl_vector *beta2     = gsl_vector_alloc(3);
    gsl_vector_set(beta2, 0, 7);
    gsl_vector_set(beta2, 1, 7);
    gsl_vector_set(beta2, 2, 7);
  apop_data *contrasts2  = apop_data_calloc(3,3,3);
    apop_data_set(contrasts2, 0, 0, -1);
    apop_data_set(contrasts2, 1, 1, -1);
    apop_data_set(contrasts2, 2, 2, -1);
    Diff (apop_linear_constraint(beta2, contrasts2, 0) , sqrt(3*49) , tol5);
    assert(apop_vector_sum(beta2)==0);
    //sharp corner: go to one wall.
    gsl_vector_set(beta2, 0, 7);
    gsl_vector_set(beta2, 1, 7);
    gsl_vector_set(beta2, 2, 7);
    apop_data_set(contrasts2, 0, 1, 1);
    Diff(apop_linear_constraint(beta2, contrasts2, 0) , sqrt(2*49), tol5);
    assert(gsl_vector_get(beta2,0)==7);
    assert(gsl_vector_get(beta2,1)==0);
    assert(gsl_vector_get(beta2,2)==0);
}

void test_data_sort(){
  int       i;
  apop_data *d = apop_text_to_data("test_data2", 0 ,1);
    apop_data_sort(d, 0, 'a');
    for (i=1; i< d->matrix->size2; i++){
        assert(apop_data_get(d, i,0) >= apop_data_get(d, i-1,0));
    }
    assert(apop_data_get(d, 0,1)== 32 || apop_data_get(d, 0,1)== 9);

    apop_data_sort(d, 0, 'd');
    for (i=1; i< d->matrix->size2; i++){
        assert(apop_data_get(d, i,0) <= apop_data_get(d, i-1,0));
    }
    assert(apop_data_get(d, 0,1)== 55);
}

void test_histograms(gsl_rng *r){
//create a million draws 
    int n = 5e5;
  double    mu = 2.8, sigmasq = 1.34;
  apop_data *d = apop_data_alloc(0,n,1);
  gsl_matrix *out   = gsl_matrix_alloc(n,1);
  int       i;
    for (i=0; i< n; i++)
        //gsl_matrix_set(d->matrix, i, 0, gsl_rng_uniform(r));
        gsl_matrix_set(d->matrix, i, 0, gsl_ran_gaussian(r, sqrt(sigmasq))+mu);
    apop_model   *hp = apop_estimate(d, apop_histogram);
    for (i=0; i< n; i++){
        apop_histogram.draw(gsl_matrix_ptr(out, i,0), r, hp);
        assert(gsl_finite(gsl_matrix_get(out, i,0)));
    }
    apop_model *outparams   = apop_estimate(apop_matrix_to_data(out), apop_normal);
    assert(fabs(outparams->parameters->vector->data[0]-mu) < 1e2);
    assert(fabs(outparams->parameters->vector->data[1]-sigmasq) < 1e2);
    apop_model_free(hp);
    apop_model_free(outparams);

    apop_model   *hp2 = apop_model_copy(apop_histogram);
    Apop_model_add_group(hp2, apop_histogram, .data=d, .bins_in=100);
    for (i=0; i< n; i++){
        apop_draw(gsl_matrix_ptr(out, i,0), r, hp);
        assert(gsl_finite(gsl_matrix_get(out, i,0)));
    }
    apop_model *outparams2   = apop_estimate(apop_matrix_to_data(out), apop_normal);
    assert(fabs(outparams2->parameters->vector->data[0]-mu) < 1e2);
    assert(fabs(outparams2->parameters->vector->data[1]-sigmasq) < 1e2);
    apop_model_free(hp2);
    apop_model_free(outparams2);
//    apop_plot_histogram(out, 100, NULL); 
}


void test_jackknife(){
  double      pv[] = {3.09,2.8762};
  apop_model m = apop_normal;
  APOP_DATA_ALLOC(d, len, 1);
  APOP_RNG_ALLOC(r, 8);
  apop_data  *p  = apop_data_alloc();
  p->vector      = apop_array_to_vector(pv, 2);
  apop_model*pp = apop_model_copy(m);
      pp->parameters    = p;
  size_t      i;
    for (i =0; i< len; i++)
        m.draw(apop_data_ptr(d, i, 0), r, pp); 
    apop_data *out = apop_jackknife_cov(d, *pp);
    //apop_data_show(out);
    //printf("%g\n",  2*gsl_pow_2(pv[1])/(len-1));
    //fflush(NULL);
    //Notice that the jackknife just ain't a great estimator here.
assert ((fabs(apop_data_get(out, 0,0) - gsl_pow_2(pv[1])/len)) < tol2 
            && fabs(apop_data_get(out, 1,1) - gsl_pow_2(pv[1])/(2*len)) < tol2*100);
    apop_data *out2 = apop_bootstrap_cov(d, m);
    assert (fabs(apop_data_get(out2, 0,0) - gsl_pow_2(pv[1])/len) < tol2
                && fabs(apop_data_get(out2, 1,1) - gsl_pow_2(pv[1])/(2*len)) < tol2);
    apop_data_free(d);
}

//In my inattention, I wrote two jackknife tests. So you get double the checks.
int test_jack(gsl_rng *r){
  int i, draws     = 2000;
  apop_data *d  =apop_data_alloc(0, draws, 1);
  apop_model  m   = apop_normal;
  double      pv[] = {1., 3.};
    m.parameters = apop_line_to_data(pv, 2,0,0);
    for (i =0; i< draws; i++)
        m.draw(apop_data_ptr(d, i, 0), r, &m); 
    apop_data *out = apop_jackknife_cov(d, m);
    double error = fabs(apop_data_get(out, 0,0)-gsl_pow_2(pv[1])/draws) //var(mu)
                + fabs(apop_data_get(out, 1,1)-gsl_pow_2(pv[1])/(2*draws))//var(sigma)
                +fabs(apop_data_get(out, 0,1)) +fabs(apop_data_get(out, 1,0));//cov(mu,sigma); should be 0.
    return (error < 1e-3);
}

void test_lognormal(gsl_rng *r){
    int j;
    apop_model *source = apop_model_copy(apop_normal);
    apop_model_clear(NULL, source);
    double mu    = gsl_ran_flat(r, -1, 1);
    double sigma = gsl_ran_flat(r, .01, 1);
    int n     = gsl_ran_flat(r,1,8e5);
    apop_data *data = apop_data_alloc(0,1,n);
    gsl_vector_set(source->parameters->vector, 0, mu);
    gsl_vector_set(source->parameters->vector, 1, sigma);
    for (j=0; j< n; j++){
        double *k   = gsl_matrix_ptr(data->matrix, 0, j);
        apop_draw(k, r, source);
        *k = exp(*k);
    }
    apop_model *out = apop_estimate(data, apop_lognormal);
    double muhat = apop_data_get(out->parameters, 0,-1);
    double sigmahat = apop_data_get(out->parameters, 1,-1);
    if (verbose) printf("mu: %g, muhat: %g, var: %g, varhat: %g\n", mu, muhat,  sigma,sigmahat);
    assert(fabs(mu-muhat)<1e-2);
    assert(fabs(sigma-sigmahat)<1e-2);
}

void test_multivariate_normal(gsl_rng *r){
  int       len = 4e5;
  size_t    i;
  apop_data *rdraws = apop_data_alloc(0, len, 2);
  double params[] = {1, 3, 0,
                    2, 0, 1};
    apop_data *p = apop_line_to_data(params, 2,2,2);
    apop_model *mv = apop_model_copy(apop_multivariate_normal);
    mv->parameters=p;
    for(i=0; i < len; i ++){
        APOP_ROW(rdraws, i, ping);
        apop_draw(ping->data, r, mv);
    }
    apop_model *est =apop_estimate(rdraws, apop_multivariate_normal);
    double error = fabs(est->parameters->vector->data[0] - p->vector->data[0])
                  +fabs(est->parameters->vector->data[1] - p->vector->data[1])
                  +fabs(est->parameters->matrix->data[0] - p->matrix->data[0])
                  +fabs(est->parameters->matrix->data[1] - p->matrix->data[1])
                  +fabs(est->parameters->matrix->data[2] - p->matrix->data[2])
                  +fabs(est->parameters->matrix->data[3] - p->matrix->data[3]);
    assert (error < 3e-2); //yes, unimpressive, but we don't wanna be here all day.
}

static void common_binomial_bit(apop_model *out, int n, double p){
    /*double phat = apop_data_get(out->parameters, 1,-1);
    double nhat = apop_data_get(out->parameters, 0,-1);
    if (verbose) printf("n: %i, p: %g, nhat: %g, phat: %g\n", n, p, phat, nhat);*/
    assert(apop_data_get(out->parameters, 0,-1) == n);
    assert(apop_data_get(out->parameters, 1,-1) - p < 1e-2);
}

void test_binomial(gsl_rng *r){
    size_t  i;
    double p = gsl_rng_uniform(r);
    int n     = gsl_ran_flat(r,1,1e5);
    apop_data *d = apop_data_alloc(1,n);
    for(i=0; i < n; i ++)
        apop_data_set(d, 0,i,(gsl_rng_uniform(r) < p));
    apop_model *out = apop_estimate(d, apop_binomial);
    apop_model *outm = apop_estimate(d, apop_multinomial);
    common_binomial_bit(out, n, p);
    common_binomial_bit(outm, n, p);
    apop_data_free(d);

    /*
    p = gsl_rng_uniform(r);
    n  = gsl_ran_flat(r,1,4e4);
    d = apop_data_calloc(0,n,2);
    for(i=0; i < n; i ++){
        if (gsl_rng_uniform(r) < p)
            apop_matrix_increment(d->matrix, i,1);
            //apop_matrix_increment(d->matrix, (int)gsl_ran_flat(r,0, n),1,1);
        else
            apop_matrix_increment(d->matrix, i,0);
            //apop_matrix_increment(d->matrix, (int)gsl_ran_flat(r,0, n),0,1);
    }
    apop_model *bint = apop_model_copy(apop_binomial);
    Apop_settings_add_group(bint, apop_rank, NULL);
    out = apop_estimate(d, *bint);
    common_binomial_bit(out, n, p);
    apop_model *mint = apop_model_copy(apop_multinomial);
    Apop_settings_add_group(mint, apop_rank, NULL);
    out = apop_estimate(d, *mint);
    common_binomial_bit(out, n, p);
    */
    apop_data_free(d);
}

void db_to_text(){
    apop_db_open(NULL);
    if (!apop_table_exists("d"))
        apop_text_to_db("data-mixed", "d", 0, 1, NULL);
    apop_data *d = apop_query_to_mixed_data ("tmttmmmt", "select * from d");
    int b_allele_col = apop_name_find(d->names, "b_all.*", 't');
    assert(!strcmp("T",  d->text[3][b_allele_col]));
    int rsid_col = apop_name_find(d->names, "rsid", 't');
    assert(!strcmp("rs2977656",  d->text[4][rsid_col]));
    assert(apop_data_get_it(d, 5, "ab")==201);

    apop_data *dcc = apop_data_copy(d); //test apop_data_copy
    assert(!strcmp("T",  dcc->text[3][b_allele_col]));
    assert(!strcmp("rs2977656",  dcc->text[4][rsid_col]));
    assert(apop_data_get_it(dcc, 5, "ab")==201);

    apop_data *dd = apop_query_to_text ("select * from d");
    b_allele_col = apop_name_find(dd->names, "b_all.*", 't');
    assert(!strcmp("T",  dd->text[3][b_allele_col]));
    rsid_col = apop_name_find(dd->names, "rsid", 't');
    assert(!strcmp("rs2977656",  dd->text[4][rsid_col]));
    
    apop_data *dc = apop_data_copy(d);
    b_allele_col = apop_name_find(dc->names, "b_all.*", 't');
    assert(!strcmp("T",  dc->text[3][b_allele_col]));
    rsid_col = apop_name_find(dc->names, "rsid", 't');
    assert(!strcmp("rs2977656",  dc->text[4][rsid_col]));
    assert(apop_data_get_it(dc, 5, "ab")==201);

    char oldtype = apop_opts.output_type;
    apop_opts.output_type = 'd';
    apop_data_print(dc, "mixedtest");
    apop_data *de = apop_query_to_mixed_data("mmmmtttt","select * from mixedtest");
    b_allele_col = apop_name_find(de->names, "b_all.*", 't');
    assert(!strcmp("T",  de->text[3][b_allele_col]));
    rsid_col = apop_name_find(de->names, "rsid", 't');
    assert(!strcmp("rs2977656",  de->text[4][rsid_col]));
    assert(apop_data_get_it(de, 5, "ab")==201);
    apop_opts.output_type = oldtype;
}

void test_blank_db_queries(){
    apop_db_open(NULL);
    apop_query("create table t (a, b, c)");
    apop_data *d = apop_query_to_data("select * from t");
    apop_data *e = apop_query_to_text("select * from t");
    gsl_matrix *f = apop_query_to_matrix("select * from t");
    gsl_vector *g = apop_query_to_vector("select * from t");
    double h = apop_query_to_float("select * from t");
    assert(d==NULL);
    assert(e==NULL);
    assert(f==NULL);
    assert(g==NULL);
    assert(gsl_isnan(h));
}

int get_factor_index(apop_data *flist, char *findme){
  int i;
    for (i=0; i< flist->textsize[0]; i++)
        if (apop_strcmp(flist->text[i][0], findme))
            return i;
    return -2;
}

//If the dummies are a separate matrix, offset=0;
//If the dummies are an addendum to main, offset=original_data->matrix->size2;
static void check_for_dummies(apop_data *d, apop_data *dum, int offset){
  int i, j, n;
    apop_data *factorlist = apop_data_get_page(d, "categor");
    for(i=0; i < d->textsize[0]; i ++)
        if ((n = get_factor_index(factorlist, d->text[i][0]))>0){
            for(j=0; j < factorlist->textsize[0]-1; j ++)
                if (j==n-1)
                    assert(apop_data_get(dum, i, j+offset));
                else
                    assert(!apop_data_get(dum, i, j+offset));
        } else
            for(j=0; j < factorlist->textsize[0]-1; j ++)
                assert(!apop_data_get(dum, i, j+offset));
}

void dummies_and_factors(){
  int i;
    apop_text_to_db("data-mixed", "genes");
    apop_data *d = apop_query_to_mixed_data("mmmt", "select aa, bb, 1, a_allele from genes");
    apop_data *dum = apop_data_to_dummies(d, 0, 't', 0);
    check_for_dummies(d, dum, 0);
    apop_text_to_factors(d, 0, 2);
    for(i=0; i < d->textsize[0]; i ++) //the set is only As and Cs.
        if (!strcmp(d->text[i][0], "A"))
            assert(apop_data_get(d, i, 2) == 0);
        else
            assert(apop_data_get(d, i, 2) == 1);
    //test combination routines
    apop_data *d2 = apop_query_to_mixed_data("mmmt", "select aa, bb, 1, a_allele from genes");
    apop_data_to_dummies(d2, 0,  't', .append='y');
    check_for_dummies(d2, d2, 3);
}

void test_vector_moving_average(){
  int   i;
  gsl_vector *v = apop_vector_realloc(NULL, 100); //using realloc as an alloc
    for(i=0; i < 100; i ++)
        gsl_vector_set(v, i, i);
    gsl_vector *unsmooth = apop_vector_moving_average(v, 1);
    for(i=0; i < 100; i ++)
        assert(gsl_vector_get(v, i) == gsl_vector_get(unsmooth, i));
    gsl_vector *slightly_smooth = apop_vector_moving_average(v, 2);
    //With evenly-spaced data, a moving average returns the original,
    //with tails missing:
    for(i=0; i < 98; i ++)
        assert(gsl_vector_get(v, i+1) == gsl_vector_get(slightly_smooth, i));
}

void test_transpose(){
    apop_data *t = apop_text_to_data("test_data", 0, 1);
    apop_data *tt = apop_data_transpose(t);
    assert(apop_data_get(tt, 0, 3) == 9);
    assert(apop_data_get(tt, 1, 0) == 4);
    assert(!strcmp(tt->names->row[2], "c"));
    assert(!tt->names->colct);
}

apop_data *generate_probit_logit_sample (gsl_vector* true_params, gsl_rng *r, apop_model *method){
  int i, j;
  double val;
  int samples = 1e5;
  apop_data *data = apop_data_alloc(samples, true_params->size);
        //generate a random vector of data X, then set the outcome to one if the probit/logit condition holds
        for (i = 0; i < samples; i++){
            apop_data_set(data, i, 0, 1);
            for (j = 1; j < true_params->size; j++)
                apop_data_set(data, i, j, (gsl_rng_uniform(r)-0.5) *2);
            APOP_ROW(data, i, asample);
            gsl_blas_ddot(asample, true_params, &val);
            if (method == &apop_probit)
                apop_data_set(data, i, 0, (gsl_ran_gaussian(r, 1) > -val));
            else   //Logit:   p(act) = e(xb) / (1+e(xb));
                apop_data_set(data, i, 0, gsl_rng_uniform(r) < exp(val)/(1+exp(val)));
        }
    return data;
}

void test_unique_elements(){
    double d[] = {0, -3, 1.2, 2.4, -2, -3, 0.1, -0.1, 1.2, -2};
    gsl_vector *dv = apop_array_to_vector(d, sizeof(d)/sizeof(double));
    gsl_vector *distinct = apop_vector_unique_elements(dv);
    assert(distinct->size == 7);
    assert(gsl_vector_get(distinct, 2) == -.1);
    assert(gsl_vector_get(distinct, 3) == 0);
    assert(gsl_vector_get(distinct, 4) == .1);

    apop_data *t = apop_text_alloc(NULL, 9, 7);
    apop_text_add(t, 0, 0, "Hi,");
    apop_text_add(t, 1, 0, "there");
    apop_text_add(t, 2, 0, ".");
    apop_text_add(t, 3, 0, "This");
    apop_text_add(t, 4, 0, "there");
    apop_text_add(t, 5, 0, "is");
    apop_text_add(t, 6, 0, "dummy");
    apop_text_add(t, 7, 0, "text");
    apop_text_add(t, 8, 0, ".");
    apop_data *dt = apop_text_unique_elements(t, 0);
    assert(dt->textsize[0] == 7);
    assert(!strcmp(".", dt->text[0][0]));
    assert(!strcmp("Hi,", dt->text[1][0]));
    assert(!strcmp("text", dt->text[5][0]));
}

void test_probit_and_logit(gsl_rng *r){
  int i;
  for (i=0; i < 3; i++){
    int param_ct = gsl_rng_uniform(r)*7 + 1; //up to seven params.
    int j;
    gsl_vector *true_params = gsl_vector_alloc(param_ct);
    for (j = 0; j < param_ct; j++)
        gsl_vector_set(true_params, j, (gsl_rng_uniform(r)-0.5)*2);
    //apop_vector_show(true_params);

    //Logit
    apop_data* data = generate_probit_logit_sample(true_params, r, &apop_logit);
    Apop_model_add_group(&apop_logit, apop_mle, .want_cov='n');
    apop_model *m = apop_estimate(data, apop_logit);
    APOP_COL(m->parameters, 0, logit_params);
    assert(apop_vector_distance(logit_params, true_params) < 0.07);
    apop_data_free(data);
    apop_model_free(m);

    //Probit
    apop_data* data2 = generate_probit_logit_sample(true_params, r, &apop_probit);
    Apop_model_add_group(&apop_probit, apop_mle, .want_cov='n');
    m = apop_estimate(data2, apop_probit);
    APOP_COL(m->parameters, 0, probit_params);
    assert(apop_vector_distance(probit_params, true_params) < 0.07);
    apop_model_free(m);
    apop_data_free(data2);
  }
}

void test_resize(){
    int i;
    //This is the multiplication table from _Modeling with Data_
    //with a +.1 to distinguish columns from rows.
  gsl_matrix *m = apop_matrix_realloc(NULL, 20,15);//check using realloc as an alloc
    gsl_matrix_set_all(m, 1);
    for (i=0; i< m->size1; i++){
        Apop_matrix_row(m, i, one_row);
        gsl_vector_scale(one_row, i+1);
    }
    for (i=0; i< m->size2; i++){
        Apop_matrix_col(m, i, one_col);
        gsl_vector_scale(one_col, i+1);
        gsl_vector_add_constant(one_col, (i+1)/10.);
    }
    apop_matrix_realloc(m, 11, 17);
    assert(gsl_matrix_get(m, 3, 5) == 4*6+.6);
    apop_matrix_realloc(m, 10, 10);
    Diff (apop_matrix_sum(m) , 55 * 56 , tol6);
    gsl_vector *v = gsl_vector_alloc(20);
    for (i=0; i< 20; i++)
        gsl_vector_set(v, i, i);
    apop_vector_realloc(v, 38);
    for (i=0; i< 20; i++)
        assert(gsl_vector_get(v, i) == i);
    apop_vector_realloc(v, 10);
    assert(apop_vector_sum(v) == 45);
}

void test_crosstabbing() {
    apop_db_close(); //gotta test it somewhere
    if (!apop_table_exists("snps"))
        apop_text_to_db("test_data_mixed", "snps", 0, 1);
    apop_query("create table snp_ct as "
                 " select a_allele, b_allele, count(*) as ct "
                 " from snps group by a_allele, b_allele ");
    apop_data *d = apop_db_to_crosstab("snp_ct", "a_allele", "b_allele", "ct");
    assert(apop_data_get_tt(d, "A", "G")==5);
    assert(apop_data_get_tt(d, "C", "G")==1);
}

void test_mvn_gamma(){
    assert(apop_multivariate_gamma(10, 1)==gsl_sf_gamma(10));
    assert(apop_multivariate_lngamma(10, 1)==gsl_sf_lngamma(10));
}

void test_data_to_db() {
  int i, j;
    if (!apop_table_exists("snps"))
        apop_text_to_db("test_data_mixed", "snps");
    apop_data *d = apop_query_to_mixed_data("tvttmmmt", "select * from snps");
    apop_data_to_db(d, "snps2");
    apop_data *d2 = apop_query_to_mixed_data("vmmmtttt", "select * from snps2");
    for (i=0; i< d2->vector->size; i++)
        assert(d->vector->data[i] == d2->vector->data[i]);
    for (i=0; i< d2->matrix->size1; i++)
        for (j=0; j< d2->matrix->size2; j++)
            assert(gsl_matrix_get(d->matrix, i, j) ==  gsl_matrix_get(d2->matrix, i, j));
    for (i=0; i< d2->textsize[0]; i++)
        for (j=0; j< d2->textsize[1]; j++)
            assert(!strcmp(d->text[i][j],d2->text[i][j]));  
}

void test_default_rng(gsl_rng *r) {
    size_t i;
    gsl_vector *o = gsl_vector_alloc(2e5);
    apop_model *ncut = apop_model_set_parameters(apop_normal, 1.1, 1.23);
    ncut->draw = NULL; //forced to use the default.
    for(i=0; i < 2e5; i ++)
        apop_draw(o->data+i, r, ncut);
    apop_model *back_out = apop_estimate(apop_vector_to_data(o), apop_normal);
    Diff(back_out->parameters->vector->data[0] , 1.1 , tol2);
    Diff(back_out->parameters->vector->data[1] , 1.23 , tol2);
    gsl_vector_free(o);
}

double ran_uniform(double in, void *r){ return gsl_rng_uniform(r);}

void test_posdef(gsl_rng *r){
    r_global = r;
    size_t j;
    for(j=0; j < 30; j ++){
        int size = gsl_rng_uniform(r) *10+1;
        apop_data *d = apop_data_alloc( size, size);
        apop_map(d, .fn_dp=ran_uniform, .param=r, .inplace=1, .part='m');
        apop_matrix_to_positive_semidefinite(d->matrix);
        assert(apop_matrix_is_positive_semidefinite(d->matrix));
    }
}

static double set_to_index(double in, int index){ return index;}
static double is_even(double in){ return !((int)in%2);}
static double is_odd(double in){ return (int)in%2;}
static double nan_even(double in){ return is_even(in) ? GSL_NAN : in; }

void row_manipulations(){
    apop_data *test= apop_data_alloc(10);
    apop_map(test, .fn_di=set_to_index, .part='v', .inplace='y');
    int rm[10] = {0,1,0,1,0,1,0,1,0,1};
    apop_data_rm_rows(test, rm);
    assert (!apop_map_sum(test, .fn_d=is_odd, .part='v'));

    apop_data *test2= apop_data_alloc(10,0,0);
    apop_map(test2, .fn_di=set_to_index, .part='v', .inplace='y');
    apop_map(test2, .fn_d=nan_even, .part='v', .inplace='y');
    apop_data_listwise_delete(test2, 'y');
    assert (5== apop_map_sum(test2, .fn_d=is_odd, .part='v'));
    assert (!apop_map_sum(test2, .fn_d=is_even, .part='v'));
}

void estimate_model(apop_data *data, apop_model *dist, int method, apop_data *true_params){
  double                  starting_pt[] =  {1.6, 1.4};//{3.2, 1.4}; //{1.82,2.1};

    Apop_model_add_group(dist, apop_mle, 
        .parent       = dist,  .starting_pt = starting_pt,
        .method       = method, .verbose   =0,
        .step_size    = 1e-1,
        .tolerance    = 1e-2,   .k         = 1.8,
        .t_initial    = 1,      .t_min     = .5,
        .use_score    = 1,      .want_cov  = 'n'
        );
    Apop_model_add_group(dist, apop_parts_wanted);
    apop_model *e    = apop_estimate(data,*dist);
    Diff(0.0, apop_vector_distance(apop_data_pack(true_params),apop_data_pack(e->parameters)), 1e-1); 
}

/*Produce random data, then try to recover the original params */
  void test_one_distribution(gsl_rng *r, apop_model *model, apop_model *true_params){
  long int        runsize             = 1e5;
  apop_data      *data; 
  size_t          i;
    //generate.
    data = apop_data_calloc(runsize,model->dsize);
    if (!strcmp(model->name, "Wishart distribution")){
        data = apop_data_calloc(runsize,4);
        true_params->parameters->vector->data[0] = runsize-4;
        for (i=0; i< runsize; i++){
            Apop_row(data, i, v)
            true_params->draw(v->data, r, true_params);
            assert(!apop_vector_map_sum(v, nan_map));
        }
    } else {
        data = apop_data_calloc(runsize,model->dsize);
        for (i=0; i< runsize; i++){
            Apop_row(data, i, v)
            true_params->draw(v->data, r, true_params);
            assert(!apop_vector_map_sum(v, nan_map));
        }
    }
    //model.score =NULL;
    estimate_model(data, model,APOP_SIMPLEX_NM, true_params->parameters);
    estimate_model(data, model,APOP_CG_PR, true_params->parameters);
    apop_data_free(data);
}

void test_cdf(gsl_rng *r, apop_model *m){//m is parameterized
    //Make random draws from the dist, then find the CDF at that draw
    //That should generate a uniform distribution.
    if (!m->cdf || apop_strcmp(m->name, "Bernoulli distribution")
                || apop_strcmp(m->name, "Binomial distribution"))
        return;
    int i, drawct = 1e5;
    apop_data *draws = apop_data_alloc(drawct, m->dsize);
    apop_data *cdfs = apop_data_alloc(drawct);
    for (i=0; i< drawct; i++){
        Apop_row(draws, i, onerow);
        apop_draw(onerow->data, r, m);
        Apop_data_row(draws, i, one_data_pt);
        apop_data_set(cdfs, i, -1, apop_cdf(one_data_pt, m));
    }
    apop_model *h = apop_estimate(cdfs, apop_histogram);
    gsl_histogram *histo = Apop_settings_get(h, apop_histogram, pdf);
    for (i=0; i< histo->n; i++)
        if(fabs(histo->bin[i]- (drawct+0.0)/histo->n)> drawct/1000.)
        printf("%g %g %g\n", histo->bin[i], (drawct+0.0)/histo->n, drawct/1000.);
        //Diff(histo->bin[i], (drawct+0.0)/histo->n, drawct/100.);
}


double  true_parameter_v[]    = {1.82,2.1};

void test_distributions(gsl_rng *r){
  if (verbose) printf("\n");
  int         i;
  apop_model* true_params;
  apop_model  null_model      = {"the null model"};
  apop_model  dist[]          = {
                apop_beta, apop_bernoulli, apop_binomial, 
               /* apop_chi_squared,*/
                apop_dirichlet, apop_exponential, 
               /* apop_f_distribution,*/
                apop_gamma, 
                apop_lognormal, apop_multinomial, apop_multivariate_normal,
                apop_normal, apop_poisson,
                /*apop_t_distribution, */ apop_uniform,
                apop_waring, /*apop_wishart,*/
                apop_yule, apop_zipf, null_model};

    for (i=0; !apop_strcmp(dist[i].name, "the null model"); i++){
        if (verbose) {printf("\t%s: ", dist[i].name); fflush(NULL);}
        true_params   = apop_model_copy(dist[i]);
        true_params->parameters = apop_line_to_data(true_parameter_v, dist[i].vbase==1 ? 1 : 2,0,0);
        if (apop_strcmp(dist[i].name, "Dirichlet distribution"))
            dist[i].dsize=2;
        if (apop_strcmp(dist[i].name, "Beta distribution"))
            true_params->parameters = apop_line_to_data((double[]){.5, .2} , 2,0,0);
        if (apop_strcmp(dist[i].name, "Bernoulli distribution"))
            true_params->parameters = apop_line_to_data((double[]){.1} , 1,0,0);
        if (apop_strcmp(dist[i].name, "Binomial distribution")){
            true_params->parameters = apop_line_to_data((double[]){15, .2} , 2,0,0);
            dist[i].dsize=15;
        }
        if (apop_strcmp(dist[i].name, "Multivariate normal distribution")){
            true_params->parameters = apop_line_to_data((double[]){15, .5, .2,
                                                                    3, .2, .5} , 2,2,2);
            dist[i].dsize=2;
        }
        if (apop_strcmp(dist[i].name, "Multinomial distribution")){
            true_params->parameters = apop_line_to_data((double[]){15, .5, .2, .1} , 4,0,0);
            dist[i].dsize=15;
        }
        if (apop_regex(dist[i].name, "gamma distribution"))
            true_params->parameters = apop_line_to_data((double[]){1.5, 2.5} , 2,0,0);
        if (apop_strcmp(dist[i].name, "t distribution"))
            true_params->parameters = apop_line_to_data((double[]){1, 3, 996} , 3,0,0);
        if (apop_strcmp(dist[i].name, "Wishart distribution")){
            true_params->parameters = apop_line_to_data((double[]){996, .2, .1,
                                                                     0, .1, .2}, 2,2,2);
            apop_vector_realloc(true_params->parameters->vector, 1);
        }
        test_one_distribution(r, dist+i, true_params);
        test_cdf(r, true_params);
        printf("Passed.\n");
    }
}

void test_pmf(){
   double x[] = {0, 0.2, 0 , 0.4, 1, .7, 0 , 0, 0};
   size_t i;
	apop_data *d= apop_data_alloc();
 	double out;
	gsl_rng *r = apop_rng_alloc(1234);
	d->vector = apop_array_to_vector(x, 9);
	apop_model *m = apop_crosstab_to_pmf(d);
	gsl_vector *v = gsl_vector_alloc(d->vector->size);
	for (i=0; i< 1e5; i++){
		apop_draw(&out, r, m);
		apop_vector_increment(v, out);
	}
    apop_vector_normalize(d->vector);
    apop_vector_normalize(v);
    for(i=0; i < v->size; i ++)
        Diff(d->vector->data[i], v->data[i], tol2);
}

void test_arms(gsl_rng *r){
    size_t i;
    gsl_vector *o = gsl_vector_alloc(3e5);
    apop_model *ncut = apop_model_set_parameters(apop_normal, 1.1, 1.23);
    ncut->draw = NULL; //testing the default.
    for(i=0; i < 3e5; i ++)
        apop_draw(o->data+i, r, ncut);
    apop_model *back_out = apop_estimate(apop_vector_to_data(o), apop_normal);
    Diff(back_out->parameters->vector->data[0] , 1.1, 1e-2)
    Diff(back_out->parameters->vector->data[1] , 1.23, 1e-2)

    apop_opts.verbose ++;
    apop_model *bcut = apop_model_set_parameters(apop_beta, 0.4, 0.43); //bimodal
    Apop_model_add_group(bcut, apop_arms, .model=bcut, .xl=1e-5, .xr=1-1e-5);
    bcut->draw = NULL; //testing the default.
    for(i=0; i < 3e5; i ++)
        apop_draw((o->data)+i, r, bcut);
    apop_model *back_outb = apop_estimate(apop_vector_to_data(o), apop_beta);
    Diff(back_outb->parameters->vector->data[0] , 0.4, 1e-2)
    Diff(back_outb->parameters->vector->data[1] , 0.43, 1e-2)
    apop_opts.verbose --;
    apop_model *test_copying = apop_model_copy(*back_outb);
    apop_model_free(ncut);
    apop_model_free(back_out);
    apop_model_free(back_outb);
    apop_model_free(test_copying);
}


void test_weighted_regression(apop_data *d,apop_model *e){
    //pretty rudimentary: set all weights to equal and see if we get the same result.
    apop_data *cp = apop_data_copy(d);
    cp->weights = gsl_vector_alloc(d->matrix->size1);
    gsl_vector_set_all(cp->weights, 3);
    apop_model *e2 = apop_estimate(cp, apop_ols);
    assert( apop_vector_distance(e2->parameters->vector, e->parameters->vector) < tol5);
}

void test_ols_offset(gsl_rng *r){
    //A thing we know about OLS: an offset in the variables should make almost no
    //difference. Here we fit OLS with data that is of the form Y = 3*Y+eps
    //and then fit with the same data but offset: Y=3*(Y+20)+eps
    //The coefficient of the constant term moves to absorb the shift;
    //The p-value for that coefficient changes;
    //everything else should remain constant.
    int i, size1 = 1000;
    apop_data *useme = apop_data_alloc(size1, 2);
    for (i =0; i< size1; i++){
        apop_data_set(useme, i, 0, i);
        apop_data_set(useme, i, 1, 3*i+gsl_ran_gaussian(r, 2));
    }
    apop_data *cp = apop_data_copy(useme);
    apop_model *zero_off = apop_estimate(useme, apop_ols);
    Apop_col(cp, 1, off);
    gsl_vector_add_constant(off, 20);
    apop_model *way_off = apop_estimate(cp, apop_ols);
    Apop_col(zero_off->info, 0, zinfo);
    Apop_col(way_off->info, 0, winfo);
    assert(apop_vector_distance(zinfo, winfo) < 1e-4);
    gsl_vector *zcov = apop_data_pack(apop_data_get_page(zero_off->parameters, "<covariance>"));
    gsl_vector *wcov = apop_data_pack(apop_data_get_page( way_off->parameters, "<covariance>"));
    assert(apop_vector_distance(zcov, wcov) < 1e-4);
    assert(apop_data_get(zero_off->parameters, 0, -1) - (apop_data_get(way_off->parameters, 0, -1)+20./3)  < 1e-3);
    assert(apop_data_get(zero_off->parameters, 1, -1) - apop_data_get(way_off->parameters, 1, -1)  < 1e-4);
}

#define do_test(text, fn)   if (verbose)    \
                                printf("%s:", text);  \
                            else printf(".");   \
                            fflush(NULL);   \
                            fn;\
                            {if (verbose) printf(" passed.\n");} 

int main(int argc, char **argv){
  int  slow_tests = 0;
  char c, opts[]  = "sqt:";
    if (argc==1)
        printf("Tests for Apophenia.\nRunning relatively faster tests.  To run slower tests (primarily simulated annealing), use -s.\nFor quieter output, use -q. For multiple threads, use -t2, -t3, ...\nRunning...\n");
    while((c = getopt(argc, argv, opts))!=-1)
        if (c == 's')
            slow_tests  ++;
        else if (c == 'q')
            verbose  --;
        else if (c == 't')
            apop_opts.thread_count  = atoi(optarg);

    //set up some global or common variables
    gsl_rng       *r              = apop_rng_alloc(8); 
    apop_data     *d  = apop_text_to_data("test_data2",0,1);
    apop_model *an_ols_model = apop_model_copy(apop_ols);
    Apop_model_add_group(an_ols_model, apop_lm, .want_cov=1, .want_expected_value= 1);
    apop_model *e  = apop_estimate(d, *an_ols_model);

    do_test("test apop_update", test_updating(r));
    do_test("weighted regression", test_weighted_regression(d,e));
    do_test("offset OLS", test_ols_offset(r));
    do_test("default RNG", test_default_rng(r));
    do_test("log and exponent", log_and_exp(r));
    do_test("test printing", test_printing());
    do_test("test rank expand/compress", rank_round_trip(r));
    do_test("test row set and remove", row_manipulations(r));
    do_test("test apop_map on apop_data_rows", test_apop_map_row());
    do_test("test optimization of multi-page parameters", pack_test());
    do_test("Kullback-Leibler divergence test", test_kl_divergence(r));
    do_test("apop_distance test", test_distances());
    do_test("test column pruning", test_prune_cols());
    do_test("test PMF", test_pmf());
    do_test("apop_pack/unpack test", apop_pack_test(r));
    do_test("test regex", test_regex());
    do_test("test adaptive rejection sampling", test_arms(r));
    do_test("test listwise delete", test_listwise_delete());
    //do_test("test fix params", test_model_fix_parameters(r));
    do_test("positive definiteness", test_posdef(r));
    //do_test("test binomial estimations", test_binomial(r));
    do_test("test data to db", test_data_to_db());
    do_test("test db to crosstab", test_crosstabbing());
    do_test("dummies and factors", dummies_and_factors());
    do_test("test vector/matrix realloc", test_resize());
    do_test("test Fisher exact test", test_fisher());
    do_test("db_to_text", db_to_text());
    do_test("test_vector_moving_average", test_vector_moving_average());
    do_test("apop_estimate->dependent test", test_predicted_and_residual(e));
    do_test("apop_f_test and apop_coefficient_of_determination test", test_f(e));
    do_test("OLS test", test_OLS(r));
    do_test("test lognormal estimations", test_lognormal(r));
    do_test("test queries returning empty tables", test_blank_db_queries());
    do_test("test jackknife covariance", test_jack(r));
    do_test("test apop_histogram model", test_histograms(r));
    do_test("test apop_data sort", test_data_sort());
    do_test("test multivariate_normal", test_multivariate_normal(r));
    do_test("nist_tests", nist_tests());
    do_test("NaN handling", test_nan_data());
    do_test("database skew and kurtosis", test_skew_and_kurt());
    do_test("test_percentiles", test_percentiles());
    do_test("weighted moments", test_weigted_moments());
    do_test("split and stack test", test_split_and_stack());
    do_test("multivariate gamma", test_mvn_gamma());
    do_test("apop_dot test", test_dot());
    do_test("apop_generalized_harmonic test", test_harmonic());
    do_test("apop_strcmp test", test_apop_strcmp());
    do_test("apop_strip_dots test", test_strip_dots());
    do_test("Inversion test", test_inversion(r));
    do_test("apop_jackknife test", test_jackknife());
    do_test("apop_matrix_summarize test", test_summarize());
    do_test("apop_linear_constraint test", test_linear_constraint());
    do_test("transposition test", test_transpose());
    do_test("test unique elements", test_unique_elements());
    if (slow_tests){
        do_test("test distributions", test_distributions(r));
        if (verbose) printf("\tSlower tests:\n");
        do_test("test ML imputation", test_ml_imputation(r));
        do_test("test probit and logit", test_probit_and_logit(r));
        do_test("Test score (dlog likelihood) calculation", test_score());
    }
    printf("\nApophenia has passed all of its tests. Yay.\n");
    return 0;
}
