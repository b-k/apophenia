/*
Here are assorted unit tests, some mechanical and some much more computation-intensive.

The mechanical tests often do things in a convoluted manner for the purpose of touching
as much of the code base as possible. If you want examples for using Apophenia, please 
don't look in the tests---try the eg directory in this distribution or online. 

For an example of the statistical tests, let us say that we wish to verify the results
of a regression. Some systems have a canned screenshot of the 'correct' regression
results that ships with the test suite, and compare a screenshot of the run to the
canned version. I don't get much confidence from this---what if the canned screenshot is
wrong? Better would be to know something about the regression results (like the relation
between the common F statistic, SSR, and SSE) and check that the fact always holds.

Those claims are true as N goes to infinity; for finite N the routines have to strike
a balance. How many draws should I make, and how much user time should I waste, before
measuring the error, and what error tolerance should I set? This is a difficult balance,
and is to some extent the key problem behind all of numeric computing.

There are two types of error bounds here. One is tighter, and therefore more prone
to false alarms, but really forces us to write better numeric code. The other is
much more permissive, and just tells us whether the computation failed to go in the
right direction. Users who run 'make check' will be running the second type of test,
because I (BK) got sick of people sending me bug reports that a test failed
because it reported an error of 1e-5 when it should have been 1e-8. There is always
room for better numeric precision; we all know this without reminders from the
post-install tests.  */

#define _GNU_SOURCE
#include <apop.h>
#include <unistd.h>

#ifdef FULL_TOLERANCE
double tol6 = 1e-6;
double tol5 = 1e-5;
double tol3 = 1e-3;
double tol2 = 1e-2;
double tol1 = 1e-1;
#else
double tol6 = 1e-3;
double tol5 = 1e-3;
double tol3 = 1e-2;
double tol2 = 1e-1;
double tol1 = 1e-1;
#endif

//assertions never return a value.
#undef Apop_assert
#define Apop_assert(expr, ...) {if (!(expr)) {fprintf(stderr, __VA_ARGS__); abort();}}

#define Diff(L, R, eps) {double left=(L), right=(R); Apop_stopif(isnan(left-right) || fabs((left)-(right))>(eps), abort(), 0, "%g is too different from %g (abitrary limit=%g).", (double)(left), (double)(right), eps);}

//A NULL-tolerant strcmp, which used to be a fn and has been deleted.
#define apop_strcmp(a, b) (((a)&&(b) && !strcmp((a), (b))) || (!(a) && !(b)))

int len = 8000;
int verbose = 1;

void test_nan_data();
void db_to_text();

void apop_data_scale (apop_data *d, double scale){
    if (d->vector) gsl_vector_scale(d->vector, scale);
    if (d->matrix) gsl_matrix_scale(d->matrix, scale);
}


void v_pow10(double *in){ *in = pow(10,*in);}
double log_for_map(gsl_vector *v){apop_vector_log(v); return apop_sum(v);}
double log_by_val(double x){return x;}

static void log_and_exp(gsl_rng *r){
    apop_data *d = apop_data_alloc(100,2);
    apop_data_add_names(d, 'c', "10", "e");
    for (int i=0; i< 100; i++){
        apop_data_set(d, i, .colname="10", .val=gsl_rng_uniform(r)*10);
        double *eval = apop_data_ptr(d, .row=i, .colname="e");
        *eval = gsl_rng_uniform(r)*10;
    }
    apop_data *d2 = apop_data_copy(d);
    Apop_col_tv(d, "10", tencol);
    apop_vector_log10(tencol);
    apop_vector_apply(tencol, v_pow10);
    Apop_col_v(d2, 0, o_tencol);
    assert(apop_vector_distance(tencol, o_tencol) < 1e-3);
    
    Apop_col_tv(d, "e", ecol);
    apop_vector_log(ecol);
    apop_vector_exp(ecol);
    Apop_col_v(d2, 1, o_ecol);
    assert(apop_vector_distance(ecol, o_ecol) < 1e-3);

    apop_data *d5 = apop_data_alloc(5,5);
    for (int i=0; i< 5; i++)
        for (int j=0; j< 5; j++)
            apop_data_set(d5, i, j, gsl_rng_uniform(r));
    double log_sum_by_v = apop_matrix_map_sum(d5->matrix, log_for_map);
    //d5 is now all logs.
    double log_sum = apop_matrix_map_all_sum(d5->matrix, log_by_val);
    Diff(log_sum_by_v, log_sum, 1e-5);
    Diff(apop_matrix_mean(d5->matrix)*25., log_sum, 1e-5);
    apop_data_free(d); apop_data_free(d2);
    apop_data_free(d5);
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
    apop_model_draws(mvn, .draws=fillme);
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
    apop_data_free(fillme);
}

void test_percentiles(){
    gsl_vector *v = gsl_vector_alloc(307);
    for (size_t i=0; i< 307; i++)
        gsl_vector_set(v, i, i);
    double *pcts_up    = apop_vector_percentiles(v, 'u');
    double *pcts_down  = apop_vector_percentiles(v, 'd');
    double *pcts_avg   = apop_vector_percentiles(v, 'a');
    for (size_t i=0; i< 101; i++){
        assert(pcts_up[i] >= pcts_down[i]);
        assert(pcts_up[i] >= pcts_avg[i]);
        assert(pcts_avg[i] >= pcts_down[i]);
    }
    assert(pcts_up[100] == pcts_down[100] && pcts_avg[100] == pcts_down[100]);
    assert(pcts_up[0] == pcts_down[0] && pcts_avg[0] == pcts_down[0]);
    assert(pcts_avg[50] == (pcts_down[50] + pcts_up[50])/2);
}

void test_score(){
    int len = 1e5;
    gsl_rng *r = apop_rng_alloc(123);
    apop_data *data = apop_data_alloc(len,1);
    apop_model *source = apop_model_set_parameters(apop_normal, 
                                gsl_ran_flat(r, -5, 5), gsl_ran_flat(r, .01, 5));
    for (size_t j=0; j< len; j++)
        apop_draw(gsl_matrix_ptr(data->matrix, j, 0), r, source);
    apop_model *estme = apop_model_copy(apop_normal);
    Apop_model_add_group(estme, apop_mle, .method= "annealing");
    apop_prep(data, estme);
    apop_maximum_likelihood(data, estme);

    apop_model *straight_est = apop_estimate(data, apop_normal);
    Diff (straight_est->parameters->vector->data[0], source->parameters->vector->data[0], tol1);
    Diff (straight_est->parameters->vector->data[1], source->parameters->vector->data[1], tol1);
    apop_model_free(straight_est); 

    double sigsqn = gsl_pow_2(estme->parameters->vector->data[1])/len;
    apop_data *cov = apop_data_get_page(estme->parameters, "cov", 'r');
    Diff (apop_data_get(cov, 0, 0), sigsqn , tol3);
    Diff (apop_data_get(cov, 1, 1), sigsqn/2 , tol3);
    double *cov1 = apop_data_ptr(estme->parameters, .page="<covariance>", .row=1, .col=1);
    Diff (*cov1 ,sigsqn/2 , tol3);
    Diff(apop_data_get(cov, 0,1) + apop_data_get(cov, 0,1), 0, tol3);
    apop_model_free(estme);
    apop_model_free(source); 
    apop_data_free(data);
}

void test_normalizations(gsl_vector *v){
    //let's check out normalizations, while we have a vector for it
    gsl_vector_scale(v, 23);
    gsl_vector_add_constant(v, 8);
    apop_data dv = (apop_data){.matrix=apop_vector_to_matrix(v)};
    apop_data_transpose(&dv);
    apop_matrix_normalize(dv.matrix, 'r', 's');
    apop_data *dvagain = apop_data_transpose(&dv, .inplace='n');
    apop_data *sum = apop_data_summarize(dvagain);
    apop_data_free(dvagain);
    Diff(apop_data_get(sum, .colname="mean"), 0, 1e-5);
    Diff(apop_data_get(sum, .colname="std dev"), 1, 1e-5);
    Diff(apop_data_get(sum, .colname="variance"), 1, 1e-5);
    apop_data_free(sum);
    gsl_matrix_free(dv.matrix);
}

//This tests the database-side functions. Only for sqlite.
void test_skew_and_kurt(gsl_rng *r){
    gsl_vector *v;
    if (apop_opts.db_engine=='s'){
        apop_table_exists(.remove='d', .name="t");
        apop_query("create table t(vals)");
        for(int i=0; i<1e4; i++)
            apop_query("insert into t values(%g)", gsl_rng_uniform(r));
        v = apop_query_to_vector("select * from t");
        Diff (apop_var(v), apop_query_to_float("select var(vals) from t"), tol6);
        Diff (apop_vector_skew(v), apop_query_to_float("select skew(vals) from t"), tol6);
        Diff (apop_vector_kurtosis(v), apop_query_to_float("select kurt(vals) from t"), tol5);
        apop_table_exists("t", 'd');
    } else {
        v = gsl_vector_alloc(1e4);
        for(int i=0; i<1e4; i++) gsl_vector_set(v, i,  gsl_rng_uniform(r));
    }
    test_normalizations(v);
    gsl_vector_free(v);
}

void test_listwise_delete(){
  apop_data *t1 = apop_data_calloc(10,10);
  apop_text_alloc(t1, 10, 10);
  for (int i=0; i< 10; i++)
      for (int j=0; j< 10; j++)
          apop_text_add(t1, i, j, "%i", i*j);
  //no NaNs yet
  apop_data *t1c = apop_data_listwise_delete(t1);
  assert(t1c->matrix->size1==10);
  assert(t1c->matrix->size2==10);
  assert(atoi(t1c->text[3][4])==12);
  //Now kill two rows
  asprintf(&(t1c->text[3][4]), "nan");
  asprintf(&(t1c->text[9][9]), "NaN");
  t1c = apop_data_listwise_delete(t1c);
  assert(t1c->matrix->size1==8);
  assert(t1c->matrix->size2==10);
  assert(atoi(t1c->text[3][4])==16);
  //check the vector
  apop_data *t2 = apop_data_calloc(10); //check on this form of calloc.
  t1->vector = t2->vector;
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
  //return NULL if every row has missing data.
  Apop_col_v(t1, 7, v)
  gsl_vector_set_all(v, GSL_NAN);
  assert(!apop_data_listwise_delete(t1));

    //btw, check non-square transpose with blanks.
    t1 = apop_data_calloc();
    apop_text_alloc(t1, 12, 10);
    for (int i=0; i< 10; i++)
        for (int j=0; j< 9; j++)
            apop_text_add(t1, i, j, "%i", i*j);
    t2 = apop_data_copy(t1);
    apop_data_transpose(t1);
    assert(!strlen(t1->text[7][11]));

    //come back
    apop_data_transpose(t1);
    assert(!strlen(t1->text[11][7]));
    for (int i=0; i< 10; i++)
        for (int j=0; j< 9; j++)
            assert(!strcmp(t1->text[i][j], t2->text[i][j]));

    apop_data_free(t1);
    t1 = apop_text_alloc(NULL, 10, 12);
    for (int i=0; i< 9; i++)
        for (int j=0; j< 10; j++)
            apop_text_add(t1, i, j, "%i", i*j);
    apop_data *t4 = apop_data_transpose(t1, .inplace='n');
    assert(!strlen(t4->text[11][7]));
    assert(atoi(t4->text[9][8])==72);
    assert(t4->textsize[0]==12);
    assert(t4->textsize[1]==10);
}


static void wmt(gsl_vector *v, gsl_vector *v2, gsl_vector *w, gsl_vector *av, gsl_vector *av2, double mean){
    assert(apop_vector_mean(av) == apop_vector_mean(v,w));
    assert(apop_vector_mean(v, w) == mean);
    Diff (apop_vector_var(av), apop_vector_var(v, w), tol5);
    Diff (apop_vector_cov(av, av2), apop_vector_cov(v, v2, w), tol5);
    Diff (apop_vector_skew_pop(av), apop_vector_skew_pop(v, w), tol5);
    Diff (apop_vector_kurtosis_pop(av), apop_vector_kurtosis_pop(v, w), tol5);
}

void test_weigted_moments(){
  //double        data[]      = {1,2,3};//checking vector_fill with this and w2
  double        alldata[]   = {1,2,3};
  double        data3[]     = {3,2,1};
  double        alldata3[]  = {3,2,1};
  double        weights[]   = {1,1,1};
  gsl_vector    *v          = apop_vector_fill(gsl_vector_alloc(3), 1, 2, 3);
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
    gsl_vector *w2 = gsl_vector_alloc(5);
    apop_vector_fill(w2, 4, 3, 2, 1, 0);
    wmt(v, v2, w2, av, av2, 1);
}

void test_split_and_stack(gsl_rng *r){
    apop_data *d1 = apop_data_alloc(10,10,10);
    int i,j, tr, tc;
    apop_data **splits, *dv2;
    for(i=-1; i< 10; i++)
        for(j=0; j< 10; j++)
            apop_data_set(d1, j, i, gsl_rng_uniform(r));

    d1->weights = gsl_vector_alloc(10);
    for(j=0; j< 10; j++)
        gsl_vector_set(d1->weights, j, gsl_rng_uniform(r));

    //vector_stacking NULLs:
    gsl_vector *orig=apop_vector_copy(d1->vector);
    gsl_vector *cp=apop_vector_stack(NULL, d1->vector);
    apop_vector_stack(d1->vector, NULL);
    assert(d1->vector->size==10);
        for(j=0; j< 10; j++){
            assert(gsl_vector_get(d1->vector, j)==gsl_vector_get(orig, j));
            assert(gsl_vector_get(cp, j)==gsl_vector_get(orig, j));
        }
    assert(!apop_vector_stack(NULL, NULL));
    assert(!apop_data_stack(NULL, NULL));
    gsl_vector_free(orig);
    gsl_vector_free(cp);

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
            assert(gsl_vector_get(dv2->weights, tr) == gsl_vector_get(d1->weights,tr));
            if(tr < i){
                assert(apop_data_get(splits[0], tr, tc)==apop_data_get(d1, tr, tc));
                assert(apop_data_get(splits[0], tr, tc)==apop_data_get(d1, tr, tc));
                assert(gsl_vector_get(splits[0]->weights, tr) == gsl_vector_get(d1->weights,tr));
            } else{
                int start=i < 0 ? 0 : i;
                assert(apop_data_get(splits[1], tr-start, tc)==apop_data_get(d1, tr, tc));
                assert(apop_data_get(splits[1], tr-start, tc)==apop_data_get(d1, tr, tc));
                assert(gsl_vector_get(splits[1]->weights, tr-start) == gsl_vector_get(d1->weights,tr));
            }
        }
        apop_data_free(dv2);
        apop_data_free(splits[0]);
        apop_data_free(splits[1]);
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
            assert(gsl_vector_get(dv2->weights, tr) == gsl_vector_get(d1->weights,tr));
            if (splits[0]) assert(gsl_vector_get(splits[0]->weights, tr) == gsl_vector_get(d1->weights,tr));
            if (splits[1]) assert(gsl_vector_get(splits[1]->weights, tr) == gsl_vector_get(d1->weights,tr));
            if(tc < i){
                assert(apop_data_get(splits[0], tr, tc)==apop_data_get(d1, tr, tc));
                assert(apop_data_get(splits[0], tr, tc)==apop_data_get(d1, tr, tc));
            } else{
                int start=i < 0 ? 0 : i;
                assert(apop_data_get(splits[1], tr, tc-start)==apop_data_get(d1, tr, tc));
                assert(apop_data_get(splits[1], tr, tc-start)==apop_data_get(d1, tr, tc));
            }
        }
        apop_data_free(splits[0]);
        apop_data_free(splits[1]);
        apop_data_free(dv2);
    }

    apop_data **notsplits = apop_data_split(d1, 800, 'c');
    assert(notsplits[0]->matrix->size1 == d1->matrix->size1);
    assert(notsplits[0]->matrix->size2 == d1->matrix->size2);
    assert(notsplits[0]->vector->size == d1->vector->size);
    assert(notsplits[1] == NULL);
    //let's try a NULL data set 
    apop_data *nulldata = apop_data_alloc();
    apop_data *onespot = apop_data_alloc(1,1);
    apop_data_set(onespot, 0, 0, 12);
    apop_data_stack(nulldata, onespot, .posn='c', .inplace='y');
    assert(nulldata->matrix->size1==1);
    assert(nulldata->matrix->size2==1);
    assert(apop_data_get(nulldata, 0, 0) == 12);
    apop_data_free(nulldata);
    apop_data_free(onespot);

    //text
    apop_data *txt = apop_text_alloc(NULL, 3,3);
    apop_data *txt2 = apop_text_alloc(NULL, 3,3);
    for (int i=0; i< 3; i++)
        for (int j=0; j< 3; j++){
            apop_text_add(txt, i, j, "(%i, %i)", i, j);
            apop_text_add(txt2, i, j, "[%i, %i]", i, j);
        }

    apop_data *rbound = apop_data_stack(txt, txt2, .posn='r');
    assert(rbound->textsize[0] == 6);
    assert(rbound->textsize[1] == 3);
    assert(!strcmp(rbound->text[2][2], "(2, 2)"));
    assert(!strcmp(rbound->text[5][2], "[2, 2]"));

    apop_data_stack(txt, txt2, .posn='c', .inplace='y');
    assert(txt->textsize[0] == 3);
    assert(txt->textsize[1] == 6);
    assert(!strcmp(txt->text[2][2], "(2, 2)"));
    assert(!strcmp(txt->text[2][5], "[2, 2]"));
    apop_data_free(txt);
    apop_data_free(txt2);

    apop_data **txtsplits = apop_data_split(rbound, 3, 'r');
    apop_data_free(rbound);
    assert(txtsplits[0]->textsize[0] ==txtsplits[0]->textsize[1]);
    assert(txtsplits[0]->textsize[0] == 3);
    assert(txtsplits[1]->textsize[0] ==txtsplits[1]->textsize[1]);
    assert(txtsplits[1]->textsize[0] == 3);
    assert(!strcmp(txtsplits[0]->text[1][1], "(1, 1)"));
    assert(!strcmp(txtsplits[1]->text[1][1], "[1, 1]"));
    apop_data_free(d1);
    apop_data_free(onespot);
    apop_data_free(nulldata);
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

    apop_data *predict_tab = apop_data_get_page(est->info, "predict", 'r');
    v   = gsl_matrix_column(predict_tab->matrix, apop_name_find(predict_tab->names, "residual", 'c')).vector;
    assert(fabs(apop_mean(&v)) < tol5);

    Apop_col_tv(predict_tab, "predicted", vv);
    gsl_blas_dgemv(CblasNoTrans, 1, m, est->parameters->vector, 0, prediction);
    gsl_vector_sub(prediction, vv);
    assert(fabs(apop_vector_sum(prediction)) < tol5);
}

/** I claim that the F test calculated via apop_F_test(est, NULL, NULL)
 equals a transformation of R^2 (after a normalization step).
*/
void test_f(apop_model *est){
    apop_data *rsq  = apop_estimate_coefficient_of_determination(est);
    apop_data *constr= apop_data_calloc(est->parameters->vector->size-1, est->parameters->vector->size);
    int i;
    for (i=1; i< est->parameters->vector->size; i++)
        apop_data_set(constr, i-1, i, 1);
    apop_data *ftab = apop_F_test(est, constr);
    apop_data *ftab2 = apop_F_test(est, NULL);
    //apop_data_show(ftab);
    //apop_data_show(ftab2);
    double n = est->data->matrix->size1;
    double K = est->parameters->vector->size-1;
    double r = apop_data_get(rsq, .rowname="R squared");
    double f = apop_data_get(ftab, .rowname="F statistic");
    double f2 = apop_data_get(ftab2, .rowname="F statistic");
    Diff (f , r*(n-K)/((1-r)*K) , tol5);
    Diff (f2 , r*(n-K)/((1-r)*K) , tol5);
}

void test_OLS(gsl_rng *r){
    apop_data *set = apop_data_alloc(0, len, 2);
    for(int i=0; i< len; i++){
        apop_data_set(set, i, 1, 100*(gsl_rng_uniform(r)-0.5));
        apop_data_set(set, i, 0, -1.4 + apop_data_get(set,i,1)*2.3);
    }
    apop_data *bkup = apop_data_copy(set);
    apop_model *out = apop_estimate(set, apop_ols);
    Diff (apop_data_get(out->parameters, 0,-1) , -1.4 , tol5);
    Diff (apop_data_get(out->parameters, 1,-1) , 2.3 , tol5);
    apop_model_free(out);

    gsl_vector *w = gsl_vector_alloc(set->matrix->size1);
    gsl_vector_set_all(w, 14);
    bkup->weights  = w;
    out = apop_estimate(bkup, apop_ols);
    Diff (apop_data_get(out->parameters, 0,-1) , -1.4 , tol5);
    Diff (apop_data_get(out->parameters, 1,-1) , 2.3 , tol5);
    apop_data_free(bkup);
    apop_data_free(set);
    apop_model_free(out);
}

#define INVERTSIZE 100
void test_inversion(gsl_rng *r){
    gsl_matrix *invme = gsl_matrix_alloc(INVERTSIZE, INVERTSIZE);
    gsl_matrix *inved;
    gsl_matrix *inved_back;
    apop_data *four = apop_data_alloc(1);
    apop_data_set(four, 0, -1, 4);
    apop_model *fourp = apop_model_copy(apop_zipf);
    fourp->parameters = four;
    for(int i=0; i<INVERTSIZE; i++)
        for(int j=0; j<INVERTSIZE; j++)
            apop_zipf->draw(gsl_matrix_ptr(invme, i,j),r,  fourp);
    apop_det_and_inv(invme, &inved, 0, 1);
    apop_det_and_inv(inved, &inved_back, 0, 1);
    apop_model_free(fourp);
    double error = 0;
    for(int i=0; i<INVERTSIZE; i++)
        for(int j=0; j<INVERTSIZE; j++)
            error += gsl_matrix_get(invme, i,j) - gsl_matrix_get(inved_back, i,j);
    assert (error < 1e-5);
    gsl_matrix_free(invme);
    gsl_matrix_free(inved);
    gsl_matrix_free(inved_back);
}

void test_summarize(){
    apop_table_exists("td", 'd');
    apop_text_to_db("test_data", .has_row_names= 0,1, .tabname = "td");
    gsl_matrix *m = apop_query_to_matrix("select * from td");
    apop_data *s = apop_data_summarize(apop_matrix_to_data(m));
    gsl_matrix_free(m);
    double t = gsl_matrix_get(s->matrix, 1,0);
    assert (t ==3);
    t = gsl_matrix_get(s->matrix, 2, 1);
    double v = sqrt((2*2 +3*3 +3*3 +4.*4.)/3.);
    assert (t == v);
    apop_data_free(s);
}

void test_dot(){
apop_data *d1   = apop_text_to_data(.text_file="test_data2",0,1); // 55 x 2
apop_data *d2   = apop_text_to_data("test_data2"); // 55 x 2
apop_data *d3   = apop_dot(d1, d2, .form2='t');
gsl_vector v1 = gsl_matrix_row(d1->matrix, 0).vector;
apop_data *dv   = apop_vector_to_data(&v1); // 2 x 1
apop_data *d7   = apop_dot(dv, dv);
    assert(apop_data_get(d7, 0, -1) == apop_data_get(d3, 0,0));
apop_data *d8   = apop_dot(d1, dv);
apop_data *d9   = apop_dot(dv, d1, .form2='t');
    gsl_vector_sub(d8->vector, d9->vector);
    assert(!apop_vector_sum(d8->vector));

    int verbosity = apop_opts.verbose;
    apop_opts.verbose = -1;
    apop_data *d10 = apop_dot(d1, d2);
    assert(d10->error == 'd');
    apop_data *d11 = apop_dot(dv, d1);
    assert(d11->error == 'd');
    apop_opts.verbose = verbosity;
    apop_data_free(d1); apop_data_free(d2); apop_data_free(d3); 
    dv->vector=NULL; apop_data_free(dv); 
    apop_data_free(d7); apop_data_free(d8);
    apop_data_free(d9); apop_data_free(d10); apop_data_free(d11);
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
    size_t ct = 1000;
    apop_data *d = apop_data_alloc(0,ct,2);
    double draw[2];
    apop_multivariate_normal->vsize =
    apop_multivariate_normal->msize1 =
    apop_multivariate_normal->msize2 = 2;
    apop_model *pp = apop_model_set_parameters(apop_multivariate_normal,
                                        8, 1, 0.5,
                                        2, 0.5, 1);
    for(int i=0; i< ct; i++){
        apop_multivariate_normal->draw(draw, r, pp);
        apop_data_set(d, i, 0, draw[0]);
        apop_data_set(d, i, 1, draw[1]);
    }

    apop_data *pcopy = apop_data_copy(pp->parameters);
    gsl_matrix_set_all(pp->parameters->matrix, GSL_NAN);
    apop_model *mep1  = apop_model_fix_params(pp);
    Apop_settings_add(mep1, apop_mle, starting_pt, ((double[]){1.5, .25, .25, 1.5}));
    apop_model *e1 = apop_estimate(d, mep1);
    gsl_vector_sub(e1->parameters->vector, pcopy->vector);
    assert(apop_vector_sum(e1->parameters->vector) < 1e-1);
    apop_model_free(e1);

    double start2[] = {7,3};
    pp->parameters = apop_data_copy(pcopy);
    gsl_vector_set_all(pp->parameters->vector, GSL_NAN);
    apop_model *mep2  = apop_model_fix_params(pp);
    apop_model_free(pp);
    Apop_settings_add(mep2, apop_mle, starting_pt, start2);
    Apop_settings_add(mep2, apop_mle, method, "PR cg");

    apop_model *e2 = apop_estimate(d, mep2);
    apop_model_free(mep2);
    gsl_matrix_sub(e2->parameters->matrix, pcopy->matrix);
    assert(apop_matrix_sum(e2->parameters->matrix) < 1e-2);
    apop_model_free(e2);
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

static void broken_est(apop_data *d, apop_model *m){
    static gsl_rng *r; if (!r) r = apop_rng_alloc(1);
    if (gsl_rng_uniform(r) < 1./100.) {
        gsl_vector_set_all(m->parameters->vector, GSL_NAN);
        return;
    }
    apop_normal->estimate(d, m);
}

static void super_broken_est(apop_data *d, apop_model *m){
    static gsl_rng *r; if (!r) r = apop_rng_alloc(1);
    if (gsl_rng_uniform(r) < 3./4.) {
        gsl_vector_set_all(m->parameters->vector, GSL_NAN);
        return;
    }
    apop_normal->estimate(d, m);
}

void test_jackknife(gsl_rng *r){
    double pv[] = {3.09,2.8762};
    int len = 2000;
    apop_model *m = apop_model_copy(apop_normal);
    apop_data *d = apop_data_alloc(0, len, 1);
    apop_data *p = apop_data_alloc();
    p->vector = apop_array_to_vector(pv, 2);
    apop_model*pp = apop_model_copy(m);
    pp->parameters = p;
    for (size_t i =0; i< len; i++)
        m->draw(apop_data_ptr(d, i, 0), r, pp); 
    apop_data *out = apop_jackknife_cov(d, m);
    //Notice that the jackknife just ain't a great estimator here.
assert ((fabs(apop_data_get(out, 0,0) - gsl_pow_2(pv[1])/len)) < tol2 
            && fabs(apop_data_get(out, 1,1) - gsl_pow_2(pv[1])/(2*len)) < tol2*100);
    apop_data *out2 = apop_bootstrap_cov(d, m, .keep_boots='y');
    assert (fabs(apop_data_get(out2) - gsl_pow_2(pv[1])/len) < tol2
                && fabs(apop_data_get(out2, 1,1) - gsl_pow_2(pv[1])/(2*len)) < tol2);
    apop_data_free(out2);

    //bootstrap should recover gracefully from a small number of NaNs...
    m->estimate = broken_est;
    out2 = apop_bootstrap_cov(d, m, .ignore_nans='y');
    assert (fabs(apop_data_get(out2) - gsl_pow_2(pv[1])/len) < tol2
                && fabs(apop_data_get(out2, 1,1) - gsl_pow_2(pv[1])/(2*len)) < tol2);


    //...but not from a large number of NaNs.
    int vvv= apop_opts.verbose;
    apop_opts.verbose = -1;
    m->estimate= super_broken_est;
    apop_data *out3 = apop_bootstrap_cov(d, m, .ignore_nans='y');
    assert(out3->error);
    apop_opts.verbose = vvv;

    apop_data_free(d);
    apop_data_free(out);
    apop_data_free(out2);
    apop_data_free(out3);
    apop_model_free(m);
}

//In my inattention, I wrote two jackknife tests. So you get double the checks.
int test_jack(gsl_rng *r){
  int i, draws     = 1000;
  apop_data *d = apop_data_alloc(draws, 1);
  apop_model *m = apop_normal;
  double pv[] = {1., 3.};
    m->parameters = apop_data_fill_base(apop_data_alloc(2), pv);
    for (i =0; i< draws; i++)
        m->draw(apop_data_ptr(d, i, 0), r, m); 
    apop_data *out = apop_jackknife_cov(d, m);
    double error = fabs(apop_data_get(out, 0,0)-gsl_pow_2(pv[1])/draws) //var(mu)
                + fabs(apop_data_get(out, 1,1)-gsl_pow_2(pv[1])/(2*draws))//var(sigma)
                +fabs(apop_data_get(out, 0,1)) +fabs(apop_data_get(out, 1,0));//cov(mu,sigma); should be 0.
    apop_data_free(d);
    apop_data_free(out);
    return (error < 1e-2);//still not very accurate.
}

void test_lognormal(gsl_rng *r){
    apop_model *source = apop_model_copy(apop_normal);
    apop_model_clear(NULL, source);
    double mu = gsl_ran_flat(r, -1, 1);
    double sigma = gsl_ran_flat(r, .01, 1);
    int n = gsl_ran_flat(r,1,8e5);
    apop_data *data = apop_data_alloc(0,1,n);
    gsl_vector_set(source->parameters->vector, 0, mu);
    gsl_vector_set(source->parameters->vector, 1, sigma);
    for (int j=0; j< n; j++){
        double *k   = gsl_matrix_ptr(data->matrix, 0, j);
        apop_draw(k, r, source);
        *k = exp(*k);
    }
    apop_model_free(source);
    apop_model *out = apop_estimate(data, apop_lognormal);
    double muhat = apop_data_get(out->parameters, 0,-1);
    double sigmahat = apop_data_get(out->parameters, 1,-1);
    //if (verbose) printf("mu: %g, muhat: %g, var: %g, varhat: %g\n", mu, muhat,  sigma,sigmahat);
    Diff(mu, muhat, 1e-2);
    Diff(sigma, sigmahat, 1e-2);
    apop_model_free(out);

    apop_model *for_mle= apop_model_copy(apop_lognormal);
    for_mle->estimate=NULL;
    apop_model *out2 = apop_estimate(data, for_mle);
    apop_model_free(for_mle);
    muhat = apop_data_get(out2->parameters, 0,-1);
    sigmahat = apop_data_get(out2->parameters, 1,-1);
    Diff(mu, muhat, 1e-2);
    Diff(sigma, sigmahat, 1e-2);
    apop_model_free(out2);
    apop_data_free(data);
}

void test_multivariate_normal(){
    int len = 5e5;
    double params[] = {1, 3, 0,
                       2, 0, 1};
    apop_data *p = apop_data_fill_base(apop_data_alloc(2, 2, 2), params);
    apop_model *mv = apop_model_copy(apop_multivariate_normal);
    mv->parameters=p;
    mv->dsize=2;
    apop_data *rdraws = apop_model_draws(mv, .count=len);
    mv->parameters=NULL;
    apop_model_free(mv);
    apop_model *est =apop_estimate(rdraws, apop_multivariate_normal);
    double error = fabs(est->parameters->vector->data[0] - p->vector->data[0])
                  +fabs(est->parameters->vector->data[1] - p->vector->data[1])
                  +fabs(est->parameters->matrix->data[0] - p->matrix->data[0])
                  +fabs(est->parameters->matrix->data[1] - p->matrix->data[1])
                  +fabs(est->parameters->matrix->data[2] - p->matrix->data[2])
                  +fabs(est->parameters->matrix->data[3] - p->matrix->data[3]);
    Diff(error, 0, 4e-2); //yes, unimpressive, but we don't wanna be here all day.
    apop_model_free(est);
    apop_data_free(rdraws);
}

static void common_binomial_bit(apop_model *out, int n, double p){
    /*double phat = apop_data_get(out->parameters, 1,-1);
    double nhat = apop_data_get(out->parameters, 0,-1);
    if (verbose) printf("n: %i, p: %g, nhat: %g, phat: %g\n", n, p, phat, nhat);*/
    assert(apop_data_get(out->parameters, 0,-1) == n);
    assert(apop_data_get(out->parameters, 1,-1) - p < 1e-2);
}

void test_binomial(gsl_rng *r){
    double p = gsl_rng_uniform(r);
    int n     = gsl_ran_flat(r,1,1e5);
    apop_data *d = apop_data_falloc((1,2), n*(1-p), n*p);
    apop_model *out = apop_estimate(d, apop_binomial);
    apop_model *outm = apop_estimate(d, apop_multinomial);
    common_binomial_bit(out, n, p);
    common_binomial_bit(outm, n, p);
    apop_model_free(outm);
    apop_data_free(d);
    apop_model_free(out);
}

void test_rownames(){
    apop_data *d = apop_data_falloc((2, 2), 0, 1, 
                                            2, 3);
    apop_data_add_names(d, 'r', "zero", "one");
    apop_data_add_names(d, 'c', "C zero", "C one");
    assert(apop_data_get(d, .rowname="zero", .col=0) == 0);
    assert(apop_data_get(d, .rowname="zero", .colname="C zero") == 0);
    assert(apop_data_get(d, .rowname="one", .col=0) == 2);

    double *oneone = apop_data_ptr(d, .rowname="one", .col=1);
    *oneone= 27;
    assert(apop_data_get(d, .rowname="one", .colname="C one")==27);

    apop_data_set(d, .rowname="one", .colname="C zero", .val=33);
    double *onezero = apop_data_ptr(d, .rowname="one", .col=0);
    assert(*onezero == 33);

    apop_data_set(d, .rowname="zero", .col=1, .val=10);
    double *zeroone = apop_data_ptr(d, .rowname="zero", .colname="C one");
    assert(*zeroone == 10);
}

int get_factor_index(apop_data *flist, char *findme){
    for (int i=0; i< flist->textsize[0]; i++)
        if (apop_strcmp(flist->text[i][0], findme))
            return i;
    return -2;
}

//If the dummies are a separate matrix, offset=0;
//If the dummies are an addendum to main, offset=original_data->matrix->size2;
static void check_for_dummies(apop_data *d, apop_data *dum, int offset){
  int n;
    apop_data *factorlist = apop_data_get_factor_names(d, 0, 't');
    for(int i=0; i < d->textsize[0]; i ++)
        if ((n = get_factor_index(factorlist, d->text[i][0]))>0){
            for(int j=0; j < factorlist->textsize[0]-1; j ++)
                if (j==n-1)
                    assert(apop_data_get(dum, i, j+offset));
                else
                    assert(!apop_data_get(dum, i, j+offset));
        } else
            for(int j=0; j < factorlist->textsize[0]-1; j ++)
                assert(!apop_data_get(dum, i, j+offset));
}

void dummies_and_factors(){
    apop_text_to_db("data-mixed", "genes");
    apop_data *d = apop_query_to_mixed_data("mmmt", "select aa, bb, 1, a_allele from genes");
    apop_data *dum = apop_data_to_dummies(d, 0, 't', 0);
    check_for_dummies(d, dum, 0);
    apop_data_to_factors(d, 't', 0, 2);
    for(int i=0; i < d->textsize[0]; i ++) //the set is only As and Cs.
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
    apop_data *tt = apop_data_transpose(t, .inplace='n');
    assert(apop_data_get(tt, 0, 3) == 9);
    assert(apop_data_get(tt, 1, 0) == 4);
    assert(!strcmp(tt->names->row[2], "c"));
    assert(!strcmp(tt->names->row[3], "d"));
    assert(!tt->names->colct);
    apop_data_transpose(t);
    assert(apop_data_get(t, 0, 3) == 9);
    assert(apop_data_get(t, 1, 0) == 4);
    assert(!strcmp(t->names->row[2], "c"));
    assert(!strcmp(t->names->row[3], "d"));
    assert(!t->names->colct);
}

apop_data *generate_probit_logit_sample (gsl_vector* true_params, gsl_rng *r, apop_model *method){
  int i, j;
  double val;
  int samples = 8e4;
  apop_data *data = apop_data_alloc(samples, true_params->size);
        //generate a random vector of data X, then set the outcome to one if the probit/logit condition holds
        for (i = 0; i < samples; i++){
            apop_data_set(data, i, 0, 1);
            for (j = 1; j < true_params->size; j++)
                apop_data_set(data, i, j, (gsl_rng_uniform(r)-0.5) *2);
            Apop_row_v(data, i, asample);
            gsl_blas_ddot(asample, true_params, &val);
            if (method == apop_probit)
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
    int param_ct = gsl_rng_uniform(r)*7 + 1; //up to seven params.
    gsl_vector *true_params = gsl_vector_alloc(param_ct);
    for (int j = 0; j < param_ct; j++)
        gsl_vector_set(true_params, j, (gsl_rng_uniform(r)-0.5)*2);

    //Logit
    apop_data* data = generate_probit_logit_sample(true_params, r, apop_logit);
    Apop_model_add_group(apop_logit, apop_mle, .tolerance=1e-5);
    Apop_model_add_group(apop_logit, apop_parts_wanted);
    apop_model *m = apop_estimate(data, apop_logit);
    Apop_col_v(m->parameters, 0, logit_params);
    assert(apop_vector_distance(logit_params, true_params) < 0.07);
    apop_data_free(data);
    apop_model_free(m);

    //Probit
    apop_data* data2 = generate_probit_logit_sample(true_params, r, apop_probit);
    Apop_model_add_group(apop_probit, apop_mle);
    Apop_model_add_group(apop_logit, apop_parts_wanted);
    m = apop_estimate(data2, apop_probit);
    Apop_col_v(m->parameters, 0, probit_params);
    assert(apop_vector_distance(probit_params, true_params) < 0.07);
    gsl_vector_free(true_params);
    apop_model_free(m);
    apop_data_free(data2);
}

void test_resize(){
    //This is the multiplication table from _Modeling with Data_
    //with a +.1 to distinguish columns from rows.
    gsl_matrix *m = apop_matrix_realloc(NULL, 20,15);//check using realloc as an alloc
    gsl_matrix_set_all(m, 1);
    for (int i=0; i< m->size1; i++){
        Apop_matrix_row(m, i, one_row);
        gsl_vector_scale(one_row, i+1);
    }
    for (int i=0; i< m->size2; i++){
        Apop_matrix_col(m, i, one_col);
        gsl_vector_scale(one_col, i+1);
        gsl_vector_add_constant(one_col, (i+1)/10.);
    }
    apop_matrix_realloc(m, 11, 17);
    assert(gsl_matrix_get(m, 3, 5) == 4*6+.6);
    apop_matrix_realloc(m, 10, 10);
    Diff (apop_matrix_sum(m) , 55 * 56 , tol6);
    gsl_vector *v = gsl_vector_alloc(20);
    for (int i=0; i< 20; i++)
        gsl_vector_set(v, i, i);
    apop_vector_realloc(v, 38);
    for (int i=0; i< 20; i++)
        assert(gsl_vector_get(v, i) == i);
    apop_vector_realloc(v, 10);
    assert(apop_vector_sum(v) == 45);
}

void test_mvn_gamma(){
    assert(apop_multivariate_gamma(10, 1)==gsl_sf_gamma(10));
    assert(apop_multivariate_lngamma(10, 1)==gsl_sf_lngamma(10));
}

void test_default_rng(gsl_rng *r) {
    gsl_vector *o = gsl_vector_alloc(2e5);
    apop_model *ncut = apop_model_set_parameters(apop_normal, 1.1, 1.23);
    ncut->draw = NULL; //forced to use the default.
    for(size_t i=0; i < 2e5; i ++)
        apop_draw(o->data+i, r, ncut);
    apop_model *back_out = apop_estimate(apop_vector_to_data(o), apop_normal);
    Diff(back_out->parameters->vector->data[0] , 1.1 , tol2);
    Diff(back_out->parameters->vector->data[1] , 1.23 , tol2);
    gsl_vector_free(o);
}

double ran_uniform(double in, void *r){ return gsl_rng_uniform(r);}
double negate(double in){ return -in;}

void test_posdef(gsl_rng *r){
    for(size_t j=0; j < 30; j ++){
        int size = gsl_rng_uniform(r) *10+1;
        apop_data *d = apop_data_alloc(size, size);
        apop_map(d, .fn_dp=ran_uniform, .param=r, .inplace='y', .part='m');
        apop_matrix_to_positive_semidefinite(d->matrix); 
        assert(apop_matrix_is_positive_semidefinite(d->matrix));

        //start over, from the negation of where you just were 
        //(guaranteed neg definite, no?)
        gsl_matrix * neg = apop_matrix_map_all(d->matrix, negate);
        apop_matrix_to_positive_semidefinite(neg);
        assert(apop_matrix_is_positive_semidefinite(neg));
        gsl_matrix_free(neg);
        apop_data_free(d);
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


void test_pmf(){
    double x[] = {0, 0.2, 0 , 0.4, 1, .7, 0 , 0, 0};
    gsl_rng *r = apop_rng_alloc(1234);
    apop_data *d = apop_data_alloc();
    d->weights = apop_array_to_vector(x, 9);
    apop_model *mc = apop_model_copy(apop_pmf);
    Apop_model_add_group(mc, apop_pmf, .draw_index= 'y');
    mc->dsize=0;
    apop_model *m = apop_estimate(d, mc);
    gsl_vector *v = gsl_vector_calloc(d->weights->size);
    for (size_t i=0; i< 1e5; i++){
        double out;
        apop_draw(&out, r, m);
        (*gsl_vector_ptr(v, out))++;
    }
    apop_vector_normalize(d->weights);
    apop_vector_normalize(v);
    for (size_t i=0; i < v->size; i ++)
        Diff(d->weights->data[i], v->data[i], 1e-2);
    apop_model_free(m);
    apop_data_free(d);
    gsl_vector_free(v);
}

void test_arms(gsl_rng *r){
    gsl_vector *o = gsl_vector_alloc(3e5);
    apop_model *ncut = apop_model_set_parameters(apop_normal, 1.1, 1.23);
    ncut->draw = NULL; //testing the default.
    for(size_t i=0; i < 3e5; i ++)
        apop_draw(o->data+i, r, ncut);
    apop_data *ov = apop_vector_to_data(o);
    apop_model *back_out = apop_estimate(ov, apop_normal);
    Diff(back_out->parameters->vector->data[0] , 1.1, 1e-2)
    Diff(back_out->parameters->vector->data[1] , 1.23, 1e-2)

    apop_opts.verbose ++;
    apop_model *bcut = apop_model_set_parameters(apop_beta, 0.4, 0.43); //bimodal
    Apop_model_add_group(bcut, apop_arms, .model=bcut, .xl=1e-5, .xr=1-1e-5);
    bcut->draw = NULL; //testing the default.
    for(size_t i=0; i < 3e5; i ++)
        apop_draw((o->data)+i, r, bcut);
    ov = apop_vector_to_data(o);
    apop_model *back_outb = apop_estimate(ov, apop_beta);
    apop_data_free(ov);
    Diff(back_outb->parameters->vector->data[0] , 0.4, 1e-2)
    Diff(back_outb->parameters->vector->data[1] , 0.43, 1e-2)
    apop_opts.verbose --;
    apop_model *test_copying = apop_model_copy(back_outb);
    apop_model_free(ncut);
    apop_model_free(back_out);
    apop_model_free(back_outb);
    apop_model_free(test_copying);
}

void test_pmf_compress(gsl_rng *r){
    apop_data *d = apop_data_alloc();
    apop_text_alloc(d, 9, 1);
    d->vector = apop_array_to_vector((double []){12., 1., 2., 2., 1., 1., 2., 2., NAN}, 9);
    apop_text_fill(d, "Dozen", "Single", "Pair", "Pair",
                      "Single", "Single", "Pair", "Pair",
                      "Nada");


    apop_data_pmf_compress(d);

    assert(d->vector->data[0]==12);
    assert(d->vector->data[1]==1);
    assert(d->vector->data[2]==2);
    assert(gsl_isnan(d->vector->data[3]));
    assert(d->weights->data[0]==1);
    assert(d->weights->data[1]==3);
    assert(d->weights->data[2]==4);
    assert(d->weights->data[3]==1);
    assert(apop_strcmp(d->text[0][0], "Dozen"));
    assert(apop_strcmp(d->text[1][0], "Single"));
    assert(apop_strcmp(d->text[2][0], "Pair"));
    assert(apop_strcmp(d->text[3][0], "Nada"));

    apop_data *b = apop_data_alloc();
    b->vector = apop_array_to_vector((double []){1.1, 2.1, 2, 1, 1}, 5);
    apop_text_alloc(b, 5, 1);
    apop_text_add(b, 0, 0, "Type 1");
    apop_text_add(b, 1, 0, "Type 1");
    apop_text_add(b, 2, 0, "Type 1");
    apop_text_add(b, 3, 0, "Type 1");
    apop_text_add(b, 4, 0, "Type 2");
    Apop_row(b, 0, arow);
    apop_data *spec = apop_data_copy(arow);
    gsl_vector_set_all(spec->vector, 1);
    apop_data *c = apop_data_to_bins(b, .binspec=spec);
    apop_data_free(b);
    assert(apop_strcmp(c->text[0][0], "Type 1"));
    assert(apop_strcmp(c->text[1][0], "Type 1"));
    assert(apop_strcmp(c->text[2][0], "Type 2"));
    assert(c->weights->data[0]==2);
    assert(c->weights->data[1]==2);
    assert(c->weights->data[2]==1);
    apop_data_free(c);

    //I assert that if I use the default binspec returned by a call to apop_data_to_bins,
    //then re-binning with the binspec explicitly stated will give identical results.
    int i, dcount = 10000;
    apop_data *draws = apop_data_alloc(dcount);
    apop_model *norm = apop_model_set_parameters(apop_normal, 0, 1);
    for (i=0 ; i<dcount; i++)
        apop_draw(draws->vector->data+i, r, norm);
    apop_data_sort(draws);
    apop_data *drawcopy = apop_data_copy(draws);
    apop_data *binned = apop_data_to_bins(draws);
    apop_data *binnedc = apop_data_to_bins(drawcopy, .binspec=apop_data_get_page(draws, "<binspec>"));
    assert(binned->weights->size == binnedc->weights->size);
    for (i=0; i< binned->weights->size; i++){
        assert(binned->vector->data[i] == binnedc->vector->data[i]);
        assert(binned->weights->data[i] == binnedc->weights->data[i]);
    }
}

void test_vtables(){
    //run an updating to make sure that the vtable has been generated.
    apop_model *n = apop_model_set_parameters(apop_normal, 0, 1);
    apop_data *d = apop_model_draws(n);
    apop_model *out = apop_update(d, n, apop_normal);
    //did it use the vtable to see this is a Normal distribution?
    assert(out->log_likelihood == apop_normal->log_likelihood);
    //updating a distribution with data from itself:
    Diff(apop_data_get(n->parameters), apop_data_get(out->parameters), 1e-1);

    assert(apop_update_vtable_drop(apop_normal, apop_normal)==0);
    assert(apop_update_vtable_drop(apop_beta, apop_binomial)==0);
}

void test_weighted_regression(apop_data *d, apop_model *e){
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
    Apop_col_v(cp, 1, off);
    gsl_vector_add_constant(off, 20);
    apop_model *way_off = apop_estimate(cp, apop_ols);
    assert(apop_vector_distance(zero_off->info->vector, way_off->info->vector) < 1e-4);
    gsl_vector *zcov = apop_data_pack(apop_data_get_page(zero_off->parameters, "<covariance>"), .use_info_pages='y');
    gsl_vector *wcov = apop_data_pack(apop_data_get_page( way_off->parameters, "<covariance>"), .use_info_pages='y');
    assert(apop_vector_distance(zcov, wcov) < 1e-4);
    assert(apop_data_get(zero_off->parameters, 0, -1) - (apop_data_get(way_off->parameters, 0, -1)+20./3)  < 5e-3);
    assert(apop_data_get(zero_off->parameters, 1, -1) - apop_data_get(way_off->parameters, 1, -1)  < 1e-4);
    apop_model_free(zero_off);
    apop_model_free(way_off);
    apop_data_free(cp);
    apop_data_free(useme);
}

#define do_test(text, fn) {if (verbose) printf("%s:", text); \
                          fflush(NULL);                      \
                          fn;                                \
                          if (verbose) printf("\nPASS.  ");} 

int main(int argc, char **argv){
    char *srcdir = getenv("srcdir");
    if (srcdir){ //if defined, this is probably via an automake rule.
        char buf[10000]; 
        if (strcmp(srcdir, ".") && strcmp(srcdir, getcwd(buf, 10000)))
            apop_system("cp %s/*dat* %s/printing_sample .", srcdir, srcdir); //needed for make distcheck. No-op in many cases.
    }
    int  slow_tests = 0;
    apop_opts.thread_count = 2;
    char c, opts[]  = "sqt:";
    if (argc==1)
        printf("Sundry tests for Apophenia.\nRunning relatively faster tests.  To run slower tests (primarily simulated annealing), use -s.\nFor quieter output, use -q. To change thread count (default=2), use -t1, -t2, -t3, ...\n");
    while((c = getopt(argc, argv, opts))!=-1)
        if (c == 's')       slow_tests++;
        else if (c == 'q')  verbose--;
        else if (c == 't')  apop_opts.thread_count = atoi(optarg);

    //set up some global or common variables
    gsl_rng *r = apop_rng_alloc(8); 
    apop_data *d = apop_text_to_data("test_data2",0,1);
    apop_model *an_ols_model = apop_model_copy(apop_ols);
    Apop_model_add_group(an_ols_model, apop_lm, .want_expected_value= 1);
    apop_model *e  = apop_estimate(d, an_ols_model);

    do_test("vtables", test_vtables());
    do_test("test listwise delete", test_listwise_delete());
    do_test("rownames", test_rownames());
    do_test("apop_dot", test_dot());
    do_test("apop_jackknife", test_jackknife(r));
    do_test("test multivariate_normal", test_multivariate_normal());
    do_test("log and exponent", log_and_exp(r));
    do_test("split and stack test", test_split_and_stack(r));
    do_test("test probit and logit", test_probit_and_logit(r));
    do_test("test probit and logit again", test_probit_and_logit(r));
    do_test("test ML imputation", test_ml_imputation(r));
    do_test("test data compressing", test_pmf_compress(r));
    do_test("weighted regression", test_weighted_regression(d,e));
    do_test("offset OLS", test_ols_offset(r));
    do_test("default RNG", test_default_rng(r));
    do_test("test row set and remove", row_manipulations());
    do_test("test PMF", test_pmf());
    do_test("apop_pack/unpack test", apop_pack_test(r));
    do_test("test adaptive rejection sampling", test_arms(r));
    //do_test("test fix params", test_model_fix_parameters(r));
    do_test("positive definiteness", test_posdef(r));
    do_test("test binomial estimations", test_binomial(r));
    do_test("dummies and factors", dummies_and_factors());
    do_test("test vector/matrix realloc", test_resize());
    do_test("test_vector_moving_average", test_vector_moving_average());
    do_test("apop_estimate->dependent test", test_predicted_and_residual(e));
    do_test("apop_f_test and apop_coefficient_of_determination test", test_f(e));
    do_test("OLS test", test_OLS(r));
    do_test("test lognormal estimations", test_lognormal(r));
    do_test("test jackknife covariance", test_jack(r));
    do_test("database skew, kurtosis, normalization", test_skew_and_kurt(r));
    do_test("test_percentiles", test_percentiles());
    do_test("weighted moments", test_weigted_moments());
    do_test("multivariate gamma", test_mvn_gamma());
    do_test("Inversion", test_inversion(r));
    do_test("apop_matrix_summarize", test_summarize());
    do_test("apop_linear_constraint", test_linear_constraint());
    do_test("transposition", test_transpose());
    do_test("test unique elements", test_unique_elements());
    if (slow_tests){
        if (verbose) printf("\tSlower tests:\n");
        do_test("Test score (dlog likelihood) calculation", test_score());
    }
    printf("\nApophenia has passed all of the sundry tests. Yay.\n");
    apop_db_close();
}
