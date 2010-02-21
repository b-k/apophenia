#include <apop.h>
#include "nist_tests.c"

#define Diff(L, R, eps) apop_assert_void(fabs((L)-(R)<(eps)), 0, 's', "%g is too different from %g.", L, R);

//I'm using the test script an experiment to see if 
//these macros add any value.
#define APOP_VECTOR_ALLOC(name, r) gsl_vector *name = gsl_vector_alloc(r)
#define APOP_DATA_ALLOC(name, r, c) apop_data *name = apop_data_alloc(0, (r),(c))
#define APOP_RNG_ALLOC(name, seed) gsl_rng *name = apop_rng_alloc(seed)

//These test functions are also displayed in the documentation as examples.
#include "../eg/test_pruning.c"     // test_prune_cols()
#include "../eg/test_distances.c"
#include "../eg/test_kl_divergence.c"
#include "../eg/test_strip_dots.c"
#include "../eg/test_harmonic.c"
#include "../eg/test_updating.c" 
#include "../eg/test_fisher.c" 
#include "../eg/test_regex.c" 
#include "../eg/test_strcmp.c" 

//One-liners for mapply:
gsl_rng *r_global;
//void random_draw(double *in) { *in = gsl_rng_uniform(r_global);}
double nan_map(double in){return gsl_isnan(in);}

double  true_parameter_v[]    = {1.82,2.1};
//double  true_parameter_v[]    = {0.8,0.2};
double  tolerance           = 1e-5;
double  lite_tolerance      = 1e-2;
int     len                 = 8000;
int     verbose             = 1;

static void compare_mvn_estimates(apop_model *L, apop_model *R, double tolerance){
    gsl_vector_sub(L->parameters->vector, R->parameters->vector);
    gsl_matrix_sub(L->parameters->matrix, R->parameters->matrix);
    assert(fabs(apop_sum(L->parameters->vector)) + fabs (apop_matrix_sum(L->parameters->matrix)) < tolerance);
}

void test_ml_imputation(gsl_rng *r){
    size_t len = 2e4;
    int i,j;
    apop_data *fillme = apop_data_alloc(0, len, 3);
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
  apop_data   *data     = apop_data_alloc(0, len,1);
    for (i=0; i<10; i++){
        apop_model  *source   = apop_model_set_parameters(apop_normal, 
                                    gsl_ran_flat(r, -5, 5), gsl_ran_flat(r, .01, 5));
        for (j=0; j< len; j++)
            apop_draw(gsl_matrix_ptr(data->matrix, j, 0), r, source);
        apop_model *estme = apop_model_copy(apop_normal);
        Apop_model_add_group(estme, apop_mle, .method= APOP_SIMAN,.parent= estme);
        apop_model *out = apop_maximum_likelihood(data, estme);

        apop_model *straight_est = apop_estimate(data, apop_normal);
        assert(fabs(straight_est->parameters->vector->data[0]- source->parameters->vector->data[0])<1e-1);//rough, I know.
        assert(fabs(straight_est->parameters->vector->data[1]- source->parameters->vector->data[1])<1e-1);

        double sigsqn = gsl_pow_2(out->parameters->vector->data[1])/len;
        apop_data *cov = apop_data_get_page(out->parameters, "cov");
        assert(fabs(apop_data_get(cov, 0,0)-sigsqn) < 1e-3);
        assert(fabs(apop_data_get(cov, 1,1)-sigsqn/2) < 1e-3);
        assert(apop_data_get(cov, 0,1) + apop_data_get(cov, 0,1) < 1e-3);
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
/*    printf ("var %g %g\n", apop_var(v) ,apop_query_to_float("select var(vals) from t"));
    printf ("skew %g %g\n", apop_vector_skew(v) ,apop_query_to_float("select skew(vals) from t"));
    printf ("kurt %g %g\n", apop_vector_kurt(v) ,apop_query_to_float("select kurt(vals) from t"));*/
    assert (fabs(apop_var(v) -apop_query_to_float("select var(vals) from t"))<1e-6);
    assert (fabs(apop_vector_skew(v) -apop_query_to_float("select skew(vals) from t"))<1e-6);
    assert (fabs(apop_vector_kurt(v) -apop_query_to_float("select kurt(vals) from t"))<1e-5);
    apop_table_exists("t",1);
}

void test_listwise_delete(){
  int i, j;
  apop_data *t1 = apop_data_calloc(0, 10,10);
  apop_text_alloc(t1, 10, 10);
  for (i=0; i< 10; i++)
      for (j=0; j< 10; j++)
          apop_text_add(t1, i, j, "%i", i*j);
  apop_data *t1c = apop_data_listwise_delete(t1);
  assert(t1c->matrix->size1==10);
  assert(t1c->matrix->size2==10);
  assert(atoi(t1c->text[3][4])==12);

  t1->vector    = gsl_vector_calloc(10);
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
    strcpy(apop_opts.db_nan, "\\.");
  apop_data *d  = apop_query_to_data("select * from nandata");
    apop_opts.output_type ='d';//check that rownames come in OK, and NaNs written right.
    apop_data_print(d, "nantest");
    apop_data_free(d);
  apop_data *d2  = apop_query_to_data("select * from nantest");
    assert(gsl_isnan(apop_data_get_tt(d2,"second", "c")));
    assert(gsl_isnan(apop_data_get_tt(d2,"third", "b")));
    assert(!apop_data_get_tt(d2,"fourth", "b"));
    apop_data_free(d2);
    strcpy(apop_opts.db_nan, "XX");
}

static void wmt(gsl_vector *v, gsl_vector *v2, gsl_vector *w, gsl_vector *av, gsl_vector *av2, double mean){
    assert(apop_vector_mean(v) == apop_vector_weighted_mean(v,NULL));
    assert(apop_vector_mean(av) == apop_vector_weighted_mean(v,w));
    assert(apop_vector_weighted_mean(v,w) == mean);
    assert(fabs(apop_vector_var(v) - apop_vector_weighted_var(v,NULL))<1e-5);
    assert(fabs(apop_vector_cov(v,v2) - apop_vector_weighted_cov(v,v2,NULL))<1e-5);
    assert(fabs(apop_vector_var(av) - apop_vector_weighted_var(v,w))<1e-5);
    assert(fabs(apop_vector_cov(av,av2) - apop_vector_weighted_cov(v,v2,w))<1e-5);
    assert(fabs(apop_vector_skew_pop(av) - apop_vector_weighted_skew(v,w))<1e-5);
    assert(fabs(apop_vector_kurtosis_pop(av) - apop_vector_weighted_kurt(v,w))<1e-5);
}

void test_weigted_moments(){
  double        data[]      = {1,2,3};
  double        alldata[]   = {1,2,3};
  double        data3[]     = {3,2,1};
  double        alldata3[]  = {3,2,1};
  double        weights[]   = {1,1,1};
  gsl_vector    *v          = apop_array_to_vector(data, 3);
  gsl_vector    *v2         = apop_array_to_vector(data3, 3);
  gsl_vector    *w          = apop_array_to_vector(weights, 3);
  gsl_vector    *av         = apop_array_to_vector(alldata, 3);
  gsl_vector    *av2        = apop_array_to_vector(alldata3, 3);
    wmt(v,v2,w,av,av2,2);
  double data2[]       = {0,1,2,3,4};
  double alldata2[]    = {0,0,0,0,1,1,1,2,2,3};
  double data4[]       = {0,1,3,2,4};
  double alldata4[]    = {0,0,0,0,1,1,1,3,3,2};
  double weights2[]    = {4,3,2,1,0};
    v             = apop_array_to_vector(data2, 5);
    v2            = apop_array_to_vector(data4, 5);
    w             = apop_array_to_vector(weights2, 5);
    av            = apop_array_to_vector(alldata2, 10);
    av2           = apop_array_to_vector(alldata4, 10);
    wmt(v,v2,w,av,av2,1);
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
    assert(fabs(apop_mean(&v)) < tolerance);

    Apop_col_t(predict_tab, "pred", vv);
    gsl_blas_dgemv(CblasNoTrans, 1, m, est->parameters->vector, 0, prediction);
    gsl_vector_sub(prediction, vv);
    assert(fabs(apop_vector_sum(prediction)) < tolerance);
}

/** I claim that the F test calculated via apop_F_test(est, NULL, NULL)
 equals a transformation of R^2.
*/
void test_f(apop_model *est){
apop_data   *rsq    = apop_estimate_coefficient_of_determination(est);
apop_data   *ftab   = apop_F_test(est, NULL);
double      n       = est->data->matrix->size1;
double      K       = est->parameters->vector->size;
double      r       = apop_data_get_ti(rsq,"R.squared", 0);
double      f       = apop_data_get(ftab, .rowname="F.stat");
    assert(fabs(f - r*(n-K)/ ((1-r)*K)) < tolerance);
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
    assert(fabs(apop_data_get(out->parameters, 0,-1) - -1.4) < tolerance);
    assert(fabs(apop_data_get(out->parameters, 1,-1) - 2.3) < tolerance);

    gsl_vector *w = gsl_vector_alloc(set->matrix->size1);
    gsl_vector_set_all(w, 14);
    bkup->weights  = w;
    out = apop_estimate(bkup, apop_ols);
    assert(fabs(apop_data_get(out->parameters, 0,-1) - -1.4) < tolerance);
    assert(fabs(apop_data_get(out->parameters, 1,-1) - 2.3) < tolerance);
}

#define INVERTSIZE 100
void test_inversion(gsl_rng *r){
  gsl_matrix  *invme   = gsl_matrix_alloc(INVERTSIZE, INVERTSIZE);
  gsl_matrix  *inved;
  gsl_matrix  *inved_back;
  int         i,j;
  double      error   = 0;
  apop_data *four     = apop_data_alloc(1,0,0);
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
 
static void fill_p(apop_data *d, int v, int m1, int m2, int w, gsl_rng *r){
    int j, k;
    if (v)
        for (j=0; j< v; j++)
            gsl_vector_set(d->vector, j, gsl_rng_uniform(r));
    if (m2)
        for (j=0; j< m1; j++)
            for (k=0; k< m2; k++)
                gsl_matrix_set(d->matrix, j, k, gsl_rng_uniform(r));
    if (w)
        for (j=0; j< w; j++)
            gsl_vector_set(d->weights, j, gsl_rng_uniform(r));
}

static void check_p(apop_data *d, apop_data *dout, int v, int m1, int m2, int w){
    int j, k;
    if (v)
        for (j=0; j< v; j++)
            assert(dout->vector->data[j] == d->vector->data[j]);
    if (w)
        for (j=0; j< w; j++)
            assert(dout->weights->data[j] == d->weights->data[j]);
    if (m2)
        for (j=0; j< m1; j++)
            for (k=0; k< m2; k++)
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
        fill_p(d, v, m1, m2, w, r);
        int second_p = i %2;
        if (second_p){
            p2 = apop_data_add_page(d, apop_data_alloc(v, m1, m2), "second p");
            outp2 = apop_data_add_page(dout, apop_data_alloc(v, m1, m2), "second p");
            if (w) {p2->weights = gsl_vector_alloc(w);
                    outp2->weights = gsl_vector_alloc(w);}
            fill_p(p2, v, m1, m2, w, r);
        }
        mid     = apop_data_pack(d, .all_pages= second_p ? 'y' : 'n');
        apop_data_unpack(mid, dout, .all_pages= second_p ? 'y' : 'n');
        check_p(d, dout, v, m1, m2, w);
        if (second_p)
            check_p(d->more, dout->more, v, m1, m2, w);
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
    assert(fabs(apop_linear_constraint(beta, contrasts, 0) - sqrt(2*49)) < tolerance);
    assert(!apop_vector_sum(beta));
    gsl_vector_set(beta, 0, 0);
    gsl_vector_set(beta, 1, 7);
    assert(fabs(apop_linear_constraint(beta, contrasts, 0) - sqrt(49/2.)) < tolerance);
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
    assert(fabs(apop_linear_constraint(beta2, contrasts2, 0) - sqrt(3*49)) < tolerance);
    assert(apop_vector_sum(beta2)==0);
    //sharp corner: go to one wall.
    gsl_vector_set(beta2, 0, 7);
    gsl_vector_set(beta2, 1, 7);
    gsl_vector_set(beta2, 2, 7);
    apop_data_set(contrasts2, 0, 1, 1);
    assert(fabs(apop_linear_constraint(beta2, contrasts2, 0) - sqrt(2*49)) < tolerance);
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
  apop_data  *p  = apop_data_alloc(0,0,0);
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
assert ((fabs(apop_data_get(out, 0,0) - gsl_pow_2(pv[1])/len)) < lite_tolerance 
            && fabs(apop_data_get(out, 1,1) - gsl_pow_2(pv[1])/(2*len)) < lite_tolerance*100);
    apop_data *out2 = apop_bootstrap_cov(d, m);
    assert (fabs(apop_data_get(out2, 0,0) - gsl_pow_2(pv[1])/len) < lite_tolerance
                && fabs(apop_data_get(out2, 1,1) - gsl_pow_2(pv[1])/(2*len)) < lite_tolerance);
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
  int       len = 3e5;
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
    double phat = apop_data_get(out->parameters, 1,-1);
    double nhat = apop_data_get(out->parameters, 0,-1);
    if (verbose) printf("n: %i, p: %g, nhat: %g, phat: %g\n", n, p, phat, nhat);
    assert(apop_data_get(out->parameters, 0,-1) == n);
    assert(apop_data_get(out->parameters, 1,-1) - p < 1e-2);
}

void test_binomial(gsl_rng *r){
    size_t  i;
    double p = gsl_rng_uniform(r);
    int n     = gsl_ran_flat(r,1,1e5);
    apop_data *d = apop_data_alloc(0,1,n);
    for(i=0; i < n; i ++)
        apop_data_set(d, 0,i,(gsl_rng_uniform(r) < p));
    apop_model *out = apop_estimate(d, apop_binomial);
    apop_model *outm = apop_estimate(d, apop_multinomial);
    common_binomial_bit(out, n, p);
    common_binomial_bit(outm, n, p);
    apop_data_free(d);

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
    apop_data_free(d);
}

void db_to_text(){
    apop_db_open(NULL);
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
  gsl_vector *v = gsl_vector_alloc(100);
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
  int samples = 1e6;
  apop_data *data = apop_data_alloc(0, samples, true_params->size);
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
  gsl_matrix *m = gsl_matrix_alloc(20,15);
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
    assert (fabs(apop_matrix_sum(m) - 55 * 56 )<1e-6);
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
    gsl_vector *o = gsl_vector_alloc(2e6);
    apop_model *ncut = apop_model_set_parameters(apop_normal, 1.1, 1.23);
    ncut->draw = NULL; //forced to use the default.
    for(i=0; i < 2e6; i ++)
        apop_draw(o->data+i, r, ncut);
    apop_model *back_out = apop_estimate(apop_vector_to_data(o), apop_normal);
    assert(fabs(back_out->parameters->vector->data[0] - 1.1) < 1e-3);
    assert(fabs(back_out->parameters->vector->data[1] - 1.23) < 1e-3);
    gsl_vector_free(o);
}

double ran_uniform(double in, void *r){ return gsl_rng_uniform(r);}

void test_posdef(gsl_rng *r){
    r_global = r;
    size_t j;
    for(j=0; j < 30; j ++){
        int size = gsl_rng_uniform(r) *10+1;
        apop_data *d = apop_data_alloc(0, size, size);
        apop_map(d, .fn_dp=ran_uniform, .param=r, .inplace=1, .part='m');
        apop_matrix_to_positive_semidefinite(d->matrix);
        assert(apop_matrix_is_positive_semidefinite(d->matrix));
    }
}

void estimate_model(apop_data *data, apop_model dist, int method){
  int                     i;
  double                  starting_pt[] =  {0.6, 0.4};//{3.2, 1.4}; //{1.82,2.1};

    Apop_model_add_group(&dist, apop_mle, 
        .parent       = &dist,  .starting_pt = starting_pt,
        .method       = method, .verbose   =0,
        .step_size    = 1e-1,
        .tolerance    = 1e-2,   .k         = 1.8,
        .t_initial    = 1,      .t_min     = .5,
        .use_score    = 1,      .want_cov  = 'n'
        );
    Apop_model_add_group(&dist, apop_lm,  .want_cov = 'n');
    apop_model *e    = apop_estimate(data,dist);
    for (i=0; i < e->parameters->vector->size; i++)
        if (verbose) printf("parameter estimate, which should be %g: %g\n",
                                true_parameter_v[i], gsl_vector_get(e->parameters->vector,i));
    for (i=0; i < e->parameters->vector->size; i++)
         assert(fabs(gsl_vector_get(e->parameters->vector,i) - true_parameter_v[i]) <= lite_tolerance *5);
}
#if 0
void estimate_model(apop_data *data, apop_model dist, int method){
  int                     i;
  double                  starting_pt[] = {3.2, 1.4},
                          wishart_start[] = {8, 2, 1, 1, 9};
                          //wishart_start[] = {1,1,2,4,3};
//apop_mle_settings  *params = apop_mle_settings_alloc(apop_matrix_to_data(data), dist);
  /*  Apop_model_add_group(&dist, apop_mle, 
        .parent       = &dist,       .starting_pt = starting_pt,
        .method       = method,      .step_size = 1e-1,
        .tolerance    = 1e-4,        .k         = 1.8,
        .t_initial    = 1,           .t_min     = .5,
        .use_score    = 1,           .want_cov  = 0);*/
    Apop_model_add_group(&dist, apop_mle, 
        .parent       = &dist,       .starting_pt = starting_pt,
        .method       = method, .verbose =1);
        /*.tolerance    = 1e-4,        .k         = 1.8,
        .t_initial    = 1,           .t_min     = .5,
        .use_score    = 1,           .want_cov  = 0);*/
    if (!strcmp(dist.name, "Wishart distribution")){
        Apop_settings_add(&dist, apop_mle, method, APOP_CG_PR);
        Apop_settings_add(&dist, apop_mle, verbose, 1);
        //Apop_settings_add(&dist, apop_mle, method, APOP_SIMAN);
        Apop_settings_add(&dist, apop_mle, starting_pt, wishart_start);
    }
    apop_model *e    = apop_estimate(data,dist);
    /*
    if (Apop_settings_get_group(e, apop_mle) && !(!strcmp(dist.name,"poisson") || !strcmp(dist.name, "Uniform distribution"))){  //then it's an MLE
        apop_model *compare_me = apop_model_copy(dist);
        apop_mle_settings *p = Apop_settings_get_group(compare_me, apop_mle);
        p->method = p->method == APOP_SIMPLEX_NM ? APOP_CG_FR : APOP_SIMPLEX_NM;
        e   = apop_estimate_restart(e, compare_me);
    }
    */
    if (verbose)
        for (i=0; i < e->parameters->vector->size; i++)
            printf("parameter estimate, which should be %g: %g\n",
            true_parameter_v[i], gsl_vector_get(e->parameters->vector,i));
    for (i=0; i < e->parameters->vector->size; i++)
         assert(fabs(gsl_vector_get(e->parameters->vector,i) -
         true_parameter_v[i]) <= lite_tolerance *5);
}
#endif

/*Produce random data, then try to recover the original params */
  void test_one_distribution(gsl_rng *r, apop_model model, apop_model *true_params){
  long int        runsize             = 1e5,
                  rowsize             = 2;
  apop_data      *data; 
  size_t          i,j;
    //generate.
    if (!strcmp(model.name, "Wishart distribution")){
        data = apop_data_calloc(0,runsize,4);
        true_params->parameters->vector->data[0] = runsize-4;
        for (i=0; i< runsize; i++){
            Apop_row(data, i, v)
            true_params->draw(v->data, r, true_params);
            assert(!apop_vector_map_sum(v, nan_map));
        }
    } else if (!strcmp(model.name, "Dirichlet distribution")){
        data = apop_data_calloc(0,runsize,2);
        for (i=0; i< runsize; i++){
            Apop_row(data, i, v)
            true_params->draw(v->data, r, true_params);
            assert(!apop_vector_map_sum(v, nan_map));
        }
    } else {
        data = apop_data_calloc(0,runsize,rowsize);
        for (i=0; i< runsize; i++)
            for (j=0; j< rowsize; j++)
                true_params->draw(apop_data_ptr(data,i,j), r, true_params);
    }
    estimate_model(data, model,APOP_SIMPLEX_NM);
    estimate_model(data, model,APOP_CG_PR);
//    estimate_model(data, model,5); //works, but takes forever.
}

void test_distributions(gsl_rng *r){
  if (verbose) printf("\n");
  int         i;
  apop_model* true_params;
  apop_model  null_model      = {"the null model"};
  apop_model  dist[]          = {apop_dirichlet, /*apop_wishart,*/apop_poisson,  /*apop_gamma,*/ apop_exponential, apop_normal, apop_t_distribution, apop_f_distribution, 
                                    /*apop_zipf, apop_yule,*/ apop_uniform, null_model};

    for (i=0; strcmp(dist[i].name, "the null model"); i++){
        if (verbose) {printf("%s: ", dist[i].name); fflush(NULL);}
        true_params   = apop_model_copy(dist[i]);
        true_params->parameters = apop_line_to_data(true_parameter_v, 2,0,0);
        if (!strcmp(dist[i].name, "Wishart distribution")){
            true_params->parameters->matrix = gsl_matrix_alloc(2,2);
            true_params->parameters->vector->data[0] = 996;
            true_params->parameters->matrix->data[0] = .2;
            true_params->parameters->matrix->data[1] = .1;
            true_params->parameters->matrix->data[2] = .1;
            true_params->parameters->matrix->data[3] = .9;
        }
        test_one_distribution(r, dist[i], true_params);
        printf("Passed.\n");
    }
}

void test_rank_distributions(gsl_rng *r){
  int         i, j, k;
  long int        runsize             = 1000,
                  rowsize             = 100;
  double          index;
  apop_model* true_params;
  apop_model  null_model      = {"the null model"};
  apop_model  dist[]          = {/*apop_gamma,*/ apop_zipf, apop_yule, null_model};

    for (i=0; strcmp(dist[i].name, "the null model"); i++){
        if (verbose) {printf("%s: ", dist[i].name); fflush(NULL);}
        true_params   = apop_model_copy(dist[i]);
        true_params->parameters = apop_line_to_data(true_parameter_v, 2,0,0);
        //generate.
        apop_data *data = apop_data_calloc(0,runsize,rowsize);
            for (k=0; k< runsize; k++)
                for (j=0; j< 1e3; j++){ 
                    apop_draw(&index, r, true_params);
                    if (index <= rowsize)
                        apop_matrix_increment(data->matrix,k,index);
                }
        //estimate.
        Apop_model_add_group(dist+i, apop_rank)
        estimate_model(data, dist[i],APOP_SIMPLEX_NM);
        estimate_model(data, dist[i],APOP_CG_PR);
        apop_data_free(data);
        Apop_settings_rm_group(dist+i, apop_rank);
        //    estimate_model(data, model,APOP_SIMAN); //works, but takes forever.
        printf("Passed.\n");
    }
}

void test_pmf(){
   double x[] = {0, 0.2, 0 , 0.4, 1, .7, 0 , 0, 0};
   size_t i;
	apop_data *d= apop_data_alloc(0,0,0);
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
        assert(fabs(d->vector->data[i] - v->data[i]) < 1e-2);
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
}


#include "pmf_test.c"

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


    do_test("Kullback-Leibler divergence test", pack_test());
    exit(8);

    do_test("Kullback-Leibler divergence test", test_kl_divergence(r));
    do_test("apop_distance test", test_distances());
    do_test("test column pruning", test_prune_cols());
    do_test("test PMF", test_pmf());
    do_test("apop_pack/unpack test", apop_pack_test(r));
    do_test("test regex", test_regex());
    do_test("test adaptive rejection sampling", test_arms(r));
    do_test("test listwise delete", test_listwise_delete());
    do_test("test distributions", test_distributions(r));
    //do_test("test rank distributions", test_rank_distributions(r));
    do_test("test ML imputation", test_ml_imputation(r));
    do_test("test apop_update", test_updating(r));
    do_test("test fix params", test_model_fix_parameters(r));
    do_test("positive definiteness", test_posdef(r));
    do_test("test binomial estimations", test_binomial(r));
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
        if (verbose) printf("Slower tests:\n");
        do_test("Test score (dlog likelihood) calculation", test_score());
        do_test("test probit and logit", test_probit_and_logit(r));
    }
    printf("\nApophenia has passed all of its tests. Yay.\n");
    return 0;
}
