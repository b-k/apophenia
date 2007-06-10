
/* Here are some ad hoc tests to verify that things are basically OK. If
you'd like more thorough tests, feel free to write them.  

Part of the incompleteness of the tests, by the way, is that most of
Apophenia was written for immediate use in certain projects, so there's
a great deal of real-world testing that didn't make it into this file.

*/

#include <apop.h>
#include "nist_tests.c"

//I'm using the test script an experiment to see if 
//these macros add any value.
#define APOP_ep_ALLOC(name) apop_ep *name = apop_ep_alloc()
#define APOP_MATRIX_ALLOC(name, r, c) gsl_matrix *name = gsl_matrix_alloc((r),(c))
#define APOP_VECTOR_ALLOC(name, r) gsl_vector *name = gsl_vector_alloc(r)
#define APOP_DATA_ALLOC(name, r, c) apop_data *name = apop_data_alloc(0, (r),(c))
#define APOP_RNG_ALLOC(name, seed) gsl_rng *name = apop_rng_alloc(seed)


double  true_parameter_v[]    = {1.82,2.1};

apop_data *true_parameter;
apop_model *true_params;

double  tolerance           = 1e-5;
double  lite_tolerance      = 1e-2;
int     len                 = 8000;
int     verbose             = 0;

int test_percentiles(){
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
    return 0;
}

//This tests the database-side functions.
int test_skew_and_kurt(){
  gsl_rng *r  = apop_rng_alloc(time(0));
  int     i;
    apop_table_exists("t",1);
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
    assert (fabs(apop_vector_kurt(v) -apop_query_to_float("select kurt(vals) from t"))<1e-6);
    apop_table_exists("t",1);
    return 0;
}

int test_listwise_delete(){
  apop_data *t1 = apop_data_calloc(0, 10,10);
  apop_data *t1c = apop_data_listwise_delete(t1);
  assert(t1c->matrix->size1==10);
  assert(t1c->matrix->size2==10);
  t1->vector    = gsl_vector_calloc(10);
  apop_data_set(t1, 4,-1, GSL_NAN);
  apop_data *t2c = apop_data_listwise_delete(t1);
  assert(t2c->matrix->size1==9);
  apop_data_set(t1, 4,-1, GSL_NAN);
  apop_data_set(t1, 7,-1, GSL_NAN);
  apop_data *t3c = apop_data_listwise_delete(t1);
  assert(t3c->matrix->size1==8);
  APOP_COL(t1, 7, v)
  gsl_vector_set_all(v, GSL_NAN);
  assert(!apop_data_listwise_delete(t1));
    return 0;
}

int test_nan_data(){

    apop_text_to_db("test_data_nans", "nandata", 0,1, NULL);
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


    return 0;
}


static void wmt(gsl_vector *v, gsl_vector *v2, gsl_vector *w, gsl_vector *av, gsl_vector *av2, double mean){
    assert(apop_vector_mean(v) == apop_vector_weighted_mean(v,NULL));
    assert(apop_vector_mean(av) == apop_vector_weighted_mean(v,w));
    assert(apop_vector_weighted_mean(v,w) == mean);
    assert(fabs(apop_vector_var(v) - apop_vector_weighted_var(v,NULL))<1e-5);
    assert(fabs(apop_vector_covar(v,v2) - apop_vector_weighted_cov(v,v2,NULL))<1e-5);
    assert(fabs(apop_vector_var(av) - apop_vector_weighted_var(v,w))<1e-5);
    assert(fabs(apop_vector_covar(av,av2) - apop_vector_weighted_cov(v,v2,w))<1e-5);
    assert(fabs(apop_vector_skew(av) - apop_vector_weighted_skew(v,w))<1e-5);
    assert(fabs(apop_vector_kurt(av) - apop_vector_weighted_kurt(v,w))<1e-5);
}

int test_weigted_moments(){
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
    return 0;
}

int test_split_and_stack(){
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
    }
    return 0;
}



/** I claim that the mean residual is near zero, and that the predicted
  value is \f$X'\beta\f$.

  */
int test_predicted_and_residual(apop_model *est){
gsl_vector  v,
            *prediction = gsl_vector_alloc(est->data->matrix->size1);
gsl_matrix  *m          = gsl_matrix_alloc(est->data->matrix->size1,est->data->matrix->size2);
    //prep an affine data matrix.
    gsl_matrix_memcpy(m, est->data->matrix);
    v   = gsl_matrix_column(m, 0).vector;
    gsl_vector_set_all(&v, 1);

    v   = gsl_matrix_column(est->expected->matrix, apop_name_find(est->expected->names, "residual", 'c')).vector;
    assert(fabs(apop_mean(&v)) < tolerance);

    v   = gsl_matrix_column(est->expected->matrix, apop_name_find(est->expected->names, "pred", 'c')).vector;
    gsl_blas_dgemv(CblasNoTrans, 1, m, est->parameters->vector, 0, prediction);
    gsl_vector_sub(prediction, &v);
    assert(fabs(apop_vector_sum(prediction)) < tolerance);
    return 0;
}

/** I claim that the F test calculated via apop_F_test(est, NULL, NULL)
 equals a transformation of R^2.

*/
int test_f(apop_model *est){
apop_data   *rsq    = apop_estimate_correlation_coefficient(est);
apop_data   *contr  = apop_data_alloc(0,0,0);
apop_data   *ftab   = apop_F_test(est, contr);
double      n       = est->data->matrix->size1;
double      K       = est->parameters->vector->size;
double      r       = apop_data_get_ti(rsq,"R.squared",-1);
double      f       = apop_data_get_ti(ftab,"F.stat",-1);
    assert(fabs(f - r*(n-K)/ ((1-r)*K)) < tolerance);
    return 0;
}


int test_OLS(){
int             i;
apop_model   *out;
////gsl_rng         *r  =  apop_rng_alloc(12);
apop_data       *bkup;
APOP_DATA_ALLOC(set, len, 2);
APOP_RNG_ALLOC(r, 23);


for(i=0; i< len; i++){
    apop_data_set(set, i, 1, 100*(gsl_rng_uniform(r)-0.5));
    apop_data_set(set, i, 0, -1.4 + apop_data_get(set,i,1)*2.3);
}
    bkup    = apop_data_copy(set);

    out = apop_OLS.estimate(set, NULL);
//    apop_estimate_print(out);
    assert(fabs(apop_data_get(out->parameters, 0,-1) - -1.4) < tolerance);
    assert(fabs(apop_data_get(out->parameters, 1,-1) - 2.3) < tolerance);

APOP_VECTOR_ALLOC(w, set->matrix->size1);
    gsl_vector_set_all(w, 14);
    bkup->weights  = w;
    out = apop_OLS.estimate(bkup, NULL);
    assert(fabs(apop_data_get(out->parameters, 0,-1) - -1.4) < tolerance);
    assert(fabs(apop_data_get(out->parameters, 1,-1) - 2.3) < tolerance);
    return 0;
}


int is_neg(double in){
    return in < 0;
}

int test_replaces(void){
gsl_vector *v   = gsl_vector_calloc(3);
    gsl_vector_set(v, 2, 2.);
    assert(apop_vector_sum(v) == 2.);

    apop_vector_replace(v, apop_double_is_zero, -2);
    assert(apop_vector_sum(v) == -2.);

    apop_vector_replace(v, is_neg, GSL_POSINF);
    assert(apop_vector_sum(v) == GSL_POSINF);

    apop_vector_replace(v, gsl_isinf, 0);
    assert(apop_vector_sum(v) == 2);

    gsl_matrix *m   = gsl_matrix_calloc(3,2);
    gsl_matrix_set(m, 2, 1, 2.);
    assert(apop_matrix_sum(m) == 2.);

    apop_matrix_replace(m, apop_double_is_zero, -2);
    assert(apop_matrix_sum(m) == -2.*5 +2);

    apop_matrix_replace(m, is_neg, GSL_POSINF);
    assert(apop_matrix_sum(m) == GSL_POSINF);

    apop_matrix_replace(m, gsl_isinf, 0);
    assert(apop_matrix_sum(m) == 2);
return 0;
}

int test_strip_dots(void){
    /* 0: replace all dots with _
      1: everything before the last dot.
      2: everything after the first dot.
      */
char teapot[]   = "tea.pot";
char many_dots[]   = "tea.pot.csv";
char *out;
    out = apop_strip_dots(teapot, 0);
    assert(!strcmp(out, "tea_pot"));
    out = apop_strip_dots(teapot, 1);
    assert(!strcmp(out, "tea"));
    out = apop_strip_dots(teapot, 2);
    assert(!strcmp(out, "pot"));
    out = apop_strip_dots(many_dots, 0);
    assert(!strcmp(out, "tea_pot_csv"));
    out = apop_strip_dots(many_dots, 1);
    assert(!strcmp(out, "tea.pot"));
    out = apop_strip_dots(many_dots, 2);
    assert(!strcmp(out, "pot.csv"));
    return 0;
}

/** test the generalized harmonic summing thing, \ref apop_generalized_harmonic.
*/
int test_harmonic(){
double		out;
int		count = 0;
	out	= apop_generalized_harmonic(270, 0.0);
	if(out !=270){
		printf("Generalized harmonic(270,0) should be 270, but it's %g. Fail.\n", out);
		count++;
	}
	out	= apop_generalized_harmonic(370, -1.0);
	if(out !=370*371/2){
		printf("Generalized harmonic(370,-1) should be 370*371/2, but it's %g. Fail.\n", out);
		count++;
	}
	out	= apop_generalized_harmonic(12, -1.0);
	if(out !=12*13/2){
		printf("Generalized harmonic(12,-1) should be 12*13/2, but it's %g. Fail.\n", out);
		count++;
	}
	return count;
}

void estimate_model(gsl_matrix *data, apop_model dist, int method){
int                     i;
double                  starting_pt[] = {3.2, 1.4};
apop_mle_params  *params = apop_mle_params_alloc(apop_matrix_to_data(data), dist);
    params->method           = method;
    params->step_size        = 1e-1;
    params->starting_pt      = starting_pt;
    params->tolerance        = 1e-4;
    params->verbose          = 0;
    params->t_initial       = 1;
    params->t_min           = .5;
    params->k               = 1.8;
    params->use_score       = 0;
    params->want_cov        = 0;
    apop_model *e    = apop_estimate(apop_matrix_to_data(data),*params->model);
    e   = apop_estimate_restart(e, method ? 0 : 1, 1);
    if (verbose)
        for (i=0; i < e->parameters->vector->size; i++)
            printf("parameter estimate, which should be %g: %g\n",
            true_parameter_v[i], gsl_vector_get(e->parameters->vector,i));
    for (i=0; i < e->parameters->vector->size; i++)
         assert(fabs(gsl_vector_get(e->parameters->vector,i) -
         true_parameter_v[i]) <= lite_tolerance *5);
}


void test_distribution(gsl_rng *r, apop_model model){
long int        runsize             = 1000,
                rowsize             = 100;
gsl_matrix      *data               = gsl_matrix_calloc(runsize,rowsize);
size_t          i,j;
     printf("\n");
    //generate.
    for (i=0; i< runsize; i++)
        for (j=0; j< rowsize; j++)
            model.draw(gsl_matrix_ptr(data,i,j), r, true_params);
    estimate_model(data, model,0);
    estimate_model(data, model,1);
//    estimate_model(data, model,5); //works, but takes forever.
}

#define INVERTSIZE 100
int test_inversion(gsl_rng *r){
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
    if (error > 1e-5) {
        printf("inversion error is too big: %g. Fail.\n", error); return 1;
    }
    return 0;
}

int test_summarize(){
gsl_matrix      *m;
apop_data       *s;
double          t, v;
    apop_text_to_db("test_data", "td", 0,1,NULL );
    m    = apop_query_to_matrix("select * from td");
    s    = apop_matrix_summarize(m);
    //apop_matrix_print(s,"\t", NULL);
    t    = gsl_matrix_get(s->matrix, 1,0);
    if (t !=3) {
        printf("apop_summarize failed to take a simple mean: %g should be three. Fail.\n", t); return 1;
        }
    t    = gsl_matrix_get(s->matrix, 2, 1);
    v    = sqrt((2*2 +3*3 +3*3 +4.*4.)/3.);
    if (t != v) {
        printf("apop_summarize failed to calcuate a std deviation: %g should be %g. Fail.\n", t,v); return 1;
        }
    return 0;
}

/** Just a 3-4-5 triangle */
int test_distances(){
    int v = apop_opts.verbose;
    apop_opts.verbose = -1;
    gsl_vector *v1 = gsl_vector_alloc(2);
    gsl_vector *v2 = gsl_vector_alloc(2);
    gsl_vector *v3 = gsl_vector_calloc(3);
    gsl_vector_set(v1, 0,2);
    gsl_vector_set(v1, 1,2);
    gsl_vector_set(v2, 0,5);
    gsl_vector_set(v2, 1,6);
    assert(apop_vector_distance(v1,v3) == 0);
    assert(apop_vector_distance(v1,v2) == 5.);
    assert(apop_vector_grid_distance(v2,v3) == 0);
    assert(apop_vector_grid_distance(v1,v2) == 7.);
    apop_opts.verbose = v;
    return 0;
}

void jtest(apop_model m, double *pv){
  APOP_DATA_ALLOC(d, len, 1);
  APOP_RNG_ALLOC(r, 8);
  apop_data  *p  = apop_data_alloc(0,0,0);
  p->vector      = apop_array_to_vector(pv, 2);
  apop_model*pp = apop_model_copy(m);
      pp->parameters    = p;
  size_t      i;
    for (i =0; i< len; i++)
        m.draw(apop_data_ptr(d, i, 0), r, pp); 
    apop_data *out = apop_jackknife_cov(d, m);
    //apop_data_show(out);
    //printf("%g\n",  2*gsl_pow_2(pv[1])/(len-1));
    //fflush(NULL);
    //Notice that the jackknife just ain't a great estimator here.
    assert (fabs(apop_data_get(out, 0,0) - pv[1]/len) < lite_tolerance
                && fabs(apop_data_get(out, 1,1) - 2*gsl_pow_2(pv[1])/(len-1)) < lite_tolerance*100);
    apop_data *out2 = apop_bootstrap_cov(out, m, r, 0);
    assert (fabs(apop_data_get(out2, 0,0) - pv[1]/len) < lite_tolerance
                && fabs(apop_data_get(out2, 1,1) - 2*gsl_pow_2(pv[1])/(len-1)) < lite_tolerance);
    apop_data_free(d);
}

int test_jackknife(){
  double      pv[] = {3.09,2.8762};
    jtest(apop_normal, pv);
    return 0;
}


int test_dot(){
apop_data *d1   = apop_text_to_data("test_data2",0,1); // 100 x 2
apop_data *d2   = apop_text_to_data("test_data2",0,1); // 100 x 2
apop_data *d3   = apop_dot(d1, d2, 0, 1);
//apop_data *d4   = apop_dot(d1, d2, 'p', 0);
gsl_vector  v1 = gsl_matrix_row(d1->matrix, 0).vector;
apop_data *d5   = apop_vector_to_data(&v1); // 2 x 1
apop_data *d7   = apop_dot(d5, d5, 0, 0);
    assert(apop_data_get(d7, 0, -1) == apop_data_get(d3, 0,0));
apop_data *d8   = apop_dot(d1, d5, 0, 0);
apop_data *d9   = apop_dot(d5, d1, 1);
    gsl_vector_sub(d8->vector, d9->vector);
    assert(!apop_vector_sum(d8->vector));
    return 0;
}

/*
static void one_matrix_for_split_to_vector(int rows, int cols, gsl_rng *r){
gsl_matrix      *m   = gsl_matrix_alloc(rows,cols),
                *m2;
int             i,j;
    for (i=0; i< rows; i++)
        for (j=0; j< cols; j++)
            gsl_matrix_set(m, i,j, gsl_rng_uniform(r));
    m2   = apop_vector_split_to_matrix(apop_matrix_stack_to_vector(m),cols);
    for (i=0; i< rows; i++)
        for (j=0; j< cols; j++)
            assert(gsl_matrix_get(m,i,j) ==gsl_matrix_get(m2,i,j));
}

int test_matrix_split_to_vector(){
gsl_rng         *r  = apop_rng_alloc(107);
    one_matrix_for_split_to_vector(10,17, r);
    one_matrix_for_split_to_vector(17,10, r);
    one_matrix_for_split_to_vector(1,17, r);
    one_matrix_for_split_to_vector(17,1, r);
    one_matrix_for_split_to_vector(10,10, r);
    return 0;
}
    do_int_test("split and stack to vector test:", test_matrix_split_to_vector());
*/



void apop_pack_test(gsl_rng *r){
  int i, j, k, v, m1,m2;
  apop_data *d, *dout;
  gsl_vector *mid;
    for (i=0; i< 10; i++){
        v   = gsl_rng_uniform(r) > 0.5 ? gsl_rng_uniform(r)*100 : 0;
        m1  = gsl_rng_uniform(r)*100;
        m2  = gsl_rng_uniform(r) > 0.5 ? gsl_rng_uniform(r)*100 : 0;
        d   = apop_data_alloc(v, m1, m2);
        if (v)
            for (j=0; j< v; j++)
                gsl_vector_set(d->vector, j, gsl_rng_uniform(r));
        if (m2)
            for (j=0; j< m1; j++)
                for (k=0; k< m2; k++)
                    gsl_matrix_set(d->matrix, j, k, gsl_rng_uniform(r));
        mid     = apop_data_pack(d);
        dout    = apop_data_unpack(mid, v, m1, m2);
        if (v)
            for (j=0; j< v; j++)
                assert(dout->vector->data[j] == d->vector->data[j]);
        if (m2)
            for (j=0; j< m1; j++)
                for (k=0; k< m2; k++)
                    assert(gsl_matrix_get(d->matrix, j, k) == gsl_matrix_get(dout->matrix, j, k));
        if (mid) gsl_vector_free(mid); 
        apop_data_free(d); apop_data_free(dout);
        }
        
}


void test_model_fix_parameters(){
  apop_data *pv   = apop_data_alloc(2,2,2);
  apop_data *mask = apop_data_calloc(2,2,2);
  size_t    i, ct = 1000;
  apop_data *d  = apop_data_alloc(0,ct,2);
  //apop_data *v  = apop_data_alloc(2,0,0);
  //apop_data *v2 = apop_data_alloc(2,0,0);
  //apop_mle_params *ep   = apop_mle_params_alloc(d, apop_multivariate_normal, NULL);
  gsl_rng *r    = apop_rng_alloc(10);
  double    draw[2];
    apop_data_set(pv, 0, -1, 8);
    apop_data_set(pv, 1, -1, 2);
    gsl_matrix_set(pv->matrix, 0,0,1);
    gsl_matrix_set(pv->matrix, 1,1,1);
    gsl_matrix_set(pv->matrix, 1,0,0.5);
    gsl_matrix_set(pv->matrix, 0,1,0.5);
  apop_model *pp = apop_model_copy(apop_multivariate_normal);
    pp->parameters  = pv;
    for(i=0; i< ct; i++){
        apop_multivariate_normal.draw(draw, r, pp);
        apop_data_set(d, i, 0, draw[0]);
        apop_data_set(d, i, 1, draw[1]);
        //*apop_data_ptr(d,i,1)    +=3;
    }

    gsl_matrix_set_all(mask->matrix, 1);
    apop_mle_params *mep1   = apop_model_fix_params(d, pv, mask, apop_multivariate_normal);
    mep1->method      = 0;
    apop_model   *e1  = mep1->model->estimate(NULL, mep1->model);
    gsl_vector_sub(e1->parameters->vector, pv->vector);
    assert(apop_vector_sum(e1->parameters->vector) < 1e-2);

  double    start2[] = {1,0,0,1};
    gsl_matrix_set_all(mask->matrix, 0);
    gsl_vector_set_all(mask->vector, 1);
    apop_mle_params *mep2   = apop_model_fix_params(d, pv, mask, apop_multivariate_normal);
    mep2->method      = 0;
    mep2->starting_pt = start2;
    apop_model   *e2  = mep2->model->estimate(NULL, mep2->model);
    gsl_matrix_sub(e2->parameters->matrix, pv->matrix);
    assert(apop_matrix_sum(e2->parameters->matrix) < 1e-2);
}

int test_linear_constraint(){
  gsl_vector *beta      = gsl_vector_alloc(2);
  gsl_vector *betaout  = gsl_vector_alloc(2);
    gsl_vector_set(beta, 0, 7);
    gsl_vector_set(beta, 1, 7);
  apop_data *contrasts  = apop_data_calloc(1,1,2);
    apop_data_set(contrasts, 0, 0, -1);
    apop_data_set(contrasts, 0, 1, -1);
    assert(fabs(apop_linear_constraint(beta, contrasts, 0, betaout) - sqrt(2*49)) < tolerance);
    assert(!apop_vector_sum(betaout));
    gsl_vector_set(beta, 0, 0);
    assert(fabs(apop_linear_constraint(beta, contrasts, 0, betaout) - sqrt(49/2.)) < tolerance);
    assert(!apop_vector_sum(betaout));
    assert(gsl_vector_get(betaout,0)==-7/2.);
    //inside corner: find the corner
  gsl_vector *beta2     = gsl_vector_alloc(3);
  gsl_vector *betaout2  = gsl_vector_alloc(3);
    gsl_vector_set(beta2, 0, 7);
    gsl_vector_set(beta2, 1, 7);
    gsl_vector_set(beta2, 2, 7);
  apop_data *contrasts2  = apop_data_calloc(3,3,3);
    apop_data_set(contrasts2, 0, 0, -1);
    apop_data_set(contrasts2, 1, 1, -1);
    apop_data_set(contrasts2, 2, 2, -1);
    assert(fabs(apop_linear_constraint(beta2, contrasts2, 0, betaout2) - sqrt(3*49)) < tolerance);
    assert(apop_vector_sum(betaout2)==0);
    //sharp corner: go to one wall.
    apop_data_set(contrasts2, 0, 1, 1);
    assert(fabs(apop_linear_constraint(beta2, contrasts2, 0, betaout2) - sqrt(2*49)) < tolerance);
    assert(gsl_vector_get(betaout2,0)==7);
    assert(gsl_vector_get(betaout2,1)==0);
    assert(gsl_vector_get(betaout2,2)==0);
    return 0;
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

void test_histograms(){
//create a million draws 
    int n = 5e5;
  double    mu = 2.8, sigmasq = 1.34;
  apop_data *d = apop_data_alloc(0,n,1);
  gsl_matrix *out   = gsl_matrix_alloc(n,1);
  int       i;
  gsl_rng   *r  = apop_rng_alloc(1234);
    for (i=0; i< n; i++)
        //gsl_matrix_set(d->matrix, i, 0, gsl_rng_uniform(r));
        gsl_matrix_set(d->matrix, i, 0, gsl_ran_gaussian(r, sqrt(sigmasq))+mu);
    apop_model   *hp = apop_histogram_params_alloc(d, 1000, NULL);
    for (i=0; i< n; i++){
        apop_histogram.draw(gsl_matrix_ptr(out, i,0), r, hp);
        assert(gsl_finite(gsl_matrix_get(out, i,0)));
    }
    apop_model *outparams   = apop_normal.estimate(apop_matrix_to_data(out), NULL);
    assert(fabs(outparams->parameters->vector->data[0]-mu) < 1e2);
    assert(fabs(outparams->parameters->vector->data[1]-sigmasq) < 1e2);
//    apop_plot_histogram(out, 100, NULL); 
}

int subtest_updating(apop_model prior, apop_model likelihood){
  int i, j, reps    = 10;
  int tsize         = 1e3;
  gsl_rng   *r      = apop_rng_alloc(2343);
  apop_data *p      = apop_data_alloc(2,0,0);
  apop_data *pb     = apop_data_alloc(1,0,0);
  apop_data *tdata  = apop_data_alloc(0,tsize,1);
    apop_data_set(p, 0, -1, 1);
    apop_data_set(p, 1, -1, 1);
    apop_data_set(pb, 0, -1, .55);
  apop_model   *prior_eps = apop_model_copy(prior);
    prior_eps->parameters   = p;
  apop_model   *pbp       = apop_model_copy(likelihood);
  pbp->parameters   = pb;
  apop_model   *outparamsbb;
  apop_model *almost_bern  = apop_model_copy(likelihood);
    strcpy(almost_bern->name, "Bernie's distribution");
    for (j=0; j< reps; j++){
        for (i=0; i< tsize; i++)
            likelihood.draw(apop_data_ptr(tdata,i,0), r, pbp);
        if (j==0)
            outparamsbb  = apop_update(tdata, *prior_eps, *almost_bern, NULL, r, 6e3, 0.53,1200);
        else
            outparamsbb  = apop_update(tdata, *outparamsbb, likelihood, NULL, r, 6e3, 0.53,1200);
        prior_eps =  apop_update(tdata, *prior_eps, likelihood, NULL, r, 0, 0,0);
    }
    float r1= prior_eps->parameters->vector->data[0]/ prior_eps->parameters->vector->data[1];
    printf("the alt:\n");
    for (i=0; i< tsize; i++)
        apop_histogram.draw(apop_data_ptr(tdata, i, 0), r, outparamsbb);
    apop_model *final= prior_eps->estimate(tdata, prior_eps);
    float r2= final->parameters->vector->data[0]/ final->parameters->vector->data[1];
    printf("%g\n", r1/r2);
    return 0;
}


int test_jack(){
  int i, draws     = 2000;
  apop_data *d  =apop_data_alloc(0, draws, 1);
  gsl_rng   *r  = apop_rng_alloc(2);
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

void test_updating (){
    subtest_updating(apop_gamma, apop_exponential);
    subtest_updating(apop_beta, apop_bernoulli);
}

#define do_int_test(text, fn)   if (verbose)    \
                                printf(text);  \
                            else printf(".");   \
                            fflush(NULL);   \
                            if (fn==0)    \
                                {if (verbose) printf(" passed.\n");} \
                           else             \
                                {printf("%s  failed.\n", text);exit(0);}

#define do_test(text, fn)   if (verbose)    \
                                printf(text);  \
                            else printf(".");   \
                            fflush(NULL);   \
                            fn;\
                            {if (verbose) printf(" passed.\n");} 

int main(){
    true_parameter  = apop_data_alloc(0,0,0);
    true_parameter->vector  = apop_array_to_vector(true_parameter_v, 2);
    true_params             = apop_model_copy(apop_gamma);//irrelevant.
    true_params->parameters = true_parameter;

  gsl_rng       *r              = apop_rng_alloc(8); 
  apop_model    dist[]          = {apop_gamma, apop_exponential, apop_normal, 
                                    apop_poisson, apop_zipf,apop_yule, apop_uniform};
  int           dist_ct         = 7,
                i;
  apop_data     *d  = apop_text_to_data("test_data2",0,1);
  apop_OLS_params *olp  = apop_OLS_params_alloc(d, apop_OLS);
    olp->want_expected_value    = 1;

  apop_model *e  = apop_estimate(d,*olp->model);
//    apop_opts.thread_count  = 2;
    do_test("test jackknife covariance", test_jack());
    do_test("test apop_update", test_updating());
    do_test("test apop_histogram model", test_histograms());
    do_test("test apop_data sort", test_data_sort());
    do_int_test("nist_tests:", nist_tests());
    do_test("apop_model_fix_parameters test:", test_model_fix_parameters());
    do_int_test("listwise delete", test_listwise_delete());
    do_int_test("NaN handling", test_nan_data());
    do_int_test("database skew and kurtosis", test_skew_and_kurt());
    do_int_test("test_percentiles:", test_percentiles());
    do_int_test("weighted moments:", test_weigted_moments());
    do_int_test("split and stack test:", test_split_and_stack());
    do_int_test("apop_dot test:", test_dot());
    do_int_test("OLS test:", test_OLS());
    do_int_test("apop_estimate->dependent test:", test_predicted_and_residual(e));
    do_int_test("apop_f_test and apop_coefficient_of_determination test:", test_f(e));
    do_int_test("apop_vector_replace test:", test_replaces());
    do_int_test("apop_generalized_harmonic test:", test_harmonic());
    do_int_test("apop_strip_dots test:", test_strip_dots());
    do_int_test("apop_distance test:", test_distances());
    do_int_test("Inversion test: ", test_inversion(r));
    do_int_test("apop_jackknife test:", test_jackknife());
    do_int_test("apop_matrix_summarize test:", test_summarize());
    do_int_test("apop_linear_constraint test:", test_linear_constraint());
    do_test("apop_pack/unpack test:", apop_pack_test(r));
    for (i=0; i< dist_ct; i++){
        do_test(dist[i].name, test_distribution(r, dist[i]));
    }
    printf("\nApophenia has passed all of its tests. Yay.\n");
    return 0;
}
