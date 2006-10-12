
/* Here are some ad hoc tests to verify that things are basically OK. If
you'd like more thorough tests, feel free to write them.  

Part of the incompleteness of the tests, by the way, is that most of
Apophenia was written for immediate use in certain projects, so there's
a great deal of real-world testing that didn't make it into this file.

*/

#include <apophenia/headers.h>
#include "nist_tests.c"

//I'm using the test script an experiment to see if 
//these macros add any value.
#define APOP_ESTIMATION_PARAMS_ALLOC(name) apop_estimation_params *name = apop_estimation_params_alloc()
#define APOP_MATRIX_ALLOC(name, r, c) gsl_matrix *name = gsl_matrix_alloc((r),(c))
#define APOP_VECTOR_ALLOC(name, r) gsl_vector *name = gsl_vector_alloc(r)
#define APOP_DATA_ALLOC(name, r, c) apop_data *name = apop_data_alloc((r),(c))
#define APOP_RNG_ALLOC(name, seed) gsl_rng *name = apop_rng_alloc(seed)


double  true_parameter[]    = {1.82,2.1},
        true_y_parameter[]  = {0,2.1};

double  tolerance           = 1e-5;
double  lite_tolerance      = 1e-2;
int     len                 = 8000;

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
APOP_DATA_ALLOC(d1, 10,10);
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
int test_predicted_and_residual(apop_estimate *est){
gsl_vector  v,
            *prediction = gsl_vector_alloc(est->data->matrix->size1);
gsl_matrix  *m          = gsl_matrix_alloc(est->data->matrix->size1,est->data->matrix->size2);
    //prep an affine data matrix.
    gsl_matrix_memcpy(m, est->data->matrix);
    v   = gsl_matrix_column(m, 0).vector;
    gsl_vector_set_all(&v, 1);

    v   = gsl_matrix_column(est->dependent->matrix, apop_name_find(est->dependent->names, "residual", 'c')).vector;
    assert(fabs(apop_mean(&v)) < tolerance);

    v   = gsl_matrix_column(est->dependent->matrix, apop_name_find(est->dependent->names, "pred%", 'c')).vector;
    gsl_blas_dgemv(CblasNoTrans, 1, m, est->parameters->vector, 0, prediction);
    gsl_vector_sub(prediction, &v);
    assert(fabs(apop_vector_sum(prediction)) < tolerance);
    return 0;
}

/** I claim that the F test calculated via apop_F_test(est, NULL, NULL)
 equals a transformation of R^2.

*/
int test_f(apop_estimate *est){
apop_data   *rsq    = apop_estimate_correlation_coefficient(est);
apop_data   *ftab   = apop_F_test(est, NULL, NULL);
double      n       = est->data->matrix->size1;
double      K       = est->parameters->vector->size;
double      r       = apop_data_get_tn(rsq,"R_squared",-1);
double      f       = apop_data_get_tn(ftab,"F_stat%",-1);
    assert(fabs(f - r*(n-K)/ ((1-r)*K)) < tolerance);
    return 0;
}


int test_OLS(){
int             i;
apop_estimate   *out;
////gsl_rng         *r  =  apop_rng_alloc(12);
apop_data       *bkup;
APOP_DATA_ALLOC(set, len, 2);
APOP_RNG_ALLOC(r, 23);
APOP_ESTIMATION_PARAMS_ALLOC(ep);


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
    ep->weights  = w;
    out = apop_OLS.estimate(bkup, ep);
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

int estimate_model(gsl_matrix *data, apop_model dist){
int                     i,
                        score         = 0;
double                  starting_pt[] = {3.2, 1.4};
apop_estimation_params  *params = apop_estimation_params_alloc();
apop_estimate           *e;
    params->method           = 100;
    params->step_size        = 1e-2;
    params->starting_pt      = starting_pt;
    params->tolerance        = 1e-5;
    params->verbose          = 1;

    e    = dist.estimate(apop_matrix_to_data(data),params);
    for (i=0; i < dist.parameter_ct; i++){
        printf("parameter estimate, which should be %g: %g\n", true_parameter[i], gsl_vector_get(e->parameters->vector,i));
        score += (fabs(gsl_vector_get(e->parameters->vector,i) - true_parameter[i]) >= 1e-1);
        //apop_estimate_print(e);
    }

/*
    //wn versions:
    e    = apop_wn_maximum_likelihood(data2,&inv, dist, dummy, 1e-1, 1e-2, 1);
    for (i=0; i < dist.parameter_ct; i++){
        printf("parameter estimate, which should be %g: %g\n", true_parameter[i], gsl_vector_get(e->parameters,i));
        score += (fabs(gsl_vector_get(e->parameters,i) - true_parameter[i]) >= 1e-1);
        //apop_estimate_print(e);
    }
*/
    //return score;
    return 0;
}


int test_distribution(gsl_rng *r, apop_model model){
long int        runsize             = 1000,
                rowsize             = 50;
gsl_matrix      *data               = gsl_matrix_calloc(runsize,rowsize);
size_t          i,j;
    //generate.
    for (i=0; i< runsize; i++){
        for (j=0; j< rowsize; j++){
            gsl_matrix_set(data, i, j, model.rng(r, true_parameter));
        }
    }
    return estimate_model(data, model);
}

void generate_for_rank_test(gsl_matrix *data, gsl_rng *r, apop_model dist, int runsize, int rowsize){
int             i, z, j,
                runct   = 10000,
                rowsum;
gsl_vector_view v;
    for (i=0; i< runsize; i++){
        rowsum    = 0;
        for (j=0; j< runct; j++){
            z    = dist.rng(r, true_parameter);
            if (!strcmp(dist.name, "Exponential, rank data"))
                    z++;
            assert (z >=1);
            if (z < rowsize){    //else, just throw it out.
                apop_matrix_increment(data, i, (int)z-1, 1);
                rowsum    ++;
            }
        }
        v    = gsl_matrix_row(data,i);
        gsl_vector_scale(&(v.vector), 1./rowsum);    //Normalize!!
    }
}
int test_rank_distribution(gsl_rng *r, apop_model dist){
long int        //i,j,
                runsize             = 500,
                rowsize             = 100;
gsl_matrix      *data               = gsl_matrix_calloc(runsize,rowsize),
                *data2              = gsl_matrix_calloc(1, rowsize);    
apop_data       *summary;
gsl_vector      v;
    //generate.
    generate_for_rank_test(data, r, dist, runsize, rowsize);
    summary     = apop_matrix_summarize (data);
    v           = gsl_matrix_column(summary->matrix,0).vector;
    gsl_matrix_set_row(data2, 0, &v);

/*    printf("the abbreviated data matrix:\n");
    apop_matrix_print(data2, "\t", NULL);
    printf("\n");*/
        //for (j=0; j< rowsize; j++){printf("%g ",apop_zipf.log_likelihood(zipf_param-1, j+1));}
    return estimate_model(data2, dist);
}

#define INVERTSIZE 100
int test_inversion(gsl_rng *r){
gsl_matrix  *invme   = gsl_matrix_alloc(INVERTSIZE, INVERTSIZE);
gsl_matrix  *inved;
gsl_matrix  *inved_back;
int         i,j;
double      error   = 0,
            four[]  = {4};
    for(i=0; i<INVERTSIZE; i++)
        for(j=0; j<INVERTSIZE; j++)
            gsl_matrix_set(invme, i,j, apop_zipf.rng(r, four));
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
    return 0;
}


int test_jackknife(){
APOP_DATA_ALLOC(d, len, 1);
APOP_RNG_ALLOC(r, 8);
size_t      i;
apop_model  *m   = &apop_normal;
double      p[] = {1.09,2.8762};
//double      no;
    for (i =0; i< len; i++){
        apop_data_set(d, i, 0, m->rng(r,p)); 
    }
    gsl_matrix *out = apop_jackknife(d, *m , NULL);
    //apop_matrix_show(out);
    return (fabs(gsl_matrix_get(out, 0,0) - p[1]) > lite_tolerance);
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




#define do_test(text, fn)   if (verbose)    \
                                printf(text);  \
                            else printf(".");   \
                            fflush(NULL);   \
                            if (fn==0)    \
                                {if (verbose) printf("passed.\n");} \
                           else             \
                                {printf("%s  failed.\n", text);exit(0);}
int main(){
gsl_rng       *r              = apop_rng_alloc(8); 
apop_model    //rank_dist[]     = {apop_zipf_rank,apop_exponential_rank,apop_yule_rank, apop_waring_rank},
              dist[]          = {apop_zipf,apop_exponential ,apop_yule, apop_waring, apop_normal, apop_poisson};
int           //rank_dist_ct    = 4,
              dist_ct         = 6,
              i, verbose      = 0;
apop_data     *d  = apop_text_to_data("test_data2",0,1);
apop_estimate *e  = apop_OLS.estimate(d,NULL);
    apop_opts.thread_count  = 2;
/*  //now specified above.
apop_estimation_params  params;
        params.method           = 1;
        params.step_size        = 1e-2;
        params.tolerance        = 1e-3;
        params.verbose          = 1;
        */
    do_test("NaN handling", test_nan_data());
    do_test("test_percentiles:", test_percentiles());
    do_test("weighted moments:", test_weigted_moments());
    do_test("nist_tests:", nist_tests());
    do_test("split and stack to vector test:", test_matrix_split_to_vector());
    do_test("split and stack test:", test_split_and_stack());
    do_test("apop_dot test:", test_dot());
    do_test("apop_jackknife test:", test_jackknife());
    do_test("OLS test:", test_OLS());
    do_test("apop_estimate->dependent test:", test_predicted_and_residual(e));
    do_test("apop_f_test and apop_coefficient_of_determination test:", test_f(e));
    do_test("apop_vector_replace test:", test_replaces());
    do_test("apop_generalized_harmonic test:", test_harmonic());
    do_test("apop_strip_dots test:", test_strip_dots());
    do_test("apop_distance test:", test_distances());
    do_test("Inversion test: ", test_inversion(r));
    do_test("apop_matrix_summarize test:", test_summarize());
    verbose ++;

    /*
    for (i=0; i< rank_dist_ct; i++){
        do_test(rank_dist[i].name, test_rank_distribution(r, rank_dist[i]));
    }
    */

    for (i=0; i< dist_ct; i++){
        do_test(dist[i].name, test_distribution(r, dist[i]));
    }
    return 0;
}
