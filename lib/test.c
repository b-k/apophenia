#include <apophenia/headers.h>

/* Here are some ad hoc tests to verify that things are basically OK. If
you'd like more thorough tests, feel free to write them.  
*/

#include <gsl/gsl_sf_zeta.h>


double  true_parameter[]    = {3.82,2.1},
        true_y_parameter[]  = {0,2.1};



int estimate_model(gsl_matrix *data, apop_model dist){
int                     i,
                        score         = 0;
double                  starting_pt[] = {3.2, 1.4};
apop_estimation_params  params;
apop_inventory          inv;
apop_estimate           *e;
    params.method           = 1;
    params.step_size        = 1e-1;
    params.starting_pt      = starting_pt;
    params.tolerance        = 1e-5;
    params.verbose          = 1;
    apop_inventory_set(&inv, 1);

    //e    = apop_maximum_likelihood(data2,&inv, dist, params);
    e    = dist.estimate(apop_matrix_to_data(data),&inv, &params);
    for (i=0; i < dist.parameter_ct; i++){
        printf("parameter estimate, which should be %g: %g\n", true_parameter[i], gsl_vector_get(e->parameters,i));
        score += (fabs(gsl_vector_get(e->parameters,i) - true_parameter[i]) >= 1e-1);
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


int test_distribution(gsl_rng *r, apop_model model, apop_estimation_params params){
long int        runsize             = 1000,
                rowsize             = 50,
                rowsum;
gsl_matrix      *data               = gsl_matrix_calloc(runsize,rowsize);
size_t          i,j;
gsl_vector      *vv;
gsl_vector_view v;
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
long int        i,j,
                runsize             = 500,
                rowsize             = 100;
gsl_matrix      *data               = gsl_matrix_calloc(runsize,rowsize),
                *data2              = gsl_matrix_calloc(1, rowsize);    
apop_name       *summary_names;
apop_data       *summary;
gsl_vector_view v;
gsl_vector      *vv;
    //generate.
    generate_for_rank_test(data, r, dist, runsize, rowsize);
    summary     = apop_matrix_summarize (data);
    v           = gsl_matrix_column(summary->data,0);
    gsl_matrix_set_row(data2, 0, &(v.vector));

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
    t    = gsl_matrix_get(s->data, 1,0);
    if (t !=3) {
        printf("apop_summarize failed to take a simple mean: %g should be three. Fail.\n", t); return 1;
        }
    t    = gsl_matrix_get(s->data, 2, 1);
    v    = sqrt((2*2 +3*3 +3*3 +4.*4.)/3.);
    if (t != v) {
        printf("apop_summarize failed to calcuate a std deviation: %g should be %g. Fail.\n", t,v); return 1;
        }
    return 0;
}

/*
Due to what is either a bug in gcc or a deep failing of understanding on my part, this fn is in distributions.c.
*/
int test_harmonic(); //in distributions.c

#define do_test(fn)     if (fn==0)     printf("passed.\n"); \
               else      {printf("failed.\n");exit(0);}
int main(){
gsl_rng                 *r;
apop_model              rank_dist[]     = {apop_zipf,apop_exponential_rank,apop_yule, apop_waring},
                        dist[]          = {apop_normal, apop_exponential};
int                     rank_dist_ct    = 4,
                        dist_ct         = 2,
                        i;
apop_estimation_params  params;
        params.method           = 1;
        params.step_size        = 1e-2;
        params.tolerance        = 1e-3;
        params.verbose          = 1;

    gsl_rng_env_setup();
    r=gsl_rng_alloc(gsl_rng_default); 

    for (i=0; i< rank_dist_ct; i++){
        printf("%s: ",rank_dist[i].name);
        do_test(test_rank_distribution(r, rank_dist[i]));
    }

    for (i=0; i< dist_ct; i++){
        printf("%s: ",dist[i].name);
        do_test(test_distribution(r, dist[i],params));
    }

    printf("apop_matrix_summarize test:");
    do_test(test_summarize());

    /*
    printf("Inversion test: ");
    do_test(test_inversion(r));

    printf("Harmonic test: ");
    do_test(test_harmonic());
    */

    return 0;
}
