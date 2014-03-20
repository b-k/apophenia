/* These are tests of distributions. The basic idea is to 
 --assume a true set of parameters
 --generate a fake data set via a few thousand draws from your preferred model.
 --estimate the parameters of a new model using the fake data
 --assert that the estimated parameters are within epsilon of the true parameters.
*/
#include <apop.h>
#include <unistd.h>
#define Diff(L, R, eps) Apop_assert(fabs((L)-(R))<(eps), "%g is too different from %g (abitrary limit=%g).", (double)(L), (double)(R), eps);

#define Print_dot if(verbose){printf(".");fflush(NULL);}
#define is_t(d) !strcmp((d)->name, "t distribution")
#define is_bernie(d) !strcmp((d)->name, "Bernoulli distribution")
#define is_binom(d) !strcmp((d)->name, "Binomial distribution")
#define is_beta(d) !strcmp((d)->name, "Beta distribution")
#define is_poisson(d) !strcmp((d)->name, "Poisson distribution")

int verbose = 1;

//The MLE of the t distribution may have non-integer value (why not?)
//Because we started with an integer value, we have to find the floor.
void tfloor(apop_model *dce){
    if (is_t(dce)) dce->parameters->vector->data[2] = floor(dce->parameters->vector->data[2]);
}

int estimate_model(apop_data *data, apop_model *dist, char *method, apop_data *true_params){
    double *starting_pt;
    if(is_bernie(dist))
        starting_pt = (double[]){.5};
    else starting_pt = (double[]) {1.6, 1.4, 10};

    Apop_settings_add_group(dist, apop_mle, 
        .starting_pt = starting_pt,
        .method       = method, .verbose   =0,
        .step_size    = 1e-1,
        .tolerance    = 1e-4,   .k         = 1.8,
        .t_initial    = 1,      .t_min     = .5
        );
    //Apop_model_add_group(dist, apop_parts_wanted);

    if((is_bernie(dist) || is_beta(dist))
       && !strcasecmp(method, "Newton hybrid"))
        return 0;
    apop_model *e = apop_estimate(data, dist);
    tfloor(e);
    Diff(0.0, apop_vector_distance(apop_data_pack(true_params), apop_data_pack(e->parameters)), 1e-1); 
    //if (is_poisson(dist)) Apop_settings_add(dist, apop_parts_wanted, covariance, 'y');
    Print_dot
    e = apop_estimate_restart(e);
    tfloor(e);
    Diff(0.0, apop_vector_distance(apop_data_pack(true_params),apop_data_pack(e->parameters)), 1e-1); 

        if (!strcmp(e->name, "Dirichlet distribution")
            || !strcmp(e->name, "Gamma distribution") //just doesn't work.
            ||(is_bernie(e) && !strcasecmp(method, "Newton hybrid"))
            ||(is_t(e)) //requires several restarts to work.
            ||(!strcmp(e->name, "Exponential distribution")) //imprecise
            || !strcmp(e->name, "Yule distribution")){
            //cycle takes all day.
            return 0;
        }

    apop_model *dc = apop_model_copy(dist);
    Apop_settings_add(dc, apop_mle, tolerance, 1e-4);
    Apop_settings_add(dc, apop_mle, dim_cycle_tolerance, fabs(apop_log_likelihood(data, e))/200.); //within .5%.
    Print_dot
    apop_model *dce = apop_estimate(data, dc);
    Print_dot
    Diff(0.0, apop_vector_distance(apop_data_pack(true_params),apop_data_pack(dce->parameters)), 1e-2); 
    return 0;
}

/*Produce random data, then try to recover the original params */
void test_one_distribution(gsl_rng *r, apop_model *model, apop_model *true_params){
    long int runsize = 1e5;
    //generate.
    apop_data *data = apop_data_calloc(runsize, model->dsize);
    if (!strcmp(model->name, "Wishart distribution")){
        data = apop_data_calloc(runsize,4);
        true_params->parameters->vector->data[0] = runsize-4;
        for (size_t i=0; i< runsize; i++){
            Apop_row_v(data, i, v)
            true_params->draw(v->data, r, true_params);
            assert(!isnan(apop_sum(v)));
        }
    } else {
        for (size_t i=0; i< runsize; i++){
            Apop_row_v(data, i, v)
            true_params->draw(v->data, r, true_params);
            assert(!isnan(apop_sum(v)));
        }
    }
    if (model->estimate) estimate_model(data, model, "", true_params->parameters);
    else { //try all the MLEs.
        estimate_model(data, model, "NM simplex", true_params->parameters);
        if(is_t(model)) return; //t distribution still v. slow to converge.
        estimate_model(data, model, "PR cg", true_params->parameters);
        estimate_model(data, model, "Newton Hybrid", true_params->parameters);
    }
    apop_data_free(data);
}

void test_cdf(gsl_rng *r, apop_model *m){//m is parameterized
    //Make random draws from the dist, then find the CDF at that draw
    //That should generate a uniform distribution.
    if (!m->cdf || is_bernie(m) || is_binom(m))
        return;
    int drawct = 1e4;
    apop_data *draws = apop_data_alloc(drawct, m->dsize);
    apop_data *cdfs = apop_data_alloc(drawct);
    for (int i=0; i< drawct; i++){
        Apop_row(draws, i, onerow);
        Apop_stopif(apop_draw(onerow->matrix->data, r, m), abort(), 0, "bad draw.");
        Apop_row(draws, i, one_data_pt);
        apop_data_set(cdfs, i, -1, apop_cdf(one_data_pt, m));
    }
    apop_model *cdf = apop_estimate(apop_data_sort(cdfs), apop_pmf);
    apop_model *u01 = apop_model_set_parameters(apop_uniform, 0, 1);
    apop_data *ktest = apop_test_kolmogorov(cdf, u01);
    //apop_data_show(ktest);
    double maxdist = apop_data_get(ktest, .rowname="max distance");
    assert(maxdist < .03); //the K-S test has high confidence of rejection with large N
    apop_data_free(ktest); apop_data_free(draws); apop_data_free(cdfs);
    apop_model_free(u01);  apop_model_free(cdf);
}

double true_parameter_v[] = {1.82,2.1};

void test_distributions(gsl_rng *r){
    if (verbose) printf("\n");
    apop_model* true_params;
    apop_model *null_model = &(apop_model){"the null model"};

#define model_no_est(base) \
    apop_model * base ## _no_est = apop_model_copy(apop_##base);\
    base ## _no_est->estimate=NULL;

    model_no_est(beta);
    model_no_est(bernoulli);
    model_no_est(gamma);
    model_no_est(exponential);
    model_no_est(poisson);

    apop_t_distribution->estimate=NULL; //find df by MLE, not observation count.
    apop_model *dist[] = {
            apop_bernoulli, bernoulli_no_est, 
            apop_beta, beta_no_est,
            apop_binomial, apop_dirichlet,
            apop_exponential, exponential_no_est,
            apop_gamma, gamma_no_est,
            apop_lognormal, apop_multinomial, 
            apop_multivariate_normal,
            apop_normal, apop_poisson, poisson_no_est,
            apop_t_distribution, apop_uniform,
            apop_yule, apop_zipf, /*apop_wishart,*/
            null_model};

    for (int i=0; strcmp(dist[i]->name, "the null model"); i++){
        if (verbose) {printf("%s: ", dist[i]->name); fflush(NULL);}
        true_params = apop_model_copy(dist[i]);
        true_params->parameters = apop_data_fill_base(apop_data_alloc(dist[i]->vsize==1 ? 1 : 2), true_parameter_v);
        if (is_beta(dist[i]))
            true_params->parameters = apop_data_falloc((2), .5, .2);
        if (is_bernie(dist[i]))
            true_params->parameters = apop_data_falloc((1), .1);
        if (is_binom(dist[i])){
            true_params->parameters = apop_data_falloc((2), 15, .2);
            dist[i]->dsize=2;
        }
        if (!strcmp(dist[i]->name, "Dirichlet distribution"))
            dist[i]->dsize=2;
        if (!strcmp(dist[i]->name, "Multivariate normal distribution")){
            true_params->parameters = apop_data_falloc((2, 2, 2), 15, .5, .2,
                                                                   3, .2, .5);
            dist[i]->dsize=2;
        }
        if (!strcmp(dist[i]->name, "Multinomial distribution")){
            true_params->parameters = apop_data_falloc((4), 15, .5, .2, .1);
            dist[i]->dsize=4;
        }
        if (apop_regex(dist[i]->name, "gamma distribution"))
            true_params->parameters = apop_data_falloc((2), 1.5, 2.5);
        if (is_t(dist[i]))
            true_params->parameters = apop_data_falloc((3), 1, 3, 16);
        if (!strcmp(dist[i]->name, "Wishart distribution")){
            true_params->parameters = apop_data_falloc((2, 2, 2), 996, .2, .1,
                                                                    0, .1, .2);
            apop_vector_realloc(true_params->parameters->vector, 1);
        }
        test_one_distribution(r, dist[i], true_params);
        test_cdf(r, true_params);
        if (verbose) {printf("\nPASS.   "); fflush(NULL);}
    }
}
static void got_bored(){ exit(0); }

int main(int argc, char **argv){
    apop_opts.thread_count = 2;
    char c, opts[] = "sqt:";
    if (argc==1)
        printf("\tDistribution tests. Each dot is an optimization run, including some methods known to be inefficient.\n\tFor quieter output, use -q. Default is two threads; change with -t1, -t3, ...\n");
    while((c = getopt(argc, argv, opts))!=-1)
        if (c == 'q')      verbose  --;
        else if (c == 't') apop_opts.thread_count  = atoi(optarg);

    gsl_rng *r = apop_rng_alloc(213452);
    signal(SIGINT, got_bored);
    test_distributions(r);
}
