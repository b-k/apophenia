#include <apop.h>

double row_offset;

void offset_rng(double *v){*v = gsl_rng_uniform(apop_rng_get_thread()) + row_offset;}
double find_tstat(gsl_vector *in){ return apop_mean(in)/sqrt(apop_var(in));}
double conf(double in, void *df){ return gsl_cdf_tdist_P(in, *(int *)df);}

//apop_vector_mean is a macro, so we can't point a pointer to it.
double mu(gsl_vector *in){ return apop_vector_mean(in);}

int main(){
    apop_data *d = apop_data_alloc(10, 100);
    gsl_rng *r = apop_rng_alloc(3242);
    for (int i=0; i< 10; i++){
        row_offset = gsl_rng_uniform(r)*2 -1; //declared and used above.
        apop_vector_apply(Apop_rv(d, i), offset_rng);
    }

    int df = d->matrix->size2-1;
    apop_data *means = apop_map(d, .fn_v = mu, .part ='r');
    apop_data *tstats = apop_map(d, .fn_v = find_tstat, .part ='r');
    apop_data *confidences = apop_map(tstats, .fn_dp = conf, .param = &df);

    printf("means:\n"); apop_data_show(means);
    printf("\nt stats:\n"); apop_data_show(tstats);
    printf("\nconfidences:\n"); apop_data_show(confidences);

    //Some sanity checks, for Apophenia's test suite.
    for (int i=0; i< 10; i++){
        //sign of mean == sign of t stat.
        assert(apop_data_get(means, i, -1) * apop_data_get(tstats, i, -1) >=0);

        //inverse of P-value should be the t statistic.
        assert(fabs(gsl_cdf_tdist_Pinv(apop_data_get(confidences, i, -1), 99) 
                    - apop_data_get(tstats, i, -1)) < 1e-5);
    }
}
