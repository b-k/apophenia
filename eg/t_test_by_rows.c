#include <apop.h>

double row_offset;
gsl_rng *r;

void offset_rng(double *v){*v = gsl_rng_uniform(r) + row_offset;}
double find_tstat(gsl_vector *in){ return apop_mean(in)/sqrt(apop_var(in));}
double conf(double in, void *df){ return gsl_cdf_tdist_P(in, *(int *)df);}

int main(){
    apop_data *d = apop_data_alloc(0, 10, 100);
    r = apop_rng_alloc(3242);
    for (int i=0; i< 10; i++){
        row_offset = gsl_rng_uniform(r)*2 -1;
        Apop_row(d, i, onerow);
        apop_vector_apply(onerow, offset_rng);
    }

    size_t df = d->matrix->size2-1;
    apop_data *means = apop_map(d, .fn_v = apop_vector_mean, .part ='r');
    apop_data *tstats = apop_map(d, .fn_v = find_tstat, .part ='r');
    apop_data *confidences = apop_map(tstats, .fn_dp = conf, .param = &df);

    printf("means:\t\t"); apop_data_show(means);
    printf("t stats:\t"); apop_data_show(tstats);
    printf("confidences:\t"); apop_data_show(confidences);
}
