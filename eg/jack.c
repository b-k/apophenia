#include <apop.h>

int main(){
    int draw_ct = 1000;
    apop_model *m = apop_model_set_parameters(apop_normal, 1, 3);
    double sigma = apop_data_get(m->parameters, 1);
    apop_data *d = apop_model_draws(m, draw_ct);
    apop_data *out = apop_jackknife_cov(d, m);
    double error = fabs(apop_data_get(out, 0,0)-gsl_pow_2(sigma)/draw_ct) //var(mu)
                + fabs(apop_data_get(out, 1,1)-gsl_pow_2(sigma)/(2*draw_ct))//var(sigma)
                +fabs(apop_data_get(out, 0,1)) +fabs(apop_data_get(out, 1,0));//cov(mu,sigma); should be 0.
    apop_data_free(d);
    apop_data_free(out);
    assert(error < 1e-2);//Not very accurate.
}


