#include <apophenia/headers.h>

int             agentct     = 5000;
int             periods     = 50;
int             binct       = 30;
double          pauselength = 0.6;
char            outfile[]   = "hist";

void initialize(gsl_vector *agentlist, gsl_rng *r){
  int   i;
    for(i=0; i<agentct; i++){
        gsl_vector_set(agentlist, i, gsl_rng_uniform(r)*100);
    } 
}

void grow(gsl_vector *agentlist, gsl_rng *r){
  int       i;
  double    g;
    for(i=0; i<agentct; i++){
        g   = gsl_ran_gaussian(r,0.1);
        gsl_vector_set(agentlist,i, 
                gsl_vector_get(agentlist,i)* exp(g));
    }
}

double estimate(gsl_vector *agentlist){
    return apop_vector_mean(agentlist);
}

int main(){
  gsl_vector    *agentlist  = gsl_vector_alloc(agentct);
  gsl_rng       *r          = apop_rng_alloc(39);
  int           i;
  FILE          *f;
    initialize(agentlist, r);
    remove(outfile);
    apop_plot_histogram(agentlist, binct, outfile);
    for(i=0; i< periods; i++){
        grow(agentlist, r);
        f   = fopen(outfile, "a");
        fprintf(f, "e\npause %g\n", pauselength);
        fclose(f);
        apop_plot_histogram(agentlist, binct, outfile);
    } 
    printf("the mean: %g\n",estimate(agentlist));
    return 0;
}
