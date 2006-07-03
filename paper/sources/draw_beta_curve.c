#include <apophenia/headers.h>

int main(){
int				i;
int				draws	= 1e6;
int				bins	= 100;
double		mu		= 0.492;    //also try 3./8.
double		sigma	= 0.093;    //also try 1./24.
char			outfile[]= "betaplot";
gsl_rng			*r		= apop_rng_alloc(0);
gsl_vector	*v		= gsl_vector_alloc(draws);
    for (i =0; i < draws; i++)
        gsl_vector_set(v, i, apop_random_beta(r, mu, sigma));
    remove(outfile);    //remove file to prevent overwriting.
    apop_plot_histogram(v, bins, outfile);
    return 0;
}
