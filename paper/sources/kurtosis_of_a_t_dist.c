#include <apophenia/headers.h>
int main(){
int         i, runct    = 5000000;
int         df, df_max  = 31;
gsl_vector 	*v  = gsl_vector_alloc(runct);
gsl_rng     *r  = apop_rng_alloc(0);
    printf("df\t k (est)\t k (actual)\n");
    for (df=1; df< df_max; df++){
        for (i=0; i< runct; i++)
            gsl_vector_set(v, i, gsl_ran_tdist(r, df));
        printf("%i\t %g", df, apop_vector_kurtosis(v)*(runct-1)/runct);
        if (df > 4)
            printf("\t%g", (3.*df - 6.)/(df-4.));
        printf("\n");
    }
	return 0;
}
