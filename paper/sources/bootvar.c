#include "matrix_convenience_fns.h"  //variance

//I assume calc_statistic() has been defined. mean() will work.

double find_boot_variance(gsl_vector *data,
                     int samples_to_take, int size_of_samples){
    gsl_vector *artificial_statistics = gsl_vector_alloc(samples_to_take);
    gsl_vector *our_sample = gsl_vector_alloc(draws);
    double return_me;
    int i;
    for (i = 0; i < samples_to_take; i++){
         boot_draw(data, size_of_samples, our_sample);
         gsl_vector_set(artificial_statistics, i, calc_statistic(our_sample));
    }
    return_me = std_dev(artificial_statistics);
    gsl_vector_free(artificial_statistics);
    gsl_vector_free(our_sample);
    return return_me;
}
