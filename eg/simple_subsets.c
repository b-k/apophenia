#include <apop.h>

#ifdef Datadir
#define DATAFILE Datadir "/" "data"
#else
#define DATAFILE "data"
#endif

int main(){
    apop_table_exists( DATAFILE , 'd');
    apop_data *d = apop_text_to_data( DATAFILE );
  
    //tally row zero of the data set's matrix by viewing it as a vector:
    gsl_vector *one_row = Apop_rv(d, 0);
    double sigma = apop_vector_sum(one_row);
    printf("Sum of row zero: %g\n", sigma);
    assert(sigma==14);
  
    //view column zero as a vector; take its mean
    double mu = apop_vector_mean(Apop_cv(d, 0));
    printf("Mean of col zero: %g\n", mu);
    assert(fabs(mu - 19./6)<1e-5);
  
    //get a sub-data set (with names) of two rows beginning at row 3; print to screen
    apop_data *six_elmts = Apop_rs(d, 3, 2);
    apop_data_print(six_elmts);
}
