#include <apop.h>

int main(){
    apop_table_exists("data", 'd');
    apop_data *d = apop_text_to_data("data");
  
    //tally row zero of the data set's matrix by viewing it as a vector:
    gsl_vector *one_row = Apop_rv(d, 0);
    double sigma = apop_vector_sum(one_row);
    printf("Sum of the first row: %g\n", sigma);
    assert(sigma==14);
  
    //view the first column as a vector; take its mean
    double mu = apop_vector_mean(Apop_cv(d, 0));
    printf("Mean of the first col: %g\n", mu);
    assert(fabs(mu - 19./6)<1e-5);
  
    //get a sub-data set (with names) of two rows beginning at row 3; print to screen
    apop_data *six_elmts = Apop_rs(d, 3, 2);
    apop_data_print(six_elmts);
}
