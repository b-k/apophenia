#include <apop.h>

int main(){
    apop_table_exists("data", 'd');
    apop_data *d = apop_text_to_data("data");
  
    //tally row zero of the data set's matrix by viewing it as a vector:
    Apop_row_v(d, 0, one_row);
    double sigma = apop_vector_sum(one_row);
    printf("Sum of the first row: %g\n", sigma);
    assert(sigma==14);
  
    //view the first column as a vector; take its mean
    Apop_col_v(d, 0, one_col);
    double mu = apop_vector_mean(one_col);
    printf("Mean of the first col: %g\n", mu);
    assert(fabs(mu - 19./6)<1e-5);
  
    //get a sub-data set (with names) of two rows beginning at row 3; print to screen
    Apop_rows(d, 3, 2, six_elmts);
    apop_data_print(six_elmts);
}
