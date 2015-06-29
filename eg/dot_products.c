/* A demonstration of dot products and various useful 
   transformations among types. */

#include <apop.h>

double eps=1e-3;//slow to converge series-->large tolerance.
#define Diff(L, R) Apop_assert(fabs((L)-(R)<(eps)), "%g is too different from %g (abitrary limit=%g).", (double)(L), (double)(R), eps);

int main(){
    int len = 3000;
    gsl_vector *v = gsl_vector_alloc(len);
    for (double i=0; i< len; i++) gsl_vector_set(v, i, 1./(i+1));
    double square;
    gsl_blas_ddot(v, v, &square);
    printf("1 + (1/2)^2 + (1/3)^2 + ...= %g\n", square);

    double pi_over_six = gsl_pow_2(M_PI)/6.;
    Diff(square, pi_over_six);

    /* Now using apop_dot, in a few forms.
       First, vector-as-data dot itself.
       If one of the inputs is a vector,
       apop_dot puts the output in a vector-as-data:*/
    apop_data *v_as_data = &(apop_data){.vector=v};
    apop_data *vdotv = apop_dot(v_as_data, v_as_data);
    Diff(gsl_vector_get(vdotv->vector, 0), pi_over_six);

    /* Wrap matrix in an apop_data set. */
    gsl_matrix *v_as_matrix = apop_vector_to_matrix(v);
    apop_data dm = (apop_data){.matrix=v_as_matrix};

    // (1 X len) vector dot (len X 1) matrix --- produce a scalar (one item vector).
    apop_data *mdotv = apop_dot(v_as_data, &dm);
    double scalarval = apop_data_get(mdotv);
    Diff(scalarval, pi_over_six);

    //(len X 1) dot (len X 1) --- bad dimensions.
    apop_opts.verbose=-1; //don't print an error.
    apop_data *mdotv2 = apop_dot(&dm, v_as_data);
    apop_opts.verbose=0; //back to safety.
    assert(mdotv2->error);

    // If we want (len X 1) dot (1 X len) --> (len X len),
    // use apop_vector_to_matrix.
    apop_data dmr = (apop_data){.matrix=apop_vector_to_matrix(v, .row_col='r')};
    apop_data *product_matrix = apop_dot(&dm, &dmr);
    //The trace is the sum of squares:
    gsl_vector_view trace = gsl_matrix_diagonal(product_matrix->matrix);
    double tracesum = apop_sum(&trace.vector);
    Diff(tracesum, pi_over_six);

    apop_data_free(product_matrix);
    gsl_matrix_free(dmr.matrix);
}
