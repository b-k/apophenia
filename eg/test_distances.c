#include <apop.h>

/* Test distance calculations using a 3-4-5 triangle */
int main(){
    gsl_vector *v1 = gsl_vector_alloc(2);
    gsl_vector *v2 = gsl_vector_alloc(2);
    apop_vector_fill(v1, 2, 2);
    apop_vector_fill(v2, 5, 6);

    assert(apop_vector_distance(v1, v1, 'd') == 0);
    assert(apop_vector_distance(v1, v2, 'd') == 1);
    assert(apop_vector_distance(v1, .metric='m') == 4);
    assert(apop_vector_distance(v2, .metric='s') == 6);
    assert(apop_vector_distance(v1,v2) == 5.); //the hypotenuse of the 3-4-5 triangle
    assert(apop_vector_distance(v1,v2, 'm') == 7.);
    assert(apop_vector_distance(v1,v2, 'L', 2) == 5.);  //L_2 norm == Euclidean
}
