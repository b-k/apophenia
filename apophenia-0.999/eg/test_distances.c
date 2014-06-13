#include <apop.h>

/* Test distance calculations using a 3-4-5 triangle */
int main(){
    gsl_vector *v1 = gsl_vector_alloc(2);
    gsl_vector *v2 = gsl_vector_alloc(2);
    gsl_vector_set(v1, 0,2);
    gsl_vector_set(v1, 1,2);
    gsl_vector_set(v2, 0,5);
    gsl_vector_set(v2, 1,6);
    assert(apop_vector_distance(v1, v1, 'd') == 0);     //discrete: if vectors are equal d==0;
    assert(apop_vector_distance(v1, v2, 'd') == 1);     //          if vectors differ d ==1
    assert(apop_vector_distance(v1, NULL, 'm') == 4.);  //length of v1, Manhattan metric
    assert(apop_vector_distance(v1,v2) == 5.);          //the hypotenuse of the 3-4-5 triangle
    assert(apop_vector_distance(v2, NULL, 's') == 6);   //length of v2, sup norm
    assert(apop_vector_distance(v1,v2, 'm') == 7.);     //distance v1 to v2 in Manhattan
    assert(apop_vector_distance(v1,v2, 'L', 2) == 5.);  //L_2 norm == Euclidean
}
