#include <apop.h>

int main(void){
gsl_vector  *in, *out;

in = gsl_vector_calloc(3);
apop_vector_fill(in, 0, 1, 2);

printf("The original vector:\n");
apop_vector_print(in);

apop_vector_normalize(in, &out, 's');
printf("Standardized with mean zero and variance one:\n");
apop_vector_print(out);
assert(apop_vector_sum(out)<1e-5);
//assert(fabs((apop_vector_var(out))- 1)<1e-5);

apop_vector_normalize(in, &out, 'r');
printf("Normalized range with max one and min zero:\n");
apop_vector_print(out);
assert(gsl_vector_max(out)==1);
assert(gsl_vector_min(out)==0);

apop_vector_normalize(in, NULL, 'p');
printf("Normalized into percentages:\n");
apop_vector_print(in);
}
