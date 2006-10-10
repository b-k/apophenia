#include <math.h>
#include <stdio.h>

int main(){
    int i;
    for(i=0; i< 1400; i+=100)
        printf("%i\t %lg \t %lg\n", i, ldexp(1,i), ldexp(1,-i));
    return 0;
}
