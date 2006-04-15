#include <stdio.h>  //printf
int globe=1;          //a global variable.

int doubling (int a_c, int b_c){
     a_c = b_c * 2;
     globe = b_c * 2;
     return a_c;
}

int main(void){
    int a = 1;
    int b = 2;
    printf("doubling returns %i\n", doubling(a,b));
    printf("a= %i\n", a);
    printf("globe= %i\n", globe);
    return 0;
}
