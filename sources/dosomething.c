#include <stdio.h>           //printf
float do_something (float a, float b, int c) {
    float d;
    d = b + c;
    return a/d;
}

int main(void){
    float output;
    float first = 7.3;
    float second = 3.3;
    int third = 2;
    output = do_something(first, second, third);
    printf("the output: %g\n", output);
    return 0;
}
