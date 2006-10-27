#include <stdlib.h> //environment handling
#include <stdio.h>

int main(){
    printf("You are: %s\n", getenv("USER"));
    setenv("USER", "Steven", 1);
    printf("But now you are: %s\n", getenv("USER"));
    return 0;
}
