#include <stdio.h>  //printf
#include <malloc.h> //malloc

int doubling (int * a_c, int b_c){
	*a_c = b_c * 2;
	return *a_c;
}

int main(void){
	int *a = malloc(sizeof(int));
	int b = 2;
	printf("%i\n", doubling(a,b));
	printf("%i\n", *a);
	free(a);
}
