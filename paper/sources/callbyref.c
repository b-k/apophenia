#include <stdio.h>  //printf
#include <malloc.h> //malloc

int doubling (int * a_c, int b_c){
	*a_c = b_c * 2;
	return *a_c;
}

int main(void){
	int *a = malloc(sizeof(int));
	int b = 2;
	printf("doubling() returns: %i\n", doubling(a,b));
	printf("a now holds: %i\n", *a);
	free(a);
}
