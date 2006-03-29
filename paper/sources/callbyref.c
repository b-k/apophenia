#include <stdio.h>  //printf
#include <malloc.h> //malloc

int doubling (int * k_c, int b_c){
	*k_c = b_c * 2;
	return *k_c;
}

int main(void){
	int *k = malloc(sizeof(int));
	int b = 2;
	printf("doubling() returns: %i\n", doubling(k,b));
	printf("a now holds: %i\n", *k);
	free(k);
}
