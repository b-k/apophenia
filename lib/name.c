#include "name.h"
#include <stdio.h>
#include <malloc.h>

apop_name * apop_name_alloc(void){
apop_name	* init_me;
	init_me	= malloc(sizeof(apop_name));
	init_me->colnames	= malloc(1);
	init_me->rownames	= malloc(1);
	init_me->colnamect	= 
	init_me->rownamect	= 0;
	return init_me;
}

int apop_name_add(apop_name * n, char *add_me, int is_column){
	if (is_column){
		(n->colnamect)++;
		n->colnames	= realloc(n->colnames, sizeof(char*) * n->colnamect);
		n->colnames[n->colnamect -1]	= malloc(sizeof(char) * (strlen(add_me) + 1));
		strcpy(n->colnames[n->colnamect -1], add_me);
		return n->colnamect;
	} //else :
		(n->rownamect)++;
		n->rownames	= realloc(n->rownames, sizeof(char*) * n->rownamect);
		n->rownames[n->rownamect -1]	= malloc(sizeof(char) * (strlen(add_me) + 1));
		strcpy(n->rownames[n->rownamect -1], add_me);
		return n->rownamect;
}

void  apop_name_print(apop_name * n){
int		i;
	if (n->colnamect > 0){
		printf("\t\t\t");
		for (i=0; i < n->colnamect; i++)
			printf("\t%s", n->colnames[i]);
		printf("\n");

	}
	if (n->rownamect > 0){
		printf("\t\t\t");
		for (i=0; i < n->rownamect; i++)
			printf("\t%s", n->rownames[i]);
		printf("\n");
	}
}
	
void  apop_name_free(apop_name * free_me){
int		i;
	for (i=0; i < free_me->colnamect; i++)
		free(free_me->colnames[i]);
	for (i=0; i < free_me->rownamect; i++)
		free(free_me->rownames[i]);
	free(free_me->colnames);
	free(free_me->rownames);
	free(free_me);
}
