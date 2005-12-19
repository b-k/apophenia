/** \file apop_name.c

Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
*/

/** Allocates a name structure
\return	An allocated, empty name structure.
\ingroup names
*/

#include "db.h" //just for apop_verbose.
#include "types.h"
#include <stdio.h>
#include <malloc.h>

apop_name * apop_name_alloc(void){
apop_name	* init_me;
	init_me	= malloc(sizeof(apop_name));
	init_me->colnames	= malloc(1);
	init_me->rownames	= malloc(1);
	init_me->depnames	= malloc(1);
	init_me->colnamect	= 
	init_me->rownamect	=
	init_me->depnamect	= 0;
	return init_me;
}

/** Adds a name to the \ref apop_name structure. Puts it at the end of the given list.

\param n 	An existing, allocated \ref apop_name structure.
\param add_me 	A string.
\param type 	If adding a dependent variable, use <tt>'d'</tt>; if adding a row name, use <tt>'r'</tt>;
If adding a (independent) column name, use <tt>'c'</tt>.
\return 	Returns the number of rows/cols/depvars after you have added the new one.
\ingroup names
*/
int apop_name_add(apop_name * n, char *add_me, char type){
	if (type == 'c'){
		(n->colnamect)++;
		n->colnames	= realloc(n->colnames, sizeof(char*) * n->colnamect);
		n->colnames[n->colnamect -1]	= malloc(sizeof(char) * (strlen(add_me) + 1));
		strcpy(n->colnames[n->colnamect -1], add_me);
		return n->colnamect;
	} 
	if (type == 'r'){
		(n->rownamect)++;
		n->rownames	= realloc(n->rownames, sizeof(char*) * n->rownamect);
		n->rownames[n->rownamect -1]	= malloc(sizeof(char) * (strlen(add_me) + 1));
		strcpy(n->rownames[n->rownamect -1], add_me);
		return n->rownamect;
	} //else:  type == 'd'
		(n->depnamect)++;
		n->depnames	= realloc(n->depnames, sizeof(char*) * n->depnamect);
		n->depnames[n->depnamect -1]	= malloc(sizeof(char) * (strlen(add_me) + 1));
		strcpy(n->depnames[n->depnamect -1], add_me);
		return n->depnamect;
}

/** Prints the given list of names to STDOUT
\param n	the \ref apop_name structure
\ingroup names
*/
void  apop_name_print(apop_name * n){
int		i;
	if (n->depnamect > 0){
		printf("\t\t\t");
		for (i=0; i < n->depnamect; i++)
			printf("\t%s", n->depnames[i]);
		printf("\n");
	}
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
	
/** Erases an \ref apop_name structure.
\ingroup names 	*/
void  apop_name_free(apop_name * free_me){
int		i;
	for (i=0; i < free_me->colnamect; i++)
		free(free_me->colnames[i]);
	for (i=0; i < free_me->rownamect; i++)
		free(free_me->rownames[i]);
	for (i=0; i < free_me->depnamect; i++)
		free(free_me->depnames[i]);
	free(free_me->colnames);
	free(free_me->rownames);
	free(free_me->depnames);
	free(free_me);
}

/** Append one list of names to another.

\param  n1  The first set of names
\param  n2  The second set of names, which will be appended after the first.
\params type    Either 'd', 'c', or 'r', stating whether you are merging the dependent var names, columns, or rows. [Default: cols]
\ingroup names */
void  apop_name_stack(apop_name * n1, apop_name *n2, char type){
size_t  n1ct, n2ct;
char **n2names;

    if (type == 'r'){
        n1ct    = n1->rownamect;
        n2ct    = n2->rownamect;
        n2names = n2->rownames;
        n1->rownames     = realloc(n1->rownames, sizeof(char *)*(n1ct+n2ct));
        memcpy((n1->rownames)+n1ct, n2names, sizeof(char *)*n2ct);
        n1->rownamect   += n2ct;
        }
    if (type == 'd'){
        n1ct    = n1->depnamect;
        n2ct    = n2->depnamect;
        n2names = n2->depnames;
        n1->depnames     = realloc(n1->depnames, sizeof(char *)*(n1ct+n2ct));
        memcpy((n1->depnames)+n1ct, n2names, sizeof(char *)*n2ct);
        n1->depnamect   += n2ct;
        }
    else {
        if (type != 'c' && apop_verbose)
            printf ("You gave me >%c<, I'm assuming you meant c; copying column names.\n",type);
        n1ct    = n1->colnamect;
        n2ct    = n2->colnamect;
        n2names = n2->colnames;
        n1->colnames     = realloc(n1->colnames, sizeof(char *)*(n1ct+n2ct));
        memcpy((n1->colnames)+n1ct, n2names, sizeof(char *)*n2ct);
        n1->colnamect   += n2ct;
        }
}

/** Remove the columns set to one in the \c drop vector.
\param n the \ref apop_name structure to be pared down
\param drop  a vector with n->colnamect elements, mostly zero, with a one marking those columns to be removed.
 */
void apop_name_rm_columns(apop_name *n, int *drop){
apop_name   *newname    = apop_name_alloc();
int         i;
    for (i=0; i< n->colnamect; i++){
        if (drop[i]==0)
            apop_name_add(newname, n->colnames[i],'c');
        else
            n->colnamect    --;
    }
    free(n->colnames);
    n->colnames = newname->colnames;
    //we need to free the newname struct, but leave the colnames intact.
    newname->colnames   = malloc(1);
    newname->colnamect  = 0;
    apop_name_free(newname);
}
