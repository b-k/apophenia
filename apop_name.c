/** \file apop_name.c

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "db.h" //just for apop_opts.verbose.
#include "types.h"
#include <stdio.h>
#include <regex.h>

/** \defgroup names apop_names: a structure and functions to handle column names for matrices. 
\ingroup types */

/** Allocates a name structure
\return	An allocated, empty name structure.
\ingroup names
*/
apop_name * apop_name_alloc(void){
apop_name	* init_me;
	init_me	= malloc(sizeof(apop_name));
	init_me->vecname	= NULL;
	init_me->colnames	= NULL;
	init_me->textnames	= NULL;
	init_me->rownames	= NULL;
	init_me->colnamect	= 
	init_me->textnamect	= 
	init_me->rownamect	= 0;
	init_me->title[0]   = '\0';
	return init_me;
}

/** Adds a name to the \ref apop_name structure. Puts it at the end of the given list.

\param n 	An existing, allocated \ref apop_name structure.
\param add_me 	A string. If NULL, do nothing; return -1.
\param type 	'r': add a row name<br>
'c': add a column name<br>
't': add a text category name<br>
'h': add a title (or a header. 't' is taken).<br>
'v': add (or overwrite) the vector name<br>
\return 	Returns the number of rows/cols/depvars after you have added the new one.
\ingroup names
*/
int apop_name_add(apop_name * n, char *add_me, char type){
    if (!add_me)
        return -1;
	if (type == 'h'){
        if (add_me){
            snprintf(n->title, 100, add_me);
		    return 1;
        } else return 0;
	} 
	if (type == 'v'){
        if (add_me){
		n->vecname	= realloc(n->vecname,  strlen(add_me) + 1);
		strcpy(n->vecname, add_me);
		return 1;
        } else return 0;
	} 
	if (type == 'r'){
		(n->rownamect)++;
		n->rownames	= realloc(n->rownames, sizeof(char*) * n->rownamect);
		n->rownames[n->rownamect -1]	= malloc(strlen(add_me) + 1);
		strcpy(n->rownames[n->rownamect -1], add_me);
		return n->rownamect;
	} 
	if (type == 't'){
		(n->textnamect)++;
		n->textnames	= realloc(n->textnames, sizeof(char*) * n->textnamect);
		n->textnames[n->textnamect -1]	= malloc(strlen(add_me) + 1);
		strcpy(n->textnames[n->textnamect -1], add_me);
		return n->textnamect;
	}
	//else assume (type == 'c'){
        if (type != 'c' && apop_opts.verbose)
            apop_error(2,'c',"You gave me >%c<, I'm assuming you meant c; copying column names.\n",type);
		(n->colnamect)++;
		n->colnames	= realloc(n->colnames, sizeof(char*) * n->colnamect);
		n->colnames[n->colnamect -1]	= malloc(strlen(add_me) + 1);
		strcpy(n->colnames[n->colnamect -1], add_me);
		return n->colnamect;
	//} 
}

/** Prints the given list of names to STDOUT
\param n	the \ref apop_name structure
\ingroup names
*/
void  apop_name_print(apop_name * n){
int		i;
	if (n->vecname){
		printf("\t\t\t");
			printf("\t%s", n->vecname);
		printf("\n");
	}
	if (n->colnamect > 0){
		printf("\t\t\t");
		for (i=0; i < n->colnamect; i++)
			printf("\t%s", n->colnames[i]);
		printf("\n");
	}
	if (n->textnamect > 0){
		printf("\t\t\t");
		for (i=0; i < n->textnamect; i++)
			printf("\t%s", n->textnames[i]);
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
	for (i=0; i < free_me->textnamect; i++)
		free(free_me->textnames[i]);
	for (i=0; i < free_me->rownamect; i++)
		free(free_me->rownames[i]);
    if (free_me->vecname);
        free(free_me->vecname);
	free(free_me->colnames);
	free(free_me->textnames);
	free(free_me->rownames);
	free(free_me);
}

/** Append one list of names to another.

Notice that if the first list is NULL, then this is a copy function. If the second is NULL, it is a no-op.

If you are copying row names to columns or vice versa, use \ref apop_name_cross_stack.

\param  n1      The first set of names
\param  n2      The second set of names, which will be appended after the first.
\param type     Either 'c', 'r', 't', or 'v' stating whether you are merging the columns, rows, or text. If 'v', then overwrite the target with the source vector name.
\ingroup names */
void  apop_name_stack(apop_name * n1, apop_name *n2, char type){
int     i;
    if (!n2)
        return;
    if (type == 'v'){
        apop_name_add(n1, n2->vecname, 'v');
        return;
    }
    if (type == 'r'){
        for (i=0; i< n2->rownamect; i++)
            apop_name_add(n1, n2->rownames[i], 'r');
        return;
    }
    if (type == 't'){
        for (i=0; i< n2->textnamect; i++)
            apop_name_add(n1, n2->textnames[i], 't');
        return;
    }
    if (type == 'c'){
        for (i=0; i< n2->colnamect; i++)
            apop_name_add(n1, n2->colnames[i], 'c');
        return;
    }
    if (apop_opts.verbose)
         printf (">%c< sent to apop_name_stack, but the only valid options are r t c. Doing nothing.\n",type);
}

/** Append one list of names to another; the source and dest list need
not be the same type. If they are, just use \ref apop_name_stack.



Notice that if the first list in NULL, then this is a copy function.

\param  n1      The first set of names
\param  n2      The second set of names, which will be appended after the first.
\param type1     Either 'c', 'r', or 't' stating whether you are merging from the columns, rows, or text. [Default: cols]
\param type2     Either 'c', 'r', or 't' stating whether you are merging to the columns, rows, or text. [Default: cols]
\ingroup names */
void  apop_name_cross_stack(apop_name * n1, apop_name *n2, char type1, char type2){
int     i;
    if (type1 == 'r'){
        for (i=0; i< n2->rownamect; i++)
            apop_name_add(n1, n2->rownames[i], type2);
        }
    else if (type1 == 't'){
        for (i=0; i< n2->textnamect; i++)
            apop_name_add(n1, n2->textnames[i], type2);
        }
    else {
        if (type1 != 'c' && apop_opts.verbose)
            printf ("You gave me >%c<, I'm assuming you meant c; copying column names.\n",type1);
        for (i=0; i< n2->colnamect; i++)
            apop_name_add(n1, n2->colnames[i], type2);
        }
}

/** Copy one \ref apop_name structure to another. That is, all data is duplicated. Usage: 

Much like \ref apop_name_memcpy, but a slightly different calling form. Usage:
\code
apop_name *out;
apop_name_memcpy(&out, in);
\endcode
 
  \param out    a structure that this function will allocate and fill
  \param in    the input names

 \ingroup names
  */
void apop_name_memcpy(apop_name **out, apop_name *in){
    *out = apop_name_alloc();
    apop_name_stack(*out, in, 'v');
    apop_name_stack(*out, in, 'c');
    apop_name_stack(*out, in, 'r');
    apop_name_stack(*out, in, 't');
    snprintf((*out)->title, 100, in->title);
}

/** Copy one \ref apop_name structure to another. That is, all data is duplicated. 

Much like \ref apop_name_memcpy, but a slightly different calling form. Usage:

\code
apop_name *out  = apop_name_copy(in);
\endcode
 
    \param in    the input names
    \return       a structure that this function will allocate and fill
    \ingroup names
  */
apop_name * apop_name_copy(apop_name *in){
apop_name *out = apop_name_alloc();
    apop_name_stack(out, in, 'v');
    apop_name_stack(out, in, 'c');
    apop_name_stack(out, in, 'r');
    apop_name_stack(out, in, 't');
    snprintf(out->title, 100, in->title);
    return out;
}


/** Remove the columns set to one in the \c drop vector.
\param n the \ref apop_name structure to be pared down
\param drop  a vector with n->colnamect elements, mostly zero, with a one marking those columns to be removed.
\ingroup names
 */
void apop_name_rm_columns(apop_name *n, int *drop){
apop_name   *newname    = apop_name_alloc();
int         i, max      = n->colnamect;
    for (i=0; i< max; i++){
        if (drop[i]==0)
            apop_name_add(newname, n->colnames[i],'c');
        else
            n->colnamect    --;
    }
    free(n->colnames);
    n->colnames = newname->colnames;
    //we need to free the newname struct, but leave the colnames intact.
    //A one-byte memory leak.
    newname->colnames   = malloc(1);
    newname->colnamect  = 0;
    apop_name_free(newname);
}



/** Finds the position of an element in a list of names.

The function uses case-insensitive regular expressions to search. 

For example, "p.val.*" will match "P value", "p.value", and "p values".

\param n        the \ref apop_name object to search.
\param in       the name you seek; see above.
\param type     'c', 'r', or 't'.
\return         The position of \c findme. If not found, returns -1.
\ingroup names
  */
size_t  apop_name_find(apop_name *n, char *in, char type){
  regex_t   re;
  char      **list;
  int       i, listct;
    if (type == 'r'){
        list    = n->rownames;
        listct  = n->rownamect;
    }
    else if (type == 't'){
        list    = n->textnames;
        listct  = n->textnamect;
    }
    else { // default: (type == 'c')
        list    = n->colnames;
        listct  = n->colnamect;
    }
    regcomp(&re, in, REG_ICASE);
    for (i = 0; i < listct; i++){
        if (!regexec(&re, list[i], 0, NULL, 0)){
            regfree(&re);
            return i;
        }
    }
    regfree(&re);
    return -1;
}
