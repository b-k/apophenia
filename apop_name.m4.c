/** \file apop_name.c */
/* Copyright (c) 2006--2009 by Ben Klemens.  Licensed under the GPLv2; see COPYING.  */

#include "apop_internal.h"
#include <stdio.h>
#include <regex.h>

/** Allocates a name structure
\return	An allocated, empty name structure.  In the very unlikely event that \c malloc fails, return \c NULL.

Because \ref apop_data_alloc uses this to set up its output, you will rarely if ever
need to call this function explicitly. You may want to use it if wrapping a \c gsl_matrix into an \ref apop_data set. For example, to put a title on a vector:

\code
apop_data *d = &(apop_data){.vector=your_vector, .names=apop_name_alloc()};
apop_name_add(d->names, "A column of numbers", 'v');
apop_data_print(d);

...
apop_name_free(d->names); //but d itself is auto-allocated; no need to free it.
\endcode
*/
apop_name * apop_name_alloc(void){
    apop_name * init_me = malloc(sizeof(apop_name));
    Apop_stopif(!init_me, return NULL, 0, "malloc failed. Probably out of memory.");
    *init_me = (apop_name){ };
	return init_me;
}

/** Adds a name to the \ref apop_name structure. Puts it at the end of the given list.

\param n 	An existing, allocated \ref apop_name structure.
\param add_me 	A string. If \c NULL, do nothing; return -1.
\param type 	'r': add a row name<br>
'c': add a matrix column name<br>
't': add a text column name<br>
'h': add a title (i.e., a header).<br>
'v': add (or overwrite) the vector name<br>
\return 	Returns the number of rows/cols/depvars after you have added the new one. But if \c add_me is \c NULL, return -1.
*/
int apop_name_add(apop_name * n, char const *add_me, char type){
    if (!add_me)
        return -1;
	if (type == 'h'){
        free(n->title);
        Asprintf(&n->title, "%s", add_me);
        return 1;
	} 
	if (type == 'v'){
		n->vector	= realloc(n->vector,  strlen(add_me) + 1);
		strcpy(n->vector, add_me);
		return 1;
	} 
	if (type == 'r'){
		n->rowct++;
		n->row	= realloc(n->row, sizeof(char*) * n->rowct);
		n->row[n->rowct -1]	= malloc(strlen(add_me) + 1);
		strcpy(n->row[n->rowct -1], add_me);
		return n->rowct;
	} 
	if (type == 't'){
		n->textct++;
		n->text	= realloc(n->text, sizeof(char*) * n->textct);
		n->text[n->textct -1]	= malloc(strlen(add_me) + 1);
		strcpy(n->text[n->textct -1], add_me);
		return n->textct;
	}
	//else assume (type == 'c')
        Apop_stopif(type != 'c', /*keep going.*/, 
            2,"You gave me >%c<, I'm assuming you meant c; "
                             " copying column names.", type);
		n->colct++;
		n->col = realloc(n->col, sizeof(char*) * n->colct);
		n->col[n->colct -1]	= malloc(strlen(add_me) + 1);
		strcpy(n->col[n->colct -1], add_me);
		return n->colct;
}

/** Prints the given list of names to stdout. Useful for debugging.

\param n  The \ref apop_name structure
*/
void apop_name_print(apop_name * n){
    if (!n) {
        printf("NULL");
        return;
    }
	if (n->title) printf("title: %s\n", n->title);
	if (n->vector){
		printf("vector:");
        printf("\t%s\n", n->vector);
	}
	if (n->colct > 0){
		printf("column:");
		for (int i=0; i < n->colct; i++)
			printf("\t%s", n->col[i]);
		printf("\n");
	}
	if (n->textct > 0){
		printf("text:");
		for (int i=0; i < n->textct; i++)
			printf("\t%s", n->text[i]);
		printf("\n");
	}
	if (n->rowct > 0){
		printf("row:");
		for (int i=0; i < n->rowct; i++)
			printf("\t%s", n->row[i]);
		printf("\n");
	}
}
	
/** Free the memory used by an \ref apop_name structure. */
void  apop_name_free(apop_name * free_me){
    if (!free_me) return; //only needed if users are doing tricky things like newdata = (apop_data){.matrix=...};
	for (size_t i=0; i < free_me->colct; i++)  free(free_me->col[i]);
	for (size_t i=0; i < free_me->textct; i++) free(free_me->text[i]);
	for (size_t i=0; i < free_me->rowct; i++)  free(free_me->row[i]);
    if (free_me->vector) free(free_me->vector);
	free(free_me->col);
	free(free_me->text);
	free(free_me->row);
	free(free_me);
}

/** Append one list of names to another.

If the first list is empty, then this is a copy function.

\param  n1      The first set of names (no default, must not be \c NULL)
\param  nadd      The second set of names, which will be appended after the first. (no default. If \c NULL, a no-op.)
\param type1     Either 'c', 'r', 't', or 'v' stating whether you are merging the
columns, rows, text, or vector. If 'v', then ignore \c typeadd and just overwrite the
target vector name with the source name. (default: 'r')
\param typeadd     Either 'c', 'r', 't', or 'v' stating whether you are merging the columns, rows, or text. If 'v', then overwrite the target with the source vector name. (default: type1)
*/
APOP_VAR_HEAD void  apop_name_stack(apop_name * n1, apop_name *nadd, char type1, char typeadd){
    apop_name * apop_varad_var(nadd, NULL); 
    if (!nadd) return;
    apop_name * apop_varad_var(n1, NULL);
    Apop_stopif(!n1, return, 0, "Can't stack onto a NULL set of names (which n1 is).");
    char apop_varad_var(type1, 'r');
    char apop_varad_var(typeadd, type1);
APOP_VAR_ENDHEAD
    int i;
    apop_name counts = (apop_name){.rowct=nadd->rowct, .textct = nadd->textct, .colct = nadd->colct};//Necessary when stacking onto self.;
    if (typeadd == 'v')
        apop_name_add(n1, nadd->vector, 'v');
    else if (typeadd == 'r')
        for (i=0; i< counts.rowct; i++)
            apop_name_add(n1, nadd->row[i], type1);
    else if (typeadd == 't')
        for (i=0; i< counts.textct; i++)
            apop_name_add(n1, nadd->text[i], type1);
    else if (typeadd == 'c')
        for (i=0; i< counts.colct; i++)
            apop_name_add(n1, nadd->col[i], type1);
    else Apop_notify(1, "'%c' sent to apop_name_stack, but the only "
                        "valid options are r t c v. Doing nothing.", typeadd);
}

/** Copy one \ref apop_name structure to another. That is, all data is duplicated.

Used internally by \ref apop_data_copy, but sometimes useful by itself. For example,
say that we have an \ref apop_data struct named \c d and a \ref gsl_matrix of the same
dimensions named \c m; we could give \c m the labels from \c d for printing:
\code
apop_data *wrapped = &(apop_data){.matrix=m, .names=apop_name_copy(d)};
apop_data_print(wrapped);
apop_name_free(wrapped->names); //wrapped itself is auto-allocated; do not free.
\endcode
 
\param in The input names
\return   A \ref apop_name struct with copies of all input names.
*/
apop_name * apop_name_copy(apop_name *in){
    apop_name *out = apop_name_alloc();
    apop_name_stack(out, in, 'v');
    apop_name_stack(out, in, 'c');
    apop_name_stack(out, in, 'r');
    apop_name_stack(out, in, 't');
    Asprintf(&out->title, "%s", in->title);
    return out;
}

/** Finds the position of an element in a list of names.

The function uses POSIX's \c strcasecmp, and so does case-insensitive search the way that function does.

\param n        the \ref apop_name object to search.
\param name     the name you seek; see above.
\param type     \c 'c' (=column), \c 'r' (=row), or \c 't' (=text). Default is \c 'c'.
\return         The position of \c findme. If \c 'c', then this may be -1, meaning the vector name. If not found, returns -2.  On error, e.g. <tt>name==NULL</tt>, returns -2.
*/
int apop_name_find(const apop_name *n, const char *name, const char type){
    Apop_stopif(!name, return -2, 0, "You asked me to search for NULL.");
    char **list;
    int listct;
    if (type == 'r' || type == 'R'){
        list = n->row;
        listct = n->rowct;
    }
    else if (type == 't' || type == 'T'){
        list = n->text;
        listct = n->textct;
    }
    else { // default type == 'c'
        list = n->col;
        listct = n->colct;
    }
    for (int i = 0; i < listct; i++)
        if (!strcasecmp(name, list[i])) return i;

    if ((type=='c' || type == 'C') && n->vector && !strcasecmp(name, n->vector)) return -1;
    return -2;
}
