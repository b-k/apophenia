#ifndef __apop_name__
#define __apop_name__

#include <string.h>

/** \defgroup names apop_names: a structure and functions to handle column names for matrices. 
\ingroup types */


/** A data set is assumed to be a matrix where each row is a single
observation and each column is a variable. Usually there is only one
dependent variable (the value to be predicted), which is the first column;
the independent variables (the predictors) follow thereafter.

This structure holds the names of these variables. You can fill it quickly
with \ref apop_db_get_names after running a query, or add names manually
with \ref apop_name_add .

Typically, the row names are not used, but they are there for your convenience.  
\ingroup names
*/
typedef struct apop_name{
	char ** colnames;
	char ** rownames;
	char ** depnames;
	int colnamect, depnamect, rownamect;
} apop_name;

apop_name * apop_name_alloc(void);
int apop_name_add(apop_name * n, char *add_me, char type);
void  apop_name_free(apop_name * free_me);
void  apop_name_print(apop_name * n);

#endif
