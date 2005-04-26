//db.c			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
#include <sqlite3.h>
#include <gsl/gsl_matrix.h>
#define ERRCHECK {if (err!=NULL) {printf("%s\n",err);  return 0;}}


int apop_table_exists(const char *q, int whattodo);
	//whattodo==1	==>kill table so it can be recreated in the main.
	//whattodo==0	==>just return the status: 1=exists 0=doesn't. 

int apop_count_cols(const char *name);
	//give me the name of a table, I'll check sqlite_master for the 
	//create statement that made it, and will use that to tell you
	//how many columns are in the table.

int apop_open_db(char *filename);
	//If filename==NULL, it'll open a database in memory

int apop_close_db(int vacuum);
	//vacuum==1: do some cleanup to minimize hard disk space
	//vacuum==0: just close the thing.

int apop_query_db(const char *q, ...);
	//Run a query but output nothing outside the DB.
	//It's fastest to compound as many queries as possible here;
	//q can be a several-line string of the form:
	//	Begin;
	//	create table blah ...;
	//	insert into blah ...;
	//	commit;
	//Also, you may use printf-type queries:
	//apop_query_db("select * from %s where date > %i", tabname, earliest);

gsl_matrix * apop_query_to_matrix(const char * fmt, ...);
	//dump a query to a matrix. 
	//do not preallocate *output.
	//	gsl_matrix * outmatrix;
	//	outmatrix = apop_query_to_matrix("select a, b, c from some_table");

int apop_matrix_to_db(gsl_matrix *data,char *tabname, char **headers);
	//dump a matrix to a database table named tabname.
	//At the moment, the headers are ignored. 
	//With no headers specified, you get columns C0, C1, C2...
