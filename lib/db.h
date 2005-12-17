//db.h			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
#ifndef apop_db_included
#define apop_db_included
#include <sqlite3.h>
#include <apophenia/types.h>
#include <gsl/gsl_matrix.h>
#define ERRCHECK {if (err!=NULL) {printf("%s\n",err);  return 0;}}

/** Set this to zero for silent mode, one for errors and warnings.

  \ingroup global_vars */
extern int apop_verbose;

int apop_table_exists(char *q, int whattodo);
	//whattodo==1	==>kill table so it can be recreated in the main.
	//whattodo==0	==>just return the status: 1=exists 0=doesn't. 

int apop_count_cols(const char *name);
	//give me the name of a table, I'll check sqlite_master for the 
	//create statement that made it, and will use that to tell you
	//how many columns are in the table.

int apop_db_open(char *filename);
int apop_open_db(char *filename);
	//If filename==NULL, it'll open a database in memory

int apop_db_close(int vacuum);
int apop_close_db(int vacuum);
	//vacuum==1: do some cleanup to minimize hard disk space
	//vacuum==0: just close the thing.

int apop_query(const char *q, ...);
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

char *** apop_query_to_chars(const char * fmt, ...);

apop_data * apop_query_to_data(const char * fmt, ...);

float apop_query_to_float(const char * fmt, ...);
	//like query_to_matrix, but returns a single number.

int apop_matrix_to_db(gsl_matrix *data,char *tabname, char **headers);
	//dump a matrix to a database table named tabname.
	//At the moment, the headers are ignored. 
	//With no headers specified, you get columns C0, C1, C2...

apop_name * apop_db_get_names(void);
int apop_db_get_cols(void);
int apop_db_get_rows(void);
	//give the column names and counts from the last query.

void apop_db_merge(char *infile);
	//copy all of the tables in the database at the given file into
	//the database apophenia has already opened. If there are
	//duplicate names, append them.

void apop_db_merge_table(char *infile, char *tabname);
	//This is just like apop_db_merge, but will
	//only pull the single table specified.

double apop_db_t_test(char * tab1, char *col1, char *tab2, char *col2);
double apop_db_paired_t_test(char * tab1, char *col1, char *col2);
	//Runs a t-test entirely inside the database. Uses the nifty
	//var() aggregator function defined by Apophenia.
#endif
