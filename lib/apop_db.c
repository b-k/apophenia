/** \file apop_db.c	An easy front end to SQLite. Includes a few nice
features like a variance, skew, and kurtosis aggregator for SQL.

 Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
*/

/** \defgroup db Database utilities 

These are convenience functions to handle interaction with SQLite. They
open one and only one database, and handle most of the interaction
therewith for you.

You will probably first use \ref apop_text_to_db to pull data into
the database, then \ref apop_query to clean the data in the database,
and finally \ref apop_query_to_data to pull some subset of the data
out for analysis.

\par Querying 
\li \ref apop_query_db: Manipulate the database, return nothing (e.g., input data).

\li \ref apop_query_to_data: Pull data into an apop_data set.

\li \ref apop_query_to_matrix: Pull data into a \c gsl_matrix.

\li \ref apop_query_to_float: Pull out a single number.

\li \ref apop_query_to_chars: Pull out columns of not-numbers.

\par Maintenance 
\li \ref apop_open_db: Optional, for when you want to use a database on disk.

\li \ref apop_close_db: If you used \ref apop_open_db, you will need to use this too.

\li \ref apop_table_exists: Check to make sure you aren't reinventing or destroying data. Also, the clean way to drop a table.

\li \ref apop_count_cols: Count the columns in a table.

\li \ref apop_db_merge: Import or merge the whole of another database into the currently open db.

\li \ref apop_db_merge_table: Import/merge just one table.

\par See also
 \li The \ref conversions, including \ref apop_convert_text_to_db and \ref apop_matrix_to_db.

 \li The \ref command_line "Command-line utilities".

\par P.S.
Apophenia reserves the right to insert temp tables into the opened database. They will all have names beginning with "apop_", so the reader is advised to not use tables with such names, and is free to ignore or delete any such tables that turn up.
 */
#include <math.h> 	                //sqrt
#include <string.h>
#include <stdarg.h>
#include <apophenia/types.h>
#include <gsl/gsl_math.h>           //GSL_NAN
#include <apophenia/db.h>
#include <apophenia/linear_algebra.h>
#include <apophenia/stats.h>	    //t_dist
#include <apophenia/regression.h>	//two_tailify

#include <apophenia/vasprintf.h>

sqlite3	*db=NULL;	                //There's only one database handle. Here it is.

apop_name *last_names = NULL;	    //The column names from the last query to matrix

int	total_rows, total_cols;		    //the counts from the last query.

/** This variable turns on some notifications. */
int apop_verbose	= 0;

                                                                                                                               
////////////////////////////////////////////////
// Part one: additional aggregate functions for calculating higher moments
////////////////////////////////////////////////

/** \page db_moments Database moments (plus pow()!)
\verbatim
select count(x), stddev(x), avg(x), var(x), variance(x), skew(x), kurt(x), kurtosis(x)
from table
group by whatever
\endverbatim

\verbatim
select pow(x,0.5), exp(x), log(x)
from table
\endverbatim

The SQL standard includes the <tt>count(x)</tt> and <tt>avg(x)</tt> aggregators,
but statisticians are usually interested in higher moments as well---at
least the variance. Therefore, SQL queries using the Apophenia library
may include any of the moments above.

<tt>var</tt> and <tt>variance</tt>; <tt>kurt</tt> and <tt>kurtosis</tt> do the same
thing. Choose the one that sounds better to you.

The  var/skew/kurtosis functions calculate ''sample'' moments, so if you want the population moment, multiply the result by (n-1)/n .

For bonus points, there are the <tt>pow(x,y)</tt>, <tt>exp(x)</tt>, and <tt>log(x)</tt> functions. They call the standard math library function of the same name to calculate \f$x^y\f$, \f$e^x\f$, and \f$\ln(x)\f$.
*/


typedef struct StdDevCtx StdDevCtx;
struct StdDevCtx {
  double avg;     /* avg of terms */
  double avg2;    /* avg of the squares of terms */
  double avg3;    /* avg of the cube of terms */
  double avg4;    /* avg of the fourth-power of terms */
  int cnt;        /* Number of terms counted */
};

static void twoStep(sqlite3_context *context, int argc, sqlite3_value **argv){
  StdDevCtx *p;
double 		x, ratio;
  if( argc<1 ) return;
  p = sqlite3_aggregate_context(context, sizeof(*p));
  if( p && argv[0] ){
    x = sqlite3_value_double(argv[0]);
    ratio	=  p->cnt/(p->cnt+1.0);
    p->cnt++;
    p->avg	*= ratio;
    p->avg2	*= ratio;
    p->avg += x/(p->cnt +0.0);
    p->avg2 += gsl_pow_2(x)/(p->cnt +0.0);
  }
}

static void threeStep(sqlite3_context *context, int argc, sqlite3_value **argv){
StdDevCtx 	*p;
double 		x, ratio;
  if( argc<1 ) return;
  p = sqlite3_aggregate_context(context, sizeof(*p));
  if( p && argv[0] ){
    x = sqlite3_value_double(argv[0]);
    ratio	=  p->cnt/(p->cnt+1.0);
    p->cnt++;
    p->avg	*= ratio;
    p->avg2	*= ratio;
    p->avg3	*= ratio;
    p->avg += x/p->cnt;
    p->avg2 += gsl_pow_2(x)/p->cnt;
    p->avg3 += gsl_pow_3(x)/p->cnt;
  }
}

static void fourStep(sqlite3_context *context, int argc, sqlite3_value **argv){
StdDevCtx 	*p;
double 		x,ratio;
  if( argc<1 ) return;
  p = sqlite3_aggregate_context(context, sizeof(*p));
  if( p && argv[0] ){
    x = sqlite3_value_double(argv[0]);
    ratio	=  p->cnt/(p->cnt+1.0);
    p->cnt++;
    p->avg	*= ratio;
    p->avg2	*= ratio;
    p->avg3	*= ratio;
    p->avg4	*= ratio;
    p->avg += x/p->cnt;
    p->avg2 += gsl_pow_2(x)/p->cnt;
    p->avg3 += gsl_pow_3(x)/p->cnt;
    p->avg4 += gsl_pow_4(x)/p->cnt;
  }
}

static void stdDevFinalize(sqlite3_context *context){
  StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
  if( p && p->cnt>1 ){
    double rCnt = p->cnt;
    sqlite3_result_double(context,
       sqrt((p->avg2*rCnt - p->avg*p->avg)/(rCnt-1.0)));
  } else if (p->cnt == 1)
    	sqlite3_result_double(context, 0);
}

static void varFinalize(sqlite3_context *context){
  StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
  if( p && p->cnt>1 ){
    double rCnt = p->cnt;
    sqlite3_result_double(context,
       (p->avg2 - gsl_pow_2(p->avg))*rCnt/(rCnt-1.0));
  } else if (p->cnt == 1)
    	sqlite3_result_double(context, 0);
}

static void skewFinalize(sqlite3_context *context){
  StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
  if( p && p->cnt>1 ){
    double rCnt = p->cnt;
    sqlite3_result_double(context,
       (p->avg3*rCnt - 3*p->avg2*p->avg*rCnt + (3*rCnt-1) * gsl_pow_3(p->avg)) / (rCnt-1.0));
  } else if (p->cnt == 1)
    	sqlite3_result_double(context, 0);
}

static void kurtFinalize(sqlite3_context *context){
  StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
  if( p && p->cnt>1 ){
    double rCnt = p->cnt;
    sqlite3_result_double(context,
       (p->avg4*rCnt - 4*p->avg3*p->avg*rCnt + 6 * gsl_pow_2(p->avg2)*gsl_pow_2(p->avg)*rCnt
						- (4*rCnt+1)* gsl_pow_4(p->avg))/(rCnt-1.0));
  } else if (p->cnt == 1)
    	sqlite3_result_double(context, 0);
}

static void powFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, pow(sqlite3_value_double(argv[0]), sqlite3_value_double(argv[1])));
}

static void expFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, exp(sqlite3_value_double(argv[0])));
}

static void logFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, log(sqlite3_value_double(argv[0])));
}


////////////////////////////////////////////////
// Part two: database querying functions, so the user doesn't have to
// touch sqlite3.
////////////////////////////////////////////////



/**
If you want to use a database on the hard drive instead of memory,
then call this once and only once before using any other database
utilities. 

If you want a disposable database which you won't use after the program
ends, don't bother with this function.

The trade-offs between an on-disk database and an in-memory db are as
one would expect: memory is faster, but is destroyed when the program
exits. SQLite includes a command line utility (<tt>sqlite3</tt>) which
let you ask queries of a database on disk, which may
be useful for debugging. There are also some graphical
front-ends; just ask your favorite search engine for <a
href="http://www.google.com/search?&q=sqlite+gui">SQLite GUI</a>.

The Apophenia package assumes you are only using a single
SQLite database at a time; if not, the \ref apop_db_merge and \ref
apop_db_merge_table functions may help.

When you are done doing your database manipulations, be sure to call \ref apop_close_db .

\param filename
The name of a file on the hard drive on which to store the database. If
<tt>NULL</tt>, then the database will be kept in memory (in which case,
the other database functions will call this function for you and you
don't need to bother).

\ingroup db
*/
int apop_db_open(char *filename){
//char	*err;
	//if (filename==NULL) 	db	=sqlite_open(":memory:",0,&err);
	//else			db	=sqlite_open(filename,0,&err);
	if (filename==NULL) 	sqlite3_open(":memory:",&db);
	else			sqlite3_open(filename,&db);
	if (db == NULL)	
		{printf("Not sure why, but the database didn't open.\n");
		return 1; }
	sqlite3_create_function(db, "stddev", 1, SQLITE_ANY, NULL, NULL, &twoStep, &stdDevFinalize);
	sqlite3_create_function(db, "var", 1, SQLITE_ANY, NULL, NULL, &twoStep, &varFinalize);
	sqlite3_create_function(db, "variance", 1, SQLITE_ANY, NULL, NULL, &twoStep, &varFinalize);
	sqlite3_create_function(db, "skew", 1, SQLITE_ANY, NULL, NULL, &threeStep, &skewFinalize);
	sqlite3_create_function(db, "kurt", 1, SQLITE_ANY, NULL, NULL, &fourStep, &kurtFinalize);
	sqlite3_create_function(db, "kurtosis", 1, SQLITE_ANY, NULL, NULL, &fourStep, &kurtFinalize);
	sqlite3_create_function(db, "pow", 2, SQLITE_ANY, NULL, &powFn, NULL, NULL);
	sqlite3_create_function(db, "exp", 2, SQLITE_ANY, NULL, &expFn, NULL, NULL);
	sqlite3_create_function(db, "log", 2, SQLITE_ANY, NULL, &logFn, NULL, NULL);
	apop_query("pragma short_column_names");
	return 0;
}

/** An alias for \ref apop_db_open. Use that one.*/ 
int apop_open_db(char *filename){ return apop_db_open(filename); }

/** Send a query to the database, return nothing 
\param fmt
An \ref sql "SQL" query.

\param ...
Your query may be in <tt>printf</tt> form. For example:
\code
char tabname[] = "demographics";
char colname[] = "heights";
int min_height = 175;
apop_query("select %s from %s where %s > %i", colname, tabname, colname, min_height);
\endcode
\ingroup db
*/
int apop_query(const char *fmt, ...){
char 		*err, *q;
va_list		argp;
	if (db==NULL) apop_db_open(NULL);
	va_start(argp, fmt);
	vasprintf(&q, fmt, argp);
	va_end(argp);
	if (apop_verbose) {printf("\n%s\n",q);}
	sqlite3_exec(db, q, NULL,NULL, &err);
	free(q);
	ERRCHECK
	return 1;
}

/** Just like \ref apop_query. Use that one. 
\ingroup db
*/
int apop_query_db(const char *fmt, ...){
char 		*err, *q;
va_list		argp;
	if (db==NULL) apop_db_open(NULL);
	va_start(argp, fmt);
	vasprintf(&q, fmt, argp);
	va_end(argp);
	if (apop_verbose) {printf("\n%s\n",q);}
	sqlite3_exec(db, q, NULL,NULL, &err);
	free(q);
	ERRCHECK
	return 1;
}

int 		isthere;//used for the next two fns only.

int tab_exists_callback(void *in, int argc, char **argv, char **whatever){
char *q	= in;
	if (!strcmp(argv[argc-1],q))
		isthere=1;
	return isthere;
}

/** Check for the existence of a table, and maybe delete it.

Recreating a table which already exists can cause errors, so it is good practice to check for existence first.
Also, this is the stylish way to delete a table, since just calling <tt>"drop table"</tt> will give you an error if the table doesn't exist.

\param q 	the table name
\param whattodo 1	==>kill table so it can be recreated in main.<br>
		0	==>return error so program can continue.
\return
0 = table does not exist<br>
1 = table was found, and if whattodo==1, has been deleted
\ingroup db
*/
int apop_table_exists(char *q, int whattodo){
char 		*err, q2[10000];
	isthere=0;
	if (db==NULL) {apop_db_open(NULL); return 0;}
	sqlite3_exec(db, "select name from sqlite_master where type='table'",tab_exists_callback,q, &err); 
	ERRCHECK
	if (whattodo==1 && isthere)
		sqlite3_exec(db,strcat(strcpy(q2, "DROP TABLE "),q),NULL,NULL, &err); ERRCHECK
	return isthere;
}


/** Counts the number of columns in a table. Occasionally useful, e.g., when data is read from a file.

\param name 	The name of the table you are inquiring about.
\return 	The number of columns the table has.
*/
int apop_count_cols(const char *name){
char 		*err, q2[5000];
int		colct	= 1;

	int count_cols_callback(void *whatever, int argc, char **argv, char **andever){
	int 		i=0;
		while(argv[0][i]!='\0')
			if (argv[0][i++]==',') 
				colct	++;
		return 0;
	}

	if (db==NULL) {printf("No database open yet."); return 0;}
	sprintf(q2, "select sql from sqlite_master where type='table' and name=\"%s\"",name);
	sqlite3_exec(db, q2, count_cols_callback,NULL, &err); 
	ERRCHECK
	return colct;
}

/**
Closes the database on disk. If you opened the database with
<tt>apop_open_db(NULL)</tt>, then this is basically optional.

\param vacuum 
1: Do clean-up to minimize the size of the database on disk.<br>
0: Don't bother; just close the database.
*/
int apop_db_close(int vacuum){
char		*err;
	if (vacuum) sqlite3_exec(db, "VACUUM", NULL, NULL, &err);
//	ERRCHECK
	sqlite3_close(db);
	return 0;
	}

/** An alias for \ref apop_db_close . */
int apop_close_db(int vacuum){ return apop_db_close(vacuum); }

int names_callback(void *o,int argc, char **argv, char **whatever){
	apop_name_add(last_names, argv[1], 'c'); 
	return 0;
}

int length_callback(void *o,int argc, char **argv, char **whatever){
	total_rows=atoi(argv[0]); 
	return 0;
}

/** Dump the results of a query into an array of strings.

You will probably be curious as to the dimensions of the
array just returned. To get this information, use \ref apop_db_get_rows
and \ref apop_db_get_cols .

\param fmt 	As with \ref apop_query , a string containing a query,
which may include <tt>printf</tt>-style tags (<tt>\%i, \%s</tt>, et cetera).

\return		An array of strings. Notice that this is always a 2-D
array, even if the query returns a single column. In that case, use
<tt>returned_tab[i][0]</tt> to refer to row <tt>i</tt>.

example
The following function will list the tables in a database (much like you could do from the command line using <tt>sqlite3 dbname.db ".table"</tt>).

\verbatim
void print_table_list(char *db_file){
char            ***tab_list;
int             row_ct, i;
        apop_db_open(db_file);
        tab_list= apop_query_to_chars("select name from sqlite_master where type==\"table\";");
        row_ct  =  apop_db_get_rows();
        for(i=0; i< row_ct; i++)
                printf("%s\n", tab_list[i][0]);
}
\endverbatim
*/
char *** apop_query_to_chars(const char * fmt, ...){
char		***output;
int		currentrow=0;
char		*q2, *err=NULL, *query;
va_list		argp;

	int db_to_chars(void *o,int argc, char **argv, char **whatever){
	int		jj;
	char ****	output = (char ****) o;
		if (*argv !=NULL){
			(*output)[currentrow]	= malloc(sizeof(char**) * argc);
			for (jj=0;jj<argc;jj++){
				if (argv[jj]==NULL){
					(*output)[currentrow][jj]	= malloc(sizeof(char*));
					strcpy((*output)[currentrow][jj], "");
				}
				else{
					(*output)[currentrow][jj]	= malloc(sizeof(char*) * strlen(argv[jj]));
					strcpy((*output)[currentrow][jj], argv[jj]);
				}
			}
			currentrow++;
		}
		return 0;
	}

	if (db==NULL) apop_db_open(NULL);
	va_start(argp, fmt);
	vasprintf(&query, fmt, argp);
	if (apop_verbose) {printf("\n%s\n",query);}
	va_end(argp);

	total_rows	= 0;
	q2		= malloc(sizeof(char)*(strlen(query)+300));
	apop_table_exists("apop_temp_table",1);
	sqlite3_exec(db,strcat(strcpy(q2,
		"CREATE TABLE apop_temp_table AS "),query),NULL,NULL, &err); ERRCHECK
	sqlite3_exec(db,"SELECT count(*) FROM apop_temp_table",length_callback,NULL, &err);
	free(query);
	ERRCHECK
	if (total_rows==0){
		output	= NULL;
	} else {
		total_cols	= apop_count_cols("apop_temp_table");
		output		= malloc(sizeof(char***) * total_rows);
		sqlite3_exec(db,"SELECT * FROM apop_temp_table",db_to_chars,&output, &err); ERRCHECK
		if (last_names !=NULL) 
			apop_name_free(last_names); 
		last_names = apop_name_alloc();
		sqlite3_exec(db,"pragma table_info(apop_temp_table)",names_callback, NULL, &err); ERRCHECK
	}
	sqlite3_exec(db,"DROP TABLE apop_temp_table",NULL,NULL, &err);  ERRCHECK
	free(q2);
	return output;
}

/** Queries the database, and dumps the result into a matrix.


\param fmt 	A string holding an \ref sql "SQL" query.
Your query may be in <tt>printf</tt> form. See \ref apop_query for an example.

\return
A <tt>gsl_matrix</tt>, which you passed in declared but not allocated.

Blanks in the database are filled with <tt>GSL_NAN</tt>s in the matrix.
*/
gsl_matrix * apop_query_to_matrix(const char * fmt, ...){
gsl_matrix	*output;
int		currentrow=0;
char		*q2, *err=NULL, *query;
va_list		argp;

	int db_to_table(void *o,int argc, char **argv, char **whatever){
	int		jj;
	gsl_matrix * 	output = (gsl_matrix *) o;
		if (*argv !=NULL){
			for (jj=0;jj<argc;jj++){
				if (argv[jj]==NULL)
					gsl_matrix_set(output,currentrow,jj, GSL_NAN);
				else
					gsl_matrix_set(output,currentrow,jj, atof(argv[jj]));
			}
			currentrow++;
		}
		return 0;
	}

	if (db==NULL) apop_db_open(NULL);
	va_start(argp, fmt);
	vasprintf(&query, fmt, argp);
	va_end(argp);
	if (apop_verbose)	printf("%s\n", query);
	total_rows= 0;
	q2	 = malloc(sizeof(char)*(strlen(query)+300));
	apop_table_exists("apop_temp_table",1);
	sqlite3_exec(db,strcat(strcpy(q2,
		"CREATE TABLE apop_temp_table AS "),query),NULL,NULL, &err); ERRCHECK
	sqlite3_exec(db,"SELECT count(*) FROM apop_temp_table",length_callback,NULL, &err);
	free(query);
	ERRCHECK
	if (total_rows==0){
		output	= NULL;
	} else {
		total_cols	= apop_count_cols("apop_temp_table");
		output		= gsl_matrix_alloc(total_rows, total_cols);
		sqlite3_exec(db,"SELECT * FROM apop_temp_table",db_to_table,output, &err); ERRCHECK
		if (last_names !=NULL) 
			apop_name_free(last_names); 
		last_names = apop_name_alloc();
		sqlite3_exec(db,"pragma table_info(apop_temp_table)",names_callback, NULL, &err); ERRCHECK
	}
	sqlite3_exec(db,"DROP TABLE apop_temp_table",NULL,NULL, &err);  ERRCHECK
	free(q2);
	return output;
}

/** Queries the database, and dumps the result into a single floating point number.
\param fmt	A string holding an \ref sql "SQL" query.
\param ...	Your query may be in <tt>printf</tt> form. See \ref apop_query for an example.
\return		A float. This calls \ref apop_query_to_matrix and returns
the (0,0)th element of the returned matrix. Thus, if your query returns
multiple lines, you will get no warning, and the function will return
the first in the list (which is not always well-defined).

If the query returns no rows at all, the function returns <tt>GSL_NAN</tt>.
*/
float apop_query_to_float(const char * fmt, ...){
gsl_matrix	*m=NULL;
va_list		argp;
char		*query;
float		out;
	va_start(argp, fmt);
	vasprintf(&query, fmt, argp);
	va_end(argp);
	m	= apop_query_to_matrix(query);
	if (m==NULL){
        if (apop_verbose)
		    printf("apop, %s, %i: Query turned up a blank table. Returning GSL_NAN.\n", __FILE__, __LINE__);
		return GSL_NAN;
	} //else
	out	= gsl_matrix_get(m, 0, 0);
	gsl_matrix_free(m);
	return out;

}

/** Queries the database, and dumps the result into an \ref apop_data set.

\param fmt 	A string holding an \ref sql "SQL" query.
Your query may be in <tt>printf</tt> form. See \ref apop_query for an example.

\return
An \ref apop_data set, which you passed in declared but not allocated.

Blanks in the database are filled with <tt>GSL_NAN</tt>s in the matrix.
\bug Currently, this is but a wrapper for \ref apop_query_to_matrix,
meaning that only numerical results are returned. If you want
non-numeric data, try \code mydata->categories  = apop_query_to_chars("select ...");\endcode. 
*/ 
apop_data * apop_query_to_data(const char * fmt, ...){
gsl_matrix	*m=NULL;
va_list		argp;
char		*query;
apop_data	*out;
	va_start(argp, fmt);
	vasprintf(&query, fmt, argp);
	va_end(argp);
	m	        = apop_query_to_matrix(query);
    out         = apop_matrix_to_data(m);
    //replace name struct allocated in apop_matrix_to_data with the
    //actual names.
    apop_name_free(out->names);
    out->names  = apop_db_get_names();
	return out;
}

/** This function returns an \ref apop_name structure with the column
names from the last <tt>apop_query_...</tt> . Since only the names from
the last query are saved, you will want to use this immediately
after your query.  */
apop_name * apop_db_get_names(void){ return last_names; }

/** This function returns the column count from the last query run. */
int apop_db_get_cols(void){ return total_cols; }

/** This function returns the row count from the last query run. */
int apop_db_get_rows(void){ return total_rows; }

/** Dump a <tt>gsl_matrix</tt> into the database.

\param data 	The name of the matrix
\param tabname	The name of the db table to be created
\param headers	A list of column names. If <tt>NULL</tt>, then the columns will be named <tt>c1</tt>, <tt>c2</tt>, <tt>c3</tt>, &c.
*/
int apop_matrix_to_db(gsl_matrix *data, char *tabname, char **headers){
int		i,j; 
int		ctr		= 0;
int		batch_size	= 100;
char		*q 		= malloc(sizeof(char)*1000);
	if (db==NULL) apop_db_open(NULL);
	sprintf(q, "create table %s (", tabname);
	for(i=0;i< data->size2; i++){
		q	=realloc(q,sizeof(char)*(strlen(q)+1000));
		if(headers == NULL) 	sprintf(q, "%s\n c%i", q,i);
		else			sprintf(q, "%s\n %s", q,headers[i]);
		if (i< data->size2-1) 	sprintf(q, "%s,",q);
		else			sprintf(q,"%s);  begin;",q);
	}
	for(i=0;i< data->size1; i++){
		q	=realloc(q,sizeof(char)*(strlen(q)+(1+data->size2)*1000));
		sprintf(q,"%s \n insert into %s values(",q,tabname);
		for(j=0;j< data->size2; j++)
			if(j< data->size2 -1 && ctr<batch_size) {
				sprintf(q,"%s %g, ",q,gsl_matrix_get(data,i,j));
				ctr++;
			} else	{
				sprintf(q,"%s %g);",q,gsl_matrix_get(data,i,j));
				sprintf(q, "%s; end;", q);
				apop_query(q);
				ctr = 0;
				sprintf(q,"begin; \n insert into %s values(",tabname);
			}
		sprintf(q,"%s )",q);
	}
	free(q);
	return 0;
}

void free_tab_list(char ****tab, int row_ct, int col_ct){
int		i,j;
	for(i=0; i< row_ct; i++){
		for(j=0; j< col_ct; j++)
			free((*tab)[i][j]); 
		free((*tab)[i]);
	}
	free(*tab);
}

/** Merge a single table from a database on the hard drive with the database currently open.

\param db_file	The name of a file on disk.
\param tabname	The name of the table in that database to be merged in.

If the table exists in the new database but not in the currently open one,
then it is simply copied over. If there is a table with the same name in
the currently open database, then the data from the new table is inserted
into the main database's table with the same name. [The function just
calls <tt>insert into main.tab select * from merge_me.tab</tt>.]

\ingroup db
\todo fix the tab_list bug.
*/
void apop_db_merge_table(char *db_file, char *tabname){
//char		***tab_list;
int		row_ct;
	if (db_file !=NULL)
		apop_query("attach database \"%s\" as merge_me;", db_file);
	apop_query_to_chars("select name from sqlite_master where name == \"%s\";", tabname);
	row_ct	= apop_db_get_rows();
	if (row_ct==0){	//just import table
		if (apop_verbose)	printf("adding in %s\n", tabname);
		apop_query("create table main.%s as select * from merge_me.%s;", tabname, tabname);
	}
	else	{			//merge tables.
		if (apop_verbose)	printf("merging in %s\n", tabname);
		apop_query("insert into main.%s select * from merge_me.%s;", tabname, tabname);
	}
	if (db_file !=NULL)
		apop_query("detach database merge_me;");
	/*if (*tab_list !=NULL)
		free_tab_list(&tab_list, row_ct, 1);*/
}

/** Merge a database on the hard drive with the database currently open.

\param db_file	The name of a file on disk.

If a table exists in the new database but not in the currently open one,
then it is simply copied over. If there are  tables with the same name
in both databases, then the data from the new table is inserted into
the main database's table with the same name. [The function just calls
<tt>insert into main.tab select * from merge_me.tab</tt>.]

\ingroup db
*/
void apop_db_merge(char *db_file){
char		***tab_list;
int		row_ct, i;
	apop_query("attach database \"%s\" as merge_me;", db_file);
	tab_list= apop_query_to_chars("select name from merge_me.sqlite_master where type==\"table\";");
	row_ct	=  apop_db_get_rows();
	for(i=0; i< row_ct; i++)
		apop_db_merge_table(NULL, tab_list[i][0]);
	apop_query("detach database merge_me;");
	free_tab_list(&tab_list, row_ct, 1);
}
                                                                                                                               
////////////////////////////////////////////////
// Part three: some stats wrappers
////////////////////////////////////////////////

/** Do a t-test entirely inside the database.
\ingroup ttest
*/
double apop_db_t_test(char * tab1, char *col1, char *tab2, char *col2){
gsl_matrix	*result1, *result2;
	result1	= apop_query_to_matrix("select avg(%s), var(%s), count(*) from %s", col1, col1, tab1);
	result2	= apop_query_to_matrix("select avg(%s), var(%s), count(*) from %s", col2, col2, tab2);
double		a_avg	= gsl_matrix_get(result1, 0, 0),
		a_var	= gsl_matrix_get(result1, 0, 1),
		a_count	= gsl_matrix_get(result1, 0, 2),
		b_avg	= gsl_matrix_get(result2, 0, 0),
		b_var	= gsl_matrix_get(result2, 0, 1),
		b_count	= gsl_matrix_get(result2, 0, 2),
		stat	= (a_avg - b_avg)/ sqrt(b_var/(b_count-1) + a_var/(a_count-1));
	return two_tailify(gsl_cdf_tdist_P(stat, a_count+b_count-2));
}

/** Do a paired t-test entirely inside the database.
\ingroup ttest
*/
double	apop_db_paired_t_test(char * tab1, char *col1, char *col2){
gsl_matrix	*result;
	result	= apop_query_to_matrix("select avg(%s - %s), var(%s - %s), count(*) from %s tab1", 
						   col1,col2,   col1, col2,          tab1);
double		avg	= gsl_matrix_get(result, 0, 0),
		var	= gsl_matrix_get(result, 0, 1),
		count	= gsl_matrix_get(result, 0, 2),
		stat	= avg/ sqrt(var/(count-1));
	return two_tailify(gsl_cdf_tdist_P(stat, count-1));
}
