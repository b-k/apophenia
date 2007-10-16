/** \file apop_db.c	An easy front end to SQLite. Includes a few nice
features like a variance, skew, and kurtosis aggregator for SQL.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

/** \defgroup db Database utilities 

These are convenience functions to handle interaction with SQLite. They
open one and only one database, and handle most of the interaction
therewith for you.

You will probably first use \ref apop_text_to_db to pull data into
the database, then \ref apop_query to clean the data in the database,
and finally \ref apop_query_to_data to pull some subset of the data
out for analysis.

\par Querying 
\li \ref apop_query: Manipulate the database, return nothing (e.g., input data).

\li \ref apop_query_to_data: Pull data into an apop_data set.

\li \ref apop_query_to_matrix: Pull data into a \c gsl_matrix.

\li \ref apop_query_to_float: Pull out a single number.

\li \ref apop_query_to_text: Pull out columns of not-numbers.

\li \ref apop_query_to_multi: Pull data into an apop_data set, but with mixed numeric/text types.

\par Maintenance 
\li \ref apop_db_open: Optional, for when you want to use a database on disk.

\li \ref apop_db_close: If you used \ref apop_db_open, you will need to use this too.

\li \ref apop_table_exists: Check to make sure you aren't reinventing or destroying data. Also, the clean way to drop a table.

\li \ref apop_count_cols: Count the columns in a table.

\li \ref apop_db_merge: Import or merge the whole of another database into the currently open db.

\li \ref apop_db_merge_table: Import/merge just one table.

\par See also
 \li The \ref conversions, including \ref apop_text_to_db and \ref apop_matrix_to_db.

 \li The \ref command_line "Command-line utilities".

\par P.S.
Apophenia reserves the right to insert temp tables into the opened database. They will all have names beginning with "apop_", so the reader is advised to not use tables with such names, and is free to ignore or delete any such tables that turn up.
 */
#include <string.h>
#include <stdarg.h>
#include <apophenia/types.h>
#include <gsl/gsl_math.h>           //GSL_NAN
#include <apophenia/db.h>
#include <apophenia/linear_algebra.h>
#include <apophenia/stats.h>	    //t_dist
#include <apophenia/conversions.h>	//apop_strip_dots
#include <apophenia/regression.h>	//two_tailify
#include <apophenia/bootstrap.h>	//apop_rng_alloc

#include <apophenia/vasprintf.h>
#include "config.h"
#include <math.h> 	                //sqrt

/** Here are where the options are initially set. */
apop_opts_type apop_opts	= { 0,              //verbose
                                'f',            //output type
                                NULL,            //output pipe
                                "\t",           //output delimiter
                                1,              //output append
                                "| ,\t",        //input delimiters
                                "row_names",    //db_name_column
                                "NaN", //db_nan
                                '\0',            //db_engine
                                1               //threadct
};


#ifdef HAVE_LIBMYSQLCLIENT
#include "apop_db_mysql.c"
#endif

apop_name *last_names = NULL;	    //The column names from the last query to matrix
static gsl_rng* db_rng  = NULL;     //the RNG for the RNG function.

#ifdef HAVE_LIBSQLITE3 
#include "apop_db_sqlite.c"
#endif
                                                                                                                               

/** Random numbers are generated inside the database using a separate
 RNG. This will initialize it for you, just like \ref apop_rng_alloc,
 except the RNG it produces is kept for internal use. If you don't call
 it, then it will be called at first use, with seed zero.

\param  seed    The seed. No need to get funny with it: 0, 1, and 2 will produce wholly different streams.
\return The RNG ready for your use.
\ingroup db
*/
void apop_db_rng_init(int seed){
    db_rng  = apop_rng_alloc(seed);
}



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

MySQL users: either set the environment variable APOP_DB_ENGINE=mysql or set \c apop_opts.db_engine = 'm'.

The Apophenia package assumes you are only using a single
SQLite database at a time; if not, the \ref apop_db_merge and \ref
apop_db_merge_table functions may help.

When you are done doing your database manipulations, be sure to call \ref apop_db_close .

\param filename
The name of a file on the hard drive on which to store the database. If
<tt>NULL</tt>, then the database will be kept in memory (in which case,
the other database functions will call this function for you and you
don't need to bother).

\return 0: everything OK<br>
        1: database did not open.

\ingroup db
*/
int apop_db_open(char *filename){
#ifdef HAVE_LIBSQLITE3
    if (!db) //check the environment.
#endif
#ifdef HAVE_LIBMYSQLCLIENT
       if(!mysql_db)  
#endif
        if (getenv("APOP_DB_ENGINE") && !strcasecmp(getenv("APOP_DB_ENGINE"), "mysql"))
            apop_opts.db_engine = 'm';

    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_db_open(filename);
#else
        {apop_error(0, 'c', "apop_db_open: Apophenia was compiled without mysql support.\n");
        return 0;
        }
#endif
#ifdef HAVE_LIBSQLITE3
        return apop_sqlite_db_open(filename);
#else
        {apop_error(0, 'c', "apop_db_open: Apophenia was compiled without sqlite support.\n");
        return 0;
        }
#endif
}


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
  va_list   argp;
	va_start(argp, fmt);
	vasprintf(&q, fmt, argp);
	va_end(argp);
	if (apop_opts.verbose) {printf("\n%s\n",q);}
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        {if (!mysql_db) {apop_error(0, 'c', "No database is open.");
            return 0;}
        apop_mysql_query(q);}
#else
        {apop_error(0, 'c', "apop_query: Apophenia was compiled without mysql support.\n");
        return 0;
        }
#endif
    else 
#ifdef HAVE_LIBSQLITE3
        {if (!db) apop_db_open(NULL);
        sqlite3_exec(db, q, NULL,NULL, &err);
	    ERRCHECK
        }
#else
        {apop_error(0, 'c', "apop_query: Apophenia was compiled without SQLite support.\n");
        return 0;
        }
#endif
	free(q);
	return 1;
}

int 		isthere;//used for the next two fns only.

static int tab_exists_callback(void *in, int argc, char **argv, char **whatever){
  char *q	= in;
	if (!strcmp(argv[argc-1],q))
		isthere=1;
	return isthere;
}

/** Check for the existence of a table, and maybe delete it.

Recreating a table which already exists can cause errors, so it is good practice to check for existence first.
Also, this is the stylish way to delete a table, since just calling <tt>"drop table"</tt> will give you an error if the table doesn't exist.

\param q 	the table name
\param whattodo 'd'	==>delete table so it can be recreated in main.<br>
		'n'	==>no action. return error so program can continue.
\return
0 = table does not exist<br>
1 = table was found, and if whattodo==1, has been deleted
\ingroup db
*/
int apop_table_exists(char *q, char whattodo){
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_table_exists(q, whattodo);
#else
        {apop_error(0, 'c', "apop_table_exists: Apophenia was compiled without mysql support.\n");
        return 0; }
#endif
#ifdef HAVE_LIBSQLITE3
  char 		*err, q2[10000];
	isthere=0;
	if (db==NULL) return 0;
	sqlite3_exec(db, "select name from sqlite_master where type='table'",tab_exists_callback,q, &err); 
	ERRCHECK
	if ((whattodo==1|| whattodo=='d') && isthere){
		sqlite3_exec(db,strcat(strcpy(q2, "DROP TABLE "),q),NULL,NULL, &err); 
        ERRCHECK
    }
	return isthere;
#endif
}

//colct is global for the count_cols callback.
int		colct;

static int count_cols_callback(void *whatever, int argc, char **argv, char **andever){
  int 		i=0;
	while(argv[0][i]!='\0')
		if (argv[0][i++]==',') 
			colct	++;
	return 0;
}

/** Counts the number of columns in a table. Occasionally useful, e.g., when data is read from a file.

\param name 	The name of the table you are inquiring about.
\return 	The number of columns the table has.
*/
int apop_count_cols(const char *name){
  char 		*err, q2[5000];
    colct   = 1;
	if (db==NULL) {printf("No database open yet."); return 0;}
	sprintf(q2, "select sql from sqlite_master where type='table' and name=\"%s\"",name);
	sqlite3_exec(db, q2, count_cols_callback,NULL, &err); 
	ERRCHECK
	return colct;
}

/**
Closes the database on disk. If you opened the database with
\c apop_db_open(NULL), then this is basically optional.

\param vacuum 
'v': vacuum---do clean-up to minimize the size of the database on disk.<br>
'q': Don't bother; just close the database.
*/
int apop_db_close(char vacuum){
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        {apop_mysql_db_close(0);
        return 0;}
#else
        {apop_error(0, 'c', "apop_db_close: Apophenia was compiled without mysql support.\n");
        return 0; }
#endif
    else{
#ifdef HAVE_LIBSQLITE3
  char		*err;
	if (vacuum==1 || vacuum=='v') 
        sqlite3_exec(db, "VACUUM", NULL, NULL, &err);
//	ERRCHECK
	sqlite3_close(db);
    db  = NULL;
	return 0;
#else
        {apop_error(0, 'c', "apop_db_close: Apophenia was compiled without SQLite support.\n");
        return 0; }
#endif
    }
}


/** Dump the results of a query into an array of strings.

\param fmt 	As with \ref apop_query , a string containing a query,
which may include <tt>printf</tt>-style tags (<tt>\%i, \%s</tt>, et cetera).

\return		An \ci apop_data structure with the <tt>text</tt>
element filled. Notice that this is always a 2-D array, even if the query
returns a single column. In that case, use <tt>returned_tab->text[i][0]</tt>
to refer to row <tt>i</tt>.


example
The following function will list the tables in a database (much like you could do from the command line using <tt>sqlite3 dbname.db ".table"</tt>).

\verbatim
#include <apophenia/headers.h>

void print_table_list(char *db_file){
apop_data   *tab_list;
int         i;
        apop_db_open(db_file);
        tab_list= apop_query_to_text("select name from sqlite_master where type==\"table\";");
        for(i=0; i< tab_list->textsize[0]; i++)
                printf("%s\n", tab_list->text[i][0]);
}

int main(int argc, char **argv){
    if (argc == 1){
        printf("Give me a database name, and I will print out the list of tables contained therein.\n");
        return 0; 
    }
    print_table_list(argv[1]);
}

\endverbatim
*/
apop_data * apop_query_to_text(const char * fmt, ...){
  va_list	argp;
  char		*query;
	va_start(argp, fmt);
	vasprintf(&query, fmt, argp);
	if (apop_opts.verbose) {printf("\n%s\n",query);}
	va_end(argp);
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_to_text(query);
#else
        {apop_error(0, 'c', "apop_query_to_text: Apophenia was compiled without mysql support.\n");
        return 0;}
#endif
#ifdef HAVE_LIBSQLITE3
        return apop_sqlite_query_to_text(query);
#else
        {apop_error(0, 'c', "apop_query_to_text: Apophenia was compiled without SQLite support.\n");
        return NULL; }
#endif
}

static int     firstcall;
static char        data_or_matrix  = 'm';
static char        full_divider[103];
static regmatch_t  result[3];
static regex_t     *regex;


size_t tr;
//apop_query_to_matrix callback.
static int db_to_table(void *o,int argc, char **argv, char **column){
  int		jj, i, ncfound = 0;
  gsl_matrix ** 	output = (gsl_matrix **) o;
    if (firstcall){
        firstcall   --;
        namecol     = -1;
        for(i=0; i<argc; i++)
            if (!strcmp(column[i], apop_opts.db_name_column)){
                namecol = i;
                ncfound = 1;
                break;
            }
	    *output		= gsl_matrix_alloc(tr, argc-ncfound);
        if (data_or_matrix == 'd')
            for(i=0; i<argc; i++)
                if (namecol != i)
                    apop_name_add(last_names, column[i], 'c');
    }
	if (argv !=NULL){
        ncfound =0;
		for (jj=0;jj<argc;jj++)
            if (jj != namecol){
                if (!argv[jj] || !regexec(regex, argv[jj], 1, result, 0))
				    gsl_matrix_set(*output,currentrow,jj-ncfound, GSL_NAN);
			    else
				    gsl_matrix_set(*output,currentrow,jj-ncfound, atof(argv[jj]));
            } else {
                apop_name_add(last_names, argv[jj], 'r');
                ncfound = 1;
            }
		currentrow++;
	}
	return 0;
}



/** Queries the database, and dumps the result into a matrix.


\param fmt 	A string holding an \ref sql "SQL" query.
Your query may be in <tt>printf</tt> form. See \ref apop_query for an example.

\return
A <tt>gsl_matrix</tt>, which you passed in declared but not allocated.

Blanks in the database (i.e., <tt> NULL</tt>s) and elements that match \ref apop_opts.db_nan
are filled with <tt>NAN</tt>s in the matrix.
*/
gsl_matrix * apop_query_to_matrix(const char * fmt, ...){
  char      *query;
  va_list	argp;
	va_start(argp, fmt);
	vasprintf(&query, fmt, argp);
	va_end(argp);
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_to_matrix(query);
#else
        {apop_error(0, 'c', "apop_query_to_matrix: Apophenia was compiled without mysql support.\n");
        return 0;}
#endif
  gsl_matrix	*output = NULL;
  char		    *q2, *err=NULL;
  size_t        total_rows  = 0;
    firstcall   = 
    currentrow  = 0;
	if (db==NULL) apop_db_open(NULL);
	if (apop_opts.verbose)	printf("%s\n", query);
	q2	 = malloc(strlen(query)+300);
	apop_table_exists("apop_temp_table",1);
	sqlite3_exec(db,strcat(strcpy(q2,
		"CREATE TABLE apop_temp_table AS "),query),NULL,NULL, &err); ERRCHECK
	sqlite3_exec(db,"SELECT count(*) FROM apop_temp_table",length_callback,&total_rows, &err);
	free(query);
	ERRCHECK
    sprintf(full_divider, "^%s$", apop_opts.db_nan);
    regex           = malloc(sizeof(regex_t));
    regcomp(regex, full_divider, REG_EXTENDED+REG_ICASE);
	if (total_rows>0){
        firstcall   = 1;
		if (last_names !=NULL) 
			apop_name_free(last_names); 
		last_names = apop_name_alloc();
        tr = total_rows;  //globalize
		sqlite3_exec(db,"SELECT * FROM apop_temp_table",db_to_table,&output, &err); ERRCHECK
	}
	assert(apop_table_exists("apop_temp_table",0));
	sqlite3_exec(db,"DROP TABLE apop_temp_table",NULL,NULL, &err);  ERRCHECK
    regfree(regex);
    free(regex);
	free(q2);
	return output;
}

/** Queries the database, and dumps the first column of the result into a gsl_vector.
\param fmt	A string holding an \ref sql "SQL" query.
\param ...	Your query may be in <tt>printf</tt> form. See \ref apop_query for an example.
\return		A <tt>gsl_vector</tt> holding the first column of the returned matrix. Thus, if your query returns
multiple lines, you will get no warning, and the function will return
the first in the list.

If the query returns no columns at all, the function returns <tt>NULL</tt>.
*/
gsl_vector * apop_query_to_vector(const char * fmt, ...){
  char		*query;
  va_list	argp;
	va_start(argp, fmt);
	vasprintf(&query, fmt, argp);
	va_end(argp);

    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_to_vector(query);
#else
        {apop_error(0, 'c', "apop_query_to_vector: Apophenia was compiled without mysql support.\n");
        return 0;}
#endif
  gsl_matrix	*m=NULL;
  gsl_vector  *out;
	if (db==NULL) apop_db_open(NULL);
	m	= apop_query_to_matrix(query);
	if (m==NULL){
        if (apop_opts.verbose)
		    printf("apop, %s, %i: Query turned up a blank table. Returning NULL.\n", __FILE__, __LINE__);
		return NULL;
	} //else
    out = gsl_vector_alloc(m->size1);
	gsl_matrix_get_col(out, m, 0);
	gsl_matrix_free(m);
	return out;

}

/** Queries the database, and dumps the result into a single double-precision floating point number.
\param fmt	A string holding an \ref sql "SQL" query.
\param ...	Your query may be in <tt>printf</tt> form. See \ref apop_query for an example.
\return		A double, actually. This calls \ref apop_query_to_matrix and returns
the (0,0)th element of the returned matrix. Thus, if your query returns
multiple lines, you will get no warning, and the function will return
the first in the list (which is not always well-defined).

If the query returns no rows at all, the function returns <tt>NAN</tt>.
*/
double apop_query_to_float(const char * fmt, ...){
  char		*query;
  va_list	argp;
	va_start(argp, fmt);
	vasprintf(&query, fmt, argp);
	va_end(argp);
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_to_float(query);
#else
        {apop_error(0, 'c', "apop_query_to_float: Apophenia was compiled without mysql support.\n");
        return 0;}
#endif
#ifdef HAVE_LIBSQLITE3
  gsl_matrix	*m=NULL;
  double		out;
	if (db==NULL) apop_db_open(NULL);
	m	= apop_query_to_matrix(query);
	if (m==NULL){
        if (apop_opts.verbose)
		    printf("apop, %s, %i: Query turned up a blank table. Returning NAN.\n", __FILE__, __LINE__);
		return GSL_NAN;
	} //else
	out	= gsl_matrix_get(m, 0, 0);
	gsl_matrix_free(m);
	return out;
#endif
}

/** This function returns an \ref apop_name structure with the column
names from the last <tt>apop_query_...</tt> . Since only the names from
the last query are saved, you will want to use this immediately
after your query.  */
apop_name * apop_db_get_names(void){
  apop_name   *out    = apop_name_copy(last_names);
    return out; 
}

/** Queries the database, and dumps the result into an \ref apop_data set.

\param fmt 	A string holding an \ref sql "SQL" query.
Your query may be in <tt>printf</tt> form. See \ref apop_query for an example.

\return
An \ref apop_data set, which you passed in declared but not allocated.

Blanks in the database (i.e., <tt> NULL</tt>s) and elements that match \ref apop_opts.db_nan
are filled with <tt>NAN</tt>s in the matrix.

If \ref apop_opts.db_name_column is set (it defaults to being "row_name"),
and the name of a column matches the name, then the row names are read from that column.

\bug Currently, this is but a wrapper for \ref apop_query_to_matrix,
meaning that only numerical results are returned. If you want
non-numeric data, try \code mydata->text  = apop_query_to_text("select ...");\endcode. 
*/ 
apop_data * apop_query_to_data(const char * fmt, ...){
  gsl_matrix	*m=NULL;
  va_list		argp;
  char		    *query;
  apop_data	    *out;
	va_start(argp, fmt);
	vasprintf(&query, fmt, argp);
	va_end(argp);
	if (apop_opts.verbose)	
        printf("%s\n", query);
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_to_data(query);
#else
        {apop_error(0, 'c', "apop_query_to_data: Apophenia was compiled without mysql support.\n");
        return 0;
        }
#endif

#ifdef HAVE_LIBSQLITE3
    data_or_matrix  = 'd';
	m	        = apop_query_to_matrix(query);
    if (!m) return NULL;
    data_or_matrix  = 'm';
    out         = apop_matrix_to_data(m);
    //replace name struct allocated in apop_matrix_to_data with the
    //actual names.
    apop_name_free(out->names);
    out->names  = apop_db_get_names();
	return out;
#endif
}


/** Query data to an \c apop_data set, but a mix of names, vectors, matrix elements, and text.

If you are querying to a matrix and maybe a name, use \c
apop_query_to_data (and set \c apop_opts.db_name_column if desired). But
if your data is a mix of text and numbers, use this.

The first argument is a character string consisting of the letters \c
nvmt, one for each column of the SQL output, indicating whether the
column is a name, vector, matrix colum, or text column. You can have only
one n and v. 

If the query produces more columns than there are elements in the column
specification, then the remainder are dumped into the text section. If
there are fewer columns produced than given in the spec, the additional
elements will be allocated but not filled (i.e., they are uninitialized
and will have garbage).

The 'n' character indicates row, meaning that \c apop_opts.db_name_column is ignored).

As with the other \c apop_query_to_... functions, the query can include printf-style format specifiers.
*/
apop_data * apop_query_to_mixed_data(const char *typelist, const char * fmt, ...){
  va_list   argp;
  char      *query;
    va_start(argp, fmt);
    vasprintf(&query, fmt, argp);
    if (apop_opts.verbose) {printf("\n%s\n",query);}
    va_end(argp);
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        {apop_error(0, 'c', "%s: Sorry, this function has only been written for SQLITE so far.\n", __func__);
        return 0;}
        //return apop_mysql_query_to_text(query);
#else
        {apop_error(0, 'c', "%s: Apophenia was compiled without mysql support.\n", __func__);
        apop_error(0, 'c', "%s: Also, this function has only been written for SQLITE so far.\n", __func__);
        return 0;}
#endif
#ifdef HAVE_LIBSQLITE3
        return apop_sqlite_multiquery(typelist, query);
#else
        {apop_error(0, 'c', "%s: Apophenia was compiled without SQLite support.\n", __func__);
        return NULL; }
#endif
}

/** Dump a <tt>gsl_matrix</tt> into the database. This function is
basically preempted by \ref apop_matrix_print. Use that one; this may soon no longer be available.

\param data 	The name of the matrix
\param tabname	The name of the db table to be created
\param headers	A list of column names. If <tt>NULL</tt>, then the columns will be named <tt>c1</tt>, <tt>c2</tt>, <tt>c3</tt>, &c.
 \ingroup conversions
*/
int apop_matrix_to_db(gsl_matrix *data, char *tabname, char **headers){
  int		i,j; 
  double    v;
  int		ctr		= 0;
  int		batch_size	= 100;
  char		*q 		= malloc(1000);
	if (db==NULL) apop_db_open(NULL);
	sprintf(q, "create table %s (", tabname);
	for(i=0;i< data->size2; i++){
		q	=realloc(q,strlen(q)+1000);
		if(headers == NULL) 	sprintf(q, "%s\n c%i", q,i);
		else			sprintf(q, "%s\n %s ", q,headers[i]);
		if (i< data->size2-1) 	sprintf(q, "%s,",q);
		else			sprintf(q,"%s);  begin;",q);
	}
	for(i=0;i< data->size1; i++){
		q	=realloc(q,strlen(q)+(1+data->size2)*1000);
		sprintf(q,"%s \n insert into %s values(",q,tabname);
		for(j=0;j< data->size2; j++){
            v   =gsl_matrix_get(data,i,j);
            if (gsl_isnan(v))
			    sprintf(q,"%s NULL%s ",
                    q, 
                    j < data->size2 -1 ? "," : ");");
            else
			    sprintf(q,"%s %g%s ",
                    q, v,
                    j < data->size2 -1 ? "," : ");");
        }
	ctr++;
	if(ctr==batch_size) {
			apop_query("%s commit;",q);
			ctr = 0;
			sprintf(q,"begin; \n insert into %s values(",tabname);
		}
	}
    if (ctr>0) 
        apop_query("%s commit;",q);
	free(q);
	return 0;
}


/** Dump an \ref apop_data set into the database.

This function is basically preempted by \ref apop_data_print. Use that
one; this may soon no longer be available.

Column names are inserted if there are any. If there are, all dots
are converted to underscores. 
Otherwise, the columns will be named <tt>c1</tt>, <tt>c2</tt>, <tt>c3</tt>, &c.

If \ref apop_opts.db_name_column is not blank (the default is "row_name"),
then a so-named column is created, and the row names are placed there.

\param set 	    The name of the matrix
\param tabname	The name of the db table to be created
\ingroup apop_data
\todo add text names.
 \ingroup conversions
*/
int apop_data_to_db(apop_data *set, char *tabname){
  int		i,j; 
  int		ctr		    = 0;
  int		batch_size	= 100;
  double    v;
  char		*q 		    = malloc(1000);
  int       use_row= strlen(apop_opts.db_name_column) &&set->names->rowct == set->matrix->size1;

    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
    {
        sprintf(q, "create table %s (", tabname);
        if (use_row) {
            sprintf(q, "%s\n %s varchar(1000), \n", q, apop_opts.db_name_column);
            q	=realloc(q,strlen(q)+1000);
        }
        for(i=0;i< set->matrix->size2; i++){
            q	=realloc(q,strlen(q)+1000);
            if(set->names->colct <= i) 	sprintf(q, "%s\n c%i double", q,i);
            else			sprintf(q, "%s\n %s  double", q,apop_strip_dots(set->names->column[i],'d'));
            if (i< set->matrix->size2-1 || set->textsize[1]) 	
                strcat(q, ", ");
        }
        for(i=0;i< set->textsize[1]; i++){
            q	=realloc(q,strlen(q)+1000);
            if(set->names->textct <= i) 	sprintf(q, "%s\n tc%i varchar(1000)", q,i);
            else			sprintf(q, "%s\n %s  varchar(1000)", q,apop_strip_dots(set->names->text[i],'d'));
            if (i< set->textsize[1]-1) 	strcat(q, ", ");
        }
        strcat(q,"); ");
        apop_query(q);
        sprintf(q, " ");
    }
#else 
        {apop_error(0, 'c', "apop_data_to_db: Apophenia was compiled without mysql support.\n");
        return 1;
        }
#endif
    else {
#ifdef HAVE_LIBSQLITE3
        if (db==NULL) apop_db_open(NULL);
        sprintf(q, "create table %s (", tabname);
        if (use_row) {
            sprintf(q, "%s\n %s, \n", q, apop_opts.db_name_column);
            q	=realloc(q,strlen(q)+1000);
        }
        for(i=0;i< set->matrix->size2; i++){
            q	=realloc(q,strlen(q)+1000);
            if(set->names->colct <= i) 	sprintf(q, "%s\n c%i", q,i);
            else			sprintf(q, "%s\n \"%s\" ", q,apop_strip_dots(set->names->column[i],'d'));
            if (i< set->matrix->size2-1 || set->textsize[1]) 	
                sprintf(q, "%s,",q);
        }
        for(i=0; i< set->textsize[1]; i++){
            q	=realloc(q,strlen(q)+1000);
            if(set->names->textct <= i) 	sprintf(q, "%s\n tc%i ", q,i);
            else			sprintf(q, "%s\n %s ", q,apop_strip_dots(set->names->text[i],'d'));
            if (i< set->textsize[1]-1) 	strcat(q, ", ");
        }
        sprintf(q,"%s);  begin;",q);
#else
        {apop_error(0, 'c', "apop_data_to_db: Apophenia was compiled without SQLite support.\n");
        return 1;}
#endif
    }
	for(i=0;i< set->matrix->size1; i++){
		q	=realloc(q,strlen(q)+(1+set->matrix->size2)*batch_size*1000);
		sprintf(q, "%s \n insert into %s values(",q, tabname);
        if (use_row)
			sprintf(q,"%s \'%s\', ",q, set->names->row[i]);
		for(j=0;j< set->matrix->size2; j++){
            v   =gsl_matrix_get(set->matrix,i,j);
            if (gsl_isnan(v))
			    sprintf(q,"%s NULL%s ",
                    q, 
                    (j < set->matrix->size2 -1 || set->textsize[1]) ? "," : " ");
            else
			    sprintf(q,"%s %g%s ",
                    q, v,
                    (j < set->matrix->size2 -1 || set->textsize[1]) ? "," : " ");
        }
		for(j=0;j< set->textsize[1]; j++)
			sprintf(q,"%s \'%s\' %c ",q, set->text[i][j], j < set->textsize[1]-1 ? ',' : ' ');
        sprintf(q,"%s);",q);
		ctr++;
		if(ctr==batch_size || apop_opts.db_engine == 'm') {
		    if(apop_opts.db_engine == 'm')   apop_query(q);
            else                                apop_query("%s commit;",q);
            ctr = 0;
		    if(apop_opts.db_engine == 'm')   
                sprintf(q," ");
            else                            
                sprintf(q,"begin; \n insert into %s values(",tabname);
			}
	}
    if ( !(apop_opts.db_engine == 'm') && ctr>0) 
        apop_query("%s commit;",q);
	free(q);
	return 0;
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
	if (db_file !=NULL)
		apop_query("attach database \"%s\" as merge_me;", db_file);
	apop_data *d = apop_query_to_text("select name from sqlite_master where name == \"%s\";", tabname);
	if (!d){	//just import table
		if (apop_opts.verbose)	printf("adding in %s\n", tabname);
		apop_query("create table main.%s as select * from merge_me.%s;", tabname, tabname);
	}
	else	{			//merge tables.
        apop_data_free(d);
		if (apop_opts.verbose)	printf("merging in %s\n", tabname);
		apop_query("insert into main.%s select * from merge_me.%s;", tabname, tabname);
	}
	if (db_file !=NULL)
		apop_query("detach database merge_me;");
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
  apop_data	*tab_list;
  int		i;
	apop_query("attach database \"%s\" as merge_me;", db_file);
	tab_list= apop_query_to_text("select name from merge_me.sqlite_master where type==\"table\";");
	for(i=0; i< tab_list->textsize[0]; i++)
		apop_db_merge_table(NULL, (tab_list->text)[i][0]);
	apop_query("detach database merge_me;");
	apop_data_free(tab_list);
}
                                                                                                                               
////////////////////////////////////////////////
// Part three: some stats wrappers
////////////////////////////////////////////////

/** Do a t-test entirely inside the database.
  Returns only the two-tailed p-value.
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
	return apop_two_tailify(gsl_cdf_tdist_P(stat, a_count+b_count-2));
}

/** Do a paired t-test entirely inside the database.
  Returns only the two-tailed p-value.
\ingroup ttest
*/
double	apop_db_paired_t_test(char * tab1, char *col1, char *col2){
  gsl_matrix	*result;
	result	= apop_query_to_matrix("select avg(%s - %s), var(%s - %s), count(*) from %s tab1", 
						   col1,col2,   col1, col2,          tab1);
  double		avg	    = gsl_matrix_get(result, 0, 0),
		        var	    = gsl_matrix_get(result, 0, 1),
		        count	= gsl_matrix_get(result, 0, 2),
		        stat	= avg/ sqrt(var/(count-1));
	return 2*GSL_MIN(gsl_cdf_tdist_P(stat, count-1),gsl_cdf_tdist_Q(stat, count-1));
}
