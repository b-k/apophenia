/** \file apop_db.c	An easy front end to SQLite. Includes a few nice
features like a variance, skew, and kurtosis aggregator for SQL. */
/* Copyright (c) 2006--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "db.h"
#include "stats.h"	    //t_dist
#include "conversions.h"	//apop_strip_dots
#include "variadic.h"
#include <config.h>

/** Here are where the options are initially set. */
apop_opts_type apop_opts	= 
          { .verbose=0,                    .output_type = 'f',
            .output_pipe = NULL,           .output_delimiter ="\t", 
            .output_append = 0,            .input_delimiters = "| ,\t", 
            .db_name_column = "row_names", .db_nan = "NaN", 
            .db_engine = '\0',             .db_user = "\0", 
            .db_pass = "\0",               .thread_count = 1,
            .rng_seed = 479901,            .version = X.XX };

#ifdef HAVE_LIBMYSQLCLIENT
#include "apop_db_mysql.c"
#endif

static gsl_rng* db_rng  = NULL;     //the RNG for the RNG function.

#ifdef HAVE_LIBSQLITE3 
#include "apop_db_sqlite.c"
#endif

//This macro declares the query string and fills it from the printf part of the call.
#define Fillin(query, fmt)          \
  char		*query;                 \
  va_list   argp;                   \
	va_start(argp, fmt);            \
	vasprintf(&query, fmt, argp);   \
	va_end(argp);                   \
	if (apop_opts.verbose) {printf("\n%s\n",query);} 

typedef struct {
    int         firstcall;
    size_t      currentrow;
    regmatch_t  result[3];
    regex_t     *regex;
    apop_data   *outdata;
} callback_t;

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

When you are done doing your database manipulations, be sure to call \ref apop_db_close if writing to disk.

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

/*
int apop_db_open(char *filename){
    int check_env = 0;
#ifdef HAVE_LIBSQLITE3
    if (!db) //check the environment.
        check_env = 1;
#endif
#ifdef HAVE_LIBMYSQLCLIENT
   if(!mysql_db)  
       check_env = 1;
#endif


    if (check_env && getenv("APOP_DB_ENGINE") && !strcasecmp(getenv("APOP_DB_ENGINE"), "mysql"))
        apop_opts.db_engine = 'm';

    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_db_open(filename);
#else
        apop_assert(0, 1, 0, 'c', "Apophenia was compiled without mysql support.\n");
#endif
#ifdef HAVE_LIBSQLITE3
    return apop_sqlite_db_open(filename);
#else
    apop_assert(0, 1, 0, 'c', "Apophenia was compiled without SQLite support.\n");
#endif
}
*/


int 		isthere;//used for the next two fns only.

static int tab_exists_callback(void *in, int argc, char **argv, char **whatever){
  char *q	= in;
	if (!strcmp(argv[argc-1],q))
		isthere=1;
	return 0;
}

/** Check for the existence of a table, and maybe delete it.

Recreating a table which already exists can cause errors, so it is good practice to check for existence first.
Also, this is the stylish way to delete a table, since just calling <tt>"drop table"</tt> will give you an error if the table doesn't exist.

\param name 	the table name (no default)
\param remove 'd'	==>delete table so it can be recreated in main.<br>
		'n'	==>no action. return error so program can continue. (default)
\return
0 = table does not exist<br>
1 = table was found, and if whattodo==1, has been deleted

This function uses the \ref designated syntax for inputs.
\ingroup db
*/
APOP_VAR_HEAD int apop_table_exists(char *name, char remove){
    char *apop_varad_var(name, NULL)
    apop_assert(name, 0, 0, 's', "You gave me a NULL table name.");
    char apop_varad_var(remove, 'n')
    return apop_table_exists_base(name, remove);
APOP_VAR_END_HEAD
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_table_exists(name, remove);
#else
        apop_assert(0, 0, 0, 'c', "Apophenia was compiled without mysql support.\n");
#endif
#ifdef HAVE_LIBSQLITE3
  char 		*err, q2[10000];
	isthere=0;
	if (db==NULL) return 0;
	sqlite3_exec(db, "select name from sqlite_master where type='table'",tab_exists_callback,name, &err); 
	ERRCHECK
	if ((remove==1|| remove=='d') && isthere){
		sqlite3_exec(db,strcat(strcpy(q2, "DROP TABLE "),name),NULL,NULL, &err); 
        ERRCHECK
    }
	return isthere;
#endif
}

/**
Closes the database on disk. If you opened the database with
\c apop_db_open(NULL), then this is basically optional.

\param vacuum 
'v': vacuum---do clean-up to minimize the size of the database on disk.<br>
'q': Don't bother; just close the database. (default = 'q')

This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD int apop_db_close(char vacuum){
    char apop_varad_var(vacuum, 'q')
    return apop_db_close_base(vacuum);
APOP_VAR_END_HEAD
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        {apop_mysql_db_close(0);
        return 0;}
#else
        apop_assert(0, 0, 0, 'c', "Apophenia was compiled without mysql support.\n");
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
        apop_assert(0, 0, 0, 'c', "Apophenia was compiled without SQLite support.\n");
#endif
    }
}

/** \defgroup queries Queries
 
  These functions query the database, and most return a value for use on the C-side.

In all cases, your query may be in <tt>printf</tt> form. For example:
\code
char tabname[] = "demographics";
char colname[] = "heights";
int min_height = 175;
apop_query("select %s from %s where %s > %i", colname, tabname, colname, min_height);
\endcode

\li Blanks in the database (i.e., <tt> NULL</tt>s) and elements that match \ref apop_opts_type "apop_opts.db_nan"
are filled with <tt>NAN</tt>s in the matrix.

  \{
  */

/** Send a query to the database, return nothing 
\param fmt A <tt>printf</tt>-style \ref sql "SQL" query.
*/
int apop_query(const char *fmt, ...){
  char 		*err;
  Fillin(q, fmt)
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        {apop_assert(mysql_db, 1, 0, 'c', "No mySQL database is open.");
        apop_mysql_query(q);}
#else
        apop_assert(0, 1, 0, 'c', "Apophenia was compiled without mysql support.\n");
#endif
    else 
#ifdef HAVE_LIBSQLITE3
        {if (!db) apop_db_open(NULL);
        sqlite3_exec(db, q, NULL,NULL, &err);
	    ERRCHECK
        }
#else
        apop_assert(0, 1, 0, 'c', "Apophenia was compiled without SQLite support.\n");
#endif
	free(q);
	return 1;
}

/** Dump the results of a query into an array of strings.

\return		An \ref apop_data structure with the <tt>text</tt>
element filled. Notice that this is always a 2-D array, even if the query
returns a single column. In that case, use <tt>returned_tab->text[i][0]</tt>
to refer to row <tt>i</tt>.

For example,
the following function will list the tables in a database (much like you could do from the command line using <tt>sqlite3 dbname.db ".table"</tt>).

\include ls_tables.c
*/
apop_data * apop_query_to_text(const char * fmt, ...){
  Fillin(query, fmt)
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_core(query, process_result_set_chars);
#else
        apop_assert(0, NULL, 0, 'c', "Apophenia was compiled without mysql support.\n");
#endif
#ifdef HAVE_LIBSQLITE3
    return apop_sqlite_query_to_text(query);
#else
    apop_assert(0, NULL, 0, 'c', "Apophenia was compiled without SQLite support.\n");
#endif
}

//apop_query_to_matrix callback.
static int db_to_table(void *qinfo, int argc, char **argv, char **column){
  int		jj, i, ncfound = 0;
  callback_t *qi= qinfo ;
    if (qi->firstcall){
        qi->firstcall   --;
        namecol     = -1;
        for(i=0; i<argc; i++)
            if (!strcmp(column[i], apop_opts.db_name_column)){
                namecol = i;
                ncfound = 1;
                break;
            }
	    qi->outdata		= apop_data_alloc(0, 1, argc-ncfound);
        for(i=0; i<argc; i++)
            if (namecol != i)
                apop_name_add(qi->outdata->names, column[i], 'c');
    } else 
        apop_matrix_realloc(qi->outdata->matrix, qi->currentrow+1, qi->outdata->matrix->size2);
	if (argv !=NULL){
        ncfound =0;
		for (jj=0;jj<argc;jj++)
            if (jj != namecol){
                if (!argv[jj] || !regexec(qi->regex, argv[jj], 1, qi->result, 0))
				    gsl_matrix_set(qi->outdata->matrix,qi->currentrow,jj-ncfound, GSL_NAN);
			    else
				    gsl_matrix_set(qi->outdata->matrix,qi->currentrow,jj-ncfound, atof(argv[jj]));
            } else {
                apop_name_add(qi->outdata->names, argv[jj], 'r');
                ncfound = 1;
            }
		(qi->currentrow)++;
	}
	return 0;
}

/** Queries the database, and dumps the result into an \ref apop_data set.

If \ref apop_opts_type "apop_opts.db_name_column" is set (it defaults to being "row_names"),
and the name of a column matches the name, then the row names are read from that column.
*/ 
apop_data * apop_query_to_data(const char * fmt, ...){
  Fillin(query, fmt)
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_core(query, process_result_set_data);
#else
        apop_assert(0, 0, 0, 'c', "Apophenia was compiled without mysql support.\n");
#endif

#ifdef HAVE_LIBSQLITE3
  char		    *err=NULL;
  char        full_divider[103];
  callback_t  qinfo;
    qinfo.firstcall   = 1;
    qinfo.currentrow  = 0;
    qinfo.outdata     = NULL;
	if (db==NULL) apop_db_open(NULL);
	if (apop_opts.verbose)	printf("%s\n", query);
    sprintf(full_divider, "^%s$", apop_opts.db_nan);
    qinfo.regex           = malloc(sizeof(regex_t));
    regcomp(qinfo.regex, full_divider, REG_EXTENDED+REG_ICASE);
    sqlite3_exec(db, query,db_to_table,&qinfo, &err); ERRCHECK
    regfree(qinfo.regex);
    free(qinfo.regex);
    free (query);
	return qinfo.outdata;
#endif
}

/** Queries the database, and dumps the result into a matrix.  */
gsl_matrix * apop_query_to_matrix(const char * fmt, ...){
    Fillin(query, fmt)
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_core(query, process_result_set_matrix);
#else
        apop_assert(0, 0, 0, 'c', "Apophenia was compiled without mysql support.\n");
#endif
    apop_data * outd = apop_query_to_data(query);
    gsl_matrix *outm = NULL;
    if (outd){
        outm = outd->matrix;
        outd->matrix = NULL;
        apop_data_free(outd);
    }
    free(query);
    return outm;
}

/** Queries the database, and dumps the first column of the result into a gsl_vector.

\return		A <tt>gsl_vector</tt> holding the first column of the returned matrix. Thus, if your query returns
multiple lines, you will get no warning, and the function will return
the first in the list.

If the query returns no columns at all, the function returns \c NULL.  */
gsl_vector * apop_query_to_vector(const char * fmt, ...){
    Fillin(query, fmt)
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_core(query, process_result_set_vector);
#else
        apop_assert(0, 0, 0, 'c', "Apophenia was compiled without mysql support.\n");
#endif
  apop_data	*d=NULL;
  gsl_vector  *out;
	if (db==NULL) apop_db_open(NULL);
	d	= apop_query_to_data(query);
    apop_assert(d, NULL, 1, 'c', "Query turned up a blank table. Returning NULL.");
    //else:
    out = gsl_vector_alloc(d->matrix->size1);
	gsl_matrix_get_col(out, d->matrix, 0);
	apop_data_free(d);
    free(query);
	return out;

}

/** Queries the database, and dumps the result into a single double-precision floating point number.
\return		A double, actually. This calls \ref apop_query_to_matrix and returns
the (0,0)th element of the returned matrix. Thus, if your query returns
multiple lines, you will get no warning, and the function will return
the first in the list (which is not always well-defined).

If the query returns no rows at all, the function returns <tt>NAN</tt>.  */
double apop_query_to_float(const char * fmt, ...){
    Fillin(query, fmt)
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_to_float(query);
#else
        apop_assert(0, 0, 0, 'c', "Apophenia was compiled without mysql support.\n");
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

/** Query data to an \c apop_data set, but a mix of names, vectors, matrix elements, and text.

If you are querying to a matrix and maybe a name, use \c
apop_query_to_data (and set \ref apop_opts_type "apop_opts.db_name_column" if desired). But
if your data is a mix of text and numbers, use this.

The first argument is a character string consisting of the letters \c nvmtw, one for each column of the SQL output, indicating whether the column is a name, vector, matrix colum, text column, or weight vector. You can have only one n, v, and w. 

If the query produces more columns than there are elements in the column specification, then the remainder are dumped into the text section. If there are fewer columns produced than given in the spec, the additional elements will be allocated but not filled (i.e., they are uninitialized and will have garbage).

The 'n' character indicates row, meaning that \ref apop_opts_type "apop_opts.db_name_column" is ignored).

As with the other \c apop_query_to_... functions, the query can include printf-style format specifiers.
*/
apop_data * apop_query_to_mixed_data(const char *typelist, const char * fmt, ...){
    Fillin(query, fmt)
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
    apop_assert(0, 0, 0, 'c', "Apophenia was compiled without SQLite support.\n");
#endif
}

/** \} end query group. */

/* Convenience function for extending a string. 
 asprintf(%q, "%s and stuff", q);
 gives you a memory leak. This takes care of that.
 */
void qxprintf(char **q, char *format, ...){
    va_list ap;
    char *r = *q;
    va_start(ap, format);
    vasprintf(q, format, ap);
    va_end(ap);
    free(r);
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
	asprintf(&q, "create table %s (", tabname);
	for(i=0;i< data->size2; i++){
		if(headers == NULL) 	qxprintf(&q, "%s\n c%i", q,i);
		else			qxprintf(&q, "%s\n %s ", q,headers[i]);
		if (i< data->size2-1) 	qxprintf(&q, "%s,",q);
		else			qxprintf(&q,"%s);  begin;",q);
	}
	for(i=0;i< data->size1; i++){
		qxprintf(&q,"%s \n insert into %s values(",q,tabname);
		for(j=0;j< data->size2; j++){
            v   =gsl_matrix_get(data,i,j);
            if (gsl_isnan(v))
			    qxprintf(&q,"%s NULL%s ",
                    q, 
                    j < data->size2 -1 ? "," : ");");
            else
			    qxprintf(&q,"%s %g%s ",
                    q, v,
                    j < data->size2 -1 ? "," : ");");
        }
	ctr++;
	if(ctr==batch_size) {
			apop_query("%s commit;",q);
			ctr = 0;
			qxprintf(&q,"begin; \n insert into %s values(",tabname);
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
Otherwise, the columns will be named \c c1, \c c2, \c c3, &c.

If \ref apop_opts_type "apop_opts.db_name_column" is not blank (the default is "row_name"),
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
  int       use_row= strlen(apop_opts.db_name_column) 
                && ((set->matrix && set->names->rowct == set->matrix->size1)
                    || (set->vector && set->names->rowct == set->vector->size));

    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
    {
        int comma;
        asprintf(&q, "create table %s (", tabname);
        if (use_row) 
            qxprintf(&q, "%s\n %s varchar(1000), \n", q, apop_opts.db_name_column);
        if (set->vector){
            comma = (set->matrix || set->textsize[1]);
            if(!set->names->vector) 
                qxprintf(&q, "%s\n vector double%s", q, comma ? ", " : " ");
            else
                qxprintf(&q, "%s\n \"%s\" double%s", q,apop_strip_dots(set->names->vector,'d'), comma ? ", " : " ");
        }
        if (set->matrix)
            for(i=0;i< set->matrix->size2; i++){
                comma = (i< set->matrix->size2-1 || set->textsize[1]);
                if(set->names->colct <= i) 
                    qxprintf(&q, "%s\n c%i double%s", q,i, comma ? ", " : " ");
                 else
                    qxprintf(&q, "%s\n %s  double%s", q,apop_strip_dots(set->names->column[i],'d'), comma ? ", " : " ");
            }
        for(i=0;i< set->textsize[1]; i++){
            comma = (i< set->textsize[1]-1);
            if (set->names->textct <= i)
                qxprintf(&q, "%s\n tc%i varchar(1000)%s", q,i, comma ? ", " : " ");
            else
                qxprintf(&q, "%s\n %s  varchar(1000)%s", q,apop_strip_dots(set->names->text[i],'d'), comma ? ", " : " ");
        }
        apop_query("%s); ", q);
        sprintf(q, " ");
    }
#else 
        apop_assert(0, 1, 0, 'c', "Apophenia was compiled without mysql support.\n");
#endif
    else
#ifdef HAVE_LIBSQLITE3
    {
        if (db==NULL) apop_db_open(NULL);
        asprintf(&q, "create table %s (", tabname);
        if (use_row) {
            qxprintf(&q, "%s\n %s, \n", q, apop_opts.db_name_column);
        }
        if (set->vector){
            if(!set->names->vector) 	qxprintf(&q, "%s\n vector", q);
            else			qxprintf(&q, "%s\n \"%s\" ", q,apop_strip_dots(set->names->vector,'d'));
            if (set->matrix || set->textsize[1]) 	
                qxprintf(&q, "%s,",q);
        }
        if (set->matrix)
            for(i=0;i< set->matrix->size2; i++){
                if(set->names->colct <= i) 	qxprintf(&q, "%s\n c%i", q,i);
                else			qxprintf(&q, "%s\n \"%s\" ", q,apop_strip_dots(set->names->column[i],'d'));
                if (i< set->matrix->size2-1 || set->textsize[1]) 	
                    qxprintf(&q, "%s,",q);
            }
        for(i=0; i< set->textsize[1]; i++){
            if(set->names->textct <= i) 	qxprintf(&q, "%s\n tc%i ", q,i);
            else			qxprintf(&q, "%s\n %s ", q, apop_strip_dots(set->names->text[i],'d'));
            if (i< set->textsize[1]-1) 	qxprintf(&q, "%s, ", q);
        }
        qxprintf(&q,"%s);  begin;",q);
    }
#else
        apop_assert(1, 0, 0, 'c', "Apophenia was compiled without SQLite support.\n");
#endif
    int lim = set->vector ? set->vector->size : set->matrix->size1;
	for(i=0; i< lim; i++){
		qxprintf(&q, "%s \n insert into %s values(",q, tabname);
        if (use_row)
			qxprintf(&q,"%s \'%s\', ",q, set->names->row[i]);
        if (set->vector){
            v   =gsl_vector_get(set->vector,i);
            if (gsl_isnan(v))
                qxprintf(&q,"%s NULL%s ",
                    q, 
                    (set->matrix || set->textsize[1]) ? "," : " ");
            else
                qxprintf(&q,"%s %g%s ",
                    q, v,
                    (set->matrix || set->textsize[1]) ? "," : " ");
        }
        if (set->matrix)
            for(j=0;j< set->matrix->size2; j++){
                v   =gsl_matrix_get(set->matrix,i,j);
                if (gsl_isnan(v))
                    qxprintf(&q,"%s NULL%s ",
                        q, 
                        (j < set->matrix->size2 -1 || set->textsize[1]) ? "," : " ");
                else if (isinf(v)==1)
                    qxprintf(&q,"%s 'inf'%s ",
                        q, 
                        (j < set->matrix->size2 -1 || set->textsize[1]) ? "," : " ");
                else if (isinf(v)==-1)
                    qxprintf(&q,"%s '-inf'%s ",
                        q, 
                        (j < set->matrix->size2 -1 || set->textsize[1]) ? "," : " ");
                else
                    qxprintf(&q,"%s %g%s ",
                        q, v,
                        (j < set->matrix->size2 -1 || set->textsize[1]) ? "," : " ");
            }
		for(j=0;j< set->textsize[1]; j++)
			qxprintf(&q,"%s \'%s\' %c ",q, set->text[i][j], j < set->textsize[1]-1 ? ',' : ' ');
        qxprintf(&q,"%s);",q);
		ctr++;
		if(ctr==batch_size || apop_opts.db_engine == 'm') {
		    if(apop_opts.db_engine == 'm')   apop_query(q);
            else                                apop_query("%s commit;",q);
            ctr = 0;
		    /*if(apop_opts.db_engine == 'm')   
                qxprintf(&q,"%s ", q);
            else                            */
		    if(apop_opts.db_engine != 'm')   
                qxprintf(&q,"begin; \n");
			}
	}
    if ( !(apop_opts.db_engine == 'm') && ctr>0) 
        apop_query("%s commit;",q);
	free(q);
	return 0;
}

/** Merge a single table from a database on the hard drive with the database currently open.

\param db_file	The name of a file on disk. [default = \c NULL]
\param tabname	The name of the table in that database to be merged in. [No default]
\param inout  Do we copy data in to the currently-open main db [\c 'i'] or out to the specified auxiliary db[\c 'o']?  [default = 'i']

If the table exists in the new database but not in the currently open one,
then it is simply copied over. If there is a table with the same name in
the currently open database, then the data from the new table is inserted
into the main database's table with the same name. [The function just
calls <tt>insert into main.tab select * from merge_me.tab</tt>.]

\ingroup db
This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD void apop_db_merge_table(char *db_file, char *tabname, char inout){
    char * apop_varad_var(tabname, NULL);
    Apop_assert_void(tabname, 0,'s', "I need a non-NULL tabname");
    char * apop_varad_var(db_file, NULL);
    char apop_varad_var(inout, 'i');
    apop_db_merge_table_base(db_file, tabname,inout);
APOP_VAR_ENDHEAD
    char maine[] = "main";
    char merge_me[] = "merge_me";
    char *from = inout == 'i' ? merge_me : maine;
    char *to = inout == 'i' ? maine : merge_me ;
	if (db_file !=NULL)
		apop_query("attach database \"%s\" as merge_me;", db_file);
	int d = apop_query_to_float("select count(*) from %s.sqlite_master where name == \"%s\";", to, tabname);
	if (!d){	//just import table
		if (apop_opts.verbose)	printf("adding in %s\n", tabname);
		apop_query("create table %s.%s as select * from %s.%s;", to, tabname, from, tabname);
	}
	else	{			//merge tables.
		if (apop_opts.verbose)	printf("merging in %s\n", tabname);
		apop_query("insert into %s.%s select * from %s.%s;", to, tabname, from, tabname);
	}
	if (db_file !=NULL)
		apop_query("detach database merge_me;");
}

/** Merge a database on the hard drive with the database currently open.

\param db_file	The name of a file on disk. [No default; can't be \c NULL]
\param inout  Do we copy data in to the currently-open main db [\c 'i'] or out to the specified auxiliary db[\c 'o']?  [default = 'i']

If a table exists in the new database but not in the currently open one,
then it is simply copied over. If there are  tables with the same name
in both databases, then the data from the new table is inserted into
the main database's table with the same name. [The function just calls
<tt>insert into main.tab select * from merge_me.tab</tt>.]

\todo This is sqlite-only right now. 

This function uses the \ref designated syntax for inputs.
\ingroup db
*/
APOP_VAR_HEAD void apop_db_merge(char *db_file, char inout){
    char * apop_varad_var(db_file, NULL);
    Apop_assert_void(db_file, 0, 's', "This function copies from a named database file to the currently in-memory database. You need to give me the name of that named db.")
    char apop_varad_var(inout, 'i');
    apop_db_merge_base(db_file, inout);
APOP_VAR_ENDHEAD
  apop_data	*tab_list;
	apop_query("attach database \"%s\" as filedb;", db_file);
	tab_list= apop_query_to_text("select name from %s.sqlite_master where type==\"table\";"
              , inout == 'i' ? "filedb" : "main" );
    if(!tab_list) return; //No tables to merge.
	for(int i=0; i< tab_list->textsize[0]; i++)
		apop_db_merge_table(db_file, (tab_list->text)[i][0], inout);
	apop_query("detach database filedb;");
	apop_data_free(tab_list);
}
                                                                                                                               
// Some stats wrappers

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
        return	fabs(1 - (1 - gsl_cdf_tdist_P(stat, a_count+b_count-2))*2); //two-tailify a one-tailed lookup.
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
