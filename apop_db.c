/** \file apop_db.c	An easy front end to SQLite. Includes a few nice
features like a variance, skew, and kurtosis aggregator for SQL. */
/* Copyright (c) 2006--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
#include <regex.h>

/** Here are where the options are initially set. */
apop_opts_type apop_opts	= 
          { .verbose=1,                    .output_type = 'f',
            .output_pipe = NULL,           .output_delimiter ="\t", 
            .output_append = 0,            .input_delimiters = "|,\t", 
            .db_name_column = "row_names", .db_nan = "NaN", 
            .db_engine = '\0',             .db_user = "\0", 
            .db_pass = "\0",               .thread_count = 1,
            .log_file = NULL,
            .rng_seed = 479901,            .version = X.XX };

#ifdef HAVE_LIBMYSQLCLIENT
#include "apop_db_mysql.c"
#endif

#define ERRCHECK {Apop_assert_c(err==NULL, 1, 0, "%s: %s",query, err);}
#define ERRCHECK_NR {Apop_assert_c(err==NULL, NULL, 0, "%s: %s",query, err);}

static gsl_rng* db_rng  = NULL;     //the RNG for the RNG function.

#include "apop_db_sqlite.c"

//This macro declares the query string and fills it from the printf part of the call.
#define Fillin(query, fmt)          \
  char		*query;                 \
  va_list   argp;                   \
	va_start(argp, fmt);            \
	vasprintf(&query, fmt, argp);   \
	va_end(argp);                   \
	Apop_notify(2, "%s", query);

typedef struct {
    int         firstcall;
    size_t      currentrow;
    regmatch_t  result[3];
    regex_t     *regex;
    apop_data   *outdata;
} callback_t;

/** Random numbers are generated inside the database using a separate RNG. This will initialize it for you, just like \ref apop_rng_alloc, except the RNG it produces is kept for internal use. If you don't call it, then it will be called at first use, with seed zero.

\param  seed    The seed. No need to get funny with it: 0, 1, and 2 will produce wholly different streams.
\return The RNG ready for your use.
\ingroup db
*/
void apop_db_rng_init(int seed){ db_rng  = apop_rng_alloc(seed); }

/**
If you want to use a database on the hard drive instead of memory, then call this once and only once before using any other database utilities. 

If you want a disposable database which you won't use after the program ends, don't bother with this function.

The trade-offs between an on-disk database and an in-memory db are as one would expect: memory is faster, but is destroyed when the program exits. SQLite includes a command line utility (<tt>sqlite3</tt>) which let you ask queries of a database on disk, which may be useful for debugging. There are also some graphical front-ends; just ask your favorite search engine for <a href="http://www.google.com/search?&q=sqlite+gui">SQLite GUI</a>.

MySQL users: either set the environment variable APOP_DB_ENGINE=mysql or set \c apop_opts.db_engine = 'm'.

The Apophenia package assumes you are only using a single SQLite database at a time; if not, the \ref apop_db_merge and \ref apop_db_merge_table functions may help.

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
int apop_db_open(char const *filename){
    if (!db) //check the environment.
#ifdef HAVE_LIBMYSQLCLIENT
       if(!mysql_db)  
#endif
        if (getenv("APOP_DB_ENGINE") && !strcasecmp(getenv("APOP_DB_ENGINE"), "mysql"))
            apop_opts.db_engine = 'm';

    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_db_open(filename);
#else
        {Apop_assert_c(0, 0, 0, "apop_db_open: Apophenia was compiled without mysql support.");}
#endif
        return apop_sqlite_db_open(filename);
}

typedef struct {
    char const *name;
    int isthere;
} tab_exists_t;

static int tab_exists_callback(void *in, int argc, char **argv, char **whatever){
    tab_exists_t *te = in;
	if (!strcmp(argv[argc-1], te->name))
		te->isthere=1;
	return 0;
}

/** Check for the existence of a table, and maybe delete it.

Recreating a table which already exists can cause errors, so it is good practice to check for existence first.  Also, this is the stylish way to delete a table, since just calling <tt>"drop table"</tt> will give you an error if the table doesn't exist.

\param name 	the table name (no default)
\param remove 'd'	==>delete table so it can be recreated in main.<br>
		'n'	==>no action. Return result so program can continue. (default)
\return
0 = table does not exist<br>
1 = table was found, and if remove=='d', has been deleted

\li In the SQLite engine, this function considers table views to be tables.

This function uses the \ref designated syntax for inputs.
\ingroup db
*/
APOP_VAR_HEAD int apop_table_exists(char const *name, char remove){
    char const *apop_varad_var(name, NULL)
    Apop_assert(name, "You gave me a NULL table name.");
    char apop_varad_var(remove, 'n')
APOP_VAR_END_HEAD
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_table_exists(name, remove);
#else
        Apop_assert(0, "Apophenia was compiled without mysql support.");
#endif
  char 		*err=NULL, *q2;
  tab_exists_t te = { .name = name };
  tab_exists_t tev = { .name = name };
	if (db==NULL) return 0;
	sqlite3_exec(db, "select name from sqlite_master where type='table'", tab_exists_callback, &te, &err); 
	sqlite3_exec(db, "select name from sqlite_master where type='view'", tab_exists_callback, &tev, &err); 
    char query[]="Selecting names from sqlite_master";//for ERRCHECK.
	ERRCHECK
	if ((remove==1|| remove=='d') && (te.isthere||tev.isthere)){
        if (te.isthere)
            asprintf(&q2, "drop table %s;", name);
        else
            asprintf(&q2, "drop view %s;", name);
		sqlite3_exec(db, q2, NULL, NULL, &err); 
        ERRCHECK
        free(q2);
    }
	return (te.isthere||tev.isthere);
}

/**
Closes the database on disk. If you opened the database with \c apop_db_open(NULL), then this is basically optional.

\param vacuum 
'v': vacuum---do clean-up to minimize the size of the database on disk.<br>
'q': Don't bother; just close the database. (default = 'q')

This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD int apop_db_close(char vacuum){
    char apop_varad_var(vacuum, 'q')
APOP_VAR_END_HEAD
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        {apop_mysql_db_close(0);
        return 0;}
#else
        Apop_assert(0, "Apophenia was compiled without mysql support.")
#endif
    else{
      char		*err;
        if (vacuum==1 || vacuum=='v') 
            sqlite3_exec(db, "VACUUM", NULL, NULL, &err);
    //	ERRCHECK
        sqlite3_close(db);
        db  = NULL;
        return 0;
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

/** Send a query to the database that returns no data.
\param fmt A <tt>printf</tt>-style SQL query.
\return 0 on success, 1 on failure.
*/
int apop_query(const char *fmt, ...){
  char 		*err=NULL;
  Fillin(query, fmt)
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        {Apop_assert_c(mysql_db, 1, 0, "No mySQL database is open.");
        apop_mysql_query(query);}
#else
        Apop_assert_c(0, 1, 0, "Apophenia was compiled without mysql support.")
#endif
    else 
        {if (!db) apop_db_open(NULL);
        sqlite3_exec(db, query, NULL,NULL, &err);
	    ERRCHECK
        }
	free(query);
	return 0;
}

/** Dump the results of a query into an array of strings.

\return		An \ref apop_data structure with the <tt>text</tt> element filled. Notice that this is always a 2-D array, even if the query returns a single column. In that case, use <tt>returned_tab->text[i][0]</tt> to refer to row <tt>i</tt>.

For example,
the following function will list the tables in a database (much like you could do from the command line using <tt>sqlite3 dbname.db ".table"</tt>).

\include ls_tables.c
*/
apop_data * apop_query_to_text(const char * fmt, ...){
    apop_data *out = NULL;
    Fillin(query, fmt)
    if (apop_opts.db_engine == 'm'){
#ifdef HAVE_LIBMYSQLCLIENT
        out = apop_mysql_query_core(query, process_result_set_chars);
#else
        Apop_assert_c(0, NULL, 0, "Apophenia was compiled without mysql support.");
#endif
    } else
        out = apop_sqlite_query_to_text(query);
    free(query);
    return out;
}

//apop_query_to_data callback.
static int db_to_table(void *qinfo, int argc, char **argv, char **column){
  int		i, ncfound = 0;
  callback_t *qi= qinfo ;
    if (qi->firstcall){
        qi->firstcall   --;
        namecol     = -1;
        for(i=0; i<argc; i++)
            if (!strcasecmp(column[i], apop_opts.db_name_column)){
                namecol = i;
                ncfound = 1;
                break;
            }
	    qi->outdata		= argc-ncfound ? apop_data_alloc(1, argc-ncfound) : apop_data_alloc( );
        for(i=0; i<argc; i++)
            if (namecol != i)
                apop_name_add(qi->outdata->names, column[i], 'c');
    } else 
        if (qi->outdata->matrix)
            apop_matrix_realloc(qi->outdata->matrix, qi->currentrow+1, qi->outdata->matrix->size2);
	if (argv !=NULL){
        ncfound =0;
		for (int jj=0;jj<argc;jj++)
            if (jj != namecol){
                double valor = 
                    !argv[jj] || apop_strcmp(argv[jj], "NULL")|| !regexec(qi->regex, argv[jj], 1, qi->result, 0)
				     ? GSL_NAN : atof(argv[jj]);
                gsl_matrix_set(qi->outdata->matrix,qi->currentrow,jj-ncfound, valor);
            } else {
                apop_name_add(qi->outdata->names, argv[jj], 'r');
                ncfound = 1;
            }
		(qi->currentrow)++;
	}
	return 0;
}

/** Queries the database, and dumps the result into an \ref apop_data set.

If \ref apop_opts_type "apop_opts.db_name_column" is set (it defaults to being "row_names"), and the name of a column matches the name, then the row names are read from that column.
*/ 
apop_data * apop_query_to_data(const char * fmt, ...){
  Fillin(query, fmt)
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_core(query, process_result_set_data);
#else
        Apop_assert_c(0, 0, 0, "Apophenia was compiled without mysql support.")
#endif

    //else
  char		    *err=NULL;
  char        full_divider[103];
  callback_t  qinfo = {.firstcall = 1,
                       .regex = malloc(sizeof(regex_t))};
	if (db==NULL) apop_db_open(NULL);
    sprintf(full_divider, "^%s$", apop_opts.db_nan);
    regcomp(qinfo.regex, full_divider, REG_EXTENDED+REG_ICASE+REG_NOSUB);
    sqlite3_exec(db, query,db_to_table,&qinfo, &err); ERRCHECK_NR
    regfree(qinfo.regex);
    free(qinfo.regex);
    free (query);
	return qinfo.outdata;
}

//These used to do more, but I'll leave them as a macro anyway in case of future expansion.
#define Store_settings  \
    int v = apop_opts.verbose; apop_opts.verbose=0;/*hack to prevent double-printing.*/ \

#define Restore_settings  \
    apop_opts.verbose=v;

/** Queries the database, and dumps the result into a matrix.

  Uses \ref apop_query_to_data and returns just the matrix part; see that function for notes.

  \li If \c apop_opts.db_name_column is set, then I'll ignore that column. It gets put into the names of the \ref apop_data set, and then thrown away when I return only the \c gsl_matrix part of that set. 
 
  \return A \c gsl_matrix.

 */
gsl_matrix * apop_query_to_matrix(const char * fmt, ...){
    Fillin(query, fmt)
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_core(query, process_result_set_matrix);
#else
        Apop_assert_c(0, 0, 0, "Apophenia was compiled without mysql support.")
#endif
    Store_settings
    apop_data * outd = apop_query_to_data("%s", query);
    Restore_settings
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

\return		A <tt>gsl_vector</tt> holding the first column of the returned matrix. Thus, if your query returns multiple lines, you will get no warning, and the function will return the first in the list.

\li Uses \ref apop_query_to_data internally, then throws away all but the first column of the matrix.
  \li If \c apop_opts.db_name_column is set, then I'll ignore that column. It gets put into the names of the \ref apop_data set, and then thrown away when I look at only the \c gsl_matrix part of that set. 

If the query returns no columns at all, the function returns \c NULL.  */
gsl_vector * apop_query_to_vector(const char * fmt, ...){
    Fillin(query, fmt)
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        return apop_mysql_query_core(query, process_result_set_vector);
#else
        Apop_assert_c(0, 0, 0, "Apophenia was compiled without mysql support.")
#endif
  apop_data	*d=NULL;
  gsl_vector  *out;
	if (db==NULL) apop_db_open(NULL);
    Store_settings
	d	= apop_query_to_data("%s", query);
    Restore_settings
    Apop_assert_c(d, NULL, 2, "Query [%s] turned up a blank table. Returning NULL.", query);
    //else:
    out = gsl_vector_alloc(d->matrix->size1);
	gsl_matrix_get_col(out, d->matrix, 0);
	apop_data_free(d);
    free(query);
	return out;

}

/** Queries the database, and dumps the result into a single double-precision floating point number.
\return		A double, actually.

\li This calls \ref apop_query_to_data and returns the (0,0)th element of the returned matrix. Thus, if your query returns multiple lines, you will get no warning, and the function will return the first in the list (which is not always well-defined; maybe use an {\em order by} clause in your query if you expect multiple lines).

  \li If \c apop_opts.db_name_column is set, then I'll ignore that column. It gets put into the names of the \ref apop_data set, and then thrown away when I look at only the \c gsl_matrix element of that set. 

If the query returns no rows at all, the function returns <tt>NAN</tt>.  */
double apop_query_to_float(const char * fmt, ...){
    double out;
    Fillin(query, fmt)
    if (apop_opts.db_engine == 'm'){
#ifdef HAVE_LIBMYSQLCLIENT
        out = apop_mysql_query_to_float(query);
#else
        apop_assert_c(0, 0, 0, "Apophenia was compiled without mysql support.")
#endif
    } else {
        apop_data	*d=NULL;
        if (db==NULL) apop_db_open(NULL);
        Store_settings
        d	= apop_query_to_data("%s", query);
        Restore_settings
        Apop_assert_c(d, GSL_NAN, 2, "Query [%s] turned up a blank table. Returning NaN.", query);
        out	= apop_data_get(d, 0, 0);
        apop_data_free(d);
    }
    free(query);
	return out;
}

/** Query data to an \c apop_data set, but a mix of names, vectors, matrix elements, and text.

If you are querying to a matrix and maybe a name, use \c
apop_query_to_data (and set \ref apop_opts_type "apop_opts.db_name_column" if desired). But
if your data is a mix of text and numbers, use this.

The first argument is a character string consisting of the letters \c nvmtw, one for each column of the SQL output, indicating whether the column is a name, vector, matrix column, text column, or weight vector. You can have only one n, v, and w. 

If the query produces more columns than there are elements in the column specification, then the remainder are dumped into the text section. If there are fewer columns produced than given in the spec, the additional elements will be allocated but not filled (i.e., they are uninitialized and will have garbage).

The 'n' character indicates row, meaning that \ref apop_opts_type "apop_opts.db_name_column" is ignored).

As with the other \c apop_query_to_... functions, the query can include printf-style format specifiers.
*/
apop_data * apop_query_to_mixed_data(const char *typelist, const char * fmt, ...){
    Fillin(query, fmt)
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
        {Apop_notify(0, "Sorry, this function has only been written for SQLITE so far.");
        return 0;}
        //return apop_mysql_query_to_text(query);
#else
        {Apop_notify(0, "Apophenia was compiled without mysql support.");
        Apop_notify(0, "Also, this function has only been written for SQLITE so far.");
        return 0;}
#endif
    //else
    apop_data *out = apop_sqlite_multiquery(typelist, query);
    free(query);
    return out;
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

static void add_a_number (char **q, char *comma, double v){
    if (gsl_isnan(v))
        qxprintf(q,"%s%c NULL ", *q, *comma);
    else if (isinf(v)==1)
        qxprintf(q,"%s%c  'inf'", *q, *comma);
    else if (isinf(v)==-1)
        qxprintf(q,"%s%c  '-inf' ", *q, *comma);
    else
        qxprintf(q,"%s%c %g ",*q ,*comma, v);
    *comma = ',';
}

/** Dump an \ref apop_data set into the database.

This function is basically preempted by \ref apop_data_print. Use that one; this may soon no longer be available.

Column names are inserted if there are any. If there are, all dots are converted to underscores.  Otherwise, the columns will be named \c c1, \c c2, \c c3, &c.

\li If \ref apop_opts_type "apop_opts.db_name_column" is not blank (the default is "row_name"), then a so-named column is created, and the row names are placed there.

\li If there are weights, they will be the last column of the table, and the column will be named "weights".

\li If the table exists; append to. If the table does not exist, create. So perhaps call \ref apop_table_exists <tt>("tabname", 'd')</tt> to ensure that the table is removed ahead of time.

\li You can also call this via \ref apop_data_print <tt>(data, "tabname", .output_type='d', .output_append='w')</tt> to overwrite a new table or with <tt>.output_append='a'</tt> to append.

\param set 	    The name of the matrix
\param tabname	The name of the db table to be created
\ingroup apop_data
\ingroup conversions
*/
void apop_data_to_db(const apop_data *set, const char *tabname, const char output_append){
    Apop_assert_c(set, , 1, "you sent me a NULL data set. Database table %s will not be created.", tabname);
    int	i,j; 
    int	ctr		    = 0;
    int	batch_size	= 100;
    char	*q 		    = malloc(1000);
    char  comma       = ' ';
    int   use_row= strlen(apop_opts.db_name_column) 
                && ((set->matrix && set->names->rowct == set->matrix->size1)
                    || (set->vector && set->names->rowct == set->vector->size));

    int tab_exists = apop_table_exists(tabname);

//To start, q will either be "begin;" if the table exists or "create table ... ; begin;" if it doesn't.
//Except mysql doesn't like transactions like this, so elide the "begin;" in that case.
    if (tab_exists && !apop_opts.db_engine == 'm')
        asprintf(&q, "begin;");
    else if (tab_exists && apop_opts.db_engine == 'm')
        asprintf(&q, " ");
    else if (apop_opts.db_engine == 'm')
#ifdef HAVE_LIBMYSQLCLIENT
    if (((output_append =='a' || output_append =='A') && apop_table_exists(tabname))){
        asprintf(&q, " ");
    else {
        asprintf(&q, "create table %s (", tabname);
        if (use_row) {
            qxprintf(&q, "%s\n %s varchar(1000)", q, apop_opts.db_name_column);
            comma = ',';
        }
        if (set->vector){
            if(!set->names->vector) 
                qxprintf(&q, "%s%c\n vector double ", q, comma);
            else
                qxprintf(&q, "%s%c\n \"%s\" double ", q,comma, apop_strip_dots(set->names->vector,'d'));
            comma = ',';
        }
        if (set->matrix)
            for(i=0;i< set->matrix->size2; i++){
                if(set->names->colct <= i) 
                    qxprintf(&q, "%s%c\n c%i double ", q, comma,i);
                 else
                    qxprintf(&q, "%s%c\n %s  double ", q, comma,apop_strip_dots(set->names->column[i],'d'));
                comma = ',';
            }
        for(i=0;i< set->textsize[1]; i++){
            if (set->names->textct <= i)
                qxprintf(&q, "%s%c\n tc%i varchar(1000) ", q, comma,i);
            else
                qxprintf(&q, "%s%c\n %s  varchar(1000) ", q, comma,apop_strip_dots(set->names->text[i],'d'));
            comma = ',';
        }
        apop_query("%s); ", q);
        sprintf(q, " ");
    }
#else 
        Apop_assert_c(0, , 0, "Apophenia was compiled without mysql support.")
#endif
    else {
        if (db==NULL) apop_db_open(NULL);
        if (((output_append =='a' || output_append =='A') && apop_table_exists(tabname)) )
            asprintf(&q, " ");
        else {
            asprintf(&q, "create table %s (", tabname);
            if (use_row) {
                qxprintf(&q, "%s\n %s", q, apop_opts.db_name_column);
                comma = ',';
            }
            if (set->vector){
                if(!set->names->vector) 	qxprintf(&q, "%s%c\n vector numeric", q, comma);
                else			qxprintf(&q, "%s%c\n \"%s\"", q, comma,apop_strip_dots(set->names->vector,'d'));
                comma = ',';
            }
            if (set->matrix)
                for(i=0;i< set->matrix->size2; i++){
                    if(set->names->colct <= i) 	
                        qxprintf(&q, "%s%c\n c%i numeric", q, comma,i);
                    else			
                        qxprintf(&q, "%s%c\n \"%s\" numeric", q, comma,apop_strip_dots(set->names->column[i],'d'));
                    comma = ',';
                }
            for(i=0; i< set->textsize[1]; i++){
                if(set->names->textct <= i) 	qxprintf(&q, "%s%c\n tc%i ", q,comma,i);
                else			qxprintf(&q, "%s%c\n %s ", q, comma, apop_strip_dots(set->names->text[i],'d'));
                comma = ',';
            }
            if (set->weights)
                qxprintf(&q, "%s%c\n \"weights\" numeric", q, comma);
            qxprintf(&q,"%s);  begin;",q);
        }
    }

    int lim = GSL_MAX(set->vector ? set->vector->size : 0,
                GSL_MAX(set->matrix ? set->matrix->size1 : 0, 
                        set->textsize[0]));
	for(i=0; i< lim; i++){
        comma = ' ';
		qxprintf(&q, "%s \n insert into %s values(",q, tabname);
        if (use_row){
            char *fixed= prep_string_for_sqlite(0, set->names->row[i]);
			qxprintf(&q, "%s %s ",q, fixed);
            free(fixed);
            comma = ',';
        }
        if (set->vector)
           add_a_number (&q, &comma, gsl_vector_get(set->vector,i));
        if (set->matrix)
            for(j=0; j< set->matrix->size2; j++)
               add_a_number (&q, &comma, gsl_matrix_get(set->matrix,i,j));
		for(j=0; j< set->textsize[1]; j++){
            char *fixed= prep_string_for_sqlite(0, set->text[i][j]);
			qxprintf(&q, "%s%c %s ",q, comma,fixed ? fixed : "''");
            free(fixed);
            comma = ',';
        }
        if (set->weights)
           add_a_number (&q, &comma, gsl_vector_get(set->weights,i));
        qxprintf(&q,"%s);",q);
		ctr++;
        apop_query("%s", q); 
        q[0]='\0';
		if(ctr==batch_size && apop_opts.db_engine != 'm') {
            ctr = 0;
            apop_query("commit;");
            qxprintf(&q,"begin; \n");
        }
	}
    if ( !(apop_opts.db_engine == 'm') && ctr>0) 
        apop_query("%s commit;",q);
	free(q);
}

/** Merge a single table from a database on the hard drive with the database currently open.

\param db_file	The name of a file on disk. [default = \c NULL]
\param tabname	The name of the table in that database to be merged in. [No default]
\param inout  Do we copy data in to the currently-open main db [\c 'i'] or out to the specified auxiliary db[\c 'o']?  [default = 'i']

If the table exists in the new database but not in the currently open one, then it is simply copied over. If there is a table with the same name in the currently open database, then the data from the new table is inserted into the main database's table with the same name. [The function just calls <tt>insert into main.tab select * from merge_me.tab</tt>.]

\ingroup db
This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD void apop_db_merge_table(char *db_file, char *tabname, char inout){
    char * apop_varad_var(tabname, NULL);
    Apop_assert_s(tabname, "I need a non-NULL tabname");
    char * apop_varad_var(db_file, NULL);
    char apop_varad_var(inout, 'i');
APOP_VAR_ENDHEAD
    char maine[] = "main";
    char merge_me[] = "merge_me";
    char *from = inout == 'i' ? merge_me : maine;
    char *to = inout == 'i' ? maine : merge_me ;
	if (db_file !=NULL)
		apop_query("attach database \"%s\" as merge_me;", db_file);
	int d = apop_query_to_float("select count(*) from %s.sqlite_master where name == \"%s\";", to, tabname);
	if (!d){	//just import table
		Apop_notify(2, "adding in %s", tabname);
		apop_query("create table %s.%s as select * from %s.%s;", to, tabname, from, tabname);
	}
	else	{			//merge tables.
		Apop_notify(2, "merging in %s", tabname);
		apop_query("insert into %s.%s select * from %s.%s;", to, tabname, from, tabname);
	}
	if (db_file !=NULL)
		apop_query("detach database merge_me;");
}

/** Merge a database on the hard drive with the database currently open.

\param db_file	The name of a file on disk. [No default; can't be \c NULL]
\param inout  Do we copy data in to the currently-open main db [\c 'i'] or out to the specified auxiliary db[\c 'o']?  [default = 'i']

If a table exists in the new database but not in the currently open one, then it is simply copied over. If there are  tables with the same name in both databases, then the data from the new table is inserted into the main database's table with the same name. [The function just calls <tt>insert into main.tab select * from merge_me.tab</tt>.]

\li This is sqlite-only; I'm not sure if it really makes much sense for mySQL.

This function uses the \ref designated syntax for inputs.
\ingroup db
*/
APOP_VAR_HEAD void apop_db_merge(char *db_file, char inout){
    char * apop_varad_var(db_file, NULL);
    Apop_assert_s(db_file, "This function copies from a named database file to the currently in-memory database. You need to give me the name of that named db.")
    char apop_varad_var(inout, 'i');
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

/** Do a \f$t\f$-test entirely inside the database.
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

/** Do a paired \f$t\f$-test entirely inside the database.
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
