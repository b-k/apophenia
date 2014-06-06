
/** \file apop_db.c	An easy front end to SQLite. Includes a few nice
features like a variance, skew, and kurtosis aggregator for SQL. */
/* Copyright (c) 2006--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"

/** Here are where the options are initially set. See the \ref apop_opts_type
    documentation for details.*/
apop_opts_type apop_opts	= 
          { .verbose=1,
            .output_delimiter ="\t",       .input_delimiters = "|,\t", 
            .db_name_column = "row_names", .nan_string = "NaN", 
            .db_engine = '\0',             .db_user = "\0", 
            .db_pass = "\0",               .thread_count = 1,
            .log_file = NULL,
            .rng_seed = 479901,            .version = 0.999 };

#define ERRCHECK {Apop_stopif(err, return 1, 0, "%s: %s",query, err); }
#define ERRCHECK_NR {Apop_stopif(err, return NULL, 0, "%s: %s",query, err); }
#define ERRCHECK_SET_ERROR(outdata) {Apop_stopif(err, if (!(outdata)) (outdata)=apop_data_alloc(); (outdata)->error='q'; sqlite3_free(err); return outdata, 0, "%s: %s",query, err); }

#include "apop_db_sqlite.c" // callback_t is defined here, btw.


#ifdef HAVE_MYSQL
//Let mysql have these.
#undef VERSION
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef PACKAGE_BUGREPORT
#include "apop_db_mysql.c"
#endif

//if !apop_opts.db_engine, run this to assign a value.
static void get_db_type(){
    if (getenv("APOP_DB_ENGINE") && (!strcasecmp(getenv("APOP_DB_ENGINE"), "mysql") || !strcasecmp(getenv("APOP_DB_ENGINE"), "mariadb")))
        apop_opts.db_engine = 'm';
    else
        apop_opts.db_engine = 's';
}

//This macro declares the query string and fills it from the printf part of the call.
#define Fillin(query, fmt)        \
    char *query;                  \
    va_list argp;                 \
	va_start(argp, fmt);          \
	Apop_stopif(vasprintf(&query, fmt, argp)==-1, , 0, "Trouble writing to a string."); \
	va_end(argp);                 \
	Apop_notify(2, "%s", query);

/** If you want to use a database on the hard drive instead of memory, then call this
once and only once before using any other database utilities.

If you want a disposable database which you won't use after the program ends, don't bother with this function.

The trade-offs between an on-disk database and an in-memory db are as one would expect: memory is faster, but is destroyed when the program exits. SQLite includes a command line utility (<tt>sqlite3</tt>) which let you ask queries of a database on disk, which may be useful for debugging. There are also some graphical front-ends; just ask your favorite search engine for <a href="http://www.google.com/search?&q=sqlite+gui">SQLite GUI</a>.

MySQL users: either set the environment variable APOP_DB_ENGINE=mysql or set \c apop_opts.db_engine = 'm'.

The Apophenia package assumes you are only using a single SQLite database at a time. You can use the SQL <tt>attach</tt> function to load other databases, or see <a href="http://modelingwithdata.org/arch/00000142.htm">this blog post</a> for further suggestions and sample code.

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
    if (!apop_opts.db_engine) get_db_type();
    if (!db) //check the environment.
#ifdef HAVE_MYSQL
       if(!mysql_db)  
#endif

    if (apop_opts.db_engine == 'm')
#ifdef HAVE_MYSQL
        return apop_mysql_db_open(filename);
#else
        {Apop_stopif(1, return -1, 0, "Apophenia was compiled without mysql support.");}
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

\li If <tt>apop_opts.stop_on_warn='n'</tt>, returns -1 on errors.

\li This function uses the \ref designated syntax for inputs.
\ingroup db
*/
#ifdef APOP_NO_VARIADIC
int apop_table_exists(char const *name, char remove){
#else
apop_varad_head(int, apop_table_exists){
    char const *apop_varad_var(name, NULL)
    Apop_stopif(!name, return -1, 0, "You gave me a NULL table name.");
    char apop_varad_var(remove, 'n')
    return apop_table_exists_base(name, remove);
}

 int apop_table_exists_base(char const *name, char remove){
#endif
    if (!apop_opts.db_engine) get_db_type();
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_MYSQL
        return apop_mysql_table_exists(name, remove);
#else
        Apop_stopif(1, return -1, 0, "Apophenia was compiled without mysql support.");
#endif
    char *err=NULL, *q2;
    tab_exists_t te = { .name = name };
    tab_exists_t tev = { .name = name };
	if (db==NULL) return 0;
	sqlite3_exec(db, "select name from sqlite_master where type='table'", tab_exists_callback, &te, &err); 
	sqlite3_exec(db, "select name from sqlite_master where type='view'", tab_exists_callback, &tev, &err); 
    char query[]="Selecting names from sqlite_master";//for ERRCHECK.
	ERRCHECK
	if ((remove==1|| remove=='d') && (te.isthere||tev.isthere)){
        if (te.isthere)
            Asprintf(&q2, "drop table %s;", name);
        else
            Asprintf(&q2, "drop view %s;", name);
		sqlite3_exec(db, q2, NULL, NULL, &err); 
        free(q2);
        ERRCHECK
    }
	return (te.isthere||tev.isthere);
}

/**
Closes the database on disk. If you opened the database with \c apop_db_open(NULL), then this is basically optional.

\param vacuum 
'v': vacuum---do clean-up to minimize the size of the database on disk.<br>
'q': Don't bother; just close the database. (default = 'q')

\return 0 on OK, nonzero on error.
\li This function uses the \ref designated syntax for inputs.
*/
#ifdef APOP_NO_VARIADIC
int apop_db_close(char vacuum){
#else
apop_varad_head(int, apop_db_close){
    char apop_varad_var(vacuum, 'q')
    return apop_db_close_base(vacuum);
}

 int apop_db_close_base(char vacuum){
#endif
    if (apop_opts.db_engine == 'm') //assume this is set by now...
#ifdef HAVE_MYSQL
        {apop_mysql_db_close(0);
        return 0;}
#else
        {Apop_stopif(1, return -1, 0, "Apophenia was compiled without mysql support.");}
#endif
    else {
        char *err, *query = "db close";//for errcheck.
        if (vacuum==1 || vacuum=='v') {
            sqlite3_exec(db, "VACUUM", NULL, NULL, &err);
            ERRCHECK
        }
        sqlite3_close(db);
    	//ERRCHECK
        db  = NULL;
    }
    return 0;
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

\li Blanks in the database (i.e., <tt> NULL</tt>s) and elements that match \ref apop_opts_type "apop_opts.nan_string" are filled with <tt>NAN</tt>s in the matrix.

  \{
  */

/** Send a query to the database that returns no data.

\li As with the \c apop_query_to_... functions, the query can include printf-style format specifiers, such as <tt>apop_query("create table %s(id, name, age);", tablename)</tt>.

\param fmt A <tt>printf</tt>-style SQL query.
\return 0 on success, 1 on failure.
*/
int apop_query(const char *fmt, ...){
    char *err=NULL;
    Fillin(query, fmt)
    if (!apop_opts.db_engine) get_db_type();
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_MYSQL
        {Apop_stopif(!mysql_db, return 1, 0, "No mySQL database is open.");
        return apop_mysql_query(query);}
#else
        Apop_stopif(1, return 1, 0, "Apophenia was compiled without mysql support.");
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

\return	 An \ref apop_data structure with the <tt>text</tt> element filled.

\param fmt A <tt>printf</tt>-style SQL query.

\li If <tt>apop_opts.db_name_column</tt> matches a column of the output table, then that column is used for row names, and therefore will not be included in the <tt>text</tt>.

\li <tt>query_output->text</tt> is always a 2-D array of strings, even if the query returns a single column. In that case, use <tt>returned_tab->text[i][0]</tt> (or equivalently, <tt>*returned_tab->text[i]</tt>) to refer to row <tt>i</tt>.

\li If an element in the database is \c NULL, the corresponding cell in the output table will be filled with the text given by \c apop_opts.nan_string. The default is \c "NaN", but you can set <tt>apop_opts.nan_string = "whatever you like"</tt> to change the text to whatever you like.

\li Returns \c NULL if your query is valid but returns zero rows.

\li As with the other \c apop_query_to_... functions, the query can include printf-style format specifiers, such as <tt>apop_query_to_text("select name from %s where id=%i;", tablename, id_number)</tt>.

For example, the following function will list the tables in an SQLite database (much like you
could do from the command line using <tt>sqlite3 dbname.db ".table"</tt>).

\exception out->error=='q' The database engine was unable to run the query (e.g.,  invalid SQL syntax). Again, a valid query that returns zero rows is not an error, and \c NULL is returned.
\exception out->error=='d' Database error.
\include ls_tables.c
*/
apop_data * apop_query_to_text(const char * fmt, ...){
    apop_data *out = NULL;
    Fillin(query, fmt)
    if (!apop_opts.db_engine) get_db_type();
    if (apop_opts.db_engine == 'm'){
#ifdef HAVE_MYSQL
        out = apop_mysql_query_core(query, process_result_set_chars);
#else
        Apop_stopif(1, apop_return_data_error('d'), 0, "Apophenia was compiled without mysql support.");
#endif
    } else out = apop_sqlite_query_to_text(query);
    free(query);
    return out;
}

//apop_query_to_data callback.
static int db_to_table(void *qinfo, int argc, char **argv, char **column){
    Apop_stopif(!argv, return -1, apop_errorlevel, "Got NULL data from SQLite.");
    int i, ncfound = 0;
    callback_t *qi= qinfo;
    if (qi->firstcall){
        qi->firstcall--;
        for(i=0; i<argc; i++)
            if (!strcasecmp(column[i], apop_opts.db_name_column)){
                qi->namecol = i;
                ncfound = 1;
                break;
            }
	    qi->outdata = argc-ncfound ? apop_data_alloc(1, argc-ncfound) : apop_data_alloc( );
        for(i=0; i<argc; i++)
            if (qi->namecol != i)
                apop_name_add(qi->outdata->names, column[i], 'c');
    } else 
        if (qi->outdata->matrix)
            apop_matrix_realloc(qi->outdata->matrix, qi->currentrow+1, qi->outdata->matrix->size2);
    ncfound =0;
    for (int jj=0;jj<argc;jj++)
        if (jj != qi->namecol){
            double valor = 
                !argv[jj] || !strcmp(argv[jj], "NULL")|| 
                (apop_opts.nan_string && !strcasecmp(apop_opts.nan_string, argv[jj]))
                 ? GSL_NAN : atof(argv[jj]);
            gsl_matrix_set(qi->outdata->matrix,qi->currentrow,jj-ncfound, valor);
        } else {
            apop_name_add(qi->outdata->names, argv[jj], 'r');
            ncfound = 1;
        }
    (qi->currentrow)++;
	return 0;
}

/** Queries the database, and dumps the result into an \ref apop_data set.

\li If \ref apop_opts_type "apop_opts.db_name_column" is set (it defaults to being "row_names"), and the name of a column matches the name, then the row names are read from that column.

\li As with the other \c apop_query_to_... functions, the query can include printf-style format specifiers, such as <tt>apop_query_to_data("select age from %s where id=%i;", tablename, id_number)</tt>.

\return If no rows are returned, \c NULL; else an \ref apop_data set with the data in place. Most data will be in the \c matrix element of the output. Column names are appropriately placed. If <tt>apop_opts.db_name_column</tt> matches one of the fields in your query's output, then that column will be used for row names (and therefore will not appear in the \c matrix).

\param fmt A <tt>printf</tt>-style SQL query.
\exception out->error=='q' Query error. A valid query that returns no rows is not an error; in that case, you get \c NULL.
*/ 
apop_data * apop_query_to_data(const char * fmt, ...){
    Fillin(query, fmt)
    if (!apop_opts.db_engine) get_db_type();
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_MYSQL
        return apop_mysql_query_core(query, process_result_set_data);
#else
        Apop_stopif(1, apop_return_data_error('d'), 0, "Apophenia was compiled without mysql support.");
#endif

    //else
    char *err=NULL;
    callback_t qinfo = {.firstcall = 1, .namecol=-1};
	if (db==NULL) apop_db_open(NULL);
    sqlite3_exec(db, query,db_to_table,&qinfo, &err); 
    free (query);
    ERRCHECK_SET_ERROR(qinfo.outdata)
	return qinfo.outdata;
}


    /** \cond doxy_ignore */
//These used to do more, but I'll leave them as a macro anyway in case of future expansion.
#define Store_settings  \
    int v = apop_opts.verbose; apop_opts.verbose=0;/*hack to prevent double-printing.*/ \

#define Restore_settings  \
    apop_opts.verbose=v;
    /** \endcond */

/** Queries the database, and dumps the result into a matrix.

  Uses \ref apop_query_to_data and returns just the matrix part; see that function for notes.

\li If \c apop_opts.db_name_column is set, then I'll ignore that column. It gets put into the names of the \ref apop_data set, and then thrown away when I return only the \c gsl_matrix part of that set. 

\li As with the other \c apop_query_to_... functions, the query can include printf-style format specifiers, such as <tt>apop_query_to_matrix("select age from %s where id=%i;", tablename, id_number)</tt>.
 
\deprecated Use \ref apop_query_to_data
\param fmt A <tt>printf</tt>-style SQL query.
\return A \c gsl_matrix.
\exception out->error=='q' Query error. A valid query that returns no rows is not an error; in that case, you get \c NULL.
 */
gsl_matrix * apop_query_to_matrix(const char * fmt, ...){
    Fillin(query, fmt)
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

/** Queries the database, and dumps the first column of the result into a \c gsl_vector.

\li Uses \ref apop_query_to_data internally, then throws away all but the first column of the matrix.

\li If \c apop_opts.db_name_column is set, then I'll ignore that column. It gets put into the names of the \ref apop_data set, and then thrown away when I look at only the \c gsl_matrix part of that set. 

\li If the query returns zero rows of data or no columns, the function returns \c NULL.  

\li As with the other \c apop_query_to_... functions, the query can include printf-style format specifiers, such as <tt>apop_query_to_vector("select age from %s where id=%i;", tablename, id_number)</tt>.

\return	 A <tt>gsl_vector</tt> holding the first column of the returned matrix. Thus, if your query returns multiple lines, you will get no warning, and the function will return the first in the list.

\param fmt A <tt>printf</tt>-style SQL query.
\exception out->error=='q' Query error. A valid query that returns no rows is not an error; in that case, you get \c NULL. */
gsl_vector * apop_query_to_vector(const char * fmt, ...){
    Fillin(query, fmt)
    if (!apop_opts.db_engine) get_db_type();
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_MYSQL
        return apop_mysql_query_core(query, process_result_set_vector);
#else
        Apop_stopif(1, return NULL, 0, "Apophenia was compiled without mysql support.");
#endif
    apop_data *d=NULL;
    gsl_vector *out;
	if (db==NULL) apop_db_open(NULL);
    Store_settings
	d	= apop_query_to_data("%s", query);
    Restore_settings
    Apop_stopif(!d, return NULL, 2, "Query [%s] turned up a blank table. Returning NULL.", query);
    //else:
    out = gsl_vector_alloc(d->matrix->size1);
	gsl_matrix_get_col(out, d->matrix, 0);
	apop_data_free(d);
    free(query);
	return out;
}

/** Queries the database, and dumps the result into a single double-precision floating point number.

\li This calls \ref apop_query_to_data and returns the (0,0)th element of the returned matrix. Thus, if your query returns multiple lines, you will get no warning, and the function will return the first in the list (which is not always well-defined; maybe use an <tt>order by</tt> clause in your query if you expect multiple lines).

\li If \c apop_opts.db_name_column is set, then I'll ignore that column. It gets put into the names of the \ref apop_data set, and then thrown away when I look at only the \c gsl_matrix element of that set. 

\li If the query returns no rows at all, the function returns <tt>NAN</tt>.

\li If the query produces a blank table, returns \c NAN, and if <tt>apop_opts.verbose>=2</tt>, prints an error.

\li As with the other \c apop_query_to_... functions, the query can include printf-style format specifiers, such as <tt>apop_query_to_float("select age from %s where id=%i;", tablename, id_number)</tt>.

\li If the query produces an error, returns \c NAN, and if <tt>apop_opts.verbose>=0</tt>, prints an error. If you need to distinguish between blank tables, NaNs in the data, and query errors, use \ref apop_query_to_data.

\param fmt A <tt>printf</tt>-style SQL query.
\return		A \c double, actually.*/
double apop_query_to_float(const char * fmt, ...){
    double out;
    Fillin(query, fmt)
    if (!apop_opts.db_engine) get_db_type();
    if (apop_opts.db_engine == 'm'){
#ifdef HAVE_MYSQL
        out = apop_mysql_query_to_float(query);
#else
        Apop_stopif(1, return NAN, 0, "Apophenia was compiled without mysql support.");
#endif
    } else {
        apop_data *d=NULL;
        if (db==NULL) apop_db_open(NULL);
        Store_settings
        d = apop_query_to_data("%s", query);
        Restore_settings
        Apop_stopif(!d, return GSL_NAN, 2, "Query [%s] turned up a blank table. Returning NaN.", query);
        Apop_stopif(d->error, return GSL_NAN, 0, "Query [%s] failed. Returning NaN.", query);
        out	= apop_data_get(d);
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

\li \ref apop_opts_type "apop_opts.db_name_column" is ignored.  Use the \c 'n' character to indicate the output column with row names.

\li As with the other \c apop_query_to_... functions, the query can include printf-style format specifiers, such as <tt>apop_query_to_mixed_data("tv", "select name, age from %s where id=%i", tablename, id_number)</tt>.

\param typelist A string consisting of the letters \c nvmtw. For example, if your query columns should go into a text column, the vector, the weights, and two matrix columns, this would be "tvwmm".
\param fmt A <tt>printf</tt>-style SQL query.
\exception out->error=='d' Dimension error. Your count of matrix parts didn't match what the query returned.
\exception out->error=='q' Query error. A valid query that returns no rows is not an error; in that case, you get \c NULL.
*/
apop_data * apop_query_to_mixed_data(const char *typelist, const char * fmt, ...){
    Fillin(query, fmt)
    if (!apop_opts.db_engine) get_db_type();
    if (apop_opts.db_engine == 'm')
#ifdef HAVE_MYSQL
        {apop_data* out = apop_mysql_mixed_query(typelist, query);
        free(query);
        return out;}
#else
        {Apop_notify(0, "Apophenia was compiled without mysql support.");
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
    Apop_stopif(vasprintf(q, format, ap)==-1, , 0, "Trouble writing to a string.");
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

static int run_prepared_statements(apop_data const *set, sqlite3_stmt *p_stmt){
#if SQLITE_VERSION_NUMBER < 3003009
     Apop_stopif(1, return -1, 0, "Attempting to use prepared statements, but using a version of SQLite that doesn't support them.");
#else
    Get_vmsizes(set) //firstcol, msize1, maxsize
    for (size_t row=0; row < maxsize; row++){
        size_t field =1;
        if (set->names->rowct>row){
            if (!strlen(set->names->row[row])) field++; //leave NULL and cleared
            Apop_stopif(sqlite3_bind_text(p_stmt, field++, set->names->row[row], -1, SQLITE_TRANSIENT),
                    return -1, apop_errorlevel, 
                    "Something wrong with the row name for line %zu, [%s].\n" , row, set->names->row[row]);
        }
        if (set->vector && set->vector->size > row)
                Apop_stopif(sqlite3_bind_double(p_stmt, field++, apop_data_get(set, row, -1)),
                    return -1, apop_errorlevel, 
                    "Something wrong with the vector element on line %zu, [%g].\n" ,row,  apop_data_get(set, row, -1));
        if (msize1 > row)
            for (size_t col=0; col < msize2; col++)
                Apop_stopif(sqlite3_bind_double(p_stmt, field++, apop_data_get(set, row, col)),
                    return -1, apop_errorlevel, 
                    "Something wrong with the matrix element %zu on line %zu, [%g].\n" ,col, row,  apop_data_get(set, row, col));
        if (*set->textsize > row)
            for (size_t col=0; col < set->textsize[1]; col++){
                if (!strlen(set->text[row][col])) field++; //leave NULL and cleared
                Apop_stopif(sqlite3_bind_text(p_stmt, field++, set->text[row][col], -1, SQLITE_TRANSIENT),
                    return -1, apop_errorlevel, 
                    "Something wrong with the row name for line %zu, [%s].\n" , row, set->text[row][col]);
            }
        if (set->weights && set->weights->size > row)
                Apop_stopif(sqlite3_bind_double(p_stmt, field++, gsl_vector_get(set->weights, row)),
                    return -1, apop_errorlevel, 
                    "Something wrong with the weight element on line %zu, [%g].\n" ,row,  gsl_vector_get(set->weights, row));
        int err = sqlite3_step(p_stmt);
        Apop_stopif(err!=0 && err != 101 //0=ok, 101=done
                    , , 0, "prepared sqlite insert query gave error code %i.\n", err);
        Apop_stopif(sqlite3_reset(p_stmt), return -1, apop_errorlevel, "SQLite error.");
        Apop_stopif(sqlite3_clear_bindings(p_stmt), return -1, apop_errorlevel, "SQLite error."); //needed for NULLs
    }
    Apop_stopif(sqlite3_finalize(p_stmt)!=SQLITE_OK, return -1, apop_errorlevel, "SQLite error.");
    return 0;
#endif
}

/** Dump an \ref apop_data set into the database.

This function is basically preempted by \ref apop_data_print. Use that one; this may soon no longer be available.

Column names are inserted if there are any. If there are, all dots are converted to underscores.  Otherwise, the columns will be named \c c1, \c c2, \c c3, &c.

\li If \ref apop_opts_type "apop_opts.db_name_column" is not blank (the default is "row_name"), then a so-named column is created, and the row names are placed there.

\li If there are weights, they will be the last column of the table, and the column will be named "weights".

\li If the table exists; append to. If the table does not exist, create. So perhaps call \ref apop_table_exists <tt>("tabname", 'd')</tt> to ensure that the table is removed ahead of time.

\li You can also call this via \ref apop_data_print <tt>(data, "tabname", .output_type='d', .output_append='w')</tt> to overwrite a new table or with <tt>.output_append='a'</tt> to append.

\li If your data set has zero data (i.e., is just a list of column names or is entirely blank), I return -1 without creating anything in the database.

\li Especially if you are using a pre-2007 version of SQLite, there may be a speed gain to wrapping the call to this function in a begin/commit pair:

\code
apop_query("begin;");
apop_data_print(dataset, .output_name="dbtab", .output_type='d');
apop_query("commit;");
\endcode


\param set 	         The name of the matrix
\param tabname	     The name of the db table to be created
\param output_append See \ref apop_prep_output.
\return 0=OK, -1=error
\ingroup apop_data
\ingroup conversions
*/
int apop_data_to_db(const apop_data *set, const char *tabname, const char output_append){
    Apop_stopif(!set, return -1, 1, "you sent me a NULL data set. Database table %s will not be created.", tabname);
    int	i,j; 
    char *q;
    char comma = ' ';
    int use_row = strlen(apop_opts.db_name_column) 
                && ((set->matrix && set->names->rowct == set->matrix->size1)
                    || (set->vector && set->names->rowct == set->vector->size));

    if (!apop_opts.db_engine) get_db_type();
    if (apop_table_exists(tabname))
        Asprintf(&q, " ");
    else if (apop_opts.db_engine == 'm')
#ifdef HAVE_MYSQL
        if (((output_append =='a' || output_append =='A') && apop_table_exists(tabname)))
            Asprintf(&q, " ");
        else {
            Asprintf(&q, "create table %s (", tabname);
            if (use_row) {
                qxprintf(&q, "%s\n %s varchar(1000)", q, apop_opts.db_name_column);
                comma = ',';
            }
            if (set->vector){
                if(!set->names->vector) 
                    qxprintf(&q, "%s%c\n vector double ", q, comma);
                else
                    qxprintf(&q, "%s%c\n %s double ", q,comma, set->names->vector);
                comma = ',';
            }
            if (set->matrix)
                for(i=0;i< set->matrix->size2; i++){
                    if(set->names->colct <= i) 
                        qxprintf(&q, "%s%c\n c%i double ", q, comma,i);
                     else
                        qxprintf(&q, "%s%c\n %s  double ", q, comma, set->names->col[i]);
                    comma = ',';
                }
            for(i=0;i< set->textsize[1]; i++){
                if (set->names->textct <= i)
                    qxprintf(&q, "%s%c\n tc%i varchar(1000) ", q, comma,i);
                else
                    qxprintf(&q, "%s%c\n %s  varchar(1000) ", q, comma, set->names->text[i]);
                comma = ',';
            }
            apop_query("%s); ", q);
            sprintf(q, " ");
        }
#else 
        Apop_stopif(1, return -1, apop_errorlevel, "Apophenia was compiled without mysql support.");
#endif
    else {
        if (db==NULL) apop_db_open(NULL);
        if (((output_append =='a' || output_append =='A') && apop_table_exists(tabname)) )
            Asprintf(&q, " ");
        else {
            Asprintf(&q, "create table %s (", tabname);
            if (use_row) {
                qxprintf(&q, "%s\n %s", q, apop_opts.db_name_column);
                comma = ',';
            }
            if (set->vector){
                if (!set->names->vector) qxprintf(&q, "%s%c\n vector numeric", q, comma);
                else qxprintf(&q, "%s%c\n \"%s\"", q, comma, set->names->vector);
                comma = ',';
            }
            if (set->matrix)
                for(i=0;i< set->matrix->size2; i++){
                    if(set->names->colct <= i) 	
                        qxprintf(&q, "%s%c\n c%i numeric", q, comma,i);
                    else			
                        qxprintf(&q, "%s%c\n \"%s\" numeric", q, comma, set->names->col[i]);
                    comma = ',';
                }
            for(i=0; i< set->textsize[1]; i++){
                if(set->names->textct <= i) qxprintf(&q, "%s%c\n tc%i ", q, comma, i);
                else qxprintf(&q, "%s%c\n %s ", q, comma, set->names->text[i]);
                comma = ',';
            }
            if (set->weights) qxprintf(&q, "%s%c\n \"weights\" numeric", q, comma);
            qxprintf(&q,"%s);",q);
            apop_query("%s", q);
            qxprintf(&q," ");
        }
    }

    Get_vmsizes(set) //firstcol, msize2, maxsize
    int col_ct = !!set->names->rowct + set->textsize[1] + msize2 - firstcol + !!set->weights;
    Apop_stopif(!col_ct, return -1, 0, "Input data set has zero columns of data (no rownames, text, matrix, vector, or weights). I can't create a table like that, sorry.");
    if(apop_use_sqlite_prepared_statements(col_ct)){
        sqlite3_stmt *statement;
        Apop_stopif(
            apop_prepare_prepared_statements(tabname, col_ct, &statement), 
            return -1, 0, "Trouble preparing prepared statements.");
        Apop_stopif(
            run_prepared_statements(set, statement), 
            return -1, 0, "error in insertions.");
    } else {
        for(i=0; i< maxsize; i++){
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
            apop_query("%s", q); 
            q[0]='\0';
        }
    }
	free(q);
    return 0;
}

// Some stats wrappers

/** Do a \f$t\f$-test entirely inside the database.
  Returns only the two-tailed p-value.
\ingroup ttest
*/
double apop_db_t_test(char * tab1, char *col1, char *tab2, char *col2){
    gsl_matrix *result1, *result2;
	result1	= apop_query_to_matrix("select avg(%s), var(%s), count(*) from %s", col1, col1, tab1);
	result2	= apop_query_to_matrix("select avg(%s), var(%s), count(*) from %s", col2, col2, tab2);
    double a_avg = gsl_matrix_get(result1, 0, 0),
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
    gsl_matrix	*result=
	        apop_query_to_matrix("select avg(%s - %s), var(%s - %s), count(*) from %s tab1",
                                            col1,col2,   col1, col2,             tab1);
    double avg	 = gsl_matrix_get(result, 0, 0),
		   var	 = gsl_matrix_get(result, 0, 1),
		   count = gsl_matrix_get(result, 0, 2),
		   stat	 = avg/ sqrt(var/(count-1));
	return 2*GSL_MIN(gsl_cdf_tdist_P(stat, count-1),gsl_cdf_tdist_Q(stat, count-1));
}
