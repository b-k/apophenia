/** \file apop_db_mysql.c
This file is included directly into \ref apop_db.c. It is read only if APOP_USE_MYSQL is defined.*/

/* Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include <my_global.h>
#include <my_sys.h>
#include <mysql.h>
#include <math.h>

static MYSQL *mysql_db; 

#define Areweconected(retval) Apop_stopif(!mysql_db, return retval, 0,  \
        "No connection to a mySQL/mariadb database. apop_db_open() failure?");

static char *opt_host_name = NULL;      /* server host (default=localhost) */
static unsigned int opt_port_num = 0;   /* port number (use built-in value) */
static char *opt_socket_name = NULL;    /* socket name (use built-in value) */
static unsigned int opt_flags = 0;      /* connection flags (none) */

#define Apop_mstopif(cond, returnop, str) \
    Apop_stopif(cond, returnop, 0,         \
         str "\n mySQL/mariadb error %u: %s\n", mysql_errno (mysql_db), mysql_error (mysql_db));

static int apop_mysql_db_open(char const *in){
    Apop_stopif(!in, return 2, 0, "MySQL needs a non-NULL db name.");
    mysql_db = mysql_init (NULL);
    Apop_stopif(!mysql_db, return 1, 0, "mysql_init() failed (probably out of memory)");
    Apop_mstopif (!mysql_real_connect (mysql_db, opt_host_name, apop_opts.db_user, apop_opts.db_pass,
                        in, opt_port_num, opt_socket_name, CLIENT_MULTI_STATEMENTS+opt_flags),
                mysql_close (mysql_db); return 1, 
                "mysql_real_connect() failed");
    return 0;
}

static void apop_mysql_db_close(int ignoreme){
    if (mysql_db) mysql_close (mysql_db);
}

/*
    //Cut & pasted & cleaned from the mysql manual.
static void process_results(void){
    else                // mysql_store_result() returned nothing; should it have?
        Apop_stopif(mysql_field_count(mysql_db) == 0,  , 0, "apop_query error");
    //else query wasn't a select & just didn't return data.
}
        */

static double apop_mysql_query(char *query){
    Apop_mstopif(mysql_query(mysql_db, query), return 1, "apop_mysql_query failed");
    MYSQL_RES *result = mysql_store_result(mysql_db);
    if (result) mysql_free_result(result);
    return 0;
}

static double apop_mysql_table_exists(char const *table, int delme){
    Areweconected(GSL_NAN);
    MYSQL_RES *res_set = mysql_list_tables(mysql_db, table);
    Apop_mstopif(!mysql_list_tables(mysql_db, table), return GSL_NAN,
          "show tables query failed.");
    int is_found = mysql_num_rows(res_set);
    mysql_free_result(res_set);
    if (!is_found) return 0;

    if (delme =='d' || delme=='D'){
       char *a_query;
       Asprintf(&a_query, "drop table %s", table);
       Apop_mstopif(mysql_query (mysql_db, a_query), GSL_NAN, 
            "table exists, but table dropping failed");
    }
    return 1;
}

#define check_and_clean(do_if_failure) \
    Apop_mstopif( mysql_errno (conn),   \
         if (out) do_if_failure; return NULL, \
         "mysql_fetch_row() failed"); \
    return out; \

static int get_name_row(unsigned int *num_fields, MYSQL_FIELD *fields){
    for(size_t i = 0; i < *num_fields; i++)
        if (!strcasecmp(fields[i].name, apop_opts.db_name_column)){
            (*num_fields)--;
            return i;
        }
    return -1;
}

static void * process_result_set_data (MYSQL *conn, MYSQL_RES *res_set) {
    MYSQL_ROW row;
    unsigned int num_fields = mysql_num_fields(res_set);
    unsigned int num_rows = mysql_num_rows (res_set);
    if (!num_fields || !num_rows) return NULL;

    MYSQL_FIELD *fields = mysql_fetch_fields(res_set);
    int name_row = get_name_row(&num_fields, fields);

    apop_data *out = apop_data_alloc(0, num_rows, num_fields);

    for(size_t i = 0; i < num_fields+ (name_row>=0); i++)
        if (i!=name_row) apop_name_add(out->names, fields[i].name, 'c');

    for (int i=0; (row = mysql_fetch_row (res_set)); i++) {
        int passed_name = 0;
        for (size_t j = 0; j < mysql_num_fields (res_set); j++){
            if (j==name_row){
                apop_name_add(out->names, row[j], 'r');
                passed_name = 1;
                continue;
            }
            if (!row[j]) apop_data_set(out, i , j-passed_name, NAN);
            else {
                char *end = NULL;
                double num = strtod(row[j], &end);
                apop_data_set(out, i , j-passed_name, *end ? NAN : num);
            }
       }
    }
    check_and_clean(apop_data_free(out))
}

static void * process_result_set_vector (MYSQL *conn, MYSQL_RES *res_set) {
    MYSQL_ROW row;
    unsigned int num_fields = mysql_num_fields(res_set);
    unsigned int num_rows = mysql_num_rows (res_set);
    if (num_fields == 0 || num_rows == 0) return NULL;
    gsl_vector *out = gsl_vector_alloc(num_rows);
    for (int j=0; (row = mysql_fetch_row (res_set)); j++){
        double valor = (!row[0] || !strcmp(row[0], "NULL"))
                           ? GSL_NAN : atof(row[0]);
        gsl_vector_set(out, j, valor);
    }
    check_and_clean(gsl_vector_free(out))
}

static void * process_result_set_chars (MYSQL *conn, MYSQL_RES *res_set) {
    MYSQL_ROW row;
    unsigned int total_cols = mysql_num_fields(res_set);
    unsigned int total_rows = mysql_num_rows(res_set);

    MYSQL_FIELD *fields = mysql_fetch_fields(res_set);
    int name_row = get_name_row(&total_cols, fields);
    apop_data *out = apop_text_alloc(NULL, total_rows, total_cols);

    for (size_t i = 0; i < total_cols + (name_row>=0); i++)
        if (i!=name_row) apop_name_add(out->names, fields[i].name, 't');

    for (int i=0; (row = mysql_fetch_row (res_set)); i++){
        int passed_name = 0;
		for (size_t jj=0; jj<total_cols; jj++){
            if (jj==name_row){
                apop_name_add(out->names, row[jj], 'r');
                passed_name = 1;
                continue;
            }
            apop_text_add(out, i, jj-passed_name, "%s", (row[jj]==NULL)?  apop_opts.nan_string : row[jj]);
		}
    }
    check_and_clean(;)
}

static void * apop_mysql_query_core(char *query, void *(*callback)(MYSQL*, MYSQL_RES*)){
    Areweconected(NULL);
    apop_data *output = NULL;
    Apop_mstopif(mysql_query (mysql_db, query), return NULL, "mysql_query() failed");
    MYSQL_RES *res_set = mysql_store_result (mysql_db);
    Apop_mstopif(!res_set, 
        if (callback == process_result_set_data || callback==process_result_set_data) apop_return_data_error('q') 
            else return NULL, 
            "mysql_store_result() failed");
    if (!res_set->row_count) goto done; //just a blank table.
    output = callback(mysql_db, res_set);

    done:
    mysql_free_result (res_set);
    return output;
}

static double apop_mysql_query_to_float(char *query){
    Areweconected(GSL_NAN);
    Apop_mstopif(mysql_query (mysql_db, query) != 0, return GSL_NAN,
          "mysql_query() failed");
    MYSQL_RES *res_set = mysql_store_result (mysql_db);
    Apop_mstopif(!res_set, return GSL_NAN, "mysql_store_result() failed");
    if (mysql_num_rows(res_set)==0) return GSL_NAN;
    MYSQL_ROW row = mysql_fetch_row (res_set);
    Apop_mstopif(mysql_errno (mysql_db),
        mysql_free_result (res_set); return GSL_NAN,
        "mysql_fetch_row() failed");
    double out = atof(row[0]);
    mysql_free_result (res_set);
    return out;
}

apop_data* apop_mysql_mixed_query(char const *intypes, char const *query){
    Areweconected(NULL);
    apop_data *out = NULL;
    Apop_mstopif(mysql_query (mysql_db, query), return NULL, "mysql_query() failed");
    MYSQL_RES *res_set = mysql_store_result(mysql_db);
    MYSQL_ROW row;
    Apop_mstopif(!res_set, return NULL, "mysql_store_result() failed");
    if (!res_set->row_count) goto done; //just a blank table.

    unsigned int total_cols = mysql_num_fields(res_set);
    unsigned int total_rows = mysql_num_rows(res_set);
    if (!total_cols || !total_rows) goto done;

    apop_qt info = { };
    count_types(&info, intypes); //in apop_db_sqlite.c
    //intypes[5] === names, vectors, mcols, textcols, weights.

    out = apop_data_alloc(info.intypes[1] ? total_rows : 0, 
                           info.intypes[2] ? total_rows : 0,  
                           info.intypes[2]);

    int requested = info.intypes[0]+info.intypes[1]+info.intypes[2]+info.intypes[3]+info.intypes[4];
    int excess = requested - total_cols;
    Apop_stopif(excess > 0, out->error='d' /*and continue.*/, 1, 
      "you asked for %i columns in your list of types(%s), but your query produced %u columns. "
      "The remainder will be placed in the text section. Output data set's ->error element set to 'd'." , requested, intypes, total_cols);
    Apop_stopif(excess < 0, out->error='d' /*and continue.*/, 1, 
      "you asked for %i columns in your list of types(%s), but your query produced %u columns. "
      "Ignoring the last %i type(s) in your list. Output data set's ->error element set to 'd'." , requested, intypes, total_cols, -excess);

    if (info.intypes[3]||excess>0) apop_text_alloc(out, total_rows, info.intypes[3] + ((excess > 0) ? excess : 0));
    if (info.intypes[4]) out->weights = gsl_vector_alloc(total_rows);

    MYSQL_FIELD *fields = mysql_fetch_fields(res_set);
    for (size_t i=0; i<total_cols; i++){
        char c = (i < requested) ? intypes[i] : 't';
        if (c == 't'|| c=='T')
            apop_name_add(out->names, fields[i].name, 't');
        else if (c == 'v'|| c=='V')
            apop_name_add(out->names, fields[i].name, 'v');
        else if (c == 'm'|| c=='M')
            apop_name_add(out->names, fields[i].name, 'c');
    }

    for (int i=0; (row = mysql_fetch_row (res_set)); i++) {
        int thism=0, thist=0;
		for (size_t j=0; j<total_cols; j++){
            char c = (j < requested) ? intypes[j] : 't';
            if (c == 'n' || c =='N')
                apop_name_add(out->names, row[j], 'r');
            else if (c == 't'|| c=='T')
                apop_text_add(out, i, thist++, "%s", (row[j]==NULL)?  apop_opts.nan_string : row[j]);
            else if (c == 'v'|| c=='V'){
                double valor = (!row[j] || !strcmp(row[j], "NULL")) ? NAN : atof(row[j]);
                gsl_vector_set(out->vector, i, valor);
            } else if (c == 'w'|| c=='W'){
                double valor = (!row[j] || !strcmp(row[j], "NULL")) ? NAN : atof(row[j]);
                gsl_vector_set(out->weights, i, valor);
            } else if (c == 'm'|| c=='M')
                gsl_matrix_set(out->matrix, i , thism++, row[j] ? atof(row[j]): GSL_NAN);
		}
    }

    done:
    mysql_free_result (res_set);
    return out;
}
