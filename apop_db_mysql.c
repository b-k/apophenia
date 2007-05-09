/** \file apop_mysql.c


This file is included directly into \ref apop_db.c. It is read only if APOP_USE_MYSQL is defined.
*/
#include <my_global.h>
#include <my_sys.h>
#include <mysql.h>


static MYSQL *mysql_db; 

static char *opt_host_name = NULL;      /* server host (default=localhost) */
static char *opt_user_name = NULL;      /* username (default=login name) */
static char *opt_password = NULL;       /* password (default=none) */
static unsigned int opt_port_num = 0;   /* port number (use built-in value) */
static char *opt_socket_name = NULL;    /* socket name (use built-in value) */
static unsigned int opt_flags = 0;      /* connection flags (none) */


/** This function and the kernel of a few other routines cut and pasted from 
_MySQL (Third Edition)_, Paul DuBois, Sams Developer's Library, March, 2005
*/
static void print_error (MYSQL *conn, char *message) {
    fprintf (stderr, "%s\n", message);
    if (conn != NULL)
    {
#if MYSQL_VERSION_ID >= 40101
        fprintf (stderr, "Error %u (%s): %s\n",
            mysql_errno (conn), mysql_sqlstate(conn), mysql_error (conn));
#else
        fprintf (stderr, "Error %u: %s\n",
            mysql_errno (conn), mysql_error (conn));
#endif
    }
}

    //Cut & pasted from the mysql manual.
static void process_results(void){
  MYSQL_RES *result;
  unsigned int num_fields;
  unsigned int num_rows;
    result = mysql_store_result(mysql_db);
    if (result){  // there are rows
        num_fields = mysql_num_fields(result);
        mysql_free_result(result);
        // retrieve rows, then call mysql_free_result(result)
    } else {  // mysql_store_result() returned nothing; should it have?
        if(mysql_field_count(mysql_db) == 0) {
            // query does not return data
            // (it was not a SELECT)
            num_rows = mysql_affected_rows(mysql_db);
        } else  // mysql_store_result() should have returned data
            fprintf(stderr, "apop_query error: %s\n", mysql_error(mysql_db));
    }
}

static int apop_mysql_db_open(char *in){
    if (!in){
        fprintf(stderr, "MySQL needs a non-NULL db name.");
        return 1;
    } else{
        mysql_db = mysql_init (NULL);
        if (!mysql_db) {
            fprintf (stderr, "mysql_init() failed (probably out of memory)\n");
            return 1;
        }
        /* connect to server */
        if (!mysql_real_connect (mysql_db, opt_host_name, opt_user_name, opt_password,
                in, opt_port_num, opt_socket_name, CLIENT_MULTI_STATEMENTS+opt_flags)) {
                    fprintf (stderr, "mysql_real_connect() to %s failed\n", in);
                    mysql_close (mysql_db);
                    return 1;
                }
    }
    return 0;
}


static double apop_mysql_query(char *query){
    if (mysql_query(mysql_db,query)) {
        print_error (mysql_db, "apop_mysql_query failed");
        return 1;
    } 
    process_results();
    return 0;
}

static void apop_mysql_db_close(int ignoreme){
        mysql_close (mysql_db);
}


static double apop_mysql_table_exists(char *table, int delme){
  MYSQL_RES         *res_set;
  MYSQL_ROW         row;
  int               is_found    =0;
    if (mysql_query (mysql_db, "show tables")){
         print_error (mysql_db, "show tables query failed.");
         return GSL_NAN;
    }
    res_set = mysql_store_result (mysql_db);     // generate result set
    if (res_set == NULL){
        print_error (mysql_db, "mysql_store_result() failed");
        return GSL_NAN;
    }
    while ((row = mysql_fetch_row (res_set)) ) {
        if(!strcmp(row[0], table)){
            is_found   ++;
            break;
        }
    }
    mysql_free_result(res_set);
    if (!is_found)
       return 0;
    if (delme){
       int len         = 100+strlen(table);
       char *a_query   = malloc(len);
       snprintf(a_query, len, "drop table %s", table);
       if (mysql_query (mysql_db, a_query)) 
           print_error (mysql_db, "table dropping failed");
       } else process_results();
    return 1;
}



static apop_data * process_result_set_data (MYSQL *conn, MYSQL_RES *res_set) {
  MYSQL_ROW        row;
  unsigned int     i, j=0;
  unsigned int num_fields = mysql_num_fields(res_set);
  apop_data *out   =apop_data_alloc(0, mysql_num_rows (res_set), num_fields);
     while ((row = mysql_fetch_row (res_set)) ) {
             for (i = 0; i < mysql_num_fields (res_set); i++) {
                 apop_data_set(out, j , i, atof(row[i]));
             }
             j++;
        }
     MYSQL_FIELD *fields = mysql_fetch_fields(res_set);
     for(i = 0; i < num_fields; i++)
         apop_name_add(out->names, fields[i].name, 'c');
        if (mysql_errno (conn)){
             print_error (conn, "mysql_fetch_row() failed");
             return NULL;
        } else
            return out;
}

static apop_data* apop_mysql_query_to_data(char *query){
  MYSQL_RES *res_set;
  apop_data *out;
    if (mysql_query (mysql_db, query) != 0){
         print_error (mysql_db, "mysql_query() failed");
         return NULL;
    }
     res_set = mysql_store_result (mysql_db);     // generate result set
     if (res_set == NULL){
          print_error (mysql_db, "mysql_store_result() failed");
        return NULL;
     }
      // process result set, and then deallocate it 
      out = process_result_set_data (mysql_db, res_set);
      mysql_free_result (res_set);
      return out;
}

double apop_mysql_query_to_float(char *query){
  MYSQL_RES *res_set;
  double out;
  MYSQL_ROW        row;
    if (mysql_query (mysql_db, query) != 0){
         print_error (mysql_db, "mysql_query() failed");
         return GSL_NAN;
    } 
    res_set = mysql_store_result (mysql_db);     // generate result set
    if (res_set == NULL){
       print_error (mysql_db, "mysql_store_result() failed");
       return GSL_NAN;
    } 
    row = mysql_fetch_row (res_set);
    out    = atof(row[0]);
    if (mysql_errno (mysql_db)){
         print_error (mysql_db, "mysql_fetch_row() failed");
         return GSL_NAN;
    } 
    mysql_free_result (res_set);
    return out;
}


gsl_vector * process_result_set_vector (MYSQL *conn, MYSQL_RES *res_set) {
  MYSQL_ROW        row;
  unsigned int     j=0;
  gsl_vector *out   =gsl_vector_alloc( mysql_num_rows (res_set));
     while ((row = mysql_fetch_row (res_set)) ) {
         gsl_vector_set(out, j,  atof(row[0]));
         j++;
    }
    if (mysql_errno (conn)){
         print_error (conn, "mysql_fetch_row() failed");
         return NULL;
    } 
    return out;
}

gsl_vector* apop_mysql_query_to_vector(char *query){
  MYSQL_RES *res_set;
  gsl_vector *out;
    if (mysql_query (mysql_db, query) != 0){
         print_error (mysql_db, "mysql_query() failed");
         return NULL;
    }
    res_set = mysql_store_result (mysql_db);     // generate result set
    if (res_set == NULL){
         print_error (mysql_db, "mysql_store_result() failed");
       return NULL;
    }
    // process result set, and then deallocate it 
    out = process_result_set_vector (mysql_db, res_set);
    mysql_free_result (res_set);
    return out;
}


gsl_matrix * process_result_set_matrix (MYSQL *conn, MYSQL_RES *res_set) {
  MYSQL_ROW        row;
  unsigned int     i, j=0;
  unsigned int num_fields = mysql_num_fields(res_set);
  gsl_matrix *out   =gsl_matrix_alloc( mysql_num_rows (res_set), num_fields);
     while ((row = mysql_fetch_row (res_set)) ) {
         for (i = 0; i < mysql_num_fields (res_set); i++) {
             gsl_matrix_set(out, j , i, atof(row[i]));
         }
         j++;
    }
    if (mysql_errno (conn)){
         print_error (conn, "mysql_fetch_row() failed");
         return NULL;
    } else
        return out;
}

gsl_matrix* apop_mysql_query_to_matrix(char *query){
  MYSQL_RES *res_set;
  gsl_matrix *out;
    if (mysql_query (mysql_db, query) != 0){
         print_error (mysql_db, "mysql_query() failed");
         return NULL;
    }else {
         res_set = mysql_store_result (mysql_db);     // generate result set
         if (res_set == NULL){
              print_error (mysql_db, "mysql_store_result() failed");
            return NULL;
         }else {
              // process result set, and then deallocate it 
              out = process_result_set_matrix (mysql_db, res_set);
              mysql_free_result (res_set);
              return out;
         }
    }
}

apop_data * process_result_set_chars (MYSQL *conn, MYSQL_RES *res_set) {
  MYSQL_ROW        row;
  unsigned int     jj, currentrow = 0;
    total_cols  = mysql_num_fields(res_set);
    total_rows  = mysql_num_rows(res_set);
  char ***out   = malloc(sizeof(char**) * total_rows );
  apop_data *out= apop_data_alloc(0,0,0);
    while ((row = mysql_fetch_row (res_set)) ) {
		out[currentrow]	= malloc(sizeof(char*) * total_cols);
		for (jj=0;jj<total_cols;jj++){
			if (row[jj]==NULL){
				out[currentrow][jj]	= malloc(sizeof(char));
				strcpy(out[currentrow][jj], "\0");
			} else {
				out[currentrow][jj]	= malloc(1+strlen(row[jj]));
				strcpy(out[currentrow][jj], row[jj]);
			}
		}
		currentrow++;
    }
    output->textegories  = out;
    output->catsize[0]  = total_rows;
    output->catsize[1]  = total_cols;
    output->textsize[0]  = total_rows;
    output->textsize[1]  = total_cols;
    if (mysql_errno (conn)){
         print_error (conn, "mysql_fetch_row() failed");
         return NULL;
    } 
    return out;
}

apop_data * apop_mysql_query_to_text(char *query){
  MYSQL_RES *res_set;
  apop_data *output;
    if (mysql_query (mysql_db, query)){
        print_error (mysql_db, "mysql_query() failed");
        return NULL;
    }
    res_set = mysql_store_result (mysql_db);     // generate result set
    if (!res_set){
       print_error (mysql_db, "mysql_store_result() failed");
       return NULL;
    }
    // process result set, and then deallocate it 
    output = process_result_set_chars (mysql_db, res_set);
    mysql_free_result (res_set);
    return output;
}
