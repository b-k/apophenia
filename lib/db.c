//db.c  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
#include <math.h> 	//sqrt
//#include "gnulib/vasprintf.h"
#include <string.h>
#include <stdarg.h>
#include <apophenia/name.h>
#include <gsl/gsl_math.h> //GSL_NAN
#include <apophenia/db.h>

#include <apophenia/vasprintf.h>

sqlite3	*db=NULL;	//There's only one database handle. Here it is.

apop_name *last_names = NULL;	//The column names from the last query to matrix

int	total_rows, total_cols;		//the counts from the last query.

int apop_verbose	= 0;

                                                                                                                               
////////////////////////////////////////////////
// Part one: additional aggregate functions for calculating higher moments
////////////////////////////////////////////////

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
    ratio	= (p->cnt/(p->cnt+1));
    p->avg	/= ratio;
    p->avg2	/= ratio;
    p->cnt++;
    p->avg += x/p->cnt;
    p->avg2 += gsl_pow_2(x)/p->cnt;
  }
}

static void threeStep(sqlite3_context *context, int argc, sqlite3_value **argv){
StdDevCtx 	*p;
double 		x, ratio;
  if( argc<1 ) return;
  p = sqlite3_aggregate_context(context, sizeof(*p));
  if( p && argv[0] ){
    x = sqlite3_value_double(argv[0]);
    ratio	= (p->cnt/(p->cnt+1));
    p->avg	/= ratio;
    p->avg2	/= ratio;
    p->avg3	/= ratio;
    p->cnt++;
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
    ratio	= (p->cnt/(p->cnt+1));
    p->avg	/= ratio;
    p->avg2	/= ratio;
    p->avg3	/= ratio;
    p->avg4	/= ratio;
    p->avg += x/p->cnt;
    p->cnt++;
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
  }
}

static void varFinalize(sqlite3_context *context){
  StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
  if( p && p->cnt>1 ){
    double rCnt = p->cnt;
    sqlite3_result_double(context,
       (p->avg2*rCnt - p->avg*p->avg)/(rCnt-1.0));
  }
}

static void skewFinalize(sqlite3_context *context){
  StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
  if( p && p->cnt>1 ){
    double rCnt = p->cnt;
    sqlite3_result_double(context,
       (p->avg3*rCnt - 3*p->avg2*p->avg*rCnt + (3*rCnt-1) * gsl_pow_3(p->avg)) / (rCnt-1.0));
  }
}

static void kurtFinalize(sqlite3_context *context){
  StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
  if( p && p->cnt>1 ){
    double rCnt = p->cnt;
    sqlite3_result_double(context,
       (p->avg4*rCnt - 4*p->avg3*p->avg*rCnt + 6 * gsl_pow_2(p->avg2)*gsl_pow_2(p->avg)*rCnt
						- (4*rCnt+1)* gsl_pow_4(p->avg))/(rCnt-1.0));
  }
}


////////////////////////////////////////////////
// Part two: database querying functions, so the user doesn't have to
// touch sqlite3.
////////////////////////////////////////////////


int apop_open_db(char *filename){
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
	apop_query_db("pragma short_column_names");
	return 0;
}

int apop_db_open(char *filename){
	return apop_open_db(filename); }

int apop_query(const char *fmt, ...){
char 		*err, *q;
va_list		argp;
	if (db==NULL) apop_open_db(NULL);
	va_start(argp, fmt);
	vasprintf(&q, fmt, argp);
	va_end(argp);
	sqlite3_exec(db, q, NULL,NULL, &err);
	free(q);
	ERRCHECK
	return 1;
}

//an identical alias:
int apop_query_db(const char *fmt, ...){
char 		*err, *q;
va_list		argp;
	if (db==NULL) apop_open_db(NULL);
	va_start(argp, fmt);
	vasprintf(&q, fmt, argp);
	va_end(argp);
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

int apop_table_exists(char *q, int whattodo){
	//whattodo==1	==>kill table so it can be recreated in main.
	//whattodo==0	==>return error so program can continue.
char 		*err, q2[10000];
	isthere=0;
	if (db==NULL) {apop_open_db(NULL); return 0;}
	sqlite3_exec(db, "select name from sqlite_master where type='table'",tab_exists_callback,q, &err); 
	ERRCHECK
	if (whattodo==1 && isthere)
		sqlite3_exec(db,strcat(strcpy(q2, "DROP TABLE "),q),NULL,NULL, &err); ERRCHECK
	return isthere;
}


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

int apop_close_db(int vacuum){
char		*err;
	if (vacuum) sqlite3_exec(db, "VACUUM", NULL, NULL, &err);
//	ERRCHECK
	sqlite3_close(db);
	return 0;
	}

int names_callback(void *o,int argc, char **argv, char **whatever){
	apop_name_add(last_names, argv[1], 'c'); 
	return 0;
}

int length_callback(void *o,int argc, char **argv, char **whatever){
	total_rows=atoi(argv[0]); 
	return 0;
}

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

	if (db==NULL) apop_open_db(NULL);
	va_start(argp, fmt);
	vasprintf(&query, fmt, argp);
	va_end(argp);

	total_rows	= 0;
	q2		= malloc(sizeof(char)*(strlen(query)+300));
	apop_table_exists("completely_temporary_table",1);
	sqlite3_exec(db,strcat(strcpy(q2,
		"CREATE TABLE completely_temporary_table AS "),query),NULL,NULL, &err); ERRCHECK
	sqlite3_exec(db,"SELECT count(*) FROM completely_temporary_table",length_callback,NULL, &err);
	free(query);
	ERRCHECK
	if (total_rows==0){
		output	= NULL;
	} else {
		total_cols	= apop_count_cols("completely_temporary_table");
		output		= malloc(sizeof(char***) * total_rows);
		sqlite3_exec(db,"SELECT * FROM completely_temporary_table",db_to_chars,&output, &err); ERRCHECK
		if (last_names !=NULL) 
			apop_name_free(last_names); 
		last_names = apop_name_alloc();
		sqlite3_exec(db,"pragma table_info(completely_temporary_table)",names_callback, NULL, &err); ERRCHECK
	}
	sqlite3_exec(db,"DROP TABLE completely_temporary_table",NULL,NULL, &err);  ERRCHECK
	free(q2);
	return output;
}

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

	if (db==NULL) apop_open_db(NULL);
	va_start(argp, fmt);
	vasprintf(&query, fmt, argp);
	va_end(argp);
	if (apop_verbose)	printf("%s\n", query);
	total_rows= 0;
	q2	 = malloc(sizeof(char)*(strlen(query)+300));
	apop_table_exists("completely_temporary_table",1);
	sqlite3_exec(db,strcat(strcpy(q2,
		"CREATE TABLE completely_temporary_table AS "),query),NULL,NULL, &err); ERRCHECK
	sqlite3_exec(db,"SELECT count(*) FROM completely_temporary_table",length_callback,NULL, &err);
	free(query);
	ERRCHECK
	if (total_rows==0){
		output	= NULL;
	} else {
		total_cols	= apop_count_cols("completely_temporary_table");
		output		= gsl_matrix_alloc(total_rows, total_cols);
		sqlite3_exec(db,"SELECT * FROM completely_temporary_table",db_to_table,output, &err); ERRCHECK
		if (last_names !=NULL) 
			apop_name_free(last_names); 
		last_names = apop_name_alloc();
		sqlite3_exec(db,"pragma table_info(completely_temporary_table)",names_callback, NULL, &err); ERRCHECK
	}
	sqlite3_exec(db,"DROP TABLE completely_temporary_table",NULL,NULL, &err);  ERRCHECK
	free(q2);
	return output;
}

float apop_query_to_float(const char * fmt, ...){
//just like the above, but returns a single number, which will be the
//(0,0)th entry in the matrix.
gsl_matrix	*m=NULL;
va_list		argp;
char		*query;
float		out;
	va_start(argp, fmt);
	vasprintf(&query, fmt, argp);
	va_end(argp);
	m	= apop_query_to_matrix(query);
	if (m==NULL){
		printf("apop, %s, %i: query turned up a blank table. Returning zero.\n", __FILE__, __LINE__);
		return 0;
	} //else
	out	= gsl_matrix_get(m, 0, 0);
	gsl_matrix_free(m);
	return out;

}

apop_name * apop_db_get_names(void){ return last_names; }
int apop_db_get_cols(void){ return total_cols; }
int apop_db_get_rows(void){ return total_rows; }

int apop_matrix_to_db(gsl_matrix *data, char *tabname, char **headers){
int		i,j; 
int		ctr		= 0;
int		batch_size	= 100;
char		*q 		= malloc(sizeof(char)*1000);
	if (db==NULL) apop_open_db(NULL);
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
				apop_query_db(q);
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

void apop_db_merge_table(char *db_file, char *tabname){
char		***tab_list;
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
	if (*tab_list !=NULL)
		free_tab_list(&tab_list, row_ct, 1);
}

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
