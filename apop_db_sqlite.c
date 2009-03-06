/** \file apop_db_sqlite.c
This file is included directly into \ref apop_db.c. It is read only if APOP_USE_SQLITE is defined.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
 */
#include <sqlite3.h>

sqlite3	*db=NULL;	                //There's only one SQLite database handle. Here it is.


/** \page db_moments Database moments (plus pow()!)
\verbatim
select count(x), stddev(x), avg(x), var(x), variance(x), skew(x), kurt(x), kurtosis(x),
std(x), stddev_samp(x), stddev_pop(x), var_samp(x), var_pop(x)
from table
group by whatever
\endverbatim

\verbatim
select sqrt(x), pow(x,0.5), exp(x), log(x), 
    sin(x), cos(x), tan(x), asin(x), acos(x), atan(x)
from table
\endverbatim

The SQL standard includes the <tt>count(x)</tt> and <tt>avg(x)</tt> aggregators,
but statisticians are usually interested in higher moments as well---at
least the variance. Therefore, SQL queries using the Apophenia library
may include any of the moments above.

<tt>var</tt> and <tt>variance</tt>; <tt>kurt</tt> and <tt>kurtosis</tt> do the same
thing. Choose the one that sounds better to you. <tt>var</tt>, <tt>var_samp</tt>, <tt>stddev</tt> and <tt>stddev_samp</tt> give sample variance/standard deviation; <tt>variance</tt>, <tt>var_pop</tt> <tt>std</tt> and <tt>stddev_pop</tt> give population standard deviation. 
The plethora of variants are for mySQL compatibility.

The  var/skew/kurtosis functions calculate sample moments, so if you want the population moment, multiply the result by (n-1)/n .

For bonus points, there are the <tt>sqrt(x)</tt>, <tt>pow(x,y)</tt>,
<tt>exp(x)</tt>, <tt>log(x)</tt>, and trig functions. They call the standard
math library function of the same name to calculate \f$\sqrt{x}\f$,
\f$x^y\f$, \f$e^x\f$, \f$\ln(x)\f$, \f$\sin(x)\f$, \f$\arcsin(x)\f$, et cetera.

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

static void stdDevFinalizePop(sqlite3_context *context){
    StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
    if( p && p->cnt>1 ){
      sqlite3_result_double(context,
         sqrt((p->avg2 - gsl_pow_2(p->avg))));
    } else if (p->cnt == 1)
      	sqlite3_result_double(context, 0);
}

static void varFinalizePop(sqlite3_context *context){
  StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
    if( p && p->cnt>1 ){
      sqlite3_result_double(context,
         (p->avg2 - gsl_pow_2(p->avg)));
    } else if (p->cnt == 1)
      	sqlite3_result_double(context, 0);
}

static void stdDevFinalize(sqlite3_context *context){
    StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
    if( p && p->cnt>1 ){
      double rCnt = p->cnt;
      sqlite3_result_double(context,
         sqrt((p->avg2 - gsl_pow_2(p->avg))*rCnt/(rCnt-1.0)));
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
         (p->avg3*rCnt - 3*p->avg2*p->avg*rCnt 
                        + 2*rCnt * gsl_pow_3(p->avg)) * rCnt/(gsl_pow_2(rCnt)-1.0));
    } else if (p->cnt == 1)
      	sqlite3_result_double(context, 0);
}

static void kurtFinalize(sqlite3_context *context){
  StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
    if( p && p->cnt>1 ){
      double rCnt = p->cnt;
      long double scale =  gsl_pow_2(rCnt)/(gsl_pow_3(rCnt)-1.0);
      sqlite3_result_double(context, 
                 scale * ((p->avg4*rCnt - 4*p->avg3*p->avg*rCnt 
                            + 6 * p->avg2*gsl_pow_2(p->avg)*rCnt
                            - 3*rCnt* gsl_pow_4(p->avg))
                          + 6./rCnt *(p->avg2 - gsl_pow_2(p->avg)))
                );
    } else if (p->cnt == 1)
      	sqlite3_result_double(context, 0);
}

static void sqrtFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, sqrt(sqlite3_value_double(argv[0]))); }

static void powFn(sqlite3_context *context, int argc, sqlite3_value **argv){
  double  base    = sqlite3_value_double(argv[0]);
  double  exp     = sqlite3_value_double(argv[1]);
    sqlite3_result_double(context, pow(base, exp));
}

static void expFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, exp(sqlite3_value_double(argv[0]))); }

static void logFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, log(sqlite3_value_double(argv[0]))); }

static void log10Fn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, log10(sqlite3_value_double(argv[0]))); }

static void sinFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, sin(sqlite3_value_double(argv[0]))); }

static void cosFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, cos(sqlite3_value_double(argv[0]))); }

static void tanFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, tan(sqlite3_value_double(argv[0]))); }

static void asinFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, asin(sqlite3_value_double(argv[0]))); }

static void acosFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, acos(sqlite3_value_double(argv[0]))); }

static void atanFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    sqlite3_result_double(context, atan(sqlite3_value_double(argv[0]))); }

static void rngFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    if (!db_rng)
        apop_db_rng_init(0);
    sqlite3_result_double(context, gsl_rng_uniform(db_rng));
}


static int apop_sqlite_db_open(char *filename){
	if (!filename) 	sqlite3_open(":memory:",&db);
	else			sqlite3_open(filename,&db);
	if (!db)	
		{if (apop_opts.verbose)
            printf("Not sure why, but the database didn't open.\n");
		return 1; }
	sqlite3_create_function(db, "stddev", 1, SQLITE_ANY, NULL, NULL, &twoStep, &stdDevFinalize);
	sqlite3_create_function(db, "std", 1, SQLITE_ANY, NULL, NULL, &twoStep, &stdDevFinalizePop);
	sqlite3_create_function(db, "stddev_samp", 1, SQLITE_ANY, NULL, NULL, &twoStep, &stdDevFinalize);
	sqlite3_create_function(db, "stddev_pop", 1, SQLITE_ANY, NULL, NULL, &twoStep, &stdDevFinalizePop);
	sqlite3_create_function(db, "var", 1, SQLITE_ANY, NULL, NULL, &twoStep, &varFinalize);
	sqlite3_create_function(db, "var_samp", 1, SQLITE_ANY, NULL, NULL, &twoStep, &varFinalize);
	sqlite3_create_function(db, "var_pop", 1, SQLITE_ANY, NULL, NULL, &twoStep, &varFinalizePop);
	sqlite3_create_function(db, "variance", 1, SQLITE_ANY, NULL, NULL, &twoStep, &varFinalizePop);
	sqlite3_create_function(db, "skew", 1, SQLITE_ANY, NULL, NULL, &threeStep, &skewFinalize);
	sqlite3_create_function(db, "kurt", 1, SQLITE_ANY, NULL, NULL, &fourStep, &kurtFinalize);
	sqlite3_create_function(db, "kurtosis", 1, SQLITE_ANY, NULL, NULL, &fourStep, &kurtFinalize);
	sqlite3_create_function(db, "sqrt", 1, SQLITE_ANY, NULL, &sqrtFn, NULL, NULL);
	sqlite3_create_function(db, "pow", 2, SQLITE_ANY, NULL, &powFn, NULL, NULL);
	sqlite3_create_function(db, "exp", 1, SQLITE_ANY, NULL, &expFn, NULL, NULL);
	sqlite3_create_function(db, "sin", 1, SQLITE_ANY, NULL, &sinFn, NULL, NULL);
	sqlite3_create_function(db, "cos", 1, SQLITE_ANY, NULL, &cosFn, NULL, NULL);
	sqlite3_create_function(db, "tan", 1, SQLITE_ANY, NULL, &tanFn, NULL, NULL);
	sqlite3_create_function(db, "asin", 1, SQLITE_ANY, NULL, &asinFn, NULL, NULL);
	sqlite3_create_function(db, "acos", 1, SQLITE_ANY, NULL, &acosFn, NULL, NULL);
	sqlite3_create_function(db, "atan", 1, SQLITE_ANY, NULL, &atanFn, NULL, NULL);
	sqlite3_create_function(db, "log", 1, SQLITE_ANY, NULL, &logFn, NULL, NULL);
	sqlite3_create_function(db, "ln", 1, SQLITE_ANY, NULL, &logFn, NULL, NULL);
	sqlite3_create_function(db, "log10", 1, SQLITE_ANY, NULL, &log10Fn, NULL, NULL);
	sqlite3_create_function(db, "ran", 0, SQLITE_ANY, NULL, &rngFn, NULL, NULL);
	sqlite3_create_function(db, "rand", 0, SQLITE_ANY, NULL, &rngFn, NULL, NULL);
	apop_query("pragma short_column_names");
    return 0;
}

static int length_callback(void *o,int argc, char **argv, char **whatever){
    size_t *total_rows = o;
	*total_rows=atoi(argv[0]); 
	return 0;
}

//these are global for the apop_db_to_... callbacks.
int		currentrow;
int     namecol;
static int firstcall;

//This is the callback for apop_query_to_text.
static int db_to_chars(void *o,int argc, char **argv, char **column){
  int		i; 
  static int	ncfound;
  int		jj, addnames = 0, ncshift=0;
  apop_data* d  = o;
    if (!d->names->textct)
        addnames    ++;
    if (firstcall){
        namecol   = -1;
        ncfound   = 0;
        firstcall = 0;
        for(i=0; i<argc; i++)
            if (!strcmp(column[i], apop_opts.db_name_column)){
                namecol = i;
                ncfound = 1;
                break;
            }
    }
    if (!d->textsize[1])
        d->textsize[1] = argc - ncfound;
    d->text  = realloc(d->text, sizeof(char*) * ++(d->textsize[0]));
    d->text[currentrow]	= malloc(sizeof(char**) * argc);
    for (jj=0; jj<argc; jj++)
        if (jj == namecol){
            apop_name_add(d->names, argv[jj], 'r'); 
            ncshift ++;
        } else {
            asprintf(&(d->text[currentrow][jj-ncshift]), (argv[jj]==NULL)? "NaN": argv[jj]);
            if(addnames)
                apop_name_add(d->names, column[jj], 't'); 
        }
    currentrow++;
    return 0;
}

apop_data * apop_sqlite_query_to_text(char *query){
  char		*err        = NULL;
  apop_data *out        = apop_data_alloc(0, 0, 0);
    currentrow  =0;
    firstcall = 1;
    if (db==NULL) apop_db_open(NULL);
    sqlite3_exec(db, query, db_to_chars, out, &err); ERRCHECK
    if (out->textsize[0]==0){
        apop_data_free(out);
        return NULL;
    }
    return out;
}

typedef struct {
    apop_data   *d;
    int         intypes[5];//names, vectors, mcols, textcols, weights.
    int         current, thisrow;
    const char  *instring;
} apop_qt;

static void count_types(apop_qt *in, const char *intypes){
  int   i   = 0;
  char  c;
    in->instring    = intypes;
    while ((c=intypes[i++]))
        if (c=='n'||c=='N')
            in->intypes[0]++;
        else if (c=='v'||c=='V')
            in->intypes[1]++;
        else if (c=='m'||c=='M')
            in->intypes[2]++;
        else if (c=='t'||c=='T')
            in->intypes[3]++;
        else if (c=='w'||c=='W')
            in->intypes[4]++;
    if (in->intypes[0]>1)
        apop_error(1, 'c', "You asked apop_query_to_mixed data for multiple row names. I'll ignore all but the last one.\n");
    if (in->intypes[1]>1)
        apop_error(1, 'c', "You asked apop_query_to_mixed for multiple vectors. I'll ignore all but the last one.\n");
    if (in->intypes[4]>1)
        apop_error(1, 'c', "You asked apop_query_to_mixed for multiple weighting vectors. I'll ignore all but the last one.\n");
}

static int multiquery_callback(void *instruct, int argc, char **argv, char **column){
  apop_qt   *in         = instruct;
  char      c;
  int       thistcol    = 0, 
            thismcol    = 0,
            i, addnames = 0;
    if (!(in->d->names->colct + in->d->names->textct + (in->d->names->vector!=NULL)))
        addnames    ++;
    if (in->d->textsize[1])
        in->d->text[in->thisrow]   = malloc(sizeof(char**) * in->d->textsize[1]);

    for (i=in->current=0; i< argc; i++){
        c   = in->instring[in->current++];
        if (c=='n'||c=='N'){
            apop_name_add(in->d->names, (argv[i]? argv[i] : "NaN")  , 'r'); 
            if(addnames)
                apop_name_add(in->d->names, column[i], 'h'); 
        } else if (c=='v'||c=='V'){
            apop_data_set(in->d, in->thisrow, -1, 
                                    argv[i] ? atof(argv[i]) : GSL_NAN);
            if(addnames)
                apop_name_add(in->d->names, column[i], 'v'); 
        } else if (c=='m'||c=='M'){
            apop_data_set(in->d, in->thisrow, thismcol++, 
                                    argv[i] ? atof(argv[i]) : GSL_NAN);
            if(addnames)
                apop_name_add(in->d->names, column[i], 'c'); 
        } else if (c=='t'||c=='T'){
            asprintf(&(in->d->text[in->thisrow][thistcol++]),
			                        argv[i] ? argv[i] : "NaN");
            if(addnames)
                apop_name_add(in->d->names, column[i], 't'); 
        } else if (c=='w'||c=='W')
                gsl_vector_set(in->d->weights, in->thisrow, 
                                    argv[i] ? atof(argv[i]) : GSL_NAN);
    }
    in->thisrow ++;
    return 0;
}

apop_data *apop_sqlite_multiquery(const char *intypes, char *query){
  //info.intypes      =  count_types(intypes);      //names, vectors, mcols, textcols.
  char		*q2,*err    = NULL;
  apop_qt   info        = {.thisrow = 0};
  size_t	total_rows 	= 0;
    count_types(&info, intypes);
	if (db==NULL) apop_db_open(NULL);
	q2		    = malloc(strlen(query) + 300);
	apop_table_exists("apop_temp_table",1);
	sqlite3_exec(db,strcat(strcpy(q2,
		"CREATE TABLE apop_temp_table AS "), query), NULL, NULL, &err); ERRCHECK
	sqlite3_exec(db, "SELECT count(*) FROM apop_temp_table", length_callback, &total_rows, &err);
	free(query);
	ERRCHECK
	if (total_rows){
        info.d               = apop_data_alloc(info.intypes[1]*total_rows, total_rows, info.intypes[2]);
        if (info.intypes[4])
            info.d->weights  = gsl_vector_alloc(total_rows);
        info.d->textsize[0]  = total_rows;
        size_t total_cols    = apop_count_cols("apop_temp_table");
        info.d->textsize[1]  = total_cols - (info.intypes[0]+info.intypes[1]+info.intypes[2]+info.intypes[4]);
        info.d->text 	     = malloc(sizeof(char***) * total_rows);
        apop_assert(info.d->textsize[1] == info.intypes[3], NULL, 0, 'c', 
            "you asked for %i rows in your list of types, but your query produced %u columns. \
            The remainder will be placed in the text section."
            , info.intypes[0]+info.intypes[1]+info.intypes[2]+info.intypes[3]+info.intypes[4], total_cols);
		sqlite3_exec(db, "SELECT * FROM apop_temp_table", multiquery_callback, &info, &err); ERRCHECK
    }
	sqlite3_exec(db,"DROP TABLE apop_temp_table",NULL,NULL, &err);  ERRCHECK
	free(q2);
	return info.d;
}
