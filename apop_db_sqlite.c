/** \file apop_db_sqlite.c
This file is included directly into \ref apop_db.c.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
 */
#include <sqlite3.h>
#include <string.h>

sqlite3	*db=NULL;	                //There's only one SQLite database handle. Here it is.


/** \page db_moments Database moments (plus pow()!)

SQLite lets users define new functions for use in queries, and Apophenia uses this facility to define a few common functions.

\li <tt>select ran() from table</tt> will produce a new random number between zero and one for every row of the input table, using \c gsl_rng_uniform. 

\li The SQL standard includes the <tt>count(x)</tt> and <tt>avg(x)</tt> aggregators,
but statisticians are usually interested in higher moments as well---at least the
variance. Therefore, SQL queries using the Apophenia library may include any of these moments:

\code
select count(x), stddev(x), avg(x), var(x), variance(x), skew(x), kurt(x), kurtosis(x),
std(x), stddev_samp(x), stddev_pop(x), var_samp(x), var_pop(x)
from table
group by whatever
\endcode

<tt>var</tt> and <tt>variance</tt>; <tt>kurt</tt> and <tt>kurtosis</tt> do the same thing. Choose the one that sounds better to you. <tt>var</tt>, <tt>var_samp</tt>, <tt>stddev</tt> and <tt>stddev_samp</tt> give sample variance/standard deviation; <tt>variance</tt>, <tt>var_pop</tt> <tt>std</tt> and <tt>stddev_pop</tt> give population standard deviation.  The plethora of variants are for mySQL compatibility.

\li The  var/skew/kurtosis functions calculate sample moments, so if you want the population moment, multiply the result by (n-1)/n .

\li Also provided: wrapper functions for standard math library
functions---<tt>sqrt(x)</tt>, <tt>pow(x,y)</tt>, <tt>exp(x)</tt>, <tt>log(x)</tt>,
and trig functions. They call the standard math library function of the same name
to calculate \f$\sqrt{x}\f$, \f$x^y\f$, \f$e^x\f$, \f$\ln(x)\f$, \f$\sin(x)\f$,
\f$\arcsin(x)\f$, et cetera.

\li The <tt>ran()</tt> function calls <tt>gsl_rng_uniform</tt> to produce a uniform
draw between zero and one. It keeps its own <tt>gsl_rng</tt>, which is intialized on
first call using the value of <tt>apop_ots.rng_seed</tt> (which is then incremented,
so the next function to use it will get a different seed).

\code
select sqrt(x), pow(x,0.5), exp(x), log(x), 
    sin(x), cos(x), tan(x), asin(x), acos(x), atan(x)
from table
\endcode

Here is a test script using many of the above.

\include db_fns.c

Here is some more realistic sample code:

\include normalizations.c
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
    if (argc<1) return;
    StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
    if (p && argv[0]){
        double x = sqlite3_value_double(argv[0]);
        double ratio = p->cnt/(p->cnt+1.0);
        p->cnt++;
        p->avg	*= ratio;
        p->avg2	*= ratio;
        p->avg += x/(p->cnt +0.0);
        p->avg2 += gsl_pow_2(x)/(p->cnt +0.0);
    }
}

static void threeStep(sqlite3_context *context, int argc, sqlite3_value **argv){
    if (argc<1) return;
    StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
    if (p && argv[0]){
        double x = sqlite3_value_double(argv[0]);
        double ratio =  p->cnt/(p->cnt+1.0);
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
    if( argc<1 ) return;
    StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
    if (p && argv[0]){
        double x = sqlite3_value_double(argv[0]);
        p->cnt++;
        p->avg = (x + p->avg * (p->cnt-1.))/p->cnt;
        p->avg2 = (gsl_pow_2(x)+ p->avg2 * (p->cnt-1.))/p->cnt;
        p->avg3 = (gsl_pow_3(x)+ p->avg3 * (p->cnt-1.))/p->cnt;
        p->avg4 = (gsl_pow_4(x)+ p->avg4 * (p->cnt-1.))/p->cnt;
    }
}

static void stdDevFinalizePop(sqlite3_context *context){
    StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
    if (p && p->cnt>1)
      sqlite3_result_double(context, sqrt((p->avg2 - gsl_pow_2(p->avg))));
    else if (p->cnt == 1)
      	sqlite3_result_double(context, 0);
}

static void varFinalizePop(sqlite3_context *context){
  StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
    if( p && p->cnt>1 )
        sqlite3_result_double(context, (p->avg2 - gsl_pow_2(p->avg)));
    else if (p->cnt == 1)
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
                        + 2*rCnt * gsl_pow_3(p->avg)) * rCnt/((rCnt-1.0)*(rCnt-2.0)));
    } else if (p->cnt == 1)
      	sqlite3_result_double(context, 0);
}

static void kurtFinalize(sqlite3_context *context){
    StdDevCtx *p = sqlite3_aggregate_context(context, sizeof(*p));
    if( p && p->cnt>1 ){
      double n = p->cnt;
      double kurtovern = p->avg4 - 4*p->avg3*p->avg
                        + 6 * p->avg2*gsl_pow_2(p->avg)
                        - 3* gsl_pow_4(p->avg);
      double var = p->avg2 - gsl_pow_2(p->avg);
      long double coeff0= n*n/(gsl_pow_3(n)*(gsl_pow_2(n)-3*n+3));
      long double coeff1= n*gsl_pow_2(n-1)+ (6*n-9);
      long double coeff2= n*(6*n-9);
      sqlite3_result_double(context, coeff0*(coeff1 * kurtovern + coeff2 * gsl_pow_2(var)));
    } else if (p->cnt == 1)
      sqlite3_result_double(context, 0);
}

static void powFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    double base = sqlite3_value_double(argv[0]);
    double exp  = sqlite3_value_double(argv[1]);
    sqlite3_result_double(context, pow(base, exp));
}

static void rngFn(sqlite3_context *context, int argc, sqlite3_value **argv){
    Staticdef(gsl_rng *, rng, apop_rng_alloc(apop_opts.rng_seed++));
    sqlite3_result_double(context, gsl_rng_uniform(rng));
}

#define sqfn(name) static void name##Fn(sqlite3_context *context, int argc, sqlite3_value **argv){ \
    sqlite3_result_double(context, name(sqlite3_value_double(argv[0]))); }

sqfn(sqrt) sqfn(exp) sqfn(log) sqfn(log10) sqfn(sin) 
sqfn(cos) sqfn(tan) sqfn(asin) sqfn(acos) sqfn(atan)


static int apop_sqlite_db_open(char const *filename){
    int status = sqlite3_open(filename ? filename : ":memory:", &db);
    Apop_stopif(status, db=NULL; return status,
            0, "The database %s didn't open.", filename ? filename : "in memory");
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
	sqlite3_create_function(db, "ln", 1, SQLITE_ANY, NULL, &logFn, NULL, NULL);
	sqlite3_create_function(db, "ran", 0, SQLITE_ANY, NULL, &rngFn, NULL, NULL);
	sqlite3_create_function(db, "pow", 2, SQLITE_ANY, NULL, &powFn, NULL, NULL);

#define sqlink(name) sqlite3_create_function(db, #name , 1, SQLITE_ANY, NULL, &name##Fn, NULL, NULL);
    sqlink(sqrt) sqlink(exp) sqlink(sin) sqlink(cos)
    sqlink(tan) sqlink(asin) sqlink(acos) sqlink(atan) sqlink(log) sqlink(log10)
	apop_query("pragma short_column_names");
    return 0;
}

typedef struct {    //for the apop_query_to_... functions.
    int       firstcall, namecol;
    size_t    currentrow;
    apop_data *outdata;
} callback_t;

//This is the callback for apop_query_to_text.
static int db_to_chars(void *qinfo,int argc, char **argv, char **column){
    callback_t *qi= qinfo;
    apop_data* d  = qi->outdata; //alias. Allocated in calling fn.
    int	addnames = 0, ncshift=0;
    if (!d->names->textct) addnames++;
    if (qi->firstcall){
        qi->firstcall = 0;
        for(int i=0; i<argc; i++)
            if (!strcasecmp(column[i], apop_opts.db_name_column)){
                qi->namecol = i;
                break;
            }
    }
    int rows = d->textsize[0];
    int cols = argc - (qi->namecol >= 0);
    apop_text_alloc(d, rows+1, cols);//doesn't move d.
    for (size_t jj=0; jj<argc; jj++)
        if (jj == qi->namecol){
            apop_name_add(d->names, argv[jj], 'r'); 
            ncshift ++;
        } else {
            apop_text_add(d, rows, jj-ncshift, (argv[jj]==NULL)? apop_opts.nan_string: argv[jj]);
            //Asprintf(&(d->text[rows][jj-ncshift]), "%s", (argv[jj]==NULL)? "NaN": argv[jj]);
            if(addnames)
                apop_name_add(d->names, column[jj], 't'); 
        }
    return 0;
}

apop_data * apop_sqlite_query_to_text(char *query){
    char *err = NULL;
    callback_t qinfo = {.outdata=apop_data_alloc(), .namecol=-1, .firstcall=1};
    if (db==NULL) apop_db_open(NULL);
    sqlite3_exec(db, query, db_to_chars, &qinfo, &err); ERRCHECK_SET_ERROR(qinfo.outdata)
    if (qinfo.outdata->textsize[0]==0){
        apop_data_free(qinfo.outdata);
        return NULL;
    }
    return qinfo.outdata;
}

typedef struct {
    apop_data  *d;
    int        intypes[5];//names, vectors, mcols, textcols, weights.
    int        current, thisrow, error_thrown;
    const char *instring;
} apop_qt;

static void count_types(apop_qt *in, const char *intypes){
    int i = 0;
    char c;
    in->instring = intypes;
    while ((c=intypes[i++]))
        if (c=='n'||c=='N')      in->intypes[0]++;
        else if (c=='v'||c=='V') in->intypes[1]++;
        else if (c=='m'||c=='M') in->intypes[2]++;
        else if (c=='t'||c=='T') in->intypes[3]++;
        else if (c=='w'||c=='W') in->intypes[4]++;
    if (in->intypes[0]>1)
        Apop_notify(1, "You asked apop_query_to_mixed data for multiple row names. I'll ignore all but the last one.");
    if (in->intypes[1]>1)
        Apop_notify(1, "You asked apop_query_to_mixed for multiple vectors. I'll ignore all but the last one.");
    if (in->intypes[4]>1)
        Apop_notify(1, "You asked apop_query_to_mixed for multiple weighting vectors. I'll ignore all but the last one.");
}

static int multiquery_callback(void *instruct, int argc, char **argv, char **column){
    apop_qt *in = instruct;
    char c;
    int thistcol    = 0, 
        thismcol    = 0,
        colct       = 0,
        i, addnames = 0;
    in->thisrow ++;
    if (!in->d) {
        in->d = in->intypes[2]
                ? apop_data_alloc(in->intypes[1], 1, in->intypes[2])
                : apop_data_alloc(in->intypes[1]);
        if (in->intypes[4])
            in->d->weights  = gsl_vector_alloc(1);
        if (in->intypes[3]){
            in->d->textsize[0]  = 1;
            in->d->textsize[1]  = in->intypes[3];
            in->d->text         = malloc(sizeof(char***));
        }
    }
    if (!(in->d->names->colct + in->d->names->textct + (in->d->names->vector!=NULL)))
        addnames++;
    if (in->d->textsize[1]){
        in->d->textsize[0]         = in->thisrow;
        in->d->text                = realloc(in->d->text, sizeof(char ***)*in->thisrow);
        in->d->text[in->thisrow-1] = malloc(sizeof(char**) * in->d->textsize[1]);
    }
    if (in->intypes[2])
        apop_matrix_realloc(in->d->matrix, in->thisrow, in->intypes[2]);
    for (i=in->current=0; i< argc; i++){
        c   = in->instring[in->current++];
        if (c=='n'||c=='N'){
            apop_name_add(in->d->names, (argv[i]? argv[i] : "NaN")  , 'r'); 
            if(addnames)
                apop_name_add(in->d->names, column[i], 'h'); 
        } else if (c=='v'||c=='V'){
            apop_vector_realloc(in->d->vector, in->thisrow);
            apop_data_set(in->d, in->thisrow-1, -1, 
                                    argv[i] ? atof(argv[i]) : GSL_NAN);
            if(addnames)
                apop_name_add(in->d->names, column[i], 'v'); 
        } else if (c=='m'||c=='M'){
            apop_data_set(in->d, in->thisrow-1, thismcol++, 
                                    argv[i] ? atof(argv[i]) : GSL_NAN);
            if(addnames)
                apop_name_add(in->d->names, column[i], 'c'); 
        } else if (c=='t'||c=='T'){
            Asprintf(&(in->d->text[in->thisrow-1][thistcol++]), "%s", 
			                        argv[i] ? argv[i] : "NaN");
            if(addnames)
                apop_name_add(in->d->names, column[i], 't'); 
        } else if (c=='w'||c=='W'){
            apop_vector_realloc(in->d->weights, in->thisrow);
            gsl_vector_set(in->d->weights, in->thisrow-1, 
                                    argv[i] ? atof(argv[i]) : GSL_NAN);
        }
        colct++;
    }
    int requested = in->intypes[0]+in->intypes[1]+in->intypes[2]+in->intypes[3]+in->intypes[4];
      Apop_stopif(colct != requested, in->error_thrown='d'; return 1, 1, 
      "you asked for %i columns in your list of types(%s), but your query produced %u columns. "
      "The remainder will be placed in the text section. Output data set's ->error element set to 'd'." , requested, in->instring, colct);
    return 0;
}

apop_data *apop_sqlite_multiquery(const char *intypes, char *query){
    Apop_stopif(!intypes, apop_return_data_error('t'), 0, "You gave me NULL for the list of input types. I can't work with that.");
    Apop_stopif(!query, apop_return_data_error('q'), 0, "You gave me a NULL query. I can't work with that.");
    char *err = NULL;
    apop_qt info = { };
    count_types(&info, intypes);
	if (!db) apop_db_open(NULL);
    sqlite3_exec(db, query, multiquery_callback, &info, &err); 
    Apop_stopif(info.error_thrown, if (!info.d) apop_data_alloc(); info.d->error='d'; return info.d,
            0, "dimension error");
    ERRCHECK_SET_ERROR(info.d)
	return info.d;
}
