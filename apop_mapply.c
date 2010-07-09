/** \file apop_mapply.c vector/matrix map/apply.  */
/* Copyright (c) 2007, 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2. 
 
   This file is a tour de force of the if statement. There are several possibilities:

   --user wants the vector, matrix, rows, columns, all items in the data set
   --user wants output in a new location, or written to the old.
   --user has extra parameters
   --user needs to know the index of the function
   --user wants the sum of the result (e.g., to find how many elements are NAN, or a sum of log-likelihoods).

   Further, Apophenia v0.22 introduced variadic, optional arguments, so
   we have a somewhat more robust syntax post-22, and also the prior syntax,
   which some may find useful.

   We thus have a lot of functions that all feed in to mapply_core, which, after several if statements,
   then dispatches segments do different threads, and either the vectorloop or forloop that does all the
   actual math.

 */
#include "stats.h"
#include "mapply.h"
#include "types.h"
#include "internal.h"
#include "conversions.h"
static gsl_vector*mapply_core(gsl_matrix *m, gsl_vector *vin, void *fn, gsl_vector *vout, int use_index, int use_param,void *param, char post_22);

typedef double apop_fn_v(gsl_vector*);
typedef void apop_fn_vtov(gsl_vector*);
typedef double apop_fn_d(double);
typedef void apop_fn_dtov(double*);
typedef double apop_fn_r(apop_data*);
typedef double apop_fn_vp(gsl_vector*, void *);
typedef double apop_fn_dp(double, void *);
typedef double apop_fn_rp(apop_data*, void *);
typedef double apop_fn_vpi(gsl_vector*, void *, int);
typedef double apop_fn_dpi(double, void *, int);
typedef double apop_fn_rpi(apop_data*, void *, int);
typedef double apop_fn_vi(gsl_vector*, int);
typedef double apop_fn_di(double, int);
typedef double apop_fn_ri(apop_data*, int);


/**
  Apply a function to every element of a data set, matrix or vector; or, apply a vector-taking function to every row or column of a matrix.

There are a lot of options: your function could take any combination of a \c gsl_vector/\c double, a parameter set, and the position of the element in the vector or matrix. As such, the function takes eight function inputs, one for each combination of vector/matrix, params/no params, index/no index. Fortunately, because 
this function uses the \ref designated syntax for inputs, you need specify only one.

For example, here is a function that will cut off each element of the input data to between \f$(-1, +1)\f$.
\code
double cutoff(double in, void *limit_in){ 
    double *limit = limit_in;
    return GSL_MAX(-*limit, GSL_MIN(*limit, in)); 
}

double param = 1;
apop_map(your_data, .fn_di=cutoff, .param=&param, .inplace=1);
\endcode

\param fn_v A function of the form <tt>double your_fn(gsl_vector *in)</tt>
\param fn_d A function of the form <tt>double your_fn(double in)</tt>
\param fn_r A function of the form <tt>double your_fn(apop_data *in)</tt>
\param fn_vp A function of the form <tt>double your_fn(gsl_vector *in, void *param)</tt>
\param fn_dp A function of the form <tt>double your_fn(double in, void *param)</tt>
\param fn_rp A function of the form <tt>double your_fn(apop_data *in, void *param)</tt>
\param fn_vpi A function of the form <tt>double your_fn(gsl_vector *in, void *param, int index)</tt>
\param fn_dpi A function of the form <tt>double your_fn(double in, void *param, int index)</tt>
\param fn_rpi A function of the form <tt>double your_fn(apop_data *in, void *param, int index)</tt>
\param fn_vi A function of the form <tt>double your_fn(gsl_vector *in, int index)</tt>
\param fn_di A function of the form <tt>double your_fn(double in, int index)</tt>
\param fn_ri A function of the form <tt>double your_fn(apop_data *in, int index)</tt>

\param in   The input data set. If \c NULL, I'll return \c NULL immediately.
\param param   A pointer to the parameters to be passed to those function forms taking a \c *param.

\param part Which part of the \c apop_data struct should I use?<br>
'v'==Just the vector<br>
'm'==Every element of the matrix, in turn<br>
'a'==Both 'v' and 'm'<br>
'r'==Apply a function \c gsl_vector \f$\to\f$ \c double to each row of the  matrix<br>
'c'==Apply a function \c gsl_vector \f$\to\f$ \c double to each column of the  matrix<br>
Default is 'a', but notice that I'll ignore a \c NULL vector or matrix, so if your data set has only a vector (for example), that's what I'll use.

\li The function forms with <tt>r</tt> in them, like \c fn_ri, are row-by-row. I'll use
\ref Apop_data_row to get each row in turn, and send it to the function. The first
implication is that your function should be expecting a \ref apop_data set with exactly
one row in it. The second is that \c part is ignored: it only makes sense to go row-by-row.

\param inplace  If zero, generate a new \ref apop_data set for output, which will contain the mapped values (and the names from the original set). If one, modify in place. The \c double \f$\to\f$ \c double versions, \c 'v', \c 'm', and \c 'a', write to exactly the same location as before. The \c gsl_vector \f$\to\f$ \c double versions, \c 'r', and \c 'c', will write to the vector. Be careful: if you are writing in place and there is already a vector there, then the original vector is lost. (Default = 0)
\param all_pages If \c 'y', then I follow the \c more pointer to subsequent pages, else I
handle only the first page of data. Default: \c 'n'.

\ingroup mapply
*/
APOP_VAR_HEAD apop_data* apop_map(apop_data *in, apop_fn_d *fn_d, apop_fn_v *fn_v, apop_fn_r *fn_r, apop_fn_dp *fn_dp, apop_fn_vp *fn_vp, apop_fn_rp *fn_rp,  apop_fn_dpi *fn_dpi, apop_fn_vpi *fn_vpi, apop_fn_rpi *fn_rpi, apop_fn_di *fn_di,  apop_fn_vi *fn_vi, apop_fn_ri *fn_ri, void *param, int inplace, char part, int all_pages){ 
    apop_data * apop_varad_var(in, NULL)
    if (!in) return NULL;
    apop_fn_v * apop_varad_var(fn_v, NULL)
    apop_fn_d * apop_varad_var(fn_d, NULL)
    apop_fn_r * apop_varad_var(fn_r, NULL)
    apop_fn_vp * apop_varad_var(fn_vp, NULL)
    apop_fn_dp * apop_varad_var(fn_dp, NULL)
    apop_fn_rp * apop_varad_var(fn_rp, NULL)
    apop_fn_vpi * apop_varad_var(fn_vpi, NULL)
    apop_fn_dpi * apop_varad_var(fn_dpi, NULL)
    apop_fn_rpi * apop_varad_var(fn_rpi, NULL)
    apop_fn_vi * apop_varad_var(fn_vi, NULL)
    apop_fn_di * apop_varad_var(fn_di, NULL)
    apop_fn_ri * apop_varad_var(fn_ri, NULL)
    int apop_varad_var(inplace, 0)
    void * apop_varad_var(param, NULL)
    char apop_varad_var(part, 'a')
    int apop_varad_var(all_pages, 'n')
APOP_VAR_ENDHEAD
    int use_param = (fn_vp || fn_dp || fn_rp || fn_vpi || fn_rpi || fn_dpi);
    int use_index  = (fn_vi || fn_di || fn_ri || fn_vpi || fn_rpi|| fn_dpi);
    //Give me the first non-null input function.
    void *fn = fn_v ? (void *)fn_v : fn_d ? (void *)fn_d : fn_r ? (void *)fn_r : fn_vp ? (void *)fn_vp : fn_dp ? (void *)fn_dp :fn_rp ? (void *)fn_rp : fn_vpi ? (void *)fn_vpi : fn_rpi ? (void *)fn_rpi: fn_dpi ? (void *)fn_dpi : fn_vi ? (void *)fn_vi : fn_di ? (void *)fn_di : fn_ri ? (void *)fn_ri : NULL;

    int by_apop_rows = fn_r || fn_rp || fn_rpi || fn_ri;

    apop_assert_s(!((part=='c' || part=='r') && (fn_d || fn_dp || fn_dpi || fn_di)), 
        "You asked for a vector-oriented operation (.part='r' or .part='c'), but gave me a "
        "scalar-oriented function. Did you mean part=='a'?");

    //Allocate output
    Get_vmsizes(in); //vsize, msize1, msize2
    apop_data *out;
    if (inplace)
       out = in;
    else 
         out = by_apop_rows ? apop_data_copy(in)
             : part == 'v' || (in->vector && ! in->matrix) ? apop_data_alloc(vsize, 0, 0)
             : part == 'm' ? apop_data_alloc(0, msize1, msize2)
             : part == 'a' ? apop_data_alloc(vsize, msize1, msize2)
             : part == 'r' ? apop_data_alloc(msize1, 0, 0)
             : part == 'c' ?  apop_data_alloc(msize2, 0, 0) : NULL;
    if (in->names){
        if (part == 'v'  || (in->vector && ! in->matrix)) {
             apop_name_stack(out->names, in->names, 'v');
             apop_name_stack(out->names, in->names, 'r');
        }
        else if (part == 'm'){
             apop_name_stack(out->names, in->names, 'r');
             apop_name_stack(out->names, in->names, 'c');
        }
        else if (by_apop_rows || part == 'a')
             out->names = apop_name_copy(in->names);
        else if (part == 'r')
             apop_name_stack(out->names, in->names, 'r');
        else if (part == 'c')
            apop_name_stack(in->names, out->names, 'r', 'c');
    }

    //Call mapply_core.
    if (by_apop_rows){
        for (int i=0; i< vsize? vsize : msize1; i++){
            Apop_data_row(out, i, the_row);
            if (fn_r) fn_r(the_row);
            else if (fn_rp)
                fn_rp(the_row, param);
            else if (fn_rpi)
                fn_rpi(the_row, param, i);
            else if (fn_ri)
                fn_ri(the_row, i);
        }
    } else {
        if (in->vector && (part == 'v' || part=='a'))
            mapply_core(NULL, in->vector, fn, out->vector, use_index, use_param, param, 'r');
        if (in->matrix && (part == 'm' || part=='a')){
            int smaller_dim = GSL_MIN(in->matrix->size1, in->matrix->size2);
            for (int i=0; i< smaller_dim; i++){
                if (smaller_dim == in->matrix->size1){
                    Apop_row(in, i, onevector);
                    Apop_row(out, i, twovector);
                    mapply_core(NULL, onevector, fn, twovector, use_index, use_param, param, 'r');
                }else{
                    Apop_col(in, i, onevector);
                    Apop_col(out, i, twovector);
                    mapply_core(NULL, onevector, fn, twovector, use_index, use_param, param, 'c');
                }
            }
        }
        if (part == 'r' || part == 'c')
            mapply_core(in->matrix, NULL, fn, out->vector, use_index, use_param, param, part);
    }
    if (all_pages && in->more)
        out->more = apop_map_base(in->more, fn_d, fn_v, fn_r, fn_dp, fn_vp, fn_rp, fn_dpi, fn_vpi, fn_rpi, fn_di, fn_vi, fn_ri, param, inplace, part, all_pages);
    return out;
}

/** A convenience function to call \ref apop_map, and return the sum of the resulting elements. Thus, this function returns a single \c double. See the \ref apop_map page for details of the inputs, which are the same here, except that \c inplace doesn't make sense---this function will always internally allocate a temp data set and free it before returning.

  See also the \ref mapply "map/apply page" for details.
 \ingroup mapply
 */
APOP_VAR_HEAD double apop_map_sum(apop_data *in, apop_fn_d *fn_d, apop_fn_v *fn_v, apop_fn_r *fn_r, apop_fn_dp *fn_dp, apop_fn_vp *fn_vp, apop_fn_rp *fn_rp,   apop_fn_dpi *fn_dpi,  apop_fn_vpi *fn_vpi, apop_fn_rpi *fn_rpi, apop_fn_di *fn_di,  apop_fn_vi *fn_vi, apop_fn_ri *fn_ri,    void *param, char part, int all_pages){ 
    apop_data * apop_varad_var(in, NULL)
    apop_fn_v * apop_varad_var(fn_v, NULL)
    apop_fn_d * apop_varad_var(fn_d, NULL)
    apop_fn_r * apop_varad_var(fn_r, NULL)
    apop_fn_vp * apop_varad_var(fn_vp, NULL)
    apop_fn_dp * apop_varad_var(fn_dp, NULL)
    apop_fn_rp * apop_varad_var(fn_rp, NULL)
    apop_fn_vpi * apop_varad_var(fn_vpi, NULL)
    apop_fn_dpi * apop_varad_var(fn_dpi, NULL)
    apop_fn_rpi * apop_varad_var(fn_rpi, NULL)
    apop_fn_vi * apop_varad_var(fn_vi, NULL)
    apop_fn_di * apop_varad_var(fn_di, NULL)
    apop_fn_ri * apop_varad_var(fn_ri, NULL)
    void * apop_varad_var(param, NULL)
    char apop_varad_var(part, 'a')
    int apop_varad_var(all_pages, 'n')
APOP_VAR_ENDHEAD 
    if (fn_r || fn_rp){
        Get_vmsizes(in);
        apop_data *copy = apop_data_copy(in);
        double outsum = 0;
        for (int i = 0; i< (vsize ? vsize : msize1); i++){
            Apop_data_row(copy, i, r);
            outsum += fn_r ? fn_r(r) : fn_rp(r, param);
        }
        return outsum;
    }
    apop_data *out = apop_map(in, fn_d, fn_v, fn_r, fn_dp, fn_vp, fn_rp, fn_dpi, fn_vpi, fn_rpi, fn_di, fn_vi, fn_ri, param, 0, part);
    double outsum = (out->vector ? apop_sum(out->vector) : 0)
                     + (out->matrix ? apop_matrix_sum(out->matrix) : 0);
    apop_data_free(out);
    return outsum + 
                ((all_pages && in->more) ? apop_map_sum_base(in->more, fn_d, fn_v, fn_r, fn_dp, fn_vp, fn_rp, fn_dpi, fn_vpi, fn_rpi, fn_di, fn_vi, fn_ri, param, part, all_pages) : 0);
}

typedef struct {
    size_t      *limlist;
    void        *fn;
    gsl_matrix  *m;
    gsl_vector  *v, *vin;
    int use_index, use_param;
    char rc;
    void *param;
} threadpass;

static void *forloop(void *t){
  threadpass      *tc = t;
  apop_fn_v    *vtod=tc->fn;
  apop_fn_vp  *fn_vp=tc->fn;
  apop_fn_vpi *fn_vpi=tc->fn;
  apop_fn_vi  *fn_vi=tc->fn;
  gsl_vector    view;
  double        val;
    for (int i= tc->limlist[0]; i< tc->limlist[1]; i++){
        view    = tc->rc == 'r' ? gsl_matrix_row(tc->m, i).vector : gsl_matrix_column(tc->m, i).vector;
        val     = 
        tc->use_param ? (tc->use_index ? fn_vpi(&view, tc->param, i) : 
                                     fn_vp(&view, tc->param) )
                      : (tc->use_index ? fn_vi(&view, i) : 
                                     vtod(&view) );
        gsl_vector_set(tc->v, i, val);
    }
  return NULL;
}

static void *oldforloop(void *t){
  threadpass      *tc = t;
  apop_fn_vtov    *vtov=tc->fn;
    if (tc->v){
        tc->rc = 'r';
        return forloop(t);
    }
    for (int i= tc->limlist[0]; i< tc->limlist[1]; i++){
        Apop_matrix_row(tc->m, i, v);
        vtov(v);
    }
  return NULL;
}

//if mapping to self, then set tc.v = in_v
static void *vectorloop(void *t){
  threadpass      *tc = t;
  double          inval, outval;
  apop_fn_d    *dtod=tc->fn;
  apop_fn_dp  *fn_dp=tc->fn;
  apop_fn_dpi *fn_dpi=tc->fn;
  apop_fn_di  *fn_di=tc->fn;
    for (int i= tc->limlist[0]; i< tc->limlist[1]; i++){
        inval   = gsl_vector_get(tc->vin, i);
        outval =
        tc->use_param ? (tc->use_index ? fn_dpi(inval, tc->param, i) : 
                                     fn_dp(inval, tc->param))
                     : (tc->use_index ? fn_di(inval, i) : 
                                     dtod(inval));
        gsl_vector_set(tc->v, i, outval);
    }
  return NULL;
}

static void *oldvectorloop(void *t){
  threadpass      *tc = t;
  double          *inval;
  apop_fn_dtov    *dtov=tc->fn;
    if (tc->v)
        return vectorloop(t);
    for (int i= tc->limlist[0]; i< tc->limlist[1]; i++){
        inval   = gsl_vector_ptr(tc->vin, i);
        dtov(inval);
    }
  return NULL;
}

static size_t *threadminmax(const int threadno, const int totalct, const int threadct){
  int       segment_size        = totalct/threadct;
  size_t    *out                = malloc(sizeof(int)*3);
        out[0]  = threadno*segment_size;
        out[1]  = (threadno==threadct-1) ? totalct : 1+(threadno+1)*segment_size;
        out[2]  = threadno;
        return out;
}

static gsl_vector*mapply_core(gsl_matrix *m, gsl_vector *vin, void *fn, gsl_vector *vout, int use_index, int use_param,void *param, char post_22){
  int           threadct    = apop_opts.thread_count;
  pthread_t     thread_id[threadct];
  int           i;
  threadpass    tp[threadct];
    for (i=0 ; i<threadct; i++)
        tp[i] = (threadpass) {
            .limlist   = threadminmax(i, 
                    m? ((!post_22 || post_22 == 'r') ? m->size1 : m->size2) : vin->size ,threadct),
            .fn = fn,   .m = m, 
            .vin = vin, .v = vout,
            .use_index = use_index, .use_param= use_param,
            .param = param, .rc = post_22
        };
    if (threadct==1){ //don't thread.
        if (m)  post_22 ? forloop(tp) : oldforloop(tp);
        else    post_22 ? vectorloop(tp) : oldvectorloop(tp);
    } else {
        for (i=0 ; i<threadct; i++){
            if (m)  pthread_create(&thread_id[i], NULL,post_22 ? forloop : oldforloop,(tp+i));
            else    pthread_create(&thread_id[i], NULL,post_22 ? vectorloop : oldvectorloop,(tp+i));
        }
        for (i=0 ; i<threadct; i++)
            pthread_join(thread_id[i], NULL);
    }
    for (i=0 ; i<threadct; i++)
        free(tp[i].limlist);
    return vout;
}

/** \defgroup mapply Map or apply a function to a vector or matrix

These functions will pull each element of a vector or matrix, or each row of a matrix, and apply a function to the given element. See the data->map/apply section of the \ref outline_mapply "outline" for many examples. 

Here are a few technical details of usage:

\li If \c apop_opts.thread_count is greater than one, then the matrix will be broken into chunks and each sent to a different thread. Notice that the GSL is generally threadsafe, and SQLite is threadsafe conditional on several commonsense caveats that you'll find in the SQLite documentation.

\li See also \ref apop_map and \ref apop_map_sum, which take in \ref apop_data sets and have options for sending in the index of the item or additional parameters.

\li The \c ...sum functions are convenience functions that just call \c ...map and then add up the contents. Thus, you will need to have adequate memory for the allocation of the temp matrix/vector.
\{ */

/** Map a function onto every row of a matrix.  The function that you input takes in a gsl_vector and returns a \c double. \c apop_matrix_map will produce a vector view of each row, and send each row to your function. It will output a \c gsl_vector holding your function's output for each row.


  \param m  The matrix
  \param fn A function of the form <tt>double fn(gsl_vector* in)</tt>

  \return A \c gsl_vector with the corresponding value for each row.

  See \ref mapply "the map/apply page" for details.
 */
gsl_vector *apop_matrix_map(const gsl_matrix *m, double (*fn)(gsl_vector*)){
  gsl_vector    *out        = gsl_vector_alloc(m->size1);
    return mapply_core((gsl_matrix*) m, NULL, fn, out, 0, 0, NULL, 0);
}

/** Apply a function to every row of a matrix.  The function that you input takes in a gsl_vector and returns nothing. \c apop_matrix_apply will produce a vector view of each row, and send each row to your function.

  \param m  The matrix
  \param fn A function of the form <tt>void fn(gsl_vector* in)</tt>

  See \ref mapply "the map/apply page" for details.
 */
void apop_matrix_apply(gsl_matrix *m, void (*fn)(gsl_vector*)){
    mapply_core(m, NULL, fn, NULL, 0, 0, NULL, 0);
}

/** Map a function onto every element of a vector.  The function that you input takes in a \c double and returns a \c double. \c apop_apply will send each element to your function, and will output a \c gsl_vector holding your function's output for each row.

  \param v  The input vector
  \param fn A function of the form <tt>double fn(double in)</tt>

  \return A \c gsl_vector (allocated by this function) with the corresponding value for each row.

  See \ref mapply "the map/apply page" for details.
 */
gsl_vector *apop_vector_map(const gsl_vector *v, double (*fn)(double)){
  gsl_vector    *out        = gsl_vector_alloc(v->size);
    return mapply_core(NULL, (gsl_vector*) v, fn, out, 0, 0, NULL, 0);
}

/** Apply a function to every row of a matrix.  The function that you input takes in a gsl_vector and returns nothing. \c apop_apply will
 send a pointer to each element of your vector to your function.

  \param v  The input vector
  \param fn A function of the form <tt>void fn(double in)</tt>

  See \ref mapply "the map/apply page" for details.
 */
void apop_vector_apply(gsl_vector *v, void (*fn)(double*)){
    mapply_core(NULL, v, fn, NULL, 0, 0, NULL, 0); }

static void apop_matrix_map_all_vector_subfn(const gsl_vector *in, gsl_vector *outv, double (*fn)(double)){
    mapply_core(NULL, (gsl_vector *) in, fn, outv, 0, 0, NULL, 0); }

/** Maps a function to every element in a matrix (as opposed to every row)

  \param in The matrix whose elements will be inputs to the function
  \param fn A function with a form like <tt>double f(double in)</tt>.
  \return a matrix of the same size as the original, with the function applied.

  See \ref mapply "the map/apply page" for details.
  */

gsl_matrix * apop_matrix_map_all(const gsl_matrix *in, double (*fn)(double)){
  gsl_matrix *out = gsl_matrix_alloc(in->size1, in->size2);
    for (size_t i=0; i< in->size1; i++){
        gsl_vector_const_view inv = gsl_matrix_const_row(in, i);
        Apop_matrix_row(out, i, v);
        apop_matrix_map_all_vector_subfn(&inv.vector, v, fn);
    }
    return out;
}

/** Applies a function to every element in a matrix (as opposed to every row)

  \param in The matrix whose elements will be inputs to the function
  \param fn A function with a form like <tt>void f(double *in)</tt> which will modify the data at the in-pointer in place.
  See \ref mapply "the map/apply page" for details.
  */
void apop_matrix_apply_all(gsl_matrix *in, void (*fn)(double *)){
    for (size_t i=0; i< in->size1; i++){
        Apop_matrix_row(in, i, v);
        apop_vector_apply(v, fn);
    }
}

/** Like \c apop_vector_map, but returns the sum of the resulting mapped function. For example, <tt>apop_vector_map_sum(v, isnan)</tt> returns the number of elements of <tt>v</tt> that are \c NaN.

  See \ref mapply "the map/apply page" for details.  */
double apop_vector_map_sum(const gsl_vector *in, double(*fn)(double)){
    gsl_vector *m = apop_vector_map (in, fn);
    double out = apop_vector_sum(m);
    gsl_vector_free(m);
    return out;
}

/** Like \c apop_matrix_map_all, but returns the sum of the resulting mapped function. For example, <tt>apop_matrix_map_all_sum(v, isnan)</tt> returns the number of elements of <tt>m</tt> that are \c NaN.

  See \ref mapply "the map/apply page" for details.  */
double apop_matrix_map_all_sum(const gsl_matrix *in, double (*fn)(double)){
    gsl_matrix *m = apop_matrix_map_all (in, fn);
    double out = apop_matrix_sum(m);
    gsl_matrix_free(m);
    return out;
}

/** Like \c apop_matrix_map, but returns the sum of the resulting mapped vector. For example, let \c log_like be a function that returns the log likelihood of an input vector; then <tt>apop_matrix_map_sum(m, log_like)</tt> returns the total log likelihood of the rows of \c m.

  See \ref mapply "the map/apply page" for details.  */
double apop_matrix_map_sum(const gsl_matrix *in, double (*fn)(gsl_vector*)){
    gsl_vector *v = apop_matrix_map (in, fn);
    double out = apop_vector_sum(v);
    gsl_vector_free(v);
    return out;
}
/** \} */
