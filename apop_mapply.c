
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
#include "apop_internal.h"
#include <stdbool.h>
static gsl_vector*mapply_core(apop_data *d, gsl_matrix *m, gsl_vector *vin, void *fn, gsl_vector *vout, bool use_index, bool use_param,void *param, char post_22, bool by_apop_rows);

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

There are a lot of options: your function could take any combination of a \c gsl_vector, a \c double, an \ref apop_data, a parameter set, and the position of the element in the vector or matrix. As such, the function takes twelve function inputs, one for each combination of vector/matrix, params/no params, index/no index. Fortunately, because 
this function uses the \ref designated syntax for inputs, you will specify only one.

For example, here is a function that will cut off each element of the input data to between \f$(-1, +1)\f$.
\code
double cutoff(double in, void *limit_in){ 
    double *limit = limit_in;
    return GSL_MAX(-*limit, GSL_MIN(*limit, in)); 
}

double param = 1;
apop_map(your_data, .fn_dp=cutoff, .param=&param, .inplace='y');
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

\param all_pages If \c 'y', then I follow the \c more pointer to subsequent pages, else I
handle only the first page of data. [I abuse this for an internal semaphore, by the way, so your input must always be nonnegative and  less than 1,000. Of course, 'y' and 'n' fit these rules fine.]
Default: \c 'n'. 

\param inplace  If 'n' (the default), generate a new \ref apop_data set for output,
which will contain the mapped values (and the names from the original set).<br>
If 'y',
modify in place. The \c double \f$\to\f$ \c double versions, \c 'v', \c 'm', and \c
'a', write to exactly the same location as before. The \c gsl_vector \f$\to\f$ \c
double versions, \c 'r', and \c 'c', will write to the vector. Be careful: if you
are writing in place and there is already a vector there, then the original vector is
lost.<br>
If 'v' (as in void), return \c NULL.  (Default = 'n')

\li The function forms with <tt>r</tt> in them, like \c fn_ri, are row-by-row. I'll use
\ref Apop_row to get each row in turn, and send it to the function. The first
implication is that your function should be expecting a \ref apop_data set with
exactly one row in it. The second is that \c part is ignored: it only makes sense to go
row-by-row. 

\li If you set \c inplace='y', then you will be modifying your input data set, row by row;
if you set \c inplace='n', then I will return an \ref apop_data set whose \c vector
element is as long as your data set (i.e., as long as the longest of your text, vector,
or matrix parts).

\li If you set <tt>apop_opts.thread_count</tt> to a value greater than one, I will split the data set into as many chunks as you specify, and process them simultaneously. You need to watch out for the usual hang-ups about multithreaded programming, but if your data is iid, and each row's processing is independent of the others, you should have no problems. Bear in mind that generating threads takes some small overhead, so simple cases like adding a few hundred numbers will actually be slower when threading.

\exception out->error='p' missing or mismatched parts error, such as \c NULL matrix when you sent a function acting on the matrix element.

\ingroup mapply
*/
#ifdef APOP_NO_VARIADIC
apop_data* apop_map(apop_data *in, apop_fn_d *fn_d, apop_fn_v *fn_v, apop_fn_r *fn_r, apop_fn_dp *fn_dp, apop_fn_vp *fn_vp, apop_fn_rp *fn_rp,  apop_fn_dpi *fn_dpi, apop_fn_vpi *fn_vpi, apop_fn_rpi *fn_rpi, apop_fn_di *fn_di,  apop_fn_vi *fn_vi, apop_fn_ri *fn_ri, void *param, int inplace, char part, int all_pages){
#else
apop_varad_head(apop_data*, apop_map){ 
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
    int apop_varad_var(inplace, 'n')
    void * apop_varad_var(param, NULL)
    int by_vectors = fn_v || fn_vp || fn_vpi || fn_vi;
    char apop_varad_var(part, by_vectors ? 'r' : 'a')
    int apop_varad_var(all_pages, 'n')
    return apop_map_base(in, fn_d, fn_v, fn_r, fn_dp, fn_vp, fn_rp, fn_dpi, fn_vpi, fn_rpi, fn_di, fn_vi, fn_ri, param, inplace, part, all_pages);
}

 apop_data* apop_map_base(apop_data *in, apop_fn_d *fn_d, apop_fn_v *fn_v, apop_fn_r *fn_r, apop_fn_dp *fn_dp, apop_fn_vp *fn_vp, apop_fn_rp *fn_rp,  apop_fn_dpi *fn_dpi, apop_fn_vpi *fn_vpi, apop_fn_rpi *fn_rpi, apop_fn_di *fn_di,  apop_fn_vi *fn_vi, apop_fn_ri *fn_ri, void *param, int inplace, char part, int all_pages){
#endif
    int use_param = (fn_vp || fn_dp || fn_rp || fn_vpi || fn_rpi || fn_dpi);
    int use_index  = (fn_vi || fn_di || fn_ri || fn_vpi || fn_rpi|| fn_dpi);
    //Give me the first non-null input function.
    void *fn = fn_v ? (void *)fn_v : fn_d ? (void *)fn_d : fn_r ? (void *)fn_r : fn_vp ? (void *)fn_vp : fn_dp ? (void *)fn_dp :fn_rp ? (void *)fn_rp : fn_vpi ? (void *)fn_vpi : fn_rpi ? (void *)fn_rpi: fn_dpi ? (void *)fn_dpi : fn_vi ? (void *)fn_vi : fn_di ? (void *)fn_di : fn_ri ? (void *)fn_ri : NULL;

    int by_apop_rows = fn_r || fn_rp || fn_rpi || fn_ri;

    Apop_stopif((part=='c' || part=='r') && (fn_d || fn_dp || fn_dpi || fn_di), 
                        apop_return_data_error(p),
                        0, "You asked for a vector-oriented operation (.part='r' or .part='c'), but "
                        "gave me a scalar-oriented function. Did you mean part=='a'?");

    //Allocate output
    Get_vmsizes(in); //vsize, msize1, msize2, maxsize
    apop_data *out =   (inplace=='y') ? in
                     : (inplace=='v') ? NULL
                     : by_apop_rows ? apop_data_alloc(GSL_MAX(in->textsize[0], maxsize))
                     : part == 'v' || (in->vector && ! in->matrix) ? apop_data_alloc(vsize)
                     : part == 'm' ? apop_data_alloc(msize1, msize2)
                     : part == 'a' ? apop_data_alloc(vsize, msize1, msize2)
                     : part == 'r' ? apop_data_alloc(maxsize)
                     : part == 'c' ?  apop_data_alloc(msize2) : NULL;
    Apop_stopif(inplace=='y' && (part=='r'||part=='c') && !in->vector, in->vector=gsl_vector_alloc(maxsize), 2, 
                            "No vector in your input data set for me to write outputs to; "
                            "allocating one for you of size %i", maxsize);
    if (in->names && out && !(inplace=='y')){
        if (part == 'v'  || (in->vector && ! in->matrix)) {
             apop_name_stack(out->names, in->names, 'v');
             apop_name_stack(out->names, in->names, 'r');
        }
        else if (part == 'm'){
             apop_name_stack(out->names, in->names, 'r');
             apop_name_stack(out->names, in->names, 'c');
        }
        else if (!by_apop_rows && part == 'a'){
             apop_name_free(out->names);
             out->names = apop_name_copy(in->names);
        } else if (by_apop_rows || part == 'r')
             apop_name_stack(out->names, in->names, 'r');
        else if (part == 'c')
            apop_name_stack(in->names, out->names, 'r', 'c');
    }

    if (by_apop_rows) mapply_core(in, NULL, NULL, fn, out ? out->vector : NULL, use_index, use_param, param, 'r', by_apop_rows);
    else {
        if (in->vector && (part == 'v' || part=='a'))
            mapply_core(NULL, NULL, in->vector, fn, out ? out->vector : NULL, use_index, use_param, param, 'r', by_apop_rows);
        if (in->matrix && (part == 'm' || part=='a')){
            int smaller_dim = GSL_MIN(in->matrix->size1, in->matrix->size2);
            for (int i=0; i< smaller_dim; i++){
                if (smaller_dim == in->matrix->size1){
                    Apop_row_v(in, i, onevector);
                    if (inplace=='v')
                        mapply_core(NULL, NULL, onevector, fn, NULL, use_index, use_param, param, 'r', by_apop_rows);
                    else {
                        Apop_row_v(out, i, twovector);
                        mapply_core(NULL, NULL, onevector, fn, twovector, use_index, use_param, param, 'r', by_apop_rows);
                    }
                } else {
                    Apop_col_v(in, i, onevector);
                    if (inplace=='v')
                        mapply_core(NULL, NULL, onevector, fn, NULL, use_index, use_param, param, 'c', by_apop_rows);
                    else {
                        Apop_col_v(out, i, twovector);
                        mapply_core(NULL, NULL, onevector, fn, twovector, use_index, use_param, param, 'c', by_apop_rows);
                    }
                }
            }
        }
        if (part == 'r' || part == 'c'){
            Apop_stopif(!in->matrix, if (!out) out=apop_data_alloc(); out->error='p'; return out,
                           0, "You asked for me to operate on the %cs of the matrix, but the matrix is NULL.", part);
            mapply_core(NULL, in->matrix, NULL, fn, out ? out->vector : NULL, use_index, use_param, param, part, by_apop_rows);
        }
    }
    if ((all_pages=='y' || all_pages=='Y') && in->more){
        out->more = apop_map_base(in->more, fn_d, fn_v, fn_r, fn_dp, fn_vp, fn_rp, fn_dpi, fn_vpi, fn_rpi, fn_di, fn_vi, fn_ri, param, inplace, part, all_pages);
        Apop_stopif(out->more->error, out->error=out->more->error, 1, "Error in subpage; marked parent page with same error code.");
    }
    return out;
}

typedef struct {
    void *fn;
    gsl_matrix  *m;
    gsl_vector  *v, *vin;
    apop_data *d;
    bool use_index, use_param;
    char rc;
    void *param;
} threadpass;

/* Mapply_core splits the database into an array of threadpass structs, then one of the following
  ...loop functions gets called, which does the actual for loop to step through the rows/columns/elements.  */

static void rowloop(threadpass *tc){
    apop_fn_r   *rtod=tc->fn;
    apop_fn_rp  *fn_rp=tc->fn;
    apop_fn_rpi *fn_rpi=tc->fn;
    apop_fn_ri  *fn_ri=tc->fn;
    Get_vmsizes(tc->d); //maxsize
    OMP_for (int i=0; i< maxsize; i++){
        Apop_row(tc->d, i, onerow);
        double val = 
        tc->use_param ? (tc->use_index ? fn_rpi(onerow, tc->param, i) : fn_rp(onerow, tc->param) )
                      : (tc->use_index ? fn_ri(onerow, i) : rtod(onerow) );
        if (tc->v) gsl_vector_set(tc->v, i, val);
    }
}

static void forloop(threadpass *tc){
    apop_fn_v   *vtod=tc->fn;
    apop_fn_vp  *fn_vp=tc->fn;
    apop_fn_vpi *fn_vpi=tc->fn;
    apop_fn_vi  *fn_vi=tc->fn;
    int max = tc->rc == 'r' ? tc->m->size1 : tc->m->size2;
    OMP_for (int i= 0; i< max; i++){
        gsl_vector view = tc->rc == 'r' ? gsl_matrix_row(tc->m, i).vector : gsl_matrix_column(tc->m, i).vector;
        double val  = 
            tc->use_param ? (tc->use_index ? fn_vpi(&view, tc->param, i) : fn_vp(&view, tc->param) )
                      : (tc->use_index ? fn_vi(&view, i) : vtod(&view) );
        if (tc->v) gsl_vector_set(tc->v, i, val);
    }
}

static void oldforloop(threadpass *tc){
    apop_fn_vtov *vtov=tc->fn;
    if (tc->v){
        tc->rc = 'r';
        return forloop(tc);
    }
    OMP_for (int i=0; i< tc->m->size1; i++){
        Apop_matrix_row(tc->m, i, v);
        vtov(v);
    }
}

//if mapping to self, then set tc.v = in_v
static void vectorloop(threadpass *tc){
    apop_fn_d   *dtod=tc->fn;
    apop_fn_dp  *fn_dp=tc->fn;
    apop_fn_dpi *fn_dpi=tc->fn;
    apop_fn_di  *fn_di=tc->fn;
    OMP_for (int i= 0; i< tc->vin->size; i++){
        double inval = gsl_vector_get(tc->vin, i);
        double outval =
        tc->use_param ? (tc->use_index ? fn_dpi(inval, tc->param, i) : 
                                     fn_dp(inval, tc->param))
                     : (tc->use_index ? fn_di(inval, i) : 
                                     dtod(inval));
        if (tc->v) gsl_vector_set(tc->v, i, outval);
    }
}

static void oldvectorloop(threadpass *tc){
    apop_fn_dtov *dtov=tc->fn;
    if (tc->v) return vectorloop(tc);
    OMP_for (int i= 0; i< tc->vin->size; i++){
        double *inval = gsl_vector_ptr(tc->vin, i);
        dtov(inval);
    }
}

static gsl_vector*mapply_core(apop_data *d, gsl_matrix *m, gsl_vector *vin, void *fn, gsl_vector *vout, bool use_index, bool use_param, void *param, char post_22, bool by_apop_rows){
    Get_vmsizes(d); //maxsize
    threadpass tp =
         (threadpass) {
            .fn = fn, .m = m, .d = d,
            .vin = vin, .v = vout,
            .use_index = use_index, .use_param= use_param,
            .param = param, .rc = post_22
        };
    if (by_apop_rows) rowloop(&tp);
    else if (m) post_22 ? forloop(&tp) : oldforloop(&tp);
    else        post_22 ? vectorloop(&tp) : oldvectorloop(&tp);
    return vout;
}

/** \defgroup mapply Map or apply a function to a vector or matrix

These functions will pull each element of a vector or matrix, or each row of a matrix, and apply a function to the given element. See the data->map/apply section of the \ref outline_mapply "outline" for many examples. 

There are two types, which were developed at different times. The \ref apop_map and \ref apop_map_sum functions use variadic function inputs to cover a lot of different types of process depending on the inputs. Other functions with types in their names, like \ref apop_matrix_map and \ref apop_vector_apply, may be easier to use in some cases. With one exception, they use the same guts, so use whichever type is convenient.

Here are a few technical details of usage:

\li If \c apop_opts.thread_count is greater than one, then the matrix will be broken into chunks and each sent to a different thread. Notice that the GSL is generally threadsafe, and SQLite is threadsafe conditional on several commonsense caveats that you'll find in the SQLite documentation.

\li Apart from \ref apop_map_sum (which does minimal internal allocation), the \c ...sum functions are convenience functions that just call \c ...map and then add up the contents. Thus, you will need to have adequate memory for the allocation of the temp matrix/vector.
\{ */

/** Map a function onto every row of a matrix.  The function that you input takes in a gsl_vector and returns a \c double. \c apop_matrix_map will produce a vector view of each row, and send each row to your function. It will output a \c gsl_vector holding your function's output for each row.


  \param m  The matrix
  \param fn A function of the form <tt>double fn(gsl_vector* in)</tt>

  \return A \c gsl_vector with the corresponding value for each row.

  \li If you input a \c NULL matrix, I return \c NULL.
  \li See \ref mapply "the map/apply page" for details.
 */
gsl_vector *apop_matrix_map(const gsl_matrix *m, double (*fn)(gsl_vector*)){
    if (!m) return NULL;
    gsl_vector *out = gsl_vector_alloc(m->size1);
    return mapply_core(NULL, (gsl_matrix*) m, NULL, fn, out, 0, 0, NULL, 0, false);
}

/** Apply a function to every row of a matrix.  The function that you input takes in a gsl_vector and returns nothing. \c apop_matrix_apply will produce a vector view of each row, and send each row to your function.

  \param m  The matrix
  \param fn A function of the form <tt>void fn(gsl_vector* in)</tt>

  \li If the matrix is \c NULL, this is a no-op and returns immediately.
  \li See \ref mapply "the map/apply page" for details.
 */
void apop_matrix_apply(gsl_matrix *m, void (*fn)(gsl_vector*)){
    if (!m) return;
    mapply_core(NULL, m, NULL, fn, NULL, 0, 0, NULL, 0, false);
}

/** Map a function onto every element of a vector.  The function that you input takes in a \c double and returns a \c double. \c apop_apply will send each element to your function, and will output a \c gsl_vector holding your function's output for each row.

  \param v  The input vector
  \param fn A function of the form <tt>double fn(double in)</tt>

  \return A \c gsl_vector (allocated by this function) with the corresponding value for each row.

  \li If you input a \c NULL vector, I return \c NULL.
  \li See \ref mapply "the map/apply page" for details.
 */
gsl_vector *apop_vector_map(const gsl_vector *v, double (*fn)(double)){
    if (!v) return NULL;
    gsl_vector *out = gsl_vector_alloc(v->size);
    return mapply_core(NULL, NULL, (gsl_vector*) v, fn, out, 0, 0, NULL, 0, false);
}

/** Apply a function to every row of a matrix.  The function that you input takes in a gsl_vector and returns nothing. \c apop_apply will
 send a pointer to each element of your vector to your function.

  \param v  The input vector
  \param fn A function of the form <tt>void fn(double in)</tt>

  \li If the vector is \c NULL, this is a no-op and returns immediately.
  \li See \ref mapply "the map/apply page" for details.
 */
void apop_vector_apply(gsl_vector *v, void (*fn)(double*)){
    if (!v) return;
    mapply_core(NULL, NULL, v, fn, NULL, 0, 0, NULL, 0, false); }

static void apop_matrix_map_all_vector_subfn(const gsl_vector *in, gsl_vector *outv, double (*fn)(double)){
    mapply_core(NULL, NULL, (gsl_vector *) in, fn, outv, 0, 0, NULL, 0, false); }

/** Maps a function to every element in a matrix (as opposed to every row)

  \param in The matrix whose elements will be inputs to the function
  \param fn A function with a form like <tt>double f(double in)</tt>.
  \return a matrix of the same size as the original, with the function applied.

  \li If you input a \c NULL matrix, I return \c NULL.
  \li See \ref mapply "the map/apply page" for details.
  */

gsl_matrix * apop_matrix_map_all(const gsl_matrix *in, double (*fn)(double)){
    if (!in) return NULL;
    gsl_matrix *out = gsl_matrix_alloc(in->size1, in->size2);
    OMP_for (size_t i=0; i< in->size1; i++){
        gsl_vector_const_view inv = gsl_matrix_const_row(in, i);
        Apop_matrix_row(out, i, v);
        apop_matrix_map_all_vector_subfn(&inv.vector, v, fn);
    }
    return out;
}

/** Applies a function to every element in a matrix (as opposed to every row)

  \param in The matrix whose elements will be inputs to the function
  \param fn A function with a form like <tt>void f(double *in)</tt> which will modify the data at the in-pointer in place.

  \li If the matrix is \c NULL, this is a no-op and returns immediately.
  \li See \ref mapply "the map/apply page" for details.
  */
void apop_matrix_apply_all(gsl_matrix *in, void (*fn)(double *)){
    if (!in) return;
    OMP_for (size_t i=0; i< in->size1; i++){
        Apop_matrix_row(in, i, v);
        apop_vector_apply(v, fn);
    }
}

/** Like \c apop_vector_map, but returns the sum of the resulting mapped function. For example, <tt>apop_vector_map_sum(v, isnan)</tt> returns the number of elements of <tt>v</tt> that are \c NaN.

  \li If you input a \c NULL vector, I return the sum of zero items: zero.
  \li See \ref mapply "the map/apply page" for details.  */
double apop_vector_map_sum(const gsl_vector *in, double(*fn)(double)){
    if (!in) return 0;
    gsl_vector *m = apop_vector_map (in, fn);
    double out = apop_vector_sum(m);
    gsl_vector_free(m);
    return out;
}

/** Like \c apop_matrix_map_all, but returns the sum of the resulting mapped function. For example, <tt>apop_matrix_map_all_sum(v, isnan)</tt> returns the number of elements of <tt>m</tt> that are \c NaN.

  \li If you input a \c NULL matrix, I return the sum of zero items: zero.
  \li See \ref mapply "the map/apply page" for details.  */
double apop_matrix_map_all_sum(const gsl_matrix *in, double (*fn)(double)){
    if (!in) return 0;
    gsl_matrix *m = apop_matrix_map_all (in, fn);
    double out = apop_matrix_sum(m);
    gsl_matrix_free(m);
    return out;
}

/** Like \c apop_matrix_map, but returns the sum of the resulting mapped vector. For example, let \c log_like be a function that returns the log likelihood of an input vector; then <tt>apop_matrix_map_sum(m, log_like)</tt> returns the total log likelihood of the rows of \c m.

  \li If you input a \c NULL matrix, I return the sum of zero items: zero.
  \li See \ref mapply "the map/apply page" for details.  */
double apop_matrix_map_sum(const gsl_matrix *in, double (*fn)(gsl_vector*)){
    if (!in) return 0;
    gsl_vector *v = apop_matrix_map (in, fn);
    double out = apop_vector_sum(v);
    gsl_vector_free(v);
    return out;
}

/** A function that effectively calls \ref apop_map and returns the sum of the resulting elements. Thus, this function returns a single \c double. See the \ref apop_map page for details of the inputs, which are the same here, except that \c inplace doesn't make sense---this function will always just add up the input function outputs.

  See also the \ref mapply "map/apply page" for details.

\li I don't copy the input data to send to your input function. Therefore, if your function modifies its inputs as a side-effect, your data set will be modified as this function runs.

\li The sum of zero elements is zero, so that is what is returned if the input \ref apop_data set is \c NULL. If <tt>apop_opts.verbose >= 2</tt> I print a warning.
 \ingroup mapply
 */
#ifdef APOP_NO_VARIADIC
double apop_map_sum(apop_data *in, apop_fn_d *fn_d, apop_fn_v *fn_v, apop_fn_r *fn_r, apop_fn_dp *fn_dp, apop_fn_vp *fn_vp, apop_fn_rp *fn_rp, apop_fn_dpi *fn_dpi,  apop_fn_vpi *fn_vpi, apop_fn_rpi *fn_rpi, apop_fn_di *fn_di, apop_fn_vi *fn_vi, apop_fn_ri *fn_ri, void *param, char part, int all_pages){
#else
apop_varad_head(double, apop_map_sum){ 
    apop_data * apop_varad_var(in, NULL)
    Apop_stopif(!in, return 0, 2, "NULL input. Returning zero.");
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
    char apop_varad_var(part, ((fn_v||fn_vp||fn_vpi||fn_vi) ? 'r' : 'a'));
    int apop_varad_var(all_pages, 'n')
    return apop_map_sum_base(in, fn_d, fn_v, fn_r, fn_dp, fn_vp, fn_rp, fn_dpi, fn_vpi, fn_rpi, fn_di, fn_vi, fn_ri, param, part, all_pages);
}

 double apop_map_sum_base(apop_data *in, apop_fn_d *fn_d, apop_fn_v *fn_v, apop_fn_r *fn_r, apop_fn_dp *fn_dp, apop_fn_vp *fn_vp, apop_fn_rp *fn_rp, apop_fn_dpi *fn_dpi,  apop_fn_vpi *fn_vpi, apop_fn_rpi *fn_rpi, apop_fn_di *fn_di, apop_fn_vi *fn_vi, apop_fn_ri *fn_ri, void *param, char part, int all_pages){
#endif 
    apop_data *mapped = apop_map(in, .fn_d=fn_d, .fn_v=fn_v, .fn_r=fn_r, 
                        .fn_dp=fn_dp, .fn_vp=fn_vp, .fn_rp=fn_rp, 
                        .fn_dpi=fn_dpi,  .fn_vpi=fn_vpi, .fn_rpi=fn_rpi, 
                        .fn_di=fn_di, .fn_vi=fn_vi, .fn_ri=fn_ri, 
                        .param=param, .part=part, .inplace='n', .all_pages='n');
    double outsum =   (mapped->vector ? apop_sum(mapped->vector) : 0)
                    + (mapped->matrix ? apop_matrix_sum(mapped->matrix) : 0);
    apop_data_free(mapped);
    return outsum + 
                    (((all_pages=='y' || all_pages=='Y') && in->more) ? 
                        apop_map_sum_base(in->more, fn_d, fn_v, fn_r, fn_dp, 
                        fn_vp, fn_rp, fn_dpi, fn_vpi, fn_rpi, fn_di, fn_vi, 
                        fn_ri, param, part, all_pages) : 0);
}
/** \} */
