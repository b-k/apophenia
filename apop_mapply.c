#include <apophenia/stats.h>
#include <apophenia/mapply.h>
#include <apophenia/types.h>

/** \file apop_mapply.c vector/matrix map/apply. 
 
Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */


//Support for apop_matrix_apply:

typedef double apop_fn_vtod(gsl_vector*);
typedef void apop_fn_vtov(gsl_vector*);
typedef double apop_fn_dtod(double);
typedef void apop_fn_dtov(double*);

typedef struct {
    size_t      *limlist;
    void        *fn;
    gsl_matrix  *m;
    gsl_vector  *vin;
    gsl_vector  *v;
} threadpass;

static void *forloop(void *t){
  threadpass    *tc = t;
  apop_fn_vtod  *vtod=tc->fn;
  apop_fn_vtov  *vtov=tc->fn;
  int           i;
  gsl_vector    view;
  double        val;
    for (i= tc->limlist[0]; i< tc->limlist[1]; i++){
        view    = gsl_matrix_row(tc->m, i).vector;
        if (tc->v){
            val     = vtod(&view);
            gsl_vector_set(tc->v, i, val);
        } else
            vtov(&view);
    }
  return NULL;
}

static void *vectorloop(void *t){
  threadpass    *tc = t;
  int           i;
  double        *inval, outval;
  apop_fn_dtod  *dtod=tc->fn;
  apop_fn_dtov  *dtov=tc->fn;
    for (i= tc->limlist[0]; i< tc->limlist[1]; i++){
        inval   = gsl_vector_ptr(tc->vin, i);
        if (tc->v){
            outval  = dtod(*inval);
            gsl_vector_set(tc->v, i, outval);
        } else
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

static gsl_vector*mapply_core(gsl_matrix *m,gsl_vector *vin, void *fn, gsl_vector *vout){
  int           threadct    = apop_opts.thread_count;
  pthread_t     thread_id[threadct];
  int           i;
  threadpass    tp[threadct];
    for (i=0 ; i<threadct; i++){
        tp[i].limlist   = threadminmax(i, m? m->size1: vin->size ,threadct);
        tp[i].fn        = fn;
        tp[i].m         = m;
        tp[i].vin       = vin;
        tp[i].v         = vout;
    }
    if (threadct==1){ //don't thread.
        if (m)  forloop(tp);
        else    vectorloop(tp);
    } else {
        for (i=0 ; i<threadct; i++){
            if (m)  pthread_create(&thread_id[i], NULL,forloop,(tp+i));
            else    pthread_create(&thread_id[i], NULL,vectorloop,(tp+i));
        }
        for (i=0 ; i<threadct; i++)
            pthread_join(thread_id[i], NULL);
    }
    for (i=0 ; i<threadct; i++)
        free(tp[i].limlist);
    return vout;
}


/** Map a function onto every row of a matrix.  The function that you
 input takes in a gsl_vector and returns a \c double. \c apop_matrix_map will
 produce a vector view of each row, and send each row to your function. It
 will output a \c gsl_vector holding your function's output for each row.

If \c apop_opts.thread_count is greater than one, then the matrix will be
broken into chunks and each sent to a different thread. Notice that the
GSL is generally threadsafe, and SQLite is 100\% not threadsafe. If
your function calls SQLite behind your back, you will find out when your
program crashes.

  \param m  The matrix
  \param fn A function of the form <tt>double fn(gsl_vector* in)</tt>

  \return A \c gsl_vector with the corresponding value for each row.

  See also \ref apop_matrix_apply, which works like this function but does not return a value.

  \ingroup convenience_fns
 */
gsl_vector *apop_matrix_map(gsl_matrix *m, double (*fn)(gsl_vector*)){
  gsl_vector    *out        = gsl_vector_alloc(m->size1);
    return mapply_core(m, NULL, fn, out);
}

/** Apply a function to every row of a matrix.  The function that you
 input takes in a gsl_vector and returns \c NULL. \c apop_matrix_apply will
 produce a vector view of each row, and send each row to your function.

If \c apop_opts.thread_count is greater than one, then the matrix will be
broken into chunks and each sent to a different thread. Notice that the
GSL is generally threadsafe, and SQLite is 100\% not threadsafe. If
your function calls SQLite behind your back, you will find out when your
program crashes.

  \param m  The matrix
  \param fn A function of the form <tt>void fn(gsl_vector* in)</tt>

  See also \ref apop_matrix_map, which works like this function but returns a value for each row.

  \ingroup convenience_fns
 */
void apop_matrix_apply(gsl_matrix *m, void (*fn)(gsl_vector*)){
    mapply_core(m, NULL, fn, NULL);
}

/** Map a function onto every element of a vector.  The function that you
 input takes in a \c double and returns a \c double. \c apop_apply will
 send each element to your function, and
 will output a \c gsl_vector holding your function's output for each row.

If \c apop_opts.thread_count is greater than one, then the matrix will be
broken into chunks and each sent to a different thread. Notice that the
GSL is generally threadsafe, and SQLite is 100\% not threadsafe. If
your function calls SQLite behind your back, you will find out when your
program crashes.

  \param v  The input vector
  \param fn A function of the form <tt>double fn(double in)</tt>

  \return A \c gsl_vector (allocated by this function) with the corresponding value for each row.

  See also \ref apop_vector_apply, which works like this function but does not return a value.

  \ingroup convenience_fns
 */
gsl_vector *apop_vector_map(gsl_vector *v, double (*fn)(double)){
  gsl_vector    *out        = gsl_vector_alloc(v->size);
    return mapply_core(NULL, v, fn, out);
}

/** Apply a function to every row of a matrix.  The function that you
 input takes in a gsl_vector and returns \c NULL. \c apop_apply will
 send a pointer to each element of your vector to your function.

If \c apop_opts.thread_count is greater than one, then the matrix will be
broken into chunks and each sent to a different thread. Notice that the
GSL is generally threadsafe, and SQLite is 100\% not threadsafe. If
your function calls SQLite behind your back, you will find out when your
program crashes.

  \param v  The input vector
  \param fn A function of the form <tt>void fn(double in)</tt>

  See also \ref apop_vector_map, which works like this function but returns a value.

  \ingroup convenience_fns
 */
void apop_vector_apply(gsl_vector *v, void (*fn)(double*)){
    mapply_core(NULL, v, fn, NULL);
}




void apop_matrix_map_all_vector_subfn(const gsl_vector *in, gsl_vector *outv, double (*fn)(double)){
    mapply_core(NULL, in, fn, outv);
}

void  apop_matrix_apply_all_vector_subfn(const gsl_vector *in, void (*fn)(double *)){
    apop_vector_apply(in, fn);
}


gsl_matrix * apop_matrix_map_all(const gsl_matrix *in, double (*fn)(double)){
  int i;
  gsl_matrix *out = gsl_matrix_alloc(in->size1, in->size2);
  gsl_vector_view v, inv;
    for (i=0; i< in->size1; i++){
        inv = gsl_matrix_row(in, i);
        v = gsl_matrix_row(out, i);
        apop_matrix_map_all_vector_subfn(&inv.vector, &v.vector, fn);
    }
    return out;
}

void apop_matrix_apply_all(gsl_matrix *in, void (*fn)(double *)){
  int   i;
  gsl_vector_view v;
    for (i=0; i< in->size1; i++){
        v   = gsl_matrix_row(in, i);
        apop_vector_apply(&v.vector, fn);
    }
}

/** like \c apop_vector_map, but returns the sum of the resulting mapped
 function. For example, <tt>apop_vector_map_sum(v, isnan)</tt> returns the number of elements of <tt>v</tt> that are \c NaN.

 Calls \c apop_vector_map internally, meaning that you will need adequate free memory.
 */
double apop_vector_map_sum(const gsl_vector *in, double(*fn)(double)){
    gsl_vector *m = apop_vector_map (in, fn);
    double out = apop_vector_sum(m);
    gsl_vector_free(m);
    return out;
}

/** like \c apop_matrix_map_all, but returns the sum of the resulting mapped
 function. For example, <tt>apop_matrix_map_all_sum(v, isnan)</tt> returns the number of elements of <tt>m</tt> that are \c NaN.

 Calls \c apop_matrix_map_all, meaning that you will need adequate free memory.
 */
double apop_matrix_map_all_sum(const gsl_matrix *in, double (*fn)(double)){
    gsl_matrix *m = apop_matrix_map_all (in, fn);
    double out = apop_matrix_sum(m);
    gsl_matrix_free(m);
    return out;
}

/** Like \c apop_matrix_map, but returns the sum of the resulting mapped
 vector. For example, let \c log_like be a function that returns the
 log likelihood of an input vector; then <tt>apop_matrix_map_sum(m,
 log_like)</tt> returns the total log likelihood of the rows of \c m.

 Calls \c apop_matrix_map, meaning that you will need adequate free memory.
 */
double apop_matrix_map_sum(const gsl_matrix *in, double (*fn)(gsl_vector*)){
    gsl_vector *v = apop_matrix_map (in, fn);
    double out = apop_vector_sum(v);
    gsl_vector_free(v);
    return out;
}
