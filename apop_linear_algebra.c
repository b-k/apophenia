
/** \file apop_linear_algebra.c	Assorted things to do with matrices,
such as take determinants or do singular value decompositions.  Includes
many convenience functions that don't actually do math but add/delete
columns, check bounds, et cetera.
*/ 
/* Copyright (c) 2006--2007, 2012 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

/** \defgroup linear_algebra 	Singular value decompositions, determinants, et cetera.  

This page describes some standard bits of linear algebra that Apophenia facilitates.

See also the printing functions, \ref apop_print, and the
\ref convenience_fns "Convenience functions".
*/

/** \defgroup convenience_fns 	Things to make life easier with the GSL
 */

#include "apop_internal.h"

void apop_gsl_error(const char *reason, const char *file, int line, int gsl_errno){
    Apop_notify(1, "%s: %s", file, reason);
    Apop_maybe_abort(1);
}
#define Checkgsl(...) if (__VA_ARGS__) {goto done;}
#define Check_gsl_with_out(...) if (__VA_ARGS__) {out->error='m'; goto done;}
#define Check_gsl_with_outmp(...) if (__VA_ARGS__) {gsl_matrix_free(*out); *out=NULL; goto done;}
#define Set_gsl_handler gsl_error_handler_t *prior_handler = gsl_set_error_handler(apop_gsl_error);
#define Unset_gsl_handler gsl_set_error_handler(prior_handler);

/**
Calculate the determinant of a matrix, its inverse, or both, via LU decomposition. The \c in matrix is not destroyed in the process.

\see apop_matrix_determinant,  apop_matrix_inverse

\param in The matrix to be inverted/determined. 
\param out If you want an inverse, this is where to place the matrix to be filled with the inverse. Will be allocated by the function. 

\param calc_det 
0: Do not calculate the determinant.\\
1: Do.

\param calc_inv
0: Do not calculate the inverse.\\
1: Do.

\return If <tt>calc_det == 1</tt>, then return the determinant. Otherwise, just returns zero.  If <tt>calc_inv!=0</tt>, 
then \c *out is pointed to the matrix inverse. In case of difficulty, I will set <tt>*out=NULL</tt> and return \c NaN.

\ingroup linear_algebra
*/

double apop_det_and_inv(const gsl_matrix *in, gsl_matrix **out, int calc_det, int calc_inv) {
    Set_gsl_handler
    Apop_stopif(in->size1 != in->size2, *out=NULL; return GSL_NAN, 0, "You asked me to invert a %zu X %zu matrix, "
            "but inversion requires a square matrix.", in->size1, in->size2);
    int sign;
    double the_determinant = GSL_NAN;
	gsl_matrix *invert_me = gsl_matrix_alloc(in->size1, in->size1);
	gsl_permutation * perm = gsl_permutation_alloc(in->size1);
	gsl_matrix_memcpy (invert_me, in);
	Checkgsl(gsl_linalg_LU_decomp(invert_me, perm, &sign))
	if (calc_inv){
		*out = gsl_matrix_alloc(in->size1, in->size1); //square.
		Check_gsl_with_outmp(gsl_linalg_LU_invert(invert_me, perm, *out))
    }
	if (calc_det)
		the_determinant	= gsl_linalg_LU_det(invert_me, sign);
    done:
	gsl_matrix_free(invert_me);
	gsl_permutation_free(perm);
    Unset_gsl_handler
	return the_determinant;
}

/**
Inverts a matrix. The \c in matrix is not destroyed in the process.
You may want to call \ref apop_matrix_determinant first to check that your input is invertible, or use \ref apop_det_and_inv to do both at once.

\param in The matrix to be inverted.
\return Its inverse.
\ingroup linear_algebra
*/
gsl_matrix * apop_matrix_inverse(const gsl_matrix *in) {
    gsl_matrix *out;
    apop_det_and_inv(in, &out, 0, 1);
    return out;
}

/**
Find the determinant of a matrix. The \c in matrix is not destroyed in the process.

See also \ref apop_matrix_inverse ,  or \ref apop_det_and_inv to do both at once.

\param in The matrix to be determined.
\return     The determinant.
\ingroup linear_algebra
*/
double apop_matrix_determinant(const gsl_matrix *in) {
    return apop_det_and_inv(in, NULL, 1, 0);
}

/** Principal component analysis: hand in a matrix and (optionally) a number of desired dimensions, and I'll return a data set where each column of the matrix is an eigenvector. The columns are sorted, so column zero has the greatest weight. The vector element of the data set gives the weights.

You also specify the number of elements your principal component space should have. If this is equal to the rank of the space in which the input data lives, then the sum of weights will be one. If the dimensions desired is less than that (probably so you can prepare a plot), then the weights will be accordingly smaller, giving you an indication of how much variation these dimensions explain. 

\param data The input matrix. (No default. If \c NULL, return \c NULL and print a warning iff <tt>apop_opts.verbose >= 1</tt>.)
I modify int in place so that each column has mean zero.

\param dimensions_we_want  (default: the size of the covariance matrix, i.e. <tt>data->size2</tt>)
The singular value decomposition will return this many of the eigenvectors with the largest eigenvalues.

\return     Returns a \ref apop_data set whose matrix is the principal component space. Each column of the returned matrix will be another eigenvector; the columns will be ordered by the eigenvalues. 
The data set's vector will be the largest eigenvalues, scaled by the total of all eigenvalues (including those that were thrown out). The sum of these returned values will give you the percentage of variance explained by the factor analysis.
\exception out->error=='a'  Allocation error.
\ingroup linear_algebra */
#ifdef APOP_NO_VARIADIC
apop_data * apop_matrix_pca(gsl_matrix *data, int const dimensions_we_want){
#else
apop_varad_head(apop_data *, apop_matrix_pca) {
    gsl_matrix * apop_varad_var(data, NULL);
    Apop_stopif(!data, return NULL, 1, "NULL data input");
    int const apop_varad_var(dimensions_we_want, data->size2);
    return apop_matrix_pca_base(data, dimensions_we_want);
}

 apop_data * apop_matrix_pca_base(gsl_matrix *data, int const dimensions_we_want){
#endif
    Set_gsl_handler
    apop_data *pc_space	= apop_data_alloc(0, data->size2, dimensions_we_want);
    Apop_stopif(pc_space->error, return pc_space, 0, "Allocation error.");
	pc_space->vector = gsl_vector_alloc(dimensions_we_want);
    Apop_stopif(!pc_space->vector, pc_space->error='a'; return pc_space, 
                0, "Allocation error setting up a %i vector.", dimensions_we_want);
    gsl_matrix *eigenvectors = gsl_matrix_alloc(data->size2, data->size2);
    gsl_vector *dummy_v 	 = gsl_vector_alloc(data->size2);
    gsl_vector *all_evalues  = gsl_vector_alloc(data->size2);
    gsl_matrix *square  	 = gsl_matrix_calloc(data->size2, data->size2);
    Apop_stopif(!eigenvectors || !dummy_v || !all_evalues || !square, pc_space->error='a'; return pc_space, 
                0, "Allocation error setting up workspace for %zu dimensions.", data->size2);
    double eigentotals	= 0;
    apop_matrix_normalize(data, 'c', 'm');
	Checkgsl(gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, data, data, 0, square))
	Checkgsl(gsl_linalg_SV_decomp(square, eigenvectors, all_evalues, dummy_v))
	for (int i=0; i< all_evalues->size; i++)
		eigentotals	+= gsl_vector_get(all_evalues, i);
	for (int i=0; i<dimensions_we_want; i++){
		Apop_col_v(&(apop_data){.matrix=eigenvectors}, i, v);
		gsl_matrix_set_col(pc_space->matrix, i, v);
		gsl_vector_set(pc_space->vector, i, gsl_vector_get(all_evalues, i)/eigentotals);
	}
    done:
	gsl_vector_free(dummy_v); 	gsl_vector_free(all_evalues);
	gsl_matrix_free(square); 	gsl_matrix_free(eigenvectors);
    Unset_gsl_handler
    return pc_space;
}

/** Take the log (base ten) of every element in a vector.

\li If the input vector is \c NULL, do nothing. 
\ingroup convenience_fns
 */
void apop_vector_log10(gsl_vector *v){
    if (!v) return;
    for (size_t i=0; i< v->size; i++){
        double *d = gsl_vector_ptr(v, i);
	    *d = log10(*d);
    }
}

/** Take the natural log of every element in a vector.

\li If the input vector is \c NULL, do nothing. 
\ingroup convenience_fns
 */
void apop_vector_log(gsl_vector *v){
    if (!v) return;
    for (size_t i=0; i< v->size; i++){
        double *d  = gsl_vector_ptr(v, i);
	    *d = gsl_sf_log(*d);
    }
}

/** Replace every vector element \f$v_i\f$ with exp(\f$v_i\f$).

\li If the input vector is \c NULL, do nothing. 
\ingroup convenience_fns
 */
void apop_vector_exp(gsl_vector *v){
    if (!v) return;
    for (size_t i=0; i< v->size; i++){
        double *d = gsl_vector_ptr(v, i);
        *d = exp(*d);
    }
}

/** Put the first vector on top of the second vector.

\param  v1  the upper vector (default=\c NULL, in which case this basically copies \c v2)
\param  v2  the second vector (default=\c NULL, in which case nothing is added)
\param  inplace If 'y', use \ref apop_vector_realloc to modify \c v1 in place; see the caveats on that function. Otherwise, allocate a new vector, leaving \c v1 unmolested. (default='n')
\return     the stacked data, either in a new vector or a pointer to \c v1.

\li This function uses the \ref designated syntax for inputs.
\ingroup convenience_fns
*/
#ifdef APOP_NO_VARIADIC
gsl_vector * apop_vector_stack(gsl_vector *v1, gsl_vector * v2, char inplace){
#else
apop_varad_head(gsl_vector *, apop_vector_stack){
    gsl_vector * apop_varad_var(v1, NULL);
    gsl_vector * apop_varad_var(v2, NULL);
    char apop_varad_var(inplace, 'n');
    return apop_vector_stack_base(v1, v2, inplace);
}

 gsl_vector * apop_vector_stack_base(gsl_vector *v1, gsl_vector * v2, char inplace){
#endif
    gsl_vector *out;
    gsl_vector t;
    if (!v1  && v2){
        out = gsl_vector_alloc(v2->size);
        gsl_vector_memcpy(out, v2);
        return out;
    } else if (!v2  && v1){
        if (inplace == 'y')
            return v1;
        out = gsl_vector_alloc(v1->size);
        gsl_vector_memcpy(out, v1);
        return out;
    } else if (!v1 && !v2)
        return NULL;
    //else:
    size_t v1size = v1->size; //save in case of reallocing.
    if (inplace == 'y' )
        out = apop_vector_realloc(v1, v1->size+v2->size);
    else {
        out = gsl_vector_alloc(v1->size + v2->size);
        t   = gsl_vector_subvector(out, 0, v1size).vector;
        gsl_vector_memcpy(&t, v1);
    }
    t   = gsl_vector_subvector(out, v1size, v2->size).vector;
    gsl_vector_memcpy(&t, v2);
    return out;
}

/** Put the first matrix either on top of or to the right of the second matrix.
  The fn returns a new matrix, meaning that at the end of this function, until you gsl_matrix_free() the original matrices, you will be taking up twice as much memory. Plan accordingly.

\param  m1  the upper/rightmost matrix (default=\c NULL, in which case this basically copies \c m2)
\param  m2  the second matrix (default = \c NULL, in which case \c m1 is returned)
\param  posn    if 'r', stack rows on top of other rows, else, e.g. 'c' stack  columns next to columns. (default ='r')
\param  inplace If 'y', use \ref apop_matrix_realloc to modify \c m1 in place; see the caveats on that function. Otherwise, allocate a new matrix, leaving \c m1 unmolested. (default='n')
\return     the stacked data, either in a new matrix or a pointer to \c m1.

\ingroup convenience_fns

For example, here is a little function to merge four matrices into a single two-part-by-two-part matrix. The original matrices are unchanged.
\code
gsl_matrix *apop_stack_two_by_two(gsl_matrix *ul, gsl_matrix *ur, gsl_matrix *dl, gsl_matrix *dr){
  gsl_matrix *output, *t;
    output = apop_matrix_stack(ul, ur, 'c');
    t = apop_matrix_stack(dl, dr, 'c');
    apop_matrix_stack(output, t, 'r', .inplace='y');
    gsl_matrix_free(t);
    return output;
}
\endcode

\li This function uses the \ref designated syntax for inputs.
*/
#ifdef APOP_NO_VARIADIC
gsl_matrix * apop_matrix_stack(gsl_matrix *m1, gsl_matrix * m2, char posn, char inplace){
#else
apop_varad_head(gsl_matrix *, apop_matrix_stack){
    gsl_matrix *apop_varad_var(m1, NULL);
    gsl_matrix *apop_varad_var(m2, NULL);
    char apop_varad_var(posn, 'r');
    char apop_varad_var(inplace, 'n');
    return apop_matrix_stack_base(m1, m2, posn, inplace);
}

 gsl_matrix * apop_matrix_stack_base(gsl_matrix *m1, gsl_matrix * m2, char posn, char inplace){
#endif
    gsl_matrix      *out;
    gsl_vector_view tmp_vector;
    if (!m1 && m2){
        out = gsl_matrix_alloc(m2->size1, m2->size2);
        gsl_matrix_memcpy(out, m2);
        return out;
    } else if (!m2 && m1) {
        if (inplace =='y')
            return m1;
        out = gsl_matrix_alloc(m1->size1, m1->size2);
        gsl_matrix_memcpy(out, m1);
        return out;
    } else if (!m2  && !m1) 
        return NULL;

    if (posn == 'r'){
        Apop_stopif(m1->size2 != m2->size2, return NULL, 0, "When stacking matrices on top of each other, they have to have the same number of columns, but  m1->size2==%zu and m2->size2==%zu. Returning NULL.", m1->size2, m2->size2);
        int m1size = m1->size1;
        if (inplace =='y')
            out = apop_matrix_realloc(m1, m1->size1 + m2->size1, m1->size2);
        else {
            out     = gsl_matrix_alloc(m1->size1 + m2->size1, m1->size2);
            for (int i=0; i< m1size; i++){
                    tmp_vector  = gsl_matrix_row(m1, i);
                    gsl_matrix_set_row(out, i, &(tmp_vector.vector));
            }
        }
        for (int i=m1size; i< m1size + m2->size1; i++){
            tmp_vector  = gsl_matrix_row(m2, i- m1size);
            gsl_matrix_set_row(out, i, &(tmp_vector.vector));
        }
        return out;
    } else {
        Apop_stopif(m1->size1 != m2->size1, return NULL, 0, "When stacking matrices side by side, "
                "they have to have the same number of rows, but m1->size1==%zu and m2->size1==%zu. Returning NULL."
                , m1->size1, m2->size1);
        int m1size = m1->size2;
        if (inplace =='y')
            out = apop_matrix_realloc(m1, m1->size1, m1->size2 + m2->size2);
        else {
            out     = gsl_matrix_alloc(m1->size1, m1->size2 + m2->size2);
            for (int i=0; i< m1size; i++){
                tmp_vector  = gsl_matrix_column(m1, i);
                gsl_matrix_set_col(out, i, &(tmp_vector.vector));
            }
        }
        for (int i=0; i< m2->size2; i++){
            tmp_vector  = gsl_matrix_column(m2, i);
            gsl_matrix_set_col(out, i+ m1size, &(tmp_vector.vector));
        }
        return out;
    } 
}

/** Delete columns from a matrix. 

  This is done via copying, so if you have an exceptionally large
  data set, you're better off producing the matrix in the perfect form
  directly.

\param in   the \c gsl_matrix to be subsetted
\return     a \c gsl_matrix with the specified columns removed. If you ask me to remove no columns, I'll return a copy of the original. If you ask me to remove all columns, I'll return \c NULL.
\param drop an array of <tt>int</tt>s. If use[7]==1, then column seven will be cut from the output. 
*/
gsl_matrix *apop_matrix_rm_columns(gsl_matrix *in, int *drop){
    int ct  = 0,  //how many columns will not be dropped?
        j   = 0;
    for (size_t i=0; i < in->size2; i++)
        if (drop[i]==0)
            ct++;
    if (ct == in->size2)
        return apop_matrix_copy(in);
    if (ct == 0)
        return NULL;
    gsl_matrix *out = gsl_matrix_alloc(in->size1, ct);
    for (size_t i=0; i < in->size2; i++){
        if (drop[i]==0){
            Apop_col_v(&(apop_data){.matrix=in}, i, v);
            gsl_matrix_set_col(out, j, v);
            j   ++;
        }
    }
    return out;
}

/** Test for a situation when a vector is diverging,
so you can preempt a procedure that is about to break on infinite values.

Alternatively, set \c max to \c INFINITY (or \c GSL_INF) to just test whether all of the matrix's elements are finite.

 \param in  A <tt>gsl_vector</tt>
 \param max An upper and lower bound to the elements of the vector. (default: GSL_POSINF)
 \return    1 if everything is bounded: not Inf, -Inf, or NaN, and \f$-\max < x < \max\f$; zero otherwise. 
 
\li A \c NULL vector has no unbounded elements, so \c NULL input returns 1. You get a warning if <tt>apop_opts.verbosity >=1</tt>.
\li This function uses the \ref designated syntax for inputs.
\ingroup convenience_fns
 */
#ifdef APOP_NO_VARIADIC
int apop_vector_bounded(const gsl_vector *in, long double max){
#else
apop_varad_head(int, apop_vector_bounded){
    const gsl_vector * apop_varad_var(in, NULL)
    Apop_stopif(!in, return 1, 1, "You sent in a NULL vector; returning 1.");
    long double apop_varad_var(max, GSL_POSINF)
    return apop_vector_bounded_base(in, max);
}

 int apop_vector_bounded_base(const gsl_vector *in, long double max){
#endif
    double x;
    for (size_t i=0; i< in->size; i++){
        x   = gsl_vector_get(in, i);
        if (!gsl_finite(x) || x> max || x< -max)
            return 0;
    }
    return 1;
}


static gsl_vector* dot_for_apop_dot(const gsl_matrix *m, const gsl_vector *v, 
                             const CBLAS_TRANSPOSE_t flip){
    #define Check_gslv(...) if (__VA_ARGS__) {gsl_vector_free(out); out=NULL;}
    gsl_vector *out = (flip ==CblasNoTrans)
                        ? gsl_vector_calloc(m->size1)
                        : gsl_vector_calloc(m->size2);
    Check_gslv(gsl_blas_dgemv (flip, 1.0, m, v, 0.0, out))
    return out;
}

/** A convenience function for dot products, which requires less prep and typing than the <tt>gsl_cblas_dgexx</tt> functions.

Second, it makes some use of the semi-overloading of the \ref apop_data structure. \c d1 may be a vector or a matrix, and the same for \c d2, so this function can do vector dot matrix, matrix dot matrix, and so on. If \c d1 includes both a vector and a matrix, then later parameters will indicate which to use.

\li This function uses the \ref designated syntax for inputs.

\param d1 the left part of \f$ d1 \cdot d2\f$
\param d2 the right part of \f$ d1 \cdot d2\f$
\param form1 't' or 'p': transpose or prime \c d1->matrix, or, if \c d1->matrix is \c NULL, read \c d1->vector as a row vector.<br>
                    'n' or 0: no transpose. (the default)<br>
                    'v': ignore the matrix and use the vector.

\param form2 As above, with \c d2.
\return     an \ref apop_data set. If two matrices come in, the vector element is \c NULL and the 
            matrix has the dot product; if either or both are vectors,
            the vector has the output and the matrix is \c NULL.

\exception out->error='a'  Allocation error.
\exception out->error='d'  dimension-matching error.
\exception out->error='m'  GSL math error.
\exception NULL If you ask me to take the dot product of NULL, I return NULL. [May some day change.]

\li Some systems auto-transpose non-conforming matrices. You input a \f$3 \times 5\f$ and
a \f$3 \times 5\f$ matrix, and the system assumes that you meant to transpose the second,
producing a \f$3 \times 5 \cdot 5 \times 3 \rightarrow 3 \times 3\f$ output. Apophenia
does not do this. First, it's ambiguous whether the output should be \f$3 \times 3\f$
or \f$5 \times 5\f$. Second, your next run might have three observations, and two \f$3 \times 3\f$ 
matrices don't require transposition; auto-transposition thus creates situations where
bugs can pop up on only some iterations of a loop.

\li For a vector \f$cdot\f$ a matrix, the vector is always treated as a row vector,
meaning that a \f$3\times 1\f$ dot a \f$3\times 4\f$ matrix is correct, and produces a
\f$1 \times 4\f$ vector.  For a vector \f$cdot\f$ a matrix, the vector is always treated
as a column vector. Requests for transposition are ignored.  

\li As a corrollary to the above rule, a vector dot a vector always produces a scalar,
 which will be put in the zeroth element of the output vector;
see the example. 

\li If you want to multiply an \f$N \times 1\f$ vector \f$\cdot\f$ a \f$1 \times N\f$
matrix produce an \f$N \times N\f$ matrix, then use \ref apop_vector_to_matrix to turn
your vectors into matrices; see the example.


\li A note for readers of <em>Modeling with Data</em>: the awkward instructions on using
this function on p 130 are now obsolete, thanks to the designated initializer syntax
for function calls. Notably, in the case where <tt>d1</tt> is a vector and <tt>d2</tt>
a matrix, then <tt>apop_dot(d1,d2,'t')</tt> won't work, because <tt>'t'</tt> now refers
to <tt>d1</tt>. Instead use <tt>apop_dot(d1,d2,.form2='t')</tt> or  <tt>apop_dot(d1,d2,0,
't')</tt>

Sample code:
\include dot_products.c

\ingroup linear_algebra
  */
#ifdef APOP_NO_VARIADIC
apop_data * apop_dot(const apop_data *d1, const apop_data *d2, char form1, char form2){
#else
apop_varad_head(apop_data *, apop_dot){
    const apop_data * apop_varad_var(d1, NULL)
    const apop_data * apop_varad_var(d2, NULL)
    Apop_stopif(!d1, return NULL, 1, "d1 is NULL; returning NULL");
    Apop_stopif(!d2, return NULL, 1, "d2 is NULL; returning NULL");
    char apop_varad_var(form1, 0)
    char apop_varad_var(form2, 0)
    return apop_dot_base(d1, d2, form1, form2);
}

 apop_data * apop_dot_base(const apop_data *d1, const apop_data *d2, char form1, char form2){
#endif
    Set_gsl_handler
    int         uselm, userm;
    gsl_matrix  *lm = d1->matrix, 
                *rm = d2->matrix;
    gsl_vector  *lv = d1->vector, 
                *rv = d2->vector;

    if (d1->matrix && form1 != 'v') uselm = 1;
    else if (d1->vector)            uselm = 0;
    else {
        Apop_stopif(form1 == 'v', return NULL, 0,
                    "You asked for a vector from the left data set, but "
                    "its vector==NULL. Returning NULL.");
        Apop_stopif(1, return NULL, 0, "The left data set has neither non-NULL "
                                  "matrix nor vector. Returning NULL.");
    }
    if (d2->matrix && form2 != 'v') userm = 1;
    else if (d2->vector)            userm = 0;
    else {
        Apop_stopif(form2 == 'v', return NULL, 0, 
                    "You asked for a vector from the right data set, but "
                    "its vector==NULL. Returning NULL.");
        Apop_stopif(1, return NULL, 0, "The right data set has neither non-NULL "
                                  "matrix nor vector. Returning NULL.");
    }
    apop_data *out = apop_data_alloc();
    #define Dimcheck(lr, lc, rr, rc) Apop_stopif((lc)!=(rr), out->error='d'; goto done,\
        0, "mismatched dimensions: %zuX%zu dot %zuX%zu. %s", (lr), (lc), (rr), (rc),\
        ((lr)==(rr)) ? " Maybe transpose the first?" \
        : ((rc)==(lc)) ? " Maybe transpose the second?" : "");

    CBLAS_TRANSPOSE_t lt, rt;
    lt  = (form1 == 'p' || form1 == 't' || form1 == 1) 
            ? CblasTrans: CblasNoTrans;
    rt  = (form2 == 'p' || form2 == 't' || form2 == 1) 
            ? CblasTrans: CblasNoTrans;
    if (uselm && userm){
        Dimcheck((lt== CblasNoTrans) ? lm->size1:lm->size2,
                 (lt== CblasNoTrans) ? lm->size2:lm->size1,
                 (rt== CblasNoTrans) ? rm->size1:rm->size2,
                 (rt== CblasNoTrans) ? rm->size2:rm->size1)
        gsl_matrix *outm = gsl_matrix_calloc((lt== CblasTrans)? lm->size2: lm->size1, 
                                             (rt== CblasTrans)? rm->size1: rm->size2);
        Check_gsl_with_out(gsl_blas_dgemm (lt,rt, 1, lm, rm, 0, outm))
        out->matrix = outm;
    } else if (!uselm && userm){
        Dimcheck((size_t)1, lv->size,
                 (rt== CblasNoTrans) ? rm->size1:rm->size2,
                 (rt== CblasNoTrans) ? rm->size2:rm->size1)
        //dgemv is always matrix first, then vector, so reverse from vm to mv:
        // if output vector has dimension matrix->size2, send CblasTrans
        // if output vector has dimension matrix->size1, send CblasNoTrans
        out->vector = dot_for_apop_dot(rm, lv
                        , (rt == CblasNoTrans) ? CblasTrans : CblasNoTrans);
        Apop_stopif(!out->vector, out->error='m'; goto done, 0, "GSL-level math error");
    } else if (uselm && !userm){
        Dimcheck((lt== CblasNoTrans) ? lm->size1:lm->size2,
                 (lt== CblasNoTrans) ? lm->size2:lm->size1,
                  rv->size , (size_t)1)
        out->vector = dot_for_apop_dot(lm, rv , lt);
        Apop_stopif(!out->vector, out->error='m'; goto done, 0, "GSL-level math error");
    } else if (!uselm && !userm){ 
        double outd;
        Check_gsl_with_out(gsl_blas_ddot(lv, rv, &outd))
        out->vector = gsl_vector_alloc(1);
        gsl_vector_set(out->vector, 0, outd);
    }

    //If using the vector, there's no meaningful name to assign.
    if (d1->names && uselm){
        if (lt == CblasTrans) apop_name_stack(out->names, d1->names, 'r', 'c');
        else                  apop_name_stack(out->names, d1->names, 'r');
    }
    if (d2->names && userm){
        if (rt == CblasTrans) apop_name_stack(out->names, d2->names, 'c', 'r');
        else                  apop_name_stack(out->names, d2->names, 'c');
    }

done:
    Unset_gsl_handler
    return out;
}
