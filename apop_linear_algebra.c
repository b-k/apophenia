/** \file apop_linear_algebra.c	Assorted things to do with matrices,
such as take determinants or do singular value decompositions.  Includes
many convenience functions that don't actually do math but add/delete
columns, check bounds, et cetera.

Copyright (c) 2006 by Ben Klemens. Licensed under the GNU GPL v2.
*/

/** \defgroup linear_algebra 	Singular value decompositions, determinants, et cetera.  

This page describes some standard bits of linear algebra that Apophenia facilitates.

See also the printing functions, \ref apop_print, and the
\ref convenience_fns "Convenience functions".
*/

/** \defgroup convenience_fns 	Things to make life easier with the GSL.
 */

/** \defgroup output		Printing to the screen or a text file

Most functions print only to the screen, but the \ref apop_print "matrix
and vector printing functions" will let you print to a text file as
well. The presumption is that statistic estimates are for your own
consumption, while you are printing a matrix for import into another program.

See \ref apop_name_print.
*/
/** \defgroup apop_print 	Asst printing functions		

The <tt>apop_*_print</tt> functions will print to screen, text file,
or database, depending on how you set \ref apop_opts.output_type.
The <tt>apop_*_show</tt> functions print only to screen, and are basically
just a convenience shell to the corresponding <tt>apop_*_print</tt>
function.

\ingroup output
*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>	//popen, I think.
#include <apophenia/linear_algebra.h> 
#include <apophenia/stats.h>
#include "math.h" //pow!
#include <apophenia/vasprintf.h>


/**
Calculate the determinant of a matrix, its inverse, or both. The \c in matrix is not destroyed in the process.

\param in
The matrix to be inverted/determined.

\param out 
If you want an inverse, this is where to place the matrix to be filled with the inverse. If <tt>calc_inv == 0</tt>, then use <tt>NULL</tt>. 

\param calc_det 
0: Do not calculate the determinant.

1: Do.

\param calc_inv
0: Do not calculate the inverse.

1: Do.

\return
If <tt>calc_det == 1</tt>, then return the determinant. Otherwise, just returns zero.
\ingroup linear_algebra
*/
double apop_det_and_inv(gsl_matrix *in, gsl_matrix **out, int calc_det, int calc_inv) {
  int 		sign;
  double 	the_determinant = 0;
	gsl_matrix *invert_me = gsl_matrix_alloc(in->size1, in->size1);
	gsl_permutation * perm = gsl_permutation_alloc(in->size1);
	invert_me = gsl_matrix_alloc(in->size1, in->size1);
	gsl_matrix_memcpy (invert_me, in);
	gsl_linalg_LU_decomp(invert_me, perm, &sign);
	if (calc_inv){
		*out	= gsl_matrix_alloc(in->size1, in->size1); //square.
		gsl_linalg_LU_invert(invert_me, perm, *out);
		}
	if (calc_det)
		the_determinant	= gsl_linalg_LU_det(invert_me, sign);
	gsl_matrix_free(invert_me);
	gsl_permutation_free(perm);
	return(the_determinant);
}

/** This comes up often enough that it deserves its own convenience function.
\param x	A vector.
\param sigma 	A symmetric matrix.
\return 	The scalar X'SX
\ingroup linear_algebra
*/
double apop_x_prime_sigma_x(gsl_vector *x, gsl_matrix *sigma){
  gsl_vector * 	sigma_dot_x	= gsl_vector_calloc(x->size);
  double		the_result;
	gsl_blas_dsymv(CblasUpper, 1, sigma, x, 0, sigma_dot_x); //sigma should be symmetric
	gsl_blas_ddot(x, sigma_dot_x, &the_result);
	gsl_vector_free(sigma_dot_x);
	return(the_result);
}

void apop_normalize_for_svd(gsl_matrix *in){
//Greene (2nd ed, p 271) recommends pre- and post-multiplying by sqrt(diag(X'X)) so that X'X = I.
  gsl_vector_view	v;
  gsl_vector	*diagonal = gsl_vector_alloc(in->size1);
  int 		i;
	//Get the diagonal, take the square root
	v	= gsl_matrix_diagonal(in);
	gsl_vector_memcpy(diagonal, &(v.vector));
	for (i=0; i<diagonal->size; i++)
		gsl_vector_set(diagonal, i, pow(gsl_vector_get(diagonal,i), .5));
	//mulitply each row and column by the diagonal vector.
	for (i=0; i<diagonal->size; i++){
		v	= gsl_matrix_column(in, i);
		gsl_vector_mul(&(v.vector), diagonal);
		v	= gsl_matrix_row(in, i);
		gsl_vector_mul(&(v.vector), diagonal);
	}
	gsl_vector_free(diagonal);
}

/**
Singular value decomposition, aka principal component analysis, aka factor analysis.

\param data 
The input matrix.

\param dimensions_we_want 
The singular value decomposition will return this many of the eigenvectors with the largest eigenvalues.

\return     Returns a \ref apop_data set whose matrix is the principal component space. Each column of the returned matrix will be another eigenvector; the columns will be ordered by the eigenvalues. 
The data set's vector will be the largest eigenvalues, scaled by the total of all eigenvalues (including those that were thrown out). The sum of these returned values will give you the percentage of variance explained by the factor analysis.

\ingroup linear_algebra */
apop_data * apop_sv_decomposition(gsl_matrix *data, int dimensions_we_want) {
//Get X'X
  gsl_matrix * 	eigenvectors 	= gsl_matrix_alloc(data->size2, data->size2);
  gsl_vector * 	dummy_v 	    = gsl_vector_alloc(data->size2);
  gsl_vector * 	all_evalues 	= gsl_vector_alloc(data->size2);
  gsl_matrix * 	square  	    = gsl_matrix_calloc(data->size2, data->size2);
  int 		    i;
  double		eigentotals	= 0;
  apop_data    *pc_space	    = apop_data_alloc(0,data->size2, dimensions_we_want);
	pc_space->vector = gsl_vector_alloc(dimensions_we_want);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, data, data, 0, square);
	apop_normalize_for_svd(square);	
	gsl_linalg_SV_decomp(square, eigenvectors, all_evalues, dummy_v);
	for (i=0; i< all_evalues->size; i++)
		eigentotals	+= gsl_vector_get(all_evalues, i);
	for (i=0; i<dimensions_we_want; i++){
		APOP_MATRIX_COL(eigenvectors, i, v);
		gsl_matrix_set_col(pc_space->matrix, i, v);
		gsl_vector_set(pc_space->vector, i, gsl_vector_get(all_evalues, i)/eigentotals);
	}
	gsl_vector_free(dummy_v); 	gsl_vector_free(all_evalues);
	gsl_matrix_free(square); 	gsl_matrix_free(eigenvectors);
    return pc_space;
}


/** Just add <tt>amt</tt> to a \c gsl_vector element. Equivalent to <tt>gsl_vector_set(gsl_vector_get(v, i) + amt, i)</tt>, but more readable (and potentially faster).

\param v The \c gsl_vector in question
\param i The location in the vector to be incremented.
\param amt The amount by which to increment. Of course, one can decrement by specifying a negative amt.
\ingroup convenience_fns
 */
inline void apop_vector_increment(gsl_vector * v, int i, double amt){
	v->data[i * v->stride]	+= amt;
}

/** Just add <tt>amt</tt> to a \c gsl_matrix element. Equivalent to <tt>gsl_matrix_set(gsl_matrix_get(m, i, j) + amt, i, j)</tt>, but more readable (and potentially faster).

\param m The \c gsl_matrix in question
\param i The row of the element to be incremented.
\param j The column of the element to be incremented.
\param amt The amount by which to increment. Of course, one can decrement by specifying a negative amt.
\ingroup convenience_fns
 */
inline void apop_matrix_increment(gsl_matrix * m, int i, int j, double amt){
	m->data[i * m->tda +j]	+= amt;
}


/** Take the log of every element in a vector.
\ingroup convenience_fns
 */
void apop_vector_log(gsl_vector *v){
  int     i;
  double  d;
    for (i=0; i< v->size; i++){
	    d   = v->data[i * v->stride];
	    v->data[i * v->stride]  = gsl_sf_log(d);
    }
}

/** Take the exponent of every element in a vector.
\ingroup convenience_fns
 */
void apop_vector_exp(gsl_vector *v){
  int     i;
  double  d;
    for (i=0; i< v->size; i++){
	    d   = v->data[i * v->stride];
	    v->data[i * v->stride]  = gsl_sf_exp(d);
    }
}

/** Put the first vector on top of the second vector.
  The fn returns a new vector, meaning that at the end of this function, until you gsl_vector_free() the original vectors, you will be taking up twice as much memory. Plan accordingly.

\param  v1  the upper vector
\param  v2  the second vector
\return     a new vector with the stacked data.

\ingroup convenience_fns
*/
gsl_vector *apop_vector_stack(gsl_vector *v1, gsl_vector * v2){
  gsl_vector      *out;
  gsl_vector      t;
    if (!v1  && v2){
        out = gsl_vector_alloc(v2->size);
        gsl_vector_memcpy(out, v2);
        return out;
    } else if (!v2  && v1){
        out = gsl_vector_alloc(v1->size);
        gsl_vector_memcpy(out, v1);
        return out;
    } else if (!v2 && !v2)
        return NULL;
    //else:
    out = gsl_vector_alloc(v1->size + v2->size);
    t   = gsl_vector_subvector(out, 0, v1->size).vector;
    gsl_vector_memcpy(&t, v1);
    t   = gsl_vector_subvector(out, v1->size, v2->size).vector;
    gsl_vector_memcpy(&t, v2);
    return out;
}

/** Put the first matrix either on top of or to the right of the second matrix.
  The fn returns a new matrix, meaning that at the end of this function, until you gsl_matrix_free() the original matrices, you will be taking up twice as much memory. Plan accordingly.

\param  m1  the upper/rightmost matrix
\param  m2  the second matrix
\param  posn    if 'r', stack rows on top of other rows, else, e.g. 'c' stack  columns next to columns.
\return     a new matrix with the stacked data.
\ingroup convenience_fns

For example, here is a little function to merge four matrices into a single two-part-by-two-part matrix:
\code
gsl_matrix *apop_stack_two_by_two(gsl_matrix *ul, gsl_matrix *ur, gsl_matrix *dl, gsl_matrix *dr){
gsl_matrix  *t1, *t2, *output;
    t1   = apop_matrix_stack(ul, ur, 'c');
    gsl_matrix_free(ul);
    gsl_matrix_free(ur);
    t2   = apop_matrix_stack(dl, dr, 'c');
    gsl_matrix_free(dl);
    gsl_matrix_free(dr);
    output  = apop_matrix_stack(t1, t2, 'r');
    gsl_matrix_free(t1);
    gsl_matrix_free(t2);
    return output;
}
\endcode
*/
gsl_matrix *apop_matrix_stack(gsl_matrix *m1, gsl_matrix * m2, char posn){
gsl_matrix      *out;
gsl_vector_view tmp_vector;
int             i;

    if (!m1 && m2){
        out = gsl_matrix_alloc(m2->size1, m2->size2);
        gsl_matrix_memcpy(out, m2);
        return out;
    } else if (!m2 && m1) {
        out = gsl_matrix_alloc(m1->size1, m1->size2);
        gsl_matrix_memcpy(out, m1);
        return out;
    } else if (!m2  && !m1) 
        return NULL;

    if (posn == 'r'){
        if (m1->size2 != m2->size2){
            printf("When stacking matrices on top of each other, they have to have the same number of columns (m1->size2==m2->size2). Returning NULL.\n");
            return NULL;
        }
        out         = gsl_matrix_alloc(m1->size1 + m2->size1, m1->size2);
        for (i=0; i< m1->size1; i++){
            tmp_vector  = gsl_matrix_row(m1, i);
            gsl_matrix_set_row(out, i, &(tmp_vector.vector));
        }
        for ( ; i< m1->size1+m2->size1; i++){   //i is not reinitialized.
            tmp_vector  = gsl_matrix_row(m2, i- m1->size1);
            gsl_matrix_set_row(out, i, &(tmp_vector.vector));
        }
        return out;
    } else {
        if (m1->size1 != m2->size1){
            printf("When stacking matrices side by side, they have to have the same number of rows (m1->size1==m2->size1). Returning NULL.\n");
            return NULL;
        }
        out         = gsl_matrix_alloc(m1->size1, m1->size2 + m2->size2);
        for (i=0; i< m1->size2; i++){
            tmp_vector  = gsl_matrix_column(m1, i);
            gsl_matrix_set_col(out, i, &(tmp_vector.vector));
        }
        for ( ; i< m1->size2+m2->size2; i++){   //i is not reinitialized.
            tmp_vector  = gsl_matrix_column(m2, i- m1->size2);
            gsl_matrix_set_col(out, i, &(tmp_vector.vector));
        }
        return out;
    } 
}

/** Delete columns from a matrix. 

  This is done via copying, so if you have an exceptionally large
  data set, you're better off producing the matrix in the perfect form
  directly.

\param in   the \c gsl_matrix to be subsetted
\return     a \c gsl_matrix with the specified columns removed.
\param drop an array of ints. If use[7]==1, then column seven will be cut from the output. 
*/
gsl_matrix *apop_matrix_rm_columns(gsl_matrix *in, int *drop){
gsl_matrix      *out;
int             i, 
                ct  = 0, 
                j   = 0;
gsl_vector_view v;
    for (i=0; i < in->size2; i++)
        if (drop[i]==0)
            ct++;
    out = gsl_matrix_alloc(in->size1, ct);
    for (i=0; i < in->size2; i++){
        if (drop[i]==0){
            v   = gsl_matrix_column(in, i);
            gsl_matrix_set_col(out, j, &(v.vector));
            j   ++;
        }
    }
    return out;
}



/** test whether any of the elements of a <tt>gsl_vector</tt> are missing.

 \param in  A <tt>gsl_vector</tt>
 \return    1 if any element is NaN; zero if all elements are numbers.
 \ingroup convenience_fns
 */
int apop_vector_isnan(gsl_vector *in){
size_t  i;
    for (i=0; i< in->size; i++)
        if (gsl_isnan(gsl_vector_get(in, i)))
            return 1;
    return 0;
}

/** test whether all of the elements of a <tt>gsl_vector</tt> are
 * finite.

 \param in  A <tt>gsl_vector</tt>
 \return    1 if everything is finite (not Inf, -Inf, or NaN); zero otherwise.
 \ingroup convenience_fns
 */
int apop_vector_finite(gsl_vector *in){
size_t  i;
    for (i=0; i< in->size; i++)
        if (!gsl_finite(gsl_vector_get(in, i)))
            return 0;
    return 1;
}

/** This is a variant of \ref apop_vector_finite, qv. This lets you
test for a situation when a vector is diverging, but not actually Inf,
so you can preempt a procedure that is about to break on infinite values.

 \param in  A <tt>gsl_vector</tt>
 \param max An upper and lower bound to the elements of the vector.
 \return    1 if everything is bounded: not Inf, -Inf, or NaN, and \f$-\max < x < \max\f$; zero otherwise.
 \ingroup convenience_fns
 */
int apop_vector_bounded(gsl_vector *in, long double max){
size_t      i;
long double x;
    for (i=0; i< in->size; i++){
        x   = gsl_vector_get(in, i);
        if (!gsl_finite(x) || x> max || x< -max)
            return 0;
    }
    return 1;
}


static apop_data *dot_for_apop_dot(gsl_matrix *m, gsl_vector *v, CBLAS_TRANSPOSE_t flip){
gsl_vector *outv;
    if (flip ==CblasNoTrans)
        outv = gsl_vector_calloc(m->size1);
    else
        outv = gsl_vector_calloc(m->size2);
    gsl_blas_dgemv (flip, 1.0, m, v, 0.0, outv);
    return apop_vector_to_data(outv);
}

/** A convenience function for dot products.

  First, this requires less typing than the <tt>gsl_cblas_dgexx</tt> functions.

  Second, it makes some use of the semi-overloading of the \ref apop_data
  structure. \c d1 may be a vector or a matrix, and the same for \c d2, so
  this function can do vector dot matrix, matrix dot matrix, and so on. If
  \c d1 includes both a vector and a matrix, then the matrix only is used.

\param d1 the left part of \f$ d1 \cdot d2\f$
\param d2 the right part of \f$ d1 \cdot d2\f$
\param t1 't' or 'p': transpose or prime d1.<br>
                    'n' or 0: no transpose.
\param t2 't' or 'p': transpose or prime d2.
                    'n' or 0: no transpose.
\return     an \ref apop_data set. If two matrices come in, the vector element is \c NULL and the 
            matrix has the dot product; if either or both are vectors,
            the vector has the output and the matrix is \c NULL
\ingroup linear_algebra
  */
//apop_data * apop_dot(apop_data *d1, apop_data *d2, char t1, char t2){
apop_data * apop_dot(apop_data *d1, apop_data *d2, ...){
int         uselm, userm;
gsl_matrix  *lm = d1->matrix, 
            *rm = d2->matrix;
gsl_vector  *lv = d1->vector, 
            *rv = d2->vector;
CBLAS_TRANSPOSE_t   lt  ,//= (t1=='t' || t1=='T' || t1=='p' || t1=='P') ? CblasTrans : CblasNoTrans,
                    rt  ;//= (t2=='t' || t2=='T' || t2=='p' || t2=='P') ? CblasTrans : CblasNoTrans;
apop_data   *out    = apop_data_alloc(0,0,0);
va_list		argp;
	va_start(argp, d2);
    if (d1->matrix)
        uselm   = 1;
    else if (d1->vector)
        uselm   = 0;
    else{
        printf("apop_dot: the left data set has neither non-NULL matrix nor vector. Returning NULL.\n");
        return NULL;
    }
    if (d2->matrix)
        userm   = 1;
    else if (d2->vector)
        userm   = 0;
    else{
        printf("apop_dot: the right data set has neither non-NULL matrix nor vector. Returning NULL.\n");
        return NULL;
    }
    if (uselm && userm){
        lt  = va_arg(argp, int) ? CblasTrans: CblasNoTrans;
        rt  = va_arg(argp, int) ? CblasTrans: CblasNoTrans;
        gsl_matrix *outm    = gsl_matrix_calloc((lt== CblasTrans)? lm->size2: lm->size1, 
                                                (rt== CblasTrans)? rm->size1: rm->size2);
        gsl_blas_dgemm (lt,rt, 1, lm, rm, 0, outm);
        out->matrix         = outm;
    } else if (!uselm && userm){
        //If output vector has dimension matrix->size2, send CblasTrans
        //If output vector has dimension matrix->size1, send CblasNoTrans
        rt  = va_arg(argp, int) ? CblasTrans: CblasNoTrans;
        if (rt == CblasNoTrans)
            out = dot_for_apop_dot(rm, lv, CblasTrans);
        else
            out = dot_for_apop_dot(rm, lv, CblasNoTrans);
    } else if (uselm && !userm){
        lt  = va_arg(argp, int) ? CblasTrans: CblasNoTrans;
        if (lt == CblasNoTrans)
            out = dot_for_apop_dot(lm, rv, CblasNoTrans);
        else
            out = dot_for_apop_dot(lm, rv, CblasTrans);
    } else if (!uselm && !userm){ 
        double outd;
        gsl_blas_ddot (lv, rv, &outd);
        out->vector = gsl_vector_alloc(1);
        gsl_vector_set(out->vector, 0, outd);
    }
	va_end(argp);
    return out;
}


/** Take a matrix and stack the second column under the first, the third
 under those, et cetera, producing a single vector. Reverse the process
 with \ref apop_vector_split_to_matrix.

 This is a relatively obscure function that you'll probably never
 use. Maybe you have in mind <tt>gsl_matrix_get_col</tt> to pull a single
 column from a matrix?


 \param m An \f$m\times n\f$ matrix to stack
 \return A single vector of length \f$m\times n\f$

 \ingroup conversions

 */
gsl_vector *apop_matrix_stack_to_vector(gsl_matrix *m){
int         i;
gsl_vector  *out    = gsl_vector_alloc(m->size1 * m->size2);
gsl_vector  view;
    for(i=0; i< m->size2; i++){
        view    = gsl_vector_subvector(out, i*m->size1, m->size1).vector;
        gsl_matrix_get_col(&view, m, i);
    }
    return out;
}

/** Put the first few elements of a vector in the first column of a
  matrix, the next few in the second column, et cetera, thus producing a matrix from a vector.
  Reverse the process with \ref apop_matrix_stack_to_vector.

 \param     v   A vector of length \f$m\f$.
 \param columns The number of columns to be output.
 \return        A matrix of dimension \f$\frac{m}{\rm cols} \times {\rm cols}\f$

 \ingroup conversions
*/
gsl_matrix *apop_vector_split_to_matrix(gsl_vector *v, int columns){
int         i;
gsl_vector  view;
gsl_matrix  *out    = gsl_matrix_alloc(v->size/columns, columns);
    for(i=0; i< columns; i++){
        view    = gsl_vector_subvector(v, i*v->size/columns, v->size/columns).vector;
        gsl_matrix_set_col(out, i, &view);
    }
    return out;
}

