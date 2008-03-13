/** \file apop_data.c	 The apop_data structure joins together a
  gsl_matrix, apop_name, and a table of strings. No biggie.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include <assert.h>
#include <gsl/gsl_matrix.h>
#include "types.h"
#include "output.h"
#include "linear_algebra.h"

/** \defgroup data_struct apop_data

  The \c apop_data structure represents a data set.  It joins together
  a gsl_vector, a gsl_matrix, an apop_name, and a table of strings. It
  tries to be minimally intrusive, so you can use it everywhere you
  would use a \c gsl_matrix or a \c gsl_vector.

  For example, let us say that you are running a regression: there is
  a vector for the dependent variable, and a matrix for the dependent
  variables. Think of them as a partitioned matrix, where the vector is column -1, and the first column of the matrix is column zero. Here is some code to print the entire matrix. Notice that the column counter \c i starts counting at -1.

  \code
  for (j = 0; j< data->matrix->size1; j++){
    printf("%s\t", apop_name_get(data->names, j, 'r'));
    for (i = -1; i< data->matrix->size2; i++)
        printf("%g\t", apop_data_get(data, j, i));
    printf("\n");
    }
    \endcode

We're generally assuming that the data vector and data matrix have the
same row count: \c data->vector->size==data->matrix->size1 . This means
that the \ref apop_name structure doesn't have separate vector_names
and row_names elements: the rownames are assumed to apply for both.

  \ingroup types
  */

/** Allocate a \ref apop_data structure, to be filled with data. The
  three arguments are: vector size, matrix rows, matrix cols. Any and
  all of these may be zero, in which case the vector or matrix will not be allocated, as appropriate.

    Your best bet for allocating the categories is to produce
them elsewhere, such as \ref apop_query_to_text and then point an
apop_data structure with a zero-sized vector to your matrix of strings.

The \c weights vector is set to \c NULL. If you need it, allocate it via
\code d->weights   = gsl_vector_alloc(row_ct); \endcode.

See also \ref apop_data_calloc.

  \param vsize              vector size, if any. 
  \param msize1, msize2     Row and column size for the matrix. If \c size2>0
  this exactly mirrors the format of \c gsl_matrix_alloc. If \c size2==-1,
  then allocate a vector. 
  
  \c apop_data_alloc(0,0,0) will produce a basically blank set, with \c out->matrix==out->vector==NULL. 

 \return    The \ref apop_data structure, allocated and ready.
 \ingroup data_struct
  */
apop_data * apop_data_alloc(const size_t vsize, const size_t msize1, const int msize2){
  apop_data  *setme       = malloc(sizeof(apop_data));
    setme->vector   = NULL;
    setme->matrix   = NULL;
    setme->weights  = NULL;
    if (msize2 > 0  && msize1 > 0)
        setme->matrix   = gsl_matrix_alloc(msize1,msize2);
    if (vsize)
        setme->vector   = gsl_vector_alloc(vsize);
    else if (msize2==-1 && msize1>0)
        setme->vector   = gsl_vector_alloc(msize1);

    setme->names        = apop_name_alloc();
    setme->text         = NULL;
    setme->textsize[0]  =
    setme->textsize[1]  = 0;
    return setme;
}


/** Allocate a \ref apop_data structure, to be filled with data; set
 everything in the allocated portion to zero. See \ref apop_data_alloc for details.

  \param vsize              vector size, if any. 
  \param msize1, msize2     Row and column size for the matrix. If \c size2>0
  this exactly follows the format of \c gsl_matrix_calloc. If \c size2==-1,
  then allocate a vector. 
  
  \c apop_data_calloc(0, 0,0) will produce a basically blank set, with \c out->matrix==out->vector==NULL. 

 \return    The \ref apop_data structure, allocated and zeroed out.
 \ingroup data_struct
  */
apop_data * apop_data_calloc(const size_t vsize, const size_t msize1, const int msize2){
  apop_data  *setme       = malloc(sizeof(apop_data));
    setme->vector   = NULL;
    setme->matrix   = NULL;
    setme->weights  = NULL;
    if (msize2 >0 && msize1 > 0)
        setme->matrix   = gsl_matrix_calloc(msize1,msize2);
    if (vsize)
        setme->vector   = gsl_vector_calloc(vsize);
    else if (msize2==-1 && msize1>0)
        setme->vector   = gsl_vector_calloc(msize1);

    setme->names        = apop_name_alloc();
    setme->text         = NULL;
    setme->textsize[0]  =
    setme->textsize[1]  = 0;
    return setme;
}


/** Wrap an \ref apop_name structure around an existing \c gsl_matrix.

\param m    The existing matrix you'd like to turn into an \ref apop_data structure.
return      The \ref apop_data structure in question.
  */
apop_data * apop_matrix_to_data(gsl_matrix *m){
  apop_assert_void(m, 1, 'c',"Converting a NULL matrix to an apop_data structure.");
  apop_data  *setme   = apop_data_alloc(0,0,0);
    setme->matrix = m;
    return setme;
}

/** Wrap an \ref apop_name structure around an existing \c gsl_vector.

\param  v   The data vector
\return     an allocated, ready-to-use \ref apop_data struture.
*/
apop_data * apop_vector_to_data(gsl_vector *v){
  apop_assert_void(v, 1, 'c',"Converting a NULL vector to an apop_data structure.");
  apop_data  *setme   = apop_data_alloc(0,0,0);
    setme->vector = v;
    return setme;
}

/** free a matrix of chars* (i.e., a char***).  */
void apop_text_free(char ***freeme, int rows, int cols){
  int  i, j;
    if (rows && cols)
        for (i=0; i < rows; i++){
            for (j=0; j < cols; j++)
                free(freeme[i][j]);
            free(freeme[i]);
        }
    free(freeme);
}


/** Free an \ref apop_data structure.
 
  As with \c free(), it is safe to send in a \c NULL pointer (in which case the funtion does nothing).

 \ingroup data_struct
  */
void apop_data_free(apop_data *freeme){
    if (!freeme) return;
    if (freeme->vector)
        gsl_vector_free(freeme->vector);
    if (freeme->matrix)
        gsl_matrix_free(freeme->matrix);
    if (freeme->weights)
        gsl_vector_free(freeme->weights);
    apop_name_free(freeme->names);
    apop_text_free(freeme->text, freeme->textsize[0] , freeme->textsize[1]);
    free(freeme);
}


/** Copy one \ref apop_data structure to another. That is, all data is duplicated.

  This function does <i>not</i> allocate the output for you. If you want such behavior, usr \ref apop_data_copy.
 
  \param out    a structure that this function will fill. Must be preallocated
  \param in    the input data

 \ingroup data_struct
 \todo This doesn't copy over the category data.
  */
void apop_data_memcpy(apop_data *out, const apop_data *in){
  size_t   i, j;
    if (!out)
        apop_error(1,'c',"%s: you are copying to a NULL vector. Do you mean to use apop_data_copy instead?\n", __func__);
    if (in->matrix){
        if (in->matrix->size1 != out->matrix->size1 ||
                in->matrix->size2 != out->matrix->size2){
            if (apop_opts.verbose)
                apop_error(1, 'c',"%s: You're trying to copy a (%zu X %zu) into a (%zu X %zu) matrix. Returning w/o any copying.\n", 
                __func__, in->matrix->size1, in->matrix->size2, 
                out->matrix->size1, out->matrix->size2);
            return;
        }
        gsl_matrix_memcpy(out->matrix, in->matrix);
    }
    if (in->vector){
        if (in->vector->size != out->vector->size){
            if (apop_opts.verbose)
                apop_error(1, 'c',"%s: You're trying to copy a %zu-elmt vector into a %zu-elmt vector. Returning w/o any copying.\n", 
                __func__, in->vector->size, out->vector->size);
            return;
        }
        gsl_vector_memcpy(out->vector, in->vector);
    }
    apop_name_stack(out->names, in->names, 'r');
    apop_name_stack(out->names, in->names, 'c');
    apop_name_stack(out->names, in->names, 't');
    out->textsize[0] = in->textsize[0];
    out->textsize[1] = in->textsize[1];
    if (in->textsize[0] && in->textsize[1]){
        out->text  = malloc(sizeof(char ***) * in->textsize[0] * in->textsize[1]);
        for (i=0; i< in->textsize[0]; i++){
            out->text[i]  = malloc(sizeof(char **) * in->textsize[1]);
            for(j=0; j < in->textsize[1]; j ++)
				asprintf(&(out->text[i][j]), in->text[i][j]);
        }
    }
}

/** Copy one \ref apop_data structure to another. That is, all data is duplicated.

  Just a front-end for \ref apop_data_memcpy for those who prefer this sort of syntax.
 
  \param in    the input data
  \return       a structure that this function will allocate and fill. If input is NULL, then this will be NULL.

 \ingroup data_struct
  */
apop_data *apop_data_copy(const apop_data *in){
    if (!in)
        return NULL;
    apop_data *out  = apop_data_alloc(0, 0, 0);
    if (in->vector)
        out->vector = gsl_vector_alloc(in->vector->size);
    if (in->matrix)
        out->matrix = gsl_matrix_alloc(in->matrix->size1, in->matrix->size2);
    apop_data_memcpy(out, in);
    return out;
}

/** Put the first data set either on top of or to the right of the second matrix.

The fn returns a new data set, meaning that at the end of this function,
until you apop_data_free() the original data sets, you will be taking up
twice as much memory. Plan accordingly. 

 For the opposite operation, see \ref apop_data_split.

\param  m1      the upper/rightmost data set
\param  m2      the second data set
\param  posn    if 'r', stack rows of m1's matrix above rows of m2's<br>
if 'c', stack columns of m1's matrix to right of m2's<br>

If m1 or m2 are NULL, this returns a copy of the other element, and if
both are NULL, you get NULL back.

\return         a new \ref apop_data set with the stacked data.
Categories get dropped. If stacking rows on rows, the output vector is the input
vectors stacked accordingly. If stacking columns by columns, the output
vector is just the vector of m1 and m2->vector doesn't appear in the
output at all.  [If you don't like this behavior, send me code showing
what you'd prefer.]

\ingroup data_struct
*/
apop_data *apop_data_stack(apop_data *m1, apop_data * m2, char posn){
  gsl_matrix  *stacked= NULL;
  apop_data   *out    = NULL;
    if (m1 == NULL)
        out = apop_data_copy(m2);
    else if (m2 == NULL)
        out = apop_data_copy(m1);
    else if (posn == 'r' || posn == 'c'){
        stacked = apop_matrix_stack(m1->matrix, m2->matrix, posn);
        out     = apop_matrix_to_data(stacked);
        if (posn == 'c')
            out->vector = apop_vector_stack(m1->vector, m2->vector);
        else 
            out->vector = m1->vector;
        out->names  = apop_name_copy(m1->names);
        apop_name_stack(out->names, m2->names, posn);
    } else
        apop_error(0, 'c',"%s: valid positions are 'r' or 'c' ('v' is deprecated); you gave me >%c<. Returning NULL.", __func__, posn);
    return out;
}

/** Split one input \ref apop_data structure into two.

 For the opposite operation, see \ref apop_data_stack.

 For this function, the \ref apop_data->vector is taken to be the -1st
 element of the matrix.
 \param in  The \ref apop_data structure to split
 \param splitpoint The index of what will be the first row/column of the
 second data set.  E.g., if this is -1 and \c r_or_c=='c', then the whole
 data set will be in the second data set; if this is the length of the
 matrix then the whole data set will be in the first data set. Another
 way to put it is that \c splitpoint will equal the number of rows/columns
 in the first matrix (unless it is -1, in which case the first matrix
 will have zero rows, or it is greater than the matrix's size, in which
 case it will have as many rows as the original).
 \param r_or_c If this is 'r' or 'R', then cleave the rows; of 'c' or 'C' cleave the columns.

 \return An array of two \ref apop_data sets. If one is empty then a
 NULL pointer will be returned.

 */
apop_data ** apop_data_split(apop_data *in, int splitpoint, char r_or_c){
    //A long, dull series of contingencies. Bonus: a valid use of goto.
  apop_data   **out   = malloc(2*sizeof(apop_data *));
  gsl_vector  v1, v2;
  gsl_matrix  m1, m2;
  int       set_v1  = 1,
            set_v2  = 1,
            set_m1  = 1,
            set_m2  = 1;
     if (r_or_c == 'r' || r_or_c == 'r') {
        if (splitpoint <=0){
            out[0]  = NULL;
            out[1]  = apop_data_copy(in);
        } else if (splitpoint >= in->matrix->size1) {
            out[0]  = apop_data_copy(in);
            out[1]  = NULL;
        } else {
            if (in->vector){
                v1      = gsl_vector_subvector(in->vector, 0, splitpoint).vector;
                v2      = gsl_vector_subvector(in->vector, splitpoint, in->vector->size - splitpoint).vector;
            } else
                set_v1  = 
                set_v2  = 0;
            if (in->matrix){
                m1      = gsl_matrix_submatrix (in->matrix, 0, 0, splitpoint, in->matrix->size2).matrix;
                m2      = gsl_matrix_submatrix (in->matrix, splitpoint, 0,
                                    in->matrix->size1 - splitpoint,  in->matrix->size2).matrix;
            } else
                set_m1  = 
                set_m2  = 0;
            goto allocation;
        }
    } else if (r_or_c == 'c' || r_or_c == 'C') {
        if (splitpoint <= -1){
            out[0]  = NULL;
            out[1]  = apop_data_copy(in);
        } else if (splitpoint >= in->matrix->size2){
            out[1]  = NULL;
            out[0]  = apop_data_copy(in);
        } else if (splitpoint == 0){
            if (in->vector)
                v1      = gsl_vector_subvector(in->vector, 0, in->vector->size).vector;
            else set_v1 = 0;
            set_v2  = 0;
            set_m1  = 0;
            if (in->matrix)
                m2      = gsl_matrix_submatrix (in->matrix, 0, 0, 
                                    in->matrix->size1,  in->matrix->size2).matrix;
            else set_m2 = 0;
            goto allocation;
        } else if (splitpoint > 0 && splitpoint < in->matrix->size2){
            if (in->vector)
                v1      = gsl_vector_subvector(in->vector, 0, in->vector->size).vector;
            else set_v1 = 0;
            set_v2  = 0;
            if (in->matrix){
                m1      = gsl_matrix_submatrix (in->matrix, 0, 0, in->matrix->size1, splitpoint).matrix;
                m2      = gsl_matrix_submatrix (in->matrix, 0, splitpoint, 
                                    in->matrix->size1,  in->matrix->size2-splitpoint).matrix;
            } else
                set_m1  = 
                set_m2  = 0;
            goto allocation;
        } else { //splitpoint >= in->matrix->size2
            if (in->vector)
                v1      = gsl_vector_subvector(in->vector, 0, in->vector->size).vector;
            else set_v1 = 0;
            set_v2  = 0;
            if (in->matrix)
                m1      = gsl_matrix_submatrix (in->matrix, 0, 0, 
                            in->matrix->size1, in->matrix->size2).matrix;
            else set_m1 = 0;
            set_m2  = 0;
            goto allocation;
        }
    } else {
        apop_error(0, 'c',"%s: Please set r_or_c == 'r' or == 'c'. Returning two NULLs.\n", __func__);
        out[0]  = NULL;
        out[1]  = NULL;
    }
    return out;

allocation:
    out[0]  = apop_data_alloc(0, 0,0);
    out[1]  = apop_data_alloc(0, 0,0);
    if (set_v1){
        out[0]->vector  = gsl_vector_alloc(v1.size);
        gsl_vector_memcpy(out[0]->vector, &v1);
    }
    if (set_v2){
        out[1]->vector  = gsl_vector_alloc(v2.size);
        gsl_vector_memcpy(out[1]->vector, &v2);
    }
    if (set_m1){
        out[0]->matrix  = gsl_matrix_alloc(m1.size1, m1.size2);
        gsl_matrix_memcpy(out[0]->matrix, &m1);
    }
    if (set_m2){
        out[1]->matrix  = gsl_matrix_alloc(m2.size1, m2.size2);
        gsl_matrix_memcpy(out[1]->matrix, &m2);
    }
    return out;
}




/** Remove the columns set to one in the \c drop vector.
\param n the \ref apop_name structure to be pared down
\param drop  a vector with n->colct elements, mostly zero, with a one marking those columns to be removed.
\ingroup names
 */
static void apop_name_rm_columns(apop_name *n, int *drop){
apop_name   *newname    = apop_name_alloc();
int         i, max      = n->colct;
    for (i=0; i< max; i++){
        if (drop[i]==0)
            apop_name_add(newname, n->column[i],'c');
        else
            n->colct    --;
    }
    free(n->column);
    n->column = newname->column;
    //we need to free the newname struct, but leave the column intact.
    //A one-byte memory leak.
    newname->column   = malloc(1);
    newname->colct  = 0;
    apop_name_free(newname);
}


/** Remove the columns set to one in the \c drop vector.
The returned data structure looks like it was modified in place, but
the data matrix and the names are duplicated before being pared down,
so if your data is taking up more than half of your memory, this may
not work.

\param d the \ref apop_data structure to be pared down. 
\param drop an array of ints. If use[7]==1, then column seven will be cut from the
output. A reminder: <tt>calloc(in->size2 * sizeof(int))</tt> will fill your array with zeros on allocation, and 
<tt>memset(use, 1, in->size2 * sizeof(int))</tt> will
quickly fill an array of ints with nonzero values.
 \ingroup data_struct
 */
void apop_data_rm_columns(apop_data *d, int *drop){
  gsl_matrix  *freeme = d->matrix;
    d->matrix = apop_matrix_rm_columns(d->matrix, drop);
    gsl_matrix_free(freeme); 
    apop_name_rm_columns(d->names, drop);
}

/** Get a pointer to an element of an \c apop_data set. Column -1 is the
  vector. For symmetry with \c gsl_vector_ptr and \c gsl_matrix_ptr.

  \param data   The data set.
  \param    i   The row
  \param    j   The column
  \return       A pointer to double pointing to the appropriate element.
  */
double * apop_data_ptr(const apop_data *data, const int i, const int j){
    if (j == -1){
        apop_assert(data->vector, NULL, 0, 's', "You asked for the vector element (i=-1) but it is NULL.");
        return gsl_vector_ptr(data->vector, i);
    } else {
        apop_assert(data->matrix, NULL, 0, 's', "You asked for the matrix element (%i, %i) but the matrix is NULL.", i, j);
        return gsl_matrix_ptr(data->matrix, i,j);
    }
}

/** Get a pointer to an element from an \ref apop_data set, using the row name but
 the column number.

Uses \ref apop_name_find for the search; see notes there on the name
matching rules.

\param  in  the \ref apop_data set.
\param  row the name of the row you seek.
\param  col the number of the column. If -1, return the vector element.

 \ingroup data_struct
 */
double *apop_data_ptr_ti(const apop_data *in, char* row, int col){
  int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert(rownum != -1,  NULL, 0,'c',"Couldn't find %s amongst the row names.", row);
    if (col>=0)
        return gsl_matrix_ptr(in->matrix, rownum, col);
    else
        return gsl_vector_ptr(in->vector, rownum);
}

/** Get a pointer to an element from an \ref apop_data set, using the column name but
 the row number.

Uses \ref apop_name_find for the search; see notes there on the name
matching rules.

 \ingroup data_struct
 */
double *apop_data_ptr_it(const apop_data *in, size_t row, char* col){
  int colnum =  apop_name_find(in->names, col, 'c');
    apop_assert(colnum != -1,  NULL, 0,'c',"Couldn't find %s amongst the column names.", col);
    return gsl_matrix_ptr(in->matrix, row, colnum);
}

/** Get a pointer to an element from an \ref apop_data set, using the row and column name.

Uses \ref apop_name_find for the search; see notes there on the name
matching rules.

 \ingroup data_struct
 */
double *apop_data_ptr_tt(const apop_data *in, char *row, char* col){
  int colnum =  apop_name_find(in->names, col, 'c');
  int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert(colnum != -1,  NULL, 0,'c',"Couldn't find %s amongst the column names.", col);
    apop_assert(rownum != -1,  NULL, 0,'c',"Couldn't find %s amongst the row names.", row);
    return gsl_matrix_ptr(in->matrix, rownum, colnum);
}

/** Returns the data element at the given point.

Q: How does <tt> apop_data_get(in, r, c)</tt> differ from
 <tt>gsl_matrix_get(in->matrix, r, c)</tt>?\\
A: It's nine characters shorter. Also, if \c c==-1, then this will
return the \c apop_data's vector element.

 \ingroup data_struct
*/
double apop_data_get(const apop_data *in, size_t row, int col){
    if (col>=0){
        apop_assert(in->matrix, 0, 0, 's', "You asked for the matrix element (%u, %i) but the matrix is NULL.", row, col);
        return gsl_matrix_get(in->matrix, row, col);
    } else {
        apop_assert(in->vector, 0, 0, 's', "You asked for the vector element (col=-1) but it is NULL.");
        return gsl_vector_get(in->vector, row);
    }
}

/** Get an element from an \ref apop_data set, using the row name but
 the column number

Uses \ref apop_name_find for the search; see notes there on the name
matching rules.

\param  in  the \ref apop_data set.
\param  row the name of the row you seek.
\param  col the number of the column. If -1, return the vector element.

 \ingroup data_struct
 */
double apop_data_get_ti(const apop_data *in, char* row, int col){
  int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert(rownum != -1,  GSL_NAN, 0,'c',"Couldn't find %s amongst the row names.", row);
    if (col>=0)
        return gsl_matrix_get(in->matrix, rownum, col);
    else
        return gsl_vector_get(in->vector, rownum);
}

/** Get an element from an \ref apop_data set, using the column name but
 the row number

Uses \ref apop_name_find for the search; see notes there on the name
matching rules.

 \ingroup data_struct
 */
double apop_data_get_it(const apop_data *in, size_t row, char* col){
  int colnum =  apop_name_find(in->names, col, 'c');
    apop_assert(colnum != -1,  GSL_NAN, 0,'c',"Couldn't find %s amongst the column names.", col);
    return gsl_matrix_get(in->matrix, row, colnum);
}

/** Get an element from an \ref apop_data set, using the row and column name.

Uses \ref apop_name_find for the search; see notes there on the name
matching rules.

 \ingroup data_struct
 */
double apop_data_get_tt(const apop_data *in, char *row, char* col){
  int colnum =  apop_name_find(in->names, col, 'c');
  int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert(colnum != -1,  GSL_NAN, 0,'c',"Couldn't find %s amongst the column names.", col);
    apop_assert(rownum != -1,  GSL_NAN, 0,'c',"Couldn't find %s amongst the row names.", row);
    return gsl_matrix_get(in->matrix, rownum, colnum);
}

/** Sets the data element at the given point.

Q: How does <tt> apop_data_set(in, row, col, data)</tt> differ from
 <tt>gsl_matrix_set(in->matrix, row, col, data)</tt>?\\
A: It's seven characters shorter.

Oh, and if \c col<0, then this will set the element of \c in->vector.
 \ingroup data_struct
*/
void apop_data_set(apop_data *in, size_t row, int col, double data){
    if (col>=0){
        apop_assert_void(in->matrix, 0, 's', "You're trying to set the matrix element (%u, %i) but the matrix is NULL.", row, col);
        gsl_matrix_set(in->matrix, row, col, data);
    } else {
        apop_assert_void(in->vector, 0, 's', "You're trying to set a vector element (row=-1) but the vector is NULL.");
        gsl_vector_set(in->vector, row, data);
    }
}
/** Set an element from an \ref apop_data set, using the row name but
 the column number

Uses \ref apop_name_find for the search; see notes there on the name
matching rules.

\param  in  the \ref apop_data set.
\param  row the name of the row you seek.
\param  col the number of the column. If -1, set the vector element.
\param  data the value to insert

 \ingroup data_struct
 */
void apop_data_set_ti(apop_data *in, char* row, int col, double data){
  int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert_void(rownum != -1, 0,'c',"Couldn't find %s amongst the row names. Making no changes.", row);
    if (col>=0)
        gsl_matrix_set(in->matrix, rownum, col, data);
    else
        gsl_vector_set(in->vector, rownum, data);
}

/** Set an element from an \ref apop_data set, using the column name but
 the row number

Uses \ref apop_name_find for the search; see notes there on the name
matching rules.

 \ingroup data_struct
 */
void apop_data_set_it(apop_data *in, size_t row, char* col, double data){
  int colnum =  apop_name_find(in->names, col, 'c');
    apop_assert_void(colnum != -1, 0,'c',"Couldn't find %s amongst the column names. Making no changes.", col);
    gsl_matrix_set(in->matrix, row, colnum, data);
}

/** Set an element from an \ref apop_data set, using the row and column name.

Uses \ref apop_name_find for the search; see notes there on the name
matching rules.

 \ingroup data_struct
 */
void apop_data_set_tt(apop_data *in, char *row, char* col, double data){
  int colnum =  apop_name_find(in->names, col, 'c');
  int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert_void(colnum != -1, 0, 'c',"Couldn't find %s amongst the column names.", col);
    apop_assert_void(rownum != -1, 0, 'c',"Couldn't find %s amongst the column names.", row);
    gsl_matrix_set(in->matrix, rownum, colnum, data);
}

/** Add a named element to a data vector. For example, this is primarily
used by the testing procedures, that produce a column of named parameters.

\param d    The \ref apop_data structure.
\param name The name to add
\param val  the vector value to add.

I use the position of the name to know where to put the value. If there
are two names in the data set, then I will put the data in the third
slot in the vector. If you use this function from start to finish in
building your vector, then you'll be fine.

*/
void apop_data_add_named_elmt(apop_data *d, char *name, double val){
    gsl_vector_set(d->vector, d->names->rowct, val);
    apop_name_add(d->names, name, 'r');
}



/** Add a string to the text element of an \c apop_data set. 


   By the way, the \c text element is simply an array of arrays of
   strings, so there are no special tricks to this function. After some
   checks, it just runs:

\code
in->text[row][col]  = malloc(strlen(yourtext)+1);
strcpy(in->text[row][col], yourtext);
\endcode

\param in   The \c apop_data set, that already has an allocated \c text element.
\param row  The row
\param col  The col
\param text The text to write.
*/
void apop_text_add(apop_data *in, const size_t row, const size_t col, const char *text){
    apop_assert_void((in->textsize[0] >= (int)row-1) && (in->textsize[1] >= (int)col-1), 0, 'c', "You asked me to put the text '%s' at (%i, %i), but the text array has size (%i, %i)\n", text, row, col, in->textsize[0], in->textsize[1]);
    if (text==NULL){
        in->text[row][col]  = malloc(strlen("NaN") +1);
        strcpy(in->text[row][col], "NaN");
    } else {
        in->text[row][col]  = malloc(strlen(text)+1);
        strcpy(in->text[row][col], text);
    }
}

/** This allocates an array of strings and puts it in the \c text element
  of an \c apop_data set. 

  \param in An \c apop_data set. It's OK to send in \c NULL, in which case an apop_data set with \c NULL \c matrix and \vector elements is returned.
  \param row    the number of rows of text.
  \param col     the number of columns of text.
  \return       A pointer to the relevant \c apop_data set. If the input was not \c NULL, then this is a repeat of the input pointer.
  */
apop_data * apop_text_alloc(apop_data *in, const size_t row, const size_t col){
  int   i, j;
    if (!in)
        in  = apop_data_alloc(0,0,0);
    in->text = malloc(sizeof(char**) * row);
    for (i=0; i< row; i++){
        in->text[i] = malloc(sizeof(char*) * col);
        for (j=0; j< col; j++)
            in->text[i][j] = NULL;
    }
    in->textsize[0] = row;
    in->textsize[1] = col;
    return in;
}



/** Transpose the matrix element of the input \ref apop_data set,
 including the row/column names. The vector and text elements of the input data set are completely ignored.

 This is really just a friendly wrapper for \c gsl_matrix_transpose_memcpy; if you have a \c gsl_matrix with no names, you may prefer to just use that function.

 \param in The input \ref apop_data set.
 \return  A newly alloced \ref apop_data set, with the appropriately transposed matrix. The vector and text elements will be \c NULL.
 */ 
apop_data *apop_data_transpose(apop_data *in){
  apop_assert(in->matrix, NULL, 0, 'c', "input data set has no matrix element, so I'm returning NULL.");
  apop_data *out = apop_data_alloc(0, in->matrix->size2, in->matrix->size1);
    gsl_matrix_transpose_memcpy(out->matrix, in->matrix);
    apop_name_cross_stack(out->names, in->names, 'r', 'c');
    apop_name_cross_stack(out->names, in->names, 'c', 'r');
    return out;
}
