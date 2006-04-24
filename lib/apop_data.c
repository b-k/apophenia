/** \file apop_data.c	 The apop_data structure joins together a
  gsl_matrix, apop_name, and a table of strings. No biggie.


Copyright (c) 2006 by Ben Klemens. Licensed under the GNU GPL v2.
*/

#include <gsl/gsl_matrix.h>
#include "apophenia/types.h"
#include "apophenia/output.h"
#include "apophenia/linear_algebra.h"

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
    printf("%s\t", apop_get_name(data->names, j, 'r'));
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

/** Allocate a \ref apop_data structure, to be filled with data. If
\c size2>0, then the matrix will be allocated; else, the vector will
be allocated. Your best bet for allocating the categories is to produce
them elsewhere, such as \ref apop_query_to_chars and then point an
apop_data structure with a zero-sized vector to your matrix of strings.

  \param size1, size2   row and column size for the matrix. If \c size2>0
  this exactly mirrors the format of \c gsl_matrix_alloc. If \c size2==-1,
  then allocate a vector. \c apop_data_alloc(0,0) will produce a basically blank set, with \c out->matrix==out->vector==NULL. 

 \return    The \ref apop_data structure, allocated and ready.
 \ingroup data_struct
  */
apop_data * apop_data_alloc(int size1, int size2){
apop_data  *setme       = malloc(sizeof(apop_data));
    setme->vector   = NULL;
    setme->matrix   = NULL;
    if (size2 > 0)
        setme->matrix   = gsl_matrix_alloc(size1,size2);
    else if (size1>0)
        setme->vector   = gsl_vector_alloc(size1);

    setme->names        = apop_name_alloc();
    setme->categories   = NULL;
    setme->catsize[0]   = 
    setme->catsize[1]   = 0;
    return setme;
}

/** Wrap an \ref apop_name structure around an existing \c gsl_matrix.

\param m    The existing matrix you'd like to turn into an \ref apop_data structure.
return      The \ref apop_data structure in question.
  */
apop_data * apop_data_from_matrix(gsl_matrix *m){
apop_data  *setme   = apop_data_alloc(0,0);
    if (m==NULL && apop_opts.verbose) 
        {printf("Warning: converting a NULL matrix to an apop_data structure.\n");}
    setme->matrix       = m;
    return setme;
}

/** A synonym for \ref apop_data_from_matrix. Use that one.
 */
apop_data * apop_matrix_to_data(gsl_matrix *m){
    return apop_data_from_matrix(m);
}

/** Wrap an \ref apop_name structure around an existing \c gsl_vector.

    A synonym for \ref apop_vector_to_data.

\param  v   The data vector
\return     an allocated, ready-to-use \ref apop_data struture.
*/
apop_data * apop_data_from_vector(gsl_vector *v){
apop_data  *setme   = malloc(sizeof(apop_data));
    if (v==NULL && apop_opts.verbose) 
        {printf("Warning: converting a NULL matrix to an apop_data structure.\n");}
    setme->vector       = v;
    setme->names        = apop_name_alloc();
    setme->categories   = NULL;
    setme->catsize[0]   = 
    setme->catsize[1]   = 0;
    return setme;
}


/** Wrap an \ref apop_name structure around an existing \c gsl_vector.

\param  v   The data vector
\return     an allocated, ready-to-use \ref apop_data struture.
*/
apop_data * apop_vector_to_data(gsl_vector *v){
    return apop_data_from_vector(v);
}

/** free a matrix of chars* (i.e., a char***).

  */
void apop_cats_free(char ***freeme, int rows, int cols){
int     i,j;
    if (rows && cols){
        for (i=0; i < rows; i++){
            for (j=0; j < cols; j++)
                if(freeme[i][j])
                    free(freeme[i][j]);
            if(freeme[i])
                free(freeme[i]);
        }
    }
        if(freeme)
            free(freeme);
}


/** Free an \ref apop_name structure.

 \ingroup data_struct
  */
void apop_data_free(apop_data *freeme){
    if (freeme->vector)
        gsl_vector_free(freeme->vector);
    if (freeme->matrix)
        gsl_matrix_free(freeme->matrix);
    apop_name_free(freeme->names);
    apop_cats_free(freeme->categories, freeme->catsize[0] , freeme->catsize[1]);
    /*
    if (freeme->catsize[0] && freeme->catsize[1]){
        for (i=0; i < freeme->catsize[0]; i++)
            for (j=0; j < freeme->catsize[1]; j++)
                if(freeme->categories[i][j])
                    free(freeme->categories[i][j]);
        if(freeme->categories)
            free(freeme->categories);
    }
    */
    free(freeme);
}


/** Copy one \ref apop_data structure to another. That is, all data is duplicated.

  This function does <i>not</i> allocate the output for you. If you want such behavior, usr \ref apop_data_copy.
 
  \param out    a structure that this function will fill. Must be preallocated
  \param in    the input data

 \ingroup data_struct
 \todo This doesn't copy over the category data.
  */
void apop_data_memcpy(apop_data *out, apop_data *in){
    if (!out)
        printf("apop_data_mecpy: you are copying to a NULL vector. Do you mean to use apop_data_copy instead?\n");
    if (in->matrix){
        if (in->matrix->size1 != out->matrix->size1 ||
                in->matrix->size2 != out->matrix->size2){
            if (apop_opts.verbose)
                printf("You're trying to copy a (%i X %i) into a (%i X %i) matrix. Returning w/o any copying.\n", 
                in->matrix->size1, in->matrix->size2, 
                out->matrix->size1, out->matrix->size2);
            return;
        }
        gsl_matrix_memcpy(out->matrix, in->matrix);
    }
    if (in->vector){
        if (in->vector->size != out->vector->size){
            if (apop_opts.verbose)
                printf("You're trying to copy a %i-elmt vector into a %i-elmt vector. Returning w/o any copying.\n", 
                in->vector->size, out->vector->size);
            return;
        }
        gsl_vector_memcpy(out->vector, in->vector);
    }
    apop_name_stack(out->names, in->names, 'r');
    apop_name_stack(out->names, in->names, 'c');
    apop_name_stack(out->names, in->names, 't');
    if (in->catsize[0] && in->catsize[1]){
        out->categories  = malloc(sizeof(char **) * in->catsize[0] * in->catsize[1]);
        memcpy( out->categories, in->categories, sizeof(char **) * in->catsize[0] * in->catsize[1]);
    }
}

/** Copy one \ref apop_data structure to another. That is, all data is duplicated.

  Just a front-end for \ref apop_data_memcpy for those who prefer this sort of syntax.
 
  \param in    the input data
  \return       a structure that this function will allocate and fill

 \ingroup data_struct
  */
apop_data *apop_data_copy(apop_data *in){
apop_data *out  = apop_data_alloc(0, 0);
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
twice as much memory. Plan accordingly. If you are stacking matrices, the output vector and categories are NULL; if you are stacking vectors, the output matrix and categories are NULL. 

The dependent variable names are not copied or modified in any way by
this function, so if you mean for them to change or get stacked, you
will have to make the modifications yourself, perhaps by using either \ref apop_name_add or \ref apop_name_stack.

\param  m1      the upper/rightmost data set
\param  m2      the second data set
\param  posn    if 'r', stack rows of m1's matrix above rows of m2's<br>
if 'c', stack columns of m1's matrix to right of m2's<br>
if 'v', stack m1's vector to over of m2's<br>

\return         a new \ref apop_data set with the stacked data.
\ingroup data_struct
*/
apop_data *apop_data_stack(apop_data *m1, apop_data * m2, char posn){
apop_data   *out    = NULL;
    if (posn == 'r' || posn == 'c'){
        gsl_matrix  *stacked= apop_matrix_stack(m1->matrix, m2->matrix, posn);
        out         = apop_matrix_to_data(stacked);
        out->names  = apop_name_copy(m1->names);
        apop_name_stack(out->names, m2->names, posn);
    } else if (posn == 'v'){
        gsl_vector  *stacked= apop_vector_stack(m1->vector, m2->vector);
        out         = apop_vector_to_data(stacked);
        out->names  = apop_name_copy(m1->names);
        apop_name_stack(out->names, m2->names, posn);
    } else{
        printf("apop_data_stack: valid positions are 'r', 'c', or 'v'; you gave me >%c<. Returning NULL.", posn);
    }
    return out;
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

/** Returns the data element at the given point.

Q: How does <tt> apop_data_get(in, r,c)</tt> differ from
 <tt>gsl_matrix_get(in->matrix, row, col)</tt>?\\
A: It's seven characters shorter. Oh, and if \c c==-1, then this will
return the \c apop_data's vector element.

 \ingroup data_struct
*/
double apop_data_get(apop_data *in, size_t row, int col){
    if (col>=0)
        return gsl_matrix_get(in->matrix, row, col);
    else
        return gsl_vector_get(in->vector, row);
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
double apop_data_get_tn(apop_data *in, char* row, int col){
int rownum =  apop_name_find(in->names, row, 'r');
    if (rownum == -1){
        if(apop_opts.verbose)
            printf("couldn't find %s amongst the column names.\n",row);
        return GSL_NAN;
    }
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
double apop_data_get_nt(apop_data *in, size_t row, char* col){
int colnum =  apop_name_find(in->names, col, 'c');
    if (colnum == -1){
        if(apop_opts.verbose)
            printf("couldn't find %s amongst the column names.\n",col);
        return GSL_NAN;
    }
    return gsl_matrix_get(in->matrix, row, colnum);
}

/** Get an element from an \ref apop_data set, using the row and column name.

Uses \ref apop_name_find for the search; see notes there on the name
matching rules.

 \ingroup data_struct
 */
double apop_data_get_tt(apop_data *in, char *row, char* col){
int colnum =  apop_name_find(in->names, col, 'c');
int rownum =  apop_name_find(in->names, row, 'r');
    if (colnum == -1){
        if(apop_opts.verbose)
            printf("couldn't find %s amongst the column names.\n",col);
        return GSL_NAN;
    }
    if (rownum == -1){
        if(apop_opts.verbose)
            printf("couldn't find %s amongst the column names.\n",row);
        return GSL_NAN;
    }
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
    if (col>=0)
        gsl_matrix_set(in->matrix, row, col, data);
    else
        gsl_vector_set(in->vector, row, data);
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
void apop_data_set_tn(apop_data *in, char* row, int col, double data){
int rownum =  apop_name_find(in->names, row, 'r');
    if (rownum == -1){
        if(apop_opts.verbose)
            printf("couldn't find %s amongst the column names.\n",row);
    }
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
void apop_data_set_nt(apop_data *in, size_t row, char* col, double data){
int colnum =  apop_name_find(in->names, col, 'c');
    if (colnum == -1){
        if(apop_opts.verbose)
            printf("couldn't find %s amongst the column names.\n",col);
    }
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
    if (colnum == -1){
        if(apop_opts.verbose)
            printf("couldn't find %s amongst the column names.\n",col);
    }
    if (rownum == -1){
        if(apop_opts.verbose)
            printf("couldn't find %s amongst the column names.\n",row);
    }
    gsl_matrix_set(in->matrix, rownum, colnum, data);
}
