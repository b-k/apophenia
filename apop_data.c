/** \file 
The apop_data structure joins together a gsl_matrix, apop_name, and a table of strings. */
/* Copyright (c) 2006--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "types.h"
#include "output.h"
#include "conversions.h"
#include "linear_algebra.h"

/** \defgroup data_struct apop_data

  The \c apop_data structure represents a data set.  It joins together a gsl_vector, a gsl_matrix, an apop_name, and a table of strings. It tries to be minimally intrusive, so you can use it everywhere you would use a \c gsl_matrix or a \c gsl_vector.

 For example, let us say that you are running a regression: there is a vector for the dependent variable, and a matrix for the dependent variables. Think of them as a partitioned matrix, where the vector is column -1, and the first column of the matrix is column zero. Here is some code to print the entire matrix. Notice that the column counter \c i starts counting at -1.

  \code
  for (j = 0; j< data->matrix->size1; j++){
    printf("%s\t", data->names->row[j]);
    for (i = -1; i< data->matrix->size2; i++)
        printf("%g\t", apop_data_get(data, j, i));
    printf("\n");
    }
    \endcode

We're generally assuming that the data vector and data matrix have the same row count: \c data->vector->size==data->matrix->size1 . This means that the \ref apop_name structure doesn't have separate vector_names and row_names elements: the rownames are assumed to apply for both.

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

\see{apop_data_calloc}

  \param vsize              vector size, if any. 
  \param msize1, msize2     Row and column size for the matrix. 
  
  \c apop_data_alloc(0,0,0) will produce a basically blank set, with \c out->matrix==out->vector==NULL. 

 \return    The \ref apop_data structure, allocated and ready.
 \ingroup data_struct
  */
apop_data * apop_data_alloc(const size_t vsize, const size_t msize1, const int msize2){
  apop_data  *setme       = malloc(sizeof(apop_data));
    apop_assert(setme, NULL, 0, 's', "malloc failed. Probably out of memory.");
    *setme = (apop_data) { }; //init to zero/NULL.
    if (msize2 > 0  && msize1 > 0){
        setme->matrix   = gsl_matrix_alloc(msize1,msize2);
        apop_assert(setme->matrix, NULL, 0, 's', "malloc failed on a %zu x %i matrix. Probably out of memory.", msize1, msize2);
    }
    if (vsize){
        setme->vector   = gsl_vector_alloc(vsize);
        apop_assert(setme->vector, NULL, 0, 's', "malloc failed on a vector of size %zu. Probably out of memory.", vsize);
    }
    //I allocate a vector of msize2==-1. This is deprecated, and will one day be deleted.
    else if (msize2==-1 && msize1>0){
        setme->vector   = gsl_vector_alloc(msize1);
        apop_assert(setme->vector, NULL, 0, 's', "malloc failed on a vector of size %zu. Probably out of memory.", msize1);
    }
    setme->names        = apop_name_alloc();
    return setme;
}

/** Allocate a \ref apop_data structure, to be filled with data; set everything in the allocated portion to zero. See \ref apop_data_alloc for details.

  \param vsize              vector size, if any. 
  \param msize1, msize2     Row and column size for the matrix. If \c size2>0
  this exactly follows the format of \c gsl_matrix_calloc. If \c size2==-1,
  then allocate a vector. 
  
  \c apop_data_calloc(0, 0,0) will produce a basically blank set, with \c out->matrix==out->vector==NULL. 

\return    The \ref apop_data structure, allocated and zeroed out.
\see apop_data_alloc 
\ingroup data_struct
  */
apop_data * apop_data_calloc(const size_t vsize, const size_t msize1, const int msize2){
  apop_data  *setme       = malloc(sizeof(apop_data));
    apop_assert(setme, NULL, 0, 's', "malloc failed. Probably out of memory.");
    *setme = (apop_data) { }; //init to zero/NULL.
    if (msize2 >0 && msize1 > 0){
        setme->matrix   = gsl_matrix_calloc(msize1,msize2);
        apop_assert(setme->matrix, NULL, 0, 's', "malloc failed on a %zu x %i matrix. Probably out of memory.", msize1, msize2);
    }
    if (vsize){
        setme->vector   = gsl_vector_calloc(vsize);
        apop_assert(setme->vector, NULL, 0, 's', "malloc failed on a vector of size %zu. Probably out of memory.", vsize);
    }
    else if (msize2==-1 && msize1>0){
        setme->vector   = gsl_vector_calloc(msize1);
        apop_assert(setme->vector, NULL, 0, 's', "malloc failed on a vector of size %zu. Probably out of memory.", msize1);
    }
    setme->names        = apop_name_alloc();
    return setme;
}

/** Wrap an \ref apop_data structure around an existing \c gsl_matrix.
 The matrix is not copied, but is pointed to by the new \ref apop_data struct.

\param m    The existing matrix you'd like to turn into an \ref apop_data structure.
\return      The \ref apop_data structure whose \c matrix pointer points to the input matrix. The rest of the struct is basically blank.
  */
apop_data * apop_matrix_to_data(gsl_matrix *m){
  apop_assert_void(m, 1, 'c',"Converting a NULL matrix to an apop_data structure.");
  apop_data  *setme   = apop_data_alloc(0,0,0);
    setme->matrix = m;
    return setme;
}

/** Wrap an \ref apop_data structure around an existing \c gsl_vector.
 The vector is not copied, but is pointed to by the new \ref apop_data struct.

\param  v   The data vector
\return     an allocated, ready-to-use \ref apop_data structure.
*/
apop_data * apop_vector_to_data(gsl_vector *v){
  apop_assert(v, NULL, 1, 'c',"Converting a NULL vector to an apop_data structure.");
  apop_data  *setme   = apop_data_alloc(0,0,0);
    setme->vector = v;
    return setme;
}

/** Free a matrix of chars* (i.e., a char***). This is the form of the
 text element of the \ref apop_data set, so you can use this for:
 \code
 apop_text_free(yourdata->text, yourdata->textsize[0], yourdata->textsize[1]);
 \endcode
 This is what \c apop_data_free uses internally.
   */
void apop_text_free(char ***freeme, int rows, int cols){
    if (rows && cols)
        for (int i=0; i < rows; i++){
            for (int j=0; j < cols; j++)
                free(freeme[i][j]);
            free(freeme[i]);
        }
    free(freeme);
}

/** Free an \ref apop_data structure.
 
As with \c free(), it is safe to send in a \c NULL pointer (in which case the function does nothing).

If the \c more pointer is not \c NULL, I will free the pointed-to data set first.
If you don't want to free data sets down the chain, set <tt>more=NULL</tt> before calling this.

 \ingroup data_struct
  */
void apop_data_free(apop_data *freeme){
    if (!freeme) return;
    if (freeme->more)
        apop_data_free(freeme->more);
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

  This function does not allocate the output structure for you for the overall structure or the vector or matrix. If you want such behavior, use \ref apop_data_copy. Both functions do allocate memory for the text.

  I don't follow the \c more pointer, though \ref apop_data_copy does.
 
  \param out    a structure that this function will fill. Must be preallocated
  \param in    the input data

 \ingroup data_struct
  */
void apop_data_memcpy(apop_data *out, const apop_data *in){
    apop_assert_void(out, 1, 'c', "you are copying to a NULL vector. "
                                  "Do you mean to use apop_data_copy instead?");
    if (in->matrix){
        apop_assert_void(in->matrix->size1 == out->matrix->size1 && in->matrix->size2 == out->matrix->size2, 
                1, 'c',"You're trying to copy a (%zu X %zu) into a (%zu X %zu) matrix. Returning w/o any copying.", 
                        in->matrix->size1, in->matrix->size2, out->matrix->size1, out->matrix->size2);
        gsl_matrix_memcpy(out->matrix, in->matrix);
    }
    if (in->vector){
        apop_assert_void(in->vector->size == out->vector->size,
                1, 'c',"You're trying to copy a %zu-elmt vector into a %zu-elmt vector. Returning w/o any copying.", 
                                 in->vector->size, out->vector->size);
        gsl_vector_memcpy(out->vector, in->vector);
    }
    if (in->weights){
        apop_assert_void(in->weights->size == out->weights->size,
                1, 'c',"Weight vector sizes don't match: you're trying "
                        "to copy a %zu-elmt vector into a %zu-elmt vector. "
                        "Weights not copied.", 
                                 in->weights->size, out->weights->size);
        gsl_vector_memcpy(out->weights, in->weights);
    }
    apop_name_free(out->names);
    out->names = apop_name_alloc();
    apop_name_stack(out->names, in->names, 'r');
    apop_name_stack(out->names, in->names, 'c');
    apop_name_stack(out->names, in->names, 't');
    out->textsize[0] = in->textsize[0]; 
    out->textsize[1] = in->textsize[1]; 
    if (in->textsize[0] && in->textsize[1]){
        out->text  = malloc(sizeof(char ***) * in->textsize[0] * in->textsize[1]);
        for (size_t i=0; i< in->textsize[0]; i++){
            out->text[i]  = malloc(sizeof(char **) * in->textsize[1]);
            for(size_t j=0; j < in->textsize[1]; j ++)
				asprintf(&(out->text[i][j]), "%s", in->text[i][j]);
        }
    }
}

/** Copy one \ref apop_data structure to another. That is, all data is duplicated.

  Basically a front-end for \ref apop_data_memcpy for those who prefer this sort of syntax. 

  Unlike \ref apop_data_memcpy, I do follow the \c more pointer.
 
  \param in    the input data
  \return       a structure that this function will allocate and fill. If input is NULL, then this will be NULL.

 \ingroup data_struct
  */
apop_data *apop_data_copy(const apop_data *in){
    if (!in)
        return NULL;
    apop_data *out  = apop_data_alloc(0, 0, 0);
    if (in->more)
        out->more = apop_data_copy(in->more);
    if (in->vector)
        out->vector = gsl_vector_alloc(in->vector->size);
    if (in->matrix)
        out->matrix = gsl_matrix_alloc(in->matrix->size1, in->matrix->size2);
    if (in->weights)
        out->weights = gsl_vector_alloc(in->weights->size);
    apop_data_memcpy(out, in);
    return out;
}

/** Put the first data set either on top of or to the left of the second matrix.

The fn returns a new data set, meaning that at the end of this function,
until you apop_data_free() the original data sets, you will be taking up
twice as much memory. Plan accordingly. 

 For the opposite operation, see \ref apop_data_split.

\param  m1      the upper/rightmost data set (default = \c NULL)
\param  m2      the second data set (default = \c NULL)
\param  posn    If 'r', stack rows of m1's matrix above rows of m2's<br>
if 'c', stack columns of m1's matrix to left of m2's<br>
(default = 'r')
\param  inplace If \c 'i' \c 'y' or 1, use \ref apop_matrix_realloc and \ref apop_vector_realloc to modify \c m1 in place; see the caveats on those function. Otherwise, allocate a new vector, leaving \c m1 unmolested. (default='n')
\return         The stacked data, either in a new \ref apop_data set or \c m1

\li If m1 or m2 are NULL, this returns a copy of the other element, and if
both are NULL, you get NULL back (except if \c m1 is \c NULL and \c inplace is \c 'y', where you'll get the original \c m1 back)
\li Text is ignored.
\li \c more is ignored.
\li If stacking rows on rows, the output vector is the input
vectors stacked accordingly. If stacking columns by columns, the output
vector is just a copy of the vector of m1 and m2->vector doesn't appear in the
output at all.  
\li The same rules for dealing with the vector(s) hold for the vector(s) of weights.
\li Names are a copy of the names for \c m1, with the names for \c m2 appended to the row or column list, as appropriate.

This function uses the \ref designated syntax for inputs.
\ingroup data_struct
*/
APOP_VAR_HEAD apop_data *apop_data_stack(apop_data *m1, apop_data * m2, char posn, char inplace){
    apop_data * apop_varad_var(m1, NULL)
    apop_data * apop_varad_var(m2, NULL)
    char apop_varad_var(posn, 'r')
    char apop_varad_var(inplace, 'n')
    inplace = (inplace == 'i' || inplace == 'y' || inplace == 1 || inplace == 'I' || inplace == 'Y') ? 1 : 0;
    return apop_data_stack_base(m1, m2, posn, inplace);
APOP_VAR_ENDHEAD
  gsl_matrix  *stacked= NULL;
  apop_data   *out    = NULL;
    if (!m1)
        return apop_data_copy(m2);
    if (!m2)
        return (inplace == 'i' || inplace == 'y') ? m1 : apop_data_copy(m1);
    apop_assert((posn == 'r' || posn == 'c'), NULL, 0, 'c', "Valid positions are 'r' or 'c'"
                             " you gave me '%c'. Returning NULL.", posn);
    if (inplace){
        out = m1;
        apop_matrix_stack(m1->matrix, m2->matrix, posn, inplace);
    } else {
        stacked = apop_matrix_stack(m1->matrix, m2->matrix, posn);
        out     = apop_matrix_to_data(stacked);
    }
    if (posn == 'r'){
        out->vector = apop_vector_stack(m1->vector, m2->vector, inplace);
        out->weights = apop_vector_stack(m1->weights, m2->weights, inplace);
    } else if (!inplace){
        out->vector = apop_vector_copy(m1->vector);
        out->weights = apop_vector_copy(m1->weights);
    }
    if (!inplace)
        out->names  = apop_name_copy(m1->names);
    apop_name_stack(out->names, m2->names, posn);
    return out;
}

/** Split one input \ref apop_data structure into two.

 For the opposite operation, see \ref apop_data_stack.

 The \ref apop_data->vector is taken to be the -1st element of the matrix.  
 
\param in  The \ref apop_data structure to split 
\param splitpoint The index of what will be the first row/column of the second data set.  E.g., if this is -1 and \c r_or_c=='c', then the whole data set will be in the second data set; if this is the length of the matrix then the whole data set will be in the first data set. Another way to put it is that \c splitpoint will equal the number of rows/columns in the first matrix (unless it is -1, in which case the first matrix will have zero rows, or it is greater than the matrix's size, in which case it will have as many rows as the original).  
\param r_or_c If this is 'r' or 'R', then cleave the rows; of 'c' or 'C' cleave the columns.

 \return An array of two \ref apop_data sets. If one is empty then a
 NULL pointer will be returned in that position. For example, for a data set of 50 rows, <tt>apop_data **out = apop_data_split(data, 100, 'r')</tt> sets <tt>out[0] = apop_data_copy(data)</tt> and <tt>out[1] = NULL</tt>.

 \li \c more pointer is ignored.
 \li Weights will be preserved. If splitting by rows, then the top and bottom parts of the weights vector will be assigned to the top and bottom parts of the main data set. If splitting by columns, identical copies of the weights vector will be assigned to both parts.
 */
apop_data ** apop_data_split(apop_data *in, int splitpoint, char r_or_c){
    //A long, dull series of contingencies. Bonus: a reasonable use of goto.
  apop_data   **out   = malloc(2*sizeof(apop_data *));
  out[0] = out[1] = NULL;
  gsl_vector  v1, v2, w1, w2;
  gsl_matrix  m1, m2;
  int       set_v1  = 1,
            set_v2  = 1,
            set_m1  = 1,
            set_m2  = 1,
            set_w1  = 1,
            set_w2  = 1;
     if (r_or_c == 'r' || r_or_c == 'r') {
        if (splitpoint <=0)
            out[1]  = apop_data_copy(in);
        else if (splitpoint >= in->matrix->size1)
            out[0]  = apop_data_copy(in);
        else {
            if (in->vector){
                v1 = gsl_vector_subvector(in->vector, 0, splitpoint).vector;
                v2 = gsl_vector_subvector(in->vector, splitpoint, in->vector->size - splitpoint).vector;
            } else
                set_v1  = 
                set_v2  = 0;
            if (in->weights){
                w1 = gsl_vector_subvector(in->weights, 0, splitpoint).vector;
                w2 = gsl_vector_subvector(in->weights, splitpoint,
                        in->weights->size - splitpoint).vector;
            } else
                set_w1  = 
                set_w2  = 0;
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
        if (in->weights){
            w1      = gsl_vector_subvector(in->weights, 0, in->weights->size).vector;
            w2      = gsl_vector_subvector(in->weights, 0, in->weights->size).vector;
        } else 
            set_w1 = 
            set_w2 = 0;

        if (splitpoint <= -1)
            out[1]  = apop_data_copy(in);
        else if (splitpoint >= in->matrix->size2)
            out[0]  = apop_data_copy(in);
        else if (splitpoint == 0){
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
    } else 
        apop_error(0, 'c',"%s: Please set r_or_c == 'r' or == 'c'. Returning two NULLs.\n", __func__);
    return out;

allocation:
    out[0]  = apop_data_alloc(0, 0,0);
    out[1]  = apop_data_alloc(0, 0,0);
    if (set_v1) out[0]->vector  = apop_vector_copy(&v1);
    if (set_v2) out[1]->vector  = apop_vector_copy(&v2);
    if (set_m1) out[0]->matrix  = apop_matrix_copy(&m1);
    if (set_m2) out[1]->matrix  = apop_matrix_copy(&m2);
    if (set_w1) out[0]->weights  = apop_vector_copy(&w1);
    if (set_w2) out[1]->weights  = apop_vector_copy(&w2);
    return out;
}

/** Remove the columns set to one in the \c drop vector.
\param n the \ref apop_name structure to be pared down
\param drop  a vector with n->colct elements, mostly zero, with a one marking those columns to be removed.
\ingroup names
 */
static void apop_name_rm_columns(apop_name *n, int *drop){
  apop_name   *newname    = apop_name_alloc();
  size_t initial_colct = n->colct;
    for (size_t i=0; i< initial_colct; i++){
        if (drop[i]==0)
            apop_name_add(newname, n->column[i],'c');
        else
            n->colct    --;
        free(n->column[i]);
    }
    free(n->column);
    n->column = newname->column;

    //we need to free the newname struct, but leave the column intact.
    newname->column   = NULL;
    newname->colct  = 0;
    apop_name_free(newname);
}


/** Remove the columns set to one in the \c drop vector.
The returned data structure looks like it was modified in place, but the data matrix and the names are duplicated before being pared down, so if your data is taking up more than half of your memory, this may not work.

\param d the \ref apop_data structure to be pared down. 
\param drop an array of ints. If use[7]==1, then column seven will be cut from the
output. A reminder: <tt>calloc(in->size2 , sizeof(int))</tt> will fill your array with zeros on allocation, and 
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

/** \def apop_data_prune_columns(d, ...)
  Keep only the columns of a data set that you name.

\param keep_names A list of names to retain (i.e. the columns that shouldn't be pruned
out). For example, if you have run \ref apop_data_summarize, you have columns for several
statistics, but may care about only one or two; see the example.

For example:
\include eg/test_pruning.c 

\li I use case-insensitive regular expressions to find your column; see \ref apop_regex for details.
\li If your regex matches multiple columns, I'll only give you the first.
\li If I can't find a column matching one of your strings, I throw an error to the screen and continue.
\hideinitializer */
void apop_data_prune_columns_base(apop_data *d, char **colnames){
    /* In types.h, you'll find an alias that takes the input, wraps it in the cruft that is
    C's compound literal syntax, and appends a final "" to the list of strings. Here, I
    just find each element of the list, using that "" as a stopper, and then call apop_data_rm_columns.*/
    int rm_list[d->matrix->size1];
    int keep_count = 0;
    char **name_step = colnames;
    //to throw errors for typos (and slight efficiency gains), I need an array of whether
    //each input colname has been used.
    while (strlen(*name_step++))
        keep_count++;
    int used_field[keep_count];
    for (int i=0; i< keep_count; i++)
        used_field[i] = 0;

    for (int i=0; i< d->names->colct; i++){
        int keep = 0;
        for (int j=0; j<keep_count; j++)
            if (!used_field[j] && apop_regex(d->names->column[i], colnames[j])){
                keep ++;
                used_field[j]++;
                break;
            }
        rm_list[i] = !keep;
    }
    apop_data_rm_columns(d, rm_list);
    for (int j=0; j<keep_count; j++)
        if (!used_field[j])
            printf("You asked me to keep column \"%s\" but I couldn't find a match for it. Typo?", colnames[j]);
}

/** \defgroup data_set_get Set/get/point to the data element at the given point
  \{

Q: How does <tt> apop_data_set(in, row, col, data)</tt> differ from
 <tt>gsl_matrix_set(in->matrix, row, col, data)</tt>?<br>
A: It's seven characters shorter.

There are a few other differences: the \ref apop_data set has names, so we can get/set elements using those names, and it has both matrix and vector elements.


\li The versions that take a column/row name use  \ref apop_name_find
for the search; see notes there on the name matching rules.
\li For those that take a column number, column -1 is the vector element. 
\li For those that take a column name, I will search the vector last---if I don't find the name among the matrix columns, but the name matches the vector name, I return -1.
\li If you give me both a .row and a .rowname, I go with the name; similarly for .col and
.colname.
\li You can give me the name of a page, e.g.
\code
double AIC = apop_data_get(data, .rowname="AIC", .col=-1, page="Info");
\endcode
Each search for a page in the \ref apop_data set involves a regular expression search
across the names of the pages, which is expensive. If you're doing one lookup, this is not
an issue; for hundreds or thousands, you may be better off spending a line of code to get the page once:
\code
apop_data *d = apop_data_get_page(data, page="Indices");
for (int i=0; i< 1e7; i++)
    apop_data_set(d, i, 0, i);
\endcode

The \c _ptr functions return a pointer to the given cell. Those functions follow the lead of \c gsl_vector_ptr and \c gsl_matrix_ptr, and like those functions, return a pointer to the appropriate \c double.

These functions use the \ref designated syntax for inputs.

\example
\code
apop_data *d = apop_data_alloc(10, 10, 10);
apop_name_add(d->names, "Zeroth row", 'r');
apop_name_add(d->names, "Zeroth col", 'c');

apop_data_set(d, 8, 0, 27);
assert(apop_data_get(d, 8, .colname="Zeroth") == 27);
double *x = apop_data_ptr(d, .col=7, .rowname="Zeroth");
*x = 270;
assert(apop_data_get(d, 0, 7) == 270);
\endcode
*/

/** \deprecated  use \ref apop_data_ptr */
double *apop_data_ptr_ti(apop_data *in, const char* row, const int col){
  int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert(rownum != -2,  NULL, 0,'c',"Couldn't find %s amongst the row names.", row);
    return (col >= 0) ? gsl_matrix_ptr(in->matrix, rownum, col)
                      : gsl_vector_ptr(in->vector, rownum);
}

/** \deprecated  use \ref apop_data_ptr */
double *apop_data_ptr_it(apop_data *in, const size_t row, const char* col){
  int colnum =  apop_name_find(in->names, col, 'c');
    apop_assert(colnum != -2,  NULL, 0,'c',"Couldn't find %s amongst the column names.", col);
    return (colnum >= 0) ? gsl_matrix_ptr(in->matrix, row, colnum)
                         : gsl_vector_ptr(in->vector, row);
}

/** \deprecated  use \ref apop_data_ptr */
double *apop_data_ptr_tt(apop_data *in, const char *row, const char* col){
  int colnum =  apop_name_find(in->names, col, 'c');
  int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert(rownum != -2,  NULL, 0,'c',"Couldn't find %s amongst the row names.", row);
    apop_assert(colnum != -2,  NULL, 0,'c',"Couldn't find %s amongst the column names.", col);
    return (colnum >= 0) ? gsl_matrix_ptr(in->matrix, rownum, colnum)
                         : gsl_vector_ptr(in->vector, rownum);
}

/** Get a pointer to an element of an \ref apop_data set. 
 See \ref data_set_get "the set/get page" for details. */
APOP_VAR_HEAD double * apop_data_ptr(apop_data *data, const int row, const int col, const char *rowname, const char *colname, const char *page){
    apop_data * apop_varad_var(data, NULL);
    Apop_assert_s(data, "You sent me a NULL data set.");
    const int apop_varad_var(row, 0);
    const int apop_varad_var(col, 0);
    const char * apop_varad_var(rowname, NULL);
    const char * apop_varad_var(colname, NULL);
    const char * apop_varad_var(page, NULL);

    apop_data *d = page ? apop_data_get_page(data, page) : data;
    if (rowname && colname)
        return apop_data_ptr_tt(d, rowname,colname);
    if (rowname && !colname)
        return apop_data_ptr_ti(d, rowname,col);
    if (!rowname && colname)
        return apop_data_ptr_it(d, row,colname);

    //else: row number, column number
    if (col == -1){
        Apop_assert_s(data->vector, "You asked for the vector element (col=-1) but it is NULL.");
        return gsl_vector_ptr(data->vector, row);
    } else {
        Apop_assert_s(data->matrix, "You asked for the matrix element (%i, %i) but the matrix is NULL.", row, col);
        return gsl_matrix_ptr(data->matrix, row,col);
    }
APOP_VAR_ENDHEAD
return NULL;//the main function is blank.
}

/** \deprecated  use \ref apop_data_get */
double apop_data_get_ti(const apop_data *in, const char* row, const int col){
  int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert(rownum != -2,  GSL_NAN, 0,'c',"Couldn't find %s amongst the row names.", row);
    return (col >= 0) ? gsl_matrix_get(in->matrix, rownum, col)
                      : gsl_vector_get(in->vector, rownum);
}

/** \deprecated  use \ref apop_data_get */
double apop_data_get_it(const apop_data *in, const size_t row, const char* col){
  int colnum =  apop_name_find(in->names, col, 'c');
    apop_assert(colnum != -2,  GSL_NAN, 0,'c',"Couldn't find %s amongst the column names.", col);
    return (colnum >= 0) ? gsl_matrix_get(in->matrix, row, colnum)
                         : gsl_vector_get(in->vector, row);
}

/** \deprecated  use \ref apop_data_get */
double apop_data_get_tt(const apop_data *in, const char *row, const char* col){
  int colnum =  apop_name_find(in->names, col, 'c');
  int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert(colnum != -2,  GSL_NAN, 0,'c',"Couldn't find %s amongst the column names.", col);
    apop_assert(rownum != -2,  GSL_NAN, 0,'c',"Couldn't find %s amongst the row names.", row);
    return (colnum >= 0) ? gsl_matrix_get(in->matrix, rownum, colnum)
                         : gsl_vector_get(in->vector, rownum);
}

/** Returns the data element at the given point.
 See \ref data_set_get "the set/get page" for details. */
APOP_VAR_HEAD double apop_data_get(const apop_data *data, const size_t row, const int col, const char *rowname, const char *colname, const char *page){
    const apop_data * apop_varad_var(data, NULL);
    Apop_assert_s(data, "You sent me a NULL data set.");
    const size_t apop_varad_var(row, 0);
    const int apop_varad_var(col, 0);
    const char * apop_varad_var(rowname, NULL);
    const char * apop_varad_var(colname, NULL);
    const char * apop_varad_var(page, NULL);
    
    const apop_data *d = page ? apop_data_get_page(data, page) : data;
    if (rowname && colname)
        return apop_data_get_tt(d, rowname,colname);
    if (rowname && !colname)
        return apop_data_get_ti(d, rowname,col);
    if (!rowname && colname)
        return apop_data_get_it(d, row,colname);
    //else: row number, column number
    if (col>=0){
        apop_assert(d->matrix, 0, 0, 's', "You asked for the matrix element (%zu, %i) but the matrix is NULL.", row, col);
        return gsl_matrix_get(d->matrix, row, col);
    } else {
        apop_assert(d->vector, 0, 0, 's', "You asked for the vector element (col=-1) but it is NULL.");
        return gsl_vector_get(d->vector, row);
    }
APOP_VAR_ENDHEAD
return 0;//the main function is blank.
}

/** Set an element from an \ref apop_data set, using the row name but the column number 
  \deprecated Use \ref apop_data_set. 
 */
void apop_data_set_ti(apop_data *in, const char* row, const int col, const double data){
  int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert_void(rownum != -2, 0,'c',"Couldn't find %s amongst the row names. Making no changes.", row);
    return (col >= 0) ? gsl_matrix_set(in->matrix, rownum, col, data)
                      : gsl_vector_set(in->vector, rownum, data);
}

/** Set an element from an \ref apop_data set, using the column name but the row number
  \deprecated Use \ref apop_data_set.  */
void apop_data_set_it(apop_data *in, const size_t row, const char* col, const double data){
  int colnum =  apop_name_find(in->names, col, 'c');
    apop_assert_void(colnum != -2, 0,'c',"Couldn't find %s amongst the column names. Making no changes.", col);
    return (colnum >= 0) ? gsl_matrix_set(in->matrix, row, colnum, data)
                         : gsl_vector_set(in->vector, row, data);
}

/** Set an element from an \ref apop_data set, using the row and column name.  
  \deprecated Use \ref apop_data_set.  */
void apop_data_set_tt(apop_data *in, const char *row, const char* col, const double data){
  int colnum =  apop_name_find(in->names, col, 'c');
  int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert_void(colnum != -2, 0, 'c',"Couldn't find %s amongst the column names.", col);
    apop_assert_void(rownum != -2, 0, 'c',"Couldn't find %s amongst the column names.", row);
    return (colnum >= 0) ? gsl_matrix_set(in->matrix, rownum, colnum, data)
                         : gsl_vector_set(in->vector, rownum, data);
}

/**  Set a data element.
 
  This function uses the \ref designated syntax, so with names you don't have to worry
  about element ordering. But the ordering of elements may still be
  noteworthy. For compatibility with older editions of Apophenia, the order is (row, col,
  value, rowname, colname), so the following would all set row 3, column 8, of \c d to 5:
  \code
  apop_data_set(d, 3, 8, 5);
  apop_data_set(d, .row = 3, .col=8, .val=5);
  apop_data_set(d, .row = 3, .colname="Column 8", .val=5);
//but:
  apop_data_set(d, .row = 3, .colname="Column 8", 5);  //invalid---the value doesn't follow the colname.
  \endcode

  This differs somewhat from the deprecated \c _tt, \c _ti, and \c _tt functions, and to
  the some extent, the GSL.

 See \ref data_set_get "the set/get page" for details. */
APOP_VAR_HEAD void apop_data_set(apop_data *data, const size_t row, const int col, const double val, const char *colname, const char *rowname, const char *page){
    apop_data * apop_varad_var(data, NULL);
    Apop_assert_s(data, "You sent me a NULL data set.");
    const size_t apop_varad_var(row, 0);
    const int apop_varad_var(col, 0);
    const double apop_varad_var(val, 0);
    const char * apop_varad_var(rowname, NULL);
    const char * apop_varad_var(colname, NULL);
    const char * apop_varad_var(page, NULL);
    
    apop_data *d = page ? apop_data_get_page(data, page) : data;
    if (rowname && colname)
        return apop_data_set_tt(d, rowname, colname, val);
    if (rowname && !colname)
        return apop_data_set_ti(d, rowname, col, val);
    if (!rowname && colname)
        return apop_data_set_it(d, row, colname, val);
    //else: row number, column number
    if (col>=0){
        Apop_assert_s(d->matrix, "You're trying to set the matrix element (%zu, %i) but the matrix is NULL.", row, col);
        return gsl_matrix_set(d->matrix, row, col, val);
    } else {
        Apop_assert_s(d->vector, "You're trying to set a vector element (row=-1) but the vector is NULL.");
        return gsl_vector_set(d->vector, row, val);
    }
APOP_VAR_ENDHEAD
//the main function is blank.
}
/** \} //End data_set_get group */

/** A convenience function to add a named element to a data vector. For example, 
many of the testing procedures use this to easily produce a column of named parameters.

\param d    The \ref apop_data structure. Must not be \c NULL.
\param name The name to add
\param val  the value to add to the vector.

\li I use the position of the name to know where to put the value. If
there are two names in the data set, then I will put the new name in
the third name slot and the data in the third slot in the vector. If
you use this function from start to finish in building your vector,
then you'll be fine.

\li If the vector is too short (or \c NULL), I will call \ref apop_vector_realloc internally to make space in the vector.

*/
void apop_data_add_named_elmt(apop_data *d, char *name, double val){
    Apop_assert_s(d, "You sent me a NULL apop_data set. Maybe allocate with apop_data_alloc(0, 0,0) to start.");
    apop_name_add(d->names, name, 'r');
    if (!d->vector)
        d->vector = gsl_vector_alloc(1);
    if (d->vector->size < d->names->rowct)
        apop_vector_realloc(d->vector, d->names->rowct);
    gsl_vector_set(d->vector, d->names->rowct-1, val);
}

/** Add a string to the text element of an \ref apop_data set.  If you
 send me a \c NULL string, I will write the string <tt>"NaN"</tt> in the given slot.

\param in   The \ref apop_data set, that already has an allocated \c text element.
\param row  The row
\param col  The col
\param fmt The text to write.
\param ... You can use a printf-style fmt and follow it with the usual variables to fill in.
*/
void apop_text_add(apop_data *in, const size_t row, const size_t col, const char *fmt, ...){
  va_list   argp;
  apop_assert_void((in->textsize[0] >= (int)row-1) && (in->textsize[1] >= (int)col-1), 0, 'c', "You asked me to put the text "
                            " '%s' at (%zu, %zu), but the text array has size (%i, %i)\n", fmt, row, col, in->textsize[0], in->textsize[1]);
    if (!fmt){
        asprintf(&(in->text[row][col]), "NaN");
        return;
    }
	va_start(argp, fmt);
    vasprintf(&(in->text[row][col]), fmt, argp);
	va_end(argp);
}

/** This allocates an array of strings and puts it in the \c text element
  of an \ref apop_data set. 

  \param in An \ref apop_data set. It's OK to send in \c NULL, in which case an apop_data set with \c NULL \c matrix and \c vector elements is returned.
  \param row    the number of rows of text.
  \param col     the number of columns of text.
  \return       A pointer to the relevant \ref apop_data set. If the input was not \c NULL, then this is a repeat of the input pointer.
  */
apop_data * apop_text_alloc(apop_data *in, const size_t row, const size_t col){
    if (!in)
        in  = apop_data_alloc(0,0,0);
    in->text = malloc(sizeof(char**) * row);
    for (size_t i=0; i< row; i++){
        in->text[i] = malloc(sizeof(char*) * col);
        for (size_t j=0; j< col; j++)
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
    apop_name_stack(out->names, in->names, 'r', 'c');
    apop_name_stack(out->names, in->names, 'c', 'r');
    return out;
}

/** This function will resize a gsl_matrix to a new height or width.

Data in the matrix will be retained. If the new height or width is smaller than the old, then data in the later rows/columns will be cropped away (in a non--memory-leaking manner). If the new height or width is larger than the old, then new cells will be filled with garbage; it is your responsibility to zero out or otherwise fill new rows/columns before use.

<b>Warning I</b>: Using this function is basically bad form---especially when used in a <tt>for</tt> loop that adds a column each time. A large number of <tt>realloc</tt>s can take a noticeable amount of time. You are thus encouraged to make an effort to determine the size of your data beforehand.

<b>Warning II</b>: The <tt>gsl_matrix</tt> is a versatile struct that can represent submatrices and other cuts from parent data. I can't deal with those, and check for such situations beforehand. [Besides, resizing a portion of a parent matrix makes no sense.]

\param m The already-allocated matrix to resize.  If you give me \c NULL, this becomes equivalent to \c gsl_matrix_alloc
\param newheight, newwidth The height and width you'd like the matrix to be.
\return m, now resized
 */
gsl_matrix * apop_matrix_realloc(gsl_matrix *m, size_t newheight, size_t newwidth){
    if (!m)
        return (newheight && newwidth) ?  gsl_matrix_alloc(newheight, newwidth) : NULL;
    size_t i, oldoffset=0, newoffset=0, realloced = 0;
    apop_assert((m->block->data==m->data) && m->owner & (m->tda == m->size2),
                NULL, 0, 's', "I can't resize submatrices or other subviews.");
    m->block->size = newheight * newwidth;
    if (m->size2 > newwidth)
        for (i=1; i< GSL_MIN(m->size1, newheight); i++){
            oldoffset +=m->size2;
            newoffset +=newwidth;
            memmove(m->data+newoffset, m->data+oldoffset, sizeof(double)*newwidth);
        } 
    else if (m->size2 < newwidth){
        m->block->data = m->data = realloc(m->data, sizeof(double) * m->block->size);
        realloced = 1;
        int height = GSL_MIN(m->size1, newheight);
        for (i= height-1; i > 0; i--){
            newoffset +=newwidth;
            memmove(m->data+(height * newwidth) - newoffset, m->data+i*m->size2, sizeof(double)*m->size2);
        }
    }
    m->size1 = newheight;
    m->tda   =
    m->size2 = newwidth;
    if (!realloced)
        m->block->data = m->data = realloc(m->data, sizeof(double) * m->block->size);
    return m;
}

/** This function will resize a gsl_vector to a new length.

 Data in the vector will be retained. If the new height is
 smaller than the old, then data in the bottom of the vector will be
 cropped away (in a non--memory-leaking manner). If the new height is larger than the old,
 then new cells will be filled with garbage; it is your responsibility
 to zero out or otherwise fill them before use.

 <b>Warning I</b>: Using this function is basically bad form---especially
 when used in a <tt>for</tt> loop that adds an element each time. A large
 number of <tt>realloc</tt>s can take a noticeable amount of time. You are
thus encouraged to make an effort to determine the size of your data
beforehand.

 <b>Warning II</b>: The <tt>gsl_vector</tt> is a versatile struct that
 can represent subvectors, matrix columns and other cuts from parent data. I can't
 deal with those, and check for such situations beforehand. [Besides,
 resizing a portion of a parent matrix makes no sense.]

\param v The already-allocated vector to resize.  If you give me \c NULL, this is equivalent to \c gsl_vector_alloc
\param newheight The height you'd like the vector to be.
\return v, now resized
 */
gsl_vector * apop_vector_realloc(gsl_vector *v, size_t newheight){
    if (!v)
        return newheight ?  gsl_vector_alloc(newheight) : NULL;
    apop_assert((v->block->data==v->data) && v->owner & (v->stride == 1),
                NULL, 0, 's', "I can't resize subvectors or other views.");
    v->block->size = newheight;
    v->size     = newheight;
    v->block->data = 
    v->data        = realloc(v->data, sizeof(double) * v->block->size);
    return v;
}

/** It's good form to get a page from your data set by name, because you
  may not know the order for the pages, and the stepping through makes
  for dull code anyway (<tt>apop_data *page = dataset; while (page->more) page= page->more;</tt>).

  \param data The \ref apop_data set to use. No default; if \c NULL,
      gives a warning if <tt>apop_opts.verbose >=1</tt> and returns \c NULL.
  \param title The name of the page to retrieve. Default=\c "Info", which
      is the name of the page of additional estimation information returned
      by estimation routines (log likelihood, status, AIC, BIC, confidence intervals, ...).
      I use \ref apop_regex, so this is really a case-insensitive regular expression.

    \return The page whose title matches what you gave me. If I don't
    find a match, return \c NULL.

This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD apop_data * apop_data_get_page(const apop_data * data, const char *title){
    const apop_data * apop_varad_var(data, NULL);
    apop_assert(data, NULL, '1', 'c', "You requested a page from a NULL data set. Returning NULL");
    const char * apop_varad_var(title, "Info");
    return apop_data_get_page_base(data, title);
APOP_VAR_ENDHEAD
    while (data && (!data->names || !apop_regex(data->names->title, title)))
        data = data->more;
    return (apop_data *) data; //de-const.
}

/** Add a page to a \ref apop_data set. It gets a name so you can find it later.

  \param dataset The input data set, to which a page will be added.
  \param newpage The page to append
  \param title The name of the new page. Remember, this is truncated at 100 characters.

  \return The new page.  I post a warning if I am appending or appending to a \c NULL data set and  <tt>apop_opts.verbose >=1 </tt>.

  \example 
  A silly example that establishes a baseline data set, adds a page,
  modifies it, and then later retrieves it.
  \code
  apop_data *d = apop_data_alloc(10, 10, 10); //the base data set.
  apop_data *a_new_page = apop_data_add_page(d, apop_data_alloc(0,2,2), "new 2 x 2 page");
  gsl_vector_set_all(a_new_page->matrix, 3);

  //later:
  apop_data *retrieved = apop_data_get_page(d, "new"); //uses regexes, not literal match.
  apop_data_show(retrieved); //print a 2x2 grid of 3s.
  \endcode
  \param
*/
apop_data * apop_data_add_page(apop_data * dataset, apop_data *newpage, const char *title){
    apop_assert(newpage, NULL, '1', 'c', "You are adding a NULL page to a data set. Doing nothing; returning NULL.");
    snprintf(newpage->names->title, 100, "%s", title);
    apop_assert(dataset, newpage, '1', 'c', "You are adding a page to a NULL data set. Returning the new page as its own data set.");
    while (dataset->more)
        dataset = dataset->more;
    dataset->more = newpage;
    return newpage;
}

/** Remove the first page from an \ref apop_data set that matches a given name.

  \li I don't check the first page, so there's no concern that the head of your list of
  pages will move. Again, the intent of the <tt>-\>more</tt> pointer in the \ref apop_data
  set is not to fully implement a linked list, but primarily to allow you to staple auxiliary
  information to a main data set.

  \param data The input data set, to which a page will be added. No default. If \c
  NULL, I return silently if <tt> apop_opts.verbose \< 1 </tt>; print an
  error otherwise.
  \param title The name of the page to remove Default: \c "Info"
  \param free_p If \c 'y', then \ref apop_data_free the page. Default: \c 'y'.
*/
APOP_VAR_HEAD void apop_data_rm_page(apop_data * data, const char *title, const char free_p){
    apop_data *apop_varad_var(data, NULL);
    apop_assert_void(data, '1', 'c', "You are removing a page from a NULL a data set. Doing nothing.");
    const char *apop_varad_var(title, "Info");
    const char apop_varad_var(free_p, 'y');
    return apop_data_rm_page_base(data, title, free_p);
APOP_VAR_ENDHEAD
    while (data->more && !apop_regex(data->more->names->title, title))
        data = data->more;
    if (data->more){
        apop_data *tmp = data->more;
        data->more = data->more->more;
        if (free_p=='y')
            free(tmp);
    }
}
