/** \file 
The apop_data structure joins together a gsl_matrix, apop_name, and a table of strings. */
/* Copyright (c) 2006--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
//apop_gsl_error is in apop_linear_algebra.c
#define Set_gsl_handler gsl_error_handler_t *prior_handler = gsl_set_error_handler(apop_gsl_error);
#define Unset_gsl_handler gsl_set_error_handler(prior_handler);

/** Allocate a \ref apop_data structure, to be filled with data.
 
\li The typical case is  three arguments, like <tt>apop_data_alloc(2,3,4)</tt>: vector size, matrix rows, matrix cols. If the first argument is zero, you get a \c NULL vector.
\li Two arguments, <tt>apop_data_alloc(2,3)</tt>,  would allocate just a matrix, leaving the vector \c NULL.
\li One argument, <tt>apop_data_alloc(2)</tt>,  would allocate just a vector, leaving the matrix \c NULL.
\li Zero arguments, <tt>apop_data_alloc()</tt>,  will produce a basically blank set, with \c out->matrix==out->vector==NULL. 

For allocating the text part, see \ref apop_text_alloc.

The \c weights vector is set to \c NULL. If you need it, allocate it via
\code d->weights   = gsl_vector_alloc(row_ct); \endcode

\see apop_data_calloc

\return The \ref apop_data structure, allocated and ready.
\exception out->error=='a'  Allocation error. The matrix, vector, or names couldn't be <tt>malloc</tt>ed, which probably means that you requested a very large data set.

\li An \ref apop_data struct, by itself, is about 72 bytes. If I can't allocate that much memory, I return \c NULL.
                But if even this much fails, your computer may be on fire and you should go put it out. 

 \li This function uses the \ref designated syntax for inputs.
 \ingroup data_struct
*/
APOP_VAR_HEAD apop_data * apop_data_alloc(const size_t size1, const size_t size2, const int size3){
    const size_t apop_varad_var(size1, 0);
    const size_t apop_varad_var(size2, 0);
    const int apop_varad_var(size3, 0);
APOP_VAR_ENDHEAD
    size_t vsize=0, msize1=0; 
    int msize2=0;
    if (size3){
        vsize = size1;
        msize1 = size2;
        msize2 = size3;
    }
    else if (size2) {
        msize1 = size1;
        msize2 = size2;
    }
    else vsize = size1;
    apop_data *setme = malloc(sizeof(apop_data));
    Apop_stopif(!setme, return NULL, -5, "malloc failed. Probably out of memory.");
    *setme = (apop_data) { }; //init to zero/NULL.
    Set_gsl_handler
    if (msize2 > 0  && msize1 > 0){
        setme->matrix = gsl_matrix_alloc(msize1,msize2);
        Apop_stopif(!setme->matrix, setme->error='a'; return setme,
                0, "malloc failed on a %zu x %i matrix. Probably out of memory.", msize1, msize2);
    }
    if (vsize){
        setme->vector = gsl_vector_alloc(vsize);
        Apop_stopif(!setme->vector, setme->error='a'; return setme,
                0, "malloc failed on a vector of size %zu. Probably out of memory.", vsize);
    }
    Unset_gsl_handler
    setme->names = apop_name_alloc();
    Apop_stopif(!setme->names, setme->error='a'; return setme,
                0, "couldn't allocate names. Probably out of memory.");
    return setme;
}

/** Allocate a \ref apop_data structure, to be filled with data; set everything in the allocated portion to zero. See \ref apop_data_alloc for details.

\return    The \ref apop_data structure, allocated and zeroed out.
\see apop_data_alloc 
\ingroup data_struct
 This function uses the \ref designated syntax for inputs.
  */
APOP_VAR_HEAD apop_data * apop_data_calloc(const size_t size1, const size_t size2, const int size3){
    const size_t apop_varad_var(size1, 0);
    const size_t apop_varad_var(size2, 0);
    const int apop_varad_var(size3, 0);
APOP_VAR_ENDHEAD
    size_t vsize=0, msize1=0; 
    int msize2=0;
    if (size3){
        vsize = size1;
        msize1 = size2;
        msize2 = size3;
    }
    else if (size2) {
        msize1 = size1;
        msize2 = size2;
    }
    else vsize = size1;
    apop_data  *setme       = malloc(sizeof(apop_data));
    apop_assert(setme, "malloc failed. Probably out of memory.");
    *setme = (apop_data) { }; //init to zero/NULL.
    if (msize2 >0 && msize1 > 0){
        setme->matrix = gsl_matrix_calloc(msize1,msize2);
        apop_assert(setme->matrix, "malloc failed on a %zu x %i matrix. Probably out of memory.", msize1, msize2);
    }
    if (vsize){
        setme->vector = gsl_vector_calloc(vsize);
        apop_assert(setme->vector, "malloc failed on a vector of size %zu. Probably out of memory.", vsize);
    }
    else if (msize2==-1 && msize1>0){
        setme->vector = gsl_vector_calloc(msize1);
        apop_assert(setme->vector, "malloc failed on a vector of size %zu. Probably out of memory.", msize1);
    }
    setme->names = apop_name_alloc();
    return setme;
}

/** Wrap an \ref apop_data structure around an existing \c gsl_matrix.
 The matrix is not copied, but is pointed to by the new \ref apop_data struct.

\param m  The existing matrix you'd like to turn into an \ref apop_data structure.
\return   The \ref apop_data structure whose \c matrix pointer points to the input matrix. The rest of the struct is basically blank.
\li If you give me a \c NULL matrix, I return a blank \ref apop_data set, equivalent to <tt>apop_data_alloc()</tt>, and print 
             a warning if <tt>apop_opts.verbosity >=1</tt>
  */
apop_data * apop_matrix_to_data(gsl_matrix *m){
    Apop_assert_c(m, apop_data_alloc(), 1, "Converting a NULL matrix to a blank apop_data structure.");
    apop_data *setme = apop_data_alloc();
    setme->matrix = m;
    return setme;
}

/** Wrap an \ref apop_data structure around an existing \c gsl_vector.
 The vector is not copied, but is pointed to by the new \ref apop_data struct.

\param  v   The data vector
\return     an allocated, ready-to-use \ref apop_data structure.
\li If you give me a \c NULL vector, I return a blank \ref apop_data set, equivalent to <tt>apop_data_alloc()</tt>, and print 
             a warning if <tt>apop_opts.verbosity >=1</tt>
*/
apop_data * apop_vector_to_data(gsl_vector *v){
    Apop_assert_c(v, apop_data_alloc(), 1, "Converting a NULL vector to a blank apop_data structure.");
    apop_data  *setme = apop_data_alloc();
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

/** Free the elements of the given \ref apop_data set and then the \ref apop_data set
  itself. Intended to be used by \ref apop_data_free, a macro that calls this to free
  elements, then sets the value to \c NULL.

\li \ref apop_data_free is a macro that calls this function and, on success, sets the input pointer to \c NULL. 
For typical cases, that's slightly more useful than this function.

\exception freeme.error='c' Circular linking is against the rules. If <tt>freeme->more == freeme</tt>, then 
I set <tt>freeme.error='c'</tt> and return. If you send in a structure like A -> B ->
B, then both data sets A and B will be marked.

\return \c 0 on OK, \c 'c' on error.
  */
char apop_data_free_base(apop_data *freeme){
    if (!freeme) return 0;
    if (freeme->more){
        Apop_stopif(freeme == freeme->more, freeme->error='c'; return 'c',
                            1, "the ->more element of this data set equals the data set itself. "
                               "This is not healthy. Not freeing; marking your data set with error='c'.");
        if (apop_data_free_base(freeme->more)) 
            Apop_stopif(freeme->more->error == 'c', freeme->error='c'; return 'c', 
                                1, "Propogating error code to parent data set");
    } 
    if (freeme->vector)  
        gsl_vector_free(freeme->vector);
    if (freeme->matrix)  
        gsl_matrix_free(freeme->matrix); 
    if (freeme->weights)
        gsl_vector_free(freeme->weights);
    apop_name_free(freeme->names);
    apop_text_free(freeme->text, freeme->textsize[0] , freeme->textsize[1]);
    free(freeme);
    return 0;
}

/** Copy one \ref apop_data structure to another. That is, all data on the first page is duplicated. [To do multiple pages, call this via a \c for loop over the data set's pages.]

  This function does not allocate the output structure or the vector, matrix, text, or weights elements---I assume you have already done this and got the dimensions right. I will assert that there is at least enough room in the destination for your data, and fail if the copy would write more elements than there are bins.


  \li If you want space allocated, use \ref apop_data_copy.
  \li I don't follow the \c more pointer, though \ref apop_data_copy does.
  \li You can use the subsetting macros, \ref Apop_data_row or \ref Apop_data_rows, to copy within a data set:

\code
//Copy the contents of row i of mydata to row j.
Apop_data_row(mydata, i, fromrow);
Apop_data_row(mydata, j, torow);
apop_data_memcpy(torow, fromrow);
\endcode
 
  \param out   a structure that this function will fill. Must be preallocated with the appropriate sizes.
  \param in    the input data

\exception out.error='d'  Dimension error; couldn't copy.
\exception out.error='p'  Part missing; e.g., in->matrix exists but out->matrix doesn't; couldn't copy.

 \ingroup data_struct
  */
void apop_data_memcpy(apop_data *out, const apop_data *in){
    Apop_assert_c(out, , 1, "you are copying to a NULL matrix. Do you mean to use apop_data_copy instead?");
    if (in->matrix){
        Apop_stopif(!out->matrix, out->error='p'; return, 1, "in->matrix exists but out->matrix does not.");
        Apop_stopif(in->matrix->size1 != out->matrix->size1 || in->matrix->size2 != out->matrix->size2, 
                out->error='d'; return,
                1, "you're trying to copy a (%zu X %zu) into a (%zu X %zu) matrix.", 
                        in->matrix->size1, in->matrix->size2, out->matrix->size1, out->matrix->size2);
        gsl_matrix_memcpy(out->matrix, in->matrix);
    }
    if (in->vector){
        Apop_stopif(!out->vector, out->error='p'; return, 1, "in->vector exists but out->vector does not.");
        Apop_stopif(in->vector->size != out->vector->size,
                out->error='d'; return,
                1, "You're trying to copy a %zu-elmt "
                        "vector into a %zu-elmt vector.", in->vector->size, out->vector->size);
        gsl_vector_memcpy(out->vector, in->vector);
    }
    if (in->weights){
        Apop_stopif(!out->weights, out->error='p'; return, 1, "in->weights exists but out->weights does not.");
        Apop_stopif(in->weights->size != out->weights->size,
                    out->error='d'; return,
                    1, "Weight vector sizes don't match: "
                    "you're trying to copy a %zu-elmt vector into a %zu-elmt vector.", 
                                 in->weights->size, out->weights->size);
        gsl_vector_memcpy(out->weights, in->weights);
    }
    if (in->names){
        apop_name_free(out->names);
        out->names = apop_name_alloc();
        apop_name_stack(out->names, in->names, 'v');
        apop_name_stack(out->names, in->names, 'r');
        apop_name_stack(out->names, in->names, 'c');
        apop_name_stack(out->names, in->names, 't');
    }
    out->textsize[0] = in->textsize[0]; 
    out->textsize[1] = in->textsize[1]; 
    if (in->textsize[0] && in->textsize[1]){
        Apop_stopif(out->textsize[0] < in->textsize[0] || out->textsize[1] < in->textsize[1],
                    out->error='d'; return,
                    1, "I am trying to copy a grid of (%zu, %zu) text elements into a grid of (%zu, %zu), "
                    "and that won't work. Please use apop_text_alloc to reallocate the right amount of data, "
                    "or use apop_data_copy for automatic allocation.",
                    in->textsize[0] , in->textsize[1] , out->textsize[0] , out->textsize[1]);
        for (size_t i=0; i< in->textsize[0]; i++)
            for(size_t j=0; j < in->textsize[1]; j ++)
				apop_text_add(out, i, j, "%s", in->text[i][j]);
    }
}

/** Copy one \ref apop_data structure to another. That is, all data is duplicated.

  Basically a front-end for \ref apop_data_memcpy for those who prefer this sort of syntax. 

  Unlike \ref apop_data_memcpy, I do follow the \c more pointer.
 
  \param in    the input data
  \return       a structure that this function will allocate and fill. If input is NULL, then this will be NULL.

 \ingroup data_struct

\exception out.error='a'  Allocation error.
\exception out.error='c'  Cyclic link: <tt>D->more == D</tt> (may be later in the chain, e.g., <tt>D->more->more = D->more</tt>) You'll have only a partial copy.
\exception out.error='d'  Dimension error; should never happen.
\exception out.error='p'  Missing part error; should never happen.
\li If the input data set has an error, then I will copy it anyway, including the error flag (which might be overwritten). I print a warning if the verbosity level is <tt>>=1</tt>.
  */
apop_data *apop_data_copy(const apop_data *in){
    if (!in) return NULL;
    apop_data *out = apop_data_alloc();
    Apop_stopif(out->error, return out, 0, "Allocation error.");
    if (in->error){
        Apop_notify(1, "the data set to be copied has an error flag of %c. Copying it.", in->error);
        out->error = in->error;
    }
    if (in->more){
        Apop_stopif(in == in->more, out->error='c'; return out,
                0, "the ->more element of this data set equals the "
                                        "data set itself. This is not healthy. Made a partial copy and set out.error='c'.");
        out->more = apop_data_copy(in->more);
        Apop_stopif(out->more->error, out->error=out->more->error; return out,
                0, "propagating an error in the ->more element to the parent apop_data set. Only a partial copy made.");
    }
    if (in->vector){
        out->vector = gsl_vector_alloc(in->vector->size);
        Apop_stopif(!out->vector, out->error='a'; return out, 0, "Allocation error on vector of size %zu.", in->vector->size);
    }
    if (in->matrix){  
        out->matrix = gsl_matrix_alloc(in->matrix->size1, in->matrix->size2);
        Apop_stopif(!out->matrix, out->error='a'; return out, 0, "Allocation error on matrix "
                    "of size %zu X %zu.", in->matrix->size1, in->matrix->size2);
    }
    if (in->weights){
        out->weights = gsl_vector_alloc(in->weights->size);
        Apop_stopif(!out->weights, out->error='a'; return out, 0, "Allocation error on weights vector of size %zu.", in->weights->size);
    }
    if (in->textsize[0] && in->textsize[1]){
        apop_text_alloc(out, in->textsize[0], in->textsize[1]);
        Apop_stopif(out->error, return out, 0, "Allocation error on text grid of size %zu X %zu.", in->textsize[0], in->textsize[1]);
    }
    apop_data_memcpy(out, in);
    return out;
}

/** Put the first data set either on top of or to the left of the second data set.

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
\exception out->error=='a' Allocation error.
\exception out->error=='d'  Dimension error; couldn't copy.

\li If m1 or m2 are NULL, this returns a copy of the other element, and if
both are NULL, you get NULL back (except if \c m2 is \c NULL and \c inplace is \c 'y', where you'll get the original \c m1 pointer back)
\li Text is handled as you'd expect: If 'r', one set of text is stacked on top of the other [number of columns must match]; if 'c', one set of text is set next to the other [number of rows must match].
\li \c more is ignored.
\li If stacking rows on rows, the output vector is the input
vectors stacked accordingly. If stacking columns by columns, the output
vector is just a copy of the vector of m1 and m2->vector doesn't appear in the
output at all.  
\li The same rules for dealing with the vector(s) hold for the vector(s) of weights.
\li Names are a copy of the names for \c m1, with the names for \c m2 appended to the row or column list, as appropriate.
\li This function uses the \ref designated syntax for inputs.
\ingroup data_struct
*/
APOP_VAR_HEAD apop_data *apop_data_stack(apop_data *m1, apop_data * m2, char posn, char inplace){
    apop_data * apop_varad_var(m1, NULL)
    apop_data * apop_varad_var(m2, NULL)
    char apop_varad_var(posn, 'r')
    Apop_assert_c((posn == 'r' || posn == 'c'), NULL, 0, "Valid positions are 'r' or 'c'"
                                                         " you gave me '%c'. Returning NULL.", posn);
    char apop_varad_var(inplace, 'n')
    inplace = (inplace == 'i' || inplace == 'y' || inplace == 1 || inplace == 'I' || inplace == 'Y') ? 1 : 0;
APOP_VAR_ENDHEAD
    if (!m1) return apop_data_copy(m2);
    if (!m2) return inplace ? m1 : apop_data_copy(m1);
    apop_data *out = NULL;
    if (inplace)
        out = m1;
    else {
        apop_data *m = m1->more; //not following the more pointer.
        m1->more =NULL;
        out = apop_data_copy(m1);
        Apop_stopif(out->error, return out, 0, "initial copy failed; leaving.");
        m1->more = m;
    }
    Get_vmsizes(m1); //original sizes of vsize, msize1, msize2.
    out->matrix = apop_matrix_stack(out->matrix, m2->matrix, posn, .inplace='y');
    if (posn == 'r'){
        out->vector  = apop_vector_stack(out->vector, m2->vector, .inplace='y');
        out->weights = apop_vector_stack(out->weights, m2->weights, .inplace='y');
    } 

    if (m2->text){ //we've already copied m1->text, if any, so if m2->text is NULL, we're done.
        if (posn=='r'){
            apop_assert(!out->text || m2->textsize[1]==out->textsize[1], 
                            "The first data set has %zu columns of text and the second has %zu columns. "
                            "I can't stack that.", out->textsize[1], m2->textsize[1]);
            int basetextsize = out->textsize[0];
            apop_text_alloc(out, basetextsize+m2->textsize[0], out->textsize[1]);
            Apop_stopif(out->error, return out, 0, "Allocation error.");
            for(int i=0; i< m2->textsize[0]; i++)
                for(int j=0; j< m2->textsize[1]; j++)
                    apop_text_add(out, i+basetextsize, j, m2->text[i][j]);
        } else {
            apop_assert(!out->text || m2->textsize[0]==out->textsize[0], 
                            "The first data set has %zu rows of text and the second has %zu rows. "
                            "I can't stack that.", out->textsize[0], m2->textsize[0]);
            int basetextsize = out->textsize[1];
            apop_text_alloc(out, out->textsize[0], basetextsize+m2->textsize[1]);
            Apop_stopif(out->error, return out, 0, "Allocation error.");
            for(int i=0; i< m2->textsize[0]; i++)
                for(int j=0; j< m2->textsize[1]; j++)
                    apop_text_add(out, i, j+basetextsize, m2->text[i][j]);
            apop_name_stack(out->names, m2->names, 't');
        }
    }
    if ((posn=='r' && m2->names->rowct) || (posn=='c' && m2->names->colct)){
        int min = posn =='r' ? m1->names->rowct : m1->names->colct;
        int max = posn =='r' ? GSL_MAX(vsize, msize1) : msize2;
        for (int k = min; k < max; k++)          //pad so the name stacking is aligned (if needed)
            apop_name_add(out->names, "", posn); 
        apop_name_stack(out->names, m2->names, posn);
    }
    return out;
}

/** Split one input \ref apop_data structure into two.

 For the opposite operation, see \ref apop_data_stack.
 
\param in  The \ref apop_data structure to split 
\param splitpoint The index of what will be the first row/column of the second data set.  E.g., if this is -1 and \c r_or_c=='c', then the whole data set will be in the second data set; if this is the length of the matrix then the whole data set will be in the first data set. Another way to put it is that \c splitpoint will equal the number of rows/columns in the first matrix (unless it is -1, in which case the first matrix will have zero rows, or it is greater than the matrix's size, in which case it will have as many rows as the original).  
\param r_or_c If this is 'r' or 'R', then put some rows in the first data set and some in the second; of 'c' or 'C', split columns into first and second data sets.

 \return An array of two \ref apop_data sets. If one is empty then a
 \c NULL pointer will be returned in that position. For example, for a data set of 50 rows, <tt>apop_data **out = apop_data_split(data, 100, 'r')</tt> sets <tt>out[0] = apop_data_copy(data)</tt> and <tt>out[1] = NULL</tt>.

 \li When splitting at a row, the text is also split.
 \li \c more pointer is ignored.
 \li The <tt>apop_data->vector</tt> is taken to be the -1st element of the matrix.  
 \li Weights will be preserved. If splitting by rows, then the top and bottom parts of the weights vector will be assigned to the top and bottom parts of the main data set. If splitting by columns, identical copies of the weights vector will be assigned to both parts.
 \li Data is copied, so you may want to call <tt>apop_data_free(in)</tt> after this.
 */
apop_data ** apop_data_split(apop_data *in, int splitpoint, char r_or_c){
    //A long, dull series of contingencies. Bonus: a reasonable use of goto.
    apop_data   **out   = malloc(2*sizeof(apop_data *));
    out[0] = out[1] = NULL;
    Apop_assert_c(in, out, 1, "input was NULL; output will be an array of two NULLs.");
    gsl_vector v1, v2, w1, w2;
    gsl_matrix m1, m2;
    int set_v1 = 1, set_v2 = 1,
        set_m1 = 1, set_m2 = 1,
        set_w1 = 1, set_w2 = 1,
        namev0 = 0, namev1 = 0,
        namer0 = 0, namer1 = 0,
        namec0 = 0, namec1 = 0,
        namersplit = -1, namecsplit = -1;
     if (r_or_c == 'r' || r_or_c == 'R') {
        if (splitpoint <=0)
            out[1]  = apop_data_copy(in);
        else if (in->matrix && splitpoint >= in->matrix->size1)
            out[0]  = apop_data_copy(in);
        else {
            namev0  =
            namev1  = 
            namec0  =
            namec1  = 1;
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
            namersplit=splitpoint;
            goto allocation;
        }
    } else if (r_or_c == 'c' || r_or_c == 'C') {
        if (in->weights){
            w1      = gsl_vector_subvector(in->weights, 0, in->weights->size).vector;
            w2      = gsl_vector_subvector(in->weights, 0, in->weights->size).vector;
        } else 
            set_w1 = 
            set_w2 = 0;
        namer0 = 1;
        namer1 = 1;

        if (splitpoint <= -1)
            out[1]  = apop_data_copy(in);
        else if (in->matrix && splitpoint >= in->matrix->size2)
            out[0]  = apop_data_copy(in);
        else if (splitpoint == 0){
            if (in->vector){
                v1      = gsl_vector_subvector(in->vector, 0, in->vector->size).vector;
                namev0  = 1;
            } else 
                set_v1 = 0;
            set_v2  = 0;
            set_m1  = 0;
            if (in->matrix){
                m2      = gsl_matrix_submatrix (in->matrix, 0, 0, 
                                    in->matrix->size1,  in->matrix->size2).matrix;
                namec1  = 1;
            } else 
                set_m2 = 0;
            goto allocation;
        } else if (splitpoint > 0 && in->matrix && splitpoint < in->matrix->size2){
            if (in->vector){
                v1      = gsl_vector_subvector(in->vector, 0, in->vector->size).vector;
                namev0  = 1;
            } else 
                set_v1 = 0;
            set_v2  = 0;
            if (in->matrix){
                m1      = gsl_matrix_submatrix (in->matrix, 0, 0, in->matrix->size1, splitpoint).matrix;
                m2      = gsl_matrix_submatrix (in->matrix, 0, splitpoint, 
                                    in->matrix->size1,  in->matrix->size2-splitpoint).matrix;
                namecsplit = splitpoint;
            } else
                set_m1  = 
                set_m2  = 0;
            goto allocation;
        } else { //splitpoint >= in->matrix->size2
            if (in->vector){
                v1      = gsl_vector_subvector(in->vector, 0, in->vector->size).vector;
                namev0  = 1;
            } else 
                set_v1 = 0;
            set_v2  = 0;
            if (in->matrix){
                m1      = gsl_matrix_submatrix (in->matrix, 0, 0, 
                            in->matrix->size1, in->matrix->size2).matrix;
                namec0 = 1;
            }
            else set_m1 = 0;
            set_m2  = 0;
            goto allocation;
        }
    } else Apop_assert_c(0, out, 0, "Please set r_or_c == 'r' or == 'c'. Returning two NULLs.");
    return out;

allocation:
    out[0]  = apop_data_alloc();
    out[1]  = apop_data_alloc();
    if (set_v1) out[0]->vector  = apop_vector_copy(&v1);
    if (set_v2) out[1]->vector  = apop_vector_copy(&v2);
    if (set_m1) out[0]->matrix  = apop_matrix_copy(&m1);
    if (set_m2) out[1]->matrix  = apop_matrix_copy(&m2);
    if (set_w1) out[0]->weights  = apop_vector_copy(&w1);
    if (set_w2) out[1]->weights  = apop_vector_copy(&w2);
    if (namev0 && out[0]) apop_name_stack(out[0]->names, in->names, 'v');
    if (namev1 && out[1]) apop_name_stack(out[1]->names, in->names, 'v');
    if (namersplit >=0)
        for (int k=0; k< in->names->rowct; k++){
            int which = (k >= namersplit);
            assert(out[which]);
            apop_name_add(out[which]->names, in->names->row[k], 'r');
        }
    else {
        if (namer0 && out[0]) apop_name_stack(out[0]->names, in->names, 'r');
        if (namer1 && out[1]) apop_name_stack(out[1]->names, in->names, 'r');
    }
    if (namecsplit >=0)
        for (int k=0; k< in->names->colct; k++){
            int which = (k >= namecsplit);
            assert(out[which]);
            apop_name_add(out[which]->names, in->names->column[k], 'c');
        }
    else {
        if (namec0 && out[0]) apop_name_stack(out[0]->names, in->names, 'c');
        if (namec1 && out[1]) apop_name_stack(out[1]->names, in->names, 'c');
    }
    //finally, the text [split by rows only]
    if (r_or_c=='r' && in->textsize[0] && in->textsize[1]){
        apop_name_stack(out[1]->names, in->names, 't');
        apop_text_alloc(out[0], splitpoint, in->textsize[1]);
        Apop_stopif(out[0]->error, return out, 0, "Allocation error.");
        if (in->textsize[0] > splitpoint){
            apop_name_stack(out[0]->names, in->names, 't');
            apop_text_alloc(out[1], in->textsize[0]-splitpoint, in->textsize[1]);
            Apop_stopif(out[1]->error, return out, 0, "Allocation error.");
        }
        for (int i=0; i< in->textsize[0]; i++)
            for (int j=0; j< in->textsize[1]; j++){
                int whichtext = (i >= splitpoint);
                int row = whichtext ? i - splitpoint : i;
                asprintf(&(out[whichtext]->text[row][j]), "%s", in->text[i][j]);
            }
    }
    return out;
}


/** Remove the columns set to one in the \c drop vector.
\param n the \ref apop_name structure to be pared down
\param drop  a vector with n->colct elements, mostly zero, with a one marking those columns to be removed.
\ingroup names
 */
static void apop_name_rm_columns(apop_name *n, int *drop){
    apop_name *newname = apop_name_alloc();
    size_t initial_colct = n->colct;
    for (size_t i=0; i< initial_colct; i++){
        if (drop[i]==0) apop_name_add(newname, n->column[i],'c');
        else            n->colct--;
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
    gsl_matrix *freeme = d->matrix;
    d->matrix = apop_matrix_rm_columns(d->matrix, drop);
    gsl_matrix_free(freeme); 
    apop_name_rm_columns(d->names, drop);
}

/** \def apop_data_prune_columns(in, ...)
  Keep only the columns of a data set that you name.

\param in The data set to prune.
\param ... A list of names to retain (i.e. the columns that shouldn't be pruned
out). For example, if you have run \ref apop_data_summarize, you have columns for several
statistics, but may care about only one or two; see the example.

For example:
\include test_pruning.c 

\li I use case-insensitive regular expressions to find your column; see \ref apop_regex for details.
\li If your regex matches multiple columns, I'll only give you the first.
\li If I can't find a column matching one of your strings, I throw an error to the screen and continue.
\li This is a macro calling \ref apop_data_prune_columns_base. It packages your list of
columns into a list of strings, adds a \c NULL string at the end, and calls that function.
\hideinitializer */
 
/** Keep only the columns of a data set that you name.
  This is the function called internally by the \ref apop_data_prune_columns macro. In
  most cases, you'll want to use that macro. An example of the two uses demonstrating the
  difference:

  \code 
    apop_data_prune_columns(d, "mean", "median");

    char *list[] = {"mean", "median", NULL};
    apop_data_prune_columns_base(d, list);
  \endcode

\param d The data set to prune.
\param colnames A null-terminated list of names to retain (i.e. the columns that shouldn't be pruned
out). 
  */
void apop_data_prune_columns_base(apop_data *d, char **colnames){
    /* In types.h, you'll find an alias that takes the input, wraps it in the cruft that is
    C's compound literal syntax, and appends a final "" to the list of strings. Here, I
    find each element of the list, using that "" as a stopper, and then call apop_data_rm_columns.*/
    Apop_assert_c(d, , 1, "You're asking me to prune a NULL data set; returning NULL.");
    int rm_list[d->matrix->size1];
    int keep_count = 0;
    char **name_step = colnames;
    //to throw errors for typos (and slight efficiency gains), I need an array of whether
    //each input colname has been used.
    while (*name_step++)
        keep_count++;
    int used_field[keep_count];
    memset(used_field, 0, keep_count*sizeof(int));

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
        apop_assert_c(used_field[j], , 1, "You asked me to keep column \"%s\" but I couldn't find a match for it. Typo?", colnames[j]);
}

/** \defgroup data_set_get Set/get/point to the data element at the given point
  \{
First, some examples:

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

Q: How does <tt> apop_data_set(in, row, col, data)</tt> differ from
 <tt>gsl_matrix_set(in->matrix, row, col, data)</tt>?<br>
A: It's seven characters shorter.

OK, there are a few other differences: the \ref apop_data set has names, so we can get/set elements using those names, and it has both matrix and vector elements.

\li The versions that take a column/row name use  \ref apop_name_find
for the search; see notes there on the name matching rules.
\li For those that take a column number, column -1 is the vector element. 
\li For those that take a column name, I will search the vector last---if I don't find the name among the matrix columns, but the name matches the vector name, I return -1.
\li If you give me both a .row and a .rowname, I go with the name; similarly for .col and
.colname.
\li You can give me the name of a page, e.g.
\code
double AIC = apop_data_get(data, .rowname="AIC", .col=-1, .page="Info");
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
*/

/* \deprecated  use \ref apop_data_ptr */
double *apop_data_ptr_ti(apop_data *in, const char* row, const int col){
    int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert_c(rownum != -2,  NULL, 0,"Couldn't find '%s' amongst the row names.", row);
    return (col >= 0) ? gsl_matrix_ptr(in->matrix, rownum, col)
                      : gsl_vector_ptr(in->vector, rownum);
}

/* \deprecated  use \ref apop_data_ptr */
double *apop_data_ptr_it(apop_data *in, const size_t row, const char* col){
    int colnum =  apop_name_find(in->names, col, 'c');
    apop_assert_c(colnum != -2,  NULL, 0,"Couldn't find '%s' amongst the column names.", col);
    return (colnum >= 0) ? gsl_matrix_ptr(in->matrix, row, colnum)
                         : gsl_vector_ptr(in->vector, row);
}

/* \deprecated  use \ref apop_data_ptr */
double *apop_data_ptr_tt(apop_data *in, const char *row, const char* col){
    int colnum =  apop_name_find(in->names, col, 'c');
    int rownum =  apop_name_find(in->names, row, 'r');
    apop_assert_c(rownum != -2,  NULL, 0,"Couldn't find '%s' amongst the row names.", row);
    apop_assert_c(colnum != -2,  NULL, 0,"Couldn't find '%s' amongst the column names.", col);
    return (colnum >= 0) ? gsl_matrix_ptr(in->matrix, rownum, colnum)
                         : gsl_vector_ptr(in->vector, rownum);
}

/** Get a pointer to an element of an \ref apop_data set. 

\li If a \c NULL vector or matrix (as the case may be), stop (unless <tt>apop_opts.stop_on_warning='n'</tt>, then return \c NULL).
\li If the row/column you requested is outside the bounds of the matrix (or the name isn't found), always return \c NULL.
\li See \ref data_set_get "the set/get page" for details. 
*/
APOP_VAR_HEAD double * apop_data_ptr(apop_data *data, const int row, const int col, const char *rowname, const char *colname, const char *page){
    apop_data * apop_varad_var(data, NULL);
    Apop_assert(data, "You sent me a NULL data set.");
    const int apop_varad_var(row, 0);
    const int apop_varad_var(col, 0);
    const char * apop_varad_var(rowname, NULL);
    const char * apop_varad_var(colname, NULL);
    const char * apop_varad_var(page, NULL);

    apop_data *d;
    if (page){
        d = apop_data_get_page(data, page);
        Apop_assert(d, "I couldn't find a page with label '%s'.", page);
    } else d = data;
    if (rowname && colname)
        return apop_data_ptr_tt(d, rowname,colname);
    if (rowname && !colname)
        return apop_data_ptr_ti(d, rowname,col);
    if (!rowname && colname)
        return apop_data_ptr_it(d, row,colname);

    //else: row number, column number
    if (col == -1){
        Apop_assert(data->vector, "You asked for the vector element (col=-1) but it is NULL.");
        return gsl_vector_ptr(data->vector, row);
    } else {
        Apop_assert(data->matrix, "You asked for the matrix element (%i, %i) but the matrix is NULL.", row, col);
        return gsl_matrix_ptr(data->matrix, row,col);
    }
APOP_VAR_ENDHEAD
    return NULL;//the main function is blank.
}

/* \deprecated  use \ref apop_data_get */
double apop_data_get_ti(const apop_data *in, const char* row, const int col){
    int rownum =  apop_name_find(in->names, row, 'r');
    Apop_assert_c(rownum != -2,  GSL_NAN, 0,"Couldn't find '%s' amongst the row names.", row);
    if (col >= 0){
        Apop_assert_nan(in->matrix, "You asked me to get the (%i, %i) element of a NULL matrix.", rownum, col);
        return gsl_matrix_get(in->matrix, rownum, col);
    } else {
        Apop_assert_nan(in->vector, "You asked me to get the %ith element of a NULL vector.", rownum);
        return gsl_vector_get(in->vector, rownum);
    }
}

/* \deprecated  use \ref apop_data_get */
double apop_data_get_it(const apop_data *in, const size_t row, const char* col){
    int colnum = apop_name_find(in->names, col, 'c');
    Apop_assert_c(colnum != -2,  GSL_NAN, 0,"Couldn't find '%s' amongst the column names.", col);
    if (colnum >= 0){
        Apop_assert_nan(in->matrix, "You asked me to get the (%zu, %i) element of a NULL matrix.", row, colnum);
        return gsl_matrix_get(in->matrix, row, colnum);
    } else {
        Apop_assert_nan(in->vector, "You asked me to get the %zuth element of a NULL vector.", row);
        return gsl_vector_get(in->vector, row);
    }
}

/* \deprecated  use \ref apop_data_get */
double apop_data_get_tt(const apop_data *in, const char *row, const char* col){
    int colnum =  apop_name_find(in->names, col, 'c');
    int rownum =  apop_name_find(in->names, row, 'r');
    Apop_assert_c(colnum != -2,  GSL_NAN, 0,"Couldn't find '%s' amongst the column names.", col);
    Apop_assert_c(rownum != -2,  GSL_NAN, 0,"Couldn't find '%s' amongst the row names.", row);
    if (colnum >= 0){
        Apop_assert_nan(in->matrix, "You asked me to get the (%i, %i) element of a NULL matrix.", rownum, colnum);
        return gsl_matrix_get(in->matrix, rownum, colnum);
    } else {
        Apop_assert_nan(in->vector, "You asked me to get the %ith element of a NULL vector.", rownum);
        return gsl_vector_get(in->vector, rownum);
    }
}

/** Returns the data element at the given point.
 
  In case of error (probably that you asked for a data point out of bounds), returns \c GSL_NAN.
 See \ref data_set_get "the set/get page" for details. */
APOP_VAR_HEAD double apop_data_get(const apop_data *data, const size_t row, const int col, const char *rowname, const char *colname, const char *page){
    const apop_data * apop_varad_var(data, NULL);
    Apop_assert(data, "You sent me a NULL data set.");
    const size_t apop_varad_var(row, 0);
    const int apop_varad_var(col, 0);
    const char * apop_varad_var(rowname, NULL);
    const char * apop_varad_var(colname, NULL);
    const char * apop_varad_var(page, NULL);
    
    const apop_data *d;
    if (page){
        d = apop_data_get_page(data, page);
        Apop_assert_nan(d, "I couldn't find a page with label '%s'.", page);
    } else d = data;
    if (rowname && colname)
        return apop_data_get_tt(d, rowname,colname);
    if (rowname && !colname)
        return apop_data_get_ti(d, rowname,col);
    if (!rowname && colname)
        return apop_data_get_it(d, row,colname);
    //else: row number, column number
    if (col>=0){
        Apop_assert_nan(d->matrix, "You asked for the matrix element (%zu, %i) but the matrix is NULL.", row, col);
        return gsl_matrix_get(d->matrix, row, col);
    } else {
        Apop_assert_nan(d->vector, "You asked for the vector element (col=-1) but it is NULL.");
        return gsl_vector_get(d->vector, row);
    }
APOP_VAR_ENDHEAD
    return 0;//the main function is blank.
}

/* The only hint the GSL gives that something failed is that the error-handler is called.
   The error handling function won't let you set an output to the function. So all we
   can do is use a global variable.
*/

static threadlocal int error_for_set; //see apop_internal.h

void apop_gsl_error_for_set(const char *reason, const char *file, int line, int gsl_errno){
    Apop_notify(1, "%s: %s", file, reason);
    Apop_maybe_abort(1);
    error_for_set = -1;
}

/* Set an element from an \ref apop_data set, using the row name but the column number 
  \deprecated Use \ref apop_data_set. 
 */
int apop_data_set_ti(apop_data *in, const char* row, const int col, const double data){
    Set_gsl_handler
    int rownum =  apop_name_find(in->names, row, 'r');
    Apop_assert_c(rownum != -2, -1, 0, "Couldn't find '%s' amongst the row names. Making no changes.", row);
    if (col >= 0) gsl_matrix_set(in->matrix, rownum, col, data);
    else          gsl_vector_set(in->vector, rownum, data);
    Unset_gsl_handler
    return error_for_set;
}

/* Set an element from an \ref apop_data set, using the column name but the row number
  \deprecated Use \ref apop_data_set.  */
int apop_data_set_it(apop_data *in, const size_t row, const char* col, const double data){
    Set_gsl_handler
    int colnum =  apop_name_find(in->names, col, 'c');
    Apop_assert_c(colnum != -2, -1, 0, "Couldn't find '%s' amongst the column names. Making no changes.", col);
    if (colnum >= 0)  gsl_matrix_set(in->matrix, row, colnum, data);
    else              gsl_vector_set(in->vector, row, data);
    Unset_gsl_handler
    return error_for_set;
}

/* Set an element from an \ref apop_data set, using the row and column name.  
  \deprecated Use \ref apop_data_set.  */
int apop_data_set_tt(apop_data *in, const char *row, const char* col, const double data){
    Set_gsl_handler
    int colnum =  apop_name_find(in->names, col, 'c');
    int rownum =  apop_name_find(in->names, row, 'r');
    Apop_assert_c(colnum != -2, -1, 0, "Couldn't find '%s' amongst the column names. Making no changes.", col);
    Apop_assert_c(rownum != -2, -1, 0, "Couldn't find '%s' amongst the column names. Making no changes.", row);
    if (colnum >= 0) gsl_matrix_set(in->matrix, rownum, colnum, data);
    else             gsl_vector_set(in->vector, rownum, data);
    Unset_gsl_handler
    return error_for_set;
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
  some extent, the GSL.

  \return 0=OK, -1=error (couldn't find row/column name, or you asked for a location outside the vector/matrix bounds).

\li  The error codes for out-of-bounds errors are thread-safe iff you are have a
C11-compliant compiler (thanks to the \c _Thread_local keyword) or GCC (thanks to its \c __thread
extension).

 See \ref data_set_get "the set/get page" for details. */
APOP_VAR_HEAD int apop_data_set(apop_data *data, const size_t row, const int col, const double val, const char *colname, const char *rowname, const char *page){
    apop_data * apop_varad_var(data, NULL);
    Apop_assert_negone(data, "You sent me a NULL data set.");
    const size_t apop_varad_var(row, 0);
    const int apop_varad_var(col, 0);
    const double apop_varad_var(val, 0);
    const char * apop_varad_var(rowname, NULL);
    const char * apop_varad_var(colname, NULL);
    const char * apop_varad_var(page, NULL);
    
    apop_data *d;
    if (page){
        d = apop_data_get_page((apop_data*)data, page);
        Apop_assert_negone(d, "I couldn't find a page with label '%s'.", page);
    } else d = data;
    if (rowname && colname)  return apop_data_set_tt(d, rowname, colname, val);
    if (rowname && !colname) return apop_data_set_ti(d, rowname, col, val);
    if (!rowname && colname) return apop_data_set_it(d, row, colname, val);
    //else: row number, column number
    Set_gsl_handler
    if (col>=0){
        Apop_assert_negone(d->matrix, "You're trying to set the matrix element (%zu, %i) but the matrix is NULL.", row, col);
        gsl_matrix_set(d->matrix, row, col, val);
    } else {
        Apop_assert_negone(d->vector, "You're trying to set a vector element (row=-1) but the vector is NULL.");
        gsl_vector_set(d->vector, row, val);
    }
    Unset_gsl_handler
    return error_for_set;
APOP_VAR_ENDHEAD
    return 0; //the main function is blank.
}
/** \} //End data_set_get group */


/** Now that you've used \ref Apop_data_row to pull a row from an \ref apop_data set,
  this function lets you write that row to another position in the same data set or a
  different data set entirely.  

  The set written to must have the same form as the original: 
  \li a vector element has to be present if one existed in the original, 
  \li same for the weights vector,
  \li the matrix in the destination has to have as many columns as in the original, and
  \li the text has to have a row long enough to hold the original
  \li If the row to be written to already has a rowname, it is overwritten.
        If <tt>d->names->rowct == row_number</tt> (all rows up to \c row_number have row names), then extend the list of row names by one to add the new name. Else, don't add the row name. 
  \li Column names (of all types) aren't touched. Maybe use \c apop_data_copy or \c apop_name_copy if you need to copy these names.

  If any of the source elements are \c NULL, I won't bother to check that element in the
  destination.

  \return 0=OK, -1=error (probably a source/destination size mismatch).

  \li  The error codes for out-of-bounds errors are thread-safe iff you are have a
  C11-compliant compiler (thanks to the \c _Thread_local keyword) or GCC (thanks to its \c __thread extension).
*/
int apop_data_set_row(apop_data * d, apop_data *row, int row_number){
    Set_gsl_handler
    if (row->vector){
        Apop_assert_negone(d->vector, "You asked me to copy an apop_data_row with a vector element to "
                "an apop_data set with no vector.");
        gsl_vector_set(d->vector, row_number, row->vector->data[0]);
    }
    if (row->matrix && row->matrix->size2 > 0){
        Apop_assert_negone(d->matrix, "You asked me to copy an apop_data_row with a matrix row to "
                "an apop_data set with no matrix.");
        Apop_row(d, row_number, a_row); 
        Apop_row(row, 0, row_to_copy); 
        gsl_vector_memcpy(a_row, row_to_copy);
    }
    if (row->textsize[1]){
        Apop_assert_negone(d->textsize[1], "You asked me to copy an apop_data_row with text to "
                "an apop_data set with no text element.");
        for (int i=0; i < row->textsize[1]; i++){
            free(d->text[row_number][i]);
            d->text[row_number][i]= strdup(row->text[0][i]);
        }
    }
    if (row->weights){
        Apop_assert_negone(d->weights, "You asked me to copy an apop_data_row with a weight to "
                "an apop_data set with no weights vector.");
        gsl_vector_set(d->weights, row_number, row->weights->data[0]);
    }
    if (row->names && row->names->rowct && d->names){
        if (row_number < d->names->rowct){
            free(d->names->row[row_number]);
            d->names->row[row_number]=strdup(row->names->row[0]);
        } else if (row_number == d->names->rowct)
            apop_name_add(d->names, row->names->row[0], 'r');
    }
    Unset_gsl_handler
    return error_for_set;
}


/** A convenience function to add a named element to a data set.  Many of Apophenia's
testing procedures use this to easily produce a column of named parameters. It is public
as a convenience.

\param d    The \ref apop_data structure. Must not be \c NULL, but may be blank (as per
allocation via \ref apop_data_alloc <tt>( )</tt> ).
\param name The name to add
\param val  the value to add to the set.

\li I use the position of the name to know where to put the value. If
there are two names in the data set, then I will put the new name in
the third name slot and the data in the third slot in the zeroth column of the matrix. If
you use this function from start to finish in building your list, then you'll be fine.
\li If the matrix is too short (or \c NULL), I will call \ref apop_matrix_realloc internally to make space.
\li The list is in the zeroth column of the matrix, because that fits well with the
defaults for \ref apop_data_get. An example:

\code
apop_data *list = apop_data_alloc();
apop_data_add_named_elmt(list, "height", 165);
apop_data_add_named_elmt(list, "weight", 60);

double height = apop_data_get(list, .rowname="height");
\endcode
*/
void apop_data_add_named_elmt(apop_data *d, char *name, double val){
    Apop_assert_n(d, "You sent me a NULL apop_data set. Maybe allocate with apop_data_alloc() to start.");
    apop_name_add(d->names, name, 'r');
    if (!d->matrix)
        d->matrix = gsl_matrix_alloc(1, 1);
    if (d->matrix->size1 < d->names->rowct)
        apop_matrix_realloc(d->matrix, d->names->rowct, d->matrix->size2);
    gsl_matrix_set(d->matrix, d->names->rowct-1, 0, val);
}

//See apop_data_add_names in types.h.
void apop_data_add_names_base(apop_data *d, const char type, char const ** names){
    if (!d->names) d->names = apop_name_alloc();
    for(char const** name = names; *name !=NULL; name++)
        apop_name_add(d->names, *name, type);
}


/** Add a string to the text element of an \ref apop_data set.  If you
 send me a \c NULL string, I will write the string <tt>"NaN"</tt> in the given slot.
 If there is already something in that slot, that string is freed (preventing memory
 leaks).

\param in   The \ref apop_data set, that already has an allocated \c text element.
\param row  The row
\param col  The col
\param fmt The text to write.
\param ... You can use a printf-style fmt and follow it with the usual variables to fill in.

\return 0=OK, -1=error (probably out-of-bounds)

\li UTF-8 or ASCII text is correctly handled.
\li Apophenia follows a general rule of not reallocating behind your back: if your text
matrix is currently of size (3,3) and you try to put an item in slot (4,4), then I display
an error rather than reallocating the text matrix.
\li Resizing a text matrix is annoying in C, so note that \ref apop_text_alloc will
reallocate to a new size if you need. For example, this code will fill the diagonals of
the text array with a message, resizing as it goes:
\li The string added is a copy (via <tt>printf</tt>), not a pointer to the input(s).
\li If there had been a string at the grid point you are writing to, 
the old one is effectively lost when the new one is placed. So, I free
the old string to prevent leaks. Remember this if you had other pointers aliasing that
string, in which case you may as well avoid this function and just use <tt>
asprintf(&(your_dataset->text[row][col]), "your string")</tt>.

\code
apop_data *list = (something already allocated.);
for (int n=0; n < 10; n++){
    apop_text_alloc(list, n+1, n+1);
    apop_text_add(list, n, n, "This is cell (%i, %i)", n, n);
}
\endcode
*/
int apop_text_add(apop_data *in, const size_t row, const size_t col, const char *fmt, ...){
    Apop_assert_negone((in->textsize[0] >= (int)row+1) && (in->textsize[1] >= (int)col+1), "You asked me to put the text "
                            " '%s' at position (%zu, %zu), but the text array has size (%zu, %zu)\n", 
                               fmt,             row, col,                  in->textsize[0], in->textsize[1]);
    free (in->text[row][col]);
    if (!fmt){
        asprintf(&(in->text[row][col]), "%s", apop_opts.db_nan);
        return 0;
    }
    va_list argp;
	va_start(argp, fmt);
    vasprintf(&(in->text[row][col]), fmt, argp);
	va_end(argp);
    return 0;
}

/** This allocates an array of strings and puts it in the \c text element of an \ref apop_data set. 

  If the \c text element already exists, then this is effectively a \c realloc function,
  reshaping to the size you specify.

  \param in An \ref apop_data set. It's OK to send in \c NULL, in which case an apop_data set with \c NULL \c matrix and \c vector elements is returned.
  \param row    the number of rows of text.
  \param col     the number of columns of text.
  \return       A pointer to the relevant \ref apop_data set. If the input was not \c NULL, then this is a repeat of the input pointer.
  \exception out->error=='a'  Allocation error.
  */
apop_data * apop_text_alloc(apop_data *in, const size_t row, const size_t col){
    if (!in) in  = apop_data_alloc();
    if (!in->text){
        if (row){
            in->text = malloc(sizeof(char**) * row);
            Apop_stopif(!in->text, in->error='a'; return in, 
                    0, "malloc failed setting up %zu rows. Probably out of memory.", row);
        }
        if (row && col)
            for (size_t i=0; i< row; i++){
                in->text[i] = malloc(sizeof(char*) * col);
                Apop_stopif(!in->text[i], in->error='a'; return in, 
                        0, "malloc failed setting up row %zu (with %zu columns). Probably out of memory.", i, col);
                for (size_t j=0; j< col; j++)
                    in->text[i][j] = strdup("");
            }
    } else { //realloc
        size_t rows_now = in->textsize[0];
        size_t cols_now = in->textsize[1];
        if (rows_now > row){
            for (int i=row; i < rows_now; i++){
                for (int j=0; j < cols_now; j++)
                    free(in->text[i][j]);
                free(in->text[i]);
            }
            in->text = realloc(in->text, sizeof(char**)*row);
            Apop_stopif(row && !in->text, in->error='a'; return in,
                            0, "realloc failed shrinking down to %zu rows from %zu rows. "
                            "There may be actual bugs eating your computer.", row, rows_now);
        }
        if (rows_now < row){
            in->text = realloc(in->text, sizeof(char**)*row);
            Apop_stopif(!in->text, in->error='a'; return in,
                            0, "realloc failed setting up %zu rows. Probably out of memory.", row);
            for (size_t i=rows_now; i < row; i++){
                in->text[i] = malloc(sizeof(char*) * col);
                Apop_stopif(!in->text[i], in->error='a'; return in, 
                        0, "malloc failed setting up row %zu (with %zu columns). Probably out of memory.", i, col);
                for (int j=0; j < cols_now; j++)
                    in->text[i][j] = strdup("");
            }
        }
        if (cols_now > col)
            for (int i=0; i < row; i++)
                for (int j=col; j < cols_now; j++)
                    free(in->text[i][j]);
        if (cols_now != col)
            for (int i=0; i < row; i++){
                in->text[i] = realloc(in->text[i], sizeof(char*)*col);
                for (int j=cols_now; j < col; j++) //happens iff cols_now < col
                    in->text[i][j] = strdup("");
            }
    }
    in->textsize[0] = row;
    in->textsize[1] = col;
    return in;
}

/** Transpose the matrix element of the input \ref apop_data set,
 including the row/column names. Returns a copy of the element(s) that got transposed.
 The vector and weights elements of the input data set are completely ignored (but see also \ref apop_vector_to_matrix, which can convert a vector to a 1 X N matrix.)

This is mostly just a friendly wrapper for \c gsl_matrix_transpose_memcpy; if you have a \c gsl_matrix with no names or text, you may prefer to just use that function.

\param in The input \ref apop_data set. If \c NULL, I return \c NULL. Default is \c NULL.
\param transpose_text If \c 'y', then also transpose the text element. Default is \c 'y'.
\return  A newly alloced \ref apop_data set, with the appropriately transposed matrix and/or text. The vector and weights elements will be \c NULL. If <tt>transpose_text='n'</tt>, then the text element of the output set will also be \c NULL.

\li Row names are written to column names of the output matrix, text, or both (whichever is not empty in the input).
\li If only the matrix or only the text have names, then the one set of names is written to the row names of the output.
\li If both matrix column names and text column names are present, text column names are lost.

\li This function uses the \ref designated syntax for inputs.
 */ 
APOP_VAR_HEAD apop_data * apop_data_transpose(apop_data const *in, char transpose_text){
    apop_data const * apop_varad_var(in, NULL);
    Apop_stopif(!in, return NULL, 1, "Transposing a NULL data set; returning NULL.");
    char apop_varad_var(transpose_text, 'y');

APOP_VAR_ENDHEAD
    Apop_stopif(!in->matrix && !*in->textsize, return apop_data_alloc(), 
            1, "input data set has neither matrix nor text elements; returning an empty data set.");
    apop_data *out = apop_data_alloc(0, in->matrix ? in->matrix->size2 : 0
                                      , in->matrix ? in->matrix->size1 : 0);
    if (in->matrix) gsl_matrix_transpose_memcpy(out->matrix, in->matrix);
    apop_name_stack(out->names, in->names, 'r', 'c');
    apop_name_stack(out->names, in->names, 'c', 'r');
    if (!(transpose_text=='y')) return out;

    apop_text_alloc(out, in->textsize[1], in->textsize[0]);
    for (int r=0; r< in->textsize[0]; r++)
        for (int c=0; c< in->textsize[1]; c++)
            apop_text_add(out, c, r, in->text[r][c]);
    if (in->names->textct && !in->names->textct)
        apop_name_stack(out->names, in->names, 't', 'r');
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
                                    "I can't resize submatrices or other subviews.");
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
smaller than the old, then data at the end of the vector will be
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
    if (!v) return newheight ? gsl_vector_alloc(newheight) : NULL;
    apop_assert((v->block->data==v->data) && v->owner & (v->stride == 1),
                                    "I can't resize subvectors or other views.");
    v->block->size = newheight;
    v->size = newheight;
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
      
  \param match If \c 'c', case-insensitive match (via \c strcasecmp); if \c 'e', exact match, if \c 'r' regular expression substring search (via \ref apop_regex). Default=\c 'r'.

    \return The page whose title matches what you gave me. If I don't find a match, return \c NULL.

\li This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD apop_data * apop_data_get_page(const apop_data * data, const char *title, const char match){
    const apop_data * apop_varad_var(data, NULL);
    apop_assert_c(data, NULL, '1', "You requested a page from a NULL data set. Returning NULL");
    const char * apop_varad_var(title, "Info");
    const char apop_varad_var(match, 'r');
APOP_VAR_ENDHEAD
    while (data && (!data->names || 
                (match!='e' && match!='c' && !apop_regex(data->names->title, title))
                || (match=='c' && strcasecmp(data->names->title, title))
                || (match=='e' && strcmp(data->names->title, title))
                ))
        data = data->more;
    return (apop_data *) data; //de-const.
}

/** Add a page to a \ref apop_data set. It gets a name so you can find it later.

  \param dataset The input data set, to which a page will be added.
  \param newpage The page to append
  \param title The name of the new page. Remember, this is truncated at 100 characters.

  \return The new page.  I post a warning if I am appending or appending to a \c NULL data set and  <tt>apop_opts.verbose >=1 </tt>.

  \li Some data is fundamentally multi-page; an optimization search over multi-page
  parameters would search the space given by all pages, for example. 
  Also, pages may be appended as output or auxiliary information, such as covariances---an
  MLE would not search over these elements. Generally, any page with a name in XML-ish
  brackets, such as <tt>\<Covariance\></tt>, will be considered informational and ignored
  by search routines, missing data routines, et cetera. This is achieved by a rule in \ref
  apop_data_pack and \ref apop_data_unpack.

  Here is a toy example that establishes a baseline data set, adds a page,
  modifies it, and then later retrieves it.
  \code
  apop_data *d = apop_data_alloc(10, 10, 10); //the base data set.
  apop_data *a_new_page = apop_data_add_page(d, apop_data_alloc(2,2), "new 2 x 2 page");
  gsl_vector_set_all(a_new_page->matrix, 3);

  //later:
  apop_data *retrieved = apop_data_get_page(d, "new"); //uses regexes, not literal match.
  apop_data_show(retrieved); //print a 2x2 grid of 3s.
  \endcode
*/
apop_data * apop_data_add_page(apop_data * dataset, apop_data *newpage, const char *title){
    apop_assert_c(newpage, NULL, '1', "You are adding a NULL page to a data set. Doing nothing; returning NULL.");
    if (title && !(newpage->names->title == title))//has title, is not just pointint to existing title
        snprintf(newpage->names->title, 100, "%s", title);
    apop_assert_c(dataset, newpage, '1', "You are adding a page to a NULL data set. Returning the new page as its own data set.");
    while (dataset->more)
        dataset = dataset->more;
    dataset->more = newpage;
    return newpage;
}

/** Remove the first page from an \ref apop_data set that matches a given name.

  \param data The input data set, to which a page will be added. No default. If \c NULL, I return silently if <tt> apop_opts.verbose < 1 </tt>; print an error otherwise.
  \param title The name of the page to remove. Default: \c "Info"
  \param free_p If \c 'y', then \ref apop_data_free the page. Default: \c 'y'.

  \return If not freed, a pointer to the \c apop_data page that I just pulled out. Thus,
  you can use this to pull a single page from a data set. I set that page's \c more
  pointer to \c NULL, to minimize any confusion about more-than-linear linked list
  topologies. If <tt>free_p=='y'</tt> (the default) or the page is not found, return \c NULL.

  \li I don't check the first page, so there's no concern that the head of your list of
  pages will move. Again, the intent of the <tt>->more</tt> pointer in the \ref apop_data
  set is not to fully implement a linked list, but primarily to allow you to staple auxiliary
  information to a main data set.

  \li If I don't find the page you want, I return NULL, and print a message if
  <tt>apop_opts.verbose >= 1</tt>.
*/
APOP_VAR_HEAD apop_data* apop_data_rm_page(apop_data * data, const char *title, const char free_p){
    apop_data *apop_varad_var(data, NULL);
    apop_assert_c(data, NULL, '1', "You are removing a page from a NULL a data set. Doing nothing.");
    const char *apop_varad_var(title, "Info");
    const char apop_varad_var(free_p, 'y');
APOP_VAR_ENDHEAD
    while (data->more && !apop_regex(data->more->names->title, title))
        data = data->more;
    Apop_assert_c(data->more, NULL, 1, "You asked me to remove '%s' but I couldn't find a page matching that regex.", title);
    if (data->more){
        apop_data *tmp = data->more;
        data->more = data->more->more;
        tmp->more = NULL;
        if (free_p=='y'){
            free(tmp);
            return NULL;
        } //else:
        return tmp;
    } else return NULL;
}

typedef int (*apop_fn_ir)(apop_data*, void*);

/** Remove the columns set to one in the \c drop vector.  
  \param in the \ref apop_data structure to be pared down
  \param drop  a vector with as many elements as the max of the vector, matrix, or text
  parts of \c in, with a one marking those columns to be removed.       \ingroup names
  \param do_drop A function that returns one for rows to drop and zero for rows to not drop. A sample function:
  \code
  int your_drop_function(apop_data *onerow, void *extra_param){
    return gsl_isnan(apop_data_get(onerow)) || !strcmp(onerow->text[0][0], "Uninteresting data point");
  }
  \endcode
  \param drop_parameter If your \c do_drop function requires additional input, put it here and it will be passed through.
  \ref apop_data_rm_rows uses \ref Apop_data_row to get a subview of the input data set of height one (and since all the default arguments default to zero, you don't have to write out things like \ref apop_data_get <tt>(onerow, .row=0, .col=0)</tt>, which can help to keep things readable).

  \li If all the rows are to be removed, then you will wind up with the same \ref apop_data set, with \c NULL \c vector, \c matrix, \c weight, and text. Therefore, you may wish to check for \c NULL elements after use. I remove rownames, but leave the other names, in case you want to add new data rows.

 \li The typical use is to provide only a list or only a function. If both are \c NULL, I return without doing anything, and print a warning if <tt>apop_opts.verbose >=1</tt>. If you provide both, I will drop the row if either the vector has a one in that row's position, or if the function returns a nonzero value.
 \li This function uses the \ref designated syntax for inputs.
*/  
APOP_VAR_HEAD void apop_data_rm_rows(apop_data *in, int *drop, apop_fn_ir do_drop, void *drop_parameter ){
    apop_data* apop_varad_var(in, NULL);
    Apop_assert_c(in, , 1, "Input data set was NULL; no changes made.");
    int* apop_varad_var(drop, NULL);
    apop_fn_ir apop_varad_var(do_drop, NULL);
    void* apop_varad_var(drop_parameter, NULL);
    Apop_assert_c(drop || do_drop, , 1, "You gave me neither a list of ints "
            "indicating which rows to drop, nor a drop_fn I can use to test "
            "each row. Returning with no changes made.");
APOP_VAR_ENDHEAD
    //First, shift columns down to the nearest not-freed row.
    int outlength = 0;
    Get_vmsizes(in); //vsize, msize1, maxsize
    for (int i=0 ; i < maxsize; i++){
        int drop_row=0;
        if (drop && drop[i]) drop_row = 1;
        else if (do_drop){
            Apop_data_row(in, i, onerow); 
            drop_row = do_drop(onerow, drop_parameter);
        }
        if (!drop_row){
            if (outlength == i) outlength++;
            else {
                Apop_data_row(in, i, thisrow);
                apop_data_set_row(in, thisrow, outlength++);
            }
        }
    }
    if (!outlength){
        gsl_vector_free(in->vector);  in->vector = NULL;
        gsl_vector_free(in->weights); in->weights = NULL;
        gsl_matrix_free(in->matrix);  in->matrix = NULL;
        apop_text_alloc(in, 0, 0);
        //leave colnames intact, remove rownames below.
    }

    //now trim excess memory:
    if (in->vector)  apop_vector_realloc(in->vector, GSL_MIN(in->vector->size, outlength));
    if (in->weights) apop_vector_realloc(in->weights, GSL_MIN(in->weights->size, outlength));
    if (in->matrix)  apop_matrix_realloc(in->matrix, GSL_MIN(in->matrix->size1, outlength), in->matrix->size2);
    if (in->text)    apop_text_alloc(in, GSL_MIN(outlength, in->textsize[0]), in->textsize[1]);
    if (in->names && in->names->rowct > outlength){
        for (int k=outlength; k< in->names->rowct; k++)
            free(in->names->row[k]);
        in->names->rowct = outlength;
    }
}
