/** \file apop_data.c	 The apop_data structure joins together a
  gsl_matrix, apop_name, and a table of strings. No biggie.


Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL.
*/

#include <gsl/gsl_matrix.h>
#include "apophenia/types.h"
#include "apophenia/output.h"
#include "apophenia/linear_algebra.h"

/** \defgroup data_struct apop_data

  The \c apop_data structure represents a data set.  It joins together a
  gsl_matrix, apop_name, and a table of strings. No biggie. It tries to be minimally intrusive, so you can use it everywhere you would use a \c gsl_matrix.

  */

/** Allocate a \ref apop_name structure, to be filled with data.

  \param size1, size2   row and column size for the matrix. Notice that this exactly mirrors the format of \c gsl_matrix_alloc.
 \return    The \ref apop_data structure in question.
 \ingroup data_struct
  */
apop_data * apop_data_alloc(int size1, int size2){
apop_data  *setme   = malloc(sizeof(apop_data));
    setme->data     = gsl_matrix_alloc(size1,size2);
    setme->names    = apop_name_alloc();
    setme->categories = NULL;
    return setme;
}

/** Wrap an \ref apop_name structure around an existing \c gsl_matrix.

\param m    The existing matrix you'd like to turn into an \ref apop_data structure.
return      The \ref apop_data structure in question.
  */
apop_data * apop_matrix_to_data(gsl_matrix *m){
apop_data  *setme   = malloc(sizeof(apop_data));
    if (m==NULL && apop_verbose) {printf("Warning: converting a NULL matrix to an apop_data structure.\n");}
    setme->data     = m;
    setme->names    = apop_name_alloc();
    setme->categories = NULL;
    return setme;
}

/** Free an \ref apop_name structure.

 \ingroup data_struct
  */
void apop_data_free(apop_data *freeme){
    gsl_matrix_free(freeme->data);
    apop_name_free(freeme->names);
    if (freeme->categories !=NULL)
        free(freeme->categories);
    free(freeme);
}


/** Put the first data set either on top of or to the right of the second matrix.

The fn returns a new data set, meaning that at the end of this function,
until you apop_data_free() the original data sets, you will be taking up
twice as much memory. Plan accordingly.

The dependent variable names are not copied or modified in any way by
this function, so if you mean for them to change or get stacked, you
will have to make the modifications yourself, perhaps by using either \ref apop_name_add or \ref apop_name_stack.

\param  m1      the upper/rightmost data set
\param  m2      the second data set
\param  posn    if 'r', stack rows of m1 above rows of m2, else, e.g. 'c', stack m1's columns to right of m2's.
\return         a new \ref apop_data set with the stacked data.
\ingroup data_struct
*/
apop_data *apop_data_stack(apop_data *m1, apop_data * m2, char posn){
gsl_matrix  *stacked    = apop_matrix_stack(m1->data, m2->data, posn);
apop_data   *out        = apop_matrix_to_data(stacked);
    memcpy(out->names, m1->names, sizeof(apop_name));
    apop_name_stack(out->names, m2->names, posn);
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
 */
void apop_data_rm_columns(apop_data *d, int *drop){
gsl_matrix  *freeme = d->data;
    d->data = apop_matrix_rm_columns(d->data, drop);
    gsl_matrix_free(freeme); 
    apop_name_rm_columns(d->names, drop);
}
