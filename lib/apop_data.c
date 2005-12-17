/** \file apop_data.c	 The apop_data structure joins together a
  gsl_matrix, apop_name, and a table of strings. No biggie.


Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL.
*/

#include <gsl/gsl_matrix.h>
#include "apophenia/types.h"
#include "apophenia/output.h"
#include "apophenia/linear_algebra.h"

/** \defgroup data_struct

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

  \param m   The existing matrix you'd like to turn into an \ref apop_data structure.
 \return    The \ref apop_data structure in question.
  */
apop_data * apop_matrix_to_data(gsl_matrix *m){
apop_data  *setme   = malloc(sizeof(apop_data));
    if (m==NULL && apop_verbose) {printf("Warning: converting a NULL matrix to an apop_data structure.\n");}
    setme->data     = m;
    setme->names    = apop_name_alloc();
    setme->categories = NULL;
    return setme;
}

/** Free a \ref apop_name structure.

 \ingroup data_struct
  */
void apop_data_free(apop_data *freeme){
    gsl_matrix_free(freeme->data);
    apop_name_free(freeme->names);
    if (freeme->categories !=NULL)
        free(freeme->categories);
    free(freeme);
}
