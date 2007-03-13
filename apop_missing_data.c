
#include <apop.h>

/** If there is an NaN anywhere in the row of data (including the matrix
    and the vector) then delete the row from the data set.

    The function returns a new data set with the NaNs removed, so
    the original data set is left unmolested. You may want to \c
    apop_data_free the original immediately after this function.

    If every row has an NaN, then this returns NULL; you may want to
    check for this after the function returns.

    \param d    The data, with NaNs
    \return     A (potentially shorter) copy of the data set, without NaNs.
*/
apop_data * apop_data_listwise_delete(apop_data *d){
  int i, j, min = 0, max = 0, height=0, has_vector=0, has_matrix=0, to_rm;
    //get to know the input.
    if (d->matrix){
        height      = d->matrix->size1;
        max         = d->matrix->size2;
        has_matrix  ++;
    } 
    if (d->vector){
        height      = height ?  height : d->vector->size;
        min         = -1;
        has_vector  ++;
    } 
    if (!has_matrix && !has_vector) {
        fprintf(stderr, "You sent to apop_data_listwise_delete a data set with void matrix and vector. Confused, it is returning NULL.\n");
        return NULL;
        }
    //find out where the NaNs are
  gsl_vector *marked = gsl_vector_calloc(height);
    for (i=0; i< d->matrix->size1; i++)
        for (j=min; j <max; j++)
            if (gsl_isnan(apop_data_get(d, i, j))){
                    gsl_vector_set(marked, i, 1);
                    break;
            }
    to_rm   = apop_sum(marked);
    //copy the good data.
    if (to_rm  == height)
        return NULL;
  apop_data *out = apop_data_alloc(height-to_rm, has_matrix ? max : -1);
    out->names  = apop_name_copy(d->names);                           ///You loser!!! Fix this. And add text!!!!
    if (has_vector && has_matrix)
        out->vector = gsl_vector_alloc(height - to_rm);
    j   = 0;
    for (i=0; i< height; i++){
        if (!gsl_vector_get(marked, i)){
            if (has_vector)
                gsl_vector_set(out->vector, j, gsl_vector_get(d->vector, i));
            if (has_matrix){
                APOP_ROW(d, i, v);
                gsl_matrix_set_row(out->matrix, j, v);
            if (d->names->rownames && d->names->rownamect > i)
                apop_name_add(out->names, d->names->rownames[i], 'r');
            }
            j++;
        }
    }
    gsl_vector_free(marked);
    return out;
}
