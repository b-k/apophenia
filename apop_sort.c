/** \file apop_sort.c 

  A few functions to sort data. One sorts an \c apop_data set in place, and one returns percentiles for a sorted vector.
  The headers are in stats.h

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "likelihoods.h"
#include <gsl/gsl_sort_vector.h>

static int find_min_unsorted(size_t *sorted, size_t height, size_t min){
    while (min<height)
        if (!sorted[min])   return min;
        else                min++;
    return -1;
}

static void temp_in_out(apop_data *d, gsl_vector *tv, double *tvd, size_t i, char **namein, char inout){
    if (inout == 'i'){
        if (d->matrix){
            APOP_ROW(d, i, row);
            gsl_vector_memcpy(tv, row);
        }
        if (d->vector)
            *tvd    = apop_data_get(d, i, -1);
        if (d->names && d->names->rowct > i){
            *namein = realloc(*namein, strlen(d->names->row[i])+1);
            strcpy(*namein, d->names->row[i]);
        }
        return;
    }//else writing out
    if (d->matrix){
        APOP_ROW(d, i, row);
        gsl_vector_memcpy(row, tv);
    }
    if (d->vector){
        double *j   = apop_data_ptr(d, i, -1);
        *j          = *tvd;
    }
    if (d->names && d->names->rowct > i){
        d->names->row[i] = realloc(d->names->row[i], strlen(*namein)+1);
        strcpy(d->names->row[i], *namein);
    }
}

static void shift(apop_data *d, size_t from,  size_t to){
  double *dp    = d->vector ? apop_data_ptr(d, to, -1): 0;
  char *namep  = d->names && d->names->rowct > to ? d->names->row[to] : NULL;
    if (d->matrix){
        APOP_ROW(d, to, tov);
        temp_in_out(d, tov, dp, from, &namep, 'i');
        return;
    }
    temp_in_out(d, NULL, dp , from, &namep, 'i');
}

/** This function sorts the whole of a \c apop_data set based on one
 column. Sorts in place, with very little additional memory used.

 Uses the \c gsl_sort_vector_index function internally, and that function just ignores NaNs; therefore this function just leaves NaNs exactly where they lay.

 \param data    The input set to be modified.
 \param sortby  The column of data by which the sorting will take place. As usual, -1 indicates the vector element.
 \param asc   If 'd' or 'D', sort in descending order; else sort in ascending order.
 \return A pointer to the data set, so you can do things like \c apop_data_show(apop_data_sort(d, -1, 'a')).
 \ingroup data_struct
*/
apop_data * apop_data_sort(apop_data *data, int sortby, char asc){
  size_t            i, j;
  size_t            height  = (sortby==-1) ? data->vector->size: data->matrix->size1;
  size_t            sorted[height];
  size_t            *perm, start=0;
  gsl_permutation   *p  = gsl_permutation_alloc(height);
  gsl_vector        *tv = data->matrix ? gsl_vector_alloc(data->matrix->size2) : NULL;
  double            tvd;
  char              *nametmp = malloc(1);
    memset(sorted, 0, sizeof(size_t)*height);
    if (sortby == -1)
        gsl_sort_vector_index (p, data->vector);
    else {
        APOP_COL(data, sortby, v);
        gsl_sort_vector_index (p, v);
    }
    perm    = p->data;
    if (asc=='d' || asc=='D')
        for (j=0; j< height/2; j++){
            tvd            = perm[j];
            perm[j]        = perm[height-1-j];
            perm[height-1-j] = tvd;
        }
    while (1){
        i           =
        start       = find_min_unsorted(sorted, height, start);
        if (i==-1) goto finished;
        temp_in_out(data, tv, &tvd, start, &nametmp, 'i');
        sorted[start]++;
        while (perm[i]!=start){
            shift(data, perm[i], i);
            sorted[perm[i]]++;
            i   = perm[i];
        }
        temp_in_out(data, tv, &tvd, i, &nametmp, 'o');
    }
finished:
    if (tv) gsl_vector_free(tv);
    if (nametmp) free(nametmp);
    gsl_permutation_free(p);
    return data;
}


/** Returns a vector of size 101, where returned_vector[95] gives the
value of the 95th percentile, for example. Returned_vector[100] is always
the maximum value, and returned_vector[0] is always the min (regardless
of rounding rule).

\param data	a gsl_vector of data.
\param rounding This will either be 'u', 'd', or 'a'. Unless your data is
exactly a multiple of 101, some percentiles will be ambiguous. If 'u',
then round up (use the next highest value); if 'd' (or anything else),
round down to the next lowest value; if 'a', take the mean of the two nearest points. If 'u' or 'a', then you can say "5% or
more  of the sample is below returned_vector[5]"; if 'd' or 'a', then you can
say "5% or more of the sample is above returned_vector[5]".  

\ingroup basic_stats
*/ 
double * apop_vector_percentiles(gsl_vector *data, char rounding){
  gsl_vector	*sorted	= gsl_vector_alloc(data->size);
  double		*pctiles= malloc(sizeof(double) * 101);
  int		i, index;
	gsl_vector_memcpy(sorted,data);
	gsl_sort_vector(sorted);
	for(i=0; i<101; i++){
		index = i*(data->size-1)/100.0;
		if (rounding == 'u' && index != i*(data->size-1)/100.0)
			index ++; //index was rounded down, but should be rounded up.
		if (rounding == 'a' && index != i*(data->size-1)/100.0)
            pctiles[i]	= (gsl_vector_get(sorted, index)+gsl_vector_get(sorted, index+1))/2.;
        else pctiles[i]	= gsl_vector_get(sorted, index);
	}
	gsl_vector_free(sorted);
	return pctiles;
}
