
/** \file apop_sort.c
Copyright (c) 2013 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
#include <stdbool.h>

static double find_smallest_larger_than(apop_data const *sort_order, double *x){
    //the next column in the sort order is the one that is not NAN, greater than x, but smaller than all other candidate values.
    double candidate_col=-100, candidate_val = INFINITY, v;
    if (sort_order->vector && !isnan(v=gsl_vector_get(sort_order->vector, 0)) 
            && v > *x && v < candidate_val){
        candidate_val = v;
        candidate_col = -1;
    }
    if (sort_order->matrix)
        for (int i=0; i< sort_order->matrix->size2; i++)
            if (!isnan(v=gsl_matrix_get(sort_order->matrix, 0, i)) 
                    && v > *x && v < candidate_val){
                candidate_val = v;
                candidate_col = i;
            }
    if (*sort_order->textsize)
        for (int i=0; i< sort_order->textsize[1]; i++)
            if (apop_opts.nan_string && strcmp(sort_order->text[0][i], apop_opts.nan_string)
                    && (v=atof(sort_order->text[0][i])) > *x && v < candidate_val){
                candidate_val = v;
                candidate_col = i+0.5;
            }
    if (sort_order->weights && !isnan(v=gsl_vector_get(sort_order->weights, 0)) 
            && v > *x && v < candidate_val){
        candidate_val = v;
        candidate_col = -2;
    }
    if (sort_order->names->rowct)
        if (apop_opts.nan_string && strcmp(*sort_order->names->row, apop_opts.nan_string)
                && (v=atof(*sort_order->names->row)) > *x && v < candidate_val){
            candidate_val = v;
            candidate_col = 0.2;
        }

    *x = candidate_val;
    return candidate_col;
}

static void generate_sort_order(apop_data const *data, apop_data const *sort_order, int cols_to_sort_ct, double *so){
/* the internal rule is that the vector is -1, the weights vector is -2, the names are
 * 0.2, and the text cols are the column+0.5. How's that for arbitrary. */
    if (sort_order) {
        double x = -INFINITY;
        for (int i=0; i< cols_to_sort_ct-1; i++)
            so[i] = find_smallest_larger_than(sort_order, &x);
    } else {
        int ctr=0;
        if (data->vector) so[ctr++] = -1;
        if (data->matrix) for(int i=0; i<data->matrix->size2; i++) so[ctr++] = i;
        if (*data->textsize) for(int i=0; i< data->textsize[1]; i++) so[ctr++] = i+0.5;
        if (data->weights) so[ctr++] = -2;
        if (data->names->rowct) so[ctr++] = 0.2;
    }
    so[cols_to_sort_ct-1] = -100;
}

#include <gsl/gsl_sort_vector.h>

static int find_min_unsorted(size_t *sorted, size_t height, size_t min){
    while (min<height)
        if (!sorted[min]) return min;
        else              min++;
    return -1;
}

static threadlocal apop_data *d;  //stdlib qsort doesn't have a hook where we can put these.
static threadlocal int offset;

static int compare_strings(const void *a, const void *b) {
    const size_t *da = (const size_t *) a;
    const size_t *db = (const size_t *) b;
    return offset==-1
        ? strcasecmp(d->names->row[*da], d->names->row[*db])
        : strcasecmp(d->text[*da][offset], d->text[*db][offset]);
}

static void rearrange(apop_data *data, size_t height, size_t *perm){
    size_t i, start=0;
    size_t sorted[height];
    memset(sorted, 0, sizeof(size_t)*height);
    while (1){
        i     =
        start = find_min_unsorted(sorted, height, start);
        if (i==-1) break;
        Apop_row(data, start, firstrow);
        apop_data *first_row_storage = apop_data_copy(firstrow);
        sorted[start]++;
        while (perm[i]!=start){
            //copy from perm[i] to i
            Apop_row(data, perm[i], onerow);
            apop_data_set_row(data, onerow, i);
            sorted[perm[i]]++;
            i = perm[i];
        }
        apop_data_set_row(data, first_row_storage, i);
        apop_data_free(first_row_storage);
    }
}

/** Sort an \ref apop_data set on an arbitrary sequence of columns. 

The \c sort_order set is a one-row data set that should look like the data set being
sorted. The easiest way to generate it is to use \ref Apop_row to pull one row of the
table, then copy and fill it. For each column you want used in the sort, assign a ranking giving whether the column should be sorted first, second, .... Columns you don't want used in the sorting should be set to \c NAN. Ties are broken by the earlier element in the default order (see below).

E.g., to sort by the last column of a five-column matrix first, then the next-to-last column, then the next-to-next-to-last, then by the first text column, then by the second text column:

\code
Apop_row(data, 0, so)
apop_data *sort_order = apop_data_copy(so);
sort_order->vector = NULL; //so it will be skipped.
Apop_data_fill(sort_order, NAN, NAN, 3, 2, 1);
apop_text_add(sort_order, 0, 0, "4");
apop_text_add(sort_order, 0, 1, "5");
apop_data_sort(data, sort_order);
\endcode

I use only comparisons, not the actual numeric values, so you can use any sequence of
numbers: (1, 2, 3) and (-1.32, 0, 27) work identically.

\li Strings are sorted case-insensitively, using \c strcasecmp. [exercise for the reader: modify the source to use Glib's locale-correct string sorting.]

\li The setup generates a lexicographic sort using the columns you specify. If you would like a different sort order, such as Euclidian distance to the origin, you can generate a new column expressing your preferred metric, and then sorting on that. See the example below.

\param data The data set to be sorted. If \c NULL, this function is a no-op that returns \c NULL.
\param sort_order A \ref apop_data set describing the order in which columns are used for sorting, as above. If \c NULL, then sort by the vector, then each matrix column, then text, then weights, then row names.
\param inplace If 'n', make a copy, else sort in place. (default: 'y').
\param asc If 'a', ascending; if 'd', descending. This is applied to all columns; column-by-column application is to do. (default: 'a').
\param col_order For internal use only. In your call, it should be \c NULL; the \ref designated syntax will takes care of it for you.

\return A pointer to the sorted data set. If <tt>inplace=='y'</tt> (the default), then this is the same as the input set.


A few examples:

\include "../tests/sort_example.c"

\li This function uses the \ref designated syntax for inputs.
*/
#ifdef APOP_NO_VARIADIC
apop_data * apop_data_sort(apop_data *data, apop_data *sort_order, char asc, char inplace, double *col_order){
#else
apop_varad_head(apop_data *, apop_data_sort){
    apop_data * apop_varad_var(data, NULL);
    Apop_stopif(!data, return NULL, 1, "You gave me NULL data to sort. Returning NULL");
    apop_data * apop_varad_var(sort_order, NULL);
    char apop_varad_var(inplace, 'y');
    char apop_varad_var(asc, 'a');
    double * apop_varad_var(col_order, NULL);
    return apop_data_sort_base(data, sort_order, asc, inplace, col_order);
}

 apop_data * apop_data_sort_base(apop_data *data, apop_data *sort_order, char asc, char inplace, double *col_order){
#endif
    if (!data) return NULL;

    apop_data *out = inplace=='n' ? apop_data_copy(data) : data;

    apop_data *xx = sort_order ? sort_order : out;
    Get_vmsizes(xx); //firstcol, msize2
    int cols_to_sort_ct = msize2 - firstcol +1 + !!(xx->weights) + xx->textsize[1] + !!xx->names->rowct;
    double so[cols_to_sort_ct];
    if (!col_order){
        generate_sort_order(out, sort_order, cols_to_sort_ct, so);
        col_order = so;
    }

    bool is_text = ((int)*col_order != *col_order);
    bool is_name = (*col_order == 0.2);

    gsl_vector_view c;
    gsl_vector *cc = NULL;
    if (!is_text && *col_order>=0){
        c = gsl_matrix_column(out->matrix, *col_order);
        cc = &c.vector;
    }
    gsl_vector *thiscol =   cc               ? cc
                          : (*col_order==-2) ? out->weights
                          : (*col_order==-1) ? out->vector
                                             : NULL;

    size_t height = thiscol   ? thiscol->size
                    : is_name ? out->names->rowct
                              : *out->textsize;

    gsl_permutation *p = gsl_permutation_alloc(height);
    if (!is_text) gsl_sort_vector_index (p, thiscol);
    else {
        gsl_permutation_init(p);
        d = out;
        offset = is_name ? -1 : *col_order-0.5;        
        qsort(p->data, height, sizeof(size_t), compare_strings);
    }

    size_t *perm = p->data;
    if (asc=='d' || asc=='D') //reverse the perm matrix.
        for (size_t j=0; j< height/2; j++){
            double t         = perm[j];
            perm[j]          = perm[height-1-j];
            perm[height-1-j] = t;
        }
    rearrange(out, height, perm);
    gsl_permutation_free(p);
    if (col_order[1] == -100) return out;

    /*Second pass:
    find blocks where all are of the same value.
    After you pass a block of size > 1 row where all vals in this col are identical,
    sort that block, using the rest of the sort order. */
    int bottom=0;
    if (!is_text){
        double last_val = gsl_vector_get(thiscol, 0);
        for (int i=1; i< height+1; i++){
            double this_val=0;
            if ((i==height || (this_val=gsl_vector_get(thiscol, i)) != last_val) 
                    && bottom != i-1){
                Apop_rows(out, bottom, i-bottom, subset);
                apop_data_sort_base(subset, sort_order, 'a', 'y', col_order+1);
            }
            if (last_val != this_val) bottom = i;
            last_val = this_val;
        }
    } else {
        char *last_val =  is_name ? out->names->row[0] : out->text[0][(int)(*col_order-0.5)];
        for (int i=1; i< height+1; i++){
            char *this_val = i==height ? NULL : is_name ? out->names->row[i] : out->text[i][(int)(*col_order-0.5)];
            if ((i==height || strcasecmp(this_val, last_val)) 
                    && bottom != i-1){
                Apop_rows(out, bottom, i-bottom, subset);
                apop_data_sort_base(subset, sort_order, 'a', 'y', col_order+1);
            }
            if (this_val && strcmp(last_val, this_val)) bottom = i;
            last_val = this_val;
        }
    }
    return out;
}
