/* Probability mass functions 
Copyright (c) 2011 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  

\amodel apop_pmf A probability mass function is commonly known as a histogram, or still more commonly,
a bar chart. It indicates that at a given coordinate, there is a given mass.

The data format for the PMF is simple: each row holds the coordinates, and the
<em>weights vector</em> holds the mass at the given point. This is in contrast to the
crosstab format, where the location is simply given by the position of the data point
in the grid.

For example, here is a typical crosstab:

<table>
<tr>            <td></td><td> col 0</td><td> col 1</td><td> col 2</td></tr>
<tr><td>row 0 </td><td>0</td><td> 8.1</td><td> 3.2</td></tr>
<tr><td>row 1 </td><td>0</td><td> 0</td><td> 2.2</td></tr>
<tr><td>row 2 </td><td>0</td><td> 7.3</td><td> 1.2</td></tr>
</table>

Here it is as a sparse listing:

<table>
<tr>        <td></td> value</td></td> dimension 1<td></td> dimension 2<td></tr>
<tr> <td>8.1</td> <td>0</td> <td>1</td> </tr>
<tr> <td>3.2</td> <td>0</td> <td>2</td> </tr>
<tr> <td>2.2</td> <td>1</td> <td>2</td> </tr>
<tr> <td>7.3</td> <td>2</td> <td>1</td> </tr>
<tr> <td>1.2</td> <td>2</td> <td>2</td> </tr>
</table>

The \c apop_pmf internally represents data in this manner. The dimensions are held
in the \c matrix element of the data set, and the cell values are held in the \c weights
element (<em>not the vector</em>).

If your data is in a crosstab (with entries in the matrix element for 2-D data or the
vector for 1-D data), then use \ref apop_crosstab_to_pmf to make the conversion.

If your data is already in the sparse listing format (which is probably the case for 3-
or more dimensional data), then just point the model to your parameter set:

\code
apop_model *my_pmf = apop_model_copy(apop_pmf);
my_pmf->more = in_data;
//or equivalently:
apop_model *my_pmf = apop_estimate(in_data, apop_pmf);
\endcode

\li If the \c weights element is \c NULL, then I assume that all rows of the data set are
equally probable.

\li Be careful: the weights are in the \c weights element of the \c apop_data set, not in
the \c vector element. If you put the weights in the \c vector and have \c NULL \c
weights, then draws are equiprobable. This will be difficult to debug.

\adoc    Input_format     As above, you can input to the \c estimate
                      routine a 2-D matrix that will be converted into this form.     
\adoc    Parameter_format  None. The list of observations and their weights are in the \c data set, not the \c parameters.
\adoc    Settings   None.    
*/

#include "apop_internal.h"

/* \adoc    estimated_data  The data you sent in is linked to (not copied).
\adoc    estimated_parameters  Still \c NULL.    */
static apop_model *estim (apop_data *d, apop_model *out){
    out->data = d;
    apop_data_free(out->parameters); //may have been auto-alloced by prep.
    out->parameters = apop_data_copy(d);
    return out;
}

/* \adoc    RNG  The first time you draw from a PMF, I will generate a CMF (Cumulative
Mass Function). For an especially large data set this may take a human-noticeable amount
of time. The CMF will be stored in <tt>parameters->weights[1]</tt>, and subsequent
draws will have no computational overhead. */
static void draw (double *out, gsl_rng *r, apop_model *m){
    Nullcheck_m(m) Nullcheck_d(m->data)
    Get_vmsizes(m->data)
    size_t current; 
    if (!m->data->weights) //all rows are equiprobable
        current = gsl_rng_uniform(r)* (GSL_MAX(vsize,msize1)-1);
    else {
        size_t size = m->data->weights->size;
        if (!m->more){
            m->more = gsl_vector_alloc(size);
            gsl_vector *cdf = m->more; //alias.
            cdf->data[0] = m->data->weights->data[0];
            for (int i=1; i< size; i++)
                cdf->data[i] = m->data->weights->data[i] + cdf->data[i-1];
            //Now make sure the last entry is one.
            gsl_vector_scale(cdf, 1./cdf->data[size-1]);
        }
        double draw = gsl_rng_uniform(r);
        //do a binary search for where draw is in the CDF.
        double *cdf = ((gsl_vector*)m->more)->data; //alias.
        size_t top = size-1, bottom = 0; 
        current = (top+bottom)/2.;
        if (current==0){//array of size one or two
            if (size!=1) 
                if (cdf[0] < draw)
                    current = 1;
        } else while (!(cdf[current]>=draw && cdf[current-1] < draw)){
            if (cdf[current] < draw){ //step up
                bottom = current;
                if (current == top-1)
                    current ++;
                else
                    current = (bottom + top)/2.;
            } else if (cdf[current-1] >= draw){ //step down
                top = current;
                if (current == bottom+1)
                    current --;
                else
                    current = (bottom + top)/2.;
            }
        }
    }
    //Done searching. Current should now be the right row index.
    Apop_data_row(m->data, current, outrow);
    int i = 0;
    if (outrow->vector)
        out[i++] = outrow->vector->data[0];
    if (outrow->matrix)
        for( ; i < outrow->matrix->size2; i ++)
            out[i] = gsl_matrix_get(outrow->matrix, 0, i);
}


static int are_equal(apop_data *left, apop_data *right){
    /* Intended by use for apop_data_pmf_compress and .p, below.
      That means we aren't bothering with comparing names, and weights are likely to be
      different, because we're using those to tally data elements. If the data set has
      a longer matrix than vector, say, then one side may have the vector element and
      the other not, so we still check that there's a match in presence of each element.

     */
    if (left->vector){
        if (left->vector->data[0] != right->vector->data[0] 
             && !(gsl_isnan(left->vector->data[0]) && gsl_isnan(right->vector->data[0]))) 
            return 0;
    } else if (right->vector) return 0;

    if (left->matrix){
        if (left->matrix->size2 != right->matrix->size2) return 0;
        for (int i=0; i< left->matrix->size2; i++){
            double L = apop_data_get(left, 0, i);
            double R = apop_data_get(right, 0, i);
            if (L != R && !(gsl_isnan(L) && gsl_isnan(R))) return 0;
        }
    }
    else if (right->matrix) return 0;

    if (left->textsize[1]){
        if (left->textsize[1] != right->textsize[1]) return 0;
        for (int i=0; i< left->textsize[1]; i++)
            if (!apop_strcmp(left->text[0][i], right->text[0][i])) return 0;
    }
    else if (right->textsize[1]) return 0;
    return 1;
}

static int find_in_data(apop_data *searchme, apop_data *findme){//findme is one row tall.
    Get_vmsizes(searchme)
    for(int i=0; i < GSL_MAX(vsize, GSL_MAX(searchme->textsize[0], msize1)); i++){
        Apop_data_row(searchme, i, onerow);
        if (are_equal(findme, onerow))
            return i;
    }
    return -1;
}
double pmf_p(apop_data *d, apop_model *m){
    Nullcheck_d(d) 
    Nullcheck_m(m) 
    int model_pmf_length;
    {
        Get_vmsizes(m->data);
        model_pmf_length = GSL_MAX(vsize, GSL_MAX(m->data->textsize[0], msize1));
    }
    Get_vmsizes(d)//firstcol, vsize, vsize1, msize2
    long double p = 1;
    for(int i=0; i< GSL_MAX(msize1, d->textsize[0]); i++){
        Apop_data_row(d, i, onerow);
        int elmt = find_in_data(m->data, onerow);
        if (elmt == -1) return 0; //Can't find one observation: prob=0;
        p *= m->data->weights
                 ? m->data->weights->data[elmt]
                 : 1./model_pmf_length; //no weights means any known event is equiprobable
    }
    return p;
}

static void pmf_print(apop_model *est){ apop_data_print(est->data); }

apop_model apop_pmf = {"PDF or sparse matrix", .dsize=-1, .estimate = estim, .draw = draw, .p=pmf_p, .print=pmf_print};


/** Say that you have added a long list of observations to a single \ref apop_data set,
  meaning that each row has weight one. There are a huge number of duplicates, perhaps because there are a handful of 
  types that keep repeating:

<table frame=box>
<tr>
<td>Vector value</td><td> Text name</td><td>Weights</td>
</tr><tr valign=bottom>
<td align=center>
</td></tr>
<tr><td>12</td><td>Dozen</td><td>1</td></tr>
<tr><td>1</td><td>Single</td><td>1</td></tr>
<tr><td>2</td><td>Pair</td><td>1</td></tr>
<tr><td>2</td><td>Pair</td><td>1</td></tr>
<tr><td>1</td><td>Single</td><td>1</td></tr>
<tr><td>1</td><td>Single</td><td>1</td></tr>
<tr><td>2</td><td>Pair</td><td>1</td></tr>
<tr><td>2</td><td>Pair</td><td>1</td></tr>
</table>

You would like to reduce this to a set of distinct values, with their weights adjusted accordingly:

<table frame=box>
<tr>
<td>Vector value</td><td> Text name</td><td>Weights</td>
</tr><tr valign=bottom>
<td align=center>
</td></tr>
<tr><td>12</td><td>Dozen</td><td>1</td></tr>
<tr><td>1</td><td>Single</td><td>3</td></tr>
<tr><td>2</td><td>Pair</td><td>4</td></tr>
</table>


\param in An \ref apop_data set that may have duplicate rows. As above, the data may be in text and/or numeric formats. If there is a \c weights vector, I will add those weights together as duplicates are merged. If there is no \c weights vector, I will create one, which is initially set to one for all values, and then aggregated as above.

\return Your input is changed in place, via \ref apop_data_rm_rows, so use \ref apop_data_copy before copying this function if you need to retain the original format. For your convenience, this function returns a pointer to your original data, which has now been pruned.

*/
apop_data *apop_data_pmf_compress(apop_data *in){
    Apop_assert_c(in, NULL, 1,  "You sent me a NULL input data set; returning NULL output.");
    Get_vmsizes(in);
    size_t max = GSL_MAX(vsize, GSL_MAX(msize1, in->textsize[0]));
    if (!in->weights){
        in->weights=gsl_vector_alloc(max);
        gsl_vector_set_all(in->weights, 1);
    }
    int *cutme = calloc(max, sizeof(int));
    for (int i=0; i< max;i++){
        Apop_data_row(in, i, subject);
        for (int j=0; j< i; j++){
            Apop_data_row(in, j, compare_me);
            if (are_equal(subject, compare_me)){
                apop_vector_increment(compare_me->weights, 0,
                                    gsl_vector_get(subject->weights, 0));
                cutme[i]=1;
                break; //j-loop only
            }
        }
    }
    apop_data_rm_rows(in, cutme);
    free(cutme);
    return in;
}

