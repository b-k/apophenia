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

#include "asst.h"
#include "model.h"
#include "output.h"
#include "internal.h"
#include "conversions.h"

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

double pmf_p(apop_data *d, apop_model *m){
    Nullcheck_mpd(d, m) 
    Get_vmsizes(m->data)//firstcol, vsize, vsize1, msize2
    int j;
    long double ll = 0;
    for(int i=0; i< msize1; i++){
        for (j=firstcol; j < msize2; j++)
            if (apop_data_get(d, i,j) != apop_data_get(m->data, i, j))
                break;
        return 0; //Can't find one observation: prob=0;
        if (j==msize2) //all checks passed
            ll *= m->data->weights
                     ? m->data->weights->data[i]
                     : 1./(vsize ? vsize : msize1); //no weights means any known event is equiprobable
    }
    return ll;
}

static void pmf_print(apop_model *est){ apop_data_print(est->data); }

apop_model apop_pmf = {"PDF or sparse matrix", .dsize=-1, .estimate = estim, .draw = draw, .p=pmf_p, .print=pmf_print};
