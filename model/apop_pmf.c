#include "asst.h"
#include "model.h"
#include "output.h"
#include "internal.h"
#include "conversions.h"

extern apop_model apop_pmf;

/** For usage, see the documentation for the \ref apop_pmf model. 
  \param d An input crosstab
  \return A PMF model which has a single line for each nonzero line of the crosstab.
 */
apop_model *apop_crosstab_to_pmf (apop_data *d){
    Get_vmsizes(d) //tsize
    int use_matrix=0, use_vector = 0;
    size_t ctr = 0;
    double x;
    if (d->matrix) use_matrix++;
    else if (d->vector) use_vector++;
    else Apop_assert(0, "You gave me an input set with neither vector nor matrix data.");
    apop_model *out = apop_model_copy(apop_pmf);
    out->parameters = apop_data_alloc(0, tsize, (use_matrix ? 2: 1));
    out->parameters->weights = gsl_vector_alloc(tsize);
    out->data = out->parameters;
    if (use_matrix){
        for(size_t i=0; i < d->matrix->size1; i ++)
            for(size_t j=0; j < d->matrix->size2; j ++)
                if ((x = gsl_matrix_get(d->matrix, i, j))) {
                    apop_data_set(out->parameters, ctr, 0, i);
                    apop_data_set(out->parameters, ctr, 1, j);
                    gsl_vector_set(out->parameters->weights, ctr++, x);
                }
    }
    else if (use_vector)
        for(size_t i=0; i < d->vector->size; i++)
            if ((x = gsl_vector_get(d->vector, i))){
                apop_data_set(out->parameters, ctr, 0, i);
                gsl_vector_set(out->parameters->weights, ctr++, x);
            }
    if (ctr){
        apop_vector_realloc(out->parameters->weights, ctr);
        apop_matrix_realloc(out->parameters->matrix, ctr, (use_matrix ? 2: 1));
    }
    return out;
}

static apop_model *estim (apop_data *d, apop_model *out){
    out->data = d;
    apop_data_free(out->parameters); //may have been auto-alloced by prep.
    out->parameters = apop_data_copy(d);
    return out;
}

static void draw (double *out, gsl_rng *r, apop_model *m){
    size_t current; 
    if (!m->data->weights) //all rows are equiprobable
        current = gsl_rng_uniform(r)*
                  ((m->data->vector ? m->data->vector->size : m->data->matrix->size1)-1);
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
    //done searching. current should now be the right row index.
    Apop_data_row(m->data, current, outrow);
    int i = 0;
    if (outrow->vector)
        out[i++] = outrow->vector->data[0];
    if (outrow->matrix)
        for( ; i < outrow->matrix->size2; i ++)
            out[i] = gsl_matrix_get(outrow->matrix, 0, i);
}

double pmf_p(apop_data *d, apop_model *m){
    Nullcheck_p(m) 
    Get_vmsizes(m->data)//firstcol, vsize, vsize1, msize2
    int j;
    double ll = 0;
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
