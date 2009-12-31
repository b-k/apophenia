#include "asst.h"
#include "model.h"
#include "internal.h"

extern apop_model apop_pmf;

static apop_model *estim (apop_data *d, apop_model *m){
    Get_vmsizes(d) //tsize
    int use_matrix=0, use_vector = 0;
    size_t ctr = 0;
    double x;
    apop_model *out = apop_model_copy(*m);
    if (d->matrix) use_matrix++;
    else if (d->vector) use_vector++;
    else apop_error(0, 's', "You gave me an input set with neither vector nor matrix data.\n");
    out->parameters = apop_data_alloc(0, tsize, (use_matrix ? 2: 1));
    out->parameters->weights = gsl_vector_alloc(tsize);
    out->more = NULL;
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

/*
//generate a CDF
static apop_model *estim (apop_data *d, apop_model *m){
    apop_model *out = apop_model_copy(apop_pmf);
    apop_assert(d->weights, NULL, 0, 's', "I expect the input to estimation of "
               "an apop_pmf to have a weights vector.\n");
    int size = d->weights->size;
    apop_assert(size, NULL, 0, 's', "I expect the input to estimation of "
               "an apop_pmf to have a weights vector with positive length.\n");
    gsl_vector *cdf = gsl_vector_alloc(size);
    out->more = cdf;
    cdf->data[0] = d->weights->data[0];
    for (int i=1; i< d->weights->size; i++)
        cdf->data[i] = d->weights->data[i] + cdf->data[i-1];
    //Now make sure the last entry is one.
    gsl_vector_scale(cdf, 1./cdf->data[cdf->size-1]);
    return out;
}
*/

static void draw (double *out, gsl_rng *r, apop_model *m){
    size_t size = m->parameters->weights->size;
    if (!m->more){
        m->more = gsl_vector_alloc(size);
        gsl_vector *cdf = m->more; //alias.
        cdf->data[0] = m->parameters->weights->data[0];
        for (int i=1; i< size; i++)
            cdf->data[i] = m->parameters->weights->data[i] + cdf->data[i-1];
        //Now make sure the last entry is one.
        gsl_vector_scale(cdf, 1./cdf->data[size-1]);
    }
    double draw = gsl_rng_uniform(r);
    //do a binary search for where draw is in the CDF.
    double *cdf = ((gsl_vector*)m->more)->data; //alias.
    size_t top = size, bottom = 0, current = (top+bottom)/2.;
    if (current==0){//array of size one or two
        if (size!=0) 
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
    //done searching. current should now be the right row index.
    Apop_row(m->parameters, current, outrow);
    for(size_t i=0; i < outrow->size; i ++)
        out[i] = gsl_vector_get(outrow, i);
}

apop_model apop_pmf = {"PDF or sparse matrix", .estimate = estim, .draw = draw };
