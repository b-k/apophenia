/** \file apop_histogram.c 

This implements a one-d histogram representing an empirical
distribution. It is primarily a wrapper for the GSL's comparable functions
in the standard \c apop_model form, for easy comparison with other models.

(c) 2007 Ben Klemens. Licensed under the GNU GPL v 2. */

#include <apophenia/model.h>
#include <apophenia/types.h>
#include <gsl/gsl_histogram.h>
#include <stdio.h>
#include <assert.h>
#include <apop.h>

apop_model apop_histogram;

/** Allocate the parameters for the \c apop_histogram model.

  This produces the \c apop_model half of the \c apop_model/apop_params
  pair. 

  \param    data The input data. As with other distributions, the data
  should be in the matrix element of the \c apop_data set, and can have any dimensions
  ($1\times 10000$, $10000\times 1$, $100\times 100$...).
  \param bins How many bins should the PDF have?
  \param params_in If you want to hook the output parameters in to an existing \c apop_model set, then put that input set here.

 */
apop_model *apop_histogram_params_alloc(apop_data *data, int bins, apop_model *model_in){
    //header is in model.h.
  apop_model *pin  = model_in;
  apop_histogram_params *hp = malloc(sizeof(*hp));
    hp->model  = apop_model_copy(apop_histogram);
    hp->model->model_params        = hp;
    hp->model->model_params_size   = sizeof(*hp);
    if (!model_in){
        pin     = hp->model;
    } else
        hp->model  = pin;
    snprintf(pin->method_name,100, "Histogram");
  size_t              i, j, sum = 0;
  double              minv    = GSL_POSINF,
                      maxv    = GSL_NEGINF,
                      minm    = GSL_POSINF,
                      maxm    = GSL_NEGINF;
    if (data->vector) gsl_vector_minmax(data->vector, &minv, &maxv);
    if (data->matrix) gsl_matrix_minmax(data->matrix, &minm, &maxm);
    hp->pdf = gsl_histogram_alloc(bins);
    gsl_histogram_set_ranges_uniform(hp->pdf, GSL_MIN(minv,minm), GSL_MAX(maxv,maxm));
   //add infinity bins.
    double  *newbins = malloc(sizeof(double)* ( bins +3));
    newbins[0]                  = GSL_NEGINF;
    memcpy((newbins + 1), hp->pdf->range, sizeof(double) * (bins+1));
    newbins[bins+1]      += 2*GSL_DBL_EPSILON; //so max won't fall in the infinity bin.
    newbins[bins+2]         = GSL_POSINF;
    gsl_histogram_free(hp->pdf);
    hp->pdf = gsl_histogram_alloc(bins+2);
    gsl_histogram_set_ranges(hp->pdf, newbins,bins+3);

    if (data->vector)
        for (i=0; i< data->vector->size; i++){
            gsl_histogram_increment(hp->pdf,gsl_vector_get(data->vector,i));
            sum++;
        }
    if (data->matrix)
        for (i=0; i< data->matrix->size1; i++)
            for (j=0; j< data->matrix->size2; j++){
            gsl_histogram_increment(hp->pdf,gsl_matrix_get(data->matrix,i,j));
            sum ++;
        }
    for (i=0; i< hp->pdf->n; i++)
        hp->pdf->bin[i]    /= (sum + 0.0);
    hp->cdf =NULL;
    return hp->model;
}

gsl_histogram *gpdf;

apop_model *est(apop_data *d, apop_model *in){
    return apop_histogram_params_alloc(d, 1000, in);

}

static double one_vector(gsl_vector *in){
  size_t    i, k;
  double    product = 0;
    for (i=0; i< in->size; i++){
        gsl_histogram_find(gpdf, gsl_vector_get(in, i), &k);
        product += gsl_histogram_get (gpdf, k);
    }
    return product;
}

static double histogram_p(const apop_data *beta, apop_data *d, apop_model *parameters){
  apop_histogram_params *hp = parameters->model_params;
  long double           product = 0;
    gpdf    = hp->pdf;
    if (d->vector)
        product += one_vector(d->vector);
    if (d->matrix){
        gsl_vector *outp = apop_matrix_map(d->matrix, one_vector);
        product += apop_vector_sum(outp);
        gsl_vector_free(outp);
    }
    return product;
}

static double histogram_log_likelihood(const apop_data *beta, apop_data *d, apop_model *parameters){
	return log(histogram_p(beta,d,parameters)) ;
}

static void histogram_rng(double *out, gsl_rng *r, apop_model* eps){
  apop_histogram_params *hp   = eps->model_params;
    if (!hp->cdf){
        hp->cdf = gsl_histogram_pdf_alloc(hp->pdf->n); //darn it---this produces a CDF!
        gsl_histogram_pdf_init(hp->cdf, hp->pdf);
    }
    do {
        *out  = gsl_histogram_pdf_sample(hp->cdf, gsl_rng_uniform(r));
    } while (!gsl_finite(*out));
}

/** The histogram model.

  This is an empirical distribution. If you have a data set from which you want to make random draws, this is overkill; instead just use something like \code 
  gsl_rng *r = apop_rng_alloc(27);
  gsl_vector *my_data = [gather data here.];
  gsl_vector_get(my_data, gsl_rng_uniform(r)*my_data->size);
  \endcode

  But this can be used anywhere a model is needed, such as the inputs and outputs to \c apop_update.

  The model is unlike most other models in that there are no parameters
  of any sort (beyond the data itself), so there is no \c estimate
  method; instead all the work of producing the histogram is done in
  \c apop_histogram_params_alloc.

\ingroup models
*/
apop_model apop_histogram = {"histogram", 0,0,0, .estimate = est, 
    .p = histogram_p, .log_likelihood = histogram_log_likelihood, 
    .draw = histogram_rng};


//By the way, here's the null model.
apop_model apop_null = {"The null model", 0,0,0};

