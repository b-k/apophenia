/** \file apop_histogram.c 

This implements a one-d histogram representing an empirical distribution. It is primarily a wrapper for the GSL's comparable functions in the standard \c apop_model form, for easy comparison with other models.*/
/* Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "mapply.h"
#include "settings.h"
#include "variadic.h"
#include <gsl/gsl_math.h>

apop_model apop_histogram;

/** Allocate the parameters for the \c apop_histogram model.

  \param    data The input data. I'll use all data in both the the matrix and vector element of the \c apop_data set, and the matrix can have any dimensions (\f$1\times 10000\f$, \f$10000\times 1\f$, \f$100\times 100\f$...).
  \param bins How many bins should the PDF have?
 */
apop_histogram_settings *apop_histogram_settings_alloc(apop_data *data, int bins){
    //header is in model.h.
  Apop_assert(data, NULL, 0, 'c', "You asked me to set up a histogram with NULL data. Returning a NULL settings group.");
  apop_histogram_settings *hp  = malloc(sizeof(apop_histogram_settings));
  size_t              i, j;
  double              minv    = GSL_POSINF,
                      maxv    = GSL_NEGINF,
                      minm    = GSL_POSINF,
                      maxm    = GSL_NEGINF;
    if (data->vector) gsl_vector_minmax(data->vector, &minv, &maxv);
    if (data->matrix) gsl_matrix_minmax(data->matrix, &minm, &maxm);
    hp->data = data;
    hp->pdf = gsl_histogram_alloc(bins);
    gsl_histogram_set_ranges_uniform(hp->pdf, GSL_MIN(minv,minm), GSL_MAX(maxv,maxm));
   //add infinity bins.
    double  *newbins = malloc(sizeof(double)* (bins + 3));
    newbins[0]                  = GSL_NEGINF;
    memcpy((newbins + 1), hp->pdf->range, sizeof(double) * (bins+1));
    newbins[bins+1]      += 2*GSL_DBL_EPSILON; //so max won't fall in the infinity bin.
    newbins[bins+2]         = GSL_POSINF;
    gsl_histogram_free(hp->pdf);
    hp->pdf = gsl_histogram_alloc(bins+2);
    gsl_histogram_set_ranges(hp->pdf, newbins,bins+3);

    if (data->vector)
        for (i=0; i< data->vector->size; i++)
            gsl_histogram_increment(hp->pdf, apop_data_get(data, i, -1));
    if (data->matrix)
        for (i=0; i< data->matrix->size1; i++)
            for (j=0; j< data->matrix->size2; j++)
                gsl_histogram_increment(hp->pdf, apop_data_get(data, i, j));
    hp->cdf        = NULL;
    hp->histobase  = NULL;
    hp->kernelbase = NULL;
    return hp;
}

/*
apop_histogram_settings * apop_histogram_settings_init(apop_histogram_settings in){
    apop_assert(in.data && in.bins, NULL, 0, 's', "I need both the .data and .bins elements to be set.");
    return apop_histogram_settings_alloc(in.data, in.bins);
}
*/

void * apop_histogram_settings_copy(apop_histogram_settings *in){
    apop_histogram_settings *out = malloc(sizeof(apop_histogram_settings));
    out->pdf = gsl_histogram_clone(in->pdf);
    out->cdf = NULL; //the GSL doesn't provide a copy function, so screw it---just regenerate.
    out->histobase  = in->histobase  ? apop_model_copy(*in->histobase)  : NULL;
    out->kernelbase = in->kernelbase ? apop_model_copy(*in->kernelbase) : NULL;
    return out;
}

void apop_histogram_settings_free(apop_histogram_settings *in){
    gsl_histogram_free(in->pdf);
    gsl_histogram_pdf_free(in->cdf);
    //Assume kernel base and histogram base are freed elsewhere.
}


apop_model *est(apop_data *d, apop_model *in){
    apop_model *out = apop_model_copy(*in);
    if (!(apop_settings_get_group(in, "apop_histogram") || apop_settings_get_group(in, "apop_kernel_density")))
        Apop_settings_add_group(out, apop_histogram, d, 1000);
    return out;
}

static double one_histo_ll(double i, void *gpdf){
  size_t    k;
    gsl_histogram_find(gpdf, i, &k);
    return log(gsl_histogram_get (gpdf, k));
}

static double histogram_ll(apop_data *d, apop_model *in){
    apop_histogram_settings *hp = apop_settings_get_group(in, "apop_histogram");
    if (!hp) apop_settings_get_group(in, "apop_kernel_density");
    apop_assert(hp, 0, 0, 's', "you sent me an unparametrized model.");
    return apop_map_sum(d, .fn_dp =one_histo_ll, .param=hp->pdf);
}

static void histogram_rng(double *out, gsl_rng *r, apop_model* in){
  apop_histogram_settings *hp = apop_settings_get_group(in, "apop_histogram");
  if (!hp) apop_settings_get_group(in, "apop_kernel_density");
  apop_assert_void(hp, 0, 's', "you sent me an unparametrized model.");
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
  method; instead all the work of producing the histogram is done in \c
  apop_histogram_settings_alloc. [Actually, there is an \c estimate method,
  but it is just an alias for the histogram_alloc function.]

\hideinitializer
\ingroup models
*/
apop_model apop_histogram = {"Histogram", .estimate = est, .log_likelihood = histogram_ll, .draw = histogram_rng};


////Kernel density estimation


static void apop_set_first_params(double in, apop_model *m){
    m->parameters->vector->data[0]  = in;
}

static double apop_getmidpt(const gsl_histogram *pdf, const size_t n){
    if (!n) 
        return pdf->range[1];
    if (n== pdf->n-1) 
        return pdf->range[n-2];
    //else
        return (pdf->range[n] + pdf->range[n+1])/2.;
}

static gsl_histogram *apop_alloc_wider_range(const gsl_histogram *in, const double padding){
//put 10% more bins on either side of the original data.
  size_t  newsize =in->n * (1+2*padding);
  double  newrange[newsize];
  double  diff = in->range[2] - in->range[1];
    memcpy(&newrange[(int)(in->n * padding)], in->range, sizeof(double)*in->n);
    for (int k = in->n *padding +1; k>=1; k--)
        newrange[k] = newrange[k+1] - diff;
    for (int k = in->n *(1+padding); k <newsize-1; k++)
        newrange[k] = newrange[k-1] + diff;
    newrange[0]          = GSL_NEGINF;
    newrange[newsize -1] = GSL_POSINF;
    gsl_histogram * out = gsl_histogram_alloc(newsize-1);
    gsl_histogram_set_ranges(out, newrange, newsize);
    return out;
}

/** Allocate and fill a kernel density, which is a smoothed histogram. 

You may either provide a histogram and a \c NULL data set, or a \c NULL histogram and a real data set, in which case I will convert the data set into a histogram and use the histogram thus created.

\param data    a data set, which, if  not \c NULL and \c !histobase , will be converted to a histogram.
  \param histobase This is the preferred format for input data. It is the histogram to be smoothed.
\param kernelbase The kernel to use for smoothing, with all parameters set and a \c p method. Popular favorites are \ref apop_normal and \ref apop_uniform.
\param set_params A function that takes in a single number and the model, and sets the parameters accordingly. The function will call this for every point in the data set. Below is the default, which is used if this is \c NULL. It simply sets the first element of the model's parameter vector to the input number; this is appropriate for a Normal distribution, where we want to center the distribution on each data point in turn.

\code
void apop_set_first_params(double in, apop_model *m){
    m->parameters->vector->data[0]  = in;
}
\endcode

*/
apop_histogram_settings *apop_kernel_density_settings_alloc(apop_data *data, 
        apop_model *histobase, apop_model *kernelbase, void (*set_params)(double, apop_model*)){
  apop_data *smallset  = apop_data_alloc(0,1,1);
  double    padding    = 0.1;

  apop_model *base     = NULL;
  apop_histogram_settings *out = NULL;
    //establish and copy the base histogram
    if(apop_settings_get_group(histobase, "apop_histogram")){
        base = histobase;
    } else if (data){
        base = apop_model_copy(apop_histogram);
        Apop_settings_add_group(base, apop_histogram, data, 1000);
    } else
        apop_error(0, 's', "I need either a histobase model with a histogram or a non-NULL data set.");

    apop_histogram_settings *bh    = apop_settings_get_group(base, "apop_histogram");
    out             = apop_histogram_settings_copy(bh);
    out->pdf        = apop_alloc_wider_range(bh->pdf, padding);
    out->kernelbase = apop_model_copy(*kernelbase);
    out->histobase  = base;
    set_params      = set_params ? set_params : apop_set_first_params;

    //finally, the double-loop producing the density.
    for (size_t i=0; i< bh->pdf->n; i++)
        if (bh->pdf->bin[i]){
            set_params(apop_getmidpt(bh->pdf,i), out->kernelbase);
            for (size_t j=1; j < out->pdf->n-1; j ++){
                smallset->matrix->data[0] = apop_getmidpt(out->pdf,j);
                out->pdf->bin[j] += bh->pdf->bin[i] * 
                        out->kernelbase->p(smallset, out->kernelbase);
            }
        }

    //set end-bins. As you can see, I've commented this part out for now.
    out->pdf->bin[0]    = 0;
    out->pdf->bin[out->pdf->n-1]  = 0;
        /*
    double ratio = out->pdf->bin[1]/(out->pdf->bin[1] + out->pdf->bin[out->pdf->n-2]);
    if (gsl_isnan(ratio)){ //then both bins are zero.
        out->pdf->bin[0]    = 0;
        out->pdf->bin[out->pdf->n-1]  = 0;
    } else {
        out->pdf->bin[0]    = (1-sum)*ratio;
        out->pdf->bin[out->pdf->n-1]  = (1-sum)*(1-ratio);
    }*/

    //normalize to one.
    double sum = 0;
    for (size_t j=1; j < out->pdf->n-1; j ++)
        sum += out->pdf->bin[j];
    for (size_t j=1; j < out->pdf->n-1; j ++)
        out->pdf->bin[j]   /= sum;

    apop_data_free(smallset);
    return out;
}

static apop_model * apop_kernel_density_estimate(apop_data * data,  apop_model *parameters){
    apop_model *m   = apop_model_set_parameters(apop_normal, 0., 1.);
    apop_model *h   = NULL;
    if (!(h = apop_settings_get_group(parameters, "apop_histogram")))
        h = apop_estimate(data, apop_histogram);
    Apop_assert(h, NULL, 0, 's', "I need either a model with a histogram or a non-NULL data set.\n");
    apop_model *out = apop_model_copy(apop_kernel_density);
    Apop_settings_add_group(out, apop_kernel_density, data, h, m, apop_set_first_params);
    return out;
}

/** The apop_kernel_density model.

A Kernel density is simply a smoothing of a histogram. At each point along the histogram, put a distribution (default: Normal(0,1)) on top of the point. Sum all of these distributions to form the output histogram.

The output is a histogram that behaves exactly like the gsl_histogram, except the histobase and kernelbase elements are set.

\hideinitializer
\ingroup models
*/
apop_model apop_kernel_density = {"kernel density estimate",
	.estimate = apop_kernel_density_estimate, .log_likelihood = histogram_ll, .draw = histogram_rng};

