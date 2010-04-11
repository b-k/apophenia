/** \file apop_histogram.c 

This implements a one-d histogram representing an empirical distribution. It is primarily a wrapper for the GSL's comparable functions in the standard \c apop_model form, for easy comparison with other models.*/
/* Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "mapply.h"
#include "internal.h"
#include "settings.h"
#include "conversions.h"
#include "deprecated.h"
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

/** Initialize an  \ref apop_histogram_settings struct. You'll probably call this via
  \code
  int binct = 100; //or some other reasonable number of histogram bins
  Apop_model_add_group(your_model, apop_histogram, .data = your_data_set, .bins_in = binct);
  \endcode
   
  The \c .data input is mandatory.
*/
apop_histogram_settings * apop_histogram_settings_init(apop_histogram_settings in){
    apop_assert(in.data && in.bins_in, NULL, 0, 's', "I need both the .data and .bins_in elements to be set.");
    return apop_histogram_settings_alloc(in.data, in.bins_in);
}

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


apop_model *est(apop_data *d, apop_model *est){
    if (!(apop_settings_get_group(est, apop_histogram) || apop_settings_get_group(est, apop_kernel_density)))
        Apop_settings_add_group(est, apop_histogram, d, 1000);
    return est;
}

static double one_histo_ll(double i, void *gpdf){
  size_t    k;
    gsl_histogram_find(gpdf, i, &k);
    return log(gsl_histogram_get (gpdf, k));
}

static double histogram_ll(apop_data *d, apop_model *in){
    apop_histogram_settings *hp = apop_settings_get_group(in, apop_histogram);
    if (!hp) apop_settings_get_group(in, apop_kernel_density);
    apop_assert(hp, 0, 0, 's', "you sent me an unparametrized model.");
    return apop_map_sum(d, .fn_dp =one_histo_ll, .param=hp->pdf);
}

static void histogram_rng(double *out, gsl_rng *r, apop_model* in){
  apop_histogram_settings *hp = apop_settings_get_group(in, apop_histogram);
  if (!hp) apop_settings_get_group(in, apop_kernel_density);
  apop_assert_void(hp, 0, 's', "you sent me an unparametrized model.");
    if (!hp->cdf){
        hp->cdf = gsl_histogram_pdf_alloc(hp->pdf->n); //darn it---this produces a CDF!
        gsl_histogram_pdf_init(hp->cdf, hp->pdf);
    }
    do {
        *out  = gsl_histogram_pdf_sample(hp->cdf, gsl_rng_uniform(r));
    } while (!gsl_finite(*out));
}

apop_model apop_histogram = {"Histogram", .dsize=1, .estimate = est, .log_likelihood = histogram_ll, .draw = histogram_rng};


////Kernel density estimation

static void apop_set_first_param(apop_data_row in, apop_model *m){
    m->parameters->vector->data[0]  = in.vector_pt ?
            *in.vector_pt : gsl_vector_get(&in.matrix_row, 0);
}

/** Allocate and fill a kernel density, which is a smoothed histogram. 

You may either provide a histogram and a \c NULL data set, or a \c NULL histogram and a real data set, in which case I will convert the data set into a histogram and use the histogram thus created.

\param data    a data set, which, if  not \c NULL and \c !histobase , will be converted to a histogram.
  \param histobase This is the preferred format for input data. It is the histogram to be smoothed.
\param kernelbase The kernel to use for smoothing, with all parameters set and a \c p method. Popular favorites are \ref apop_normal and \ref apop_uniform.
\param set_params A function that takes in a single number and the model, and sets
the parameters accordingly. The function will call this for every point in the data
set. Here is the default, which is used if this is \c NULL. It simply sets the first
element of the model's parameter vector to the input number; this is appropriate for a
Normal distribution, where we want to center the distribution on each data point in turn.

\code
void apop_set_first_param(double in, apop_model *m){
    m->parameters->vector->data[0] = in;
}
\endcode

For a Uniform[0,1] recentered around each point, you'd want to put this function in your code:

\code
void set_midpoint(double in, apop_model *m){
    m->parameters->vector->data[0] = in-0.5;
    m->parameters->vector->data[1] = in+0.5;
}
\endcode

*/

apop_kernel_density_settings *apop_kernel_density_settings_init(apop_kernel_density_settings in){
    //If there's a PMF associated with the model, run with it.
    //else, generate one from the data.
    apop_kernel_density_settings *out = malloc(sizeof(apop_kernel_density_settings));
    apop_varad_setting(in, out, base_pmf, apop_estimate(in.base_data, apop_pmf));
    apop_varad_setting(in, out, kernel, apop_model_set_parameters(apop_normal, 0, 1));
    apop_varad_setting(in, out, set_fn, apop_set_first_param);
    out->own_pmf = !in.base_pmf;
    if (!out->kernel->parameters)
        apop_prep(out->base_data, out->kernel);
    return out;
}

void * apop_kernel_density_settings_copy(apop_kernel_density_settings *in){
    apop_kernel_density_settings *out = malloc(sizeof(apop_kernel_density_settings));
    *out  = *in;
    out->own_pmf    =
    out->own_kernel = 0;
    return out;
}

void apop_kernel_density_settings_free(apop_kernel_density_settings *in){
    if (in->own_pmf)
        apop_model_free(in->base_pmf);
    if (in->own_kernel)
        apop_model_free(in->kernel);
}

apop_model *apop_kernel_estimate(apop_data *d, apop_model *m){
    //Uh, nothing. Just run the init fn.
    if (!apop_settings_get_group(m, apop_kernel_density))
        apop_model_add_group(m, apop_kernel_density, .base_data=d);
    return m;
}

double kernel_p(apop_data *d, apop_model *m){
    Get_vmsizes(d);
    long double p = 0;
    apop_kernel_density_settings *ks = apop_settings_get_group(m, apop_kernel_density);
    apop_data *pmf_data = apop_settings_get(m, apop_kernel_density, base_pmf)->parameters;
    int len = pmf_data->weights ? pmf_data->weights->size
                                : pmf_data->vector ? pmf_data->vector->size
                                                   : pmf_data->matrix->size1;
    for (int k = 0; k < len; k++){
        apop_data_row r = apop_data_get_row(pmf_data, k);
        (ks->set_fn)(r, ks->kernel);
        p += apop_p(d, ks->kernel);
    }
    p /= len;
    return p;
}

apop_model apop_kernel_density = {"kernel density estimate", .dsize=1,
	.estimate = apop_kernel_estimate, .p = kernel_p};
