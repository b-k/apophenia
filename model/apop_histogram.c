/* Histograms via the GSL histogram.
This implements a one-d histogram representing an empirical distribution. It is primarily a wrapper for the GSL's comparable functions in the standard \c apop_model form, for easy comparison with other models.

Copyright (c) 2007, 2010 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  

\amodel apop_histogram A one-dimensional histogram.

  This is an empirical distribution. If you have a data set from which you want to make random draws, this is overkill; instead just use something like \code 
  gsl_rng *r = apop_rng_alloc(27);
  gsl_vector *my_data = [gather data here.];
  gsl_vector_get(my_data, gsl_rng_uniform(r)*my_data->size);
  \endcode

In fact, this is the core of the \ref apop_pmf model, which also works for multidimensional or text data; do consider using that model instead of this one.

But this can be used anywhere a model is needed, such as the inputs and outputs to \c apop_update.

\adoc    Name  <tt>Histogram</tt>
\adoc    Input_format  A free-format pile of scalars, in vector, matrix, or both.        
\adoc    Parameter_format  No parameters; all information is kept in the \ref apop_histogram_settings struct.   
\adoc    estimated_parameters None.  
\adoc    estimated_settings Really just initializes an  \ref apop_histogram_settings struct using some default values. Instead of  calling the \c estimate method, you could also do a call like:
  \code
  int binct = 100; //or some other reasonable number of histogram bins
  Apop_model_add_group(your_model, apop_histogram, .data = your_data_set, .bins_in = binct);
  \endcode
   
  When adding a settings group, the \c .data and \c .bins_in inputs are mandatory.

\adoc    RNG  The first call produces a cumulative density tally,
            and so will take several microseconds longer than later calls.  
\adoc    settings    \ref apop_histogram_settings   
*/

#include "arms.h"
#include "model.h"
#include "mapply.h"
#include "internal.h"
#include "settings.h"
#include "conversions.h"
#include "deprecated.h"
#include "variadic.h"
#include <gsl/gsl_math.h>

apop_histogram_settings *apop_histogram_settings_alloc(apop_data *data, int bins){
    //header is in model.h.
  Apop_assert_c(data, NULL, 0, "You asked me to set up a histogram with NULL data. Returning a NULL settings group.");
  apop_histogram_settings *hp  = malloc(sizeof(apop_histogram_settings));
    size_t i, j;
    double minv    = GSL_POSINF,
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
    return hp;
}

apop_histogram_settings * apop_histogram_settings_init(apop_histogram_settings in){
    Apop_assert(in.data && in.bins_in, "I need both the .data and .bins_in elements to be set.");
    return apop_histogram_settings_alloc(in.data, in.bins_in);
}

Apop_settings_copy(apop_histogram,
    out->pdf = gsl_histogram_clone(in->pdf);
    out->cdf = NULL; //the GSL doesn't provide a copy function, so screw it---just regenerate.
    out->histobase = in->histobase ? apop_model_copy(*in->histobase) : NULL;
    out->ownerbase = in->histobase ? 1 : 0;
)

Apop_settings_free(apop_histogram,
    gsl_histogram_free(in->pdf);
    gsl_histogram_pdf_free(in->cdf);
    if (in->ownerbase)
        apop_model_free(in->histobase);
)


apop_model *est(apop_data *d, apop_model *est){
    if (!(apop_settings_get_group(est, apop_histogram)))
        Apop_settings_add_group(est, apop_histogram, d, 1000);
    return est;
}

static double one_histo_ll(double i, void *gpdf){
    size_t k;
    gsl_histogram_find(gpdf, i, &k);
    return log(gsl_histogram_get (gpdf, k));
}

static double histogram_ll(apop_data *d, apop_model *in){
    Nullcheck_mpd(d, in);
    apop_histogram_settings *hp = apop_settings_get_group(in, apop_histogram);
    if (!hp) apop_settings_get_group(in, apop_kernel_density);
    Apop_assert(hp, "you sent me an unparametrized model.");
    return apop_map_sum(d, .fn_dp =one_histo_ll, .param=hp->pdf);
}

static void histogram_rng(double *out, gsl_rng *r, apop_model* in){
    Nullcheck_m(in);
    apop_histogram_settings *hp = apop_settings_get_group(in, apop_histogram);
    if (!hp) apop_settings_get_group(in, apop_kernel_density);
    Apop_assert(hp, "you sent me an unparametrized model.");
    if (!hp->cdf){
        hp->cdf = gsl_histogram_pdf_alloc(hp->pdf->n); //darn it---this produces a CDF!
        gsl_histogram_pdf_init(hp->cdf, hp->pdf);
    }
    do {
        *out  = gsl_histogram_pdf_sample(hp->cdf, gsl_rng_uniform(r));
    } while (!gsl_finite(*out));
}

static void histogram_print(apop_model *est){
    gsl_histogram *h = Apop_settings_get(est, apop_histogram, pdf);
    gsl_histogram_fprintf (apop_opts.output_pipe, h, "%g", "%g");
}

apop_model apop_histogram = {"Histogram", .dsize=1, .estimate = est, .log_likelihood = histogram_ll, 
                    .draw = histogram_rng, .print=histogram_print};

/*\amodel apop_kernel_density The kernel density smoothing of a PMF

A Kernel density is simply a smoothing of a histogram. At each point
along the histogram, put a distribution (default: Normal(0,1)) on top
of the point. Sum all of these distributions to form the output histogram.

Elements of \ref apop_kernel_density_settings that you may want to set:

\li data a data set, which, if  not \c NULL and \c !base_pmf , will be converted to an \ref apop_pmf model.
\li base_pmf This is the preferred format for input data. It is the histogram to be smoothed.
\li kernelbase The kernel to use for smoothing, with all parameters set and a \c p method. Popular favorites are \ref apop_normal and \ref apop_uniform.
\li set_params A function that takes in a single number and the model, and sets
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

\adoc    Input_format  I'll estimate a \ref apop_pmf internally, so I
                      follow that format, which is one observation (of any format) per line.        
\adoc    Parameter_format  None    
\adoc    Estimated_parameters None
\adoc    Estimated_settings  The estimate method basically just runs
                             <tt>apop_model_add_group(your_data, apop_kernel_density);</tt>    
\adoc    settings  \ref apop_kernel_density_settings.    
\adoc    Examples
This example sets up and uses KDEs based on a Normal and a Uniform distribution.

\include kernel.c
*/

static void apop_set_first_param(apop_data *in, apop_model *m){
    m->parameters->vector->data[0]  = in->vector ?
            in->vector->data[0] : gsl_matrix_get(in->matrix, 0, 0);
}

Apop_settings_init(apop_kernel_density, 
    //If there's a PMF associated with the model, run with it.
    //else, generate one from the data.
    Apop_varad_set(base_pmf, apop_estimate(in.base_data, apop_pmf));
    Apop_varad_set(kernel, apop_model_set_parameters(apop_normal, 0, 1));
    Apop_varad_set(set_fn, apop_set_first_param);
    out->own_pmf = !in.base_pmf;
    out->own_kernel = 1;
    if (!out->kernel->parameters)
        apop_prep(out->base_data, out->kernel);
)

Apop_settings_copy(apop_kernel_density,
    out->own_pmf    =
    out->own_kernel = 0;
)

Apop_settings_free(apop_kernel_density,
    if (in->own_pmf)
        apop_model_free(in->base_pmf);
    if (in->own_kernel)
        apop_model_free(in->kernel);
)

static apop_model *apop_kernel_estimate(apop_data *d, apop_model *m){
    Nullcheck_d(d);
    if (!apop_settings_get_group(m, apop_kernel_density))
        apop_model_add_group(m, apop_kernel_density, .base_data=d);
    return m;
}

static double kernel_p_cdf_base(apop_data *d, apop_model *m,
        double (*fn)(apop_data*,apop_model*)){
    Nullcheck_mpd(d, m);
    Get_vmsizes(d);
    long double total = 0;
    apop_kernel_density_settings *ks = apop_settings_get_group(m, apop_kernel_density);
    apop_data *pmf_data = apop_settings_get(m, apop_kernel_density, base_pmf)->data;
    int len = pmf_data->weights ? pmf_data->weights->size
                                : pmf_data->vector ? pmf_data->vector->size
                                                   : pmf_data->matrix->size1;
    for (int k = 0; k < len; k++){
        Apop_data_row(pmf_data, k, r);
        (ks->set_fn)(r, ks->kernel);
        total += fn(d, ks->kernel);
    }
    double weight = pmf_data->weights ? apop_sum(pmf_data->weights) : len;
    total /= weight;
    return total;
}

static double kernel_p(apop_data *d, apop_model *m){
    return kernel_p_cdf_base(d, m, apop_p);
}

/* \adoc    CDF Sums the CDF to the given point of all the sub-distributions.*/
static double kernel_cdf(apop_data *d, apop_model *m){
    return kernel_p_cdf_base(d, m, apop_cdf);
}

/* \adoc    RNG  Randomly selects a data point, then randomly draws from that sub-distribution.*/
static void kernel_draw(double *d, gsl_rng *r, apop_model *m){
    //randomly select a point, using the weights.
    apop_kernel_density_settings *ks = apop_settings_get_group(m, apop_kernel_density);
    apop_model *pmf = apop_settings_get(m, apop_kernel_density, base_pmf);
    apop_data *point = apop_data_alloc(1, pmf->dsize);
    Apop_row(point, 0, draw_here);
    apop_draw(draw_here->data, r, pmf);
    (ks->set_fn)(point, ks->kernel);
    //Now draw from the distribution around that point.
    apop_draw(d, r, ks->kernel);
    apop_data_free(point);
}

apop_model apop_kernel_density = {"kernel density estimate", .dsize=1,
	.estimate = apop_kernel_estimate, .p = kernel_p, .cdf=kernel_cdf, .draw=kernel_draw};
