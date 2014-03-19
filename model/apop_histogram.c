/* Histograms via the GSL histogram.
This implements a one-d histogram representing an empirical distribution. It is primarily a wrapper for the GSL's comparable functions in the standard \c apop_model form, for easy comparison with other models.

Copyright (c) 2007, 2010, 2013 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
*/

#include "apop_internal.h"
#include <gsl/gsl_math.h>


/*\amodel apop_kernel_density The kernel density smoothing of a PMF or histogram.

At each point along the histogram, put a distribution (default: Normal(0,1)) on top
of the point. Sum all of these distributions to form the output distribution.

Elements of \ref apop_kernel_density_settings that you may want to set:

\li \c data a data set, which, if  not \c NULL and \c !base_pmf , will be converted to an \ref apop_pmf model.
\li \c base_pmf This is the preferred format for input data. It is the histogram to be smoothed.
\li \c kernelbase The kernel to use for smoothing, with all parameters set and a \c p method. Popular favorites are \ref apop_normal and \ref apop_uniform.
\li \c set_params A function that takes in a single number and the model, and sets
the parameters accordingly. The function will call this for every point in the data
set. Here is the default, which is used if this is \c NULL. It simply sets the first
element of the model's parameter vector to the input number; this is appropriate for a
Normal distribution, where we want to center the distribution on each data point in turn.

\code
static void apop_set_first_param(apop_data *in, apop_model *m){
    m->parameters->vector->data[0]  = apop_data_get(in);
}
\endcode

For a Uniform[0,1] recentered around the first element of the PMF matrix, you could put this function in your code:

\code
void set_midpoint(apop_data * in, apop_model *m){
    apop_data_set(m->parameters, 0, -1, apop_data_get(in)+0.5);
    apop_data_set(m->parameters, 1, -1, apop_data_get(in)-0.5);
}
\endcode

\adoc Input_format  I'll estimate a \ref apop_pmf internally, so I
                   follow that format, which is one observation (of any format) per line.
\adoc Parameter_format  None
\adoc Estimated_parameters None
\adoc Estimated_settings  The estimate method basically just runs
                          <tt>apop_model_add_group(your_data, apop_kernel_density);</tt>
\adoc Settings  \ref apop_kernel_density_settings
\adoc Examples
This example sets up and uses KDEs based on a Normal and a Uniform distribution.

\include kernel.c
*/

static void apop_set_first_param(apop_data *in, apop_model *m){
    m->parameters->vector->data[0] = in->vector ? in->vector->data[0] 
                                                 : gsl_matrix_get(in->matrix, 0, 0);
}

Apop_settings_init(apop_kernel_density, 
    //If there's a PMF associated with the model, run with it.
    //else, generate one from the data.
    Apop_varad_set(base_pmf, apop_estimate(in.base_data, apop_pmf));
    Apop_varad_set(kernel, apop_model_set_parameters(apop_normal, 0, 1));
    Apop_varad_set(set_fn, apop_set_first_param);
    out->own_pmf = !in.base_pmf;
    out->own_kernel = !in.kernel;
    if (!out->kernel->parameters) apop_prep(out->base_data, out->kernel);
)

Apop_settings_copy(apop_kernel_density,
    out->own_pmf    =
    out->own_kernel = 0;
)

Apop_settings_free(apop_kernel_density,
    if (in->own_pmf)    apop_model_free(in->base_pmf);
    if (in->own_kernel) apop_model_free(in->kernel);
)

static void apop_kernel_estimate(apop_data *d, apop_model *m){
    Nullcheck_d(d, );
    if (!apop_settings_get_group(m, apop_kernel_density))
        Apop_settings_add_group(m, apop_kernel_density, .base_data=d);
}

static long double kernel_p_cdf_base(apop_data *d, apop_model *m,
        double (*fn)(apop_data*,apop_model*)){
    Nullcheck_d(d, GSL_NAN);
    Nullcheck_m(m, GSL_NAN);
    long double total = 0;
    apop_kernel_density_settings *ks = apop_settings_get_group(m, apop_kernel_density);
    apop_data *pmf_data = apop_settings_get(m, apop_kernel_density, base_pmf)->data;
    Get_vmsizes(pmf_data); //maxsize
    for (size_t k = 0; k < maxsize; k++){
        Apop_row(pmf_data, k, r);
        double wt = r->weights ? *r->weights->data : 1;
        (ks->set_fn)(r, ks->kernel);
        total += fn(d, ks->kernel)*wt;
    }
    long double weight = pmf_data->weights ? apop_sum(pmf_data->weights) : maxsize;
    total /= weight;
    return total;
}

static long double kernel_p(apop_data *d, apop_model *m){
    return kernel_p_cdf_base(d, m, apop_p);
}

/* \adoc    CDF Sums the CDF to the given point of all the sub-distributions.*/
static long double kernel_cdf(apop_data *d, apop_model *m){
    return kernel_p_cdf_base(d, m, apop_cdf);
}

/* \adoc    RNG  Randomly selects a data point, then randomly draws from that sub-distribution.
 Returns 0 on success, 1 if unable to pick a sub-distribution (meaning the weights over the distributions are somehow broken), and 2 if unable to draw from the sub-distribution.
 */
static int kernel_draw(double *d, gsl_rng *r, apop_model *m){
    //randomly select a point, using the weights.
    apop_kernel_density_settings *ks = apop_settings_get_group(m, apop_kernel_density);
    apop_model *pmf = apop_settings_get(m, apop_kernel_density, base_pmf);
    apop_data *point = apop_data_alloc(1, pmf->dsize);
    Apop_row_v(point, 0, draw_here);
    Apop_stopif(apop_draw(draw_here->data, r, pmf), return 1, 0, "Unable to use the PMF over kernels to select a kernel from which to draw.");
    (ks->set_fn)(point, ks->kernel);
    //Now draw from the distribution around that point.
    Apop_stopif(apop_draw(d, r, ks->kernel), return 2, 0, "unable to draw from a single selected kernel.");
    apop_data_free(point);
    return 0;
}

apop_model *apop_kernel_density = &(apop_model){"kernel density estimate", .dsize=1,
	.estimate = apop_kernel_estimate, .p = kernel_p, .cdf=kernel_cdf, .draw=kernel_draw};
