/** \file apop_kernel.c 

  A Kernel density is simply a smoothing of a histogram. At each point
  along the histogram, put a distribution (default: Normal(0,1)) on
  top of the point. Sum all of these distributions to form the output
  histogram.

  The output is a histogram that behaves exactly like the gsl_histogram,
  except the histobase and kernelbase elements are set.

Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.
*/

#include <apophenia/model.h>
#include <apophenia/types.h>
#include <apophenia/histogram.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <stdio.h>
#include <assert.h>


static void apop_set_first_params(double in, apop_model *m){
    m->parameters->vector->data[0]  = in;
}

void apop_histogram_plot(apop_model *in, char *outfile){
  int             i, k;
  FILE *          f;
  apop_histogram_params   *inhist  = in->model_settings;
  gsl_histogram   *h      = inhist->pdf;

    double midpoints[h->n]; //cut 'n' pasted from kernel density alloc.
    midpoints[0]            = h->range[1];
    midpoints[h->n]  = h->range[h->n-1];
    for(k=1; k < h->n-1; k ++)
        midpoints[k] = (h->range[k+1] +h->range[k])/2.;

    if (apop_opts.output_type == 'p')
        f   = apop_opts.output_pipe;
    else 
        f = (outfile ? fopen(outfile, "a") : stdout);
	fprintf(f, "set key off					                    ;\n\
                        plot '-' with lines\n");
	for (i=1; i < h->n-1; i++)
	    fprintf(f, "%4f\t %g\n", midpoints[i], gsl_histogram_get(h, i));
	fprintf(f, "e\n");
    if (apop_opts.output_type == 'p')
        fflush(f);
    else if (outfile)    
        fclose(f);
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
  int     k;
    memcpy(&newrange[(int)(in->n * padding)], in->range, sizeof(double)*in->n);
    for (k = in->n *padding +1; k>=1; k--)
        newrange[k] = newrange[k+1] - diff;
    for (k = in->n *(1+padding); k <newsize-1; k++)
        newrange[k] = newrange[k-1] + diff;
    newrange[0]          = GSL_NEGINF;
    newrange[newsize -1] = GSL_POSINF;
    gsl_histogram * out = gsl_histogram_alloc(newsize-1);
    gsl_histogram_set_ranges(out, newrange, newsize);
    return out;
}

/** Allocate and fill a kernel density, which is a smoothed histogram. 


  You may either provide a histogram and a \c NULL data set, or a \c
  NULL histogram and a real data set, in which case I will convert the
  data set into a histogram and use the histogram thus created.

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
apop_model *apop_kernel_density_params_alloc(apop_data *data, 
        apop_model *histobase, apop_model *kernelbase, void (*set_params)(double, apop_model*)){
  size_t    i, j;
  apop_data *smallset          = apop_data_alloc(0,1,1);
  apop_histogram_params *out   = malloc(sizeof(apop_histogram_params));
    set_params                 = set_params ? set_params : apop_set_first_params;
    out->model                 = apop_model_copy(apop_kernel_density);
    out->model->model_settings = out;
    out->kernelbase            = apop_model_copy(*kernelbase);
    out->histobase = data && !histobase ?
                apop_histogram_params_alloc(data, 1000)
                : apop_model_copy(*histobase);
  apop_histogram_params *bh    = out->histobase->model_settings;

    double  padding = 0.1;
    out->pdf        = apop_alloc_wider_range(bh->pdf, padding);

    //finally, the double-loop producing the density.
    for (i=0; i< bh->pdf->n; i++)
        if(bh->pdf->bin[i]){
            set_params(apop_getmidpt(bh->pdf,i), out->kernelbase);
            for(j=1; j < out->pdf->n-1; j ++){
                smallset->matrix->data[0] = apop_getmidpt(out->pdf,j);
                out->pdf->bin[j] += bh->pdf->bin[i] * 
                        out->kernelbase->p(smallset, out->kernelbase);
            }
        }
    //normalize
    double sum = 0;
    for(j=1; j < out->pdf->n-1; j ++){
        sum +=
        out->pdf->bin[j]   /= bh->pdf->n;
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
    apop_data_free(smallset);
    apop_histogram_normalize(out->model);
    return out->model;
}

static apop_model * apop_kernel_density_estimate(apop_data * data,  apop_model *parameters){
    apop_model *m   = apop_model_set_parameters(apop_normal, 0., 1.);
    apop_model *h   = apop_estimate(data, apop_histogram);
	return apop_kernel_density_params_alloc(data, h, m, apop_set_first_params);
}

static double apop_kernel_density_log_likelihood(apop_data *d, apop_model *p){
    return log(apop_histogram.p(d,p));
}

static double apop_kernel_density_p(apop_data *d, apop_model *p){
    return apop_histogram.p(d, p);
}

static void apop_kernel_density_rng( double *out, gsl_rng *r, apop_model* eps){
    apop_histogram.draw(out, r, eps);
}

/** The apop_kernel_density model.
Takes in a histogram, and smooths it out via the kernel density method.


\ingroup models
*/
apop_model apop_kernel_density = {"kernel density estimate", 1,0,0,
	.estimate = apop_kernel_density_estimate, .p = apop_kernel_density_p,
    .log_likelihood = apop_kernel_density_log_likelihood,
     .draw = apop_kernel_density_rng};
