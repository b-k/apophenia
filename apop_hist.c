/** \file apop_hist.c	PMF and CMF manipulations.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
 (Except psmirnov2x, Copyright R Project, but also licensed under the GPL.)
 \author Ben Klemens
 */


/*
The characters:
    vectors
    hisograms (PDFs)
    histogram_pdfs (CDFs)
    models

Things to do:
    --produce PDFs from vectors
    --produce PDFs from models
    --produce CDFs from PDFs
    --produce CDFs from vectors
    --produce synced PDFs from (v,v), (v,m)
    --produce synced PDFs (from scratch w/shared maxmin, from scratch w/one having infinite range)
    --type I v type II trials.

    --send two vectors to goodness-of-fit test.
    --send vector, model to goodness-of-fit test.
    --send synced pdfs to goodness-of-fit test.


    OK, here are my functions:

    apop_vector_to_pmf()    //incorporate both synced and unsynced.
    apop_model_to_pmf()     //same.
    apop_pmf_to_cmf()
    apop_pmf_goodness_of_fit()

    filters provided:
        vector          -> pmf
        (vector, vector)-> pmf
        (pmf, model)    -> pmf
        pmf             -> cmf
        (pmf, pmf)      -> gof


*/

#include "db.h"     //just for apop_opts
#include "apophenia/stats.h"
#include "apophenia/histogram.h"
#include "apophenia/bootstrap.h" //rng_alloc
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_vector.h>

/** \defgroup histograms Histograms, PMFs, and CMFs

The GSL provides a few structures that basically accumulate data into
bins. The first is the <tt>gsl_histogram</tt> structure, that produces a PMF.

To produce a PMF from a vector, use \ref apop_vector_to_histogram; to produce
a PMF via random draws from a model, use \ref apop_model_to_pmf.


The second structure from the GSL incrementally sums up the PMF's bins to
produce a CMF. The CMF can be used to map from a draw from a Uniform[0,1]
to a draw from the PMF.  Because it can be used to draw from the PMF,
the GSL calls this the <tt>gsl_histogram_pdf</tt> structure. That's right:
the the data in the <tt>gsl_histogram_pdf</tt> structure is a cumulative
sum---a CMF.

Anyway, here are some functions to deal with these various histograms and such.

 */



/** Give me an existing histogram (i.e., an \c apop_model) and I'll
 create a new histogram with the same bins, but with data from the vector you provide 

\param template An \c apop_model produced using a form like \c apop_estimate(yourdata, apop_histogram).
\param indata The new data to be binned.
\ingroup histograms
*/
apop_model *apop_histogram_refill_with_vector(apop_model *template, gsl_vector *indata){
  if (!template || strcmp(template->name, "Histogram"))
      apop_error(0, 's', "%s: The first argument needs to be an apop_histogram model.\n", __func__);
  size_t  i;
    apop_model *out = apop_model_copy(*template); 
    gsl_histogram *hout  = ((apop_histogram_params*) out->model_settings)->pdf;
    gsl_histogram_reset(hout);
    for (i=0; i< indata->size; i++)
        gsl_histogram_increment(hout, gsl_vector_get(indata, i));
    return out;
}


/** Give me an existing histogram (i.e., an \c apop_model) and I'll
 create a new histogram with the same bins, but with data from \c draws
 random draws from the parametrized model you provide.

 Unlike with most other histogram-genrating functions, this one will normalize the output to integrate to one.

\param template An \c apop_model produced using a form like \c apop_estimate(yourdata, apop_histogram).
\param m The model to be drawn from. Because this function works via random draws, the model needs to have a 
\c draw method.
\param draws The number of random draws to make.
\param r The \c gsl_rng used to make random draws.
\ingroup histograms
*/
apop_model *apop_histogram_refill_with_model(apop_model *template, apop_model *m, long int draws, gsl_rng *r){
  if (!template || strcmp(template->name, "Histogram"))
      apop_error(0, 's', "%s: The first argument needs to be an apop_histogram model.\n", __func__);
  if (!m || !m->draw)
      apop_error(0, 's', "%s: The second argument needs to be an apop_model with a function to make radom draws.\n", __func__);
  long double i;
  double d;
    apop_model *out = apop_model_copy(*template); 
    gsl_histogram *hout  = ((apop_histogram_params*) out->model_settings)->pdf;
    gsl_histogram_reset(hout);
    for (i=0; i< draws; i++){
        m->draw(&d, r, m);
        gsl_histogram_increment(hout, d);
    }
    apop_histogram_normalize(out);
    return out;
}


/** Test the goodness-of-fit between a data vector and a model.

This just produces a PDF and calls \ref apop_pdf_test_goodness_of_fit.
apop_data *apop_pdf_test_goodness_of_fit(gsl_vector *v, apop_model *m, gsl_vector *params, int bins, int draws){
gsl_histogram_pdf   *h      = apop_vector_to_pdf(v, bins);
apop_data           *out    = apop_pdf_test_goodness_of_fit(h, m, params, bins);
    gsl_histogram_pdf_free(h);
    return out;
}
*/

static apop_data *gof_output(double diff, int bins){
apop_data   *out    = apop_data_alloc(0,4,-1);
double      pval    = gsl_cdf_chisq_P(diff, bins-1);
    apop_data_add_named_elmt(out, "Chi squared statistic", diff);
    apop_data_add_named_elmt(out, "df", bins-1);
    apop_data_add_named_elmt(out, "p value", pval);
    apop_data_add_named_elmt(out, "confidence", 1 - pval);
    return out;
}

/** Test the goodness-of-fit between two histograms (in \c apop\_model form). I assume that the histograms are aligned.

  \todo It'd be nice if this could test histograms where one has the infinibins and the other doesn't.
  \ingroup histograms
*/
apop_data *apop_histograms_test_goodness_of_fit(apop_model *m0, apop_model *m1){
  gsl_histogram *h0 = ((apop_histogram_params *)m0->model_settings)->pdf;
  gsl_histogram *h1 = ((apop_histogram_params *)m1->model_settings)->pdf;
  int     i;
  double  diff    = 0,
          bins    = h0->n;
    if (h0->n == h1->n){
        for (i=0; i< bins; i++)
            if (h0->bin[i]==0){
                if(apop_opts.verbose)
                    printf ("element %i of the first vector is zero. Skipping it.\n", i);
            } else
                diff    += gsl_pow_2(h0->bin[i] - h1->bin[i])/h0->bin[i];
    } else
        apop_error(0, 's', "%s: Sorry, I haven't implemented the case where the bin counts of the two histograms are unequal.", __func__);
    return gof_output(diff, bins);
}

/*psmirnov2x is cut/pasted/trivially modified from the R project. Copyright them. */
static double psmirnov2x(double x, int m, int n) {
    double md, nd, q, *u, w, out;
    int i, j;

    if(m > n) {
        i = n; n = m; m = i;
    }
    md = (double) (m);
    nd = (double) (n);
    q = floor(x * md * nd - 1e-7) / (md * nd);
    u = (double *) malloc((n + 1)* sizeof(double));

    for(j = 0; j <= n; j++) {
        u[j] = ((j / nd) > q) ? 0 : 1;
    }
    for(i = 1; i <= m; i++) {
        w = (double)(i) / ((double)(i + n));
        if((i / md) > q)
            u[0] = 0;
        else
            u[0] = w * u[0];
        for(j = 1; j <= n; j++) {
            if(fabs(i / md - j / nd) > q)
            u[j] = 0;
            else
            u[j] = w * u[j] + u[j - 1];
        }
    }
    out = u[n];
    free(u);
    return out;
}



/** Run the Kolmogorov test to determine whether two distributions are
 identical.

 \param m1, m2  Two matching \ref apop_histograms, probably produced via \ref apop_vectors_to_histogram or \ref apop_model_to_histogram.

 \return The \f$p\f$-value from the Kolmogorov test that the two distributions are equal.

 \ingroup histograms
 */
//apop_data *apop_test_kolmogorov(gsl_histogram *h1, gsl_histogram *h2){
apop_data *apop_test_kolmogorov(apop_model *m1, apop_model *m2){
  gsl_histogram *h1   = ((apop_histogram_params *)m1)->pdf;
  gsl_histogram *h2   = ((apop_histogram_params *)m2)->pdf;
  double    cdf1      = 0,
            cdf2      = 0,
            sum1      = 0,
            sum2      = 0,
            diff      = 0;
  int       i, offset = 0;
  gsl_histogram *first=NULL, *second=NULL;
  //If empirical data, h1->n == h2->n; if one is a distribution
  //then there are -inf and +inf bins; else, bail out.
    if(h1->n == h2->n){
        offset  = 0;
        first   = h1;
        second  = h2;
    } else if (h1->n == h2->n+2){
        offset  = 1;
        first   = h1;
        second  = h2;
        cdf1    = first->bin[0];
    } else if (h2->n == h1->n+2){
        offset  = 1;
        first   = h2;
        second  = h1;
        cdf1    = first->bin[0];
    } else 
        apop_error(0, 's', "%s: needs matching histograms.  Produce them via apop_vectors_to_histogram or apop_model_to_histogram. Returning NULL.\n", __func__);
    //Scaling step. 
    for (i=0; i< first->n; i++)
        sum1    += first->bin[i];
    for (i=0; i< second->n; i++)
        sum2    += second->bin[i];
printf("sum1: %g; sum2: %g\n", sum1, sum2);
    //Find point of greatest difference
    for (i=0; i< h2->n; i++){
        cdf1    += first->bin[i+offset]/sum1;
        cdf2    += second->bin[i]/sum2;
        diff     = GSL_MAX(diff, fabs(cdf1-cdf2));
    }
    //return 1-psmirnov2x(diff-0.02, n1, n2);
   // return 1-psmirnov2x(0.2204, n1, n2);

    apop_data   *out    = apop_data_alloc(0,3,-1);
    apop_data_add_named_elmt(out, "max distance", diff);
    apop_data_add_named_elmt(out, "p value, 2 tail", 1-psmirnov2x(diff, sum1, sum2));
    apop_data_add_named_elmt(out, "confidence, 2 tail", psmirnov2x(diff, sum1, sum2));
    return out;
}

/** Scale a histogram so it integrates to one (and is thus a proper PMF). */
void apop_histogram_normalize(apop_model *m){
  gsl_histogram *h = ((apop_histogram_params *)m->model_settings)->pdf;
  int           i;
  long double   sum = 0;
    if (!h)
        apop_error(0, 's', "%s: You sent me a model which is not a histogram or which is unparametrized.\n", __func__);
    for (i=0; i< h->n; i++)
        sum       += h->bin[i];
    if (!sum){
        apop_error(0, 'c', "%s: You sent me a histogram with a total density of zero. Returning same.\n", __func__);
        return;
    }
    for (i=0; i< h->n; i++)
        h->bin[i] /= sum;
}
