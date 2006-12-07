/** \file apop_histogram.c	PMF and CMF manipulations.

 Copyright 2006 by Ben Klemens. Licensed under the GNU GPL v2.
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
        vector          -> cmf
        pmf             -> cmf
        (pmf, pmf)      -> gof


*/

#include "db.h"     //just for apop_opts
#include "apophenia/stats.h"
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
the GSL calls this the <tt>gsl_histogram_PMF</tt> structure. That's right:
the the data in the <tt>gsl_histogram_PMF</tt> structure is a cumulative
sum---a CMF.

Anyway, here are some functions to deal with these various histograms and such.

 */

/** Produce a <tt>gsl_histogram</tt> from a vector. 

  \param data   The data vector
  \param bins   The number of (evenly-spaced) bins
  \ingroup histograms
  */
gsl_histogram * apop_vector_to_histogram(gsl_vector *data, int bins){
int                 i;
gsl_histogram       *h  = gsl_histogram_alloc(bins);
double              min, max;
    gsl_vector_minmax(data, &min, &max);
    gsl_histogram_set_ranges_uniform(h, min, max+10*GSL_DBL_EPSILON);
    for (i=0; i< data->size; i++)
        gsl_histogram_increment(h,gsl_vector_get(data,i));
    return h;
}

/** produce a GSL_histogram_pdf structure.

The GSL provides a means of making random draws from a data set, or to
put it another way, to produce an artificial PDF from a data set. It
does so by taking a histogram and producing a CDF. 

This function takes the requisite steps for you, producing a histogram
from the data and then converting it to a <tt>gsl_histogram_pdf</tt>
structure from which draws can be made. Usage:

\code
    //assume data is a gsl_vector* already filled with data.
    gsl_histogram_pdf *p = apop_vector_to_cmf(data, 1000);
    gsl_rng_env_setup();
    gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);

    //Draw from the PDF:
    gsl_histogram_pdf_sample(p, gsl_rng_uniform(r));

    //Eventually, clean up:
    gsl_histogram_pdf_free(p);
\endcode

\param data a <tt>gsl_vector*</tt> with the sample data
\param bins The number of bins in the artificial PDF. It is OK if most
of the bins are empty, so feel free to set this to a thousand or even
a million, depending on the level of resoultion your data has.

\ingroup histograms
*/
gsl_histogram_pdf * apop_vector_to_cmf(gsl_vector *data, int bins){
gsl_histogram_pdf   *p  = gsl_histogram_pdf_alloc(bins);
gsl_histogram       *h  = apop_vector_to_histogram(data,bins);
    gsl_histogram_pdf_init(p, h);
    gsl_histogram_free(h);
    return p;
}

/** Make two histograms that share bin structures: same max/min, same
 bin partitions. You can then plot the two together or run goodness-of-fit tests.

\param v1       a gsl_vector
\param v2       another gsl_vector
\param bins     How many bins?
\return an array of two histograms. 

\code
gsl_histogram *out = apop_vectors_to_histograms(d1, d2, 1000);
apop_test_goodness_of_fit(out[0], out[1]);
\endcode

\ingroup histograms
*/
gsl_histogram ** apop_vectors_to_histograms(gsl_vector *v1, gsl_vector *v2, int bins){
int                 i;
gsl_histogram       **h  = malloc(sizeof(gsl_histogram)*2);
double              max0, max1, max2, min0, min1,min2;
    h[0]    = gsl_histogram_alloc(bins);
    h[1]    = gsl_histogram_alloc(bins);
    gsl_vector_minmax(v1, &min1, &max1);
    gsl_vector_minmax(v2, &min2, &max2);
    min0    = GSL_MIN(min1, min2);
    max0    = GSL_MAX(max1, max2)+1e-5;
    gsl_histogram_set_ranges_uniform(h[0], min0, max0);
    gsl_histogram_set_ranges_uniform(h[1], min0, max0);
    for (i=0; i< v1->size; i++)
        gsl_histogram_increment(h[0], gsl_vector_get(v1,i));
    for (i=0; i< v2->size; i++)
        gsl_histogram_increment(h[1], gsl_vector_get(v2,i));
    return h;
}

/** Produce a histogram for a model matching an existing histogram of data.

The method is to produce a histogram for the PDF using the RNG.

\todo The double* that gets sent in to the model RNGs is ungraceful.
\ingroup histograms
*/
gsl_histogram * apop_model_to_histogram(apop_model m, gsl_histogram *h, 
                        int draws, gsl_vector *params, gsl_rng *r){
int     i;
//int     bc      = h ? h->n + 2 : bins;
int     bc      =  h->n + 2;
    //Must be careful that the histograms match, but the range for the
    //PDF is \f$(-\infty, infty)\f$. 
            //GSL documentation says newbins should be one elemnt too big.
    double  *newbins = malloc(sizeof(double)* (bc +1));
    newbins[0]                  = GSL_NEGINF;
    memcpy((newbins + 1), h->range, sizeof(double) * h->n);
    newbins[bc]                 = GSL_POSINF;
    gsl_histogram *modelhist    = gsl_histogram_alloc(bc);
    gsl_histogram_set_ranges(modelhist, newbins, bc+1);
    for (i=0; i< draws; i++)
        gsl_histogram_increment(modelhist, m.draw(r,params, NULL));
    for (i=0; i< modelhist->n; i++)
        modelhist->bin[i] /= draws;
    return modelhist;
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
apop_data   *out    = apop_data_alloc(4,-1);
double      pval    = gsl_cdf_chisq_P(diff, bins-1);
    apop_data_add_named_elmt(out, "Chi squared statistic", diff);
    apop_data_add_named_elmt(out, "df", bins-1);
    apop_data_add_named_elmt(out, "p value", pval);
    apop_data_add_named_elmt(out, "confidence", 1 - pval);
    return out;
}

/** Test the goodness-of-fit between two histograms

  \todo Right now, I'm assuming the histograms are aligned---h0->bins == h1->bins. This needs to go away.
  \ingroup histograms
*/
apop_data *apop_histograms_test_goodness_of_fit(gsl_histogram *h0, gsl_histogram *h1){
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
        printf("Sorry, I haven't implemented the case where the bin counts of the two histograms are unequal.");
    return gof_output(diff, bins);
}

/** Test the goodness-of-fit between two vectors
 Here, I assume that the vectors are already histogram-like,
 and that the first element of v1 corresponds to the first
 element of v2, et cetera. If you have unordered data, see \ref
 apop_histograms_test_goodness_of_fit.

The denominator of the chi-squared test is the first vector, unless
the second vector is longer:

\code
gsl_vector  *denom  = (v1->size > v0->size) ? v1 : v0;
\endcode

If any elements of your denominator vector are zero, then the chi-squared
statistic will be GSL_POSINF. The claim of the test is that the test
vector was drawn from the denominator vector, and if the denominator
vector has density zero at some point while the test vector does not,
then the claim must be false.

Meanwhile, if either vector's value is GSL_NAN, then that bin is
skipped. If you really think the zeros are a fluke, then try
\code
apop_vector_replace(v1, apop_double_is_zero, GSL_NAN);
apop_vector_replace(v2, apop_double_is_zero, GSL_NAN);
apop_vectors_test_goodness_of_fit(v1,v2);
\endcode


\param  v0, v1  The first vector is assumed to be the denominator of
the test (unless the second is longer; see notes).

\return     An \ref apop_data table with the chi-squared statistic, df, p value, confidence, et cetera.

\ingroup histograms
*/
apop_data *apop_vectors_test_goodness_of_fit(gsl_vector *v0, gsl_vector *v1){
double      diff    = 0,
            d0, d1;
gsl_vector  *denom  = (v1->size > v0->size) ? v1 : v0,
            *num    = (v1->size > v0->size) ? v0 : v1;
int         i, 
            nanct   = 0;
    for (i=0; i< denom->size; i++){
        d0  = gsl_vector_get(denom, i);
        d1  = (i < num->size) ? gsl_vector_get(num, i): 0;
        if (gsl_isnan(d0) || gsl_isnan(d1))
            nanct   ++;
        else {
            if (d0 == 0){
                if(apop_opts.verbose)
                    printf ("element %i of the denominator vector is zero. Returning GSL_POSINF\n", i);
                diff     = GSL_POSINF;
                break;
            } else
                diff    += gsl_pow_2(d0 - d1)/d0;
        }
    }
    return gof_output(diff, v0->size - nanct);
}

/** Test the goodness-of-fit between a histogram and a model.

  The method is to produce a histogram for the PDF using the RNG, then do a chi-squared test on the two PDFs.

\todo The double* that gets sent in to the model RNGs is ungraceful.
\ingroup histograms
*/
apop_data *apop_model_test_goodness_of_fit(gsl_vector *v1, apop_model m,
int bins, long int draws, gsl_vector *params, gsl_rng *r){
int     i, count    = bins;
double  diff        = 0,
        sum;
gsl_histogram *h1   = apop_vector_to_histogram(v1, bins),
            *h2     = apop_model_to_histogram(m, h1, draws, params, r);
    sum      = gsl_histogram_sum(h1);//h1 is not normalized; h2 is.
    diff    += gsl_pow_2(h2->bin[0]);
    for (i=0; i< bins; i++){
        if (h1->bin[i])
            diff    += gsl_pow_2(h1->bin[i]/sum - h2->bin[i+1])/(h1->bin[i]/sum);
        else
            count   --;//ignore this bin.
        }
    diff    += gsl_pow_2(h2->bin[bins+2]);
    gsl_histogram_free(h1);
    gsl_histogram_free(h2);
    return gof_output(diff, count);
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

 \param h1, h2  Two matching histograms, probably produced via \ref apop_vectors_to_histogram or \ref apop_model_to_histogram.

 \param n1, n2 How many data points are in your set? Since this is
 not inferred from the histograms, you can set them to an arbitrarily
 large number of bins.

 \return The \f$p\f$-value from the Kolmogorov test that the two distributions are equal.

 \ingroup histograms
 */
apop_data *apop_test_kolmogorov(gsl_histogram *h1, gsl_histogram *h2){
  double    cdf1    = 0,
            cdf2    = 0,
            sum1    = 0,
            sum2    = 0,
            diff    = 0;
  int       offset, i;
  gsl_histogram *first, *second;
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
    } else {
        printf("apop_test_kolmogorov needs matching histograms.  Produce them via apop_vectors_to_histogram or apop_model_to_histogram. Returning NULL.\n");
        return NULL;
    }
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

apop_data   *out    = apop_data_alloc(3,-1);
    apop_data_add_named_elmt(out, "max distance", diff);
    apop_data_add_named_elmt(out, "p value, 2 tail", 1-psmirnov2x(diff, sum1, sum2));
    apop_data_add_named_elmt(out, "confidence, 2 tail", psmirnov2x(diff, sum1, sum2));
    return out;
}
