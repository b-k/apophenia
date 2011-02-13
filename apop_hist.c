/** \file apop_hist.c */
/* Functions that work with the GSL histogram.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
 (Except psmirnov2x, Copyright R Project, but also licensed under the GPL.)
*/

#include "db.h"     //just for apop_opts
#include "asst.h" //rng_alloc
#include "stats.h"
#include "model.h"
#include "settings.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_vector.h>

/** Give me an existing histogram (i.e., an \c apop_model) and I'll
 create a new histogram with the same bins, but with data from the vector you provide 

\param template An \c apop_model produced using a form like \c apop_estimate(yourdata, apop_histogram).
\param indata The new data to be binned.
\ingroup histograms
*/
apop_model *apop_histogram_vector_reset(apop_model *template, gsl_vector *indata){
  apop_assert_s(template && !strcmp(template->name, "Histogram"), "The first argument needs to be a model with appropriate apop_histogram settings.");
  size_t  i;
  apop_model *out = apop_model_copy(*template); 
  gsl_histogram *hout  = Apop_settings_get(out, apop_histogram, pdf);
    gsl_histogram_reset(hout);
    for (i=0; i< indata->size; i++)
        gsl_histogram_increment(hout, gsl_vector_get(indata, i));
    return out;
}

/** Give me an existing histogram (i.e., an \c apop_model) and I'll
 create a new histogram with the same bins, but with data from \c draws
 random draws from the parametrized model you provide.

Unlike with most other histogram-generating functions, this one will normalize the output to integrate to one.

\param base An \c apop_model produced using a form like \c apop_estimate(yourdata, apop_histogram). I.e. a histogram model to be used as a template. (No default)
\param m The model to be drawn from. Because this function works via random draws, the model needs to have a 
\c draw method. (No default)
\param draws The number of random draws to make. (arbitrary default = 1e5)
\param rng The \c gsl_rng used to make random draws. (default: see note on \ref autorng)

\li This function uses the \ref designated syntax for inputs.

\ingroup histograms
*/
APOP_VAR_HEAD apop_model *apop_histogram_model_reset(apop_model *base, apop_model *m, long int draws, gsl_rng *rng){
    static gsl_rng *spare = NULL;
    apop_model* apop_varad_var(base, NULL);
    apop_model* apop_varad_var(m, NULL);
    long int apop_varad_var(draws, 1e5);
  Apop_assert(base && !strcmp(base->name, "Histogram"), "The first argument needs to be a model with appropriate apop_histogram settings.");
  Apop_assert(m && m->draw, "The second argument needs to be an apop_model with a 'draw' function that I can use to make random draws.");
    gsl_rng *apop_varad_var(rng, NULL)
    if (!rng && !spare) 
        spare = apop_rng_alloc(++apop_opts.rng_seed);
    if (!rng)  rng = spare;
APOP_VAR_ENDHEAD
    double d;
    apop_model *out = apop_model_copy(*base); 
    gsl_histogram *hout  = Apop_settings_get(out, apop_histogram, pdf);
    gsl_histogram_reset(hout);
    for (long int i=0; i< draws; i++){
        m->draw(&d, rng, m);
        gsl_histogram_increment(hout, d);
    }
    apop_histogram_normalize(out);
    return out;
}

/** Test the goodness-of-fit between two histograms (in \c apop_model form). I assume that the histograms are aligned.
  \ingroup histograms
*/
apop_data *apop_histograms_test_goodness_of_fit(apop_model *m0, apop_model *m1){
  gsl_histogram *h0 = Apop_settings_get(m0, apop_histogram, pdf);
  gsl_histogram *h1 = Apop_settings_get(m1, apop_histogram, pdf);
    Apop_assert(h0, "The first model you gave me has a NULL PDF.");
    Apop_assert(h1, "The second model you gave me has a NULL PDF.");
    Apop_assert(h0->n == h1->n, "Sorry, I haven't implemented the case where the bin counts of the two histograms are unequal.");

  int     df      = h0->n,
          bins    = h0->n;
  double  diff    = 0;
    for (int i=0; i< bins; i++)
        if (h0->bin[i]==0){
            Apop_notify(1, "element %i of the first vector is zero. Skipping it.", i);
            df --;
        } else 
            diff    += gsl_pow_2(h0->bin[i] - h1->bin[i])/h0->bin[i];
    //Data gathered. Now output
    apop_data   *out    = apop_data_alloc(0,4,-1);
    double      toptail = gsl_cdf_chisq_Q(diff, df-1);
    apop_data_add_named_elmt(out, "Chi squared statistic", diff);
    apop_data_add_named_elmt(out, "df", df-1);
    apop_data_add_named_elmt(out, "p value",  toptail); 
    apop_data_add_named_elmt(out, "confidence", 1 - toptail);
    return out;
}

/*psmirnov2x is cut/pasted/trivially modified from the R project. Copyright them. */
static double psmirnov2x(double x, int m, int n) {
    double md, nd, q, *u, w, out;
    int tmp, j;

    if(m > n) {
        tmp = n; n = m; m = tmp;
    }
    md = (double) (m);
    nd = (double) (n);
    q = floor(x * md * nd - 1e-7) / (md * nd);
    u = (double *) malloc((n + 1)* sizeof(double));

    for(j = 0; j <= n; j++) 
        u[j] = ((j / nd) > q) ? 0 : 1;
    for(int i = 1; i <= m; i++) {
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

 \param m1, m2  Two matching \ref apop_histogram "apop_histograms", probably produced via \ref apop_histogram_vector_reset or \ref apop_histogram_model_reset.

 \return An \ref apop_data set including the \f$p\f$-value from the Kolmogorov test that the two distributions are equal.

 \ingroup histograms
 */
apop_data *apop_test_kolmogorov(apop_model *m1, apop_model *m2){
  gsl_histogram *h1 = Apop_settings_get(m1, apop_histogram, pdf);
  gsl_histogram *h2 = Apop_settings_get(m2, apop_histogram, pdf);
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
        Apop_assert_c(0, NULL, 0, "I need matching histograms.  Produce them via apop_histogram_vector_reset "
                       "or apop_histogram_model_reset. Returning NULL.");
    //Scaling step. 
    for (i=0; i< first->n; i++)
        sum1    += first->bin[i];
    for (i=0; i< second->n; i++)
        sum2    += second->bin[i];
    if (apop_opts.verbose)
        printf("sum1: %g; sum2: %g\n", sum1, sum2);
    //Find point of greatest difference
    for (i=0; i< h2->n; i++){
        cdf1    += first->bin[i+offset]/sum1;
        cdf2    += second->bin[i]/sum2;
        diff     = GSL_MAX(diff, fabs(cdf1-cdf2));
    }

    apop_data   *out    = apop_data_alloc();
    apop_data_add_named_elmt(out, "max distance", diff);
    apop_data_add_named_elmt(out, "p value, 2 tail", 1-psmirnov2x(diff, first->n, second->n));
    apop_data_add_named_elmt(out, "confidence, 2 tail", psmirnov2x(diff, first->n, second->n));
    return out;
}

/** Scale a histogram so it integrates to one (and is thus a proper PMF). */
void apop_histogram_normalize(apop_model *m){
  gsl_histogram *h = Apop_settings_get(m, apop_histogram, pdf);
  int           i;
  long double   sum = 0;
  Apop_assert(h, "You sent me a model which is not a histogram or which is unparametrized.");
    for (i=0; i< h->n; i++)
        sum       += h->bin[i];
    Apop_assert_c(sum, ,0, "You sent me a histogram with a total density of zero. Returning same.");
    for (i=0; i< h->n; i++)
        h->bin[i] /= sum;
}
