/** \file apop_hist.c */
/* Functions that work with the GSL histogram.

Copyright (c) 2006--2007, 2010 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
 (Except psmirnov2x, Copyright R Project, but also licensed under the GPL.)
*/

#include "apop_internal.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_vector.h>

/** \defgroup histograms The GSL's histograms and Apophenia's PMFs. */

/** Give me an existing histogram (as a filled-in copy of the \c apop_histogram model) and I'll
 create a new histogram with the same bins, but with data from the vector you provide 

\param template An \c apop_model produced using a form like \c apop_estimate(yourdata, apop_histogram).
\param indata The new data to be binned.
\ingroup histograms
*/
apop_model *apop_histogram_vector_reset(apop_model *template, gsl_vector *indata){
  Apop_assert(template && !strcmp(template->name, "Histogram"), "The first argument needs to be a model with appropriate apop_histogram settings.");
  apop_model *out = apop_model_copy(*template); 
  gsl_histogram *hout  = Apop_settings_get(out, apop_histogram, pdf);
    gsl_histogram_reset(hout);
    for (size_t i=0; i< indata->size; i++)
        gsl_histogram_increment(hout, gsl_vector_get(indata, i));
    return out;
}

/** Give me an existing histogram (i.e., an \c apop_model) and I'll
 create a new histogram with the same bins, but with data from \c draws
 random draws from the parametrized model you provide.

This function will normalize the output histogram to integrate to one.

\param base An \c apop_model produced using a form like \c apop_estimate(yourdata, apop_histogram). I.e. a histogram model to be used as a template. (No default)
\param m The model to be drawn from. Because this function works via random draws, the model needs to have a 
\c draw method. (No default)
\param draws The number of random draws to make. (arbitrary default = 100,000)
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

/** Make random draws from an \ref apop_model, and bin them using a binspec in the style
 of \ref apop_data_to_bins. If you have a data set that used the same binspec, you now have synced histograms, which you can plot or sensibly test hypotheses about.

The output is normalized to integrate to one.

\param binspec A description of the bins in which to place the draws; see \ref apop_data_to_bins. (default: as in \ref apop_data_to_bins.)
\param model The model to be drawn from. Because this function works via random draws, the model needs to have a 
\c draw method. (No default)
\param draws The number of random draws to make. (arbitrary default = 10,000)
\param bin_count If no bin spec, the number of bins to use (default: as per \ref apop_data_to_bins, \f$\sqrt(N)\f$)
\param rng The \c gsl_rng used to make random draws. (default: see note on \ref autorng)

\return An \ref apop_pmf model.

\li This function uses the \ref designated syntax for inputs.

\ingroup histograms
*/
APOP_VAR_HEAD apop_model *apop_model_to_pmf(apop_model *model, apop_data *binspec, long int draws, int bin_count, gsl_rng *rng){
    apop_model* apop_varad_var(model, NULL);
    Apop_assert(model && model->draw, "The second argument needs to be an apop_model with a 'draw' function "
                              "that I can use to make random draws.");
    apop_data* apop_varad_var(binspec, NULL);
    int apop_varad_var(bin_count, 0);
    long int apop_varad_var(draws, 1e4);
    gsl_rng *apop_varad_var(rng, NULL)
    static gsl_rng *spare = NULL;
    if (!rng && !spare) 
        spare = apop_rng_alloc(++apop_opts.rng_seed);
    if (!rng) rng = spare;
APOP_VAR_ENDHEAD
    Get_vmsizes(binspec);
    apop_data *outd = apop_data_alloc(draws, model->dsize); 
    for (long int i=0; i< draws; i++){
        Apop_row(outd, i, ach);
        apop_draw(ach->data, rng, model);
    }
    apop_data_to_bins(outd, binspec, .bin_count=bin_count);
    apop_vector_normalize(outd->weights);
    return apop_estimate(outd, apop_pmf);
} 

/** Test the goodness-of-fit between either two \ref apop_pmf models or two \ref apop_histogram models. 

I assume that the histograms are synced: for PMFs, you've used \ref apop_data_to_bins to generate two histograms using the same binspec, or you've used \ref apop_data_pmf_compress to guarantee that each observation value appears exactly once in each data set.  For histograms, you've use \ref apop_histogram_vector_reset or \ref apop_histogram_model_reset to ensure histograms in sync.
  
In any case, you are confident that all values in the \c observed set appear in the \c expected set with nonzero weight; otherwise this will return a \f$\chi^2\f$ statistic of \c GSL_POSINF, indicating that it is impossible for the \c observed data to have been drawn from the \c expected distribution.

\li If an observation row has weight zero, I skip it. if <tt>apop_opts.verbose >=1 </tt> I will show a warning.

  \ingroup histograms
*/
apop_data *apop_histograms_test_goodness_of_fit(apop_model *observed, apop_model *expected){
    int df; 
    double diff;
    if (Apop_settings_get_group(observed, apop_histogram)){
        gsl_histogram *h0 = Apop_settings_get(observed, apop_histogram, pdf);
        Apop_assert(Apop_settings_get_group(expected, apop_histogram), "If the first model has histogram settings (which it does), "
                                            "the second also has to have them (which it does not).");
        gsl_histogram *h1 = Apop_settings_get(expected, apop_histogram, pdf);
        Apop_assert(h0->n == h1->n, "Sorry, I haven't implemented the case where the bin counts of the two histograms are unequal.");

      df       = h0->n;
      diff     = 0;
      int bins = h0->n;
        for (int i=0; i< bins; i++)
            if (h0->bin[i]==0){
                Apop_notify(1, "element %i of the first vector is zero. Skipping it.", i);
                df --;
            } else 
                diff    += gsl_pow_2(h0->bin[i] - h1->bin[i])/h1->bin[i];
    } else { //PMF version
        df      = observed->data->weights->size;
        diff    = 0;
        for (int i=0; i< observed->data->weights->size; i++){
            Apop_data_row(observed->data, i, one_obs);
            double obs_val = gsl_vector_get(observed->data->weights, i);
            double exp_val = apop_p(one_obs, expected);
            if (exp_val == 0){
                diff = GSL_POSINF; 
                break;
            }
            if (obs_val==0){
                Apop_notify(1, "element %i of the observed data has weight zero. Skipping it.", i);
                df --;
            } else 
                diff += gsl_pow_2(obs_val - exp_val)/exp_val;
        }
    }
    //Data gathered. Now output
    apop_data   *out    = apop_data_alloc();
    double      toptail = gsl_cdf_chisq_Q(diff, df-1);
    sprintf(out->names->title, "Goodness-of-fit test via Chi-squared statistic");
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
        w = (double)(i) / (i + nd);
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

/** Run the Kolmogorov-Smirnov test to determine whether two distributions are identical.

 \param m1, m2  Two models, most likely of \ref apop_pmf type. I will ue the cdf method, so if your function doesn't have one, expect this to run the slow default. I run it for each row of each data set, so if your model has a \c NULL at the data, I won't know what to do. 
 
I will also give special handling to two synced \ref apop_histogram "apop_histograms" (probably produced via \ref apop_histogram_vector_reset or \ref apop_histogram_model_reset), using the histogram bins to define how the CDF is built.

 \return An \ref apop_data set including the \f$p\f$-value from the Kolmogorov test that the two distributions are equal.

 \code
 apop_data_sort(data1);
 apop_data_sort(data2);
 apop_model *m1 = apop_estimate(data1, apop_pmf);
 apop_model *m2 = apop_estimate(data2, apop_pmf);
 apop_data *ktest = apop_test_kolmogorov(m1, m2);
 //The sort of output you could pull out:
 double k_statistic = apop_data_get(ktest, .rowname="max distance");
 double pval = apop_data_get(ktest, .rowname="p value, 2 tail");
 double confidence_of_inequality = apop_data_get(ktest, .rowname="confidence, 2 tail");
 \endcode

  \li I assume that the data sets are sorted.

 \ingroup histograms
 */
apop_data *apop_test_kolmogorov(apop_model *m1, apop_model *m2){
    double largest_diff;
    int maxsize1, maxsize2;
    if (Apop_settings_get_group(m1, apop_histogram)){
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
            Apop_assert(0, "I need matching histograms. Produce them via apop_histogram_vector_reset "
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
        largest_diff = diff;
        maxsize1 = first->n;
        maxsize2 = second->n;
    } else { //no pair of histograms
        Apop_assert(m1->data, "I will test the CDF at each point in the data set, but the first model has a NULL data set. "
                              "Maybe generate, then apop_data_sort, a few thousand random draws?");
        Apop_assert(m2->data, "I will test the CDF at each point in the data set, but the second model has a NULL data set. "
                              "Maybe generate, then apop_data_sort, a few thousand random draws?");
        {Get_vmsizes(m1->data); maxsize1 = maxsize;}//copy one of the macro's variables 
        {Get_vmsizes(m2->data); maxsize2 = maxsize;}//  to the full function's scope.
        largest_diff=GSL_NEGINF;
        for (size_t i=0; i< maxsize1; i++){
            Apop_data_row(m1->data, i, arow);
            largest_diff     = GSL_MAX(largest_diff, fabs(apop_cdf(arow, m1)-apop_cdf(arow, m2)));
        }
        for (size_t i=0; i< maxsize2; i++){     //There should be matched data rows, so there is redundancy.
            Apop_data_row(m2->data, i, arow); // Feel free to submit a smarter version.
            largest_diff     = GSL_MAX(largest_diff, fabs(apop_cdf(arow, m1)-apop_cdf(arow, m2)));
        }
    }
    apop_data   *out    = apop_data_alloc();
    sprintf(out->names->title, "Kolmogorov-Smirnov test");
    apop_data_add_named_elmt(out, "max distance", largest_diff);
    apop_data_add_named_elmt(out, "p value, 2 tail", 1-psmirnov2x(largest_diff, maxsize1, maxsize2));
    apop_data_add_named_elmt(out, "confidence, 2 tail", psmirnov2x(largest_diff, maxsize1, maxsize2));
    return out;
}

/** Scale a histogram so it integrates to one (and is thus a proper PMF). */
void apop_histogram_normalize(apop_model *m){
    gsl_histogram *h = Apop_settings_get(m, apop_histogram, pdf);
    int           i;
    long double   sum = 0;
    Apop_assert_n(h, "You sent me a model which is not a histogram or which is unparametrized.");
    for (i=0; i< h->n; i++)
        sum       += h->bin[i];
    Apop_assert_c(sum, ,0, "You sent me a histogram with a total density of zero. Returning same.");
    for (i=0; i< h->n; i++)
        h->bin[i] /= sum;
}


/** Create a histogram from data by putting data into bins of fixed width. 

\param indata The input data that will be binned. This is modified in place, so make a copy if you want to retain the original data.
\param close_top_bin Normally, a bin covers the range from the point equal to its minimum to points strictly less than
the minimum plus the width.  if \c 'y', then the top bin includes points less than or equal to the upper bound. This solves the problem of displaying histograms where the top bin is just one point.
\param binspec This is an \ref apop_data set with the same number of columns as \c indata. 
If you want a fixed size for the bins, then the first row of the bin spec is the bin width for each column.
This allows you to specify a width for each dimension, or specify the same size for all with something like:

\param bin_count If you don't provide a bin spec, I'll provide this many evenly-sized bins. Default: \f$\sqrt(N)\f$.  \code
Apop_data_row(indata, 0, firstrow);
apop_data *binspec = apop_data_copy(firstrow);
gsl_matrix_set_all(binspec->matrix, 10); //bins of size 10 for all dim.s
apop_data_to_bins(indata, binspec);
\endcode
The presumption is that the first bin starts at zero in all cases. You can add a second row to the spec to give the offset for each dimension.  Default: NULL. if no binspec and no binlist, then a grid with offset equal to the min of the column, and bin size such that it takes \f$\sqrt{N}\f$ bins to cover the range to the max element. 


\return A pointer to \c indata, now properly binned.  If you didn't give me a binspec, then I attach one to your data set as a page named \c \<binspec\>, so you can snap a second data set to the same grid using 
\code
apop_data_to_bins(first_set, NULL);
apop_data_to_bins(second_set, apop_data_get_page(first_set, "<binspec>"));
\endcode


  The text segment, if any, is not binned. I use \ref apop_data_pmf_compress as the final step in the binning, 
  and that does respect the text segment. So given the standard one-unit grid, the data:

<table frame=box>
<tr>
<td>Vector value</td><td> Text name</td><td>Weights</td>
</tr><tr valign=bottom>
<td align=center>
</td></tr>
<tr><td>1.1</td><td>Type 1</td><td>1</td></tr>
<tr><td>2.1</td><td>Type 1</td><td>1</td></tr>
<tr><td>2</td><td>Type 1</td><td>1</td></tr>
<tr><td>1</td><td>Type 1</td><td>1</td></tr>
<tr><td>1</td><td>Type 2</td><td>1</td></tr>
</table>

will bin into a histogram like:

<table frame=box>
<tr>
<td>Vector value</td><td> Text name</td><td>Weights</td>
</tr><tr valign=bottom>
<td align=center>
</td></tr>
<tr><td>1</td><td>Type 1</td><td>1</td></tr>
<tr><td>2</td><td>Type 1</td><td>2</td></tr>
<tr><td>1</td><td>Type 1</td><td>1</td></tr>
<tr><td>1</td><td>Type 2</td><td>1</td></tr>
</table>
*/
APOP_VAR_HEAD apop_data *apop_data_to_bins(apop_data *indata, apop_data *binspec, int bin_count, char close_top_bin){
    apop_data *apop_varad_var(indata, NULL);
    Apop_assert_c(indata, NULL, 1, "NULL input data set, so returning NULL output data set.");
    apop_data *apop_varad_var(binspec, NULL);
    char apop_varad_var(close_top_bin, 'n');
    int apop_varad_var(bin_count, 0);
APOP_VAR_ENDHEAD
    Get_vmsizes(indata); //firstcol, vsize, msize1, msize2
    double binwidth, offset, max=0;
    apop_data *bs = binspec ? binspec
                    : apop_data_add_page(indata, 
                        apop_data_alloc(vsize? 2: 0, msize1? 2: 0, indata->matrix ? msize2: 0),
                        "<binspec>");
    for (int j= firstcol; j< msize2; j++){
        Apop_col(indata, j, onecol);
        if (binspec){
           binwidth = apop_data_get(binspec, 0, j);
           offset = ((binspec->vector && binspec->vector->size==2 )
                   ||(binspec->matrix && binspec->matrix->size1==2)) ? apop_data_get(binspec, 1, j) : 0;
        } else {
            Apop_col(bs, j, abin);
            max = gsl_vector_max(onecol);
            offset = abin->data[1] = gsl_vector_min(onecol);
            binwidth = abin->data[0] = (max - offset)/(bin_count ? bin_count : sqrt(onecol->size));
        }
        for (int i=0; i< onecol->size; i++){
            double val = gsl_vector_get(onecol, i);
            if (close_top_bin=='y' && val == max && val!=offset) 
                val -= 2*GSL_DBL_EPSILON;
            gsl_vector_set(onecol, i, (floor((val -offset)/binwidth))*binwidth+offset);
        }
    }
    apop_data_pmf_compress(indata);
    return indata;
}

/** For usage, see the documentation for the \ref apop_pmf model. 
  \param d An input crosstab
  \return A PMF model which has a single line for each nonzero line of the crosstab.

\li  The \ref apop_pmf really keeps its information in the \c data element, where the list of observations and their weights are held. When an \ref apop_model is freed, the \c data element is untouched, but the \c parameters element is freed. This function generates a new data set that is basically internal to the model and should be freed with the model, and this is achieved by pointing both the \c data and \c parameters elements to the same data set. You probably won't ever touch the \c data set inside the model returned by this function, but bear in mind that, unlike the typical case, it will disappear after you free the model.
 */
apop_model *apop_crosstab_to_pmf (apop_data *d){
    Get_vmsizes(d) //tsize
    int use_matrix=0, use_vector = 0;
    size_t ctr = 0;
    double x;
    if (d->matrix) use_matrix++;
    else if (d->vector) use_vector++;
    else Apop_assert(0, "You gave me an input set with neither vector nor matrix data.");
    apop_model *out = apop_model_copy(apop_pmf);
    out->parameters = apop_data_alloc(0, tsize, (use_matrix ? 2: 1));
    out->parameters->weights = gsl_vector_alloc(tsize);
    out->data = out->parameters;
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

