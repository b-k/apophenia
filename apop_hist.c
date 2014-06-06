
/** \file apop_hist.c */
/* Functions that work with PMFs and histograms.

Copyright (c) 2006--2007, 2010, 2013 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
 (Except psmirnov2x, Copyright R Project, but also licensed under the GPL.)
*/

#include "apop_internal.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_vector.h>
#include <stdbool.h>

/** \defgroup histograms The GSL's histograms and Apophenia's PMFs. */

/** Make random draws from an \ref apop_model, and bin them using a binspec in the style
 of \ref apop_data_to_bins. If you have a data set that used the same binspec, you now have synced histograms, which you can plot or sensibly test hypotheses about.

The output is normalized to integrate to one.

\param binspec A description of the bins in which to place the draws; see \ref apop_data_to_bins. (default: as in \ref apop_data_to_bins.)
\param model The model to be drawn from. Because this function works via random draws, the model needs to have a 
\c draw method. (No default)
\param draws The number of random draws to make. (arbitrary default = 10,000)
\param bin_count If no bin spec, the number of bins to use (default: as per \ref apop_data_to_bins, \f$\sqrt(N)\f$)
\param rng The \c gsl_rng used to make random draws. (default: an RNG from \ref apop_rng_get_thread)

\return An \ref apop_pmf model.

\li This function uses the \ref designated syntax for inputs.

\ingroup histograms
*/
#ifdef APOP_NO_VARIADIC
apop_model * apop_model_to_pmf(apop_model *model, apop_data *binspec, long int draws, int bin_count, gsl_rng *rng){
#else
apop_varad_head(apop_model *, apop_model_to_pmf){
    apop_model* apop_varad_var(model, NULL);
    Apop_assert(model && model->draw, "The second argument needs to be an apop_model with a 'draw' function "
                              "that I can use to make random draws.");
    apop_data* apop_varad_var(binspec, NULL);
    int apop_varad_var(bin_count, 0);
    long int apop_varad_var(draws, 1e4);
    gsl_rng *apop_varad_var(rng, apop_rng_get_thread())
    return apop_model_to_pmf_base(model, binspec, draws, bin_count, rng);
}

 apop_model * apop_model_to_pmf_base(apop_model *model, apop_data *binspec, long int draws, int bin_count, gsl_rng *rng){
#endif
    Get_vmsizes(binspec);
    apop_data *outd = apop_data_alloc(draws, model->dsize); 
    for (long int i=0; i< draws; i++){
        Apop_row_v(outd, i, ach);
        apop_draw(ach->data, rng, model);
    }
    apop_data *outbinned = apop_data_to_bins(outd, binspec, .bin_count=bin_count);
    apop_data_free(outd);
    apop_vector_normalize(outbinned->weights);
    return apop_estimate(outbinned, apop_pmf);
} 

/** Test the goodness-of-fit between two \ref apop_pmf models. 

If you send two histograms, I assume that the histograms are synced: for PMFs,
you've used \ref apop_data_to_bins to generate two histograms using the same binspec,
or you've used \ref apop_data_pmf_compress to guarantee that each observation value
appears exactly once in each data set.

In any case, you are confident that all values in the \c observed set appear in the \c
expected set with nonzero weight; otherwise this will return a \f$\chi^2\f$ statistic
of \c GSL_POSINF, indicating that it is impossible for the \c observed data to have
been drawn from the \c expected distribution.

\li If an observation row has weight zero, I skip it. if <tt>apop_opts.verbose >=1 </tt> I will show a warning.

  \ingroup histograms
*/
apop_data *apop_histograms_test_goodness_of_fit(apop_model *observed, apop_model *expected){
    int df = observed->data->weights->size;
    double diff = 0;
    for (int i=0; i< observed->data->weights->size; i++){
        Apop_row(observed->data, i, one_obs);
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
    //Data gathered. Now output
    apop_data   *out    = apop_data_alloc();
    double      toptail = gsl_cdf_chisq_Q(diff, df-1);
    Asprintf(&out->names->title, "Goodness-of-fit test via Chi-squared statistic");
    apop_data_add_named_elmt(out, "Chi squared statistic", diff);
    apop_data_add_named_elmt(out, "df", df-1);
    apop_data_add_named_elmt(out, "p value",  toptail); 
    apop_data_add_named_elmt(out, "confidence", 1 - toptail);
    return out;
}

/*Everything from here to psmirnov2x (inclusive) is cut/pasted/trivially modified from the R project. Copyright them. */
static void m_multiply(long double *A, long double *B, long double *C, int m) {
    /* Auxiliary routine used by K().
       Matrix multiplication.
    */
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) {
            long double s = 0;
            for (int k = 0; k < m; k++)
                s += A[i * m + k] * B[k * m + j];
            C[i * m + j] = s;
        }
}

static void m_power(long double *A, int eA, long double *V, int *eV, int m, int n) {
    /* Auxiliary routine used by K().
       Matrix power.
    */
    long double *B;
    int eB, i;

    if (n == 1) {
        for (i = 0; i < m * m; i++)
            V[i] = A[i];
        *eV = eA;
        return;
    }
    m_power(A, eA, V, eV, m, n / 2);
    B = calloc(m * m, sizeof(long double));
    m_multiply(V, V, B, m);
    eB = 2 * (*eV);
    if ((n % 2) == 0) {
        for (i = 0; i < m * m; i++)
            V[i] = B[i];
        *eV = eB;
    } else {
        m_multiply(A, B, V, m);
        *eV = eA + eB;
    }
    if (V[(m / 2) * m + (m / 2)] > 1e140) {
        for (i = 0; i < m * m; i++)
            V[i] = V[i] * 1e-140;
        *eV += 140;
    }
    free(B);
}

/* The two-sided one-sample 'exact' distribution */
static double kolmogorov_2x(int n, double d) {
    /* Compute Kolmogorov's distribution.
       Code published in
	 George Marsaglia and Wai Wan Tsang and Jingbo Wang (2003),
	 "Evaluating Kolmogorov's distribution".
	 Journal of Statistical Software, Volume 8, 2003, Issue 18.
	 URL: http://www.jstatsoft.org/v08/i18/.
    */

   int k, m, i, j, g, eH, eQ;
   long double h, s, *H, *Q;

   /* The faster right-tail approximation is omitted here.
      s = d*d*n; 
      if(s > 7.24 || (s > 3.76 && n > 99)) 
          return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
   */
   k = (n * d) + 1;
   m = 2 * k - 1;
   h = k - n * d;
   H = calloc(m * m, sizeof(long double));
   Q = calloc(m * m, sizeof(long double));
   for(i = 0; i < m; i++)
       for(j = 0; j < m; j++)
           if (i - j + 1 < 0) H[i * m + j] = 0;
           else               H[i * m + j] = 1;
   for(i = 0; i < m; i++) {
       H[i * m] -= pow(h, i + 1);
       H[(m - 1) * m + i] -= pow(h, (m - i));
   }
   H[(m - 1) * m] += ((2 * h - 1 > 0) ? pow(2 * h - 1, m) : 0);
   for(i = 0; i < m; i++)
       for(j=0; j < m; j++)
           if(i - j + 1 > 0)
               for(g = 1; g <= i - j + 1; g++)
                   H[i * m + j] /= g;
   eH = 0;
   m_power(H, eH, Q, &eQ, m, n);
   s = Q[(k - 1) * m + k - 1];
   for(i = 1; i <= n; i++) {
       s = s * i / n;
       if(s < 1e-140) {
           s *= 1e140;
           eQ -= 140;
       }
   }
   s *= pow(10., eQ);
   free(H);
   free(Q);
   return(s);
}

static double psmirnov2x(double x, int m, int n) {
    if(m > n) {
        int tmp = n; n = m; m = tmp;
    }
    double md = m;
    double nd = n;
        // q has 0.5/mn added to ensure that rounding error doesn't
        // turn an equality into an inequality, eg abs(1/2-4/5)>3/10
    long double q = (.5+floor(x * md * nd - 1e-7)) / (md * nd);
    long double u[n+1];

    for(int j = 0; j <= n; j++) 
        u[j] = ((j / nd) > q) ? 0 : 1;
    for(int i = 1; i <= m; i++) {
        long double w = i/(i + nd);
        u[0] = (i / md) > q 
                ? 0
                : w * u[0];
        for(int j = 1; j <= n; j++) 
            u[j] = fabs(i / md - j / nd) > q
                    ? 0
                    : w * u[j] + u[j - 1];
    }
    return u[n];
}

/** Run the Kolmogorov-Smirnov test to determine whether two distributions are identical.

\param m1 A sorted PMF model. I.e., a model estimated via something like 

\code
apop_model *m1 = apop_estimate(apop_data_sort(input_data), apop_pmf);
\endcode

\param m2  Another \ref apop_model. If it is a PMF, then I will use a two-sample test,
which is different from the one-sample test used if this is not a PMF.

\return An \ref apop_data set including the \f$p\f$-value from the Kolmogorov-Smirnov
test that the two distributions are equal.

\exception out->error='m'  Model error: \c m1 is not an \ref apop_pmf. I verify this
by checking whether <tt>m1->cdf == apop_pmf->cdf</tt>.

\li If you are using a \ref apop_pmf model, <b>the data set(s) must be sorted before
you call this.</b> See \ref apop_data_sort and the discussion of CDFs in the \ref
apop_pmf documentation. If you don't do this, the test will almost certainly reject
the null hypothesis that \c m1 and \c m2 are identical.

Here is an example, which tests whether a set of draws from a Normal(0, 1) matches a
sequence of Normal distributions with increasing mean.

\include ks_tests.c

\ingroup histograms
*/
apop_data *apop_test_kolmogorov(apop_model *m1, apop_model *m2){
    Apop_stopif(m1->cdf != apop_pmf->cdf, apop_return_data_error('m'), 
            0, "First model has to be a PMF. I check whether m1->cdf == apop_pmf->cdf.");
    bool m2_is_pmf = (m2->cdf == apop_pmf->cdf);

    int maxsize1, maxsize2;
    {Get_vmsizes(m1->data); maxsize1 = maxsize;} //copy one of the macro's variables 
    {Get_vmsizes(m2->data); maxsize2 = maxsize;} //to the full function's scope.
    double largest_diff = GSL_NEGINF;
    double sum = 0;
    for (size_t i=0; i< maxsize1; i++){
        Apop_row(m1->data, i, arow);
        sum += m1->data->weights ? gsl_vector_get(m1->data->weights, i) : 1./maxsize1;
        largest_diff = GSL_MAX(largest_diff, fabs(sum-apop_cdf(arow, m2)));
    }
    if (m2_is_pmf){
        double sum = 0;
        for (size_t i=0; i< maxsize2; i++){     //There could be matched data rows to m1, so there is redundancy.
            Apop_row(m2->data, i, arow);   // Feel free to submit a smarter version.
            sum += m2->data->weights ? gsl_vector_get(m2->data->weights, i) : 1./maxsize2;
            largest_diff = GSL_MAX(largest_diff, fabs(sum-apop_cdf(arow, m2)));
        }
    }
    apop_data *out = apop_data_alloc();
    Asprintf(&out->names->title, "Kolmogorov-Smirnov test");
    apop_data_add_named_elmt(out, "max distance", largest_diff);
    double ps = m2_is_pmf ? psmirnov2x(largest_diff, maxsize1, maxsize2)
                          : kolmogorov_2x(maxsize1, largest_diff);
    apop_data_add_named_elmt(out, "p value, 2 tail", 1-ps);
    apop_data_add_named_elmt(out, "confidence, 2 tail", ps);
    return out;
}

/** Create a histogram from data by putting data into bins of fixed width. 

\param indata The input data that will be binned. This is copied and the copy will be modified.
\param close_top_bin Normally, a bin covers the range from the point equal to its minimum to points strictly less than
the minimum plus the width.  if \c 'y', then the top bin includes points less than or equal to the upper bound. This solves the problem of displaying histograms where the top bin is just one point.
\param binspec This is an \ref apop_data set with the same number of columns as \c indata. 
If you want a fixed size for the bins, then the first row of the bin spec is the bin width for each column.
This allows you to specify a width for each dimension, or specify the same size for all with something like:

\param bin_count If you don't provide a bin spec, I'll provide this many evenly-sized bins. Default: \f$\sqrt(N)\f$.  \code
Apop_row(indata, 0, firstrow);
apop_data *binspec = apop_data_copy(firstrow);
gsl_matrix_set_all(binspec->matrix, 10); //bins of size 10 for all dim.s
apop_data_to_bins(indata, binspec);
\endcode
The presumption is that the first bin starts at zero in all cases. You can add a second row to the spec to give the offset for each dimension.  Default: NULL. if no binspec and no binlist, then a grid with offset equal to the min of the column, and bin size such that it takes \f$\sqrt{N}\f$ bins to cover the range to the max element. 


\return A pointer to a binned \ref apop_data set.  If you didn't give me a binspec, then I attach one to the output set as a page named \c \<binspec\>, so you can snap a second data set to the same grid using 
\code
apop_data_to_bins(first_set, NULL);
apop_data_to_bins(second_set, apop_data_get_page(first_set, "<binspec>"));
\endcode


  The text segment, if any, is not binned. I use \ref apop_data_pmf_compress as the final step in the binning, 
  and that does respect the text segment. 

Here is a sample program highlighting the difference between \ref apop_data_to_bins and \ref apop_data_pmf_compress .

\include binning.c
*/
#ifdef APOP_NO_VARIADIC
apop_data * apop_data_to_bins(apop_data *indata, apop_data *binspec, int bin_count, char close_top_bin){
#else
apop_varad_head(apop_data *, apop_data_to_bins){
    apop_data *apop_varad_var(indata, NULL);
    Apop_assert_c(indata, NULL, 1, "NULL input data set, so returning NULL output data set.");
    apop_data *apop_varad_var(binspec, NULL);
    char apop_varad_var(close_top_bin, 'n');
    int apop_varad_var(bin_count, 0);
    return apop_data_to_bins_base(indata, binspec, bin_count, close_top_bin);
}

 apop_data * apop_data_to_bins_base(apop_data *indata, apop_data *binspec, int bin_count, char close_top_bin){
#endif
    Get_vmsizes(indata); //firstcol, vsize, msize1, msize2
    double binwidth, offset, max=0;
    apop_data *out = apop_data_copy(indata);
    apop_data *bs = binspec ? binspec
                    : apop_data_add_page(out, 
                        apop_data_alloc(vsize? 2: 0, msize1? 2: 0, indata->matrix ? msize2: 0),
                        "<binspec>");
    for (int j= firstcol; j< msize2; j++){
        Apop_col_v(out, j, onecol);
        if (binspec){
           binwidth = apop_data_get(binspec, 0, j);
           offset = ((binspec->vector && binspec->vector->size==2 )
                   ||(binspec->matrix && binspec->matrix->size1==2)) ? apop_data_get(binspec, 1, j) : 0;
        } else {
            Apop_col_v(bs, j, abin);
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
    apop_data_pmf_compress(out);
    return out;
}
