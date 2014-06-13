
/** \file apop_tests.c	 */
/* Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
 
At the moment, the header for  apop_test_anova is in \c asst.h.
 */
#include "apop_internal.h"

static apop_data * produce_t_test_output(int df, double stat, double diff){
    double pval, qval, two_tail;
    if (!gsl_isnan(stat)){
        pval    = gsl_cdf_tdist_P(stat, df);
        qval    = 1-pval;
        two_tail= 2*GSL_MIN(fabs(pval-.5),fabs(qval-0.5));
    } else {
        pval    = GSL_NAN;
        qval    = GSL_NAN;
        two_tail= GSL_NAN;
    }
    apop_data *out = apop_data_alloc();
    apop_data_add_named_elmt(out, "mean left - right", diff);
    apop_data_add_named_elmt(out, "t statistic", stat);
    apop_data_add_named_elmt(out, "df", df);
    apop_data_add_named_elmt(out, "p value, 1 tail", GSL_MIN(pval,qval));
    apop_data_add_named_elmt(out, "confidence, 1 tail", 1 - GSL_MIN(pval,qval));
    apop_data_add_named_elmt(out, "p value, 2 tail", 1- two_tail);
    apop_data_add_named_elmt(out, "confidence, 2 tail", two_tail);
    return out;
}

/** Answers the question: with what confidence can I say that the means of these two columns of data are different?
<tt>apop_paired_t_test</tt> answers the question: with what confidence can I say that the mean difference between the two columns is zero?

If \c apop_opts.verbose is >=1, then display some information to stdout, like the mean/var/count for both vectors and the t statistic.

\ingroup ttest
\param {a, b} two columns of data
\return an \ref apop_data set with the following elements:
    mean left - right:    the difference in means; if positive, first vector has larger mean, and one-tailed test is testing \f$L > R\f$, else reverse if negative.<br>
    t statistic:    used for the test<br>
    df:             degrees of freedom<br>
    p value, 1 tail: the p-value for a one-tailed test that one vector mean is greater than the other.
    confidence, 1 tail: 1- p value.
    p value, 2 tail: the p-value for the two-tailed test that left mean = right mean.
    confidence, 2 tail: 1-p value
*/
apop_data *	apop_t_test(gsl_vector *a, gsl_vector *b){
    int a_count = a->size,
        b_count = b->size;
    double a_avg = apop_vector_mean(a);
    double a_var = (a_count > 1) ? apop_vector_var(a) : 0,
           b_avg = apop_vector_mean(b),
    b_var = b_count > 1 ? apop_vector_var(b): 0,
    stat = (a_avg - b_avg)/ sqrt(
                        (b_count > 1 ? b_var/(b_count-1) : 0) 
                        + (a_count > 1 ? a_var/(a_count-1) : 0) 
                        );
    if (apop_opts.verbose >=1){
        printf("1st avg: %g; 1st std dev: %g; 1st count: %i.\n", a_avg, sqrt(a_var), a_count);
        printf("2st avg: %g; 2st std dev: %g; 2nd count: %i.\n", b_avg, sqrt(b_var), b_count);
        printf("t-statistic: %g.\n", stat);
    }
    int df = a_count+b_count-2;
    return produce_t_test_output(df, stat, a_avg - b_avg);
}

/** Answers the question: with what confidence can I say that the mean difference between the two columns is zero?

If \c apop_opts.verbose is >=2, then display some information, like the mean/var/count for both vectors and the t statistic, to stderr.

\ingroup ttest
\param {a, b} two columns of data
\return an \ref apop_data set with the following elements:
    mean left - right:    the difference in means; if positive, first vector has larger mean, and one-tailed test is testing \f$L > R\f$, else reverse if negative.<br>
    t statistic:    used for the test<br>
    df:             degrees of freedom<br>
    p value, 1 tail: the p-value for a one-tailed test that one vector mean is greater than the other.
    confidence, 1 tail: 1- p value.
    p value, 2 tail: the p-value for the two-tailed test that left mean = right mean.
    confidence, 2 tail: 1-p value
*/
apop_data * apop_paired_t_test(gsl_vector *a, gsl_vector *b){
    gsl_vector *diff = gsl_vector_alloc(a->size);
    gsl_vector_memcpy(diff, a);
    gsl_vector_sub(diff, b);
    int count = a->size; 
    double avg = apop_vector_mean(diff),
    var = apop_vector_var(diff),
    stat = avg/ sqrt(var/(count-1));
    gsl_vector_free(diff);
    Apop_notify(2, "avg diff: %g; diff std dev: %g; count: %i; t-statistic: %g.\n", avg, sqrt(var), count, stat);
    return produce_t_test_output(count-1, stat, avg);
}

/** Runs an F-test specified by \c q and \c c. Your best bet is to see
 the chapter on hypothesis testing in  <a href="http://modelingwithdata.org">Modeling With Data</a>, p 309. It will tell you that:
 \f[{N-K\over q}
 {({\bf Q}'\hat\beta - {\bf c})' [{\bf Q}' ({\bf X}'{\bf X})^{-1} {\bf Q}]^{-1} ({\bf Q}' \hat\beta - {\bf c})
 \over {\bf u}' {\bf u} } \sim F_{q,N-K},\f]
 and that's what this function is based on.

 \param est     an \ref apop_model that you have already calculated. (No default)
 \param contrast       The matrix \f${\bf Q}\f$ and the vector \f${\bf c}\f$, where each row represents a hypothesis. (Defaults: if matrix is \c NULL, it is set to the identity matrix with the top row missing. If the vector is \c NULL, it is set to a zero matrix of length equal to the height of the contrast matrix. Thus, if the entire \c apop_data set is NULL or omitted, we are testing the hypothesis that all but \f$\beta_1\f$ are zero.)
 \return An \c apop_data set with a few variants on the confidence with which we can reject the joint hypothesis.
 \todo There should be a way to get OLS and GLS to store \f$(X'X)^{-1}\f$. In fact, if you did GLS, this is invalid, because you need \f$(X'\Sigma X)^{-1}\f$, and I didn't ask for \f$\Sigma\f$.

\li There are two approaches to an \f$F\f$-test: the ANOVA approach, which is typically built around the claim that all effects but the mean are zero; and the more general regression form, which allows for any set of linear claims about the data. If you send a \c NULL contrast set, I will generate the set of linear contrasts that are equivalent to the ANOVA-type approach. Readers of {\em Modeling with Data}, note that there's a bug in the book that claims that the traditional ANOVA approach also checks that the coefficient for the constant term is also zero; this is not the custom and doesn't produce the equivalence presented in that and other textbooks.

\exception out->error='a'  Allocation error.
\exception out->error='d'  dimension-matching error.
\exception out->error='i'  matrix inversion error.
\exception out->error='m'  GSL math error.
\li This function uses the \ref designated syntax for inputs.
*/
#ifdef APOP_NO_VARIADIC
apop_data * apop_f_test(apop_model *est, apop_data *contrast){
#else
apop_varad_head(apop_data *, apop_f_test){
    apop_model *apop_varad_var(est, NULL)
    Nullcheck_m(est, NULL);
    Nullcheck_d(est->data, NULL);
    apop_data * apop_varad_var(contrast, NULL);
    int free_data=0,free_matrix=0,free_vector=0;
    if (!contrast) contrast = apop_data_alloc(),free_data=1;
    if (!contrast->matrix) {
        int size = est->parameters->vector->size;
        contrast->matrix= gsl_matrix_calloc(size - 1, size);
        for (int i=1; i< size; i++)
            apop_data_set(contrast, i-1, i, 1);
    }
    if (!contrast->vector) contrast->vector = gsl_vector_calloc(contrast->matrix->size1),free_vector=1;

    apop_data *out = apop_f_test_base(est, contrast);
    if (free_data) {apop_data_free(contrast); return out;}
    if (free_matrix) gsl_matrix_free(contrast->matrix);
    if (free_vector) gsl_vector_free(contrast->vector);
    return out;
    return apop_f_test_base(est, contrast);
}

 apop_data * apop_f_test_base(apop_model *est, apop_data *contrast){
#endif
    apop_data *out = apop_data_alloc();
    Asprintf(&out->names->title, "F test");
    size_t contrast_ct = contrast->vector->size;
    Apop_stopif(contrast->matrix->size1 != contrast_ct,  out->error='d'; return out,
            0, "I counted %zu contrasts by the size of either contrast->vector or "
            "est->parameters->vector->size, but you gave me a matrix with %zu rows. Those should match."
            , contrast_ct, contrast->matrix->size1);
    double f_stat, pval;

    Get_vmsizes(est->data); //msize1, msize2
    int data_df = msize1 - contrast_ct;

    //Find (\underbar x)'(\underbar x), where (\underbar x) = the data with means removed
    long double means[msize2];
    for (int i=1; i< msize2; i++){
        Apop_col_v(est->data, i, onecol)
        means[i] = apop_vector_mean(onecol);
    }
    means[0]=0;// don't screw with the ones column.
    apop_data *xpx = apop_data_alloc(msize2, msize2);
    Apop_stopif(xpx->error, apop_data_free(xpx); out->error='a'; return out, 0, "allocation error");
    for (int i=0; i< msize2; i++)
        for (int j=0; j< msize2; j++){ //at this loop, we calculate one cell in the dot prouct
            long double total = 0;
            for (int c=0; c<msize1; c++)
                total +=  (gsl_matrix_get(est->data->matrix, c, i) - means[i])
                         *(gsl_matrix_get(est->data->matrix, c, j) - means[j]);
            apop_data_set(xpx, i, j, total);
        }

    apop_data xpxinv = (apop_data){.matrix=apop_matrix_inverse(xpx->matrix)};
    Apop_stopif(!xpxinv.matrix, out->error='i'; return out, 0, "inversion of X'X error");
    apop_data *qprimexpxinv = apop_dot(contrast, &xpxinv, 'm', 'm');
    apop_data *qprimexpxinvq = apop_dot(qprimexpxinv, contrast, 'm', 't');
    Apop_stopif(qprimexpxinvq->error || qprimexpxinv->error, out->error='m'; return out, 0, "broken dot");
    apop_data qprimexpxinvqinv = (apop_data){.matrix=apop_matrix_inverse(qprimexpxinvq->matrix)};
    Apop_stopif(!qprimexpxinvqinv.matrix, out->error='i'; return out, 0, "inversion of Q'(X'X)^{-1}Q error");
    apop_data_free(qprimexpxinvq);
    apop_data_free(qprimexpxinv);
    apop_data *qprimebeta = apop_dot(contrast, est->parameters, 'm', 'v');
    Apop_stopif(qprimebeta->error, out->error='m'; return out, 0, "broken dot");
    gsl_vector_sub(qprimebeta->vector, contrast->vector);
    apop_data *qprimebetaminusc_qprimexpxinvqinv = apop_dot(&qprimexpxinvqinv, qprimebeta, .form2='v');
    Apop_stopif(qprimebetaminusc_qprimexpxinvqinv->error, out->error='m'; return out, 0, "broken dot");
    gsl_blas_ddot(qprimebeta->vector, qprimebetaminusc_qprimexpxinvqinv->vector, &f_stat);
    apop_data_free(xpx);
    apop_data_free(qprimebeta);
    apop_data_free(qprimebetaminusc_qprimexpxinvqinv);

    apop_data *r_sq_list = apop_estimate_coefficient_of_determination (est);
    double variance = apop_data_get(r_sq_list, .rowname="sse");
    f_stat *=  data_df / (variance * contrast_ct);
    pval    = (contrast_ct > 0 && data_df > 0) ? gsl_cdf_fdist_Q(f_stat, contrast_ct, data_df): GSL_NAN; 

    apop_data_add_named_elmt(out, "F statistic", f_stat);
    apop_data_add_named_elmt(out, "p value", pval);
    apop_data_add_named_elmt(out, "confidence", 1- pval);
    apop_data_add_named_elmt(out, "df1", contrast_ct);
    apop_data_add_named_elmt(out, "df2", data_df);
    return out;
}

static double one_chi_sq(apop_data *d, int row, int col, int n){
    Apop_row_v(d, row, vr);
    Apop_col_v(d, col, vc);
    double rowexp  = apop_vector_sum(vr)/n;
    double colexp  = apop_vector_sum(vc)/n;
    double observed = apop_data_get(d, row, col);
    double expected = n * rowexp * colexp;
    return gsl_pow_2(observed - expected)/expected; 
}

/** Run a Chi-squared test on an ANOVA table, i.e., an NxN table with the null hypothesis that all cells are equally likely.

 \param d The input data, which is a crosstab of various elements. They don't have to sum to one.
 \see apop_test_fisher_exact
 \ingroup asst_tests
 */
apop_data * apop_test_anova_independence(apop_data *d){
    Apop_stopif(!d || !d->matrix, return NULL, 0, "You sent me data with no matrix element. Returning NULL.");
    double total = 0;
    //You can have a one-column or one-row matrix if you want; else df = (rows-1)*(cols-1)
    double df = d->matrix->size1==1 ? d->matrix->size2-1 : d->matrix->size2 == 1 ? d->matrix->size1 
                              : (d->matrix->size1 - 1)* (d->matrix->size2 - 1);
    Apop_stopif(!df, return NULL, 0, "You sent a degenerate matrix. Returning NULL.");
    int n = apop_matrix_sum(d->matrix);
    for (size_t row=0; row <d->matrix->size1; row++)
        for (size_t col=0; col <d->matrix->size2; col++)
            total += one_chi_sq(d, row, col, n);
    apop_data *out = apop_data_alloc();
    double chisq = gsl_cdf_chisq_Q(total, df);
    apop_data_add_named_elmt(out, "chi squared statistic", total);
    apop_data_add_named_elmt(out, "df", df);
    apop_data_add_named_elmt(out, "p value", chisq);
    return out;
}


static apop_data* apop_anova_one_way(char *table, char *data, char *grouping){
    //ANOVA has always just been a process of filling in a form, and
    //that's what this function does.
    apop_data *out = apop_data_calloc(3, 6);
    apop_name_add(out->names, "sum of squares", 'c');
    apop_name_add(out->names, "df", 'c');
    apop_name_add(out->names, "mean squares", 'c');
    apop_name_add(out->names, "F ratio", 'c');
    apop_name_add(out->names, "p value", 'c');
    apop_name_add(out->names, "confidence", 'c');
    apop_name_add(out->names, grouping, 'r');
    apop_name_add(out->names, "residual", 'r');
    apop_name_add(out->names, "total", 'r');
 
    //total sum of squares:
    apop_data* tss = apop_query_to_data("select var_pop(%s), count(*) from %s", data, table);
    Apop_stopif(!tss, apop_return_data_error('q'), 0, "Query 'select var_pop(%s), count(*) from %s' returned NULL. Does that look right to you?", data, table);
    apop_data_set(out, 2, 0, apop_data_get(tss, 0, 0)*apop_data_get(tss, 0, 1)); //total sum of squares
    double total_df = apop_data_get(tss, 0, 1);
    apop_data_set(out, 2, 1, apop_data_get(tss, 0, 1)); //total df.

    //within group sum of squares:
    apop_data* wss = apop_query_to_data("select var_pop(%s), count(*) from %s group by %s", data,  table, grouping);
    double group_df = wss->matrix->size1-1;
    apop_data_set(out, 0, 0, apop_data_get(wss, 0, 0)*group_df); //within sum of squares
    apop_data_set(out, 0, 1, group_df);

    //residuals are just total-wss
    apop_data_set(out, 1, 0, apop_data_get(out, 2, 0) - apop_data_get(out, 0,0)); //residual sum of squares
    double residual_df = total_df - group_df;
    apop_data_set(out, 1, 1, residual_df); //residual df

    apop_data_set(out, 0, 2, apop_data_get(out, 0, 0)/apop_data_get(out, 0, 1));//mean SS within
    apop_data_set(out, 1, 2, apop_data_get(out, 1, 0)/apop_data_get(out, 1, 1));//mean SS residual

    apop_data_set(out, 0, 3, apop_data_get(out, 0, 2)/apop_data_get(out, 1, 2));//F ratio
    apop_data_set(out, 0, 4, gsl_cdf_fdist_P(apop_data_get(out, 0, 3), group_df, residual_df));//pval
    apop_data_set(out, 0, 5, 1- apop_data_get(out, 0, 4));//confidence

    apop_data_free(tss);
    apop_data_free(wss);
    return out;
}

/** This function produces a traditional one- or two-way ANOVA table. It
  works from data in an SQL table, using queries of the form <tt>select
  data from table group by grouping1, grouping2</tt>.

  \param table The table to be queried. Anything that can go in an SQL <tt>from</tt> clause is OK, so this can be a plain table name or a temp table specification like <tt>(select ... )</tt>, with parens.
  \param data The name of the column holding the count or other such data
  \param grouping1 The name of the first column by which to group data
  \param grouping2 If this is \c NULL, then the function will return a one-way ANOVA. Otherwise, the name of the second column by which to group data in a two-way ANOVA.
 */
#ifdef APOP_NO_VARIADIC
apop_data* apop_anova(char *table, char *data, char *grouping1, char *grouping2){
#else
apop_varad_head(apop_data*, apop_anova){
    char *apop_varad_var(table, NULL)
    Apop_stopif(!table, return NULL, 0, "I need the name of a table in the SQL database.");
    if (!strchr(table, ')')) //if you found ()s,then it is a temp table spec.
        Apop_stopif(!apop_table_exists(table), return NULL, 0, "I couldn't find the table %s in the database.", table);
    char *apop_varad_var(data, NULL)
    Apop_stopif(!data, return NULL, 0, "I need the name of the column in the %s table with the count or other data.", table);
    char *apop_varad_var(grouping1, NULL)
    Apop_stopif(!data, return NULL, 0, "I need at least grouping1, a column in the %s table.", table);
    char *apop_varad_var(grouping2, NULL)
    return apop_anova_base(table, data, grouping1, grouping2);
}

 apop_data* apop_anova_base(char *table, char *data, char *grouping1, char *grouping2){
#endif
    apop_data *first = apop_anova_one_way(table, data, grouping1);
    Apop_stopif(first->error, return first, 0, "Error (%c) running one-way ANOVA.", first->error);
    if (!grouping2) return first;
    apop_data *second = apop_anova_one_way(table, data, grouping2);
    char *joined = NULL;
    Asprintf(&joined, "%s, %s", grouping1, grouping2);
    apop_data *interaction = apop_anova_one_way(table, data, joined);
    apop_data *out = apop_data_calloc(5, 6);
    apop_name_stack(out->names, first->names, 'c');
    apop_data_add_names(out, 'r', first->names->row[0], second->names->row[0],
                                  "interaction", "residual", "total");
    Apop_row_v(first, 0, firstrow);
    Apop_row_v(second, 0, secondrow);
    Apop_row_v(interaction, 0, interrow);
    Apop_row_v(first, 2, totalrow);
    gsl_matrix_set_row(out->matrix, 0, firstrow);
    gsl_matrix_set_row(out->matrix, 1, secondrow);
    gsl_matrix_set_row(out->matrix, 2, interrow);
    gsl_matrix_set_row(out->matrix, 4, totalrow);
    
    //residuals are just total-wss
    apop_data_set(out, 3, 0, apop_data_get(out, 4, 0) 
            - gsl_vector_get(firstrow, 0) - gsl_vector_get(secondrow, 0) - gsl_vector_get(interrow, 0)); //residual sum of squares
    double residual_df = apop_data_get(out, 4, 1) 
            - gsl_vector_get(firstrow, 1) - gsl_vector_get(secondrow, 1) - gsl_vector_get(interrow, 1); //residual df
    apop_data_set(out, 3, 1, residual_df);

    apop_data_set(out, 3, 2, apop_data_get(out, 3, 0)/apop_data_get(out, 3, 1));//mean SS residual

    apop_data_set(out, 0, 3, apop_data_get(out, 0, 2)/apop_data_get(out, 3, 2));//F ratio
    apop_data_set(out, 0, 4, gsl_cdf_fdist_P(apop_data_get(out, 0, 3), gsl_vector_get(firstrow, 1), residual_df));//pval
    apop_data_set(out, 0, 5, 1- apop_data_get(out, 0, 4));//confidence

    apop_data_set(out, 1, 3, apop_data_get(out, 1, 2)/apop_data_get(out, 3, 2));//F ratio
    apop_data_set(out, 1, 4, gsl_cdf_fdist_P(apop_data_get(out, 1, 3), gsl_vector_get(secondrow, 1), residual_df));//pval
    apop_data_set(out, 1, 5, 1- apop_data_get(out, 1, 4));//confidence

    apop_data_set(out, 2, 3, apop_data_get(out, 2, 2)/apop_data_get(out, 3, 2));//F ratio
    apop_data_set(out, 2, 4, gsl_cdf_fdist_P(apop_data_get(out, 2, 3), gsl_vector_get(interrow, 1), residual_df));//pval
    apop_data_set(out, 2, 5, 1- apop_data_get(out, 2, 4));//confidence

    free(joined);
    apop_data_free(first);
    apop_data_free(second);
    apop_data_free(interaction);
    return out;
}

/** This is a convenience function to do the lookup of a given statistic along a given
distribution. You give me a statistic, its (hypothesized) distribution, and whether
to use the upper tail, lower tail, or both. I will return the odds of a Type I error
given the model---in statistician jargon, the \f$p\f$-value.  [Type I error: odds of
rejecting the null hypothesis when it is true.]

   For example, 
   \code
   apop_test(1.3);
   \endcode

will return the density of the standard Normal distribution that is more than 1.3 from zero.  
If this function returns a small value, we can be confident that the statistic is significant. Or, 
   \code
   apop_test(1.3, "t", 10, tail='u');
   \endcode

will give the appropriate odds for an upper-tailed test using the \f$t\f$-distribution with 10 degrees of freedom (e.g., a \f$t\f$-test of the null hypothesis that the statistic is less than or equal to zero).

Several more distributions are supported; see below.

 \li For a two-tailed test (the default), this returns the density outside the range. I'll only do this for symmetric distributions.
 \li For an upper-tail test ('u'), this returns the density above the cutoff
 \li For a lower-tail test ('l'), this returns the density below the cutoff 
 
\param statistic    The scalar value to be tested.  
\param distribution  The name of the distribution; see below.
\param p1  The first parameter for the distribution; see below.
\param p2  The second parameter for the distribution; see below.
\param tail 'u' = upper tail; 'l' = lower tail; anything else = two-tailed. (default = two-tailed)

\return The odds of a Type I error given the model (the \f$p\f$-value).

Here is a list of distributions you can use, and their parameters.

\c "normal" or \c "gaussian" 
\li p1=mu, p2=sigma
\li default (0, 1)

\c "lognormal"  
\li p1=mu, p2=sigma
\li default (0, 1) 
\li Remember, mu and sigma refer to the Normal one would get after exponentiation
\li One-tailed tests only

\c "uniform"  
\li p1=lower edge, p2=upper edge
\li default (0, 1)
\li two-tailed tests are run relative to the center, (p1+p2)/2.

\c "t"  
\li p1=df
\li no default

\c "chi squared", \c "chi", \c "chisq": 
\li p1=df
\li no default
\li One-tailed tests only; default='u' (\f$p\f$-value for typical cases)

\c "f"  
\li p1=df1, p2=df2
\li no default
\li One-tailed tests only

\li This function uses the \ref designated syntax for inputs.
*/
#ifdef APOP_NO_VARIADIC
double apop_test(double statistic, char *distribution, double p1, double p2, char tail){
#else
apop_varad_head(double, apop_test){
    double  apop_varad_var(statistic, 0);
    char*  apop_varad_var(distribution, NULL);
    double apop_varad_var(p1, 0);
    double apop_varad_var(p2, 0);
    int is_chi = strcasecmp(distribution, "chi squared")|| strcasecmp(distribution, "chi")
                     || strcasecmp(distribution, "chisq");
     Apop_stopif(!strcasecmp(distribution, "f") && (!p1 || !p2), return NAN, 0, "I need both a p1 and p2 parameter specifying the degrees of freedom.");
     Apop_stopif((!strcasecmp(distribution, "t") || !strcasecmp(distribution, "f") || is_chi)
             && !p1, return NAN, 0, "I need a p1 parameter specifying the degrees of freedom.");
     if (!p2 && (!distribution || !strcasecmp(distribution, "normal") || !strcasecmp(distribution, "gaussian") ))
         p2 = 1;
     if (!p2 && p1 >= 0 && !strcasecmp(distribution, "uniform"))
         p2 = 1;

    char apop_varad_var(tail, 0);
    if (!tail) tail = is_chi ? 'u' : 'a';
    return apop_test_base(statistic, distribution, p1, p2, tail);
}

 double apop_test_base(double statistic, char *distribution, double p1, double p2, char tail){
#endif
    //This is a long and boring function. I am aware that there are
    //clever way to make it shorter.
     if (!distribution || !strcasecmp(distribution, "normal") || !strcasecmp(distribution, "gaussian") ){
         if (tail == 'u')
             return gsl_cdf_gaussian_Q(p1-statistic, p2);
         else if (tail == 'l')
             return gsl_cdf_gaussian_P(p1-statistic, p2);
         else
             return 2 * gsl_cdf_gaussian_Q(fabs(p1-statistic), p2);
     }
    else if (!strcasecmp(distribution, "lognormal")){
         if (tail == 'u')
             return gsl_cdf_lognormal_Q(statistic, p1, p2);
         else if (tail == 'l')
             return gsl_cdf_lognormal_P(statistic, p1, p2);
         else
             Apop_assert(0, "A two-tailed test doesn't really make sense for the lognormal. Please specify either tail= 'u' or tail= 'l'.");
     }
    else if (!strcasecmp(distribution, "t")){
         if (tail == 'u')
             return gsl_cdf_tdist_Q(statistic, p1);
         else if (tail == 'l')
             return gsl_cdf_tdist_P(statistic, p1);
         else
             return 2 * gsl_cdf_tdist_Q(fabs(statistic), p1);
     }
    else if (!strcasecmp(distribution, "f")){
         if (tail == 'u')
             return gsl_cdf_fdist_Q(statistic, p1, p2);
         else if (tail == 'l')
             return gsl_cdf_fdist_P(statistic, p1, p2);
         else
             Apop_assert(0, "A two-tailed test doesn't really make sense for the %s. Please specify either tail= 'u' or tail= 'l'.", distribution);
     }
    else if (!strcasecmp(distribution, "chi squared")|| !strcasecmp(distribution, "chi")
                                                || !strcasecmp(distribution, "chisq")){
         if (tail == 'u')
             return gsl_cdf_chisq_Q(statistic, p1);
         else if (tail == 'l')
             return gsl_cdf_chisq_P(statistic, p1);
         else
             Apop_assert(0, "A two-tailed test doesn't really make sense for the %s. Please specify either tail= 'u' or tail= 'l'.", distribution);
     }
    else if (!strcasecmp(distribution, "uniform")){
         if (tail == 'u')
             return gsl_cdf_flat_Q(statistic, p1, p2);
         else if (tail == 'l')
             return gsl_cdf_flat_P(statistic, p1, p2);
         else
             return 2 * gsl_cdf_flat_Q(fabs(statistic - (p1+p2)/2.), p1, p2);
     }
    Apop_assert(0, "Sorry, but I don't recognize %s as a distribution", distribution);
}
