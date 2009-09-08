/** \file apop_tests.c	 */
/* Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
 
At the moment, the header for  apop_test_anova is in \c asst.h.
 */
#include "asst.h"
#include "types.h"
#include "stats.h"
#include "conversions.h"

static double one_chi_sq(apop_data *d, int row, int col, int n){
    APOP_ROW(d, row, vr);
    APOP_COL(d, col, vc);
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
  apop_assert(d && d->matrix, NULL, 0,'c', "You sent me data with no matrix element. Returning NULL.");
  double total = 0;
  size_t row, col;
  //You can have a one-column or one-row matrix if you want; else df = (rows-1)*(cols-1)
  double df    = d->matrix->size1==1 ? d->matrix->size2-1 : d->matrix->size2 == 1 ? d->matrix->size1 
                            : (d->matrix->size1 - 1)* (d->matrix->size2 - 1);
  apop_assert(df, NULL, 0,'c', "You sent a degenerate matrix. Returning NULL.");
  int n     = apop_matrix_sum(d->matrix);
    for (row=0; row <d->matrix->size1; row++)
        for (col=0; col <d->matrix->size2; col++)
            total += one_chi_sq(d, row, col, n);
    apop_data *out = apop_data_alloc(3,0,0);
    double chisq   = gsl_cdf_chisq_Q(total, df);
    apop_data_add_named_elmt(out, "chi squared statistic", total);
    apop_data_add_named_elmt(out, "df", df);
    apop_data_add_named_elmt(out, "p value", chisq);
    return out;
}


static apop_data* apop_anova_one_way(char *table, char *data, char *grouping){
    //ANOVA has always just been a process of filling in a form, and
    //that's what this function does. I use apop_data_get instead of
    //apop_data_get_tt for efficiency reasons.
  apop_data *out = apop_data_calloc(0, 3, 6);
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
    apop_assert(tss, NULL, 0, 's', "Query 'select var_pop(%s), count(*) from %s' returned NULL. Does that look right to you?", data, table);
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

  \param table The table to be queried. 
  \param data The name of the column holding the data
  \param grouping1 The name of the first column by which to group data
  \param grouping2 If this is \c NULL, then the function will return a one-way ANOVA. Otherwise, the name of the second column by which to group data in a two-way ANOVA.
 */
apop_data* apop_anova(char *table, char *data, char *grouping1, char *grouping2){
    apop_data *first = apop_anova_one_way(table, data, grouping1);
    if (!grouping2)
        return first;
    apop_data *second = apop_anova_one_way(table, data, grouping2);
    char *joined = NULL;
    asprintf(&joined, "%s, %s", grouping1, grouping2);
    apop_data *interaction = apop_anova_one_way(table, data, joined);
    apop_data *out         = apop_data_calloc(0, 5, 6);
    apop_name_stack(out->names, first->names, 'c');
    apop_name_add(out->names, first->names->row[0], 'r');
    apop_name_add(out->names, second->names->row[0], 'r');
    apop_name_add(out->names, "interaction", 'r');
    apop_name_add(out->names, "residual", 'r');
    apop_name_add(out->names, "total", 'r');


    APOP_ROW(first, 0, firstrow);
    APOP_ROW(second, 0, secondrow);
    APOP_ROW(interaction, 0, interrow);
    APOP_ROW(first, 2, totalrow);
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


/**
This is a convenience function to do the lookup of a given statistic
along a given distribution. You give me a statistic, its (hypothesized)
distribution, and whether to use the
upper tail, lower tail, or both. I will return the odds of a Type I error given the model---in statistician jargon, the \f$p\f$-value.
[Type I error: odds of rejecting the null hypothesis when it is true.]

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
\li One-tailed tests only

\c "f"  
\li p1=df1, p2=df2
\li no default
\li One-tailed tests only

 */

APOP_VAR_HEAD double apop_test(double statistic, char *distribution, double p1, double p2, char tail){
    double  apop_varad_var(statistic, 0);
    char*  apop_varad_var(distribution, NULL);
    double apop_varad_var(p1, 0);
    double apop_varad_var(p2, 0);
     apop_assert(strcasecmp(distribution, "f") || p1, 0, 0, 's', "I need both a p1 and p2 parameter specifying the degrees of freedom.")
     apop_assert(strcasecmp(distribution, "t") || strcasecmp(distribution, "f")
             || strcasecmp(distribution, "chi squared")|| strcasecmp(distribution, "chi")
                                                || strcasecmp(distribution, "chisq")
             || p1, 0, 0, 's', "I need a p1 parameter specifying the degrees of freedom.")
     if (!p2 && (!distribution || !strcasecmp(distribution, "normal") || !strcasecmp(distribution, "gaussian") ))
         p2 = 1;
     if (!p2 && p1 >= 0 && !strcasecmp(distribution, "uniform"))
         p2 = 1;

    char apop_varad_var(tail, 'a');
    return apop_test_base(statistic, distribution, p1, p2, tail);
APOP_VAR_ENDHEAD
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
             apop_error(0, 's', "A two-tailed test doesn't really make sense for the lognormal. Please specify either tail= 'u' or tail= 'l'.");
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
             apop_error(0, 's', "A two-tailed test doesn't really make sense for the %s. Please specify either tail= 'u' or tail= 'l'.", distribution);
     }
    else if (!strcasecmp(distribution, "chi squared")|| !strcasecmp(distribution, "chi")
                                                || !strcasecmp(distribution, "chisq")){
         if (tail == 'u')
             return gsl_cdf_chisq_Q(statistic, p1);
         else if (tail == 'l')
             return gsl_cdf_chisq_P(statistic, p1);
         else
             apop_error(0, 's', "A two-tailed test doesn't really make sense for the %s. Please specify either tail= 'u' or tail= 'l'.", distribution);
     }
    else if (!strcasecmp(distribution, "uniform")){
         if (tail == 'u')
             return gsl_cdf_flat_Q(statistic, p1, p2);
         else if (tail == 'l')
             return gsl_cdf_flat_P(statistic, p1, p2);
         else
             return 2 * gsl_cdf_flat_Q(fabs(statistic - (p1+p2)/2.), p1, p2);
     }
    apop_error(0,'s', "Sorry, but I don't recognize %s as a distribution", distribution);
    return 0; //shutting up the compiler.
}
