/** \file apop_tests.c	ANOVA.\n
 \author Ben Klemens
Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
 
At the moment, the header for  apop_test_ANOVA is in \c asst.h.
 */
#include <asst.h>
#include <types.h>
#include <stats.h>

static double one_chi_sq(apop_data *d, int row, int col, int n){
    APOP_ROW(d, row, vr);
    APOP_COL(d, col, vc);
    double rowexp  = apop_vector_sum(vr)/n;
    double colexp  = apop_vector_sum(vc)/n;
    double observed = apop_data_get(d, row, col);
    double expected = n * rowexp * colexp;
    return gsl_pow_2(observed - expected)/expected; 
}

/** Run a Chi-squared test on an ANOVA table, i.e., an NxN table with
 the null hypothesis that all cells are equally likely.

 \param d The input data, which is a crosstab of various elements. They don't have to sum to one.
 \ingroup asst_tests
 */
apop_data * apop_test_ANOVA_independence(apop_data *d){
  double total = 0;
  size_t row, col;
  if (!d || !d->matrix) {
      apop_error(0,'c', "%s: You sent me data with no matrix element. Returning NULL.\n", __func__);
      return NULL;
  }
  //You can have a one-column or one-row matrix if you want; else df = (rows-1)*(cols-1)
  double df    = d->matrix->size1==1 ? d->matrix->size2-1 : d->matrix->size2 == 1 ? d->matrix->size1 
                            : (d->matrix->size1 - 1)* (d->matrix->size2 - 1);
  if (!df) {
      apop_error(0,'c', "%s: You sent a degenerate matrix. Returning NULL.\n", __func__);
      return NULL;
  }
  int    n     = apop_matrix_sum(d->matrix);
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





static apop_data* apop_ANOVA_one_way(char *table, char *data, char *grouping){
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
apop_data* apop_ANOVA(char *table, char *data, char *grouping1, char *grouping2){
    apop_data *first = apop_ANOVA_one_way(table, data, grouping1);
    if (!grouping2)
        return first;
    apop_data *second = apop_ANOVA_one_way(table, data, grouping2);
    char *joined = NULL;
    asprintf(&joined, "%s, %s", grouping1, grouping2);
    apop_data *interaction = apop_ANOVA_one_way(table, data, joined);
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
