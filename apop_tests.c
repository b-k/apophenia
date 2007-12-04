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
apop_data * apop_test_ANOVA(apop_data *d){
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
