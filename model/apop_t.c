/* apop_t.c: the t-distribution, for modeling purposes.
Copyright (c) 2009, 2013 by Ben Klemens.  Licensed under the GPLv2; see COPYING.  */
#include "apop_internal.h"

//There used to be a χ^2 and F model, but nobody used them and they were largely untested.
//They last appeared in commit 2b4715111704cee3a86fca1b16903c4408bdacb9 if you'd like to recover them.

static void apop_t_estimate(apop_data *d, apop_model *m){
    Apop_stopif(!d, m->error='d'; return, 0, "No data with which to count df. (the default estimation method)");
    Get_vmsizes(d); //vsize, msize1, msize2, tsize
    double vmu = vsize ? apop_mean(d->vector) : 0;
    double v_sum_sq = vsize ? apop_var(d->vector)*(vsize-1) : 0;
    double m_sum_sq = 0;
    double mmu = 0;
   if (msize1) {
       apop_matrix_mean_and_var(d->matrix, &mmu, &m_sum_sq);
       m_sum_sq *= msize1*msize2-1;
   }
    apop_data_add_named_elmt(m->parameters, "μ", (vmu *vsize + mmu * msize1*msize2)/tsize);
    apop_data_add_named_elmt(m->parameters, "σ", sqrt(((tsize-3.)/(tsize-1)) * (v_sum_sq + m_sum_sq)/(tsize-1))); 
    apop_data_add_named_elmt(m->parameters, "df", tsize-1);
    apop_data_add_named_elmt(m->info, "log likelihood", m->log_likelihood(d, m));
}

static double one_t(double in, void *params){ 
    double mu = ((double*)params)[0];
    double sigma = ((double*)params)[1];
    double df = ((double*)params)[2];
    return log(gsl_ran_tdist_pdf((in-mu)/sigma, df)); 
}

static long double apop_tdist_llike(apop_data *d, apop_model *m){ 
    Nullcheck_mpd(d, m, GSL_NAN);
    double *params = m->parameters->vector->data;
    double sigma = params[1];
    Get_vmsizes(d); //tsize
    return apop_map_sum(d, .fn_dp=one_t, .param=params) - tsize * log(sigma);
}

int apop_t_dist_draw(double *out, gsl_rng *r, apop_model *m){ 
    Nullcheck_mp(m, 1);
    double mu = m->parameters->vector->data[0];
    double sigma = m->parameters->vector->data[1];
    double df = m->parameters->vector->data[2];
    *out = gsl_ran_tdist(r, df) * sigma + mu;
    return 0;
}

static long double apop_t_dist_cdf(apop_data *in, apop_model *m){
    Nullcheck_mp(m, GSL_NAN);
    double val = in->vector ? apop_data_get(in, 0, -1) : apop_data_get(in, 0, 0);
    double mu = m->parameters->vector->data[0];
    double sigma = m->parameters->vector->data[1];
    double df = m->parameters->vector->data[2];
    return gsl_cdf_tdist_P ((val-mu)/sigma, df);
}

static long double apop_t_dist_constraint(apop_data *beta, apop_model *m){
    Staticdef(apop_data *, d_constr, apop_data_falloc((2,2,3),
                             0, 0, 1, 0,  //0 < sigma
                            .9, 0, 0, 1)); //.9 < df
    double out= apop_linear_constraint(m->parameters->vector, d_constr);
    return out;
}


/*\amodel apop_t_distribution The t distribution, primarily for descriptive purposes.

If you want to test a hypothesis, you probably don't need this, and should instead
use \ref apop_test.

Given a statistic of the right properties, you can find the odds that the statistic
is above or below a cutoff on the t-, F, or chi-squared distribution using the \ref
apop_test function.

In that world, those three distributions are actually parameter free. The data is
assumed to be normalized to be based on a mean zero, variance one process, you get
the degrees of freedom from the size of the data, and the distribution is fixed.

For modeling purposes, more could be done. For example, the t-distribution is a favorite
proxy for Normal-like situations where there are fat tails relative to the Normal
(i.e., high kurtosis). Or, you may just prefer not to take the step of normalizing
your data---one could easily rewrite the theorems underlying the t-distribution without
the normalizations.

In such a case, the researcher would not want to fix the \f$df\f$, because \f$df\f$
indicates the fatness of the tails, which has some optimal value given the data.
Thus, there are two modes of use for these distributions:

\li Parameterized, testing style: the degrees of freedom are determined
from the data, and all necessary normalizations are assumed. Thus, this code---

\code
apop_data *t_for_testing = apop_estimate(data, apop_t)
\endcode

---will return exactly the type of \f$t\f$-distribution one would use for testing. 

\li By removing the \c estimate method---
\code
apop_model *spare_t = apop_model_copy(apop_t);
spare_t->estimate = NULL;
apop_model *best_fitting_t = apop_estimate(your_data, spare_t);
\endcode
---I will find the best \f$df\f$ via maximum likelihood, which may be desirable for
to find the best-fitting model for descriptive purposes.

\adoc    Input_format     Unordered list of scalars in the matrix and/or vector.     
\adoc    Parameter_format  vector->data[0] = mu<br>
                            vector->data[1] = sigma<br>
                            vector->data[2] = df 
\adoc    Estimate_results  I'll just count elements and set \f$df = n-1\f$. If you set the \c estimate method to \c NULL, via MLE.
\adoc    settings   \ref apop_mle_settings, \ref apop_parts_wanted_settings   
*/

apop_model *apop_t_distribution  = &(apop_model){"t distribution", 3, .dsize=1, .estimate = apop_t_estimate, 
         .log_likelihood = apop_tdist_llike, .draw=apop_t_dist_draw, .cdf=apop_t_dist_cdf,
         .constraint=apop_t_dist_constraint };
