/** \file apop_t_f_chi.c	t, F, chi squared, and Wishart distributions. */
/* Copyright (c) 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

/** \defgroup tfchi t-, chi-squared, F-, Wishart distributions

Most of these distributions are typically used for testing purposes.  For such a situation, you don't need the models here.
Given a statistic of the right properties, you can find the odds that the statistic is above or below a cutoff on the t-, F, or chi-squared distribution using the \ref apop_test function. 

In that world, those three distributions are actually parameter free. The data is assumed to be normalized to be based on a mean zero, variance one process, you get the degrees of freedom from the size of the data, and the distribution is fixed.

For modeling purposes, more could be done. For example, the t-distribution is a favorite proxy for Normal-like situations where there are fat tails relative to the Normal (i.e., high kurtosis). Or, you may just prefer not to take the step of normalizing your data---one could easily rewrite the theorems underlying the t-distribution without the normalizations.

In such a case, the researcher would not want to fix the \f$df\f$, because \f$df\f$ indicates the fatness of the tails, which has some optimal value given the data. 
Thus, there are two modes of use for these distributions: 

\li Parameterized, testing style: the degrees of freedom are determined
from the data, and all necessary normalizations are assumed. Thus, this code---

\code
apop_data *t_for_testing = apop_estimate(data, apop_t)
\endcode

---will return exactly the type of \f$t\f$-distribution one would use for testing. 

\li Descriptive: Just add an \ref apop_mle_settings group, and I'll find the best \f$df\f$ via maximum likelihood.

\code
Apop_settings_add_group(&apop_t, apop_mle, data);
apop_data *t_for_description = apop_estimate(data, apop_t);
\endcode

\c df works for all four distributions here; \c df2 makes sense only for the \f$F\f$, 

For the Wishart, the degrees of freedom and covariance matrix are always estimated via MLE.

*/
#include "mapply.h"
#include "variadic.h"
#include "likelihoods.h"
#include "model.h"
#include "internal.h"
#include <gsl/gsl_eigen.h>

double df, df2; 
int len;

apop_data * get_df(apop_data *d){
    apop_data *out = apop_data_alloc(1,0,0);
    int df = (d->vector ? d->vector->size : 0)
            +(d->matrix ? d->matrix->size1 * d->matrix->size2  : 0)
            -1;
    apop_data_add_named_elmt(out, "df", df);
    return out;
}

apop_model* apop_t_chi_estimate(apop_data *d, apop_model *m){
    Nullcheck(d); Nullcheck_m(m);
    apop_mle_settings *s = Apop_settings_get_group(m, apop_mle);
    if (!s){
        Apop_assert(d, NULL, 0, 's', "No data with which to count df. (the default estimation method)");
        apop_model *out = apop_model_copy(*m);
        out->parameters = get_df(d);
        apop_data_add_named_elmt(out->info, "log likelihood", out->log_likelihood(d, out));
        return out;
    } else 
        return apop_maximum_likelihood(d, m);
}

apop_model* apop_fdist_estimate(apop_data *d, apop_model *m){
    Nullcheck(d); Nullcheck_m(m);
    apop_mle_settings *s = Apop_settings_get_group(m, apop_mle);
    if (!s){
        Apop_assert(d, NULL, 0, 's', "No data with which to count df. (the default estimation method)");
        apop_model *out = apop_model_copy(*m);
        out->parameters = apop_data_alloc(2,0,0);
        apop_data_add_named_elmt(out->parameters, "df", d->vector->size -1);
        apop_data_add_named_elmt(out->parameters, "df2", d->matrix->size1 * d->matrix->size2 -1);
        apop_data_add_named_elmt(out->info, "log likelihood", apop_f_distribution.log_likelihood(d, out));
        return out;
    } else
        return apop_maximum_likelihood(d, m);
}

static double one_f(double in, void *df_in){ 
    double *df = df_in; 
    return log(gsl_ran_fdist_pdf(in, df[0], df[1])); 
}

static double one_t(double in, void *df){ return log(gsl_ran_tdist_pdf(in, *(double*)df)); }
static double one_chisq(double in, void *df){ return log(gsl_ran_chisq_pdf(in, *(double*)df)); }

double apop_tdist_llike(apop_data *d, apop_model *m){ 
    Nullcheck(d); Nullcheck_m(m); Nullcheck_p(m);
    double df = m->parameters->vector->data[0];
    return apop_map_sum(d, .fn_dp=one_t, .param=&df);
}

double apop_chisq_llike(apop_data *d, apop_model *m){ 
    Nullcheck(d); Nullcheck_m(m); Nullcheck_p(m);
    double df = m->parameters->vector->data[0];
    return apop_map_sum(d, .fn_dp=one_chisq, .param =&df);
}

double apop_fdist_llike(apop_data *d, apop_model *m){ 
    Nullcheck(d); Nullcheck_m(m); Nullcheck_p(m);
    double df[2];
    df[0] = m->parameters->vector->data[0];
    df[1] = m->parameters->vector->data[1];
    return apop_map_sum(d, .fn_dp=one_f, .param =df);
}

void apop_t_dist_draw(double *out, gsl_rng *r, apop_model *m){ 
    Nullcheck_mv(m); Nullcheck_pv(m);
    *out = gsl_ran_tdist (r, m->parameters->vector->data[0]);
}

void apop_f_dist_draw(double *out, gsl_rng *r, apop_model *m){
    Nullcheck_mv(m); Nullcheck_pv(m);
    *out = gsl_ran_fdist (r, m->parameters->vector->data[0], m->parameters->vector->data[1]);
}

void apop_chisq_dist_draw(double *out, gsl_rng *r, apop_model *m){
    Nullcheck_mv(m); Nullcheck_pv(m);
    *out = gsl_ran_chisq (r, m->parameters->vector->data[0]);
}

/** The multivariate generalization of the Gamma distribution.
\f$
\Gamma_p(a)=
\pi^{p(p-1)/4}\prod_{j=1}^p
\Gamma\left[ a+(1-j)/2\right]. \f$

See also \ref apop_multivariate_lngamma, which is more numerically stable in most cases.
*/
double apop_multivariate_gamma(double a, double p){
    double out = pow(M_PI, p*(p-1.)/4.);
    double factor = 1;
    for (int i=1; i<=p; i++)
        factor *= gsl_sf_gamma(a+(1-i)/2.);
    return out * factor;
}

/** The log of the multivariate generalization of the Gamma; see also
 \ref apop_multivariate_gamma.
*/
double apop_multivariate_lngamma(double a, double p){
    double out = M_LNPI * p*(p-1.)/4.;
    double factor = 0;
    for (int i=1; i<=p; i++){
        factor += gsl_sf_lngamma(a+(1-i)/2.);
    }
    return out + factor;
}

static void find_eigens(gsl_matrix **subject, gsl_vector *eigenvals, gsl_matrix *eigenvecs){
   gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc((*subject)->size1);
   gsl_eigen_symmv(*subject, eigenvals, eigenvecs, w);
   gsl_eigen_symmv_free (w);
   gsl_matrix_free(*subject); *subject  = NULL;
}

static void diagonal_copy(gsl_vector *v, gsl_matrix *m, char in_or_out){
    gsl_vector_view dv = gsl_matrix_diagonal(m);
    if (in_or_out == 'i') 
        gsl_vector_memcpy(&(dv.vector), v);
    else  
        gsl_vector_memcpy(v, &(dv.vector));
}

static double diagonal_size(gsl_matrix *m){
    gsl_vector_view dv = gsl_matrix_diagonal(m);
    return apop_sum(&dv.vector);
}

static double biggest_elmt(gsl_matrix *d){ 
    return  GSL_MAX(fabs(gsl_matrix_max(d)), fabs(gsl_matrix_min(d)));
}

/** Test whether the input matrix is positive semidefinite.

A covariance matrix will always be PSD, so this function can tell you whether your matrix is a valid covariance matrix.

Consider the 1x1 matrix in the upper left of the input, then the 2x2 matrix in the upper left, on up to the full matrix. If the matrix is PSD, then each of these has a positive determinant. This function thus calculates \f$N\f$ determinants for an \f$N\f$x\f$N\f$ matrix.

\param m The matrix to test. If \c NULL, I will return zero---not PSD.
\param semi If anything but 's', check for positive definite, not semidefinite. (default 's')

See also \ref apop_matrix_to_positive_semidefinite, which will change the input to something PSD.

This function uses the \ref designated syntax for inputs.

    */
APOP_VAR_HEAD int apop_matrix_is_positive_semidefinite(gsl_matrix *m, char semi){
    gsl_matrix * apop_varad_var(m, NULL);
    apop_assert(m, 0, 1, 'c', "You gave me a NULL matrix. I will take this as not positive semidefinite.");
    char apop_varad_var(semi, 's');
    return apop_matrix_is_positive_semidefinite_base(m, semi);
APOP_VAR_ENDHEAD
    for (int i=1; i<= m->size1; i++){
        gsl_matrix mv =gsl_matrix_submatrix (m, 0, 0, i, i).matrix;
        double det = apop_matrix_determinant(&mv);
        if ((semi == 'd' && det <0) || det <=0)
            return 0;
    }
    return 1;
}

/**  First, this function passes tests, but is under development.
  
    It takes in a matrix and converts it to the `closest' positive
    semidefinite matrix.

    \param m On input, any matrix; on output, a positive semidefinite matrix.
    \return the distance between the original and new matrices.

    \li See also the test function \ref apop_matrix_is_positive_semidefinite.
    \li This function can be used as (the core of) a model constraint.

   Adapted from the R Matrix package's nearPD, which is 
   Copyright (2007) Jens Oehlschlägel [and is GPL].
 */
double apop_matrix_to_positive_semidefinite(gsl_matrix *m){
    if (apop_matrix_is_positive_semidefinite(m)) 
        return 0; 
    double diffsize=0, dsize;
    apop_data *qdq; 
    gsl_matrix *d = apop_matrix_copy(m);
    double orig_diag_size = diagonal_size(d);
    int size = d->size1;
    gsl_vector *diag = gsl_vector_alloc(size);
    diagonal_copy(diag, d, 'o');
    double origsize = biggest_elmt(d);
    do {
        //get eigenvals
        apop_data *eigenvecs = apop_data_alloc(0, size, size);
        gsl_vector *eigenvals = gsl_vector_calloc(size);
        gsl_matrix *junk_copy = apop_matrix_copy(d);
        find_eigens(&junk_copy, eigenvals, eigenvecs->matrix);
        
        //prune positive only
        int j=0;
        int plussize = eigenvecs->matrix->size1;
        int *mask = calloc(eigenvals->size , sizeof(int));
        for (int i=0; i< eigenvals->size; i++)
            plussize -= 
            mask[i] = (gsl_vector_get(eigenvals, i) <= 0);
        
        //construct Q = pruned eigenvals
        apop_data_rm_columns(eigenvecs, mask);
        if (!eigenvecs->matrix) break;
        
        //construct D = positive eigen diagonal
        apop_data *eigendiag = apop_data_calloc(0, plussize, plussize);
        for (int i=0; i< eigenvals->size; i++)
            if (!mask[i]) {
                apop_data_set(eigendiag, j, j, eigenvals->data[i]);
                j++;
            }

        // Our candidate is QDQ', symmetrized, with the old diagonal subbed in.
        apop_data *qd = apop_dot(eigenvecs, eigendiag);
        qdq = apop_dot(qd, eigenvecs, .form2='t');
        for (int i=0; i< qdq->matrix->size1; i++)
            for (int j=i+1; j< qdq->matrix->size1; j++){
                double avg = (apop_data_get(qdq, i, j) +apop_data_get(qdq, j, i)) /2.;
                apop_data_set(qdq, i, j, avg);
                apop_data_set(qdq, j, i, avg);
            }
        diagonal_copy(diag, qdq->matrix, 'i');
        
        // Evaluate progress, clean up.
        dsize = biggest_elmt(d);
        gsl_matrix_sub(d, qdq->matrix);
        diffsize = biggest_elmt(d);
        apop_data_free(qd); gsl_matrix_free(d);
        apop_data_free(eigendiag); free(mask);
        apop_data_free(eigenvecs); gsl_vector_free(eigenvals);
        d = qdq->matrix;
        qdq->matrix=NULL; apop_data_free(qdq);
    } while (diffsize/dsize > 1e-3);

    apop_data *eigenvecs = apop_data_alloc(0, size, size);
    gsl_vector *eigenvals = gsl_vector_calloc(size);
    gsl_matrix *junk_copy = apop_matrix_copy(d);
    find_eigens(&junk_copy, eigenvals, eigenvecs->matrix);
    //make eigenvalues more positive
    double score =0;
    for (int i=0; i< eigenvals->size; i++){
        double v = gsl_vector_get(eigenvals, i);
        if (v < 1e-1){
            gsl_vector_set(eigenvals, i, 1e-1);
            score += 1e-1 - v;
        }
    }
    if (!score){
        apop_data_free(eigenvecs); gsl_vector_free(eigenvals);
        return 0;
    }
    apop_data *eigendiag = apop_data_calloc(0, size, size);
    diagonal_copy(eigenvals, eigendiag->matrix, 'i');
    double new_diag_size = diagonal_size(eigendiag->matrix);
    gsl_matrix_scale(eigendiag->matrix, orig_diag_size/new_diag_size);
    apop_data *qd = apop_dot(eigenvecs, eigendiag);
    qdq = apop_dot(qd, eigenvecs, .form2='t');

    gsl_matrix_memcpy(m, qdq->matrix);
    assert(apop_matrix_is_positive_semidefinite(m));
    apop_data_free(qdq); gsl_vector_free(diag);
    return diffsize/origsize;
}

/*  This is junk. Please ignore it for now. Thanks.
  */
static double pos_def(apop_data *data, apop_model *candidate){
    double penalty = fabs(candidate->parameters->vector->data[0] - (data->matrix->size1 - data->matrix->size2));
    candidate->parameters->vector->data[0] = data->matrix->size1 - data->matrix->size2;
    return penalty + apop_matrix_to_positive_semidefinite(candidate->parameters->matrix);
}

typedef struct{
    double paramdet;
    gsl_matrix *wparams;
    int len;
} wishartstruct_t;

double one_wishart_row(gsl_vector *in, void *ws_in){
    wishartstruct_t *ws = ws_in;
    gsl_matrix *inv;
    gsl_matrix *inv_dot_params = gsl_matrix_alloc(ws->len, ws->len);
    apop_data *square= apop_data_alloc(0, ws->len, ws->len);
    apop_data_unpack(in, square);
    double datadet = apop_det_and_inv(square->matrix, &inv, 1, 1);
    double out = log(datadet) * ((df - len -1.)/2.);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, inv, ws->wparams, 0, inv_dot_params);   
    assert(datadet);
    gsl_vector_view diag = gsl_matrix_diagonal(inv_dot_params);
    double trace = apop_sum(&diag.vector);
    out += -0.5 * trace;
    out -= log(2) * len*df/2.;
    out -= log(ws->paramdet) * df/2.;
    out -= apop_multivariate_lngamma(df/2, len);
    gsl_matrix_free(inv);
    apop_data_free(square);
    assert(isfinite(out));
    return out;
}

static double wishart_ll(apop_data *in, apop_model *m){
    Nullcheck(in); Nullcheck_m(m); Nullcheck_p(m);
    df = m->parameters->vector->data[0];
    wishartstruct_t ws = {
            .paramdet = apop_matrix_determinant(m->parameters->matrix),
            .wparams = m->parameters->matrix,
            .len = sqrt(in->matrix->size2)
        };
    if (ws.paramdet < 1e-3) return GSL_NEGINF;
    double ll =  apop_map_sum(in, .fn_vp = one_wishart_row, .param=&ws, .part='r');
    //return apop_matrix_map_sum(in->matrix, one_wishart_row);
    printf("------------%g\n", ll);
    return ll;
}

static void apop_wishart_draw(double *out, gsl_rng *r, apop_model *m){
    /*
Translated from the Fortran by BK. Fortran comments:

C          SUBROUTINE DWSHRT(D, N, NP, NNP, SB, SA)
C
C       ALGORITHM AS 53  APPL. STATIST. (1972) VOL.21, NO.3
C
C     Wishart variate generator.  On output, SA is an upper-triangular
C     matrix of size NP * NP [...]
C     whose elements have a Wishart(N, SIGMA) distribution.
*/
    Nullcheck_mv(m); Nullcheck_pv(m);
    int DF, np = m->parameters->matrix->size1;
    int n = m->parameters->vector->data[0];
    if (!m->more) { 
        gsl_matrix *ccc = apop_matrix_copy(m->parameters->matrix);
        gsl_linalg_cholesky_decomp(ccc);
        for (int i=0; i < ccc->size1; i++) //zero out the upper diagonal
            for (int j=i+1; j < ccc->size2; j++)
                gsl_matrix_set(ccc, i, j, 0);
        m->more = apop_matrix_to_data(ccc);
    }
    apop_data *Chol = m->more;
    m->more_size = sizeof(apop_data);
    apop_data *rmatrix = apop_data_calloc(0, np, np);
    static apop_model *std_normal = NULL;
    if (!std_normal) std_normal = apop_model_set_parameters(apop_normal, 0.0, 1.0);

//C     Load diagonal elements with square root of chi-square variates
    for(int i = 0; i< np; i++){
        DF = n - i + 2.;
        apop_data_set(rmatrix, i, i, gsl_ran_chisq(r, DF));
    }
    
    double ndraw;
    for(int i = 0; i< np; i++) //off-diagonal triangles: Normals.
          for(int j = 0; j< i; j++){
            apop_draw(&ndraw, r, std_normal);
            assert (!gsl_isnan(ndraw));
            apop_data_set(rmatrix, i, j, ndraw);
            apop_data_set(rmatrix, j, i, ndraw);
          }

    //Now find C * rand * rand' * C'
    apop_data *cr = apop_dot(Chol, rmatrix);
    apop_data *crr = apop_dot(cr, rmatrix, .form2='t');
    apop_data *crrc = apop_dot(crr, Chol, .form2='t');
    memmove(out, crrc->matrix->data, sizeof(double)*np*np);
    assert(!gsl_isnan(*out));
    apop_data_free(rmatrix); apop_data_free(cr);
    apop_data_free(crrc);    apop_data_free(crr);
}

static void wishart_prep(apop_data *d, apop_model *m){
     m->parameters = apop_data_alloc(1,sqrt(d->matrix->size2),sqrt(d->matrix->size2));
 }

apop_model apop_wishart  = {"Wishart distribution", 1, -1, -1, .dsize=-1, .draw = apop_wishart_draw,
         .log_likelihood = wishart_ll, .constraint = pos_def, .prep=wishart_prep};

apop_model apop_t_distribution  = {"t distribution", 1, 0, 0, .estimate = apop_t_chi_estimate, 
         .log_likelihood = apop_tdist_llike, .draw=apop_t_dist_draw };

apop_model apop_f_distribution  = {"F distribution", 2, 0, 0, .estimate = apop_fdist_estimate, 
        .log_likelihood = apop_fdist_llike, .draw=apop_f_dist_draw };

apop_model apop_chi_squared  = {"Chi squared distribution", 1, 0, 0, .estimate = apop_t_chi_estimate,  
        .log_likelihood = apop_chisq_llike, .draw=apop_chisq_dist_draw };
