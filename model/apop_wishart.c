/* apop_wishart.c: the Wishart distribution, for modeling purposes.
Copyright (c) 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#include "apop_internal.h"

#if 0

static long double pos_def(apop_data *data, apop_model *candidate){
    return apop_matrix_to_positive_semidefinite(candidate->parameters->matrix);
}

typedef struct{
    double df;
    gsl_matrix *paraminv;
    int len;
} wishartstruct_t;

static double one_wishart_row(gsl_vector *in, void *ws_in){
    wishartstruct_t *ws = ws_in;
    gsl_matrix *invparams_dot_data = gsl_matrix_alloc(ws->len, ws->len);
    apop_data *square= apop_data_alloc(ws->len, ws->len);
    apop_data_unpack(in, square);
    double datadet = apop_matrix_determinant(square->matrix);
    assert(datadet);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, ws->paraminv, square->matrix, 0, invparams_dot_data);   
    gsl_vector_view diag = gsl_matrix_diagonal(invparams_dot_data);
    double trace = apop_sum(&diag.vector);
    gsl_matrix_free(invparams_dot_data);
    apop_data_free(square);
    double out= log(datadet) * (ws->df - ws->len -1.)/2. - trace*ws->df/2.;
    assert(isfinite(out));
    return out;
}

static long double wishart_ll(apop_data *in, apop_model *m){
    Nullcheck_mpd(in, m, GSL_NAN);
    wishartstruct_t ws = {
            .paraminv = apop_matrix_inverse(m->parameters->matrix),
            .len = sqrt(in->matrix->size2),
            .df = m->parameters->vector->data[0]
        };
    double paramdet = apop_matrix_determinant(m->parameters->matrix);
    if (paramdet < 1e-3) return GSL_NEGINF;
    double ll =  apop_map_sum(in, .fn_vp = one_wishart_row, .param=&ws, .part='r');
    double k = log(ws.df)*ws.df/2.;
    k -= M_LN2 * ws.len* ws.df/2.;
    k -= log(paramdet) * ws.df/2.;
    k -= apop_multivariate_lngamma(ws.df/2., ws.len);
    return ll + k*in->matrix->size1;
}

static int apop_wishart_draw(double *out, gsl_rng *r, apop_model *m){
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
    Nullcheck_mp(m, );
    int np = m->parameters->matrix->size1;
    int n = m->parameters->vector->data[0];
    if (!m->more) { 
        gsl_matrix *ccc = apop_matrix_copy(m->parameters->matrix);
        gsl_linalg_cholesky_decomp(ccc);
        for (int i=0; i < ccc->size1; i++) //zero out the upper diagonal
            for (int j=i+1; j < ccc->size2; j++)
                gsl_matrix_set(ccc, i, j, 0);
        m->more = apop_matrix_to_data(ccc);
        m->more_size = sizeof(apop_data);
    }
    apop_data *Chol = m->more;
    apop_data *rmatrix = apop_data_calloc(np, np);
    Staticdef(apop_model *, std_normal, apop_model_set_parameters(apop_normal, 0, 1));

//C     Load diagonal elements with square root of chi-square variates
    for(int i = 0; i< np; i++){
        int DF = n - i;
        apop_data_set(rmatrix, i, i, sqrt(gsl_ran_chisq(r, DF)));
    }
    
    for(int i = 1; i< np; i++) //off-diagonal triangles: Normals.
          for(int j = 0; j< i; j++){
            double ndraw;
            apop_draw(&ndraw, r, std_normal);
            assert (!gsl_isnan(ndraw));
            apop_data_set(rmatrix, i, j, ndraw);
          }
    //Now find C * rand * rand' * C'
    apop_data *cr = apop_dot(Chol, rmatrix);
    apop_data *crr = apop_dot(cr, rmatrix, .form2='t');
    apop_data *crrc = apop_dot(crr, Chol, .form2='t');
    memmove(out, crrc->matrix->data, sizeof(double)*np*np);
    apop_data_free(rmatrix); apop_data_free(cr);
    apop_data_free(crrc);    apop_data_free(crr);
    return 0;
}

static long double wishart_constraint(apop_data *d, apop_model *m){
    double out= apop_matrix_to_positive_semidefinite(m->parameters->matrix);
    double df_minus_dim = m->parameters->vector->data[0] - (m->parameters->matrix->size1-2)-1e-4;
    if (df_minus_dim <= 0){
        out += df_minus_dim;
        m->parameters->vector->data[0] = m->parameters->matrix->size1+2+1e-4;
    }
    return out;
}

static void wishart_prep(apop_data *d, apop_model *m){
     m->parameters = apop_data_alloc(1,sqrt(d->matrix->size2),sqrt(d->matrix->size2));
}

static long double fixed_wishart_ll(apop_data *in, apop_model *m){
    //Let the mean of the input covariances be CM.
    //We need to estimate the df via MLE.
    //However, the right value of the wishart covariance grid is CM/df.
    //So, for a value of df that we're trying, scale CM appropriately.

    gsl_matrix_scale(m->parameters->matrix, 1./m->parameters->vector->data[0]);
    double out = wishart_ll(in, m);
    gsl_matrix_scale(m->parameters->matrix, m->parameters->vector->data[0]);
    return out;
}

static void wishart_estimate(apop_data *d, apop_model *m){
    Nullcheck_m(m, );
    //apop_data_set(m->parameters, 0, -1, d->matrix->size1);
    //Start with cov matrix via mean of inputs; df=NaN
    apop_data_set(m->parameters, 0, -1, GSL_NAN);
    apop_data *summ=apop_data_summarize(d);
    Apop_col_t(summ, "mean", means);
    gsl_vector *t = m->parameters->vector; //mask this while unpacking
    m->parameters->vector=NULL;
    apop_data_unpack(means, m->parameters);
    m->parameters->vector=t;

    //Estimate a model with fixed cov matrix and blank (NaN) df.
    apop_model *modified_wish = apop_model_copy(m);
    modified_wish->log_likelihood = fixed_wishart_ll;
    apop_model *fixed_wish = apop_model_fix_params(modified_wish);
    apop_model *est_via_fix = apop_estimate(d, fixed_wish);

    //copy df from fixed version to the real thing; clean up.
    t->data[0] = apop_data_get(est_via_fix->parameters, 0, -1);
    gsl_matrix_scale(m->parameters->matrix, 1./t->data[0]);
    apop_data_free(summ);
    apop_model_free(modified_wish);
    apop_model_free(fixed_wish);
}

/*\amodel apop_wishart The Wishart distribution, which is currently somewhat untested. 

Here's the likelihood function. \f$p\f$ is the dimension of the data and covariance
matrix, \f$n\f$ is the degrees of freedom, \f$\mathbf{V}\f$ is the \f$p\times p\f$
matrix of Wishart parameters, and \f${\mathbf{W}}\f$ is the \f$p\times p\f$ matrix whose
likelihood is being evaluated.  \f$\Gamma_p(\cdot)\f$ is the \ref apop_multivariate_gamma
"multivariate gamma function".

\f[
P(\mathbf{W}, \mathbf{V}) = \frac{\left|\mathbf{W}\right|^\frac{n-p-1}{2}}
                         {2^\frac{np}{2}\left|{\mathbf V}\right|^\frac{n}{2}\Gamma_p(\frac{n}{2})} \exp\left(-\frac{1}{2}{\rm Tr}({\mathbf V}^{-1}\mathbf{W})\right)\f]

See also notes in \ref tfchi.

\adoc    Input_format     Each row of the input matrix is a single square matrix,
                      flattened; use \ref apop_data_pack to convert your
                      sequence of matrices into rows.     
\adoc    Parameter_format  \f$N\f$ (the degrees of freedom) is the zeroth element of the vector. The matrix holds the matrix of parameters.
\adoc    Estimate_results  Via MLE.    
\adoc    Prep_routine   Just allocates the parameters based on the size of the input data.       
\adoc    RNG  You can use this to generate random covariance matrices, should you need them. See example below. 
\adoc    settings   \ref apop_mle_settings, \ref apop_parts_wanted_settings    
\adoc    Examples Making some random draws:

\code
apop_model *m = apop_estimate(yr_data, apop_wishart);
gsl_matrix *rmatrix = gsl_matrix_alloc(10, 10);
gsl_rng *r = apop_rng_alloc(8765);
for (int i=0; i< 1e8; i++){
    apop_draw(rmatrix->data, r, m);
    do_math_with_matrix(rmatrix);
}
\endcode */
apop_model *apop_wishart  = &(apop_model){"Wishart distribution", 1, -1, -1, .dsize=-1, .estimate=wishart_estimate, .draw = apop_wishart_draw,
         .log_likelihood = wishart_ll, .constraint = pos_def, .prep=wishart_prep, .constraint=wishart_constraint};
#endif
