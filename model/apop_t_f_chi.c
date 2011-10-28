/* \file apop_t_f_chi.c	t, F, chi squared, and Wishart distributions.
Copyright (c) 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#include "apop_internal.h"

apop_model* apop_t_estimate(apop_data *d, apop_model *m){
    Apop_assert(d, "No data with which to count df. (the default estimation method)");
    Get_vmsizes(d); //vsize, msize1, msize2, tsize
    apop_model *out = apop_model_copy(*m);
    double vmu = vsize ? apop_mean(d->vector) : 0;
    double mmu = msize1 ? apop_matrix_mean(d->matrix) : 0;
    double v_sum_sq = vsize ? apop_var(d->vector)*(vsize-1) : 0;
    double m_sum_sq = msize1 ? apop_matrix_var_m(d->matrix, mmu)*(msize1*msize2-1) : 0;
    apop_name_add(out->parameters->names, "mean", 'r');
    apop_name_add(out->parameters->names, "standard deviation",  'r');
    apop_name_add(out->parameters->names, "df", 'r');
    apop_data_set(out->parameters, 0, -1, (vmu *vsize + mmu * msize1*msize2)/tsize);
    apop_data_set(out->parameters, 1, -1, sqrt((v_sum_sq*vsize + m_sum_sq * msize1*msize2)/(tsize-1))); 
    apop_data_set(out->parameters, 2, -1, tsize-1);
    apop_data_add_named_elmt(out->info, "log likelihood", out->log_likelihood(d, out));
    return out;
}

apop_model* apop_chi_estimate(apop_data *d, apop_model *m){
    Apop_assert(d, "No data with which to count df. (the default estimation method)");
    Get_vmsizes(d); //vsize, msize1, msize2
    apop_model *out = apop_model_copy(*m);
    apop_data_add_named_elmt(out->parameters, "df", tsize - 1);
    apop_data_add_named_elmt(out->info, "log likelihood", out->log_likelihood(d, out));
    return out;
}

apop_model* apop_fdist_estimate(apop_data *d, apop_model *m){
    Apop_assert(d, "No data with which to count df. (the default estimation method)");
    apop_name_add(m->parameters->names, "df",  'r');
    apop_name_add(m->parameters->names, "df2",  'r');
    apop_data_set(m->parameters, 0, -1, d->vector->size -1);
    apop_data_set(m->parameters, 1, -1, d->matrix->size1 * d->matrix->size2 -1);
    apop_data_add_named_elmt(m->info, "log likelihood", apop_f_distribution.log_likelihood(d, m));
    return m;
}

static double one_f(double in, void *df_in){ 
    double *df = df_in; 
    return log(gsl_ran_fdist_pdf(in, df[0], df[1])); 
}

static double one_t(double in, void *params){ 
    double mu = ((double*)params)[0];
    double sigma = ((double*)params)[1];
    double df = ((double*)params)[2];
    return log(gsl_ran_tdist_pdf((in-mu)/(sigma/sqrt(df)), df)); 
}
static double one_chisq(double in, void *df){ return log(gsl_ran_chisq_pdf(in, *(double*)df)); }

double apop_tdist_llike(apop_data *d, apop_model *m){ 
    Nullcheck_mpd(d, m);
    double *params = m->parameters->vector->data;
    return apop_map_sum(d, .fn_dp=one_t, .param=&params);
}

double apop_chisq_llike(apop_data *d, apop_model *m){ 
    Nullcheck_mpd(d, m);
    double df = m->parameters->vector->data[0];
    return apop_map_sum(d, .fn_dp=one_chisq, .param =&df);
}

double apop_fdist_llike(apop_data *d, apop_model *m){ 
    Nullcheck_mpd(d, m);
    double df[2];
    df[0] = m->parameters->vector->data[0];
    df[1] = m->parameters->vector->data[1];
    return apop_map_sum(d, .fn_dp=one_f, .param =df);
}

void apop_t_dist_draw(double *out, gsl_rng *r, apop_model *m){ 
    Nullcheck_mp(m);
    double mu = m->parameters->vector->data[0];
    double sigma = m->parameters->vector->data[1];
    double df = m->parameters->vector->data[2];
    *out = gsl_ran_tdist (r, df)*(sigma/sqrt(df))+mu;
}

double apop_t_dist_cdf(apop_data *in, apop_model *m){
    Nullcheck_mp(m);
    double val = in->vector ? apop_data_get(in, 0, -1) : apop_data_get(in, 0, 0);
    double mu = m->parameters->vector->data[0];
    double sigma = m->parameters->vector->data[1];
    double df = m->parameters->vector->data[2];
    return gsl_cdf_tdist_P ((val-mu)/sigma, df);
}

void apop_f_dist_draw(double *out, gsl_rng *r, apop_model *m){
    Nullcheck_mp(m);
    *out = gsl_ran_fdist (r, m->parameters->vector->data[0], m->parameters->vector->data[1]);
}

void apop_chisq_dist_draw(double *out, gsl_rng *r, apop_model *m){
    Nullcheck_mp(m);
    *out = gsl_ran_chisq (r, m->parameters->vector->data[0]);
}

double apop_t_dist_constraint(apop_data *beta, apop_model *m){
    Staticdef(apop_data *, d_constr, apop_data_fill(apop_data_alloc(2,2,3),
                             0, 0, 1, 0,  //0 < sigma
                            .9, 0, 0, 1)); //.9 < df
    return apop_linear_constraint(m->parameters->vector, d_constr);
}


static double pos_def(apop_data *data, apop_model *candidate){
    return apop_matrix_to_positive_semidefinite(candidate->parameters->matrix);
}

typedef struct{
    double paramdet, df;
    gsl_matrix *wparams;
    int len;
} wishartstruct_t;

static double one_wishart_row(gsl_vector *in, void *ws_in){
    wishartstruct_t *ws = ws_in;
    gsl_matrix *inv;
    gsl_matrix *inv_dot_params = gsl_matrix_alloc(ws->len, ws->len);
    apop_data *square= apop_data_alloc(0, ws->len, ws->len);
    apop_data_unpack(in, square);
    double datadet = apop_det_and_inv(square->matrix, &inv, 1, 1);
    double out = log(datadet) * ((ws->df - ws->len -1.)/2.);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, inv, ws->wparams, 0, inv_dot_params);   
    assert(datadet);
    gsl_vector_view diag = gsl_matrix_diagonal(inv_dot_params);
    double trace = apop_sum(&diag.vector);
    out += -0.5 * trace;
    out -= log(2) * ws->len* ws->df/2.;
    out -= log(ws->paramdet) * ws->df/2.;
    out -= apop_multivariate_lngamma(ws->df/2., ws->len);
    gsl_matrix_free(inv);
    apop_data_free(square);
    assert(isfinite(out));
    return out;
}

static double wishart_ll(apop_data *in, apop_model *m){
    Nullcheck_mpd(in, m);
    wishartstruct_t ws = {
            .paramdet = apop_matrix_determinant(m->parameters->matrix),
            .wparams = m->parameters->matrix,
            .len = sqrt(in->matrix->size2),
            .df = m->parameters->vector->data[0]
        };
    if (ws.paramdet < 1e-3) return GSL_NEGINF;
    double ll =  apop_map_sum(in, .fn_vp = one_wishart_row, .param=&ws, .part='r');
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
    Nullcheck_mp(m);
    int DF, np = m->parameters->matrix->size1;
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
        DF = n - i + 2.;
        apop_data_set(rmatrix, i, i, gsl_ran_chisq(r, DF));
    }
    
    for(int i = 0; i< np; i++) //off-diagonal triangles: Normals.
          for(int j = 0; j< i; j++){
            double ndraw;
            apop_draw(&ndraw, r, std_normal);
            assert (!gsl_isnan(ndraw));
            apop_data_set(rmatrix, i, j, ndraw);
//            apop_draw(&ndraw, r, std_normal);
//            apop_data_set(rmatrix, j, i, ndraw);
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

double wishart_constraint(apop_data *d, apop_model *m){
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

apop_model *wishart_estimate(apop_data *d, apop_model *m){
    Nullcheck_m(m);
    apop_data_set(m->parameters, 0, -1, d->matrix->size1);
    apop_data *summ=apop_data_summarize(d);
    Apop_col_t(summ, "mean", means);
    gsl_vector *t = m->parameters->vector; //mask this while unpacking
    m->parameters->vector=NULL;
    apop_data_unpack(means, m->parameters);
    gsl_matrix_scale(m->parameters->matrix, 1./t->data[0]);
    m->parameters->vector=t;
    apop_data_free(summ);
    return m;
}

/*\amodel apop_wishart The Wishart distribution, which is currently somewhat untested. 

Here's the likelihood function. \f$p\f$ is the dimension of the data and covariance
matrix, \f$n\f$ is the degrees of freedom, \f$\mathbf{V}\f$ is the \f$p\times p\f$
matrix of Wishart parameters, and \f${\mathbf{W}}\f$ is the \f$p\times p\f$ matrix whose
likelihood is being evaluated.  \f$\Gamma_p(\cdot)\f$ is the \ref apop_multivariate_gamma
"multivariate gamma function".

\f$
P(\mathbf{W}) = \frac{\left|\mathbf{W}\right|^\frac{n-p-1}{2}}
                         {2^\frac{np}{2}\left|{\mathbf V}\right|^\frac{n}{2}\Gamma_p(\frac{n}{2})} \exp\left(-\frac{1}{2}{\rm Tr}({\mathbf V}^{-1}\mathbf{W})\right)\f$

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
apop_model apop_wishart  = {"Wishart distribution", 1, -1, -1, .dsize=-1, .estimate=wishart_estimate, .draw = apop_wishart_draw,
         .log_likelihood = wishart_ll, .constraint = pos_def, .prep=wishart_prep, .constraint=wishart_constraint};

/*\amodel apop_t_distribution The t distribution, primarily for descriptive purposes.

If you want to test a hypothesis, you probably don't need this, and should instead use \ref apop_test.  See notes in \ref tfchi.  

\adoc    Input_format     Unordered list of scalars in the matrix and/or vector.     
\adoc    Parameter_format  vector->data[0] = mu<br>
                            vector->data[1] = sigma<br>
                            vector->data[2] = df 
\adoc    Estimate_results  I'll just count elements and set \f$df = n-1\f$. If you set the \c estimate method to \c NULL, via MLE.
\adoc    settings   \ref apop_mle_settings, \ref apop_parts_wanted_settings   
*/

apop_model apop_t_distribution  = {"t distribution", 3, 0, 0, .dsize=1, .estimate = apop_t_estimate, 
         .log_likelihood = apop_tdist_llike, .draw=apop_t_dist_draw, .cdf=apop_t_dist_cdf,
         .constraint=apop_t_dist_constraint };

/*\amodel apop_f_distribution The F distribution, for descriptive purposes.

 If you want to test a hypothesis, you probably don't need this, and should instead use \ref apop_test.  See notes in \ref tfchi.  

\adoc    Input_format     Unordered list of scalars in the matrix and/or vector.     
\adoc    Parameter_format  Zeroth and first elements of the vector are the \f$df\f$s. 
\adoc    Estimate_results  Count elements and set \f$df=\f$ vector count minus one, and 
                          \f$df2=\f$ matrix count minus one. If you set the \c estimate method to \c NULL, via MLE.    
\adoc    settings   \ref apop_mle_settings    
*/
apop_model apop_f_distribution  = {"F distribution", 2, 0, 0, .dsize=1, .estimate = apop_fdist_estimate, 
        .log_likelihood = apop_fdist_llike, .draw=apop_f_dist_draw };

/*\amodel apop_chi_squared The \f$\chi^2\f$ distribution, for descriptive purposes.

 If you want to test a hypothesis, you probably don't need this, and should instead use \ref apop_test.  See notes in \ref tfchi.  

\adoc    Input_format     Unordered list of scalars in the matrix and/or vector.     
\adoc    Parameter_format  Zeroth element of the vector is the \f$df\f$. 
\adoc    Estimate_results  If you do not set an \ref apop_mle_settings group beforehand, I'll just count elements and
                          set \f$df = n-1\f$. Else, via MLE.    
\adoc    settings   \ref apop_mle_settings    
*/

apop_model apop_chi_squared  = {"Chi squared distribution", 1, 0, 0, .dsize=1, .estimate = apop_chi_estimate,  
        .log_likelihood = apop_chisq_llike, .draw=apop_chisq_dist_draw };
