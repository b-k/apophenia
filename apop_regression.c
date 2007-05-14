/** \file apop_regression.c	Generally, if it assumes something is  Normally distributed, it's here.\n
Copyright (c) 2006 by Ben Klemens. Licensed under the GNU GPL v2.
 \author Ben Klemens
 */

/** \defgroup regression  OLS/GLS: The linear projection methods */
/** \defgroup ttest  T-tests: comparing two vectors */
/** \defgroup asst_tests  Various means of hypothesis testing.

 See also the goodness of fit tests in \ref histograms.
 */

#include "db.h"     //just for apop_opts
#include "types.h"
#include "stats.h"
#include "output.h"
#include "regression.h"
#include "conversions.h"
#include "model/model.h"
#include "linear_algebra.h"
#include <search.h> //lsearch
#include <stdlib.h> //bsearch
#include <assert.h> 
#include <gsl/gsl_blas.h>

/** Allocate an \c apop_OLS_params structure. 

 \param data the data
 \param model   The model, like \c apop_OLS or \c apop_WLS.
 \param method_params   If you are using a nonstandard method to estimate the model, put the params there.
 \return an \c apop_OLS_params 
 */
apop_OLS_params * apop_OLS_params_alloc(apop_data *data, apop_model model){
  apop_OLS_params *out  = malloc(sizeof(*out));
    out->destroy_data       =  0;
    out->want_cov           =  1;
    out->want_expected_value=  1;
    out->model                 = apop_model_copy(model);
    apop_model_clear(data, out->model);
    out->model->model_params   = out;
    out->model->model_params_size = sizeof(*out);
    return out;
}

/** GSL gives p-values for a one-tailed test; convert it to two, assuming a
 symmetric distribution.
 This function is silly and needs to go.
*/
  
double apop_two_tailify(double in){
	return	fabs(1 - (1 - in)*2);
}

static apop_data * produce_t_test_output(int df, double stat, double diff){
  apop_data *out    = apop_data_alloc(0,7,-1);
  double    pval, qval, two_tail;
  if(!gsl_isnan(stat)){
        pval    = gsl_cdf_tdist_P(stat, df);
        qval    = gsl_cdf_tdist_Q(stat, df);
        two_tail= 2*GSL_MIN(pval,qval);
  } else {
        pval    = GSL_NAN;
        qval    = GSL_NAN;
        two_tail= GSL_NAN;
}
    apop_data_add_named_elmt(out, "mean left - right", diff);
    apop_data_add_named_elmt(out, "t statistic", stat);
    apop_data_add_named_elmt(out, "df", df);
    apop_data_add_named_elmt(out, "p value, 1 tail", GSL_MIN(pval,qval));
    apop_data_add_named_elmt(out, "confidence, 1 tail", 1 - GSL_MIN(pval,qval));
    apop_data_add_named_elmt(out, "p value, 2 tail", two_tail);
    apop_data_add_named_elmt(out, "confidence, 2 tail", 1 - two_tail);
    return out;
}


/** Answers the question: with what confidence can I say that the means of these two columns of data are different?
<tt>apop_paired_t_test</tt> answers the question: with what confidence can I say that the mean difference between the two columns is zero?

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
int		a_count	= a->size,
		b_count	= b->size;
double		a_avg	= apop_vector_mean(a),
		a_var	= apop_vector_var(a),
		b_avg	= apop_vector_mean(b),
		b_var	= apop_vector_var(b),
		stat	= (a_avg - b_avg)/ sqrt(b_var/(b_count-1) + a_var/(a_count-1));
	if (apop_opts.verbose){
		printf("1st avg: %g; 1st std dev: %g; 1st count: %i.\n", a_avg, sqrt(a_var), a_count);
		printf("2st avg: %g; 2st std dev: %g; 2st count: %i.\n", b_avg, sqrt(b_var), b_count);
		printf("t-statistic: %g.\n", stat);
	}

int         df      = a_count+b_count-2;
    return produce_t_test_output(df, stat, a_avg - b_avg);
}

/** Answers the question: with what confidence can I say that the mean difference between the two columns is zero?

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
gsl_vector	*diff	= gsl_vector_alloc(a->size);
	gsl_vector_memcpy(diff, a);
	gsl_vector_sub(diff, b);
int		count	= a->size; 
double		avg	= apop_vector_mean(diff),
		var	= apop_vector_var(diff),
		stat	= avg/ sqrt(var/(count-1));
	gsl_vector_free(diff);
	if (apop_opts.verbose){
		printf("avg diff: %g; diff std dev: %g; count: %i; t-statistic: %g.\n", avg, sqrt(var), count, stat);
	}
int     df      = count-1;
    return produce_t_test_output(df, stat, avg);
}

static double two_tailed_t_test(double mean, double variance, long int ct){
double  t_stat  = fabs(mean)/sqrt(variance),
        p       = 1 - gsl_cdf_tdist_P(t_stat, ct),
        q       = 1 - gsl_cdf_tdist_Q(-t_stat, ct);
    return  p + q;
}

/** For many, it is a knee-jerk reaction to a parameter estimation to
test whether each individual parameter differs from zero. This function
does that.

\param est  The \ref apop_estimate, which includes pre-calculated
parameter estimates, var-covar matrix, and the original data set.

Returns nothing. At the end of the routine, the est->parameters->matrix
includes a set of t-test values: p value, confidence (=1-pval), t statistic, standard deviation, one-tailed Pval, one-tailed confidence.

*/
void apop_estimate_parameter_t_tests (apop_model *est){
int     i, df;
double  val, var, pval, tstat, rootn, stddev, two_tail;
    if (!est->data)
        return;
    if (!est->covariance || !est->parameters){
        apop_error(1,'c', "%s: You asked me to estimate t statistics, but I'm missing either the covariance matrix or the parameters (probably the cov matrix)\n", __func__);
        return;
    }
    est->parameters->matrix = gsl_matrix_alloc(est->parameters->vector->size, 7);
    apop_name_add(est->parameters->names, "p value", 'c');
    apop_name_add(est->parameters->names, "confidence", 'c');
    apop_name_add(est->parameters->names, "t statistic", 'c');
    apop_name_add(est->parameters->names, "standard deviation", 'c');
    apop_name_add(est->parameters->names, "p value, 1 tail", 'c');
    apop_name_add(est->parameters->names, "confidence, 1 tail", 'c');
    apop_name_add(est->parameters->names, "df", 'c');
    df       = est->data->matrix   ?
                    est->data->matrix->size1:
                    est->data->vector->size;
    df      -= est->parameters->vector->size;
    rootn    = sqrt(df);
    for (i=0; i< est->parameters->vector->size; i++){
        val     = apop_data_get(est->parameters, i, -1);
        var     = apop_data_get(est->covariance, i, i);
        stddev  = sqrt(var);
        tstat   = val/stddev;
        pval    = (df > 0)? gsl_cdf_tdist_Q(tstat, df): GSL_NAN;
        two_tail= (df > 0)? two_tailed_t_test(val, var, df): GSL_NAN;
        apop_data_set_it(est->parameters, i, "df", df);
        apop_data_set_it(est->parameters, i, "t statistic", tstat);
        apop_data_set_it(est->parameters, i, "standard deviation", stddev);
        apop_data_set_it(est->parameters, i, "p value", two_tail);
        apop_data_set_it(est->parameters, i, "confidence", 1-two_tail);
        apop_data_set_it(est->parameters, i, "p value, 1 tail", pval);
        apop_data_set_it(est->parameters, i, "confidence, 1 tail", 1-pval);
    }
}

/** Runs an F-test specified by \c q and \c c. Your best bet is to see
 the chapter on hypothesis testing in the <a href="http://apophenia.sourceforge.net/gsl_stats.pdf">PDF manual</a> (check the index for F-tests). It will tell you that:
 \f[{N-K\over q}
 {({\bf Q}'\hat\beta - {\bf c})' [{\bf Q}' ({\bf X}'{\bf X})^{-1} {\bf Q}]^{-1} ({\bf Q}' \hat\beta - {\bf c})
 \over {\bf u}' {\bf u} } \sim F_{q,N-K},\f]
 and that's what this function is based on.

 At the moment, this copies the data set. Plan accordingly.

 \param est     an \ref apop_model that you have already calculated.
 \param q       The matrix \f${\bf Q}\f$, where each row represents a hypothesis.
 \param c       The vector \f${\bf c}\f$. The PDF manual explains all of this.
 \return The confidence with which we can reject the joint hypothesis.
 \todo There should be a way to get OLS and GLS to store \f$(X'X)^{-1}\f$. In fact, if you did GLS, this is invalid, because you need \f$(X'\Sigma X)^{-1}\f$, and I didn't ask for \f$\Sigma\f$.
 */
apop_data *apop_F_test (apop_model *est, apop_data *contrast){
gsl_matrix      *set        = est->data->matrix;
gsl_matrix      *q          = contrast->matrix;
gsl_vector      *c          = contrast->vector;
gsl_matrix      *data       = gsl_matrix_alloc(set->size1, set->size2);    //potentially huge.
gsl_matrix      *xpx        = gsl_matrix_calloc(set->size2, set->size2);
gsl_matrix      *xpxinv     = gsl_matrix_calloc(set->size2, set->size2);
gsl_vector      *qprimebeta = gsl_vector_calloc(est->parameters->vector->size);
gsl_matrix      *qprimexpxinv       = gsl_matrix_calloc(est->parameters->vector->size, set->size2);
gsl_matrix      *qprimexpxinvq      = gsl_matrix_calloc(est->parameters->vector->size, est->parameters->vector->size);
gsl_matrix      *qprimexpxinvqinv   = gsl_matrix_calloc(est->parameters->vector->size, est->parameters->vector->size);
gsl_vector      *qprimebetaminusc_qprimexpxinvqinv   = gsl_vector_calloc(est->parameters->vector->size);
gsl_vector      error       = gsl_matrix_column(est->expected->matrix, apop_name_find(est->expected->names, "residual", 'c')).vector;
gsl_vector      v;
double          f_stat, variance, pval;
int             q_df,
                data_df     = set->size1 - est->parameters->vector->size;
apop_data       *out        = apop_data_alloc(0,3,-1);
    sprintf(out->names->title, "F test");
    gsl_matrix_memcpy(data, set);
    v   = gsl_matrix_column(data, 0).vector;
    gsl_vector_set_all(&v, 1);
    apop_matrix_normalize(data, 'c', 'm');
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, data, data, 0, xpx);   
    gsl_matrix_free(data);
    if (q != NULL){
        q_df    = q->size1;
	    apop_det_and_inv(xpx, &xpxinv, 0, 1);		
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, q, xpxinv, 0, qprimexpxinv);  
        gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1, qprimexpxinv, q,  0,  qprimexpxinvq);  
	    apop_det_and_inv(qprimexpxinvq, &qprimexpxinvqinv, 0, 1);		
        gsl_blas_dgemv(CblasNoTrans, 1, q, est->parameters->vector, 0, qprimebeta);
    } else {
        q_df    = est->parameters->vector->size;
        gsl_matrix_memcpy(qprimexpxinvqinv, xpx);
        gsl_vector_memcpy(qprimebeta, est->parameters->vector);
    }
    if (c !=NULL)
        gsl_vector_sub(qprimebeta, c);  //else, c=0, so this is a no-op.
    gsl_blas_dgemv(CblasNoTrans, 1, qprimexpxinvqinv, qprimebeta, 0, qprimebetaminusc_qprimexpxinvqinv);
    gsl_blas_ddot(qprimebeta, qprimebetaminusc_qprimexpxinvqinv, &f_stat);

    gsl_blas_ddot(&error, &error, &variance);
    f_stat  *=  data_df / (variance * q_df);
    pval    = (q_df > 0 && data_df > 0) ? gsl_cdf_fdist_P(f_stat, q_df, data_df): GSL_NAN; 
    gsl_matrix_free(xpx);
    gsl_matrix_free(xpxinv);
    gsl_matrix_free(qprimexpxinv);
    gsl_matrix_free(qprimexpxinvq);
    gsl_matrix_free(qprimexpxinvqinv);
    gsl_vector_free(qprimebeta);
    gsl_vector_free(qprimebetaminusc_qprimexpxinvqinv);
    apop_data_add_named_elmt(out, "F statistic", f_stat);
    apop_data_add_named_elmt(out, "p value", pval);
    apop_data_add_named_elmt(out, "confidence", 1- pval);
    return out;
}

/** a synonym for \ref apop_F_test, qv. */
apop_data * apop_f_test (apop_model *est, apop_data *contrast){
return apop_F_test(est, contrast);
}

//shift first col to depvar, rename first col "one".
static void prep_names (apop_model *e){
  int i;
  apop_OLS_params   *p = e->model_params;
	if (e->data->names->colnamect > 0) {		
		//apop_name_add(n, n->colnames[0], 'd');
        apop_name_add(e->expected->names, e->data->names->colnames[0], 'c');
        apop_name_add(e->expected->names, "predicted", 'c');
        apop_name_add(e->expected->names, "residual", 'c');
        if (e->parameters)
            snprintf(e->parameters->names->title, 100, "Regression of %s", e->data->names->colnames[0]);
		sprintf(e->data->names->colnames[0], "1");
        apop_name_add(e->parameters->names, "1", 'r');
        apop_name_add(e->parameters->names, "parameters", 'v');
        for(i=1; i< e->data->names->colnamect; i++)
            apop_name_add(e->parameters->names, e->data->names->colnames[i], 'r');
        if (p->want_cov){
            if (e->data->names){
                apop_name_stack(e->covariance->names, e->data->names, 'c');
                apop_name_cross_stack(e->covariance->names, e->data->names, 'c', 'r');
            }
		    sprintf(e->covariance->names->colnames[0], "1");
		    sprintf(e->covariance->names->rownames[0], "1");
        }
	}
}

void xpxinvxpy(gsl_matrix *data, gsl_vector *y_data, gsl_matrix *xpx, gsl_vector* xpy, apop_model *out){
  apop_OLS_params   *p = out->model_params;
	if (p->want_cov + p->want_expected_value == 0 ){	
		//then don't calculate (X'X)^{-1}
		gsl_linalg_HH_solve (xpx, xpy, out->parameters->vector);
		return;
	} //else:
  gsl_vector 	*error = gsl_vector_alloc(data->size1);
  gsl_vector 	predicted;
  gsl_matrix	*cov;
  double        s_sq;
	cov	= gsl_matrix_alloc(data->size2, data->size2);
	apop_det_and_inv(xpx, &cov, 0, 1);	    //not yet cov, just (X'X)^-1.
	gsl_blas_dgemv(CblasNoTrans, 1, cov, xpy, 0, out->parameters->vector);      // \beta=(X'X)^{-1}X'Y
	gsl_blas_dgemv(CblasNoTrans, 1, data, out->parameters->vector, 0, error);   // X'\beta ==predicted
	gsl_vector_sub(error,y_data);           //X'\beta - Y == error
    gsl_blas_ddot(error, error, &s_sq);   // e'e
    s_sq    /= data->size1 - data->size2;   //\sigma^2 = e'e / df
	gsl_matrix_scale(cov, s_sq);            //cov = \sigma^2 (X'X)^{-1}
	if (p->want_expected_value){
        gsl_matrix_set_col(out->expected->matrix, 0, y_data);
        gsl_matrix_set_col(out->expected->matrix, 2, error);
        predicted   = gsl_matrix_column(out->expected->matrix, 1).vector;
        gsl_vector_set_zero(&predicted);
        gsl_vector_add(&predicted, y_data);
        gsl_vector_sub(&predicted, error);
    }
    gsl_vector_free(error);
    out->covariance->matrix	= cov;
}

/** generalized least squares.

\ingroup regression

The first column is the dependent variable, the remaining columns are the independent variables. NB: \c data is destroyed by this function. If you want to keep it, make a copy beforehand.

\param set
The first column is the dependent variable, and the remaining columns the independent. Is destroyed in the process, so make a copy beforehand if you need.

\param sigma 
A known variance-covariance matrix, of size <tt>(data->size1, data->size1)</tt>. Survives the function intact. The first column refers to the constant unit vector, so it's always zero.

\param ep
Most notable for its <tt>uses</tt> element, 
If NULL, do everything; else, produce those \ref apop_model elements which you specify. You always get the parameters and never get the log likelihood.

\return
A pointer to an \ref apop_model structure with the appropriate elements filled. See the description in \ref apop_OLS .

\todo 
Since the first column and row of the var/covar matrix is always zero, users shouldn't have to make it.
 */
/*
apop_model * apop_estimate_GLS(apop_data *set, gsl_matrix *sigma){
apop_model      *modded_ols;
    modded_ols              = apop_model_copy(apop_GLS);
    modded_ols->parameter_ct= set->matrix->size2;
apop_estimate	*out	= apop_params_alloc(set, *modded_ols, NULL);
gsl_vector 	*y_data		= gsl_vector_alloc(set->matrix->size1);
gsl_matrix 	*temp		= gsl_matrix_calloc(set->matrix->size2, set->matrix->size1);
gsl_vector 	*xsy 		= gsl_vector_calloc(set->matrix->size2);
gsl_matrix 	*xsx 		= gsl_matrix_calloc(set->matrix->size2, set->matrix->size2);
gsl_matrix 	*sigma_inverse;	//= gsl_matrix_alloc(data->size1, data->size1);
gsl_vector_view	v 		= gsl_matrix_column(set->matrix, 0);
    prep_names(out);
	gsl_matrix_get_col(y_data, set->matrix, 0);
	gsl_vector_set_all(&(v.vector), 1);	                                            //affine: first column is ones.
	apop_det_and_inv(sigma, &sigma_inverse, 0, 1);					                //find sigma^{-1}
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, set->matrix, sigma_inverse, 0, temp); 	//temp = X' \sigma^{-1}.
	gsl_matrix_free(sigma_inverse);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, temp, set->matrix, 0, xsx);    		//(X' \sigma^{-1} X)
	gsl_blas_dgemv(CblasNoTrans, 1, temp, y_data, 0, xsy);     			            //(X' \sigma^{-1} y)
	gsl_matrix_free(temp);
	xpxinvxpy(set->matrix, y_data, xsx, xsy, out);
	gsl_vector_free(y_data); gsl_vector_free(xsy);
	return out;
}
*/

/** ordinary least squares.

\ingroup regression

\param inset The first column is the dependent variable, and the remaining columns the independent. Is destroyed in the process, so make a copy beforehand if you need.

\param epin    An \ref apop_model object. The only
thing we look at is the \c destroy_data element. If this is NULL or
\c destroy_data==0, then the entire data set is copied off, and then
mangled. If \c destroy_data==1, then this doesn't copy off the data set,
but destroys it in place.

\return
Will return an \ref apop_model <tt>*</tt>.
<tt>The_result->parameters</tt> will hold the coefficients; the first
coefficient will be the coefficient on the constant term, and the
remaining will correspond to the independent variables. It will therefore
be of size <tt>(data->size2)</tt>. Do not pre-allocate.

If you asked for it, the covariance matrix, confidence levels, and residuals will also be returned. The confidence intervals give the level of certainty with which we can reject the hypothesis that the given coefficient is zero.

\b sample 

First, you will need a file named <tt>data</tt> in comma-separated form. The first column is the dependent variable; the remaining columns are the independent. For example:
\verbatim
Y, X_1, X_2, X_3
2,3,4,5
1,2,9,3
4,7,9,0
2,4,8,16
1,4,2,9
9,8,7,6
\endverbatim

The program:
\code
#include <apop.h>

int main(void){
apop_data       *data;
apop_model   *est;
    apop_text_to_db("data","d",0,1,NULL);
    data = apop_query_to_data("select * from d");
    est  = apop_OLS.estimate(data, NULL);
    apop_params_show(est);
    return 0;
}
\endcode

If you saved this code to <tt>sample.c</tt>, then you can compile it with
\verbatim
gcc sample.c -lapophenia -lgsl -lgslcblas -lsqlite3 -o run_me
\endverbatim

and then run it with <tt>./run_me</tt>. Alternatively, you may prefer to compile the program using a \ref makefile .

Feeling lazy? The program above was good form and demonstrated useful
features, but the code below will do the same thing in four lines:

\code
#include <apop.h>
int main(){
    apop_estimate_show(apop_OLS.estimate(apop_text_to_data("data", 0, 0), NULL));
    return 0; }
\endcode


 */
apop_model * apop_estimate_OLS(apop_data *inset, apop_model *ep){
  apop_data         *set;
  apop_model       *epout;
  gsl_vector        *weights    = NULL;
  int               i;
  apop_OLS_params  *olp;
    if (!ep) {
        olp             = apop_OLS_params_alloc(inset, apop_OLS);
        epout           = olp->model;
    } else {
        olp             = ep->model_params;
        epout           = ep;
    }
    set = olp->destroy_data ? inset : apop_data_copy(inset); 
    
    //prep weights.
    if (olp->destroy_data)
        weights = epout->data->weights;  //may be NULL.
    else
        weights = apop_vector_copy(epout->data->weights); //may be NULL.
    if (weights)
        for (i =0; i< weights->size; i++)
            gsl_vector_set(weights, i, sqrt(gsl_vector_get(weights, i)));

  gsl_vector    *y_data     = gsl_vector_alloc(set->matrix->size1); 
  gsl_vector    *xpy        = gsl_vector_calloc(set->matrix->size2);
  gsl_matrix    *xpx        = gsl_matrix_calloc(set->matrix->size2, set->matrix->size2);
    if (olp->want_expected_value)
        epout->expected   = apop_data_alloc(0, set->matrix->size1, 3);
    if (olp->want_cov)
        epout->covariance = apop_data_alloc(0, set->matrix->size1, set->matrix->size1);
    prep_names(epout);
    APOP_COL(set, 0, firstcol);
    gsl_vector_memcpy(y_data,firstcol);
    gsl_vector_set_all(firstcol, 1);     //affine: first column is ones.
    if (weights){
        gsl_vector_mul(y_data, weights);
        for (i = 0; i < set->matrix->size2; i++){
            APOP_COL(set, i, v);
            gsl_vector_mul(v, weights);
        }
    }

    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, set->matrix, set->matrix, 0, xpx);   //(X'X)
    gsl_blas_dgemv(CblasTrans, 1, set->matrix, y_data, 0, xpy);       //(X'y)
    xpxinvxpy(set->matrix, y_data, xpx, xpy, epout);
    gsl_vector_free(y_data); gsl_vector_free(xpy);

    if (!olp->destroy_data)
        apop_data_free(set);
    apop_estimate_parameter_t_tests(epout);
    epout->status   = 1;
    return epout;
}

/** The partitioned regression.
        Give me two data sets, and I'll take the steps needed to produce
        OLS coefficients.  If you give me intermediate results, I'll
        use them instead of calculating them myself; that way the
        fixed-effects regression is really easy.

        The partitioned regression reduces to a GLS regression where
        the Sigma matrix (usually interpreted as the var-covar matrix)
        for the first set is \f$M_2 = I - X_2'(X_2 X_2') X_2\f$
        and the Sigma for for the second set is \f$M_1 = I - X_1'(X_1
        X_1') X_1\f$.  So, calculate \f$M_1\f$ and \f$M_2\f$, run GLS
        twice, and splice together the results.

        I assume that the Y vector (the dependent variable) is the first
        column of the _second_ matrix. The textbook custom seems to be
        to put auxiliary dummy variables on the left, and I won't 
        go against that here.

\param set1    the first half of the data set
\param set2    the second half of the data set
\param  m1      If you have these, give them to me and I won't waste time calculating them.
\param  m2      If you have these, give them to me and I won't waste time calculating them.

\bug The cross-variances are assumed to be zero, which is wholeheartedly false. It's not too big a deal because nobody ever uses them for anything.
*/
/*
apop_model * apop_partitioned_OLS(apop_data *set1, apop_data *set2, gsl_matrix *m1, gsl_matrix *m2){
apop_inventory  actual_uses    = apop_inventory_filter(uses, apop_GLS.inventory_filter);
    prep_inventory_names(set1->names);
apop_estimate	*out1, *out2,
                *out	= apop_params_alloc(set1->matrix->size1, set1->matrix->size2 + set2->matrix->size2, set1->names, actual_uses);
gsl_matrix      *t1,*t2, *augmented_first_matrix, *zero1, *zero2,
                *dependent  =gsl_matrix_alloc(set1->matrix->size1, 1);
gsl_vector_view d;
int             i;


    if(m1==NULL){
    gsl_matrix 	*xpx1 		= gsl_matrix_calloc(set1->matrix->size1, set1->matrix->size1),
                *xpxinv1, *xxpxinv1, *xxpxinvx1;
        m1        = gsl_matrix_alloc(set1->matrix->size2, set1->matrix->size2);
	    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, set1->matrix, set1->matrix, 0, xpx1);	//(X'X)
	    apop_det_and_inv(xpx1, &xpxinv1, 0, 1);		
        gsl_matrix_free(xpx1);
        xxpxinv1   = gsl_matrix_calloc(set1->matrix->size1, set1->matrix->size1);
	    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, set1->matrix, xpxinv1, 0, xxpxinv1);	//X(X'X)^{-1}
        gsl_matrix_free(xpxinv1);
        xxpxinvx1  = gsl_matrix_calloc(set1->matrix->size1, set1->matrix->size1);
	    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1, xxpxinv1, set1->matrix, 0, xxpxinvx1);	//X(X'X)^{-1}X'
        gsl_matrix_free(xxpxinv1);
        gsl_matrix_set_identity(m1);
        gsl_matrix_sub(m1, xxpxinvx1);   // M = I - X(X'X)^{-1}X'
        gsl_matrix_free(xxpxinvx1);
    }

    if(m2==NULL){
    gsl_matrix 	*xpx2 		= gsl_matrix_calloc(set2->matrix->size1, set2->matrix->size1),
                *xpxinv2, *xxpxinv2, *xxpxinvx2;
        m2        = gsl_matrix_alloc(set2->matrix->size2, set2->matrix->size2);
	    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, set2->matrix, set2->matrix, 0, xpx2);	//(X'X)
	    apop_det_and_inv(xpx2, &xpxinv2, 0, 1);		
        gsl_matrix_free(xpx2);
        xxpxinv2   = gsl_matrix_calloc(set2->matrix->size1, set2->matrix->size1), 
	    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, set2->matrix, xpxinv2, 0, xxpxinv2);	//X(X'X)^{-1}
        gsl_matrix_free(xpxinv2);
        xxpxinvx2  = gsl_matrix_calloc(set2->matrix->size1, set2->matrix->size1);
	    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1, xxpxinv2, set2->matrix, 0, xxpxinvx2);	//X(X'X)^{-1}X'
        gsl_matrix_free(xxpxinv2);
        gsl_matrix_set_identity(m2);
        gsl_matrix_sub(m2, xxpxinvx2);   // M = I - X(X'X)^{-1}X'
        gsl_matrix_free(xxpxinvx2);
    }


        //the first matrix needs a dependent variable column, then we
        //can regress.
    d                       = gsl_matrix_column(set2->matrix, 0);
    gsl_matrix_set_col(dependent, 0, &(d.vector));
    augmented_first_matrix  = apop_matrix_stack(dependent, set1->matrix, 'c');
apop_data *newfirst = malloc(sizeof(apop_data));
    newfirst->matrix  = augmented_first_matrix;
    newfirst->names = set1->names;
    out1    = apop_estimate_GLS(newfirst, &actual_uses, m2);
    out2    = apop_estimate_GLS(set2, &actual_uses, m1);
    for (i=0; i< out1->parameters->vector->size; i++)
        gsl_vector_set(out->parameters->vector, i, gsl_vector_get(out1->parameters->vector, i));
    for (   ; i< out2->parameters->vector->size; i++)
        gsl_vector_set(out->parameters->vector, i, gsl_vector_get(out2->parameters->vector, i - out1->parameters->vector->size));

    //The covariance matrix is a stacking-up of the above matrices. I
    //cheat here: the cross-variances are assumed zero, which is
    //blatantly false.
    zero1           = gsl_matrix_calloc(out1->covariance->size1, out2->covariance->size2);
    t1              = apop_matrix_stack(out1->covariance, zero1, 'c');
    gsl_matrix_free(zero1);
    zero2           = gsl_matrix_calloc(out2->covariance->size1, out1->covariance->size2);
    t2              = apop_matrix_stack(zero2, out2->covariance, 'c');
    gsl_matrix_free(zero1);
    out->covariance = apop_matrix_stack(t1,t2,'r');
    return out;
}
*/
//
//
//Cut and pasted from the GNU std library documentation:
static int compare_doubles (const void *a, const void *b) {
    const double *da = (const double *) a;
    const double *db = (const double *) b;
    return (*da > *db) - (*da < *db);
}

typedef const char * ccp;
static int strcmpwrap(const void *a, const void *b){
    const ccp *aa = a;
    const ccp *bb = b;
    return strcmp(*aa, *bb);}

/** A utility to make a matrix of dummy variables. You give me a single
vector that lists the category number for each item, and I'll return
a gsl_matrix with a single one in each row in the column specified.

\param  in The gsl_vector of categories
\param  keep_first  if zero, return 
    a matrix where each row has a one in the (column specified MINUS
    ONE). That is, the zeroth category is dropped, the first category
    has an entry in column zero, et cetera. If you don't know why this
    is useful, then this is what you need. If you know what you're doing
    and need something special, set this to one and the first category won't be dropped.
param type 'd'==data column (-1==vector), 't'==text column.
\return out
*/
apop_data * apop_data_produce_dummies(apop_data *d, int col, char type, int keep_first){
  size_t      i, index,
              prior_elmt_ctr  = 107,
              elmt_ctr        = 0;
  apop_data   *out; 
  double      *elmts          = malloc(sizeof(double)),
              val;
  char        n[1000], **tval,
              **telmts        = malloc(sizeof(char*));
  gsl_vector  *in             = NULL;
    if (type == 'd'){
        if ((col == -1) && !d->vector)
            apop_error(0, 's', "%s: You asked for the vector element (col==-1) but the data's vector element is NULL.\n", __func__);
        if (col >=0 && col >= d->matrix->size2)
            apop_error(0, 's', "%s: You asked for the matrix element %i but the data's matrix element has only %i columns.\n", 
                    __func__, col, d->matrix->size2);
        APOP_COL(d, col, in_t);
        in  = in_t;
    } else
        if (col >= d->textsize[1])
            apop_error(0, 's', "%s: You asked for the text element %i but the data's text element has only %i elements.\n", 
                    __func__, col, d->textsize[1]);
    //first, create an ordered list of unique elements.
    //Note to the reader: this may or may not be more efficient on the db side.
    int s = type == 'd' ? in->size : d->textsize[0];
    if (type == 'd')
        for (i=0; i< s; i++){
            if (prior_elmt_ctr != elmt_ctr)
                elmts   = realloc(elmts, sizeof(double)*(elmt_ctr+1));
            prior_elmt_ctr  = elmt_ctr;
            val     =  gsl_vector_get(in, i);
            lsearch (&val, elmts, &elmt_ctr, sizeof(double), compare_doubles);
            qsort(elmts, elmt_ctr, sizeof(double), compare_doubles);
        }
    else 
        for (i=0; i< s; i++){
            if (prior_elmt_ctr != elmt_ctr)
                telmts   = realloc(telmts, sizeof(char**)*(elmt_ctr+1));
            prior_elmt_ctr  = elmt_ctr;
            tval    =  &(d->text[i][col]);
            lsearch (tval, telmts, &elmt_ctr, sizeof(char**), strcmpwrap);
            qsort(telmts, elmt_ctr, sizeof(char**), strcmpwrap);
        }

    //Now go through the input vector, and for row i find the posn of the vector's
    //name in the element list created above (j), then change (i,j) in
    //the dummy matrix to one.
    if (keep_first)     out  = apop_data_calloc(0, s, elmt_ctr);
    else                out  = apop_data_calloc(0, s, elmt_ctr-1);
    for (i=0; i< s; i++){
        if (type == 'd'){
            val     = gsl_vector_get(in, i);
            index   = ((size_t)bsearch(&val, elmts, elmt_ctr, sizeof(double), compare_doubles) - (long int)elmts)/sizeof(double);
        } else 
            index   = ((size_t)bsearch(&(d->text[i][col]), telmts, elmt_ctr, sizeof(char**), strcmpwrap) - (size_t)telmts)/sizeof(char**);
        if (keep_first)
            gsl_matrix_set(out->matrix, i, index,1); 
        else if (index>0)
            gsl_matrix_set(out->matrix, i, index-1, 1); 
        //else don't keep first, and index==0; throw it out.
    }
    //Add names:
    i   = (keep_first) ? 0 : 1;
    for (  ; i< elmt_ctr; i++){
        if (type =='d')
            sprintf(n,"dummy %g", elmts[i]);
        else
            sprintf(n, telmts[i]);
        apop_name_add(out->names, n, 'c');
    }
    return out;
}

/** A fixed-effects regression.
    The input is a data matrix for a regression, plus a single vector giving the fixed effect vectors.

    The solution of a fixed-effects regression is via a partitioned
    regression. Given that the data set is divided into columns
    \f$\beta_1\f$ and \f$\beta_2\f$, then the reader may 

\todo finish this documentation. [Was in a rush today.]
*/
apop_model *apop_estimate_fixed_effects_OLS(apop_data *data,  gsl_vector *categories){
apop_data *dummies = apop_data_produce_dummies(apop_vector_to_data(categories),-1, 'd', 0);
    apop_data_stack(data, dummies, 'c');
    return apop_OLS.estimate(dummies, NULL);
}

/** Good ol' \f$R^2\f$.  Let \f$Y\f$ be the dependent variable,
\f$\epsilon\f$ the residual,  \f$n\f$ the number of data points, and \f$k\f$ the number of independent vars (including the constant).

  \f$ SST \equiv \sum (Y_i - \bar Y) ^2 \f$
  \f$ SSE \equiv \sum \epsilon ^2       \f$
  \f$ R^2 \equiv 1 - {SSE\over SST}     \f$
  \f$ R^2_{adj} \equiv R^2 - {(k-1)\over (n-k)}(1-R^2)     \f$

  Internally allocates (and frees) a vector the size of your data set.
\param  in  The estimate. I need residuals to have been calculated, and the first column of in->data needs to be the dependent variable.

\return: a \f$1 \times 5\f$ apop_data table with the following fields:
"R_squared"\\
"R_squared_adj"\\
"SSE"\\
"SST"\\
"SSR"

I need the \c apop_model sent in to have the expected value table. This
is the default, but if you explicitly set want_expected_value = 0,
then turn back and unset it.

\ingroup regression
  */
apop_data *apop_estimate_correlation_coefficient (apop_model *in){
  double          sse, sst, rsq, adjustment;
  size_t          obs     = in->data->matrix->size1;
  size_t          indep_ct= in->data->matrix->size2 - 1;
  gsl_vector      v;  
  apop_data       *out    = apop_data_alloc(0, 5,-1);
  apop_OLS_params *p      = in->model_params;
    if (!p->want_expected_value){
        apop_error(0, 'c', "I need an estimate that used want_expected_value to calculate the correlation coefficient. returning NULL.\n");
        return NULL;
    }
    v   = gsl_matrix_column(in->expected->matrix, 
                apop_name_find(in->expected->names, "residual", 'c')).vector;
    gsl_blas_ddot(&v, &v, &sse);
    v   = gsl_matrix_column(in->expected->matrix, 0).vector; //actual.
    sst = apop_vector_var(&v) * (v.size - 1);
//    gsl_blas_ddot(&v, &v, &sst);
    rsq = 1. - (sse/sst);
    adjustment  = ((indep_ct -1.) /(obs - indep_ct)) * (1.-rsq) ;
    apop_data_add_named_elmt(out, "R_squared", rsq);
    apop_data_add_named_elmt(out, "R_squared_adj", rsq - adjustment);
    apop_data_add_named_elmt(out, "SSE", sse);
    apop_data_add_named_elmt(out, "SST", sst);
    apop_data_add_named_elmt(out, "SSR", sst - sse);
    return out;
}

/** A synonym for \ref apop_estimate_correlation_coefficient, q.v. 
 \ingroup regression
 */
apop_data *apop_estimate_r_squared (apop_model *in){
    return apop_estimate_correlation_coefficient(in);
}
