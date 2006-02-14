/** \file apop_regression.c	Generally, if it assumes something is  Normally distributed, it's here.\n
 Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
 \author Ben Klemens
 */

/** \defgroup regression  OLS/GLS: The linear projection methods */
/** \defgroup ttest  T-tests: comparing two vectors */
/** \defgroup asst_tests  Various means of hypothesis testing.*/

#include "db.h"     //just for apop_opts
#include "types.h"
#include "stats.h"
#include "model/model.h"
#include "output.h"
#include "regression.h"
#include "linear_algebra.h"
#include <search.h> //lsearch
#include <stdlib.h> //bsearch
#include <gsl/gsl_blas.h>


double two_tailify(double in){
//GSL gives me a one-tailed test; convert it to two.
	return	fabs(1 - (1 - in)*2);
}

/** Answers the question: with what confidence can I say that the means of these two columns of data are different?
<tt>apop_paired_t_test</tt> answers the question: with what confidence can I say that the mean difference between the two columns is zero?

\ingroup ttest
\param {a, b} two columns of data
\return the confidence level---if it is close to one, you can reject the null, while <tt>apop_t_test(a, a)</tt> will return zero.
*/
double	apop_t_test(gsl_vector *a, gsl_vector *b){
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
	return two_tailify(gsl_cdf_tdist_P(stat, a_count+b_count-2));
}

/** Answers the question: with what confidence can I say that the mean difference between the two columns is zero?

\ingroup ttest
\param {a, b} two columns of data
\return plus or minus the confidence level---if it is close to one, you can reject the null, while <tt>apop_paired_t_test(a, a)</tt> will return zero. If the confidence level is positive, then mean(a) > mean(b), and the confidence level tells us how firmly we can make that statement for the population. If the confidence level is negative, then mean(b) < mean(a).
*/
double	apop_paired_t_test(gsl_vector *a, gsl_vector *b){
gsl_vector	*diff	= gsl_vector_alloc(a->size);
	gsl_vector_memcpy(diff, a);
	gsl_vector_sub(diff, b);
int		count	= a->size, 
		factor	= 1;
double		avg	= apop_vector_mean(diff),
		var	= apop_vector_var(diff),
		stat	= avg/ sqrt(var/(count-1));
	if (apop_vector_mean(diff) < 0) 
		factor = -1;
	gsl_vector_free(diff);
	if (apop_opts.verbose){
		printf("avg diff: %g; diff std dev: %g; count: %i; t-statistic: %g.\n", avg, sqrt(var), count, stat);
	}
	return factor * two_tailify(gsl_cdf_tdist_P(stat, count-1));
}

/** Runs an F-test specified by \c q and \c c. Your best bet is to see
 the chapter on "Gaussian Tricks" in the <a href="http://apophenia.sourceforge.net/gsl_stats.pdf">PDF manual</a> (check the index for F-tests). It will tell you that:
 \f[{N-K\over q}
 {({\bf Q}'\hat\beta - {\bf c})' [{\bf Q}' ({\bf X}'{\bf X})^{-1} {\bf Q}]^{-1} ({\bf Q}' \hat\beta - {\bf c})
 \over {\bf u}' {\bf u} } \sim F_{q,N-K},\f]
 and that's what this function is based on.
 \param est     an \ref apop_estimate that you have already calculated.
 \param set    your \ref apop_data set. If NULL, use the one included in \c est.
 \param q       The matrix \f${\bf Q}\f$, where each row represents a hypothesis.
 \param c       The vector \f${\bf c}\f$. The PDF manual explains all of this.
 \return The confidence with which we can reject the joint hypothesis.
 \todo There should be a way to get OLS and GLS to store \f$(X'X)^{-1}\f$. In fact, if you did GLS, this is invalid, because you need \f$(X'\Sigma X)^{-1}\f$, and I didn't ask for \f$\Sigma\f$.
 */
double apop_F_test(apop_estimate *est, apop_data *set, gsl_matrix *q, gsl_vector *c){
gsl_matrix      *xpx        = gsl_matrix_calloc(set->data->size2, set->data->size2);
gsl_matrix      *xpxinv     = gsl_matrix_calloc(set->data->size2, set->data->size2);
gsl_vector      *qprimebeta = gsl_vector_calloc(est->parameters->size);
gsl_matrix      *qprimexpxinv       = gsl_matrix_calloc(est->parameters->size, set->data->size2);
gsl_matrix      *qprimexpxinvq      = gsl_matrix_calloc(est->parameters->size, est->parameters->size);
gsl_matrix      *qprimexpxinvqinv   = gsl_matrix_calloc(est->parameters->size, est->parameters->size);
gsl_vector      *qprimebetaminusc_qprimexpxinvqinv   = gsl_vector_calloc(est->parameters->size);
double          f_stat;
double          variance;
int             q_df,
                data_df     = set->data->size1 - est->parameters->size;
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, set->data, set->data, 0, xpx);   
    if (q != NULL){
        q_df    = q->size1;
	    apop_det_and_inv(xpx, &xpxinv, 0, 1);		
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, q, xpxinv, 0, qprimexpxinv);  
        gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1, qprimexpxinv, q,  0,  qprimexpxinvq);  
	    apop_det_and_inv(qprimexpxinvq, &qprimexpxinvqinv, 0, 1);		
        gsl_blas_dgemv(CblasNoTrans, 1, q, est->parameters, 0, qprimebeta);
    } else {
        q_df    = est->parameters->size;
        gsl_matrix_memcpy(qprimexpxinvqinv, xpx);
        gsl_vector_memcpy(qprimebeta, est->parameters);
    }
    if (c !=NULL)
        gsl_vector_sub(qprimebeta, c);  //else, c=0, so this is a no-op.
    gsl_blas_dgemv(CblasNoTrans, 1, qprimexpxinvqinv, qprimebeta, 0, qprimebetaminusc_qprimexpxinvqinv);
    gsl_blas_ddot(qprimebetaminusc_qprimexpxinvqinv, qprimebeta, &f_stat);

    gsl_blas_ddot(est->parameters, est->parameters, &variance);
    f_stat  *=  data_df / (variance * q_df);
    printf("the f statistic: %g\n", f_stat);
    gsl_matrix_free(xpx);
    gsl_matrix_free(xpxinv);
    gsl_matrix_free(qprimexpxinv);
    gsl_matrix_free(qprimexpxinvq);
    gsl_matrix_free(qprimexpxinvqinv);
    gsl_vector_free(qprimebeta);
    gsl_vector_free(qprimebetaminusc_qprimexpxinvqinv);
    return 1 - gsl_cdf_fdist_P(f_stat, q_df, data_df);
}

/** a synonym for \ref apop_F_test, qv. */
double apop_f_test(apop_estimate *est, apop_data *set, gsl_matrix *q, gsl_vector *c){
return apop_F_test(est, set, q, c);
}

//shift first col to depvar, rename first col "one".
static void prep_inventory_names(apop_name *n){
	if (n != NULL && n->colnamect > 0) {		
		apop_name_add(n, n->colnames[0], 'd');
		sprintf(n->colnames[0], "1");
	}
}

void xpxinvxpy(gsl_matrix *data, gsl_vector *y_data, gsl_matrix *xpx, gsl_vector* xpy, apop_estimate *out){
	if (out->uses.covariance + out->uses.confidence + out->uses.residuals == 0 ){	
		//then don't calculate (X'X)^{-1}
		gsl_linalg_HH_solve (xpx, xpy, out->parameters);
		return;
	} //else:
gsl_vector 	*error;
gsl_matrix	*cov;
double		upu;
int		i;
	error	= gsl_vector_alloc(data->size1);
	cov	= gsl_matrix_alloc(data->size2, data->size2);
	apop_det_and_inv(xpx, &cov, 0, 1);		//(X'X)^{-1} (not yet cov)
	gsl_blas_dgemv(CblasNoTrans, 1, cov, xpy, 0, out->parameters);
	gsl_blas_dgemv(CblasNoTrans, 1, data, out->parameters, 0, error);
	if (out->uses.predicted)	
		gsl_vector_memcpy(out->predicted, error);
	gsl_vector_sub(y_data, error);	//until this line, 'error' is the predicted values
	gsl_blas_ddot(error, error, &upu);
	gsl_matrix_scale(cov, upu/data->size2);	//Having multiplied by the variance, it's now it's the covariance.
	if (out->uses.confidence)
		for (i=0; i < data->size2; i++)  // confidence[i] = |1 - (1-N(Mu[i],sigma[i]))*2|
			gsl_vector_set(out->confidence, i,
			   two_tailify(gsl_cdf_gaussian_P(gsl_vector_get(out->parameters, i), gsl_matrix_get(cov, i, i))));
	if (out->uses.residuals == 0) 	gsl_vector_free(error);
	else 				out->residuals	= error;
	if (out->uses.covariance == 0) 	gsl_matrix_free(cov);
	else 				out->covariance	= cov;
}

/** generalized least squares.

\ingroup regression

The first column is the dependent variable, the remaining columns are the independent variables. NB: \c data is destroyed by this function. If you want to keep it, make a copy beforehand.

\param set
The first column is the dependent variable, and the remaining columns the independent. Is destroyed in the process, so make a copy beforehand if you need.

\param sigma 
A known variance-covariance matrix, of size <tt>(data->size1, data->size1)</tt>. Survives the function intact. The first column refers to the constant unit vector, so it's always zero.

\param uses 
If NULL, do everything; else, produce those \ref apop_estimate elements which you specify. You always get the parameters and never get the log likelihood.

\return
A pointer to an \ref apop_estimate structure with the appropriate elements filled. See the description in \ref apop_OLS .

\todo 
Since the first column and row of the var/covar matrix is always zero, users shouldn't have to make it.
 */
apop_estimate * apop_estimate_GLS(apop_data *set, apop_inventory *uses, gsl_matrix *sigma){
    prep_inventory_names(set->names);
apop_model      *modded_ols;
    modded_ols              = apop_model_copy(apop_GLS);
    modded_ols->parameter_ct= set->data->size2;
apop_estimate	*out	= apop_estimate_alloc(set, *modded_ols, uses, NULL);
gsl_vector 	*y_data		= gsl_vector_alloc(set->data->size1);
gsl_matrix 	*temp		= gsl_matrix_calloc(set->data->size2, set->data->size1);
gsl_vector 	*xsy 		= gsl_vector_calloc(set->data->size2);
gsl_matrix 	*xsx 		= gsl_matrix_calloc(set->data->size2, set->data->size2);
gsl_matrix 	*sigma_inverse;	//= gsl_matrix_alloc(data->size1, data->size1);
gsl_vector_view	v 		= gsl_matrix_column(set->data, 0);
	gsl_matrix_get_col(y_data, set->data, 0);
	gsl_vector_set_all(&(v.vector), 1);	                                            //affine: first column is ones.
	apop_det_and_inv(sigma, &sigma_inverse, 0, 1);					                //find sigma^{-1}
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, set->data, sigma_inverse, 0, temp); 	//temp = X' \sigma^{-1}.
	gsl_matrix_free(sigma_inverse);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, temp, set->data, 0, xsx);    		//(X' \sigma^{-1} X)
	gsl_blas_dgemv(CblasNoTrans, 1, temp, y_data, 0, xsy);     			            //(X' \sigma^{-1} y)
	gsl_matrix_free(temp);
	xpxinvxpy(set->data, y_data, xsx, xsy, out);
	gsl_vector_free(y_data); gsl_vector_free(xsy);
	return out;
}

/** ordinary least squares.

\ingroup regression

\param inset The first column is the dependent variable, and the remaining columns the independent. Is destroyed in the process, so make a copy beforehand if you need.

\param uses If <tt>NULL</tt>, then you get everything.  If a pointer to
an \ref apop_inventory , then you get what you ask for. Log likelihood is
not calculated; you always get the parameter estimates.
\param ep    An \ref apop_estimation_params object. The only
thing we look at is the \c destroy_data element. If this is NULL or
\c destroy_data==0, then the entire data set is copied off, and then
mangled. If \c destroy_data==1, then this doesn't copy off the data set,
but destroys it in place.

\return
Will return an \ref apop_estimate <tt>*</tt>.
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
#include <apophenia/headers.h>

int main(void){
apop_data       *data;
apop_estimate   *est;
apop_inventory  invent;
    apop_text_to_db("data","d",0,1,NULL);
    apop_inventory_set(&invent, 0);
    invent.parameters   = 1;
    data = apop_query_to_data("select * from d");
    est  = apop_OLS.estimate(data, &invent, NULL);
    apop_estimate_print(est);
    return 0;
}
\endcode

If you saved this code to <tt>sample.c</tt>, then you can compile it with
\verbatim
gcc sample.c -lapophenia -lgsl -lgslcblas -lsqlite3 -o run_me
\endverbatim

and then run it with <tt>./run_me</tt>. Alternatively, you may prefer to compile the program using a \ref makefile .

Feeling lazy? The program above was good form and demonstrated useful
features, but the code below will do the same thing in five lines:

\code
#include <apophenia/headers.h>
int main(void){
    apop_text_to_db("data","d",NULL);
    apop_estimate_print(apop_OLS.estimate(apop_query_to_data("select * from d"), NULL, NULL));
    return 0; }
\endcode


 */
apop_estimate * apop_estimate_OLS(apop_data *inset, apop_inventory *uses, 
                                            apop_estimation_params *ep){
apop_model      *modded_ols;
apop_data       *set;

    //check whether we get to destroy the data set or need to copy it.
    if ((ep == NULL) || (ep->destroy_data==0))
        apop_data_memcpy(&set, inset); 
    else
        set = inset;

    modded_ols              = apop_model_copy(apop_OLS); 
    modded_ols->parameter_ct= set->data->size2;
    prep_inventory_names(set->names);
apop_estimate	*out		= apop_estimate_alloc(inset, *modded_ols, uses, ep);
gsl_vector      *y_data     = gsl_vector_alloc(set->data->size1); 
gsl_vector      *xpy        = gsl_vector_calloc(set->data->size2);
gsl_matrix      *xpx        = gsl_matrix_calloc(set->data->size2, set->data->size2);
gsl_vector_view v           = gsl_matrix_column(set->data, 0);
    gsl_matrix_get_col(y_data, set->data, 0);
    gsl_vector_set_all(&(v.vector), 1);                 //affine: first column is ones.
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, set->data, set->data, 0, xpx);   //(X'X)
    gsl_blas_dgemv(CblasTrans, 1, set->data, y_data, 0, xpy);       //(X'y)
    xpxinvxpy(set->data, y_data, xpx, xpy, out);
    gsl_vector_free(y_data); gsl_vector_free(xpy);

    if ((ep == NULL) || (ep->destroy_data==0)){
        apop_data_free(set);
    }
    return out;
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
\param uses     If NULL, do everything; else, produce those \ref apop_estimate elements which you specify. You always get the parameters and never get the log likelihood.

\bug The cross-variances are assumed to be zero, which is wholeheartedly false. It's not too big a deal because nobody ever uses them for anything.
*/
/*
apop_estimate * apop_partitioned_OLS(apop_data *set1, apop_data *set2, gsl_matrix *m1, gsl_matrix *m2, apop_inventory *uses){
apop_inventory  actual_uses    = apop_inventory_filter(uses, apop_GLS.inventory_filter);
    prep_inventory_names(set1->names);
apop_estimate	*out1, *out2,
                *out	= apop_estimate_alloc(set1->data->size1, set1->data->size2 + set2->data->size2, set1->names, actual_uses);
gsl_matrix      *t1,*t2, *augmented_first_matrix, *zero1, *zero2,
                *dependent  =gsl_matrix_alloc(set1->data->size1, 1);
gsl_vector_view d;
int             i;


    if(m1==NULL){
    gsl_matrix 	*xpx1 		= gsl_matrix_calloc(set1->data->size1, set1->data->size1),
                *xpxinv1, *xxpxinv1, *xxpxinvx1;
        m1        = gsl_matrix_alloc(set1->data->size2, set1->data->size2);
	    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, set1->data, set1->data, 0, xpx1);	//(X'X)
	    apop_det_and_inv(xpx1, &xpxinv1, 0, 1);		
        gsl_matrix_free(xpx1);
        xxpxinv1   = gsl_matrix_calloc(set1->data->size1, set1->data->size1);
	    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, set1->data, xpxinv1, 0, xxpxinv1);	//X(X'X)^{-1}
        gsl_matrix_free(xpxinv1);
        xxpxinvx1  = gsl_matrix_calloc(set1->data->size1, set1->data->size1);
	    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1, xxpxinv1, set1->data, 0, xxpxinvx1);	//X(X'X)^{-1}X'
        gsl_matrix_free(xxpxinv1);
        gsl_matrix_set_identity(m1);
        gsl_matrix_sub(m1, xxpxinvx1);   // M = I - X(X'X)^{-1}X'
        gsl_matrix_free(xxpxinvx1);
    }

    if(m2==NULL){
    gsl_matrix 	*xpx2 		= gsl_matrix_calloc(set2->data->size1, set2->data->size1),
                *xpxinv2, *xxpxinv2, *xxpxinvx2;
        m2        = gsl_matrix_alloc(set2->data->size2, set2->data->size2);
	    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, set2->data, set2->data, 0, xpx2);	//(X'X)
	    apop_det_and_inv(xpx2, &xpxinv2, 0, 1);		
        gsl_matrix_free(xpx2);
        xxpxinv2   = gsl_matrix_calloc(set2->data->size1, set2->data->size1), 
	    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, set2->data, xpxinv2, 0, xxpxinv2);	//X(X'X)^{-1}
        gsl_matrix_free(xpxinv2);
        xxpxinvx2  = gsl_matrix_calloc(set2->data->size1, set2->data->size1);
	    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1, xxpxinv2, set2->data, 0, xxpxinvx2);	//X(X'X)^{-1}X'
        gsl_matrix_free(xxpxinv2);
        gsl_matrix_set_identity(m2);
        gsl_matrix_sub(m2, xxpxinvx2);   // M = I - X(X'X)^{-1}X'
        gsl_matrix_free(xxpxinvx2);
    }


        //the first matrix needs a dependent variable column, then we
        //can regress.
    d                       = gsl_matrix_column(set2->data, 0);
    gsl_matrix_set_col(dependent, 0, &(d.vector));
    augmented_first_matrix  = apop_matrix_stack(dependent, set1->data, 'c');
apop_data *newfirst = malloc(sizeof(apop_data));
    newfirst->data  = augmented_first_matrix;
    newfirst->names = set1->names;
    out1    = apop_estimate_GLS(newfirst, &actual_uses, m2);
    out2    = apop_estimate_GLS(set2, &actual_uses, m1);
    for (i=0; i< out1->parameters->size; i++)
        gsl_vector_set(out->parameters, i, gsl_vector_get(out1->parameters, i));
    for (   ; i< out2->parameters->size; i++)
        gsl_vector_set(out->parameters, i, gsl_vector_get(out2->parameters, i - out1->parameters->size));

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

//Cut and pasted from the GNU std library documentation:
static int compare_doubles (const void *a, const void *b) {
    const double *da = (const double *) a;
    const double *db = (const double *) b;
    return (*da > *db) - (*da < *db);
}

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
\return out
*/
apop_data * apop_produce_dummies(gsl_vector *in, int keep_first){
size_t      i, index,
            prior_elmt_ctr  = 107,
            elmt_ctr        = 0;
apop_data   *out; 
double      *elmts          = malloc(sizeof(double)),
            val;
char        n[1000];

    //first, create an ordered list of unique elements.
    //Note to the reader: this may or may not be more efficient on the
    //db side.
    for (i=0; i< in->size; i++){
        if (prior_elmt_ctr != elmt_ctr)
            elmts   = realloc(elmts, sizeof(double)*(elmt_ctr+1));
        prior_elmt_ctr  = elmt_ctr;
        val             =  gsl_vector_get(in, i);
        lsearch (&val, elmts, &elmt_ctr, sizeof(double), compare_doubles);
    }
    qsort(elmts, elmt_ctr, sizeof(double), compare_doubles);

    //Now go through the input vector, and for row i find the posn of the vector's
    //name in the element list created above (j), then change (i,j) in
    //the dummy matrix to one.
    if (keep_first)     out  = apop_data_alloc(in->size, elmt_ctr+1);
    else                out  = apop_data_alloc(in->size, elmt_ctr);
    gsl_matrix_set_zero(out->data);
    for (i=0; i< in->size; i++){
        val     = gsl_vector_get(in, i);
        index   = ((long int)bsearch(&val, elmts, elmt_ctr, sizeof(double), compare_doubles) - (long int)elmts)/sizeof(double);
        if (keep_first)
            gsl_matrix_set(out->data, i, index,1); 
        else if (index>0)
            gsl_matrix_set(out->data, i, index-1, 1); 
        //else don't keep first, and index==0; throw it out.
    }
    //Add names:
    i   = (keep_first) ? 0 : 1;
    for (  ; i< elmt_ctr; i++){
        sprintf(n,"dummy %g", elmts[i]);
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
apop_estimate *apop_estimate_fixed_effects_OLS(apop_data *data, apop_inventory *uses, gsl_vector *categories){
apop_data *dummies = apop_produce_dummies(categories, 0);
    apop_data_stack(data, dummies, 'c');
    return apop_OLS.estimate(dummies,  uses, NULL);
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

\ingroup regression
  */
apop_data *apop_estimate_correlation_coefficient(apop_estimate *in){
double          sse, sst, rsq;
size_t          obs     = in->data->data->size1;
size_t          indep_ct= in->data->data->size2 - 1;
gsl_vector      *y      = gsl_vector_alloc(obs);
gsl_vector_view v;  
apop_data       *out    = apop_data_alloc(5,1);
    //notice that the sum of squares is just vector dot vector.
    gsl_blas_ddot(in->residuals, in->residuals, &sse);

    //copy the dependent var to a temp vector, subtract means, square.
    v   = gsl_matrix_column(in->data->data, 0);
    gsl_vector_memcpy(y, &(v.vector));
    gsl_vector_add_constant(y, - apop_vector_mean(y));
    gsl_blas_ddot(y, y, &sst);

    rsq = 1 - (sse/sst);

    apop_name_add(out->names, "R_squared", 'r');
    gsl_matrix_set(out->data, 0, 0, rsq);
    apop_name_add(out->names, "R_squared_adj", 'r');
    gsl_matrix_set(out->data, 1, 0, rsq - ((indep_ct -1) /(obs - indep_ct)) * (1-rsq) );
    apop_name_add(out->names, "SSE", 'r');
    gsl_matrix_set(out->data, 2, 0, sse);
    apop_name_add(out->names, "SST", 'r');
    gsl_matrix_set(out->data, 3, 0, sst);
    apop_name_add(out->names, "SSR", 'r');
    gsl_matrix_set(out->data, 4, 0, sst - sse);

    gsl_vector_free(y);
    return out;
}

/** A synonym for \ref apop_estimate_correlation_coefficient, q.v. 
 \ingroup regression
 */
apop_data *apop_estimate_r_squared(apop_estimate *in){
    return apop_estimate_correlation_coefficient(in);
}
