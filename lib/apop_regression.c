/** \file apop_regression.c	Generally, if it assumes something is  Normally distributed, it's here.\n
 Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
 \author Ben Klemens
 */

/** \defgroup regression  OLS/GLS: The linear projection methods */
/** \defgroup ttest  T-tests: comparing two vectors */
/** \defgroup asst_tests  Various means of hypothesis testing.*/

#include "regression.h"
#include <gsl/gsl_blas.h>
#include "stats.h"
#include "linear_algebra.h"
#include <apophenia/types.h>

extern int apop_verbose;

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
double		a_avg	= apop_mean(a),
		a_var	= apop_var(a),
		b_avg	= apop_mean(b),
		b_var	= apop_var(b),
		stat	= (a_avg - b_avg)/ sqrt(b_var/(b_count-1) + a_var/(a_count-1));
	if (apop_verbose){
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
double		avg	= apop_mean(diff),
		var	= apop_var(diff),
		stat	= avg/ sqrt(var/(count-1));
	if (apop_mean(diff) < 0) 
		factor = -1;
	gsl_vector_free(diff);
	if (apop_verbose){
		printf("avg diff: %g; diff std dev: %g; count: %i; t-statistic: %g.\n", avg, sqrt(var), count, stat);
	}
	return factor * two_tailify(gsl_cdf_tdist_P(stat, count-1));
}

void prep_inventory_OLS(apop_name *n, apop_inventory *in, apop_inventory *out){
//These are the rules going from what you can ask for to what you'll get.
	if (in == NULL) 	//then give the user the works.
		apop_inventory_set(out, 1);
	else {
		apop_inventory_copy(*in, out);
		out->covariance	= 1;	//always calculated.
	}
	out->log_likelihood	= 0;
	out->parameters		= 1;
	if (n == NULL)
		out->names		= 0;
	else {		//shift first col to depvar, rename first col "one".
		out->names		= 1;
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

\param data
The first column is the dependent variable, and the remaining columns the independent. Is destroyed in the process, so make a copy beforehand if you need.

\param sigma 
A known variance-covariance matrix, of size <tt>(data->size1, data->size1)</tt>. Survives the function intact. The first column refers to the constant unit vector, so it's always zero.

\param n    An \c apop_name structure, describing the columns.

\param uses 
If NULL, do everything; else, produce those \ref apop_estimate elements which you specify. You always get the parameters and never get the log likelihood.

\return
A pointer to an \ref apop_estimate structure with the appropriate elements filled. See the description in \ref apop_OLS .

\todo 
Since the first column and row of the var/covar matrix is always zero, users shouldn't have to make it.
 */
apop_estimate * apop_GLS(gsl_matrix *data, gsl_matrix *sigma, apop_name * n, apop_inventory *uses){
apop_inventory	actual_uses;
	prep_inventory_OLS(n, uses, &actual_uses);
apop_estimate	*out		= apop_estimate_alloc(data->size1, data->size2, n, actual_uses);
gsl_vector 	*y_data		= gsl_vector_alloc(data->size1);
gsl_matrix 	*temp		= gsl_matrix_calloc(data->size2, data->size1);
gsl_vector 	*xsy 		= gsl_vector_calloc(data->size2);
gsl_matrix 	*xsx 		= gsl_matrix_calloc(data->size2, data->size2);
gsl_matrix 	*sigma_inverse;	//= gsl_matrix_alloc(data->size1, data->size1);
gsl_vector_view	v 		= gsl_matrix_column(data, 0);
	gsl_matrix_get_col(y_data, data, 0);
	gsl_vector_set_all(&(v.vector), 1);	//affine: first column is ones.
	apop_det_and_inv(sigma, &sigma_inverse, 0, 1);					//find sigma^{-1}
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, data, sigma_inverse, 0, temp); 	//temp = X' \sigma^{-1}.
	gsl_matrix_free(sigma_inverse);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, temp, data, 0, xsx);    		//(X' \sigma^{-1} X)
	gsl_blas_dgemv(CblasNoTrans, 1, temp, y_data, 0, xsy);     			//(X' \sigma^{-1} y)
	gsl_matrix_free(temp);
	xpxinvxpy(data, y_data, xsx, xsy, out);
	gsl_vector_free(y_data); gsl_vector_free(xsy);
	return out;
}

/** ordinary least squares.

\ingroup regression

The first column is the dependent variable, the remaining columns are the independent variables. NB: \c data is destroyed by this function. If you want to keep it, make a copy beforehand.

\param data
The first column is the dependent variable, and the remaining columns the independent. Is destroyed in the process, so make a copy beforehand if you need.

\param n
An \c apop_name structure, specifying which outputs you want.

\param uses If <tt>NULL</tt>, then you get everything.  If a pointer to
an \ref apop_inventory , then you get what you ask for. Log likelihood is
not calculated; you always get the parameter estimates.

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
gsl_matrix      *data;
apop_estimate   *est;
apop_name       *n;
     apop_convert_text_to_db("data","d",NULL);
     data       = apop_query_to_matrix("select * from d");
     n          = apop_db_get_names();
     est  = apop_OLS(data, n, NULL);
     apop_estimate_print(est);
     return 0;
}
\endcode

If you saved this code to <tt>sample.c</tt>, then you can compile it with
\verbatim
gcc sample.c -lapophenia -lgsl -lgslcblas -lsqlite3 -o run_me
\endverbatim

and then run it with <tt>./run_me</tt>. Alternatively, you may prefer to compile the program using a \ref makefile .


 */
apop_estimate * apop_OLS(gsl_matrix *data, apop_name * n, apop_inventory *uses){
apop_inventory	actual_uses;
	prep_inventory_OLS(n, uses, &actual_uses);
apop_estimate	*out		= apop_estimate_alloc(data->size1, data->size2, n, actual_uses);
gsl_vector 	*y_data		= gsl_vector_alloc(data->size1);
gsl_vector 	*xpy 		= gsl_vector_calloc(data->size2);
gsl_matrix 	*xpx 		= gsl_matrix_calloc(data->size2, data->size2);
gsl_vector_view	v 		= gsl_matrix_column(data, 0);
	gsl_matrix_get_col(y_data, data, 0);
	gsl_vector_set_all(&(v.vector), 1);	//affine: first column is ones.
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, data, data, 0, xpx);	//(X'X)
	gsl_blas_dgemv(CblasTrans, 1, data, y_data, 0, xpy);     	//(X'y)
	xpxinvxpy(data, y_data, xpx, xpy, out);
	gsl_vector_free(y_data); gsl_vector_free(xpy);
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

\param data1    the first half of the data set
\param data2    the second half of the data set
\param  M1      If you have these, give them to me and I won't waste time calculating them.
\param  M2      If you have these, give them to me and I won't waste time calculating them.
\param n1        An \c apop_name structure, describing the columns of the first data set.
\param n2        An \c apop_name structure, describing the columns of the first data set.
\param uses     If NULL, do everything; else, produce those \ref apop_estimate elements which you specify. You always get the parameters and never get the log likelihood.

\bug The cross-variances are assumed to be zero, which is wholeheartedly false. It's not too big a deal because nobody ever uses them for anything.
*/
apop_estimate * apop_partitioned_OLS(gsl_matrix *data1, gsl_matrix *data2, gsl_matrix *m1, gsl_matrix *m2, 
                                                    apop_name * n1, apop_name * n2, apop_inventory *uses){
apop_inventory	actual_uses;
	prep_inventory_OLS(n1, uses, &actual_uses); //FIX for more names
apop_estimate	*out1, *out2,
                *out	= apop_estimate_alloc(data1->size1, data1->size2 + data2->size2, n1, actual_uses);
gsl_matrix      *t1,*t2, *augmented_first_matrix, *zero1, *zero2,
                *dependent  =gsl_matrix_alloc(data1->size1, 1);
gsl_vector_view d;
int             i;


    if(m1!=NULL){
    gsl_matrix 	*xpx1 		= gsl_matrix_calloc(data1->size1, data1->size1),
                *xpxinv1,
                *xxpxinv1   = gsl_matrix_calloc(data1->size1, data1->size1), 
                *xxpxinvx1  = gsl_matrix_calloc(data1->size1, data1->size1);
        m1        = gsl_matrix_alloc(data1->size2, data1->size2);
	    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, data1, data1, 0, xpx1);	//(X'X)
	    apop_det_and_inv(xpx1, &xpxinv1, 0, 1);		
        gsl_matrix_free(xpx1);
	    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, data1, xpxinv1, 0, xxpxinv1);	//X(X'X)^{-1}
        gsl_matrix_free(xpxinv1);
	    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1, xxpxinv1, data1, 0, xxpxinvx1);	//X(X'X)^{-1}X'
        gsl_matrix_free(xxpxinv1);
        gsl_matrix_set_identity(m1);
        gsl_matrix_sub(m1, xxpxinvx1);   // M = I - X(X'X)^{-1}X'
        gsl_matrix_free(xxpxinvx1);
    }


    if(m2!=NULL){
    gsl_matrix 	*xpx2 		= gsl_matrix_calloc(data2->size2, data2->size2),
                *xpxinv2,
                *xxpxinv2   = gsl_matrix_calloc(data2->size2, data2->size2), 
                *xxpxinvx2  = gsl_matrix_calloc(data2->size2, data2->size2);
        m2        = gsl_matrix_alloc(data2->size2, data2->size2);
	    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 2, data2, data2, 0, xpx2);	//(X'X)
	    apop_det_and_inv(xpx2, &xpxinv2, 0, 1);		
        gsl_matrix_free(xpx2);
	    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, data2, xpxinv2, 0, xxpxinv2);	//X(X'X)^{-1}
        gsl_matrix_free(xpxinv2);
	    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1, xxpxinv2, data2, 0, xxpxinvx2);	//X(X'X)^{-1}X'
        gsl_matrix_free(xxpxinv2);
        gsl_matrix_set_identity(m2);
        gsl_matrix_sub(m2, xxpxinvx2);   // M = I - X(X'X)^{-1}X'
        gsl_matrix_free(xxpxinvx2);
    }


        //the first matrix needs a dependent variable column, then we
        //can regress.
    d                       = gsl_matrix_column(data2, 0);
    gsl_matrix_set_col(dependent, 0, &(d.vector));
    augmented_first_matrix  = apop_matrix_stack(dependent, data1, 'r');
    out1    = apop_GLS(augmented_first_matrix, m2, n1, &actual_uses);
    out2    = apop_GLS(data2, m1, n2, &actual_uses);
    for (i=0; i< out1->parameters->size; i++)
        gsl_vector_set(out->parameters, i, gsl_vector_get(out1->parameters, i));
    for (   ; i< out2->parameters->size; i++)
        gsl_vector_set(out->parameters, i, gsl_vector_get(out2->parameters, i - out1->parameters->size));

    //The covariance matrix is a stacking-up of the above matrices. I
    //cheat here: the cross-variances are assumed zero, which is
    //blatantly false.
    zero1           = gsl_matrix_calloc(out1->covariance->size1, out2->covariance->size2);
    t1              = apop_matrix_stack(out1->covariance, zero1, 'r');
    gsl_matrix_free(zero1);
    zero2           = gsl_matrix_calloc(out2->covariance->size1, out1->covariance->size2);
    t2              = apop_matrix_stack(zero2, out2->covariance, 'r');
    gsl_matrix_free(zero1);
    out->covariance = apop_matrix_stack(t1,t2,'t');
    return out;
}



/** A utility to make a matrix of dummy variables. You give me a single
vector that lists the category number for each item, and I'll return
a gsl_matrix with a single one in each row in the column specified.

\param  cats The gsl_vector of categories
\param  keep_first  if zero, return 
    a matrix where each row has a one in the (column specified MINUS
    ONE). That is, the zeroth category is dropped, the first category
    has an entry in column zero, et cetera. If you don't know why this
    is useful, then this is what you need. If you know what you're doing
    and need something special, set this to one and the first category won't be dropped.
\return out
\todo This won't work well if the categories do not include observations
zero through max. That is, the function should produce a second vector
that produces neat categories from sloppy input categories.
*/
gsl_matrix * apop_produce_dummies(gsl_vector *in, int keep_first){
int         max_category    = gsl_vector_max(in),
            i, col;
gsl_matrix  *out; 
    if (keep_first){
        out  = gsl_matrix_calloc(in->size, max_category+1);
        for (i=0; i< in->size; i++)
            gsl_matrix_set(out, i, gsl_vector_get(in, i),1); 
       }
    else{
        out  = gsl_matrix_calloc(in->size, max_category);
        for (i=0; i< in->size; i++){
            col =  gsl_vector_get(in, i);
            if (col>0)
                gsl_matrix_set(out, i, col, 1); 
        }
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
apop_estimate *apop_fixed_effects_OLS(gsl_matrix *data, gsl_vector *categories, apop_name * n, apop_inventory *uses){
gsl_matrix *dummies = apop_produce_dummies(categories, 0);
    return apop_partitioned_OLS(dummies, data, NULL, NULL, NULL, n, uses);
}

