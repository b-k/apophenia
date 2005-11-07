/** \file regression.c	Generally, if it assumes something is  Normally distributed, it's here.\n
 Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
 \author Ben Klemens
 */

/** \defgroup regression  OLS/GLS: The linear projection methods */
/** \defgroup ttest  T-tests: comparing two vectors */
/** \defgroup asst_tests  Various means of hypothesis testing.*/

#include <apophenia/regression.h>
#include <gsl/gsl_blas.h>
#include <apophenia/stats.h>
#include <apophenia/linear_algebra.h>
#include <apophenia/estimate.h>

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
\return the confidence level---if it is close to one, you can reject the null, while <tt>apop_paired_t_test(a, a)</tt> will return zero.
*/
double	apop_paired_t_test(gsl_vector *a, gsl_vector *b){
gsl_vector	*diff	= gsl_vector_alloc(a->size);
	gsl_vector_memcpy(diff, a);
	gsl_vector_sub(diff, b);
int		count	= a->size;
double		avg	= apop_mean(diff),
		var	= apop_var(diff),
		stat	= avg/ sqrt(var/(count-1));
	gsl_vector_free(diff);
	if (apop_verbose){
		printf("avg diff: %g; diff std dev: %g; count: %i; t-statistic: %g.\n", avg, sqrt(var), count, stat);
	}
	return two_tailify(gsl_cdf_tdist_P(stat, count-1));
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

\param n
An \c apop_name structure, specifying which outputs you want.

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
