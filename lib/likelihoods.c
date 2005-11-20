/** \file likelihoods.c	The MLE functions. Call them with an \ref apop_model.

This file includes a number of distributions and models whose parameters
one would estimate using maximum likelihood techniques.

Each typically includes four functions: the likelihood function, the 
derivative of the likelihood function, a function that calls both of them,
and a user-usable function which takes in data and a blank vector, fills 
the vector with the most likely parameters, and returns the likelihood
of those parameters.

At the bottom are the maximum likelihood procedures themselves. There
are two: the no-derivative version and the with-derivative version.
Use the with-derivative version wherever possible---in fact, it is at
the moment entirely unused, but is just here for future use.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL.
*/
#include "likelihoods.h"
#include <assert.h>
#include <gsl/gsl_deriv.h>


		///////////////////////
		//MLE support functions
		///////////////////////

//Including numerical differentiation and a couple of functions to
//negate the likelihood fns without bothering the user.

/** If you are taking a numerical gradient, you need to set this
 variable to the function for which the gradient will be taken. 
 \ingroup global_vars */
double (*apop_fn_for_derivative) (const gsl_vector *beta, void *d);


typedef struct grad_params{
	gsl_vector	*beta;
	void		*d;
	int		dimension;
} grad_params;

static double one_d(double b, void *p){
grad_params	*params	= p;
	gsl_vector_set(params->beta, params->dimension, b);
	return apop_fn_for_derivative(params->beta, params->d);
}

/**The GSL provides one-dimensional numerical differentiation; here's
 the multidimensional extension.
 
 Before using this fn., you need to set \ref apop_fn_for_derivative to
 the function you are interested (i.e., rather than passing the function
 as an argument, it's passed as a global variable.) E.g., 
 \code
 apop_fn_for_derivative = apop_exponential.log_likelihood;
 \endcode

 \ingroup linear_algebra
 \todo This fn has a hard-coded tolerance (1e-4).
 */
void apop_numerical_gradient(const gsl_vector *beta, void *d , gsl_vector *out){
int		i;
gsl_function	F;
double		result, err;
grad_params 	gp;
	gp.beta		= gsl_vector_alloc(beta->size);
	gp.d		= d;
	F.function	= one_d;
	F.params	= &gp;
	for (i=0; i< beta->size; i++){
		gp.dimension	= i;
		gsl_vector_memcpy(gp.beta, beta);
		gsl_deriv_central(&F, gsl_vector_get(beta,i), 1e-4, &result, &err);
		gsl_vector_set(out, i, result);
	}
}


/** The MLEs don't return everything you could ask for, so this pares
down the input inventory to a modified output mle.

\todo This should really just modify the input inventory. */
void prep_inventory_mle(apop_inventory *in, apop_inventory *out){
//These are the rules going from what you can ask for to what you'll get.
	if (in == NULL){ 	//then give the user the works.
		apop_inventory_set(out, 1);
	//OK, some things are not yet implemented.
	out->residuals		= 0;
		return;
	}//else:
	apop_inventory_copy(*in, out);
	out->log_likelihood	= 1;
	out->parameters		= 1;
	//OK, some things are not yet implemented.
	out->names		= 0;
	out->confidence		= 0;
	out->residuals		= 0;
	if (out->confidence==1)
		out->covariance = 1;
}


/* They always tell you to just negate your likelihood function to turn
a minimization routine into a maximization routine---and this is the
sort of annoying little detail that Apophenia is intended to take care
of for you. The next few functions do the negation, so you have one
less sign that you have to remember. 


*/
typedef struct negshell_params{
	void *d;
	apop_model	model;
} negshell_params;

static double negshell (const gsl_vector * beta, void * params){
negshell_params	*p			= params;
/*
gsl_vector *returned_beta		= gsl_vector_alloc(beta->size);
int		i;
double		penalty, 
		total_penalty		= 0;
static gsl_vector *last_good_beta	= NULL;
static double	base_for_penalty	= 0;
	for (i=0; i< p->model.constraint_ct; i++){
		penalty	= (p->model.constraint)[i](beta, p->d, returned_beta);
		if (penalty > 0)
			total_penalty	-=  penalty * p->model.log_likelihood(returned_beta, p->d);
	}
	if (total_penalty == 0)
	*/
	return - p->model.log_likelihood(beta, p->d);
	////else
	//return total_penalty;
}

static void dnegshell (const gsl_vector * beta, void * params, gsl_vector * g){
negshell_params *p	= params;
	p->model.dlog_likelihood(beta, p->d, g);
	gsl_vector_scale(g, -1);
}

static void fdf_shell(const gsl_vector *beta, void *params, double *f, gsl_vector *df){
negshell_params *p	= params;

	if (p->model.fdf==NULL){
		*f	= negshell(beta, params);
		dnegshell(beta, params, df);
	} else	{
		p->model.fdf(beta, p->d, f, df);
		(*f) 	*= -1;
		gsl_vector_scale(df, -1);
	}
}


			//////////////////////////////////////////
			//The max likelihood functions themselves. 
			//Mostly straight out of the GSL manual.
			//////////////////////////////////////////

//The next two variables are just for the current_gradient fn.
static size_t	current_dimension;
static double 	(*likelihood_for_hess) (const gsl_vector *beta, void *d);
void		(*dlikelihood_for_hess) (const gsl_vector *beta, void *d, gsl_vector *out);

	//Turn the full vector of derivatives into a single number.
	//Is this ever inefficient for large parameter spaces.
double gradient_for_hessian(const gsl_vector *beta, void * d){
double 		(*tmp) (const gsl_vector *beta, void *d);
gsl_vector	*waste	= gsl_vector_alloc(beta->size);
double		return_me;

	if (dlikelihood_for_hess == apop_numerical_gradient){	//then we need a numerical derivative.
		tmp			= apop_fn_for_derivative;
		apop_fn_for_derivative 	= likelihood_for_hess;
		apop_numerical_gradient(beta, d , waste);
		apop_fn_for_derivative 	= tmp;
	}
	else 	dlikelihood_for_hess(beta, d , waste);
	return_me	= gsl_vector_get(waste, current_dimension);
	gsl_vector_free(waste);
	return return_me;
}

/** finds the Hessian.

Calculate the matrix of second derivatives of a function.

If the function is a log likelihood, then the information matrix is the inverse of the negation of this.


\todo This is inefficient, by an order of n, where n is the number of parameters.
\ingroup linear_algebra
 */
gsl_matrix * apop_numerical_hessian(apop_model dist, gsl_vector *beta, void * d){
gsl_vector_view	v;
gsl_matrix	*out	= gsl_matrix_alloc(beta->size, beta->size);
	likelihood_for_hess	= dist.log_likelihood;
	if (dist.dlog_likelihood)
		dlikelihood_for_hess	= dist.dlog_likelihood;
	else	dlikelihood_for_hess	= apop_numerical_gradient;
	apop_fn_for_derivative	= gradient_for_hessian;
	for (current_dimension=0; current_dimension< beta->size; current_dimension++){
		v	= gsl_matrix_row(out, current_dimension);
		apop_numerical_gradient(beta, d , &(v.vector));
	}
	return out;
}

static void calculate_var_covar_matrix(apop_model dist, apop_estimate *est, gsl_matrix *data){
int		i;
double 		(*tmp) (const gsl_vector *beta, void *d);
gsl_matrix	*hessian;
	tmp			= apop_fn_for_derivative;
	//The information matrix is the inverse of the negation of the
	//hessian.
	hessian			= apop_numerical_hessian(dist, est->parameters, data);
	gsl_matrix_scale(hessian, -1);
	apop_det_and_inv(hessian, &(est->covariance), 0, 1);
	gsl_matrix_free(hessian);

	//Confidence intervals are just a cute convenience. We're
	//assuming Normality, which only works asymptotically anyway.
	if (est->uses.confidence == 0)
		return;
	//else:
	for (i=0; i< est->parameters->size; i++) // confidence[i] = |1 - (1-N(Mu[i],sigma[i]))*2|
		gsl_vector_set(est->confidence, i,
			fabs(1 - (1 - gsl_cdf_gaussian_P(gsl_vector_get(est->parameters, i), 
			gsl_matrix_get(est->covariance, i, i)))*2));
	apop_fn_for_derivative 	= tmp;
}


/*
	//find the variance-covariance matrix, using $df/d\theta \cdot df/d\theta$
static void calculate_var_covar_matrix(apop_model dist, apop_estimate *est, int betasize, gsl_matrix *data){
int		i;
gsl_vector	*diff;
gsl_matrix	*pre_cov;
gsl_vector_view v;
	pre_cov			= gsl_matrix_alloc(betasize, betasize);
	//estimate->covariance	= gsl_matrix_alloc(betasize, betasize);
	diff			= gsl_vector_alloc(betasize);
	dist.dlog_likelihood(est->parameters, data, diff);
	for (i=0; i< betasize; i++){
		gsl_matrix_set_row(pre_cov, i, diff);
	}
	for (i=0; i< betasize; i++){
		v	= gsl_matrix_column(pre_cov, i);
		gsl_vector_scale(&(v.vector), -gsl_vector_get(diff, i));
	}
	apop_det_and_inv(pre_cov, &(est->covariance), 0, 1);
	gsl_matrix_free(pre_cov);
	gsl_vector_free(diff);

	if (est->uses.confidence == 0)
		return;
	//else:
	for (i=0; i<betasize; i++) // confidence[i] = |1 - (1-N(Mu[i],sigma[i]))*2|
		gsl_vector_set(est->confidence, i,
			fabs(1 - (1 - gsl_cdf_gaussian_P(gsl_vector_get(est->parameters, i), 
			gsl_matrix_get(est->covariance, i, i)))*2));
}
*/


/** The maximum likelihood calculations, given a derivative of the log likelihood.

If no derivative exists, will calculate a numerical gradient.

\param data	the data matrix
\param uses	an inventory, which will be pared down and folded into the output \ref apop_estimate
\param	dist	the \ref apop_model object: waring, probit, zipf, &amp;c.
\param	starting_pt	an array of doubles suggesting a starting point. If NULL, use zero.
\param step_size	the initial step size.
\param tolerance	the precision the minimizer uses. Only vaguely related to the precision of the actual var.
\param verbose		Y'know.
\return	an \ref apop_estimate with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.

  \todo readd names */
apop_estimate *	apop_maximum_likelihood_w_d(gsl_matrix * data, apop_inventory *uses,
			apop_model dist, double *starting_pt, double step_size, double tolerance, int verbose){
gsl_multimin_function_fdf 	minme;
gsl_multimin_fdfminimizer 	*s;
gsl_vector 			*x;
int				iter 	= 0, 
				status,
				betasize= dist.parameter_ct;
apop_inventory			actual_uses;
apop_estimate			*est;
negshell_params			nsp;
	if (betasize == -1)	betasize = data->size2 - 1;
	prep_inventory_mle(uses, &actual_uses);
	s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, betasize);
	//s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, betasize);
	est	= apop_estimate_alloc(data->size1, betasize, NULL, actual_uses);
	if (dist.dlog_likelihood == NULL ){
		dist.dlog_likelihood 	= apop_numerical_gradient;
		apop_fn_for_derivative	= dist.log_likelihood;
	}
	if (starting_pt==NULL){
		x	= gsl_vector_alloc(betasize);
  		gsl_vector_set_all (x,  0);
	}
	else 	apop_convert_array_to_vector(starting_pt, &x, betasize);
	nsp.d		= data;
	nsp.model	= dist;
	minme.f		= negshell;
	minme.df	= dnegshell;
	minme.fdf	= fdf_shell;
	minme.n		= betasize;
	minme.params	= &nsp;
	gsl_multimin_fdfminimizer_set (s, &minme, x, step_size, tolerance);
      	do { 	iter++;
		status 	= gsl_multimin_fdfminimizer_iterate(s);
		if (status) 	break; 
		status = gsl_multimin_test_gradient(s->gradient, tolerance);
        	if (verbose){
	        	printf ("%5i %.5f  f()=%10.5f gradient=%.3f\n", iter, gsl_vector_get (s->x, 0),  s->f, gsl_vector_get(s->gradient,0));
		}
        	if (status == GSL_SUCCESS){
			est->status	= 0;
		   	if(verbose)	printf ("Minimum found.\n");
		}
       	 }
	//while (status == GSL_CONTINUE && iter < MAX_ITERATIONS_w_d);
	while (status != GSL_SUCCESS && iter < MAX_ITERATIONS_w_d);
	if(iter==MAX_ITERATIONS_w_d) {
		est->status	= 1;
		if (verbose) printf("No min!!\n");
	}
	//Clean up, copy results to output estimate.
	gsl_vector_memcpy(est->parameters, s->x);
	gsl_multimin_fdfminimizer_free(s);
	if (starting_pt==NULL) 
		gsl_vector_free(x);
	est->log_likelihood	= dist.log_likelihood(est->parameters, data);
	if (est->uses.covariance) 
		calculate_var_covar_matrix(dist, est, data);
	return est;
}

apop_estimate *	apop_maximum_likelihood_no_d(gsl_matrix * data, apop_inventory *uses,
			apop_model dist, double *starting_pt, double step_size, double tolerance, int verbose){
int			status,
			iter 		= 0,
			betasize	= dist.parameter_ct;
size_t 			i;
gsl_multimin_function 	minme;
gsl_multimin_fminimizer *s;
gsl_vector 		*x, *ss;
double			size;
apop_inventory		actual_uses;
apop_estimate		*est;
negshell_params		nsp;
	if (betasize == -1)	betasize = data->size2 - 1;
	s	= gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, betasize);
	ss	= gsl_vector_alloc(betasize);
	prep_inventory_mle(uses, &actual_uses);
	est	= apop_estimate_alloc(data->size1, betasize, NULL, actual_uses);
	est->status	= 1;	//assume failure until we score a success.
	if (starting_pt==NULL)
  		gsl_vector_set_all (x,  0);
	else
		apop_convert_array_to_vector(starting_pt, &x, betasize);
  	gsl_vector_set_all (ss,  step_size);
	nsp.model	= dist;
	nsp.d		= data;
	minme.f		= negshell;
	minme.n		= betasize;
	minme.params	= &nsp;
	gsl_multimin_fminimizer_set (s, &minme, x,  ss);
      	do { 	iter++;
		status 	= gsl_multimin_fminimizer_iterate(s);
		if (status) 	break; 
		size	= gsl_multimin_fminimizer_size(s);
	       	status 	= gsl_multimin_test_size (size, tolerance); 
		if(verbose){
			printf ("%5d ", iter);
			for (i = 0; i < betasize; i++) {
				printf ("%8.3e ", gsl_vector_get (s->x, i)); } 
			printf ("f()=%7.3f size=%.3f\n", s->fval, size);
       			if (status == GSL_SUCCESS) {
				est->status	= 0;
	   			if(verbose){
					printf ("Minimum found at:\n");
					printf ("%5d ", iter);
					for (i = 0; i < betasize; i++) {
						printf ("%8.3e ", gsl_vector_get (s->x, i)); } 
					printf ("f()=%7.3f size=%.3f\n", s->fval, size);
				}
			}
		}
      	} while (status == GSL_CONTINUE && iter < MAX_ITERATIONS);
	if (iter == MAX_ITERATIONS && verbose)
		printf("Minimization reached maximum number of iterations.");

	gsl_vector_memcpy(est->parameters, s->x);
	gsl_multimin_fminimizer_free(s);
	est->log_likelihood	= dist.log_likelihood(est->parameters, data);
	//if (est->uses.confidence == 0)
		return est;
}


/** The maximum likelihood calculations

\param data	the data matrix
\param uses	an inventory, which will be pared down and folded into the output \ref apop_estimate
\param	dist	the \ref apop_model object: waring, probit, zipf, &amp;c.
\param params	an \ref apop_estimation_params structure, featuring:<br>
starting_pt	an array of doubles suggesting a starting point. If NULL, use zero.<br>
step_size	the initial step size.<br>
tolerance	the precision the minimizer uses. Only vaguely related to the precision of the actual var.<br>
verbose		Y'know.<br>
method		0: Nelder-Mead simplex<br>
			1: conjugate gradient. If no gradient is available, use numerical approximations. (default)
\return	an \ref apop_estimate with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.

  \todo re-add names 
 \ingroup mle */
apop_estimate *	apop_maximum_likelihood(gsl_matrix * data, apop_inventory *uses,
			apop_model dist, apop_estimation_params params){

	if (params.method==0)
		return apop_maximum_likelihood_no_d(data, uses, dist, params.starting_pt, params.step_size, params.tolerance, params.verbose);
	//else:
	return apop_maximum_likelihood_w_d(data, uses, dist, params.starting_pt, params.step_size, params.tolerance, params.verbose);
}


/** This function goes row by row through <tt>m</tt> and calculates the
likelihood of the given row, putting the result in <tt>v</tt>. 
You can use this to find the variance of the estimator if other means fail.

\param m 	A GSL matrix, exactly like those used for probit, Waring, Gamma, &c MLEs.

\param v	A vector which will hold the likelihood of each row of m. Declare but do not allocate.

\param dist	An \ref apop_model object whose log likelihood function you'd like to use.

\param fn_beta		The parameters at which you will evaluate the likelihood. If <tt>e</tt> is an \ref
			apop_estimate, then one could use <tt>e->parameters</tt>.

This functions is used in the sample code in the \ref mle section.

\ingroup mle */
void apop_make_likelihood_vector(gsl_matrix *m, gsl_vector **v, 
				apop_model dist, gsl_vector* fn_beta){
gsl_matrix_view mm;
int             i;
	*v	= gsl_vector_alloc(m->size1);
        for(i=0; i< m->size1; i++){
                mm      = gsl_matrix_submatrix(m,i,0, 1,m->size2);      //get a single row
                gsl_vector_set(*v, i, dist.log_likelihood(fn_beta, &(mm.matrix)));
        }
}
