/** \file apop_mle.c	The MLE functions. Call them with an \ref apop_model.

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

Copyright (c) 2006 by Ben Klemens. Licensed under the GNU GPL v2.
*/
#include "likelihoods.h"
#include <assert.h>
#include <gsl/gsl_deriv.h>
static apop_estimate * apop_annealing(apop_model m, apop_data *data, apop_estimation_params *ep);  //below.

/** \page trace_path Plotting the path of an ML estimation.

If \c apop_opts.mle_trace_path has a name of positive length, then every time
the MLE evaluates the function, then the value will be output to a table
in the database with the given name. You can then plot this table to
get an idea of the path the estimation routine used to arrive at its MLE.

First, set this variable and run the MLE. The begin/commit wrapper
speeds things up a touch, but this will clearly be slower than without
taking notes:

\code
    strcpy(apop_opts.mle_trace_path, "path");
    apop_query("begin;");
    e   = apop_zipf.estimate(...);
    apop_query("commit;");
\endcode


Then, plot using a function like the below. Notice that you will want
\c splot for 3-d variables and \c plot for 2-d. The change in width
and pointsize is to remind the eye that the lines connecting the points
only indicate the path the maximizer went along, not actual values of
the function.

\code
static void plotme(char *outfile){
FILE            *f;
gsl_matrix      *traced_path;
    f       = fopen(outfile, "w");  //overwrites. Use "a" to append.
    fprintf(f,"splot '-' with linespoints linewidth 0.5 pointsize 2\n");
    traced_path = apop_query_to_matrix("select * from %s", apop_opts.mle_trace_path);
    fclose(f);
    apop_matrix_print(traced_path, "\t", outfile);
}
\endcode

Finally, call gnuplot:
\code 
gnuplot -persist < plotme
(or)
gnuplot plotme -
\endcode

Below is a sample of the sort of output one would get:<br>
\image latex "search.gif" "An ML search, tracing out the surface of the function" width=\textwidth
\image html "search.gif" "An ML search, tracing out the surface of the function" 

\ingroup mle
*/


		///////////////////////
		//MLE support functions
		///////////////////////

typedef double 	(*apop_fn_with_void) (const gsl_vector *beta, void *d);
typedef	void 	(*apop_df_with_void)(const gsl_vector *beta, void *d, gsl_vector *gradient);
typedef	void 	(*apop_fdf_with_void)( const gsl_vector *beta, void *d, double *f, gsl_vector *df);

//Including numerical differentiation and a couple of functions to
//negate the likelihood fns without bothering the user.

/** If you are taking a numerical gradient, you need to set this
 variable to the function for which the gradient will be taken. 
 \ingroup global_vars */
apop_fn_with_void apop_fn_for_derivative;


typedef struct grad_params{
	gsl_vector	*beta;
	apop_data	*d;
	int		    dimension;
} grad_params;

static double one_d(double b, void *p){
grad_params	*params	= p;
    gsl_vector_set(params->beta, params->dimension, b);
	return apop_fn_for_derivative(params->beta, params->d);
}

/**The GSL provides one-dimensional numerical differentiation; here's
 the multidimensional extension.
 
 Before using this fn., you need to set \ref apop_fn_for_derivative to
 the function in which you are interested (i.e., rather than passing the function
 as an argument, it's passed as a global variable.) E.g., 
 \code
 apop_fn_for_derivative = apop_exponential.log_likelihood;
 \endcode

 \ingroup linear_algebra
 \todo This fn has a hard-coded tolerance (1e-4).
 */
void apop_numerical_gradient(const gsl_vector *beta, apop_data *d , gsl_vector *out){
int		        i;
gsl_function	F;
double		    result, err;
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
static void prep_inventory_mle(apop_inventory in){
//These are the rules going from what you can ask for to what you'll get.
	in.log_likelihood	= 1;
	in.parameters		= 1;
	//OK, some things are not yet implemented.
	in.names		    = 0;
	in.confidence		= 0;
	//OK, some things are not yet implemented.
	in.dependent		= 
	in.predicted		= 0;
	if (in.confidence==1)
		in.covariance = 1;
	if (in.covariance==1)
		in.confidence = 1;
}


/* They always tell you to just negate your likelihood function to turn
a minimization routine into a maximization routine---and this is the
sort of annoying little detail that Apophenia is intended to take care
of for you. The next few functions do the negation, so you have one
less sign that you have to remember. 

These fns also take care of checking constraints, if any.

*/

apop_model negshell_model;

/*
typedef struct negshell_params{
	void        *d;
	apop_model	model;
} negshell_params;
*/

static void insert_path_into_db(gsl_vector *beta, double out){
    if (beta->size == 1){
        if(!apop_table_exists(apop_opts.mle_trace_path, 0))
            apop_query("create table  %s (beta0, ll);", apop_opts.mle_trace_path);
        apop_query("insert into %s values (%g, %g);", apop_opts.mle_trace_path, gsl_vector_get(beta,0), out);
    } else {
        if(!apop_table_exists(apop_opts.mle_trace_path, 0))
            apop_query("create table  %s (beta0, beta1, ll);", apop_opts.mle_trace_path);
        apop_query("insert into %s values (%g, %g, %g);", apop_opts.mle_trace_path, gsl_vector_get(beta,0), gsl_vector_get(beta,1), out);
    }

}

typedef double apop_constify  (const gsl_vector * beta, apop_data * d);

static double negshell (gsl_vector * beta, apop_data * d){
gsl_vector 	    *returned_beta	= gsl_vector_alloc(beta->size);
double		    penalty, out    = 0, 
                done            = 0; 
static double	base_for_penalty= 0;
	if (negshell_model.constraint){
		penalty	= negshell_model.constraint(beta, d, returned_beta);
		if (penalty > 0){
			base_for_penalty	= negshell_model.log_likelihood(returned_beta, d);
			out = -base_for_penalty + penalty;
            done++;
		}
	}
    if (!done)
	    out = - negshell_model.log_likelihood(beta, d);
    if (strlen(apop_opts.mle_trace_path) > 0)
        insert_path_into_db(beta,-out);
	gsl_vector_free(returned_beta);
    return out;
}


/* The derivative-calculating routine.
If the constraint binds
    then: take the numerical derivative of negshell, which will be the
    numerical derivative of the penalty.
    else: just find dlog_likelihood. If the model doesn't have a
    dlog likelihood or the user asked to ignore it, then the main
    maximum likelihood fn replaced model.dlog_likelihood with
    apop_numerical_gradient anyway.
Finally, reverse the sign, since the GSL is trying to minimize instead of maximize.
*/
//static void dnegshell (const gsl_vector * beta, apop_data * d, gsl_vector * g){
static void dnegshell (gsl_vector * beta, apop_data * d, gsl_vector * g){
gsl_vector      *returned_beta  = gsl_vector_alloc(beta->size);
double 		    (*tmp) (const gsl_vector *beta, void *d);
	if (negshell_model.constraint && negshell_model.constraint(beta, d, returned_beta)){
		    tmp			            = (apop_fn_with_void) apop_fn_for_derivative;
		    apop_fn_for_derivative 	= (apop_fn_with_void) negshell;
		    apop_numerical_gradient(beta, d, g);
		    apop_fn_for_derivative 	= tmp;
	} else {
	    negshell_model.dlog_likelihood(beta, d, g);
	    gsl_vector_scale(g, -1);
    }
	gsl_vector_free(returned_beta);
}

//static void fdf_shell(const gsl_vector *beta, apop_data *d, double *f, gsl_vector *df){
static void fdf_shell(gsl_vector *beta, apop_data *d, double *f, gsl_vector *df){
	if (negshell_model.fdf==NULL){
		*f	= negshell(beta, d);
		dnegshell(beta, d, df);
	} else	{
		negshell_model.fdf(beta, d, f, df);
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
static double 	(*likelihood_for_hess) (const gsl_vector *beta, apop_data *d);
void		(*dlikelihood_for_hess) (const gsl_vector *beta, apop_data *d, gsl_vector *out);

	//Turn the full vector of derivatives into a single number.
	//Is this ever inefficient for large parameter spaces.
double gradient_for_hessian(const gsl_vector *beta, apop_data * d){
//double 		(*tmp) (const gsl_vector *beta, void *d);
apop_fn_with_void 		tmp;
gsl_vector	*waste	= gsl_vector_alloc(beta->size);
double		return_me;

	if (dlikelihood_for_hess == apop_numerical_gradient){	//then we need a numerical derivative.
		tmp			            = (apop_fn_with_void) apop_fn_for_derivative;
		apop_fn_for_derivative 	= (apop_fn_with_void) likelihood_for_hess;
		apop_numerical_gradient(beta, d , waste);
		apop_fn_for_derivative 	= tmp;
	}
	else 	dlikelihood_for_hess(beta, d , waste);
	return_me	= gsl_vector_get(waste, current_dimension);
	gsl_vector_free(waste);
	return return_me;
}

/** Calculates the matrix of second derivatives of the log likelihood function.

The information matrix is the inverse of the negation of this.

\todo This is inefficient, by an order of n, where n is the number of parameters.
\ingroup linear_algebra
 */
gsl_matrix * apop_numerical_second_derivative(apop_model dist, gsl_vector *beta, apop_data * d){
gsl_vector_view	v;
gsl_matrix	*out	= gsl_matrix_alloc(beta->size, beta->size);
	likelihood_for_hess	= dist.log_likelihood;
	if (dist.dlog_likelihood)
		dlikelihood_for_hess	= dist.dlog_likelihood;
	else	dlikelihood_for_hess	= apop_numerical_gradient;
	apop_fn_for_derivative	= (apop_fn_with_void) gradient_for_hessian;
	for (current_dimension=0; current_dimension< beta->size; current_dimension++){
		v	= gsl_matrix_row(out, current_dimension);
		apop_numerical_gradient(beta, d , &(v.vector));
	}
	return out;
}

/** Calculate the Hessian.

  This is a synonym for \ref apop_numerical_second_derivative, q.v.
*/
gsl_matrix * apop_numerical_hessian(apop_model dist, gsl_vector *beta, apop_data * d){
	return apop_numerical_second_derivative(dist, beta, d);
}

/** Feeling lazy? Rather than doing actual pencil-and-paper math to find
your variance-covariance matrix, just use the negative inverse of the Hessian.

\param dist	The model
\param est	The estimate, with the parameters already calculated. The var/covar matrix will be placed in est->covariance.
\param data	The data
\ingroup basic_stats
*/
void apop_numerical_covariance_matrix(apop_model dist, apop_estimate *est, apop_data *data){
//int		i;
apop_fn_with_void tmp;
gsl_matrix	    *hessian;
	tmp			= apop_fn_for_derivative;
	//The information matrix is the inverse of the negation of the hessian.
	hessian			= apop_numerical_hessian(dist, est->parameters->vector, data);
	gsl_matrix_scale(hessian, -1);
	apop_det_and_inv(hessian, &(est->covariance->matrix), 0, 1);
	gsl_matrix_free(hessian);

	if (est->estimation_params.uses.confidence == 0)
		return;
	//else:
        apop_estimate_parameter_t_tests(est);
	apop_fn_for_derivative 	= tmp;
}


/** An alias for \ref apop_numerical_covariance_matrix. Use that one. */
void apop_numerical_var_covar_matrix(apop_model dist, apop_estimate *est, apop_data *data){
    apop_numerical_covariance_matrix(dist, est, data);}


/**  Here's the [dx/df]'[dx/df] version of the Hessian calculation.

Feeling lazy? Rather than doing actual pencil-and-paper math to find
your variance-covariance matrix, just use the negative inverse of the Hessian.

\param dist	The model
\param est	The estimate, with the parameters already calculated. The var/covar matrix will be placed in est->covariance.
\param data	The data
\ingroup basic_stats
*/
/*
void apop_numerical_var_covar_matrix(apop_model dist, apop_estimate *est, gsl_matrix *data){
int		i;
int         betasize    = est->parameters->vector->size;
gsl_vector	*diff;
gsl_matrix	*pre_cov	= gsl_matrix_alloc(betasize, betasize);
gsl_vector  v;
	//estimate->covariance	= gsl_matrix_alloc(betasize, betasize);
	diff			= gsl_vector_alloc(betasize);
	dist.dlog_likelihood(est->parameters->vector, data, diff);
	for (i=0; i< betasize; i++){
		gsl_matrix_set_row(pre_cov, i, diff);
	}
	for (i=0; i< betasize; i++){
		v	= gsl_matrix_column(pre_cov, i).vector;
		gsl_vector_scale(&v, -gsl_vector_get(diff, i));
	}
	apop_det_and_inv(pre_cov, &(est->covariance->matrix), 0, 1);
	gsl_matrix_free(pre_cov);
	gsl_vector_free(diff);

	if (est->estimation_params.uses.confidence == 0)
		return;
	//else:
	for (i=0; i<betasize; i++) // confidence[i] = |1 - (1-N(Mu[i],sigma[i]))*2|
		gsl_vector_set(est->confidence, i,
			fabs(1 - (1 - gsl_cdf_gaussian_P(gsl_vector_get(est->parameters->vector, i), 
			gsl_matrix_get(est->covariance->matrix, i, i)))*2));
}
*/


/* The maximum likelihood calculations, given a derivative of the log likelihood.

If no derivative exists, will calculate a numerical gradient.

\param data	the data matrix
\param	dist	the \ref apop_model object: waring, probit, zipf, &amp;c.
\param	starting_pt	an array of doubles suggesting a starting point. If NULL, use a vector whose elements are all 0.1 (zero has too many pathological cases).
\param step_size	the initial step size.
\param tolerance	the precision the minimizer uses. Only vaguely related to the precision of the actual var.
\param verbose		Y'know.
\return	an \ref apop_estimate with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.

  \todo readd names */
static apop_estimate *	apop_maximum_likelihood_w_d(apop_data * data,
			apop_model dist, apop_estimation_params *est_params){
gsl_multimin_function_fdf 	minme;
gsl_multimin_fdfminimizer 	*s;
gsl_vector 			        *x;
int				            iter 	= 0, 
				            status  = 0,
				            betasize= dist.parameter_ct;
apop_estimate			    *est;
	if (betasize == -1)	{
        dist.parameter_ct   =
        betasize            = data->matrix->size2 - 1;
    }
    betasize    *=  (est_params ? est_params->params_per_column : 1);
	prep_inventory_mle(est_params->uses);
	est	= apop_estimate_alloc(data, dist, est_params);
    if (!est_params || est_params->method/100 ==2)
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, betasize);
    else if (est_params->method/100 == 3)
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_pr, betasize);
    else
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, betasize);
	if (dist.dlog_likelihood == NULL || (est_params && est_params->method % 100 ==1))
		dist.dlog_likelihood 	= apop_numerical_gradient;
	apop_fn_for_derivative	= (apop_fn_with_void) dist.log_likelihood;
	if (!est_params || est_params->starting_pt==NULL){
		x	= gsl_vector_alloc(betasize);
  		gsl_vector_set_all (x,  0.1);
	}
	else 	x   = apop_array_to_vector(est_params->starting_pt, betasize);
	minme.f		= (apop_fn_with_void) negshell;
	minme.df	= (apop_df_with_void) dnegshell;
	minme.fdf	= (apop_fdf_with_void) fdf_shell;
	minme.n		= betasize;
	minme.params	= data;
	negshell_model	= dist;
	gsl_multimin_fdfminimizer_set (s, &minme, x, 
            est_params ? est_params->step_size : 0.05, 
            est_params ? est_params->tolerance : 1e-3);
      	do { 	iter++;
		status 	= gsl_multimin_fdfminimizer_iterate(s);
		if (status) 	break; 
		status = gsl_multimin_test_gradient(s->gradient, 
                                            est_params ? est_params->tolerance : 1e-3);
        	if (!est_params || est_params->verbose){
	        	printf ("%5i %.5f  f()=%10.5f gradient=%.3f\n", iter, gsl_vector_get (s->x, 0),  s->f, gsl_vector_get(s->gradient,0));
		}
        	if (status == GSL_SUCCESS){
			est->status	= 0;
		   	if(!est_params || est_params->verbose)	printf ("Minimum found.\n");
		}
       	 }
	while (status == GSL_CONTINUE && iter < MAX_ITERATIONS_w_d);
	//while (status != GSL_SUCCESS && iter < MAX_ITERATIONS_w_d);
	if(iter==MAX_ITERATIONS_w_d) {
		est->status	= 1;
		if (!est_params || est_params->verbose) printf("No min!!\n");
	}
	//Clean up, copy results to output estimate.
	gsl_vector_memcpy(est->parameters->vector, s->x);
	gsl_multimin_fdfminimizer_free(s);
	if (est_params && est_params->starting_pt==NULL) 
		gsl_vector_free(x);
	est->log_likelihood	= dist.log_likelihood(est->parameters->vector, data);
	if (!est_params || est->estimation_params.uses.covariance) 
		apop_numerical_covariance_matrix(dist, est, data);
	return est;
}

static apop_estimate *	apop_maximum_likelihood_no_d(apop_data * data, 
			apop_model dist, apop_estimation_params * est_params){
    //, double *starting_pt, double step_size, double tolerance, int verbose){
int			            status,
			            iter 		= 0,
			            betasize	= dist.parameter_ct;
size_t 			        i;
gsl_multimin_function 	minme;
gsl_multimin_fminimizer *s;
gsl_vector 		        *x, *ss;
double			        size;
apop_estimate		    *est;
	if (betasize == -1)	{
        dist.parameter_ct   =
        betasize            = data->matrix->size2 - 1;
    }
    betasize    *=  (est_params ? est_params->params_per_column : 1);
	s	= gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, betasize);
	ss	= gsl_vector_alloc(betasize);
	x	= gsl_vector_alloc(betasize);
	prep_inventory_mle(est_params->uses);
	//est	= apop_estimate_alloc(data->size1, betasize, NULL, actual_uses);
	est	= apop_estimate_alloc(data, dist, est_params);
	est->status	= 1;	//assume failure until we score a success.
	if (est_params->starting_pt==NULL)
  		gsl_vector_set_all (x,  0);
	else
		x   = apop_array_to_vector(est_params->starting_pt, betasize);
  	gsl_vector_set_all (ss,  est_params->step_size);
	minme.f		    = (apop_fn_with_void) negshell;
	minme.n		    = betasize;
	minme.params	= data;
	negshell_model	= dist;
	gsl_multimin_fminimizer_set (s, &minme, x,  ss);
      	do { 	iter++;
		status 	= gsl_multimin_fminimizer_iterate(s);
		if (status) 	break; 
		size	= gsl_multimin_fminimizer_size(s);
	   	status 	= gsl_multimin_test_size (size, est_params->tolerance); 
		if(est_params->verbose){
			printf ("%5d ", iter);
			for (i = 0; i < betasize; i++) {
				printf ("%8.3e ", gsl_vector_get (s->x, i)); } 
			printf ("f()=%7.3f size=%.3f\n", s->fval, size);
       			if (status == GSL_SUCCESS) {
				est->status	= 0;
	   			if(est_params->verbose){
					printf ("Minimum found at:\n");
					printf ("%5d ", iter);
					for (i = 0; i < betasize; i++) {
						printf ("%8.3e ", gsl_vector_get (s->x, i)); } 
					printf ("f()=%7.3f size=%.3f\n", s->fval, size);
				}
			}
		}
      	} while (status == GSL_CONTINUE && iter < MAX_ITERATIONS);
	if (iter == MAX_ITERATIONS && est_params->verbose)
		printf("Minimization reached maximum number of iterations.");

	gsl_vector_memcpy(est->parameters->vector, s->x);
	gsl_multimin_fminimizer_free(s);
	est->log_likelihood	= dist.log_likelihood(est->parameters->vector, data);
	if (est->estimation_params.uses.covariance) 
		apop_numerical_covariance_matrix(dist, est, data);
	return est;
}


/** The maximum likelihood calculations

\param data	the data matrix
\param	dist	the \ref apop_model object: waring, probit, zipf, &amp;c.
\param params	an \ref apop_estimation_params structure, featuring:<br>
starting_pt:	an array of doubles suggesting a starting point. If NULL, use zero.<br>
an \ref apop_inventory: which will be pared down and folded into the output \ref apop_estimate<br>
step_size:	the initial step size.<br>
tolerance:	the precision the minimizer uses. Only vaguely related to the precision of the actual var.<br>
verbose:	Y'know.<br>
method:		The sum of a method and a gradient-handling rule.
\li 000: Nelder-Mead simplex (gradient handling rule is irrelevant)
\li 100: conjugate gradient (Fletcher-Reeves) (default)
\li 200: conjugate gradient (BFGS: Broyden-Fletcher-Goldfarb-Shanno)
\li 300: conjugate gradient (Polak-Ribiere)
\li 500: \ref simanneal "simulated annealing"
\li 0: If no gradient is available, use numerical approximations. (default)
\li 1: Use numerical approximations even if an explicit dlog likelihood is given. <br>
Thus, the default method is 100+0 = 100. To use the Nelder-Mead simplex algorithm, use 0, or to use the Polak_Ribiere method ignoring any analytic dlog likelihood function, use 201.
\return	an \ref apop_estimate with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.

 \ingroup mle */
apop_estimate *	apop_maximum_likelihood(apop_data * data, apop_model dist, apop_estimation_params *params){
	if (params && params->method/100==5)
        return apop_annealing(dist, data, params);  //below.
    if (params && params->method/100==0)
		return apop_maximum_likelihood_no_d(data, dist, params);
	//else:
	return apop_maximum_likelihood_w_d(data, dist, params);
}


/** This function goes row by row through <tt>m</tt> and calculates the
likelihood of the given row, putting the result in <tt>v</tt>. 
You can use this to find the variance of the estimator if other means fail.

\param m 	A GSL matrix, exactly like those used for probit, Waring, Gamma, &c MLEs.

\param v	A vector which will hold the likelihood of each row of m. Declare but do not allocate.

\param dist	An \ref apop_model object whose log likelihood function you'd like to use.

\param fn_beta		The parameters at which you will evaluate the likelihood. If <tt>e</tt> is an \ref
			apop_estimate, then one could use <tt>e->parameters->vector</tt>.

This functions is used in the sample code in the \ref mle section.

\ingroup mle */
void apop_make_likelihood_vector(gsl_matrix *m, gsl_vector **v, 
				apop_model dist, gsl_vector* fn_beta){
gsl_matrix      mm;
int             i;
	*v	= gsl_vector_alloc(m->size1);
        for(i=0; i< m->size1; i++){
                mm      = gsl_matrix_submatrix(m,i,0, 1,m->size2).matrix;      //get a single row
                gsl_vector_set(*v, i, dist.log_likelihood(fn_beta, apop_matrix_to_data(&mm)));
        }
}

/** Input an earlier estimate, and then I will re-start the MLE search
 where the last one ended. You can specify greater precision or a new
 search method.

The prior estimation parameters are copied over.  If the estimate
converged to an OK value, then use the converged value; else use the
starting point from the last estimate.

Only one estimate is returned, either the one you sent in or a new
one. The loser (which may be the one you sent in) is freed. That is,
there is no memory leak when you do
\code
est = apop_estimate_restart(est, 200, 1e-2);
\endcode

 \param e   An \ref apop_estimate that is the output from a prior MLE estimation.
 \param new_method  If -1, use the prior method; otherwise, a new method to try.
 \param scale       Something like 1e-2. The step size and tolerance
                    will both be mutliplied by this amount. Of course, if this is 1, nothing changes.

\return         At the end of this procedure, we'll have two \ref
    apop_estimate structs: the one you sent in, and the one produced using the
    new method/scale. If the new estimate includes any NaNs/Infs, then
    the old estimate is returned (even if the old estimate included
    NaNs/Infs). Otherwise, the estimate with the largest log likelihood
    is returned.

\ingroup mle
\todo The tolerance for testing boundaries are hard coded (1e4). Will need to either add another input term or a global var.
*/ 
apop_estimate * apop_estimate_restart(apop_estimate *e,int  new_method, int scale){
double                 *start_pt2;
apop_estimate          *out;
apop_estimation_params new_params;
            //copy off the old params; modify the starting pt, method, and scale
    memcpy(&new_params, &(e->estimation_params),sizeof(apop_estimation_params));
    if (apop_vector_bounded(e->parameters->vector,1e4)){
        apop_vector_to_array(e->parameters->vector, &start_pt2);
	    new_params.starting_pt	= start_pt2;
    }
    else
	    new_params.starting_pt	= e->estimation_params.starting_pt;
    new_params.tolerance   = e->estimation_params.tolerance * scale;
    new_params.step_size   = e->estimation_params.step_size * scale;
    new_params.method	    = new_method;
	out	                    = (e->model)->estimate(e->data, &new_params);
            //Now check whether the new output is better than the old
    if (apop_vector_bounded(out->parameters->vector, 1e4) && out->log_likelihood > e->log_likelihood){
        apop_estimate_free(out);
        return e;
    } //else:
    apop_estimate_free(e);
    return out;
}




//////////////////////////
// Simulated Annealing.

/** \page simanneal Notes on simulated annealing

Simulated annealing is a controlled random walk.
As with the other methods, the system tries a new point, and if it
is better, switches. Initially, the system is allowed to make large
jumps, and then with each iteration, the jumps get smaller, eventually
converging. Also, there is some decreasing probability that if the new
point is {\em less} likely, it will still be chosen. Simulated annealing
is best for situations where there may be multiple local optima. Early
in the random walk, the system can readily jump from one to another;
later it will fine-tune its way toward the optimum. The number of points
tested is basically not dependent on the function: if you give it a
4,000 step program, that is basically how many steps it will take.
If you know your function is globally convex (as are most standard
probability functions), then this method is overkill.

The GSL's simulated annealing system doesn't actually do very much. It
basically provides a for loop that calls a half-dozen functions that we
the users get to write. So, the file \ref apop_mle.c handles all of this
for you. The likelihood function is taken from the model, the metric
is the Manhattan metric, the copy/destroy functions are just the usual
vector-handling fns., et cetera. The reader who wants further control
is welcome to override these functions.

 \ingroup mle
 */


/* set up parameters for this simulated annealing run. Cut 'n' pasted
 * from GSL traveling salesman sample code */

#define N_TRIES 200             /* how many points do we try before stepping */
#define ITERS_FIXED_T 200      /* how many iterations for each T? */
#define K 1.0                   /* Boltzmann constant */
#define MU_T 1.002              /* damping factor for temperature */
#define T_MIN 5.0e-1

static apop_model   annealing_model;
static apop_data    *annealing_data;
//static negshell_params    anneal_nsp;

static double annealing_energy(void *beta)
{
  //return -annealing_model.log_likelihood(beta, annealing_data->matrix);
  return negshell(beta, annealing_data);
}

/** We use the Manhattan metric to correspond to the annealing_step fn below.
 */
static double annealing_distance(void *xp, void *yp)
{
  return apop_vector_grid_distance(xp, yp);
}


/** The algorithm: 
    --randomly pick dimension
    --shift by some amount of remaining step size
    --repeat for all dims
This will give a move \f$\leq\f$ step_size on the Manhattan metric.
*/
static void annealing_step(const gsl_rng * r, void *xp, double step_size){

gsl_vector  *beta       = xp,
            *original   = gsl_vector_alloc(beta->size),
            *dummy      = gsl_vector_alloc(beta->size);
int         dims_used[beta->size];
int         dims_left, dim, sign, first_run = 0;
double      step_left, amt;
    gsl_vector_memcpy(original, beta);
    do{
        memset(dims_used, 0, beta->size * sizeof(int));
        dims_left   = beta->size;
        step_left   = step_size;
        if (!first_run)
            gsl_vector_memcpy(beta, original);
        else
            first_run   ++; 
        while (dims_left){
            do {
                dim = apop_random_int(0, beta->size, r);
            } while (dims_used[dim]);
            dims_used[dim]  ++;
            dims_left       --;
            sign    = (gsl_rng_uniform(r) > 0.5) ? 1 : -1;
            amt     = gsl_rng_uniform(r);
//printf("ss: %g, amt: %g, vector: ", step_left, amt * step_left * sign);
            apop_vector_increment(beta, dim, amt * step_left * sign); 
            step_left   *= amt;
//apop_vector_show(beta);
//printf("x");
        }
    } while (negshell_model.constraint && negshell_model.constraint(beta, annealing_data, dummy));
    gsl_vector_free(dummy);
    gsl_vector_free(original);
}

static void annealing_print(void *xp)
{
    apop_vector_show(xp);
}

static void annealing_memcpy(void *xp, void *yp){
    gsl_vector_memcpy(yp, xp);
}

static void *annealing_copy(void *xp){
    return apop_vector_copy(xp);
}

static void annealing_free(void *xp){
    gsl_vector_free(xp);
}

apop_estimate * apop_annealing(apop_model m, apop_data *data, apop_estimation_params *ep){
    //annealing_model = m;
    //annealing_data  = data;
    if (m.parameter_ct == -1) {
        m.parameter_ct   = data->matrix->size2 - 1;
    }
apop_estimate   *est  = apop_estimate_alloc(data, m, ep);
const gsl_rng   * r ;
gsl_vector      *beta;
    r   = gsl_rng_alloc(gsl_rng_env_setup()) ; 
    if (ep && ep->starting_pt)
        beta = apop_array_to_vector(ep->starting_pt, m.parameter_ct* (ep ? ep->params_per_column : 1));
    else{
        beta  = gsl_vector_alloc(m.parameter_ct * (ep ? ep->params_per_column : 1));
        gsl_vector_set_all(beta, 1);
    }
	annealing_data      = data;
	negshell_model	    = m;

gsl_siman_print_t printing_fn   = NULL;
    if (ep && ep->verbose)
        printing_fn = annealing_print;

gsl_siman_params_t params = {N_TRIES, 
                    ITERS_FIXED_T, 
                    (ep) ? ep->step_size : 1,
                    K, 
                    -m.log_likelihood(beta, data), 
                    MU_T, 
                    T_MIN};

    gsl_siman_solve(r,        //   const gsl_rng * r
          beta,          //   void * x0_p
          annealing_energy, //   gsl_siman_Efunc_t Ef
          annealing_step,   //   gsl_siman_step_t take_step
          annealing_distance, // gsl_siman_metric_t distance
          printing_fn,      //gsl_siman_print_t print_position
          annealing_memcpy, //   gsl_siman_copy_t copyfunc
          annealing_copy,   //   gsl_siman_copy_construct_t copy_constructor
          annealing_free,   //   gsl_siman_destroy_t destructor
         m.parameter_ct,    //   size_t element_size
         params);           //   gsl_siman_params_t params

    //Clean up, copy results to output estimate.
    gsl_vector_memcpy(est->parameters->vector, beta);
    est->log_likelihood = m.log_likelihood(est->parameters->vector, data);
    if (est->estimation_params.uses.covariance)
        apop_numerical_covariance_matrix(m, est, data);
    return est;
}
