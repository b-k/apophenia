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

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL.
*/
#include "likelihoods.h"
#include <assert.h>
#include <gsl/gsl_deriv.h>

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

//Including numerical differentiation and a couple of functions to
//negate the likelihood fns without bothering the user.

/** If you are taking a numerical gradient, you need to set this
 variable to the function for which the gradient will be taken. 
 \ingroup global_vars */
double (*apop_fn_for_derivative) (const gsl_vector *beta, void *d);


typedef struct grad_params{
	gsl_vector	*beta;
	void		*d;
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
 \todo This fn has a hard-coded tolerance (1e-6).
 */
void apop_numerical_gradient(const gsl_vector *beta, void *d , gsl_vector *out){
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
typedef struct negshell_params{
	void        *d;
	apop_model	model;
} negshell_params;

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

static double negshell (const gsl_vector * beta, void * params){
negshell_params	*p			    = params;
gsl_vector 	    *returned_beta	= gsl_vector_alloc(beta->size);
double		    penalty, out    = 0, 
                done            = 0; 
static double	base_for_penalty= 0;
	if (p->model.constraint){
		penalty	= p->model.constraint(beta, p->d, returned_beta);
		if (penalty > 0){
			base_for_penalty	= p->model.log_likelihood(returned_beta, p->d);
			out = -base_for_penalty + penalty;
            done++;
		}
	}
    if (!done)
	    out = - p->model.log_likelihood(beta, p->d);
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
static void dnegshell (const gsl_vector * beta, void * params, gsl_vector * g){
negshell_params *p		        = params;
gsl_vector      *returned_beta  = gsl_vector_alloc(beta->size);
double 		    (*tmp) (const gsl_vector *beta, void *d);
	if (p->model.constraint && p->model.constraint(beta, p->d, returned_beta)){
		    tmp			            = apop_fn_for_derivative;
		    apop_fn_for_derivative 	= negshell;
		    apop_numerical_gradient(beta, params, g);
		    apop_fn_for_derivative 	= tmp;
	} else {
	    p->model.dlog_likelihood(beta, p->d, g);
	    gsl_vector_scale(g, -1);
    }
	gsl_vector_free(returned_beta);
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

/** Calculates the matrix of second derivatives of the log likelihood function.

The information matrix is the inverse of the negation of this.

\todo This is inefficient, by an order of n, where n is the number of parameters.
\ingroup linear_algebra
 */
gsl_matrix * apop_numerical_second_derivative(apop_model dist, gsl_vector *beta, void * d){
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

/** Calculate the Hessian.

  This is a synonym for \ref apop_numerical_second_derivative, q.v.
*/
gsl_matrix * apop_numerical_hessian(apop_model dist, gsl_vector *beta, void * d){
	return apop_numerical_second_derivative(dist, beta, d);
}

/** Feeling lazy? Rather than doing actual pencil-and-paper math to find
your variance-covariance matrix, just use the negative inverse of the Hessian.

\param dist	The model
\param est	The estimate, with the parameters already calculated. The var/covar matrix will be placed in est->covariance.
\param data	The data
\ingroup basic_stats
*/
void apop_numerical_var_covar_matrix(apop_model dist, apop_estimate *est, gsl_matrix *data){
int		i;
double 		(*tmp) (const gsl_vector *beta, void *d);
gsl_matrix	*hessian;
	tmp			= apop_fn_for_derivative;
	//The information matrix is the inverse of the negation of the hessian.
	hessian			= apop_numerical_hessian(dist, est->parameters, data);
	gsl_matrix_scale(hessian, -1);
	apop_det_and_inv(hessian, &(est->covariance->data), 0, 1);
	gsl_matrix_free(hessian);

	//Confidence intervals are just a cute convenience. We're
	//assuming Normality, which only works asymptotically anyway.
	if (est->estimation_params.uses.confidence == 0)
		return;
	//else:
	for (i=0; i< est->parameters->size; i++) // confidence[i] = |1 - (1-N(Mu[i],sigma[i]))*2|
		gsl_vector_set(est->confidence, i,
			fabs(1 - (1 - gsl_cdf_gaussian_P(gsl_vector_get(est->parameters, i), 
			gsl_matrix_get(est->covariance->data, i, i)))*2));
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

	if (est->estimation_params.uses.confidence == 0)
		return;
	//else:
	for (i=0; i<betasize; i++) // confidence[i] = |1 - (1-N(Mu[i],sigma[i]))*2|
		gsl_vector_set(est->confidence, i,
			fabs(1 - (1 - gsl_cdf_gaussian_P(gsl_vector_get(est->parameters, i), 
			gsl_matrix_get(est->covariance, i, i)))*2));
}
*/


/* The maximum likelihood calculations, given a derivative of the log likelihood.

If no derivative exists, will calculate a numerical gradient.

\param data	the data matrix
\param	dist	the \ref apop_model object: waring, probit, zipf, &amp;c.
\param	starting_pt	an array of doubles suggesting a starting point. If NULL, use zero.
\param step_size	the initial step size.
\param tolerance	the precision the minimizer uses. Only vaguely related to the precision of the actual var.
\param verbose		Y'know.
\return	an \ref apop_estimate with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.

  \todo readd names */
static apop_estimate *	apop_maximum_likelihood_w_d(apop_data * data,
			apop_model dist, apop_estimation_params *est_params){
            //double *starting_pt, double step_size, double tolerance, int method, int verbose){
gsl_multimin_function_fdf 	minme;
gsl_multimin_fdfminimizer 	*s;
gsl_vector 			        *x;
int				            iter 	= 0, 
				            status  = 0,
				            betasize= dist.parameter_ct;
apop_estimate			    *est;
negshell_params			    nsp;
	if (betasize == -1)	{
        dist.parameter_ct   =
        betasize            = data->data->size2 - 1;
    }
	prep_inventory_mle(est_params->uses);
	est	= apop_estimate_alloc(data, dist, est_params);
    if (est_params->method/100 ==2)
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, betasize);
    else if (est_params->method/100 == 3)
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_pr, betasize);
    else
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, betasize);
	if (dist.dlog_likelihood == NULL || (est_params->method % 100 ==1))
		dist.dlog_likelihood 	= apop_numerical_gradient;
	apop_fn_for_derivative	= dist.log_likelihood;
	if (est_params->starting_pt==NULL){
		x	= gsl_vector_alloc(betasize);
  		gsl_vector_set_all (x,  0);
	}
	else 	x   = apop_array_to_vector(est_params->starting_pt, betasize);
	nsp.d		= data->data;
	nsp.model	= dist;
	minme.f		= negshell;
	minme.df	= dnegshell;
	minme.fdf	= fdf_shell;
	minme.n		= betasize;
	minme.params	= &nsp;
	gsl_multimin_fdfminimizer_set (s, &minme, x, est_params->step_size, est_params->tolerance);
      	do { 	iter++;
		status 	= gsl_multimin_fdfminimizer_iterate(s);
		if (status) 	break; 
		status = gsl_multimin_test_gradient(s->gradient, est_params->tolerance);
        	if (est_params->verbose){
	        	printf ("%5i %.5f  f()=%10.5f gradient=%.3f\n", iter, gsl_vector_get (s->x, 0),  s->f, gsl_vector_get(s->gradient,0));
		}
        	if (status == GSL_SUCCESS){
			est->status	= 0;
		   	if(est_params->verbose)	printf ("Minimum found.\n");
		}
       	 }
	while (status == GSL_CONTINUE && iter < MAX_ITERATIONS_w_d);
	//while (status != GSL_SUCCESS && iter < MAX_ITERATIONS_w_d);
	if(iter==MAX_ITERATIONS_w_d) {
		est->status	= 1;
		if (est_params->verbose) printf("No min!!\n");
	}
	//Clean up, copy results to output estimate.
	gsl_vector_memcpy(est->parameters, s->x);
	gsl_multimin_fdfminimizer_free(s);
	if (est_params->starting_pt==NULL) 
		gsl_vector_free(x);
	est->log_likelihood	= dist.log_likelihood(est->parameters, data->data);
	if (est->estimation_params.uses.covariance) 
		apop_numerical_var_covar_matrix(dist, est, data->data);
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
negshell_params		    nsp;
	if (betasize == -1)	{
        dist.parameter_ct   =
        betasize            = data->data->size2 - 1;
    }
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
	nsp.model	    = dist;
	nsp.d		    = data->data;
	minme.f		    = negshell;
	minme.n		    = betasize;
	minme.params	= &nsp;
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

	gsl_vector_memcpy(est->parameters, s->x);
	gsl_multimin_fminimizer_free(s);
	est->log_likelihood	= dist.log_likelihood(est->parameters, data->data);
	if (est->estimation_params.uses.covariance) 
		apop_numerical_var_covar_matrix(dist, est, data->data);
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
\li 0: If no gradient is available, use numerical approximations. (default)
\li 1: Use numerical approximations even if an explicit dlog likelihood is given. <br>
Thus, the default method is 100+0 = 100. To use the Nelder-Mead simplex algorithm, use 0, or to use the Polak_Ribiere method ignoring any analytic dlog likelihood function, use 201.
\return	an \ref apop_estimate with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.

  \todo re-add names 
 \ingroup mle */
apop_estimate *	apop_maximum_likelihood(apop_data * data, apop_model dist, apop_estimation_params *params){

	if (params->method/100==0)
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
    if (apop_vector_bounded(e->parameters,1e4)){
        apop_vector_to_array(e->parameters, &start_pt2);
	    new_params.starting_pt	= start_pt2;
    }
    else
	    new_params.starting_pt	= e->estimation_params.starting_pt;
    new_params.tolerance   = e->estimation_params.tolerance * scale;
    new_params.step_size   = e->estimation_params.step_size * scale;
    new_params.method	    = new_method;
	out	                    = (e->model)->estimate(e->data, &new_params);
            //Now check whether the new output is better than the old
    if (apop_vector_bounded(out->parameters, 1e4) && out->log_likelihood > e->log_likelihood){
        apop_estimate_free(out);
        return e;
    } //else:
    apop_estimate_free(e);
    return out;
}
