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
#include "apop_findzeros.c"

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

typedef	void 	(*apop_df_with_void)(const gsl_vector *beta, void *d, gsl_vector *gradient);
typedef	void 	(*apop_fdf_with_void)( const gsl_vector *beta, void *d, double *f, gsl_vector *df);


//Including numerical differentiation and a couple of functions to
//negate the likelihood fns without bothering the user.


typedef struct {
	gsl_vector	*beta;
	apop_data	*d;
	int		    dimension;
} grad_params;

typedef struct {
    apop_model  *model;
    apop_data   *data;
    apop_fn_with_void   *f;
    apop_df_with_void   *df;
    void        *params;
    grad_params *gp;
    gsl_vector  *beta;
    int         use_constraint;
}   infostruct;

static apop_estimate * apop_annealing(infostruct*);  //below.

static double one_d(double b, void *in){
  infostruct    *i   =in;
    gsl_vector_set(i->gp->beta, i->gp->dimension, b);
    apop_data   *p  = apop_data_unpack(i->gp->beta, i->model->vsize, i->model->msize1, i->model->msize2);
	double out= (*(i->f))(p, i->gp->d, i->params);
    apop_data_free(p);
    return out;
}

/**The GSL provides one-dimensional numerical differentiation; here's
 the multidimensional extension.
 
 \code
 apop_numerical_gradient(your_model.log_likelihood, beta, params, out);
 \endcode

 \ingroup linear_algebra
 \todo This fn has a hard-coded tolerance (1e-5).
 */
void apop_numerical_gradient(apop_fn_with_void ll, const gsl_vector *beta, apop_data* d, void *v , gsl_vector *out){
  int		    j;
  gsl_function	F;
  double		result, err;
  grad_params 	gp;
  infostruct    i;
    i.params    = v;
    i.f         = &ll;
	gp.beta		= gsl_vector_alloc(beta->size);
	gp.d		= d;
    i.gp        = &gp;
	F.function	= one_d;
	F.params	= &i;
	for (j=0; j< beta->size; j++){
		gp.dimension	= j;
		gsl_vector_memcpy(gp.beta, beta);
		gsl_deriv_central(&F, gsl_vector_get(beta,j), 1e-5, &result, &err);
		gsl_vector_set(out, j, result);
	}
}


/** The MLEs don't return everything you could ask for, so this pares
down the input inventory to a modified output mle.

\todo This should really just modify the input inventory. */
static void prep_inventory_mle(apop_ep * in){
//These are the rules going from what you can ask for to what you'll get.
	in->uses.log_likelihood	= 1;
	in->uses.parameters		= 1;
	//OK, some things are not yet implemented.
	in->uses.names		    = 0;
	in->uses.confidence		= 0;
	//OK, some things are not yet implemented.
	in->uses.dependent		= 
	in->uses.predicted		= 0;
	if (in->uses.confidence==1)
		in->uses.covariance = 1;
	if (in->uses.covariance==1)
		in->uses.confidence = 1;
}


/* They always tell you to just negate your likelihood function to turn
a minimization routine into a maximization routine---and this is the
sort of annoying little detail that Apophenia is intended to take care
of for you. The next few functions do the negation, so you have one
less sign that you have to remember. 

These fns also take care of checking constraints, if any.

*/


static void insert_path_into_db(const gsl_vector *beta, double out){
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

static double negshell (const gsl_vector * betain, void * in){
  infostruct    *i              = in;
  apop_data 	*returned_beta	= NULL;
  double		penalty         = 0,
                out             = 0; 
  static double	base_for_penalty= 0;
  apop_data     *beta           = apop_data_unpack(betain, i->model->vsize, i->model->msize1, i->model->msize2);
	if (i->use_constraint && i->model->constraint){
        returned_beta	= apop_data_alloc(i->model->vsize, i->model->msize1, i->model->msize2);
		penalty	= i->model->constraint(beta, i->data, returned_beta, i->params);
		if (penalty > 0){
			base_for_penalty	= i->model->log_likelihood(returned_beta, i->data, i->params);
			out = -base_for_penalty + penalty;
		}
	}
    if (!penalty){
	    out = - i->model->log_likelihood(beta, i->data, i->params);
    }
    if (strlen(apop_opts.mle_trace_path))
        insert_path_into_db(betain,-out);
	if (returned_beta) apop_data_free(returned_beta);
    return out;
}


/* The derivative-calculating routine.
If the constraint binds
    then: take the numerical derivative of negshell, which will be the
    numerical derivative of the penalty.
    else: just find dlog_likelihood. If the model doesn't have a
    dlog likelihood or the user asked to ignore it, then the main
    maximum likelihood fn replaced model.score with
    apop_numerical_gradient anyway.
Finally, reverse the sign, since the GSL is trying to minimize instead of maximize.
*/
//static void dnegshell (const gsl_vector * beta, apop_data * d, gsl_vector * g){
static void dnegshell (gsl_vector * betain, infostruct * i, gsl_vector * g){
  apop_data    *returned_beta  = apop_data_alloc(i->model->vsize, i->model->msize1, i->model->msize2);
  apop_data     *beta           = apop_data_unpack(betain, i->model->vsize, i->model->msize1, i->model->msize2);
  gsl_vector     *usemep          = betain;
  apop_data     *useme          = beta;
	if (i->model->constraint && i->model->constraint(beta, i->data, returned_beta, i->params)){
		usemep  = apop_data_pack(returned_beta);
        useme   = returned_beta;
    }
    if (i->model->score)
        i->model->score(useme, i->data, g, i->params);
    else
        apop_numerical_gradient(i->model->log_likelihood, usemep, i->data, i->params, g);
    gsl_vector_scale(g, -1);
	apop_data_free(returned_beta);
    apop_data_free(beta);
}

//static void fdf_shell(const gsl_vector *beta, apop_data *d, double *f, gsl_vector *df){
static void fdf_shell(gsl_vector *beta, infostruct *i, double *f, gsl_vector *df){
	//if (negshell_model.fdf==NULL){
		*f	= negshell(beta, i);
		dnegshell(beta, i, df);
	/*} else	{
		negshell_model.fdf(beta, d, f, df);
		(*f) 	*= -1;
		gsl_vector_scale(df, -1);
	}*/
}



			//////////////////////////////////////////
			//The max likelihood functions themselves. 
			//Mostly straight out of the GSL manual.
			//////////////////////////////////////////

/*
static size_t	current_dimension;

	//Turn the full vector of derivatives into a single number.
	//Is this ever inefficient for large parameter spaces.
double gradient_for_hessian(const gsl_vector *beta, apop_data * d){
//double 		(*tmp) (const gsl_vector *beta, void *d);
apop_fn_with_void 		tmp;
gsl_vector	*waste	= gsl_vector_alloc(beta->size);
double		return_me;

	if (dlikelihood_for_hess == gradient_shell){	//then we need a numerical derivative.
		tmp			            = (apop_fn_with_void) apop_fn_for_derivative;
		apop_fn_for_derivative 	= (apop_fn_with_void) likelihood_for_hess;
		apop_numerical_gradient(beta, d , waste);
		apop_fn_for_derivative 	= tmp;
	}
	else 	dlikelihood_for_hess(beta, d , waste, model_params);
	return_me	= gsl_vector_get(waste, current_dimension);
	gsl_vector_free(waste);
	return return_me;
}
*/

/** Calculates the matrix of second derivatives of the log likelihood function.

\todo This is inefficient, by an order of n, where n is the number of parameters.
\ingroup linear_algebra
gsl_matrix * apop_numerical_second_derivative(apop_model dist, gsl_vector *beta, apop_data * d, void *params){
gsl_vector_view	v;
gsl_matrix	*out	= gsl_matrix_alloc(beta->size, beta->size);
	likelihood_for_hess	= dist.log_likelihood;
	if (dist.score)
		dlikelihood_for_hess	= dist.score;
	else	dlikelihood_for_hess	= gradient_shell;
	apop_fn_for_derivative	= (apop_fn_with_void) gradient_for_hessian;
	for (current_dimension=0; current_dimension< beta->size; current_dimension++){
		v	= gsl_matrix_row(out, current_dimension);
		apop_numerical_gradient(beta, d , &(v.vector));
	}
	return out;
}
 */

/** Calculate the Hessian.

  This is a synonym for \ref apop_numerical_second_derivative, q.v.
gsl_matrix * apop_numerical_hessian(apop_model dist, gsl_vector *beta, apop_data * d, void *params){
	return apop_numerical_second_derivative(dist, beta, d, params);
}
*/

/** Feeling lazy? Rather than doing actual pencil-and-paper math to find
your variance-covariance matrix, just use the negative inverse of the Hessian.

\param dist	The model
\param est	The estimate, with the parameters already calculated. The var/covar matrix will be placed in est->covariance.
\param data	The data
\ingroup basic_stats
void apop_numerical_covariance_matrix(apop_model dist, apop_estimate *est, apop_data *data){
//int		i;
apop_fn_with_void tmp;
gsl_matrix	    *hessian;
	tmp			= apop_fn_for_derivative;
	hessian			= apop_numerical_second_derivative(dist, est->parameters->vector, data);
	gsl_matrix_scale(hessian, -1);
	apop_det_and_inv(hessian, &(est->covariance->matrix), 0, 1);
	gsl_matrix_free(hessian);

	if (est->ep.uses.confidence == 0)
		return;
	//else:
        apop_estimate_parameter_t_tests(est);
	apop_fn_for_derivative 	= tmp;
}
*/
void apop_numerical_covariance_matrix(apop_model dist, apop_estimate *est, apop_data *data){
    //As you can see, this is a placeholder.
    return;
}


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
	dist.score(est->parameters->vector, data, diff);
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

	if (est->ep.uses.confidence == 0)
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
			apop_model dist, infostruct *i){
  gsl_multimin_function_fdf minme;
  gsl_multimin_fdfminimizer *s;
  gsl_vector 			    *x;
  int				        iter 	= 0, 
				            status  = 0,
				            betasize= dist.vsize + dist.msize1* dist.msize2;
  apop_estimate			    *est;
  apop_ep                   *ep     = i->params;
	if (betasize == 0 || betasize == -1)	{
        dist.vsize   =
        betasize     = data->matrix->size2 - betasize;
    }
	prep_inventory_mle(ep);
	est	= apop_estimate_alloc(data, dist, ep);
    if (!ep || ep->method/100 ==2 || ep->method == 2)
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, betasize);
    else if (ep->method/100 == 3 || ep->method == 3)
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_pr, betasize);
    else
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, betasize);
	/*if (dist.score == NULL || (ep && (ep->method % 100 ==1 || ep->method == 1)))
		i->df 	= gradient_shell;
    else
		i->df 	= dist.score;*/
	//i->f    = dist.log_likelihood;
	if (!ep || ep->starting_pt==NULL){
		x	= gsl_vector_alloc(betasize);
  		gsl_vector_set_all (x,  0.1);
	}
	else 	x   = apop_array_to_vector(ep->starting_pt, betasize);
	minme.f		= negshell;
	minme.df	= (apop_df_with_void) dnegshell;
	minme.fdf	= (apop_fdf_with_void) fdf_shell;
	minme.n		= betasize;
	minme.params	= i;
	gsl_multimin_fdfminimizer_set (s, &minme, x, 
            ep ? ep->step_size : 0.05, 
            ep ? ep->tolerance : 1e-3);
      	do { 	iter++;
		status 	= gsl_multimin_fdfminimizer_iterate(s);
		if (status) 	break; 
		status = gsl_multimin_test_gradient(s->gradient, 
                                            ep ? ep->tolerance : 1e-3);
        	if (!ep || ep->verbose){
	        	printf ("%5i %.5f  f()=%10.5f gradient=%.3f\n", iter, gsl_vector_get (s->x, 0),  s->f, gsl_vector_get(s->gradient,0));
		}
        	if (status == GSL_SUCCESS){
			est->status	= 0;
		   	if(!ep || ep->verbose)	printf ("Minimum found.\n");
		}
       	 }
	while (status == GSL_CONTINUE && iter < MAX_ITERATIONS_w_d);
	//while (status != GSL_SUCCESS && iter < MAX_ITERATIONS_w_d);
	if(iter==MAX_ITERATIONS_w_d) {
		est->status	= 1;
		if (!ep || ep->verbose) printf("No min!!\n");
	}
	//Clean up, copy results to output estimate.
    est->parameters = apop_data_unpack(s->x, dist.vsize, dist.msize1, dist.msize2);
	gsl_multimin_fdfminimizer_free(s);
	if (ep && ep->starting_pt==NULL) 
		gsl_vector_free(x);
	est->log_likelihood	= dist.log_likelihood(est->parameters, data, ep);
	if (!ep || est->ep.uses.covariance) 
		apop_numerical_covariance_matrix(dist, est, data);
	return est;
}

static apop_estimate *	apop_maximum_likelihood_no_d(apop_data * data, 
			apop_model dist, infostruct * i){
    //, double *starting_pt, double step_size, double tolerance, int verbose){
  int			            status,
			                iter 		= 0,
			                betasize	= dist.vsize + dist.msize1* dist.msize2;
  size_t 			        j;
  gsl_multimin_function 	minme;
  gsl_multimin_fminimizer   *s;
  gsl_vector 		        *x, *ss;
  double			        size;
  apop_estimate		        *est;
  apop_ep                   *ep         = i->params; 
	if (betasize == 0 || betasize == -1)	{
        dist.vsize   =
        betasize     = data->matrix->size2 - betasize;
    }
	s	= gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, betasize);
	ss	= gsl_vector_alloc(betasize);
	x	= gsl_vector_alloc(betasize);
	prep_inventory_mle(ep);
	//est	= apop_estimate_alloc(data->size1, betasize, NULL, actual_uses);
	est	= apop_estimate_alloc(data, dist, ep);
	est->status	= 1;	//assume failure until we score a success.
	if (ep->starting_pt==NULL)
  		gsl_vector_set_all (x,  0);
	else
		x   = apop_array_to_vector(ep->starting_pt, betasize);
  	gsl_vector_set_all (ss,  ep->step_size);
	minme.f		    = negshell;
	minme.n		    = betasize;
	minme.params	= i;
	gsl_multimin_fminimizer_set (s, &minme, x,  ss);
      	do { 	iter++;
		status 	= gsl_multimin_fminimizer_iterate(s);
		if (status) 	break; 
		size	= gsl_multimin_fminimizer_size(s);
	   	status 	= gsl_multimin_test_size (size, ep->tolerance); 
		if(ep->verbose){
			printf ("%5d ", iter);
			for (j = 0; j < betasize; j++) {
				printf ("%8.3e ", gsl_vector_get (s->x, j)); } 
			printf ("f()=%7.3f size=%.3f\n", s->fval, size);
       			if (status == GSL_SUCCESS) {
				est->status	= 0;
	   			if(ep->verbose){
					printf ("Minimum found at:\n");
					printf ("%5d ", iter);
					for (j = 0; j < betasize; j++) {
						printf ("%8.3e ", gsl_vector_get (s->x, j)); } 
					printf ("f()=%7.3f size=%.3f\n", s->fval, size);
				}
			}
		}
      	} while (status == GSL_CONTINUE && iter < MAX_ITERATIONS);
	if (iter == MAX_ITERATIONS && ep->verbose)
		printf("Minimization reached maximum number of iterations.");

    est->parameters = apop_data_unpack(s->x, dist.vsize, dist.msize1, dist.msize2);
	gsl_multimin_fminimizer_free(s);
	est->log_likelihood	= dist.log_likelihood(est->parameters, data, i->params);
	if (est->ep.uses.covariance) 
		apop_numerical_covariance_matrix(dist, est, data);
	return est;
}


/** The maximum likelihood calculations

\param data	the data matrix
\param	dist	the \ref apop_model object: waring, probit, zipf, &amp;c.
\param params	an \ref apop_ep structure, featuring:<br>
starting_pt:	an array of doubles suggesting a starting point. If NULL, use zero.<br>
step_size:	the initial step size.<br>
tolerance:	the precision the minimizer uses. Only vaguely related to the precision of the actual var.<br>
verbose:	Y'know.<br>
method:		The sum of a method and a gradient-handling rule.
\li 0: Nelder-Mead simplex (gradient handling rule is irrelevant)
\li 1: conjugate gradient (Fletcher-Reeves) (default)
\li 2: conjugate gradient (BFGS: Broyden-Fletcher-Goldfarb-Shanno)
\li 3: conjugate gradient (Polak-Ribiere)
\li 5: \ref simanneal "simulated annealing"
\li 10: Find a root of the derivative via Newton's method
\li 11: Find a root of the derivative via the Broyden Algorithm
\li 12: Find a root of the derivative via the Hybrid method
\li 13: Find a root of the derivative via the Hybrid method; no internal scaling
\return	an \ref apop_estimate with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.

 \ingroup mle */
apop_estimate *	apop_maximum_likelihood(apop_data * data, apop_model dist, apop_ep *params){
  infostruct    info    = {&dist, data, NULL, NULL, params, NULL, NULL, 1};
	if (params && (params->method/100==5 || params->method == 5))
        return apop_annealing(&info);  //below.
    else if (params && params->method==0)
		return apop_maximum_likelihood_no_d(data, dist, &info);
    else if (params && params->method >= 10 && params->method <= 13)
        return  find_roots (data, dist, params);
	//else:
	return apop_maximum_likelihood_w_d(data, dist, &info);
}


/** This function goes row by row through <tt>m</tt> and calculates the
likelihood of the given row, putting the result in <tt>v</tt>. 
You can use this to find the variance of the estimator if other means fail.

\param m 	A GSL matrix, exactly like those used for probit, Waring, Gamma, &c MLEs.

\param v	A vector that will hold the likelihood of each row of m. Declare but do not allocate.

\param dist	An \ref apop_model object whose log likelihood function you'd like to use.

\param fn_beta		The parameters at which you will evaluate the likelihood. If <tt>e</tt> is an \ref
			apop_estimate, then one could use <tt>e->parameters->vector</tt>.

This functions is used in the sample code in the \ref mle section.

void apop_make_likelihood_vector(gsl_matrix *m, gsl_vector **v, 
				apop_model dist, gsl_vector* fn_beta){
gsl_matrix      mm;
int             i;
	*v	= gsl_vector_alloc(m->size1);
        for(i=0; i< m->size1; i++){
                mm      = gsl_matrix_submatrix(m,i,0, 1,m->size2).matrix;      //get a single row
                gsl_vector_set(*v, i, dist.log_likelihood(fn_beta, apop_matrix_to_data(&mm), model_params));
        }
}
\ingroup mle */

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
apop_ep new_params;
            //copy off the old params; modify the starting pt, method, and scale
    memcpy(&new_params, &(e->ep),sizeof(apop_ep));
    if (apop_vector_bounded(e->parameters->vector,1e4)){
        apop_vector_to_array(e->parameters->vector, &start_pt2);
	    new_params.starting_pt	= start_pt2;
    }
    else
	    new_params.starting_pt	= e->ep.starting_pt;
    new_params.tolerance   = e->ep.tolerance * scale;
    new_params.step_size   = e->ep.step_size * scale;
    new_params.method	    = new_method;
	out	                    = (e->model)->estimate(e->data, &new_params);
            //Now check whether the new output is better than the old
    if (apop_vector_bounded(out->parameters->vector, 1e4) && out->log_likelihood > e->log_likelihood){
        apop_estimate_free(e);
        return out;
    } //else:
    apop_estimate_free(out);
    return e;
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

Verbosity: if ep->verbose==1, show likelihood,  temp, &c. in a table;
if ep->verbose>1, show that plus the vector of params.

 \ingroup mle
 */


/* set up parameters for this simulated annealing run. Cut 'n' pasted
 * from GSL traveling salesman sample code */

#define N_TRIES 200             /* how many points do we try before stepping */
#define ITERS_FIXED_T 200      /* how many iterations for each T? */
#define K 1.0                   /* Boltzmann constant */
#define MU_T 1.005              /* damping factor for temperature */
#define T_MIN 5.0e-1

static double annealing_energy(void *in)
{
  infostruct *i = in;
    return negshell(((infostruct*)in)->beta, i);
}

/** We use the Manhattan metric to correspond to the annealing_step fn below.
 */
static double annealing_distance(void *xin, void *yin)
{
  return apop_vector_grid_distance(((infostruct*)xin)->beta, ((infostruct*)yin)->beta);
}


/** The algorithm: 
    --randomly pick dimension
    --shift by some amount of remaining step size
    --repeat for all dims
This will give a move \f$\leq\f$ step_size on the Manhattan metric.
*/
static void annealing_step(const gsl_rng * r, void *in, double step_size){
  infostruct  *i          = in;
  gsl_vector  *original   = gsl_vector_alloc(i->beta->size);
  int         dims_used[i->beta->size];
  int         dims_left, dim, sign, first_run = 0;
  double      step_left, amt;
  apop_data     *testme     = NULL,
                *dummy      = apop_data_alloc(i->model->vsize, i->model->msize1, i->model->msize2);
    gsl_vector_memcpy(original, i->beta);
    do{
        memset(dims_used, 0, i->beta->size * sizeof(int));
        dims_left   = i->beta->size;
        step_left   = step_size;
        if (!first_run)
            gsl_vector_memcpy(i->beta, original);
        else
            first_run   ++; 
        while (dims_left){
            do {
                dim = apop_random_int(0, i->beta->size, r);
            } while (dims_used[dim]);
            dims_used[dim]  ++;
            dims_left       --;
            sign    = (gsl_rng_uniform(r) > 0.5) ? 1 : -1;
            amt     = gsl_rng_uniform(r);
//printf("ss: %g, amt: %g, vector: ", step_left, amt * step_left * sign);
            apop_vector_increment(i->beta, dim, amt * step_left * sign); 
            step_left   *= amt;
//apop_vector_show(beta);
//printf("x");
            if (testme) apop_data_free(testme);
            testme      = apop_data_unpack(i->beta, i->model->vsize, i->model->msize1, i->model->msize2);
        }
    } while (i->model->constraint && i->model->constraint(testme, i->data, dummy, i->params));
    apop_data_free(dummy);
    apop_data_free(testme);
    gsl_vector_free(original);
}

static void annealing_print2(void *xp) { return; }

static void annealing_print(void *xp) {
    apop_vector_show(((infostruct*)xp)->beta);
}

static void annealing_memcpy(void *xp, void *yp){
    memcpy(yp, xp, sizeof(infostruct));
    ((infostruct*)yp)->beta = gsl_vector_alloc(((infostruct*)xp)->beta->size);
    gsl_vector_memcpy(((infostruct*)yp)->beta, ((infostruct*)xp)->beta);
}

static void *annealing_copy(void *xp){
    infostruct *out = malloc(sizeof(infostruct));
    memcpy(out, xp, sizeof(infostruct));
    out->beta        = gsl_vector_alloc(((infostruct*)xp)->beta->size);
    gsl_vector_memcpy(out->beta, ((infostruct*)xp)->beta);
    return out;
}

static void annealing_free(void *xp){
    gsl_vector_free(((infostruct*)xp)->beta);
    free(xp);
}

apop_estimate * apop_annealing(infostruct *i){
  apop_data     *data   = i->data;
  apop_model     *m     = i->model;
  int           paramct = i->model->vsize + i->model->msize1* i->model->msize2;
    if (m->vsize == -1) {
        m->vsize   = i->data->matrix->size2 - 1;
    }
  apop_ep           *ep = i->params;
  apop_estimate *est    = apop_estimate_alloc(i->data, *(i->model), ep);
const gsl_rng   * r ;
gsl_vector      *beta;
    r   = gsl_rng_alloc(gsl_rng_env_setup()) ; 
    if (ep && ep->starting_pt)
        beta = apop_array_to_vector(ep->starting_pt, paramct);
    else{
        beta  = gsl_vector_alloc(paramct);
        gsl_vector_set_all(beta, 1);
    }
	i->beta             = beta;
    i->use_constraint   = 0; //negshell doesn't check it; annealing_step does.

gsl_siman_print_t printing_fn   = NULL;
    if (ep && ep->verbose>1)
        printing_fn = annealing_print;
    else if (ep && ep->verbose)
        printing_fn = annealing_print2;

gsl_siman_params_t params = {N_TRIES, 
                    ITERS_FIXED_T, 
                    (ep) ? ep->step_size : 1,
                    K, 
                    212.6, //-i->model->log_likelihood(beta, i->data, ep), 
                    MU_T, 
                    T_MIN};

    gsl_siman_solve(r,        //   const gsl_rng * r
          i,                //   void * x0_p
          annealing_energy, //   gsl_siman_Efunc_t Ef
          annealing_step,   //   gsl_siman_step_t take_step
          annealing_distance, // gsl_siman_metric_t distance
          printing_fn,      //gsl_siman_print_t print_position
          annealing_memcpy, //   gsl_siman_copy_t copyfunc
          annealing_copy,   //   gsl_siman_copy_construct_t copy_constructor
          annealing_free,   //   gsl_siman_destroy_t destructor
         paramct,    //   size_t element_size
         params);           //   gsl_siman_params_t params

    //Clean up, copy results to output estimate.
    gsl_vector_memcpy(est->parameters->vector, i->beta);
    est->log_likelihood = i->model->log_likelihood(est->parameters, data, ep);
    if (est->ep.uses.covariance)
        apop_numerical_covariance_matrix(*m,est, data);
    return est;
}



//The reader will recognize this as solving the linear equation.
static double set_one_constraint_elmt(gsl_vector *beta, apop_data *c, int constraint_no, int beta_no){
   double   denom   = gsl_matrix_get(c->matrix, constraint_no, beta_no);
   assert (denom);
  int       i;
  double    score   = gsl_vector_get(c->vector, constraint_no);
    for (i=0; i< beta->size; i++)
        if(i!=beta_no)
            score   -= gsl_matrix_get(c->matrix, constraint_no, i) * gsl_vector_get(beta, i);
    score   /= denom;
    return score;
}

/** Use this for producing models. 

  \param beta The parameter vector to be tested against the constraint.
  \param d   an apop_data set, where each row represents one constraint.
  \param returned_beta  The pre-allocated backup vector that will be modified if the given parameter vector is out of bounds.

  \bugs The tolerance is hard-coded, wich is _terrible_.

*/
double apop_linear_constraint(gsl_vector *beta, void * d, gsl_vector *returned_beta){
  double  tolerance     = 1e-3;
  double  penalty       = 0,
            p, goal;
  int       i, j, modded[beta->size];
  apop_data *c          = d;
  apop_data *v          = apop_data_alloc(0,0,0);
    v->vector = beta;
  apop_data *test  = apop_dot(c, v, 0);//c->matrix
    memset(modded, 0, sizeof(int)*beta->size);
    for (i=0; i< test->vector->size; i++){
        p   = gsl_vector_get(test->vector, i);
        goal= gsl_vector_get(c->vector, i);
        if (p> 0){
            penalty += p;
            j = i;
            while (modded[j] || !gsl_matrix_get(c->matrix, i, j)){
                j   = (j+1)% beta->size;
                if (j==i)   //then _everything_ has been modded. Oh my.
                    memset(modded, 0, sizeof(int)*beta->size);
            }
            gsl_vector_set(returned_beta, j, set_one_constraint_elmt(beta,c,i, j) + tolerance);
            modded[j]   ++;
        }
    }
    apop_data_free(test);
    return penalty;
}
