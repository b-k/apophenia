/** \file likelihoods.c	The MLE functions. Call them with an \ref apop_likelihood.

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

/*
Feed in data, the parameters (to be output), the # of parameters, and a
pointer to the three non-user functions above.  You'll get the most
likely betas back out.  The return value of the function itself is the
likelihood evaluated with the most likely betas.  */


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




////////////////////
//The max likelihood functions themselves. Mostly straight out of the GSL manual.
////////////////////



/** The maximum likelihood calculations, given a derivative of the log likelihood.

Don't bother calling this fn. Call \ref apop_maximum_likelihood, and if the derivative exists in the distribution, then this will be called.

\param data	the data matrix
\param uses	an inventory, which will be pared down and folded into the output \ref apop_estimate
\param	dist	the \ref apop_likelihood object: waring, probit, zipf, &amp;c.
\param	starting_pt	an array of doubles suggesting a starting point. If NULL, use zero.
\param step_size	the initial step size.
\param tolerance	the precision the minimizer uses. Only vaguely related to the precision of the actual var.
\param verbose		Y'know.
\return	an \ref apop_estimate with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.

  \todo readd names */
apop_estimate *	apop_maximum_likelihood_w_d(gsl_matrix * data, apop_inventory *uses,
			apop_likelihood dist, double *starting_pt, double step_size, double tolerance, int verbose){

	void fdf(const gsl_vector *beta, void *d, double *f, gsl_vector *df){
		*f	= dist.log_likelihood(beta, d);
		dist.dlog_likelihood(beta, d, df);
	}

gsl_multimin_function_fdf 	minme;
gsl_multimin_fdfminimizer 	*s;
gsl_vector 			*x, *diff;
gsl_vector_view 		v;
int				iter =0, status, i;
int				betasize	= dist.parameter_ct;
apop_inventory			actual_uses;
apop_estimate			*est;
	if (betasize == -1)	betasize = data->size2 - 1;
	prep_inventory_mle(uses, &actual_uses);
	s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, betasize);
	//s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, betasize);
	est	= apop_estimate_alloc(data->size1, betasize, NULL, actual_uses);
	if (starting_pt==NULL){
		x	= gsl_vector_alloc(betasize);
  		gsl_vector_set_all (x,  0);
	}
	else
		apop_convert_array_to_vector(starting_pt, &x, betasize);
	minme.f		= dist.log_likelihood;
	minme.df	= dist.dlog_likelihood;
	if (!dist.fdf)
		minme.fdf	= fdf;
	else	minme.fdf	= dist.fdf;
	minme.n		= betasize;
	minme.params	= data;
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


	gsl_vector_memcpy(est->parameters, s->x);
	gsl_multimin_fdfminimizer_free(s);
	if (starting_pt==NULL) gsl_vector_free(x);

	est->log_likelihood	= -dist.log_likelihood(est->parameters, data);
	//Epilogue:
	//find the variance-covariance matrix, using $df/d\theta \cdot df/d\theta$
	if (est->uses.covariance == 0) 
		return est;
	//else:
gsl_matrix	*pre_cov;
	pre_cov			= gsl_matrix_alloc(betasize, betasize);
	//estimate->covariance	= gsl_matrix_alloc(betasize, betasize);
	diff			= gsl_vector_alloc(betasize);
	dist.dlog_likelihood(est->parameters, data, diff);
	for (i=0; i< betasize; i++){
		gsl_matrix_set_row(pre_cov, i, diff);
		v	= gsl_matrix_row(pre_cov, i);
		gsl_vector_scale(&(v.vector), gsl_vector_get(diff, i));
	}
	apop_det_and_inv(pre_cov, &(est->covariance), 0,1);
	gsl_matrix_free(pre_cov);
	gsl_vector_free(diff);

	if (est->uses.confidence == 0)
		return est;
	//else:
	for (i=0; i<betasize; i++) // confidence[i] = |1 - (1-N(Mu[i],sigma[i]))*2|
		gsl_vector_set(est->confidence, i,
			fabs(1 - (1 - gsl_cdf_gaussian_P(gsl_vector_get(est->parameters, i), 
			gsl_matrix_get(est->covariance, i, i)))*2));
	return est;
}


/** The maximum likelihood calculations

\param data	the data matrix
\param uses	an inventory, which will be pared down and folded into the output \ref apop_estimate
\param	dist	the \ref apop_likelihood object: waring, probit, zipf, &amp;c.
\param	starting_pt	an array of doubles suggesting a starting point. If NULL, use zero.
\param step_size	the initial step size.
\param tolerance	the precision the minimizer uses. Only vaguely related to the precision of the actual var.
\param verbose		Y'know.
\return	an \ref apop_estimate with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.

  \todo re-add names 
 \ingroup mle */
apop_estimate *	apop_maximum_likelihood(gsl_matrix * data, apop_inventory *uses,
			apop_likelihood dist, double *starting_pt, double step_size, double tolerance, int verbose){

	//First, if the derivative is available, you're in the wrong place:
	//
	//Except this hasn't been tested enough for me to be happy with it.
	if (dist.dlog_likelihood!=NULL)
		return apop_maximum_likelihood_w_d(data, uses, dist, starting_pt, step_size, tolerance, verbose);

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

	minme.f		= dist.log_likelihood;
	minme.n		= betasize;
	minme.params	= data;
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
	est->log_likelihood	= -dist.log_likelihood(est->parameters, data);
	//if (est->uses.confidence == 0)
		return est;
}


/** This function goes row by row through <tt>m</tt> and calculates the
likelihood of the given row, putting the result in <tt>v</tt>. 
You can use this to find the variance of the estimator if other means fail.

\param m 	A GSL matrix, exactly like those used for probit, Waring, Gamma, &c MLEs.

\param v	A vector which will hold the likelihood of each row of m. Declare but do not allocate.

\param dist	An \ref apop_likelihood object whose log likelihood function you'd like to use.

\param fn_beta		The parameters at which you will evaluate the likelihood. If <tt>e</tt> is an \ref
			apop_estimate, then one could use <tt>e->parameters</tt>.

This functions is used in the sample code in the \ref mle section.

\ingroup mle */
void apop_make_likelihood_vector(gsl_matrix *m, gsl_vector **v, 
				apop_likelihood dist, gsl_vector* fn_beta){
gsl_matrix_view mm;
int             i;
	*v	= gsl_vector_alloc(m->size1);
        for(i=0; i< m->size1; i++){
                mm      = gsl_matrix_submatrix(m,i,0, 1,m->size2);      //get a single row
                gsl_vector_set(*v, i, dist.log_likelihood(fn_beta, &(mm.matrix)));
        }
}
