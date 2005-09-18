/** \file likelihoods.c	An easy front end to SQLite. Includes a few nice
features like a variance, skew, and kurtosis aggregator for SQL.

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
#include <apophenia/likelihoods.h>

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



/** The maximum likelihood calculations

\param data	the data matrix
\param uses	an inventory, which will be pared down and folded into the output \ref apop_estimate
\param	dist	the \ref apop_distribution object: waring, probit, zipf, &amp;c.
\param	starting_pt	an array of doubles suggesting a starting point. If NULL, use zero.
\param step_size	the initial step size.
\param verbose		Y'know.
\return	an \ref apop_estimate with the parameter estimates, &c.

  \todo readd names */
apop_estimate *	apop_maximum_likelihood_w_d(gsl_matrix * data, apop_inventory *uses,
			apop_distribution dist, double *starting_pt, double step_size, int verbose){

	void fdf(const gsl_vector *beta, void *d, double *f, gsl_vector *df){
		*f	= dist.log_likelihood(beta, d);
		dist.dlog_likelihood(beta, d, df);
	}

gsl_multimin_function_fdf 	minme;
gsl_multimin_fdfminimizer 	*s;
gsl_vector 			*x, *ss, *diff;
gsl_vector_view 		v;
int				iter =0, status, i;
int				betasize	= dist.parameter_ct;
apop_inventory			actual_uses;
apop_estimate			*est;
	prep_inventory_mle(uses, &actual_uses);
	s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, betasize);
	//s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, betasize);
	ss	= gsl_vector_alloc(betasize);
	est	= apop_estimate_alloc(data->size1, betasize, NULL, actual_uses);
	if (starting_pt==NULL){
		x	= gsl_vector_alloc(betasize);
  		gsl_vector_set_all (x,  0);
	}
	else
		apop_convert_array_to_vector(starting_pt, &x, betasize);
  	gsl_vector_set_all (ss,  step_size);
	minme.f		= dist.log_likelihood;
	minme.df	= dist.dlog_likelihood;
	if (!dist.fdf)
		minme.fdf	= fdf;
	else	minme.fdf	= dist.fdf;
	minme.n		= betasize;
	minme.params	= (void *)data;
	gsl_multimin_fdfminimizer_set (s, &minme, x, .001, 1e-5);

      	do { 	iter++;
		status 	= gsl_multimin_fdfminimizer_iterate(s);
//		if (status) 	break; 
        	if (verbose){
	        	printf ("%5i %.5f  f()=%10.5f gradient=%.3f\n", iter, gsl_vector_get (s->x, 0),  s->f, gsl_vector_get(s->gradient,0));
		}
		status = gsl_multimin_test_gradient(s->gradient, 1e-1);
	//	size	= gsl_multimin_fminimizer_size(s);
        //	status 	= gsl_multimin_test_size (size, 1e-5); 
        	if (verbose && status == GSL_SUCCESS){
		   		printf ("Minimum found.\n");
		}
       	 }
	while (status == GSL_CONTINUE && iter < MAX_ITERATIONS_w_d);
	if(iter==MAX_ITERATIONS_w_d) printf("No min!!\n");

	gsl_vector_memcpy(est->parameters, s->x);
	gsl_multimin_fdfminimizer_free(s);

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


/** The maximum likelihood calculations, if you don't have the
 derivative of the log likelihood fn on hand.

\param data	the data matrix
\param uses	an inventory, which will be pared down and folded into the output \ref apop_estimate
\param	dist	the \ref apop_distribution object: waring, probit, zipf, &amp;c.
\param	starting_pt	an array of doubles suggesting a starting point. If NULL, use zero.
\param step_size	the initial step size.
\param verbose		Y'know.
\return	an \ref apop_estimate with the parameter estimates, &c.

  \todo readd names */
apop_estimate *	apop_maximum_likelihood(gsl_matrix * data, apop_inventory *uses,
			apop_distribution dist, double *starting_pt, double step_size, int verbose){
int			iter =0, status;
gsl_multimin_function 	minme;
gsl_multimin_fminimizer *s;
gsl_vector 		*x, *ss;
double			size;
int				betasize	= dist.parameter_ct;
apop_inventory			actual_uses;
apop_estimate			*est;
	s	= gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, betasize);
	ss	= gsl_vector_alloc(betasize);
	prep_inventory_mle(uses, &actual_uses);
	est	= apop_estimate_alloc(data->size1, betasize, NULL, actual_uses);
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
        	status 	= gsl_multimin_test_size (size, 1e-3); 
		if(verbose){
		int i;
			printf ("%5d ", iter);
			for (i = 0; i < betasize; i++) {
				printf ("%8.3e ", gsl_vector_get (s->x, i)); } 
			printf ("f()=%7.3f size=%.3f\n", s->fval, size);
       			if (status == GSL_SUCCESS) {
	   			printf ("Minimum found at:\n");
				printf ("%5d ", iter);
				for (i = 0; i < betasize; i++) {
					printf ("%8.3e ", gsl_vector_get (s->x, i)); } 
				printf ("f()=%7.3f size=%.3f\n", s->fval, size);
			}
		}
      	} while (status == GSL_CONTINUE && iter < MAX_ITERATIONS);
	if (iter == MAX_ITERATIONS)
		printf("Minimization reached maximum number of iterations.");

	gsl_vector_memcpy(est->parameters, s->x);
	gsl_multimin_fminimizer_free(s);
	est->log_likelihood	= -dist.log_likelihood(est->parameters, data);
	//if (est->uses.confidence == 0)
		return est;
}


/** This function goes row by row through <tt>m</tt> and calculates the
likelihood of the given row, putting the result in <tt>v</tt>. Under some
setups, you will need this to find the variance of the estimator.

\param m 	A GSL matrix, exactly like those used for probit, Waring, Gamma, &c MLEs.

\param v	A vector which will hold the likelihood of each row of {{{m}}}. Declare but do not allocate.

\param likelihood_fn	The address of any of the {{{apop_xxx_likelihood}}} functions, e.g., {{{&apop_probit_likelihood}}}.

\param fn_beta	The parameters at which you will evaluate the likelihood. Probably the output of {{{apop_mle_xxx}}}.

<b>Notes, examples</b><br>
Vuong (1989) (<a
href="http://links.jstor.org/sici?sici=0012-9682%28198903%2957%3A2%3C307%3ALRTFMS%3E2.0.CO%3B2-J">Jstor
link</a>) shows that in most cases, the log likelihood ratio is normally
distributed. So to compare a Waring model to a Gamma model:

\code
void compare_gamma_and_waring(gsl_matrix *m, gsl_vector *gamma_parameters, gsl_vector *waring_parameters)
{
gsl_vector      *waring_ll, *gamma_ll;
double          likelihood, mean, std_dev,
                starting_pt_w[2]= {2.12, .40},
                starting_pt_g[2] = {0.12, .40};

        gamma_parameters   = apop_mle_gamma(m, &likelihood, starting_pt_g, .01, 0);
        printf("your most likely gamma parameters are: %g and %g w/log_likelihood: %g\n",
                            gsl_vector_get(gamma_parameters,0), gsl_vector_get(gamma_parameters,1), -likelihood);
        waring_parameters   = apop_mle_waring(m, &likelihood, starting_pt_w, .01, 0);
        printf("your most likely waring parameters are: %g and %g w/log_likelihood: %g\n",
                            gsl_vector_get(waring_parameters,0), gsl_vector_get(waring_parameters,1), -likelihood);
        apop_make_likelihood_vector(m, &gamma_ll, &apop_gamma_likelihood, gamma_parameters);
        apop_make_likelihood_vector(m, &waring_ll, &apop_waring_likelihood, waring_parameters);
        gsl_vector_scale(gamma_ll, -1);        //calculate (waring - gamma) / sqrt(n)
        gsl_vector_add(waring_ll, gamma_ll);
        gsl_vector_scale(waring_ll, 1/sqrt(m->size1));
        mean = apop_mean(waring_ll);
        std_dev  = sqrt(apop_var_m(waring_ll, mean));
        printf("The log likelihood ratio between Waring model and Gamma model is mean %g var %g\n", mean, std_dev);
        if (mean > 0)
           printf("The Waring is a better fit than the Gamma with %g certainty", gsl_cdf_gaussian_P(mean,std_dev));
        else
           printf("The Gamma is a better fit than the Waring with %g certainty", gsl_cdf_gaussian_P(-mean,std_dev));
}
\endcode
\todo This code needs to be rewritten now that everything is objectified. There are better ways to do it anyway.
\ingroup mle */
void apop_make_likelihood_vector(gsl_matrix *m, gsl_vector **v, 
				double (*likelihood_fn)(const gsl_vector *beta, void *d), gsl_vector* fn_beta){
gsl_matrix_view mm;
int             i;
	*v	= gsl_vector_alloc(m->size1);
        for(i=0; i< m->size1; i++){
                mm      = gsl_matrix_submatrix(m,i,0, 1,m->size2);      //get a single row
                gsl_vector_set(*v, i, likelihood_fn(fn_beta, (void *) &(mm.matrix)));
        }
}
