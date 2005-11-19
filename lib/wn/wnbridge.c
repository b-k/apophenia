#include "db.h"
#include "name.h"
#include "stats.h"
#include "model.h"
#include "output.h"
#include "bootstrap.h"
#include "regression.h"
#include "model.h"
#include "conversions.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <stdio.h>
#include <assert.h>
#include "wnlib.h"

apop_model 	the_dist;
gsl_matrix 	*the_data;
int		the_size;

double function(double v[]){
gsl_vector	*beta	=gsl_vector_alloc(the_size);
int		i;
double		r;
	for (i=0;i< the_size; i++)
		gsl_vector_set(beta, i, v[i]);
	r	=  the_dist.log_likelihood(beta, the_data);
	gsl_vector_free(beta);
	return r;
}
void gradient(double grad[],double v[]){
gsl_vector	*beta	=gsl_vector_alloc(the_size);
gsl_vector 	*g	=gsl_vector_alloc(the_size);
int		i;
	for (i=0;i< the_size; i++)
		gsl_vector_set(beta, i, v[i]);
	the_dist.dlog_likelihood(beta, the_data, g);
	for (i=0;i< the_size; i++)
		grad[i]	= gsl_vector_get(g, i);
}

/** the bridge to Will Naylor's max likelihood thing.

\todo input starting point for vector.
*/
apop_estimate * apop_wn_maximum_likelihood(gsl_matrix * data, apop_inventory *uses,
		                        apop_model dist, double *starting_pt, double step_size, double tolerance, int verbose){

	//some globals:
	the_data	= data;
	the_dist	= dist;
	the_size	= dist.parameter_ct;
double 		*vect;
int 		code, i;
double 		val_min;

  wn_gpmake("no_free");

  wn_make_vect(&vect,the_size);


  vect[0] = 1.3;
  if (the_size > 1) 
	  vect[1] = 0.25;

  wn_conj_direction_method(&code,&val_min, vect,
                           vect,the_size,(function),WN_IHUGE);
  /*
  wn_conj_gradient_method(&code,&val_min,
                          vect,the_size,(function),(gradient),WN_IHUGE);
  */

  wn_print_vect(vect,the_size);
  if (apop_verbose){
  	printf("final result: code = %d   ",code);
  	printf("    ob = %lf\n",val_min);
}

  wn_gpfree();

apop_inventory	inv;
	apop_inventory_set(&inv, 1);
apop_estimate	*est = apop_estimate_alloc(data->size1, data->size2, NULL, inv);
  	for (i=0;i<the_size; i++)
		gsl_vector_set(est->parameters, i, vect[i]);
	est->log_likelihood	= -dist.log_likelihood(est->parameters, data);
	est->status		= code;
	//Epilogue:
	//find the variance-covariance matrix, using $df/d\theta \cdot df/d\theta$
	if (est->uses.covariance == 0 || dist.dlog_likelihood == NULL) 
		return est;
	//else:
gsl_matrix	*pre_cov;
gsl_vector	*diff;
gsl_vector_view	v;
	pre_cov			= gsl_matrix_alloc(the_size, the_size);
	//estimate->covariance	= gsl_matrix_alloc(betasize, betasize);
	diff			= gsl_vector_alloc(the_size);
	dist.dlog_likelihood(est->parameters, data, diff);
	for (i=0; i< the_size; i++){
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
	for (i=0; i<the_size; i++) // confidence[i] = |1 - (1-N(Mu[i],sigma[i]))*2|
		gsl_vector_set(est->confidence, i,
			fabs(1 - (1 - gsl_cdf_gaussian_P(gsl_vector_get(est->parameters, i), 
			gsl_matrix_get(est->covariance, i, i)))*2));
	return est;
}

