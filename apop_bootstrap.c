/** \file apop_bootstrap.c

Copyright (c) 2006 by Ben Klemens. Licensed under the GNU GPL v2.
 */

/** \defgroup boot Bootstrapping

There is currently only one jackknifing procedure, to find the covariance
matrix of a set of parameters.  See the Monte Carlo chapter of the book
for details.

\todo It would be nice if one had a means of producing random views of the input data, rather than requiring the copying of half the data set for every run. Todo: write such a function.
*/

#include "stats.h"
#include "output.h"
#include "bootstrap.h"
#include "likelihoods.h"
#include <apophenia/model.h>
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>


static apop_ep *parameter_prep(apop_ep *in){
  apop_ep *e;
    if (in){
        e   = malloc(sizeof(*e));
        memcpy(e, in, sizeof(*e));
    } else
        e   = apop_ep_alloc();
    apop_inventory_set(&(e->uses), 0);      //our re-run will only ask parameters.
    e->uses.parameters  = 1;
    return e;
}

/** Give me a data set and a function that goes from data set to
	parameter, and I'll give you the jackknifed standard deviation (sqrt(var)) of the parameter 
 

The function returns a full matrix so you can jackknife every parameter
at once. 

\param in	    The data set. A gsl_matrix where each row is a single data point
\param model    An \ref apop_model, whose \c estimate method will be used here.
\param ep        The \ref apop_ep for your model, to be passed in directly.
\return         An \c apop_data set whose matrix element is the estimated covariance matrix of the parameters.

\ingroup boot
 */
apop_data * apop_jackknife_cov(apop_data *in, apop_model model, apop_ep *ep){
  int           i;
  apop_estimate *boot_est;
  apop_ep       *e              = parameter_prep(ep);
  int           n               = in->matrix->size1;
  apop_data     *subset         = apop_data_alloc(in->matrix->size1 - 1, in->matrix->size2);
  apop_data     *array_of_boots = NULL;
  apop_estimate *overall_est    = model.estimate(subset, e);
  int           paramct         = overall_est->parameters->vector->size;
  gsl_vector    *pseudoval      = gsl_vector_alloc(paramct);

//Allocate a matrix, get a reduced view of the original, and copy.
  gsl_matrix  mv      = gsl_matrix_submatrix(in->matrix, 1,0, in->matrix->size1-1, in->matrix->size2).matrix;
    gsl_matrix_memcpy(subset->matrix, &mv);

	array_of_boots          = apop_data_alloc(in->matrix->size1, overall_est->parameters->vector->size);
    array_of_boots->names   = in->names;
    for(i = -1; i< (int) subset->matrix->size1; i++){
        //Get a view of row i, and copy it to position i-1 in the
        //short matrix.
        if (i >= 0){
            APOP_ROW(in, i, v);
            gsl_matrix_set_row(subset->matrix, i, v);
        }
        boot_est        = model.estimate(subset, e);
        gsl_vector_memcpy(pseudoval, overall_est->parameters->vector);
        gsl_vector_scale(pseudoval, n);
        gsl_vector_scale(boot_est->parameters->vector, n-1);
        gsl_vector_sub(pseudoval, boot_est->parameters->vector);
        gsl_matrix_set_row(array_of_boots->matrix, i+1, pseudoval);
        apop_estimate_free(boot_est);
    }
    apop_data   *out    = apop_data_covar(array_of_boots);
    apop_data_free(subset);
    gsl_vector_free(pseudoval);
    apop_estimate_free(overall_est);
    return out;
}

/** Initialize an RNG.
 
  Uses the Tausworth routine.

\param  seed    The seed. No need to get funny with it: 0, 1, and 2 will produce wholly different streams.
\return The RNG ready for your use.
\ingroup convenience_fns
*/
gsl_rng *apop_rng_alloc(int seed){
  static int first_use    = 1;
    if (first_use){
       first_use --;
       gsl_rng_env_setup();
    }
  gsl_rng *setme  =  gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(setme, seed);
    return setme;
}



/*
                Take everything herein as debris.

apop_data * apop_bootstrap_cov(apop_data * data, apop_model model, apop_ep *epin, gsl_rng *r, float boot_size_d, int boot_iterations) {
    if (boot_iterations ==0)    boot_iterations	= 1000;
    if (boot_size_d ==0)        boot_size_d	= 3;
  apop_ep           *ep  = parameter_prep(epin);
  size_t	        i, j, row;
  apop_data	        *subset	= apop_data_alloc(data->matrix->size1*boot_size_d, data->matrix->size2);
  apop_data         *array_of_boots = NULL,
                    *summary;
  apop_estimate     *e;
	for (i=0; i<boot_iterations; i++){
		//create the data set
		for (j=0; j< data->matrix->size1*boot_size_d; j++){
			row	= gsl_rng_uniform_int(r, data->matrix->size1);
			APOP_ROW(data, row,v);
			gsl_matrix_set_row(subset->matrix, j, v);
		}
		//get the parameter estimates.
		e   = model.estimate(subset, ep);
		if (!e) i--;
        else {
			if (i==0){
				array_of_boots	        = apop_data_alloc(boot_iterations, e->parameters->vector->size);
                array_of_boots->names   = data->names;
            }
			gsl_matrix_set_row(array_of_boots->matrix,i,e->parameters->vector);
            apop_estimate_free(e);
		} 
	}
	summary	= apop_data_covariance_matrix(array_of_boots, 1);
    apop_data_free(array_of_boots);
    apop_data_free(subset);
	return summary;
}
*/
