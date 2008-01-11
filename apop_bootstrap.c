/** \file apop_bootstrap.c

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

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
  gsl_rng *setme  =  gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(setme, seed);
    return setme;
}

/** Give me a data set and a model, and I'll give you the jackknifed covariance matrix of the model parameters.
 
\param in	    The data set. An \c apop_data set where each row is a single data point.
\param model    An \ref apop_model, whose \c estimate method will be used here.
            Use a stripped- down version that sets \c want_cov, \c want_expected_value, and anything else to zero, because the jackknife only uses the parameter estimates. If your model calls this function to estimate the covariance and \c want_cov=1, then you've got an infinite loop.
\return         An \c apop_data set whose matrix element is the estimated covariance matrix of the parameters.

\ingroup boot
 */
apop_data * apop_jackknife_cov(apop_data *in, apop_model model){
  apop_assert(in,  NULL, 0, 's', "You sent me NULL input data.");
  apop_assert(in->matrix,  NULL, 0, 's', "At the moment, this function uses only the matrix element of the input data.");
  apop_model   *e              = apop_model_copy(model);
  apop_model_clear(in, e);
  int           i, n            = in->matrix->size1;
  apop_data     *subset         = apop_data_alloc(0, n - 1, in->matrix->size2);
  apop_data     *array_of_boots = NULL;
  apop_model *overall_est       = model.estimate(in, e);
    gsl_vector_scale(overall_est->parameters->vector, n); //do it just once.
  int           paramct         = overall_est->parameters->vector->size;
  gsl_vector    *pseudoval      = gsl_vector_alloc(paramct);

//Allocate a matrix, get a reduced view of the original, and copy.
  gsl_matrix  mv      = gsl_matrix_submatrix(in->matrix, 1,0, n-1, in->matrix->size2).matrix;
    gsl_matrix_memcpy(subset->matrix, &mv);

	array_of_boots          = apop_data_alloc(0, n, overall_est->parameters->vector->size);
    array_of_boots->names   = apop_name_copy(in->names);
    for(i = -1; i< (int) subset->matrix->size1; i++){
        //Get a view of row i, and copy it to position i-1 in the short matrix.
        if (i >= 0){
            APOP_ROW(in, i, v);
            gsl_matrix_set_row(subset->matrix, i, v);
        }
        apop_model *est = model.estimate(subset, e);
        gsl_vector_memcpy(pseudoval, overall_est->parameters->vector);// *n above.
        gsl_vector_scale(est->parameters->vector, n-1);
        gsl_vector_sub(pseudoval, est->parameters->vector);
        gsl_matrix_set_row(array_of_boots->matrix, i+1, pseudoval);
        apop_model_free(est);
    }
    apop_data   *out    = apop_data_covar(array_of_boots);
    gsl_matrix_scale(out->matrix, 1./(n-1.));
    apop_data_free(subset);
    gsl_vector_free(pseudoval);
    apop_model_free(overall_est);
    apop_model_free(e);
    return out;
}


/** Give me a data set and a model, and I'll give you the bootstrapped covariance matrix of the parameter estimates.

\param data	    The data set. An \c apop_data set where each row is a single data point.
\param model    An \ref apop_model, whose \c estimate method will be used here.
\param r        An RNG that you have initialized (probably with \c apop_rng_alloc)
\param boot_iterations How many bootstrap draws should I make? A positive integer; if you express indifference by specifying zero, I'll make 1,000 draws.
\return         An \c apop_data set whose matrix element is the estimated covariance matrix of the parameters.

\ingroup boot
 */
apop_data * apop_bootstrap_cov(apop_data * data, apop_model model, gsl_rng *r, int boot_iterations) {
    if (boot_iterations ==0)    boot_iterations	= 1000;
  apop_model        *e              = apop_model_copy(model);
  apop_model_clear(data, e);
  size_t	        i, j, row;
  apop_data	        *subset	        = apop_data_alloc(0,data->matrix->size1, data->matrix->size2);
  apop_data         *array_of_boots = NULL,
                    *summary;
	for (i=0; i<boot_iterations; i++){
		//create the data set
		for (j=0; j< data->matrix->size1; j++){
			row	= gsl_rng_uniform_int(r, data->matrix->size1);
			APOP_ROW(data, row,v);
			gsl_matrix_set_row(subset->matrix, j, v);
		}
		//get the parameter estimates.
		apop_model *est = model.estimate(subset, e);
        if (i==0){
            array_of_boots	        = apop_data_alloc(0,boot_iterations, est->parameters->vector->size);
            array_of_boots->names   = apop_name_copy(data->names);
        }
        gsl_matrix_set_row(array_of_boots->matrix, i, est->parameters->vector);
        apop_model_free(est);
	}
	summary	= apop_data_covariance(array_of_boots);
    gsl_matrix_scale(summary->matrix, 1./boot_iterations);
    apop_data_free(array_of_boots);
    apop_data_free(subset);
    apop_model_free(e);
	return summary;
}
