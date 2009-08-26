/** \file apop_bootstrap.c

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

/** \defgroup boot Bootstrapping

*/

#include "model.h"
#include "variadic.h"
#include "likelihoods.h"

/** Initialize a \c gsl_rng.
 
  Uses the Tausworth routine.

\param  seed    The seed. No need to get funny with it: 0, 1, and 2 will produce wholly different streams.
\return The RNG ready for your use.
\ingroup convenience_fns
*/
gsl_rng *apop_rng_alloc(int seed){
  static int first_use    = 1;
    if (first_use){
       first_use = 0;
       gsl_rng_env_setup();
    }
  gsl_rng *setme  =  gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(setme, seed);
    return setme;
}

/** Give me a data set and a model, and I'll give you the jackknifed covariance matrix of the model parameters.

The basic algorithm for the jackknife (with many details glossed over): create a sequence of data
sets, each with exactly one observation removed, and then produce a new set of parameter estimates 
using that slightly shortened data set. Then, find the covariance matrix of the derived parameters.

Should I use the jackknife or the bootstrap? As a broad rule of thumb, the jackknife works best on models that are closer to linear. The worse a linear approximation does (at the given data), the worse the jackknife approximates the variance.

Sample usage:
\code
apop_data_show(apop_jackknife_cov(your_data, your_model));
\endcode
 
\param in	    The data set. An \ref apop_data set where each row is a single data point.
\param model    An \ref apop_model, that will be used internally by \ref apop_estimate.
            
\return         An \c apop_data set whose matrix element is the estimated covariance matrix of the parameters.

\ingroup boot
 */
apop_data * apop_jackknife_cov(apop_data *in, apop_model model){
  apop_assert(in,  NULL, 0, 's', "You sent me NULL input data.");
  apop_model   *e              = apop_model_copy(model);
  apop_model_clear(in, e);
  int           i, n            = in->matrix->size1;
  apop_data     *subset         = apop_data_alloc(in->vector ? in->vector->size -1 : 0, n - 1, in->matrix->size2);
  apop_data     *array_of_boots = NULL;
  apop_model *overall_est       = apop_estimate(in, *e);
  gsl_vector *overall_params    = apop_data_pack(overall_est->parameters);
    gsl_vector_scale(overall_params, n); //do it just once.
  int           paramct         = overall_params->size;
  gsl_vector    *pseudoval      = gsl_vector_alloc(paramct);

//Allocate a matrix, get a reduced view of the original, and copy.
  gsl_matrix  mv      = gsl_matrix_submatrix(in->matrix, 1,0, n-1, in->matrix->size2).matrix;
    gsl_matrix_memcpy(subset->matrix, &mv);
    if (in->vector){
        gsl_vector v = gsl_vector_subvector(in->vector, 1, n-1).vector;
        gsl_vector_memcpy(subset->vector, &v);
    }
	array_of_boots          = apop_data_alloc(0, n, overall_params->size);
    array_of_boots->names   = apop_name_copy(in->names);

    for(i = -1; i< (int) subset->matrix->size1; i++){
        //Get a view of row i, and copy it to position i-1 in the short matrix.
        if (i >= 0){
            Apop_row(in, i, v);
            gsl_matrix_set_row(subset->matrix, i, v);
            if (subset->vector)
                gsl_vector_set(subset->vector, i, apop_data_get(in, i, -1));
        }
        apop_model *est = apop_estimate(subset, *e);
        gsl_vector *estp = apop_data_pack(est->parameters);
        gsl_vector_memcpy(pseudoval, overall_params);// *n above.
        gsl_vector_scale(estp, n-1);
        gsl_vector_sub(pseudoval, estp);
        gsl_matrix_set_row(array_of_boots->matrix, i+1, pseudoval);
        apop_model_free(est);
        gsl_vector_free(estp);
    }
    apop_data   *out    = apop_data_covariance(array_of_boots);
    gsl_matrix_scale(out->matrix, 1./(n-1.));
    apop_data_free(subset);
    gsl_vector_free(pseudoval);
    apop_model_free(overall_est);
    gsl_vector_free(overall_params);
    apop_model_free(e);
    return out;
}


/** Give me a data set and a model, and I'll give you the bootstrapped covariance matrix of the parameter estimates.

\param data	    The data set. An \c apop_data set where each row is a single data point. (No default)
\param model    An \ref apop_model, whose \c estimate method will be used here. (No default)
\param iterations How many bootstrap draws should I make? (default: 1,000) 
\param rng        An RNG that you have initialized, probably with \c apop_rng_alloc. (Default: see \ref autorng)
\return         An \c apop_data set whose matrix element is the estimated covariance matrix of the parameters.

This function uses the \ref designated syntax for inputs.
\ingroup boot
 */
APOP_VAR_HEAD apop_data * apop_bootstrap_cov(apop_data * data, apop_model model, gsl_rng *rng, int iterations) {
    static gsl_rng *spare = NULL;
    apop_data * apop_varad_var(data, NULL);
    int apop_varad_var(iterations, 1000);
    apop_assert(data, NULL, 0, 's', "The data element can't be NULL.");
    gsl_rng * apop_varad_var(rng, NULL);
    if (!rng && !spare) 
        spare = apop_rng_alloc(++apop_opts.rng_seed);
    if (!rng)  rng = spare;
    return apop_bootstrap_cov_base(data, varad_in.model, rng, iterations);
APOP_VAR_END_HEAD
  apop_model        *e              = apop_model_copy(model);
  apop_model_clear(data, e);
  size_t	        i, j, row;
  apop_data     *subset         = apop_data_alloc(data->vector ? data->vector->size : 0
                                                    , data->matrix->size1, data->matrix->size2);
  apop_data         *array_of_boots = NULL,
                    *summary;
	for (i=0; i<iterations; i++){
		//create the data set
		for (j=0; j< data->matrix->size1; j++){
			row	= gsl_rng_uniform_int(rng, data->matrix->size1);
			APOP_ROW(data, row,v);
			gsl_matrix_set_row(subset->matrix, j, v);
            if (subset->vector)
                gsl_vector_set(subset->vector, j, apop_data_get(data, row, -1));
		}
		//get the parameter estimates.
		apop_model *est = apop_estimate(subset, *e);
        gsl_vector *estp = apop_data_pack(est->parameters);
        if (i==0){
            array_of_boots	        = apop_data_alloc(0,iterations, estp->size);
            array_of_boots->names   = apop_name_copy(data->names);
        }
        gsl_matrix_set_row(array_of_boots->matrix, i, estp);
        apop_model_free(est);
        gsl_vector_free(estp);
	}
	summary	= apop_data_covariance(array_of_boots);
    gsl_matrix_scale(summary->matrix, 1./iterations);
    apop_data_free(array_of_boots);
    apop_data_free(subset);
    apop_model_free(e);
	return summary;
}
