
/** \file apop_bootstrap.c

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the GPLv2; see COPYING.  */

#include "apop_internal.h"

/** Initialize a \c gsl_rng.
 
  Uses the Tausworth routine.

\param  seed    The seed. No need to get funny with it: 0, 1, and 2 will produce wholly different streams.
\return The RNG ready for your use.

\li If you are confident that your code is debugged and would like a new stream of values every time your program runs (provided your runs are more than a second apart), seed with the time:

\include draw_some_normals.c
*/
gsl_rng *apop_rng_alloc(int seed){
    static int first_use = 1;
    if (first_use){
       first_use = 0;
       OMP_critical(rng_env_setup) //GSL makes vague promises about thread-safety
       gsl_rng_env_setup();
    }
    gsl_rng *setme = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(setme, seed);
    return setme;
}

/** Give me a data set and a model, and I'll give you the jackknifed covariance matrix of the model parameters.

The basic algorithm for the jackknife (glossing over the details): create a sequence of data
sets, each with exactly one observation removed, and then produce a new set of parameter estimates 
using that slightly shortened data set. Then, find the covariance matrix of the derived parameters.

\li Jackknife or bootstrap? As a broad rule of thumb, the jackknife works best on models
    that are closer to linear. The worse a linear approximation does (at the given data),
    the worse the jackknife approximates the variance.

\param in	    The data set. An \ref apop_data set where each row is a single data point.
\param model    An \ref apop_model, that will be used internally by \ref apop_estimate.
            
\exception out->error=='n'   \c NULL input data.
\return         An \c apop_data set whose matrix element is the estimated covariance matrix of the parameters.
\see apop_bootstrap_cov

For example:
\include jack.c
*/
apop_data * apop_jackknife_cov(apop_data *in, apop_model *model){
    Apop_stopif(!in, apop_return_data_error(n), 0, "The data input can't be NULL.");
    Get_vmsizes(in); //msize1, msize2, vsize
    apop_model *e = apop_model_copy(model);
    int i, n = GSL_MAX(msize1, GSL_MAX(vsize, in->textsize[0]));
    apop_model *overall_est = e->parameters ? e : apop_estimate(in, e);//if not estimated, do so
    gsl_vector *overall_params = apop_data_pack(overall_est->parameters);
    gsl_vector_scale(overall_params, n); //do it just once.
    gsl_vector *pseudoval = gsl_vector_alloc(overall_params->size);

    //Copy the original, minus the first row.
    apop_data *subset = apop_data_copy(Apop_rs(in, 1, n-1));
    apop_name *tmpnames = in->names; 
    in->names = NULL;  //save on some copying below.

    apop_data *array_of_boots = apop_data_alloc(n, overall_params->size);

    for(i = -1; i< n-1; i++){
        //Get a view of row i, and copy it to position i-1 in the short matrix.
        if (i >= 0) apop_data_memcpy(Apop_r(subset, i), Apop_r(in, i));
        apop_model *est = apop_estimate(subset, e);
        gsl_vector *estp = apop_data_pack(est->parameters);
        gsl_vector_memcpy(pseudoval, overall_params);// *n above.
        gsl_vector_scale(estp, n-1);
        gsl_vector_sub(pseudoval, estp);
        gsl_matrix_set_row(array_of_boots->matrix, i+1, pseudoval);
        apop_model_free(est);
        gsl_vector_free(estp);
    }
    in->names = tmpnames;
    apop_data *out = apop_data_covariance(array_of_boots);
    gsl_matrix_scale(out->matrix, 1./(n-1.));
    apop_data_free(subset);
    gsl_vector_free(pseudoval);
    apop_data_free(array_of_boots);
    if (e!=overall_est)
        apop_model_free(overall_est);
    apop_model_free(e);
    gsl_vector_free(overall_params);
    return out;
}

/** Give me a data set and a model, and I'll give you the bootstrapped covariance matrix of the parameter estimates.

\param data	    The data set. An \c apop_data set where each row is a single data point. (No default)
\param model    An \ref apop_model, whose \c estimate method will be used here. (No default)
\param iterations How many bootstrap draws should I make? (default: 1,000) 
\param rng        An RNG that you have initialized, probably with \c apop_rng_alloc. (Default: an RNG from \ref apop_rng_get_thread)
\param boot_store  If not \c NULL, put the list of drawn parameter values here, with one parameter set per row. Sample use: 
\code
apop_data *boots;
apop_bootstrap_cov(data, model, .boot_store=&boots);
apop_data_print(boots);
\endcode
The rows are packed via \ref apop_data_pack, so use \ref apop_data_unpack if needed. (Default: \c NULL)
\param ignore_nans If \c 'y' and any of the elements in the estimation return \c NaN, then I will throw out that draw and try again. If \c 'n', then I will write that set of statistics to the list, \c NaN and all. I keep count of throw-aways; if there are more than \c iterations elements thrown out, then I throw an error and return with estimates using data I have so far. That is, I assume that \c NaNs are rare edge cases; if they are as common as good data, you might want to rethink how you are using the bootstrap mechanism. (Default: 'n')
\return         An \c apop_data set whose matrix element is the estimated covariance matrix of the parameters.
\exception out->error=='n'   \c NULL input data.
\exception out->error=='N'   \c too many NaNs.

\li This function uses the \ref designated syntax for inputs.

This example is a sort of demonstration of the Central Limit Theorem. The model is
a simulation, where each call to the estimation routine produces the mean/std dev of
a set of draws from a Uniform Distribution. Because the simulation takes no inputs,
\ref apop_bootstrap_cov simply re-runs the simulation and calculates a sequence of
mean/std dev pairs, and reports the covariance of that generated data set.

\include boot_clt.c

\see apop_jackknife_cov
 */
#ifdef APOP_NO_VARIADIC
apop_data * apop_bootstrap_cov(apop_data * data, apop_model *model, gsl_rng *rng, int iterations, char keep_boots, char ignore_nans, apop_data **boot_store){
#else
apop_varad_head(apop_data *, apop_bootstrap_cov) {
    apop_data * apop_varad_var(data, NULL);
    apop_model *model = varad_in.model;
    int apop_varad_var(iterations, 1000);
    gsl_rng * apop_varad_var(rng, apop_rng_get_thread());
    char apop_varad_var(keep_boots, 'n');
    apop_data** apop_varad_var(boot_store, NULL);
    char apop_varad_var(ignore_nans, 'n');
    return apop_bootstrap_cov_base(data, model, rng, iterations, keep_boots, ignore_nans, boot_store);
}

 apop_data * apop_bootstrap_cov_base(apop_data * data, apop_model *model, gsl_rng *rng, int iterations, char keep_boots, char ignore_nans, apop_data **boot_store){
#endif
    Get_vmsizes(data); //vsize, msize1, msize2
    apop_model *e = apop_model_copy(model);
    apop_data *subset = apop_data_copy(data);
    apop_data *array_of_boots = NULL,
              *summary;
    //prevent and infinite regression of covariance calculation.
    Apop_model_add_group(e, apop_parts_wanted); //default wants for nothing.
    size_t i, nan_draws=0;
    apop_name *tmpnames = (data && data->names) ? data->names : NULL; //save on some copying below.
    if (data && data->names) data->names = NULL;

    int height = GSL_MAX(msize1, GSL_MAX(vsize, (data?(*data->textsize):0)));
	for (i=0; i<iterations && nan_draws < iterations; i++){
		for (size_t j=0; j< height; j++){       //create the data set
			size_t randrow	= gsl_rng_uniform_int(rng, height);
            apop_data_memcpy(Apop_r(subset, j), Apop_r(data, randrow));
		}
		//get the parameter estimates.
		apop_model *est = apop_estimate(subset, e);
        gsl_vector *estp = apop_data_pack(est->parameters);
        if (!gsl_isnan(apop_sum(estp))){
            if (i==0){
                array_of_boots	      = apop_data_alloc(iterations, estp->size);
                apop_name_stack(array_of_boots->names, est->parameters->names, 'c', 'v');
                apop_name_stack(array_of_boots->names, est->parameters->names, 'c', 'c');
                apop_name_stack(array_of_boots->names, est->parameters->names, 'c', 'r');
            }
            gsl_matrix_set_row(array_of_boots->matrix, i, estp);
        } else if (ignore_nans=='y'){
            i--; 
            nan_draws++;
        }
        apop_model_free(est);
        gsl_vector_free(estp);
	}
    if(data) data->names = tmpnames;
    apop_data_free(subset);
    apop_model_free(e);
    int set_error=0;
    Apop_stopif(i == 0 && nan_draws == iterations, apop_return_data_error(N),
                1, "I ran into %i NaNs and no not-NaN estimations, and so stopped. "
                       , iterations);
    Apop_stopif(nan_draws == iterations,  set_error++;
            apop_matrix_realloc(array_of_boots->matrix, i, array_of_boots->matrix->size2),
                1, "I ran into %i NaNs, and so stopped. Returning results based "
                       "on %zu bootstrap iterations.", iterations, i);
	summary	= apop_data_covariance(array_of_boots);
    if (boot_store) *boot_store = array_of_boots;
    else            apop_data_free(array_of_boots);
    if (set_error) summary->error = 'N';
	return summary;
}
