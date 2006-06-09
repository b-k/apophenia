/** \file apop_bootstrap.c

Copyright (c) 2006 by Ben Klemens. Licensed under the GNU GPL v2.
 */

/** \defgroup boot Bootstrapping

The jackknife procedure takes the following steps:

for (1000 iterations){      <br>
	generate subset of data     <br>
	run estimation      <br>
	write down parameter estimate       <br>
	}       <br>
calculate variance.       

[If you have a consistent estimator, the mean should be darn close to the mean of the full data set, by the way.]

So that's what the jackknife procedure does. Most of this is trivial: the
real work is in the process of generating the subsets of the data. The
assumption is that your data set is in a gsl_matrix where each row is
a data element, and no rows are special.

\todo It would be nice if one had a means of producing random views of the input data, rather than requiring the copying of half the data set for every run. Todo: write such a function.
*/

#include "stats.h"
#include "output.h"
#include "bootstrap.h"
#include "likelihoods.h"
#include <apophenia/model.h>
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>


/** Give me a data set and a function that goes from data set to
	parameter, and I'll give you the jackknifed standard deviation (sqrt(var)) of the parameter 
 

The function returns a \c gsl_matrix instead of just a double so that you
can jackknife every parameter at once. Remember that if \c e is an
\ref apop_estimate, then \c e->parameters is a \c gsl_vector.

\param in	    The data set. A gsl_matrix where each row is a single data point
\param model    An \ref apop_model, whose \c estimate method will be used here.
\param ep        The \ref apop_estimation_params for your model, to be passed in directly.

\todo I couldn't find a reference on how big the jackknife subsample should be relative to the data set. So I hard-coded the subsample size to 1/3 the original data set. If you have better, code it in.
\ingroup boot
 */
gsl_matrix * apop_jackknife(apop_data *in, apop_model model, apop_estimation_params *ep){
gsl_vector              v;
int                     i;
apop_data               *subset  = apop_data_alloc(in->matrix->size1 - 1, in->matrix->size2);
apop_data               *array_of_boots = NULL;
apop_estimation_params  *e;
apop_estimate           *boot_est;

//Allocate a matrix, get a reduced view of the original, and copy.
gsl_matrix  *reduced= subset->matrix;
gsl_matrix  mv      = gsl_matrix_submatrix(in->matrix, 1,0, in->matrix->size1-1, in->matrix->size2).matrix;
    gsl_matrix_memcpy(reduced, &mv);

    //prep the parameters.
    if (ep){
        e   = malloc(sizeof(*e));
        memcpy(e, ep, sizeof(*e));
    } else
        e   = apop_estimation_params_alloc();
    apop_inventory_set(&(e->uses), 0);      //our re-run will only ask parameters.
    e->uses.parameters  = 1;
	boot_est        = model.estimate(subset, e);
	array_of_boots  = apop_data_alloc(in->matrix->size1, boot_est->parameters->vector->size);
    i = -1;
    //printf("Fuck. %i, %i, %i\n", reduced->size1, i, (i< (int) reduced->size1));
    
    for(i = -1; i< (int) reduced->size1; i++){
        //Get a view of row i, and copy it to position i-1 in the
        //short matrix.
        if (i >= 0){
            v   = gsl_matrix_row(in->matrix, i).vector;
            gsl_matrix_set_row(reduced, i, &v);
	        boot_est        = model.estimate(subset, e);
        }
        gsl_matrix_set_row(array_of_boots->matrix, i+1, boot_est->parameters->vector);
        apop_estimate_free(boot_est);
    }
    apop_data   *out    = apop_data_covar(array_of_boots);
    gsl_matrix_scale(out->matrix, gsl_pow_2(in->matrix->size1-1));
    return out->matrix;
}

/*
                Take everything herein as debris.


gsl_matrix * apop_jackknife_multirow(apop_data *data, apop_model model, apop_estimation_params *ep){
int		        i, j, row;
int             boot_iterations	= 1000;
apop_data	    *subset	        = apop_data_alloc(data->matrix->size1, data->matrix->size2);
apop_estimation_params  *e;
apop_data       *array_of_boots = NULL;
gsl_vector_view	v;
apop_estimate   *boot_est;
gsl_rng		    *rn	            = gsl_rng_alloc(gsl_rng_default);
    if (ep){
        e   = malloc(sizeof(*e));
        memcpy(e, ep, sizeof(*e));
    } else
        e   = apop_estimation_params_alloc();
    apop_inventory_set(&(e->uses), 0);      //our re-run will only ask parameters.
    e->uses.parameters  = 1;
	for (i=0; i<boot_iterations; i++){
		//create the data set
		for (j=0; j< data->matrix->size1; j++){
			row	= (int) gsl_rng_uniform_int(rn, data->matrix->size1);
			v	= gsl_matrix_row(data->matrix, row);
			gsl_matrix_set_row(subset->matrix, j, &(v.vector));
		}
		//get the parameter estimates.
		boot_est    = model.estimate(subset, e);
		if (i==0)
			array_of_boots	= apop_data_alloc(boot_iterations, boot_est->parameters->vector->size);
        gsl_matrix_set_row(array_of_boots->matrix,i,boot_est->parameters->vector);
        apop_estimate_free(boot_est);
	}
    apop_data   *out    = apop_data_covar(array_of_boots); 
    gsl_matrix_scale(out->matrix, gsl_pow_2(in->matrix->size1-1)); 
    return out->matrix;
}

*/


/*
gsl_vector * old_bootstrap(gsl_matrix * data, gsl_vector * (*boot_fn)(gsl_matrix *, void *, void* , void*), int boot_iterations,
	void *params_1, void *params_2, void *params_3) {
int		        i, j, row;
gsl_matrix	    *subset	= gsl_matrix_alloc(data->size1/3, data->size2);
apop_data       *array_of_boots = NULL,
                *summary;
gsl_vector_view	v;
gsl_vector	    *output, *b;
gsl_rng		    *rn	=gsl_rng_alloc(gsl_rng_default);

if (boot_iterations ==0) boot_iterations	= 1000;

	for (i=0; i<boot_iterations; i++){
		//create the data set
		for (j=0; j< data->size1/3; j++){
			row	= (int) gsl_rng_uniform_int(rn, data->size1);
			v	= gsl_matrix_row(data, row);
			gsl_matrix_set_row(subset, row, &(v.vector));
		}
		//get the parameter estimates.
		b	= boot_fn(data, params_1, params_2, params_3);
		if (b!=NULL) 	{
			if (i==0)
				array_of_boots	= apop_data_alloc(boot_iterations, b->size);
			gsl_matrix_set_row(array_of_boots->matrix,i,b);
		} else		i--;

	}
	summary	= apop_data_summarize(array_of_boots);
    apop_data_free(array_of_boots);
	output	= gsl_vector_alloc(summary->matrix->size1);
	gsl_matrix_get_col(output, summary->matrix, 1);
	return output;
}
*/


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
