/** \file apop_bootstrap.c

Bootstrapping!!!

The bootstrap procedure takes the following steps:

for (1000 iterations){
	generate subset of data
	run estimation
	write down parameter estimate
	}
calculate variance. 

[If you have a consistent estimator, the mean should be darn close to the mean of the full data set, by the way.]

So that's what the bootstrap procedure does. Most of this is trivial: the
real work is in the process of generating the subsets of the data. The
assumption is that your data set is in a gsl_matrix where each row is
a data element, and no rows are special.

You need to provide a function which takes a data set as an input and spits out a single number as output.

\todo It would be nice if one had a means of producing random views of the input data, rather than requiring the copying of half the data set for every run. Todo: write such a function.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL.
*/

#include "stats.h"
#include "output.h"
#include "bootstrap.h"
#include "likelihoods.h"
#include <apophenia/model.h>
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>


/** Give me a data set and a function that goes from data set to
	parameter, and I'll give you the bootstrapped standard deviation (sqrt(var)) of the parameter 
 

The function you write will have the following header:
\code{	gsl_vector * boot_fn(gsl_matrix * data, void * param_1, void*  param_2, void* param_3);}
You have _three_ parameters that you can input. Of course, you technically
only need one, since that one can be a struct with many elements, but
having to write a struct to pass in two variables is annoying. This
form pushes the problem back to when you have four or more parameters,
at which point they're probably already in a struct.

The function returns a gsl_vector instead of just a double so that you
can bootstrap every parameter at once. Remember that if \code{e} is an
\ref apop_estimate, then \code{e->parameters} is a gsl_vector.

\param	data	The data set. A gsl_matrix where each row is a single data point
\param boot_fn	The function by which parameter estimates are found; see the notes.
\param boot_iterations	How many subsamples to draw. If you forget and set this to zero, then I assume 1000
\param params_1	a parameter to send to your boot_fn.
\param params_2	a parameter to send to your boot_fn.
\param params_3	a parameter to send to your boot_fn.

\todo I couldn't find a reference on how big the bootstrap subsample should be relative to the data set. So I hard-coded the subsample size to 1/3 the original data set. If you can find better, code it in.
 */
gsl_vector * bootstrap(gsl_matrix * data, gsl_vector * (*boot_fn)(gsl_matrix *, void *, void* , void*), int boot_iterations,
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
			gsl_matrix_set_row(array_of_boots->data,i,b);
		} else		i--;

	}
	summary	= apop_matrix_summarize(array_of_boots);
    apop_data_free(array_of_boots);
	output	= gsl_vector_alloc(summary->data->size1);
	gsl_matrix_get_col(output, summary->data, 1);
	return output;
}

