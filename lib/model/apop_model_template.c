/** apop_MODELNAME.c 
 
 (This file is not handled by Doxygen)

Base file copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
Feel free to augment this with your own copyright: modifications (c) you, today.
*/
#include "model.h"

//The default list. You probably don't need them all.
#include "types.h"
#include "conversions.h"
#include "likelihoods.h"
#include "model.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <stdio.h>
#include <assert.h>


/* Every model should have an estimate function. If you are doing an
 MLE, then it should be the one line function here, and then you will need
 to fill in the MODELNAME_log_likelihood function below. If you are not
 doing an MLE, then you won't need the MODELNAME_log_likelihood function,
 but will probably instead do some substantial math here in this function.

 a
*/
static apop_estimate * MODELNAME_estimate(apop_data * data,  void *parameters){
	return apop_maximum_likelihood(data,  apop_MODELNAME, parameters);
}


/*
Often, the input data is a gsl_matrix, and you need to go line by line
through the matrix, calculating the log likelihood. That's the sample
code here.  
*/
static double MODELNAME_log_likelihood(const gsl_vector *beta, void *d){
int		    i;
double	    loglike    	= 0;
gsl_matrix 	*data 		= d;		//type cast void to gsl_matrix.
	for(i=0;i< data->size1; i++){
		loglike      += 1;//PLACE MATH HERE;
	}
	return loglike ;
}

/* The derivative of the MODELNAME distribution, for use in likelihood
  minimization. 
  The format is often the same as above: go line by line through a gsl_matrix.
  The sample is a three-dimensional parameter vector.

You can delete this function entirely if so inclined. If so, remember
to replace this function with NULL in the model definition below.
 */
static void MODELNAME_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
int		    i,j;
double	    dtotal[3];
gsl_matrix 	*data 	= d;		//type cast
    dtotal[0]  = 0,
    dtotal[1]  = 0,
    dtotal[2]  = 0;
	for(i=0; i< data->size1; i++){
		dtotal[0]  += 0; //PLACE MATH HERE
		dtotal[1]  += 0; //PLACE MATH HERE
		dtotal[2]  += 0; //PLACE MATH HERE
	}
	for(j=0; j< beta->size; j++){
	    gsl_vector_set(gradient,j,dtotal[j]);
    }
}


/* For constrained optimizations, you will need a constraint function.
This one will check whether beta[0] > 7. 

You can delete this function entirely if so inclined. If so, remember
to replace this function with NULL in the model definition below.
 */
static double MODELNAME_constraint(gsl_vector *beta, void * d, gsl_vector *returned_beta){
double  limit       = 7,
        tolerance   = 1e-1;
double  mu          = gsl_vector_get(beta, 0);
    if (mu > limit) 
        return 0;
    //else:
    gsl_vector_memcpy(returned_beta, beta);
    gsl_vector_set(returned_beta, 0, limit + tolerance);
    return limit - mu;    
}



/** You may be able to save some time by calculating both log likelihood
and dlog likelihood at the same time. If so, fill in here. 

You can delete this function entirely if so inclined. If so, remember
to replace this function with NULL in the model definition below.
	*/
static void MODELNAME_fdf( const gsl_vector *beta, void *d, double *f, gsl_vector *df){
	*f	= MODELNAME_log_likelihood(beta, d);
	MODELNAME_dlog_likelihood(beta, d, df);
}

static double MODELNAME_rng(gsl_rng* r, gsl_vector * a){
    //place math here.
    return 0;
}

/** The MODELNAME model.
You should describe the format of the input data here.

--The second parameter is either the number of parameters your model has, or a negative number indicating a special case; see the manual.
--The inventory template's elements should be one if you will return the given item; zero if you do not. It is currently correct for an MLE.
--If you deleted any of the functions above, replace their names with NULL here.

\ingroup models
*/
apop_model apop_MODELNAME = {"MODELNAME", -1, 
	MODELNAME_estimate, MODELNAME_log_likelihood, MODELNAME_dlog_likelihood, MODELNAME_fdf, MODELNAME_constraint, MODELNAME_rng};
