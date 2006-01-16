/** \file apop_stats.c	Basic moments and some distributions.

 Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
 \author Ben Klemens
 */

#include "db.h"     //just for apop_opts
#include "apophenia/stats.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_vector.h>

/** \defgroup basic_stats Some basic statistical functions. 

Many of these are juse one-line convenience functions.

\b moments

\li \ref vector_moments: including mean, variance, kurtosis, and (for pairs of vectors) covariance.
\li \ref matrix_moments: get the covariance matrix from a data matrix.
\li \ref db_moments: var, skew, kurtosis via SQL queries.

\b normalizations

\li \ref apop_normalize_vector: scale and shift a vector.
\li \ref apop_normalize_matrix: remove the mean from each column of the data set.

\b distributions

\li \ref apop_random_beta: Give the mean and variance, and this will draw from the appropriate Beta distribution.
\li \ref apop_multivariate_normal_prob: Evalute a multivariate normal at a given point.
*/

/** \defgroup vector_moments Calculate moments (mean, var, kurtosis) for the data in a gsl_vector.
\ingroup basic_stats

These functions simply take in a GSL vector and return its mean, variance, or kurtosis; the covariance functions take two GSL vectors as inputs.

\ref apop_cov and \ref apop_covar are identical; \ref apop_kurtosis and
\ref apop_kurt are identical; pick the one which sounds better to you.



For \ref apop_var_m<tt>(vector, mean)</tt>, <tt>mean</tt> is the mean of the
vector. This saves the trouble of re-calcuating the mean if you've
already done so. E.g.,

\code
gsl_vector *v;
double mean, var;

//Allocate v and fill it with data here.
mean = apop_mean(v);
var  = apop_var_m(v, mean);
printf("Your vector has mean %g and variance %g\n", mean, var);
\endcode
\ingroup basic_stats */

/** \defgroup matrix_moments   Calculate moments such as the covariance matrix

\ingroup basic_stats */

/** Returns the sum of the data in the given vector.
\ingroup convenience_fns
*/
inline long double apop_sum(gsl_vector *in){
    if (in==NULL){
        if (apop_opts.verbose)
            printf("You just asked me to sum a NULL. Returning zero.\n");
        return 0;
    }
int     i;
double  out = 0;
    for (i=0; i< in->size; i++)
        out += gsl_vector_get(in, i);
	return out; 
}

/** Returns the mean of the data in the given vector.
\ingroup vector_moments
*/
inline double apop_mean(gsl_vector *in){
	return gsl_stats_mean(in->data,in->stride, in->size); }

/** Returns the variance of the data in the given vector.
\ingroup vector_moments
*/
inline double apop_var(gsl_vector *in){
	return gsl_stats_variance(in->data,in->stride, in->size); }

/** Returns the kurtosis of the data in the given vector.
\ingroup vector_moments
*/
inline double apop_kurtosis(gsl_vector *in){
	return gsl_stats_kurtosis(in->data,in->stride, in->size); }

/** Returns the kurtosis of the data in the given vector.
\ingroup vector_moments
*/
inline double apop_kurt(gsl_vector *in){
	return gsl_stats_kurtosis(in->data,in->stride, in->size); }

/** Returns the variance of the data in the given vector, given that you've already calculated the mean.
\param in	the vector in question
\param mean	the mean, which you've already calculated using \ref apop_mean.
\ingroup vector_moments
*/
inline double apop_var_m(gsl_vector *in, double mean){
	return gsl_stats_variance_m(in->data,in->stride, in->size, mean); }

/** returns the covariance of two vectors
\ingroup vector_moments
*/
inline double apop_covar(gsl_vector *ina, gsl_vector *inb){
	return gsl_stats_covariance(ina->data,ina->stride,inb->data,inb->stride,inb->size); }

/** returns the correllation coefficient of two vectors. It's just
\f$ {\hbox{cov}(a,b)\over \sqrt(\hbox{var}(a)) * \sqrt(\hbox{var}(b))}.\f$
\ingroup vector_moments
*/
inline double apop_correlation(gsl_vector *ina, gsl_vector *inb){
	return apop_covar(ina, inb) / sqrt(apop_var(ina) * apop_var(inb));
}

/** returns the covariance of two vectors
\ingroup vector_moments
*/
inline double apop_cov(gsl_vector *ina, gsl_vector *inb){return  apop_covar(ina,inb);}


/** This function will normalize a vector, either such that it has mean
zero and variance one, or such that it ranges between zero and one, or sums to one.

\param in 	A gsl_vector which you have already allocated and filled

\param out 	If normalizing in place, <tt>NULL</tt> (or anything else; it will be ignored).
If not, the address of a <tt>gsl_vector</tt>. Do not allocate.

\param in_place 	0: <tt>in</tt> will not be modified, <tt>out</tt> will be allocated and filled.<br> 
1: <tt>in</tt> will be modified to the appropriate normalization.

\param normalization_type 
0: normalized vector will range between zero and one. Replace each X with (X-min) / (max - min).<br>
1: normalized vector will have mean zero and variance one. Replace
each X with \f$(X-\mu) / \sigma\f$, where \f$\sigma\f$ is the sample
standard deviation.<br>
2: normalized vector will sum to one. E.g., start with a set of observations in bins, end with the percentage of observations in each bin.

\b example 
\code
#include <apophenia/headers.h>

int main(void){
gsl_vector  *in, *out;

in = gsl_vector_calloc(3);
gsl_vector_set(in, 1, 1);
gsl_vector_set(in, 2, 2);

printf("The orignal vector:\n");
apop_vector_print(in, "\t", NULL);

apop_normalize_vector(in, &out, 0, 1);
printf("Normalized with mean zero and variance one:\n");
apop_vector_print(out, "\t", NULL);

apop_normalize_vector(in, &out, 0, 0);
printf("Normalized with max one and min zero:\n");
apop_vector_print(out, "\t", NULL);

apop_normalize_vector(in, NULL, 1, 2);
printf("Normalized into percentages:\n");
apop_vector_print(in, "\t", NULL);

return 0;
}
\endcode
\ingroup basic_stats */
void apop_normalize_vector(gsl_vector *in, gsl_vector **out, int in_place, int normalization_type){
double		mu, min, max;
	if (in_place) 	
		out	= &in;
	else {
		*out 	= gsl_vector_alloc (in->size);
		gsl_vector_memcpy(*out,in);
	}

	if (normalization_type == 1){
		mu	= apop_mean(in);
		gsl_vector_add_constant(*out, -mu);			//subtract the mean
		gsl_vector_scale(*out, 1/(sqrt(apop_var_m(in, mu))));	//divide by the std dev.
	} 
	else if (normalization_type == 0){
		min	= gsl_vector_min(in);
		max	= gsl_vector_max(in);
		gsl_vector_add_constant(*out, -min);
		gsl_vector_scale(*out, 1/(max-min));	

	}
	else if (normalization_type == 2){
		mu	= apop_mean(in);
		gsl_vector_scale(*out, 1/(mu * in->size));	
	}
}

/** For each column in the given matrix, normalize so that the column
has mean zero and variance one. Replace each X with \f$(X-\mu) / \sigma\f$, where \f$\sigma\f$ is the sample standard deviation
\ingroup basic_stats */
void apop_normalize_matrix(gsl_matrix *data){
gsl_vector_view v;
int             j;
        for (j = 0; j < data->size2; j++){
                v	= gsl_matrix_column(data, j);
		gsl_vector_add_constant(&(v.vector), -apop_mean(&(v.vector)));
        }
}

/** Input: any old vector. Output: 1 - the p-value for a chi-squared test to answer the question, "with what confidence can I reject the hypothesis that the variance of my data is zero?"

\param in a gsl_vector of data.
\ingroup asst_tests
 */
inline double apop_test_chi_squared_var_not_zero(gsl_vector *in){
gsl_vector	*normed;
int		i;
double 		sum=0;
	apop_normalize_vector(in,&normed, 0, 1);
	gsl_vector_mul(normed,normed);
	for(i=0;i< normed->size; 
			sum +=gsl_vector_get(normed,i++));
	gsl_vector_free(normed);
	return gsl_cdf_chisq_P(sum,in->size); 
}


inline double apop_double_abs(double a) {if(a>0) return a; else return -a;}

/** This function needs to be replaced with something that handles
names.

\todo Add names to this fn; otherwise apop_print_matrix supercedes it. */
void apop_view_matrix(gsl_matrix *a){
int 		i,j;
	for(i=0;i< a->size1; i++){
		for(j=0;j< a->size2; j++)
			printf("%g ", gsl_matrix_get(a,i,j));
		printf("\n");
	}
}

/** The Beta distribution is useful for modeling because it is bounded
between zero and one, and can be either unimodal (if the variance is low)
or bimodal (if the variance is high), and can have either a slant toward
the bottom or top of the range (depending on the mean).

The distribution has two parameters, typically named \f$\alpha\f$ and \f$\beta\f$, which
can be difficult to interpret. However, there is a one-to-one mapping
between (alpha, beta) pairs and (mean, variance) pairs. Since we have
good intuition about the meaning of means and variances, this function
takes in a mean and variance, calculates alpha and beta behind the scenes,
and returns a random draw from the appropriate Beta distribution.

\param m
The mean the Beta distribution should have. Notice that m
is in [0,1].

\param v
The variance which the Beta distribution should have. It is in (0, 1/12),
where (1/12) is the variance of a Uniform(0,1) distribution. The closer
to 1/12, the worse off you are.

\param r
An already-declared and already-initialized {{{gsl_rng}}}.

\return
Returns one random draw from the given distribution


Example:
\verbatim
	gsl_rng  *r;
	double  a_draw;
	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);
	a_draw = apop_random_beta(.25, 1.0/24.0, r);
\endverbatim
\ingroup convenience_fns
*/
double apop_random_beta(double m, double v, gsl_rng *r) {
double 		k        = (m * (1- m)/ v) -1 ;
        return gsl_ran_beta(r, m* k ,  k*(1 - m) );
}

/**
Evalutates \f[
{\exp(-{1\over 2} (X-\mu)' \Sigma^{-1} (x-\mu))
\over
   sqrt((2 \pi)^n det(\Sigma))}\f]

\param x 
The point at which the multivariate normal should be evaluated

\param mu
The vector of means

\param sigma 
The variance-covariance matrix for the dimensions

\param first_use
I assume you will be evaluating multiple points in the same distribution,
so this saves some calcualtion (notably, the determinant of sigma).<br> 
first_use==1: re-calculate everything<br> 
first_use==0: Assume \c mu and \c sigma have not changed since the last time the function was called.
\ingroup convenience_fns
*/
double apop_multivariate_normal_prob(gsl_vector *x, gsl_vector* mu, gsl_matrix* sigma, int first_use){
static double 	determinant = 0;
static gsl_matrix* inverse = NULL;
static int	dimensions =1;
gsl_vector*	x_minus_mu = gsl_vector_alloc(x->size);
double		numerator;
	gsl_vector_memcpy(x_minus_mu, x);	//copy x to x_minus_mu, then
	gsl_vector_sub(x_minus_mu, mu);		//subtract mu from that copy of x
	if (first_use){
		if (inverse !=NULL) free(inverse);
		dimensions	= x->size;
		inverse 	= gsl_matrix_alloc(dimensions, dimensions);
		determinant	= apop_det_and_inv(sigma, &inverse, 1,1);
	}
	if (determinant == 0) {printf("x"); return(GSL_NEGINF);} //tell minimizer to look elsewhere.
	numerator	= exp(- apop_x_prime_sigma_x(x_minus_mu, inverse) / 2);
printf("(%g %g %g)", numerator, apop_x_prime_sigma_x(x_minus_mu, inverse), (numerator / pow(2 * M_PI, (float)dimensions/2) * sqrt(determinant)));
	return(numerator / pow(2 * M_PI, (float)dimensions/2) * sqrt(determinant));
}

/** give me a random double between min and max [inclusive].

\param min, max 	Y'know.
\param r		The random number generator you're using.
*/
double apop_random_double(double min, double max, gsl_rng *r){
double		base = gsl_rng_uniform(r);
	return base * (max - min) - min;
}

/** Returns a vector of size 101, where returned_vector[95] gives the
value of the 95th percentile, for example. Returned_vector[100] is always
the maximum value, and returned_vector[0] is always the min (regardless
of rounding rule).

\param data	a gsl_vector of data.
\param rounding This will either be 'u' or 'd'. Unless your data is
exactly a multiple of 100, some percentiles will be ambiguous. If 'u',
then round up (use the next highest value); if 'd' (or anything else),
round down to the next lowest value. If 'u', then you can say "5% or
more  of the sample is below returned_vector[5]"; if 'd', then you can
say "5% or more of the sample is above returned_vector[5]".  

\ingroup basic_stats
*/ 
double * apop_percentiles(gsl_vector *data, char rounding){
gsl_vector	*sorted	= gsl_vector_alloc(data->size);
double		*pctiles= malloc(sizeof(double) * 101);
int		i, index;
	gsl_vector_memcpy(sorted,data);
	gsl_sort_vector(sorted);
	for(i=0; i<101; i++){
		index = i*(data->size-1)/100.0;
		if (rounding == 'u' && index != i*(data->size-1)/100.0)
			index ++; //index was rounded down, but should be rounded up.
		pctiles[i]	= gsl_vector_get(sorted, index);
	}
	gsl_vector_free(sorted);
	return pctiles;

}

/** Returns the sum of the elements of a matrix. Occasionally convenient.

  \param m	the matrix to be summed. 
\ingroup convenience_fns*/
long double apop_matrix_sum(gsl_matrix *m){
int 		i,j;
long double	sum	= 0;
	for (j=0; j< m->size1; j++)
		for (i=0; i< m->size2; i++)
			sum     += (gsl_matrix_get(m, j, i));
	return sum;
}

/** Returns the mean of all elements of a matrix.

  \param data	the matrix to be averaged. 
\ingroup convenience_fns*/
double apop_matrix_mean(gsl_matrix *data){
double          avg     = 0;
int             i,j, cnt= 0;
double          x, ratio;
        for(i=0; i < data->size1; i++)
                for(j=0; j < data->size2; j++){
                        x       = gsl_matrix_get(data, i,j);
                        ratio   = cnt/(cnt+1.0);
                        cnt     ++;
                        avg     *= ratio;
                        avg     += x/(cnt +0.0);
                }
	return avg;
}

/** Returns the variance of all elements of a matrix, given the
mean. If you want to calculate both the mean and the variance, use \ref
apop_matrix_mean_and_var.

\param data	the matrix to be averaged. 
\param mean	the pre-calculated mean
\ingroup convenience_fns*/
double apop_matrix_var_m(gsl_matrix *data, double mean){
double          avg2    = 0;
int             i,j, cnt= 0;
double          x, ratio;
        for(i=0; i < data->size1; i++)
                for(j=0; j < data->size2; j++){
                        x       = gsl_matrix_get(data, i,j);
                        ratio   = cnt/(cnt+1.0);
                        cnt     ++;
                        avg2    *= ratio;
                        avg2    += gsl_pow_2(x)/(cnt +0.0);
                }
        return mean - gsl_pow_2(mean); //E[x^2] - E^2[x]
}

/** Returns the mean and variance of all elements of a matrix.

  \param data	the matrix to be averaged. 
\param	mean	where to put the mean to be calculated
\param	var	where to put the variance to be calculated
\ingroup convenience_fns*/
void apop_matrix_mean_and_var(gsl_matrix *data, double *mean, double *var){
double          avg     = 0,
                avg2    = 0;
int             i,j, cnt= 0;
double          x, ratio;
        for(i=0; i < data->size1; i++)
                for(j=0; j < data->size2; j++){
                        x       = gsl_matrix_get(data, i,j);
                        ratio   = cnt/(cnt+1.0);
                        cnt     ++;
                        avg     *= ratio;
                        avg2    *= ratio;
                        avg     += x/(cnt +0.0);
                        avg2    += gsl_pow_2(x)/(cnt +0.0);
                }
	*mean	= avg;
        *var	= avg2 - gsl_pow_2(avg); //E[x^2] - E^2[x]
}

/** Put summary information about the columns of a table (mean, std dev, variance) in a table.

\param indata The table to be summarized. An \ref apop_data structure.
\return     An \ref apop_data structure with one row for each column in the original table, and a column for each summary statistic.
\ingroup    output
\todo At the moment, only gives the mean, standard deviation, and variance
of the data in each column; should give more in the near future.
\todo We should probably let this summarize rows as well.
*/
apop_data * apop_data_summarize(apop_data *indata){
int		    i;
gsl_vector_view	v;
apop_data	*out	= apop_data_alloc(indata->data->size2, 3);
double		mean, stddev,var;
char		rowname[10000]; //crashes on more than 10^9995 columns.
	apop_name_add(out->names, "mean", 'c');
	apop_name_add(out->names, "std dev", 'c');
	apop_name_add(out->names, "variance", 'c');
	if (indata->names !=NULL)
		for (i=0; i< indata->names->colnamect; i++)
			apop_name_add(out->names, indata->names->colnames[i], 'r');
	else
		for (i=0; i< indata->data->size2; i++){
			sprintf(rowname, "col %i", i);
			apop_name_add(out->names, rowname, 'r');
		}
	for (i=0; i< indata->data->size2; i++){
                v       = gsl_matrix_column(indata->data, i);
		mean	= apop_mean(&(v.vector));
		var 	= apop_var_m(&(v.vector),mean);
		stddev	= sqrt(var);
		gsl_matrix_set(out->data, i, 0, mean);
		gsl_matrix_set(out->data, i, 1, stddev);
		gsl_matrix_set(out->data, i, 2, var);
	}	
	return out;
}

/** Put summary information about the columns of a table (mean, std dev, variance) in a table.

 This is just the version of \ref apop_data_summarize for when
 you have a gsl_matrix instead of an \ref apop_data set. In
 fact, here's the source code for this function: <tt>return
 apop_data_summarize(apop_data_from_matrix(m));</tt>

 */
apop_data * apop_matrix_summarize(gsl_matrix *m){
    return apop_data_summarize(apop_data_from_matrix(m));
}
