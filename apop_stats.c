/** \file apop_stats.c	Basic moments and some distributions.

 \author Ben Klemens
Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "db.h"     //just for apop_opts
#include "stats.h"
#include <gsl/gsl_rng.h>

/** \defgroup basic_stats Some basic statistical functions

Many of these are just one-line convenience functions for finding moments and normalizing matrices.
*/

/** \defgroup vector_moments Calculate moments (mean, var, kurtosis) for the data in a gsl_vector.
\ingroup basic_stats

These functions simply take in a GSL vector and return its mean, variance, or kurtosis; the covariance functions take two GSL vectors as inputs.

\ref apop_vector_kurtosis and \ref apop_vector_kurt are identical; pick the one which sounds better to you.

See also \ref db_moments.


For \ref apop_vector_var_m<tt>(vector, mean)</tt>, <tt>mean</tt> is the mean of the
vector. This saves the trouble of re-calcuating the mean if you've
already done so. E.g.,

\code
gsl_vector *v;
double mean, var;

//Allocate v and fill it with data here.
mean = apop_vector_mean(v);
var  = apop_vector_var_m(v, mean);
printf("Your vector has mean %g and variance %g\n", mean, var);
\endcode
\ingroup basic_stats */

/** \defgroup matrix_moments   Calculate moments such as the covariance matrix

\ingroup basic_stats */

/** Returns the sum of the data in the given vector.
\ingroup convenience_fns
*/
inline long double apop_vector_sum(const gsl_vector *in){
  apop_assert(in, 0, 1,'c', "You just asked me to sum a NULL. Returning zero.\n")
  int     i;
  double  out = 0;
    for (i=0; i< in->size; i++)
        out += gsl_vector_get(in, i);
	return out; 
}

/** Returns the sum of the data in the given vector.

  An alias for \ref apop_vector_sum.
\ingroup convenience_fns
*/
inline long double apop_sum(const gsl_vector *in){
    return apop_vector_sum(in);
}

/** Returns the mean of the data in the given vector.
\ingroup vector_moments
*/
inline double apop_vector_mean(const gsl_vector *in){
	return gsl_stats_mean(in->data,in->stride, in->size); }

/** Returns the mean of the data in the given vector.

  An alias for \ref apop_vector_mean.
\ingroup vector_moments
*/
inline double apop_mean(const gsl_vector *in){
	return apop_vector_mean(in); 
}

/** Returns the variance of the data in the given vector.
\ingroup vector_moments
*/
inline double apop_vector_var(const gsl_vector *in){
	return gsl_stats_variance(in->data,in->stride, in->size); }

/** Returns the variance of the data in the given vector.

  An alias for \ref apop_vector_var.
\ingroup vector_moments
*/
inline double apop_var(const gsl_vector *in){
	return apop_vector_var(in); 
}

/** Returns the sample skew (divide by \f$n-1\f$) of the data in the given vector.
\ingroup vector_moments
*/
inline double apop_vector_skew(const gsl_vector *in){
	return gsl_stats_skew(in->data,in->stride, in->size)
                *pow(apop_vector_var(in),3./2)*in->size/(in->size -1.); }

/** Returns the sample kurtosis (divide by \f$n-1\f$) of the data in the given vector.
  This does not normalize the output: the kurtosis of a \f${\cal N}(0,1)\f$ is three, not zero.
\ingroup vector_moments
*/
inline double apop_vector_kurtosis(const gsl_vector *in){
	return ((gsl_stats_kurtosis(in->data,in->stride, in->size)+3)
                *gsl_pow_2(apop_vector_var(in)))*in->size/(in->size -1.); }

/** Returns the sample kurtosis (divide by \f$n-1\f$) of the data in the given vector.
  This does not normalize the output: the kurtosis of a \f${\cal N}(0,1)\f$ is three, not zero.
\ingroup vector_moments
*/
inline double apop_vector_kurt(const gsl_vector *in){
	return apop_vector_kurtosis(in);}

/** Returns the variance of the data in the given vector, given that you've already calculated the mean.
\param in	the vector in question
\param mean	the mean, which you've already calculated using \ref apop_vector_mean.
\ingroup vector_moments
*/
inline double apop_vector_var_m(const gsl_vector *in, const double mean){
	return gsl_stats_variance_m(in->data,in->stride, in->size, mean); }

/** returns the covariance of two vectors
\ingroup vector_moments
*/
inline double apop_vector_cov(const gsl_vector *ina, const gsl_vector *inb){
	return gsl_stats_covariance(ina->data,ina->stride,inb->data,inb->stride,inb->size); }

/** returns the correllation coefficient of two vectors. It's just
\f$ {\hbox{cov}(a,b)\over \sqrt(\hbox{var}(a)) * \sqrt(\hbox{var}(b))}.\f$
\ingroup vector_moments
*/
inline double apop_vector_correlation(const gsl_vector *ina, const gsl_vector *inb){
	return apop_vector_cov(ina, inb) / sqrt(apop_vector_var(ina) * apop_vector_var(inb));
}


/** returns the scalar distance (standard Euclidian metric) between two vectors. Simply \f$\sqrt{\sum_i{(a_i - b_i)^2}},\f$
where \f$i\f$ iterates over dimensions.

\ingroup convenience_fns
*/
double apop_vector_distance(const gsl_vector *ina, const gsl_vector *inb){
  apop_assert(ina, 0, 0, 'c', "first input vector is NULL. Returning 0.\n");
  apop_assert(inb, 0, 0, 'c', "second input vector is NULL. Returning 0.\n");
  apop_assert(ina->size == inb->size, 0, 0,'c', 
                "You sent a vector of size %u and a vector of size %u. Returning zero.\n", ina->size, inb->size);
  double  dist    = 0;
  size_t  i;
    for (i=0; i< ina->size; i++){
        dist    += gsl_pow_2(gsl_vector_get(ina, i) - gsl_vector_get(inb, i));
    }
	return sqrt(dist); 
}

/** returns the scalar Manhattan metric distance  between two vectors. Simply \f$\sum_i{|a_i - b_i|},\f$
where \f$i\f$ iterates over dimensions.

\ingroup convenience_fns
*/
double apop_vector_grid_distance(const gsl_vector *ina, const gsl_vector *inb){
  apop_assert(ina, 0, 0, 'c', "first input vector is NULL. Returning 0.\n");
  apop_assert(inb, 0, 0, 'c', "second input vector is NULL. Returning 0.\n");
  apop_assert(ina->size == inb->size, 0, 0,'c', 
                "You sent a vector of size %u and a vector of size %u. Returning zero.\n", ina->size, inb->size);
  double  dist    = 0;
  size_t  i;
    for (i=0; i< ina->size; i++){
        dist    += fabs(gsl_vector_get(ina, i) - gsl_vector_get(inb, i));
    }
	return dist; 
}

/** This function will normalize a vector, either such that it has mean
zero and variance one, or such that it ranges between zero and one, or sums to one.

\param in 	A gsl_vector which you have already allocated and filled

\param out 	If normalizing in place, <tt>NULL</tt>.
If not, the address of a <tt>gsl_vector</tt>. Do not allocate.

\param normalization_type 
'r': normalized vector will range between zero and one. Replace each X with (X-min) / (max - min).<br>
's': normalized vector will have mean zero and variance one. Replace
each X with \f$(X-\mu) / \sigma\f$, where \f$\sigma\f$ is the sample
standard deviation.<br>
'p': normalized vector will sum to one. E.g., start with a set of observations in bins, end with the percentage of observations in each bin.<br>
'm': normalize to mean zero: Replace each X with \f$(X-\mu)\f$<br>

\b Example 
\code
#include <apop.h>

int main(void){
gsl_vector  *in, *out;

in = gsl_vector_calloc(3);
gsl_vector_set(in, 1, 1);
gsl_vector_set(in, 2, 2);

printf("The orignal vector:\n");
apop_vector_show(in);

apop_vector_normalize(in, &out, 's');
printf("Standardized with mean zero and variance one:\n");
apop_vector_show(out);

apop_vector_normalize(in, &out, 'r');
printf("Normalized range with max one and min zero:\n");
apop_vector_show(out);

apop_vector_normalize(in, NULL, 'p');
printf("Normalized into percentages:\n");
apop_vector_show(in);
}
\endcode
\ingroup basic_stats */
void apop_vector_normalize(gsl_vector *in, gsl_vector **out, const char normalization_type){
  if (!in){
      apop_error(1, 'c', "Input vector is NULL. Doing nothing.\n");
      return; }
  double		mu, min, max;
	if (!out) 	
		out	= &in;
	else {
		*out 	= gsl_vector_alloc (in->size);
		gsl_vector_memcpy(*out,in);
	}
        //the numbers are deprecated and will go away.
	if ((normalization_type == 's')|| (normalization_type == 1)){
		mu	= apop_vector_mean(in);
		gsl_vector_add_constant(*out, -mu);
		gsl_vector_scale(*out, 1/(sqrt(apop_vector_var_m(in, mu))));
	} 
	else if ((normalization_type == 'r')||(normalization_type == 0)){
        gsl_vector_minmax(in, &min, &max);
		gsl_vector_add_constant(*out, -min);
		gsl_vector_scale(*out, 1/(max-min));	

	}
	else if ((normalization_type == 'p') ||(normalization_type == 2)){
		mu	= apop_vector_mean(in);
		gsl_vector_scale(*out, 1/(mu * in->size));	
	}
	else if ((normalization_type == 'm')){
		mu	= apop_vector_mean(in);
		gsl_vector_add_constant(*out, -mu);
	}
}

/** Normalize  each row or column in the given matrix, one by one.

  Basically just a convenience fn to iterate through the columns and run \ref apop_vector_normalize for you.

\param data     The data set to normalize.
\param row_or_col   Either 'r' or 'c'.
\param normalization     see \ref apop_vector_normalize.

\ingroup basic_stats */
void apop_matrix_normalize(gsl_matrix *data, const char row_or_col, const char normalization){
  if (!data){
      apop_error(1, 'c', "input matrix is NULL. Doing nothing.\n");
      return; }
  gsl_vector  v;
  int         j;
        if (row_or_col == 'r')
            for (j = 0; j < data->size1; j++){
                v	= gsl_matrix_row(data, j).vector;
                apop_vector_normalize(&v, NULL, normalization);
            }
        else
            for (j = 0; j < data->size2; j++){
                v	= gsl_matrix_column(data, j).vector;
                apop_vector_normalize(&v, NULL, normalization);
            }
}

/** Input: any old vector. Output: 1 - the p-value for a chi-squared test to answer the question, "with what confidence can I reject the hypothesis that the variance of my data is zero?"

\param in a gsl_vector of data.
\ingroup asst_tests
 */
inline double apop_test_chi_squared_var_not_zero(const gsl_vector *in){
  apop_assert(in, 0, 0, 'c', "input vector is NULL. Doing nothing.\n");
  gsl_vector	*normed;
  int		    i;
  double 		sum=0;
	apop_vector_normalize((gsl_vector *)in,&normed, 1);
	gsl_vector_mul(normed,normed);
	for(i=0;i< normed->size; 
			sum +=gsl_vector_get(normed,i++));
	gsl_vector_free(normed);
	return gsl_cdf_chisq_P(sum,in->size); 
}

/** give me a random double between min and max [inclusive].

\param min, max 	Y'know.
\param r		The random number generator you're using.
*/
double apop_random_double(double min, double max, gsl_rng *r){
  double		base = gsl_rng_uniform(r);
	return base * (max - min) - min;
}

/** give me a random integer between min and max [inclusive].

\param min, max 	Y'know.
\param r		The random number generator you're using.
*/
int apop_random_int(const double min, const double max, const gsl_rng *r){
  double		base = gsl_rng_uniform(r);
	return (int) (base * (max - min + 1) - min);
}

/** Returns the sum of the elements of a matrix. Occasionally convenient.

  \param m	the matrix to be summed. 
\ingroup convenience_fns*/
  long double apop_matrix_sum(const gsl_matrix *m){
  int 		i,j;
  long double	sum	= 0;
	for (j=0; j< m->size1; j++)
		for (i=0; i< m->size2; i++)
			sum     += (gsl_matrix_get(m, j, i));
	return sum;
}

/** Returns the mean of all elements of a matrix.

  Calculated to avoid overflow errors.

  \param data	the matrix to be averaged. 
\ingroup convenience_fns*/
double apop_matrix_mean(const gsl_matrix *data){
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
double apop_matrix_var_m(const gsl_matrix *data, double mean){
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
\ingroup convenience_fns
*/
void apop_matrix_mean_and_var(const gsl_matrix *data, double *mean, double *var){
  long double avg     = 0,
              avg2    = 0;
  int         i,j, cnt= 0;
  long double x, ratio;
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
\return     An \ref apop_data structure with one row for each column in the original table, and a column for each summary statistic. May have a <tt>weights</tt> element.
\ingroup    output
\todo At the moment, only gives the mean, standard deviation, and variance
of the data in each column; should give more in the near future.
\todo We should probably let this summarize rows as well.
*/
apop_data * apop_data_summarize(apop_data *indata){
  apop_assert(indata, NULL, 0, 'c', "You sent me a NULL apop_data set. Returning NULL.\n");
  apop_assert(indata->matrix, NULL, 0, 'c', "You sent me an apop_data set with a NULL matrix. Returning NULL.\n");
  int		    i;
  apop_data	*out	= apop_data_alloc(0,indata->matrix->size2, 3);
  double		mean, stddev,var;
  char		rowname[10000]; //crashes on more than 10^9995 columns.
	apop_name_add(out->names, "mean", 'c');
	apop_name_add(out->names, "std dev", 'c');
	apop_name_add(out->names, "variance", 'c');
	if (indata->names !=NULL)
		for (i=0; i< indata->names->colct; i++)
			apop_name_add(out->names, indata->names->column[i], 'r');
	else
		for (i=0; i< indata->matrix->size2; i++){
			sprintf(rowname, "col %i", i);
			apop_name_add(out->names, rowname, 'r');
		}
	for (i=0; i< indata->matrix->size2; i++){
        APOP_MATRIX_COL(indata->matrix, i, v);
        if (indata->weights){
            mean	= apop_vector_mean(v);
            var 	= apop_vector_var_m(v,mean);
            stddev	= sqrt(var);
        } else{
            mean	= apop_vector_weighted_mean(v,indata->weights);
            var 	= apop_vector_weighted_var(v,indata->weights);
            stddev	= sqrt(var);
        } 
		gsl_matrix_set(out->matrix, i, 0, mean);
		gsl_matrix_set(out->matrix, i, 1, stddev);
		gsl_matrix_set(out->matrix, i, 2, var);
	}	
	return out;
}

/** Put summary information about the columns of a table (mean, std dev, variance) in a table.

 This is just the version of \ref apop_data_summarize for when
 you have a gsl_matrix instead of an \ref apop_data set. In
 fact, here's the source code for this function: <tt>return
 apop_data_summarize(apop_matrix_to_data(m));</tt>

 */
apop_data * apop_matrix_summarize(gsl_matrix *m){
  apop_assert(m,  NULL, 0, 'c', "You sent me a NULL gsl_matrix. Returning NULL.\n");
    return apop_data_summarize(apop_matrix_to_data(m));
}


/** Find the weighted mean. 

\param  v   The data vector
\param  w   the weight vector. If NULL, assume equal weights.
\return     The weighted mean
*/
double apop_vector_weighted_mean(const gsl_vector *v,const  gsl_vector *w){
  int           i;
  long double   sum = 0, wsum = 0;
    if (!w)
        return apop_vector_mean(v);
    apop_assert(v,  0, 0, 'c', "data vector is NULL. Returning zero.\n");
    apop_assert(v->size,  0, 1, 'c', "data vector has size 0. Returning zero.\n");
    apop_assert(w->size == v->size,  0, 0,'c', "data vector has size %u; weighting vector has size %u. Returning zero.\n", v->size, w->size);
    for (i=0; i< w->size; i++){
        sum  += gsl_vector_get(w, i) * gsl_vector_get(v,i); 
        wsum += gsl_vector_get(w, i); 
    }
    return sum/wsum;
}

/** Find the sample variance of a weighted vector.

Note: Apophenia tries to be smart about reading the weights. If weights
sum to one, then the system uses \c w->size as the number of elements,
and returns the usual sum over \f$n-1\f$. If weights > 1, then the
system uses the total weights as \f$n\f$. Thus, you can use the weights
as standard weightings or to represent elements that appear repeatedly.

\param  v   The data vector
\param  w   the weight vector. If NULL, assume equal weights.
\return     The weighted sample variance. 
*/
double apop_vector_weighted_var(const gsl_vector *v, const gsl_vector *w){
  int           i;
  long double   sum = 0, wsum = 0, sumsq = 0, vv, ww;
    if (!w)
        return apop_vector_var(v);
    apop_assert(v,  0, 0, 'c', "data vector is NULL. Returning zero.\n");
    apop_assert(v->size,  0, 0,'c', "data vector has size 0. Returning zero.\n");
    apop_assert(w->size == v->size,  0, 0,'c', "data vector has size %u; weighting vector has size %u. Returning zero.\n", v->size, w->size);
    //Using the E(x^2) - E^2(x) form.
    for (i=0; i< w->size; i++){
        vv    = gsl_vector_get(v,i);
        ww    = gsl_vector_get(w,i);
        sum  += ww * vv;
        sumsq+= ww * gsl_pow_2(vv); 
        wsum += ww; 
    }
    double len = (wsum < 1.1 ? w->size : wsum);
    return (sumsq/len - gsl_pow_2(sum/len)) * len/(len -1.);
}

static double skewkurt(const gsl_vector *v, const gsl_vector *w, const int exponent, const char *fn_name){
  int           i;
  long double   wsum = 0, sumcu = 0, vv, ww, mu;
    if (!w)
        return exponent ==3 ? apop_vector_skew(v) : apop_vector_kurtosis(v);
    apop_assert(v,  0, 0, 'c', "%s: data vector is NULL. Returning zero.\n", fn_name);
    apop_assert(v->size,  0, 1, 'c',"%s: data vector has size 0. Returning zero.\n", fn_name);
    apop_assert(w->size == v->size,  0, 1, 'c',"%s: data vector has size %i; weighting vector has size %i. Returning zero.\n", fn_name, v->size, w->size);
    //Using the E(x - \bar x)^3 form, which is lazy.
    mu  = apop_vector_weighted_mean(v, w);
    for (i=0; i< w->size; i++){
        vv    = gsl_vector_get(v,i);
        ww    = gsl_vector_get(w,i);
        sumcu+= ww * gsl_pow_int(vv - mu, exponent); 
        wsum += ww; 
    }
    double len = wsum < 1.1 ? w->size : wsum;
    return sumcu/(len -1);


}

/** Find the sample skew of a weighted vector.

Note: Apophenia tries to be smart about reading the weights. If weights
sum to one, then the system uses \c w->size as the number of elements,
and returns the usual sum over \f$n-1\f$. If weights > 1, then the
system uses the total weights as \f$n\f$. Thus, you can use the weights
as standard weightings or to represent elements that appear repeatedly.

\param  v   The data vector
\param  w   the weight vector. If NULL, assume equal weights.
\return     The weighted skew. No sample adjustment given weights.
\todo   \c apop_vector_weighted_skew and \c apop_vector_weighted_kurt are lazily written.
*/
double apop_vector_weighted_skew(const gsl_vector *v, const gsl_vector *w){
    return skewkurt(v,w,3, "apop_vector_weighted_skew");
}

/** Find the sample kurtosis of a weighted vector.

\param  v   The data vector
\param  w   the weight vector. If NULL, assume equal weights.
\return     The weighted kurtosis. No sample adjustment given weights.
\todo   \c apop_vector_weighted_skew and \c apop_vector_weighted_kurt are lazily written.
*/
double apop_vector_weighted_kurt(const gsl_vector *v, const gsl_vector *w){
    return skewkurt(v,w,4, "apop_vector_weighted_kurt");
}

/** Find the sample covariance of a pair of weighted vectors. This only
makes sense if the weightings are identical, so the function takes only one weighting vector for both.

\param  v1, v2   The data vectors
\param  w   the weight vector. If NULL, assume equal weights.
\return     The weighted sample covariance
*/
double apop_vector_weighted_cov(const gsl_vector *v1, const gsl_vector *v2, const gsl_vector *w){
  int           i;
  long double   sum1 = 0, sum2 = 0, wsum = 0, sumsq = 0, vv1, vv2, ww;
    if (!w)
        return apop_vector_cov(v1,v2);
    apop_assert(v1,  0, 0, 'c', "first data vector is NULL. Returning zero.\n");
    apop_assert(v2,  0, 0, 'c', "second data vector is NULL. Returning zero.\n");
    apop_assert(v1->size,  0, 1, 'c', "apop_vector_weighted_variance: data vector has size 0. Returning zero.\n");
    apop_assert((w->size == v1->size) && (w->size == v2->size),  0, 0, 'c', "apop_vector_weighted_variance: data vectors have sizes %i and %i; weighting vector has size %i. Returning zero.\n", v1->size, v2->size, w->size);
    //Using the E(x^2) - E^2(x) form.
    for (i=0; i< w->size; i++){
        vv1   = gsl_vector_get(v1,i);
        vv2   = gsl_vector_get(v2,i);
        ww    = gsl_vector_get(w,i);
        sum1 += ww * vv1;
        sum2 += ww * vv2;
        sumsq+= ww * vv1 * vv2;
        wsum += ww; 
    }
    double len = (wsum < 1.1 ? w->size : wsum);
    return (sumsq/len  - sum1*sum2/gsl_pow_2(len)) *(len/(len-1));
}


/** Returns the variance/covariance matrix relating each column with each other.

This is the \c gsl_matrix  version of \ref apop_data_covariance; if you have column names, use that one.

\param in 	A data matrix: rows are observations, columns are variables.
\param normalize
'n', 'N', or 1 = subtract the mean from each column, thus changing the input data but speeding up the computation.<br>
anything else (like 0)= don't modify the input data

\return Returns the variance/covariance matrix relating each column with each other. This function allocates the matrix for you.
This is the sample version---dividing by \f$n-1\f$, not \f$n\f$.
\ingroup matrix_moments */
gsl_matrix *apop_matrix_covariance(gsl_matrix *in, const char normalize){
  apop_assert(in, NULL, 0, 'c', "input matrix is NULL. Returning NULL.\n");
  gsl_matrix	*out;
  int		i,j;
  double		means[in->size2];
  gsl_vector_view	v, v1, v2;
	if (normalize == 1 || normalize=='n' || normalize == 'N'){
		out	= gsl_matrix_alloc(in->size2, in->size2);
		apop_matrix_normalize(in,'c', 0);
		gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, in, in, 0, out);
	    gsl_matrix_scale(out, 1.0/(in->size1-1));
	}
	else{
		out	= gsl_matrix_calloc(in->size2, in->size2);
		for(i=0; i< in->size2; i++){
			v		= gsl_matrix_column(in, i);
			means[i]	= apop_mean(&(v.vector));
		}
		for(i=0; i< in->size2; i++){
			v1		= gsl_matrix_column(in, i);
			for(j=i; j< in->size2; j++){
			    v2		= gsl_matrix_column(in, j);
                apop_matrix_increment(out, i, j, 
					   gsl_stats_covariance_m (v1.vector.data, v1.vector.stride, v2.vector.data,  v2.vector.stride,
                                                    v2.vector.size, means[i], means[j]));
				if (i != j)	//set the symmetric element.
					gsl_matrix_set(out, j, i, gsl_matrix_get(out, i, j));
			}
	    }
    }
	return out;
}

/** Returns the matrix of correlation coefficients (\f$\sigma^2_{xy}/(\sigma_x\sigma_y)\f$) relating each column with each other.

This is the \c gsl_matrix  version of \ref apop_data_covariance; if you have column names, use that one.

\param in 	A data matrix: rows are observations, columns are variables.
\param normalize
'n' or 'N' = subtract the mean from each column, thus changing the input data but speeding up the computation.<br>
anything else (like 0)= don't modify the input data

\return Returns the variance/covariance matrix relating each column with each other. This function allocates the matrix for you.
\ingroup matrix_moments */
gsl_matrix *apop_matrix_correlation(gsl_matrix *in, const char normalize){
  apop_assert(in,  NULL, 0, 'c', "input matrix is NULL. Returning NULL.\n");
  gsl_matrix      *out    = apop_matrix_covariance(in, normalize);
  int             i;
  double          std_dev;
    for(i=0; i< in->size2; i++){
        APOP_MATRIX_COL(in, i, cvin);
        APOP_MATRIX_COL(out, i, cvout);
        APOP_MATRIX_ROW(out, i, rvout);
        std_dev     = sqrt(apop_var(cvin));
        gsl_vector_scale(cvout, 1.0/std_dev);
        gsl_vector_scale(rvout, 1.0/std_dev);
    }
    return out;
}

/** Returns the variance/covariance matrix relating each column of the matrix to each other column.

This is the \ref apop_data version of \ref apop_matrix_covariance; if you don't have column names or weights, or would like to use the speed-saving and data-destroying normalization option, use that one.

\param in 	An \ref apop_data set. If the weights vector is set, I'll take it into account.

\return Returns a \ref apop_data set the variance/covariance matrix relating each column with each other.

\ingroup matrix_moments */
apop_data *apop_data_covariance(const apop_data *in){
  apop_assert(in,  NULL, 0, 'c', "You sent me a NULL apop_data set. Returning NULL.\n");
  apop_assert(in->matrix,  NULL, 0, 'c', "You sent me an apop_data set with a NULL matrix. Returning NULL.\n");
  apop_data   *out = apop_data_alloc(0,in->matrix->size2, in->matrix->size2);
  int         i, j;
  double      var;
    for (i=0; i < in->matrix->size2; i++){
        for (j=i; j < in->matrix->size2; j++){
            APOP_COL(in, i, v1);
            APOP_COL(in, j, v2);
            var = apop_vector_weighted_cov(v1, v2, in->weights);
            gsl_matrix_set(out->matrix, i,j, var);
            if (i!=j)
                gsl_matrix_set(out->matrix, j,i, var);
        }
    }
    apop_name_stack(out->names, in->names, 'c');
    apop_name_cross_stack(out->names, in->names, 'r', 'c');
    return out;
}

/** Returns the matrix of correlation coefficients (\f$\sigma^2_{xy}/(\sigma_x\sigma_y)\f$) relating each column with each other.

This is the \ref apop_data version of \ref apop_matrix_correlation; if you don't have column names or weights, (or want the option for the faster, data-destroying version), use that one.

\param in 	A data matrix: rows are observations, columns are variables. If you give me a weights vector, I'll use it.

\return Returns the variance/covariance matrix relating each column with each other. This function allocates the matrix for you.
\ingroup matrix_moments */
apop_data *apop_data_correlation(const apop_data *in){
  apop_data *out = apop_data_covariance(in);
  int       i;
  double    std_dev;
    for(i=0; i< in->matrix->size2; i++){
        APOP_COL(in, i, cvin);
        APOP_COL(out, i, cvout);
        APOP_ROW(out, i, rvout);
        std_dev     = sqrt(apop_vector_weighted_var(cvin,in->weights));
        gsl_vector_scale(cvout, 1.0/std_dev);
        gsl_vector_scale(rvout, 1.0/std_dev);
    }
    return out;
}
