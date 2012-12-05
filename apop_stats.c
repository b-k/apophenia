/** \file apop_stats.c	Basic moments and some distributions. */
/* Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_eigen.h>


/** \defgroup vector_moments Calculate moments (mean, var, kurtosis) for the data in a gsl_vector.

These functions simply take in a GSL vector and return its mean, variance, or kurtosis; the covariance functions take two GSL vectors as inputs.

\ref apop_vector_kurtosis and \ref apop_vector_kurt are identical; pick the one which sounds better to you.

\see db_moments

For \ref apop_vector_var_m<tt>(vector, mean)</tt>, <tt>mean</tt> is the mean of the
vector. This saves the trouble of re-calculating the mean if you've
already done so. E.g.,

\code
gsl_vector *v;
double mean, var;

//Allocate v and fill it with data here.
mean = apop_vector_mean(v);
var  = apop_vector_var_m(v, mean);
printf("Your vector has mean %g and variance %g\n", mean, var);
\endcode
*/

/** Returns the sum of the data in the given vector.
\ingroup convenience_fns
*/
long double apop_vector_sum(const gsl_vector *in){
    Apop_stopif(!in, return 0, 1, "You just asked me to sum a NULL. Returning zero.")
    long double out = 0;
    for (size_t i=0; i< in->size; i++)
        out += gsl_vector_get(in, i);
	return out; 
}

/** \def apop_sum(in)
  An alias for \ref apop_vector_sum. Returns the sum of the data in the given vector.
\ingroup convenience_fns
*/

/**  \def apop_vector_mean(in)
 
 Returns the mean of the data in the given vector.
\ingroup vector_moments
*/

/**  \def apop_mean(in)
  An alias for \ref apop_vector_mean.  Returns the mean of the data in the given vector.
\ingroup vector_moments
*/

/** \def apop_vector_var(in)
 Returns the variance of the data in the given vector.
 
  This uses (n-1) in the denominator of the sum; i.e., it corrects for the bias introduced by using \f$\bar x\f$ instead of \f$\mu\f$.

  At the moment, there is no var_pop function. Just multiply this by (n-1)/n if you need that.
\ingroup vector_moments
*/

/** \def apop_var(in)
  An alias for \ref apop_vector_var.
Returns the variance of the data in the given vector.
\ingroup vector_moments
*/

/** Returns an unbiased estimate of the sample skew (population skew times  y \f$n^2/(n^2-1)\f$) of the data in the given vector.
\ingroup vector_moments
*/
double apop_vector_skew(const gsl_vector *in){
	return apop_vector_skew_pop(in) * gsl_pow_2(in->size)/((in->size -1.)*(in->size -2.)); }

/** Returns the population skew \f$(\sum_i (x_i - \mu)^3/n))\f$ of the data in the given vector.
 
  Some people like to normalize the skew by dividing by variance\f$^{3/2}\f$; that's not done here, so you'll have to do so separately if need be.
\ingroup vector_moments
*/
double apop_vector_skew_pop(const gsl_vector *in){
  //This code is cut/pasted/modified from the GSL. 
  //I reimplement the skew calculation here without the division by var^3/2 that the GSL does. 

    size_t n = in->size;
    long double avg = 0;
    long double mean = apop_vector_mean(in);
    for (size_t i = 0; i < n; i++) {
        const long double x = gsl_vector_get(in, i) - mean; 
        avg += (x * x * x  - avg)/(i + 1);
    } 
    return avg;
}

/** Returns the population kurtosis (\f$\sum_i (x_i - \mu)^4/n)\f$) of the data in the given vector.
 
  Some people like to normalize the kurtosis by dividing by variance squared, or by subtracting three; those things are  not done here, so you'll have to do them separately if need be.
\ingroup vector_moments
*/
double apop_vector_kurtosis_pop(const gsl_vector *in){
  //This code is cut/pasted/modified from the GSL. 
  //I reimplement the kurtosis calculation here without the division by var^2 that the GSL does. 

    size_t n = in->size;
    long double avg  = 0;
    long double mean = apop_vector_mean(in);
    for (size_t i = 0; i < n; i++) {
        const long double x = gsl_vector_get(in, i) - mean; 
        avg += (x * x * x * x - avg)/(i + 1);
    } 
    return avg;
}


/** Returns the sample kurtosis (divide by \f$n-1\f$) of the data in the given vector. Corrections are made to produce an unbiased result.

  \li This does not normalize the output: the kurtosis of a \f${\cal N}(0,1)\f$ is three \f$\sigma^4\f$, not three, one, or zero.
\ingroup vector_moments
*/
double apop_vector_kurtosis(const gsl_vector *in){
    size_t n = in->size;
    long double coeff0= n*n/(gsl_pow_3(n)*(gsl_pow_2(n)-3*n+3));
    long double coeff1= n*gsl_pow_2(n-1)+ (6*n-9);
    long double coeff2= n*(6*n-9);
    return  coeff0 *(coeff1 * apop_vector_kurtosis_pop(in) + coeff2 * gsl_pow_2(apop_vector_var(in)*(n-1.)/n));
}

/** \def apop_vector_kurt(in)
  Returns the sample kurtosis (divide by \f$n-1\f$) of the data in the given vector.
  This does not normalize the output: the kurtosis of a \f${\cal N}(0,1)\f$ is three, not zero.

  An alias for \ref apop_vector_kurtosis.
\ingroup vector_moments
*/

/** Returns the variance of the data in the given vector, given that you've already calculated the mean.
\param in	the vector in question
\param mean	the mean, which you've already calculated using \ref apop_vector_mean.
\ingroup vector_moments
*/
double apop_vector_var_m(const gsl_vector *in, const double mean){
	return gsl_stats_variance_m(in->data,in->stride, in->size, mean); }

/** Returns the covariance of two vectors
\ingroup vector_moments
*/
double apop_vector_cov(const gsl_vector *ina, const gsl_vector *inb){
	return gsl_stats_covariance(ina->data,ina->stride,inb->data,inb->stride,inb->size); }

/** Returns the correlation coefficient of two vectors. It's just
\f$ {\hbox{cov}(a,b)\over \sqrt(\hbox{var}(a)) * \sqrt(\hbox{var}(b))}.\f$
\ingroup vector_moments
*/
double apop_vector_correlation(const gsl_vector *ina, const gsl_vector *inb){
	return apop_vector_cov(ina, inb) / sqrt(apop_vector_var(ina) * apop_vector_var(inb)); }


/** Returns the distance between two vectors, where distance is defined
 based on the third (optional) parameter:

 - 'e' or 'E' (the default): scalar distance (standard Euclidean metric) between two vectors. Simply \f$\sqrt{\sum_i{(a_i - b_i)^2}},\f$
where \f$i\f$ iterates over dimensions.
 - 'm' or 'M'  Returns the Manhattan metric distance  between two vectors: \f$\sum_i{|a_i - b_i|},\f$
where \f$i\f$ iterates over dimensions.
 - 'd' or 'D' The discrete norm: if \f$a = b\f$, return zero, else return one.
 - 's' or 'S' The sup norm: find the dimension where \f$|a_i - b_i|\f$ is largest, return the distance along that one dimension.
 - 'l' or 'L' The \f$L_p\f$ norm, \f$\left(\sum_i{(a_i - b_i)^2}\right)^{1/p},\f$. The value of \f$p\f$ is set by the fourth (optional) argument.

 \param ina First vector (No default, must not be \c NULL)
 \param inb Second vector (Default = zero)
 \param metric The type of metric, as above.
 \param norm  If you are using an \f$L_p\f$ norm, this is \f$p\f$. Must be strictly greater than zero. (default = 2)

 Notice that the defaults are such that
 \code
 apop_vector_distance(v);
 apop_vector_distance(v, .metric = 's');
 \endcode
 gives you the standard Euclidean length of \c v and its longest element.

\include test_distances.c

This function uses the \ref designated syntax for inputs.
\ingroup convenience_fns
*/
APOP_VAR_HEAD double apop_vector_distance(const gsl_vector *ina, const gsl_vector *inb, const char metric, const double norm){
    static gsl_vector *zero = NULL;
    const gsl_vector * apop_varad_var(ina, NULL);
    Apop_assert(ina, "The first vector has to be non-NULL.");
    const gsl_vector * apop_varad_var(inb, NULL);
    if (!inb){
        if (!zero || zero->size !=ina->size){
            if (zero) gsl_vector_free(zero);
            zero = gsl_vector_calloc(ina->size);
        }
        inb = zero;
    }
    const char apop_varad_var(metric, 'e');
    const double apop_varad_var(norm, 2);
APOP_VAR_ENDHEAD
    Apop_assert(ina->size == inb->size, "I need equal-sized vectors, but "
                "you sent a vector of size %zu and a vector of size %zu. ", ina->size, inb->size);
    double  dist    = 0;
    size_t  i;
    if (metric == 'e' || metric == 'E'){
        for (i=0; i< ina->size; i++){
            dist    += gsl_pow_2(gsl_vector_get(ina, i) - gsl_vector_get(inb, i));
        }
        return sqrt(dist); 
    }
    if (metric == 'm' || metric == 'M'){ //redundant with vector_grid_distance, below.
        for (i=0; i< ina->size; i++) 
            dist    += fabs(gsl_vector_get(ina, i) - gsl_vector_get(inb, i));
        return dist; 
    }
    if (metric == 'd' || metric == 'D'){
        for (i=0; i< ina->size; i++) 
            if (gsl_vector_get(ina, i) != gsl_vector_get(inb, i))
                return 1;
        return 0;
    }
    if (metric == 's' || metric == 'S'){
        for (i=0; i< ina->size; i++) 
            dist = GSL_MAX(dist, fabs(gsl_vector_get(ina, i) - gsl_vector_get(inb, i)));
        return dist;
    }
    if (metric == 'l' || metric == 'L'){
        for (i=0; i< ina->size; i++)
            dist += pow(fabs(gsl_vector_get(ina, i) - gsl_vector_get(inb, i)), norm);
        return pow(dist, 1./norm); 
    }
  Apop_assert(0, "I couldn't find the metric type you gave, %c, in my list of supported types.", metric);
}

/** This function will normalize a vector, either such that it has mean
zero and variance one, or such that it ranges between zero and one, or sums to one.

\param in 	A gsl_vector which you have already allocated and filled. \c NULL input gives \c NULL output. (No default)

\param out 	If normalizing in place, \c NULL.
If not, the address of a <tt>gsl_vector</tt>. Do not allocate. (default = \c NULL.)

\param normalization_type 
'p': normalized vector will sum to one. E.g., start with a set of observations in bins, end with the percentage of observations in each bin. (the default)<br>
'r': normalized vector will range between zero and one. Replace each X with (X-min) / (max - min).<br>
's': normalized vector will have mean zero and variance one. Replace
each X with \f$(X-\mu) / \sigma\f$, where \f$\sigma\f$ is the sample
standard deviation.<br>
'm': normalize to mean zero: Replace each X with \f$(X-\mu)\f$<br>

\b Example 
\code
#include <apop.h>

int main(void){
gsl_vector  *in, *out;

in = gsl_vector_calloc(3);
gsl_vector_set(in, 1, 1);
gsl_vector_set(in, 2, 2);

printf("The original vector:\n");
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

This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD void apop_vector_normalize(gsl_vector *in, gsl_vector **out, const char normalization_type){
    gsl_vector * apop_varad_var(in, NULL);
    Apop_assert_c(in,  , 1, "Input vector is NULL. Doing nothing.");
    gsl_vector ** apop_varad_var(out, NULL);
    const char apop_varad_var(normalization_type, 'p');
APOP_VAR_END_HEAD
    double mu, min, max;
	if (!out) out = &in;
	else {
		*out = gsl_vector_alloc (in->size);
		gsl_vector_memcpy(*out,in);
	}
	if ((normalization_type == 's')){
		mu	= apop_vector_mean(in);
		gsl_vector_add_constant(*out, -mu);
		gsl_vector_scale(*out, 1/(sqrt(apop_vector_var_m(in, mu))));
	} 
	else if ((normalization_type == 'r')){
        gsl_vector_minmax(in, &min, &max);
		gsl_vector_add_constant(*out, -min);
		gsl_vector_scale(*out, 1/(max-min));	

	}
	else if ((normalization_type == 'p')){
		long double sum	= apop_sum(in);
        Apop_assert_n(sum, "The vector sums to zero, so I can't normalize it to sum to one.");
		gsl_vector_scale(*out, 1/sum);	
	}
	else if ((normalization_type == 'm')){
		mu	= apop_vector_mean(in);
		gsl_vector_add_constant(*out, -mu);
	}
}

/** Normalize  each row or column in the given matrix, one by one.

  Basically just a convenience fn to iterate through the columns or rows and run \ref apop_vector_normalize for you.

\param data     The data set to normalize.
\param row_or_col   Either 'r' or 'c'.
\param normalization     see \ref apop_vector_normalize.
*/
void apop_matrix_normalize(gsl_matrix *data, const char row_or_col, const char normalization){
    Apop_assert_c(data, , 1, "input matrix is NULL. Doing nothing.");
    if (row_or_col == 'r')
        for (size_t j = 0; j < data->size1; j++){
            Apop_matrix_row(data, j, v);
            apop_vector_normalize(v, NULL, normalization);
        }
    else
        for (size_t j = 0; j < data->size2; j++){
            Apop_matrix_col(data, j, v);
            apop_vector_normalize(v, NULL, normalization);
        }
}

/** Returns the sum of the elements of a matrix. Occasionally convenient.

  \param m	the matrix to be summed. 
\ingroup convenience_fns*/
long double apop_matrix_sum(const gsl_matrix *m){
    Apop_assert_c(m, 0, 1, "You just asked me to sum a NULL. Returning zero.")
    long double	sum	= 0;
	for (size_t j=0; j< m->size1; j++)
		for (size_t i=0; i< m->size2; i++)
			sum     += gsl_matrix_get(m, j, i);
	return sum;
}

/** Returns the mean of all elements of a matrix.

  Calculated to avoid overflow errors.

  \param data	the matrix to be averaged. 
\ingroup convenience_fns*/
double apop_matrix_mean(const gsl_matrix *data){
    double  avg     = 0;
    int     cnt= 0;
    double  x, ratio;
    for(size_t i=0; i < data->size1; i++)
        for(size_t j=0; j < data->size2; j++){
            x       = gsl_matrix_get(data, i,j);
            ratio   = cnt/(cnt+1.0);
            cnt     ++;
            avg     *= ratio;
            avg     += x/(cnt +0.0);
        }
	return avg;
}

/** Returns the mean and population variance of all elements of a matrix.
 
  \li If you want sample variance, multiply the result returned by <tt>(data->size1*data->size2)/(data->size1*data->size2-1.0)</tt>.

  \param data	the matrix to be averaged. 
\param	mean	where to put the mean to be calculated
\param	var	where to put the variance to be calculated
\ingroup convenience_fns
*/
void apop_matrix_mean_and_var(const gsl_matrix *data, double *mean, double *var){
    long double avg     = 0,
                avg2    = 0;
    size_t cnt= 0;
    long double x, ratio;
    for(size_t i=0; i < data->size1; i++)
        for(size_t j=0; j < data->size2; j++){
            x     = gsl_matrix_get(data, i,j);
            ratio = cnt/(cnt+1.0);
            cnt   ++;
            avg   *= ratio;
            avg2  *= ratio;
            avg   += x/(cnt +0.0);
            avg2  += gsl_pow_2(x)/(cnt +0.0);
        }
	*mean = avg;
    *var  = avg2 - gsl_pow_2(avg); //E[x^2] - E^2[x]
}

/** Put summary information about the columns of a table (mean, std dev, variance, min, median, max) in a table.

\param indata The table to be summarized. An \ref apop_data structure.
\return     An \ref apop_data structure with one row for each column in the original table, and a column for each summary statistic. May have a <tt>weights</tt> element.
\exception out->error='a'  Allocation error.

\li This function gives more columns than you probably want; use \ref apop_data_prune_columns to pick the ones you want to see.
\todo We should probably let this summarize rows as well. 
\ingroup    output */
apop_data * apop_data_summarize(apop_data *indata){
    Apop_assert_c(indata, NULL, 0, "You sent me a NULL apop_data set. Returning NULL.");
    Apop_assert_c(indata->matrix, NULL, 0, "You sent me an apop_data set with a NULL matrix. Returning NULL.");
    apop_data *out = apop_data_alloc(indata->matrix->size2, 6);
    double mean, var;
    char rowname[10000]; //crashes on more than 10^9995 columns.
	apop_name_add(out->names, "mean", 'c');
	apop_name_add(out->names, "std dev", 'c');
	apop_name_add(out->names, "variance", 'c');
	apop_name_add(out->names, "min", 'c');
	apop_name_add(out->names, "median", 'c');
	apop_name_add(out->names, "max", 'c');
	if (indata->names !=NULL)
        apop_name_stack(out->names,indata->names, 'r', 'c');
	else
		for (size_t i=0; i< indata->matrix->size2; i++){
			sprintf(rowname, "col %zu", i);
			apop_name_add(out->names, rowname, 'r');
		}
	for (size_t i=0; i< indata->matrix->size2; i++){
        Apop_matrix_col(indata->matrix, i, v);
        if (!indata->weights){
            mean	= apop_vector_mean(v);
            var 	= apop_vector_var_m(v,mean);
        } else {
            mean	= apop_vector_weighted_mean(v,indata->weights);
            var 	= apop_vector_weighted_var(v,indata->weights);
        } 
        double *pctiles = apop_vector_percentiles(v);
		gsl_matrix_set(out->matrix, i, 0, mean);
		gsl_matrix_set(out->matrix, i, 1, sqrt(var));
		gsl_matrix_set(out->matrix, i, 2, var);
		gsl_matrix_set(out->matrix, i, 3, pctiles[0]);
		gsl_matrix_set(out->matrix, i, 4, pctiles[50]);
		gsl_matrix_set(out->matrix, i, 5, pctiles[100]);
        free(pctiles);
	}
	return out;
}

/** Find the weighted mean. 

\param  v   The data vector
\param  w   the weight vector. If NULL, assume equal weights.
\return     The weighted mean */
double apop_vector_weighted_mean(const gsl_vector *v,const  gsl_vector *w){
    long double   sum = 0, wsum = 0;
    if (!w) return apop_vector_mean(v);
    Apop_assert_c(v,  0, 1, "data vector is NULL. Returning zero.");
    Apop_assert_c(v->size,  0, 1, "data vector has size 0. Returning zero.");
    Apop_assert_c(w->size == v->size,  0, 0, "data vector has size %zu; weighting vector has size %zu. Returning zero.", v->size, w->size);
    for (size_t i=0; i< w->size; i++){
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
\return     The weighted sample variance.  */
double apop_vector_weighted_var(const gsl_vector *v, const gsl_vector *w){
    if (!w) return apop_vector_var(v);
    long double sum = 0, wsum = 0, sumsq = 0, vv, ww;
    apop_assert_c(v,  0, 1, "data vector is NULL. Returning zero.\n");
    apop_assert_c(v->size, 0, 1, "data vector has size 0. Returning zero.\n");
    apop_assert_c(w->size == v->size, GSL_NAN, 0, "data vector has size %zu; weighting vector has size %zu. Returning NaN.\n", v->size, w->size);
    //Using the E(x^2) - E^2(x) form.
    for (size_t i=0; i< w->size; i++){
        vv    = gsl_vector_get(v,i);
        ww    = gsl_vector_get(w,i);
        sum  += ww * vv;
        sumsq+= ww * gsl_pow_2(vv); 
        wsum += ww; 
    }
    double len = (wsum < 1.1 ? w->size : wsum);
    return (sumsq/len - gsl_pow_2(sum/len)) * len/(len -1.);
}

static double wskewkurt(const gsl_vector *v, const gsl_vector *w, const int exponent, const char *fn_name){
    long double   wsum = 0, sumcu = 0, vv, ww, mu;
    Apop_assert_c(v,  0, 1, "%s: data vector is NULL. Returning zero.", fn_name);
    Apop_assert_c(v->size,  0, 1,"%s: data vector has size 0. Returning zero.", fn_name);
    Apop_assert_c(w->size == v->size,  0, 1,"%s: data vector has size %zu; weighting vector has size %zu. Returning zero.", fn_name, v->size, w->size);
    //Using the E(x - \bar x)^3 form, which is lazy.
    mu  = apop_vector_weighted_mean(v, w);
    for (size_t i=0; i< w->size; i++){
        vv    = gsl_vector_get(v,i);
        ww    = gsl_vector_get(w,i);
        sumcu+= ww * gsl_pow_int(vv - mu, exponent); 
        wsum += ww; 
    }
    double len = wsum < 1.1 ? w->size : wsum;
    return sumcu/len;
}

/** Find the population skew of a weighted vector.

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
    return wskewkurt(v,w,3, "apop_vector_weighted_skew");
}

/** Find the population kurtosis of a weighted vector.

\param  v   The data vector
\param  w   the weight vector. If NULL, assume equal weights.
\return     The weighted kurtosis. No sample adjustment given weights.
\todo   \c apop_vector_weighted_skew and \c apop_vector_weighted_kurt are lazily written.
*/
double apop_vector_weighted_kurt(const gsl_vector *v, const gsl_vector *w){
    return wskewkurt(v,w,4, "apop_vector_weighted_kurt");
}

/** Find the sample covariance of a pair of weighted vectors. This only
makes sense if the weightings are identical, so the function takes only one weighting vector for both.

\param  v1, v2   The data vectors
\param  w   the weight vector. If NULL, assume equal weights.
\return     The weighted sample covariance
*/
double apop_vector_weighted_cov(const gsl_vector *v1, const gsl_vector *v2, const gsl_vector *w){
    if (!w) return apop_vector_cov(v1,v2);
    long double sum1 = 0, sum2 = 0, wsum = 0, sumsq = 0, vv1, vv2, ww;
    Apop_assert_c(v1,  0, 1, "first data vector is NULL. Returning zero.");
    Apop_assert_c(v2,  0, 1, "second data vector is NULL. Returning zero.");
    Apop_assert_c(v1->size,  0, 1, "apop_vector_weighted_variance: data vector has size 0. Returning zero.");
    Apop_assert_c((w->size == v1->size) && (w->size == v2->size), GSL_NAN, 0, "apop_vector_weighted_variance: data vectors have sizes %zu and %zu; weighting vector has size %zu. Returning NaN.", v1->size, v2->size, w->size);
    //Using the E(x^2) - E^2(x) form.
    for (size_t i=0; i< w->size; i++){
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

\param in 	A data matrix: rows are observations, columns are variables. (No default, must not be \c NULL)
\param normalize
'n', 'N', or 1 = subtract the mean from each column, thus changing the input data but speeding up the computation.<br>
anything else (like 0)= don't modify the input data (default = no modification)

\return Returns the variance/covariance matrix relating each column with each other. This function allocates the matrix for you.
This is the sample version---dividing by \f$n-1\f$, not \f$n\f$.
It uses the \ref designated syntax for inputs.
\ingroup matrix_moments */
APOP_VAR_HEAD gsl_matrix *apop_matrix_covariance(gsl_matrix *in, const char normalize){
    gsl_matrix *apop_varad_var(in, NULL)
    Apop_assert_c(in,  NULL, 0, "Input matrix is NULL. Returning same.");
    const char apop_varad_var(normalize, 0)
APOP_VAR_ENDHEAD
    gsl_matrix *out;
    double means[in->size2];
    gsl_vector_view	v1, v2;
	if (normalize == 1 || normalize=='n' || normalize == 'N'){
		out	= gsl_matrix_alloc(in->size2, in->size2);
		apop_matrix_normalize(in,'c', 0);
		gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, in, in, 0, out);
	    gsl_matrix_scale(out, 1.0/(in->size1-1));
	}
	else{
		out	= gsl_matrix_calloc(in->size2, in->size2);
		for(size_t i=0; i< in->size2; i++){
			Apop_matrix_col(in, i, v);
			means[i]	= apop_mean(v);
		}
		for(size_t i=0; i< in->size2; i++){
			v1		= gsl_matrix_column(in, i);
			for(size_t j=i; j< in->size2; j++){
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

/** Returns the matrix of correlation coefficients \f$(\sigma^2_{xy}/(\sigma_x\sigma_y))\f$ relating each column with each other.

This is the \c gsl_matrix  version of \ref apop_data_correlation; if you have column names, use that one.

\param in 	A data matrix: rows are observations, columns are variables. (No default, must not be \c NULL)
\param normalize
'n' or 'N' = subtract the mean from each column, thus changing the input data but speeding up the computation.<br>
anything else (like 0)= don't modify the input data (default = no modification)

\return Returns the variance/covariance matrix relating each column with each other. This function allocates the matrix for you.

This function uses the \ref designated syntax for inputs.
\ingroup matrix_moments */
APOP_VAR_HEAD gsl_matrix *apop_matrix_correlation(gsl_matrix *in, const char normalize){
    gsl_matrix *apop_varad_var(in, NULL)
    apop_assert_c(in,  NULL, 0, "Input matrix is NULL; returning NULL.");
    const char apop_varad_var(normalize, 0)
APOP_VAR_ENDHEAD
    gsl_matrix *out = apop_matrix_covariance(in, normalize);
    for(size_t i=0; i< in->size2; i++){
        APOP_MATRIX_COL(in, i, cvin);
        APOP_MATRIX_COL(out, i, cvout);
        APOP_MATRIX_ROW(out, i, rvout);
        double std_dev = sqrt(apop_var(cvin));
        gsl_vector_scale(cvout, 1.0/std_dev);
        gsl_vector_scale(rvout, 1.0/std_dev);
    }
    return out;
}

/** Returns the variance/covariance matrix relating each column of the matrix to each other column.

This is the \ref apop_data version of \ref apop_matrix_covariance; if you don't have column names or weights, or would like to use the speed-saving and data-destroying normalization option, use that one.

\param in 	An \ref apop_data set. If the weights vector is set, I'll take it into account.

\return Returns a \ref apop_data set the variance/covariance matrix relating each column with each other.
\exception out->error='a'  Allocation error.
\ingroup matrix_moments */
apop_data *apop_data_covariance(const apop_data *in){
    Apop_assert_c(in,  NULL, 1, "You sent me a NULL apop_data set. Returning NULL.");
    Apop_assert_c(in->matrix,  NULL, 1, "You sent me an apop_data set with a NULL matrix. Returning NULL.");
    apop_data *out = apop_data_alloc(in->matrix->size2, in->matrix->size2);
    Apop_stopif(out->error, return out, 0, "allocation error.");
    for (size_t i=0; i < in->matrix->size2; i++){
        for (size_t j=i; j < in->matrix->size2; j++){
            Apop_col(in, i, v1);
            Apop_col(in, j, v2);
            double var = apop_vector_weighted_cov(v1, v2, in->weights);
            gsl_matrix_set(out->matrix, i,j, var);
            if (i!=j) gsl_matrix_set(out->matrix, j,i, var);
        }
    }
    apop_name_stack(out->names, in->names, 'c');
    apop_name_stack(out->names, in->names, 'r', 'c');
    return out;
}

/** Returns the matrix of correlation coefficients \f$(\sigma^2_{xy}/(\sigma_x\sigma_y))\f$ relating each column with each other.

This is the \ref apop_data version of \ref apop_matrix_correlation; if you don't have column names or weights, (or want the option for the faster, data-destroying version), use that one.

\param in 	A data matrix: rows are observations, columns are variables. If you give me a weights vector, I'll use it.

\return Returns the variance/covariance matrix relating each column with each other. This function allocates the matrix for you.
\exception out->error='a'  Allocation error.
\ingroup matrix_moments */
apop_data *apop_data_correlation(const apop_data *in){
    apop_data *out = apop_data_covariance(in);
    for(size_t i=0; i< in->matrix->size2; i++){
        Apop_col(in, i, cvin);
        Apop_col(out, i, cvout);
        Apop_row(out, i, rvout);
        double std_dev = sqrt(apop_vector_weighted_var(cvin, in->weights));
        gsl_vector_scale(cvout, 1.0/std_dev);
        gsl_vector_scale(rvout, 1.0/std_dev);
    }
    return out;
}

static void get_one_row(apop_data *p, apop_data *a_row, int i, int min, int max){
    for (int j=min; j< max; j++)
        apop_data_set(a_row, 0, j-min, apop_data_get(p, i, j));
}

/** Kullback-Leibler divergence.

  This measure of the divergence of one distribution from another
  has the form \f$ D(p,q) = \sum_i \ln(p_i/q_i) p_i \f$.
  Notice that it is not a distance, because there is an asymmetry
  between \f$p\f$ and \f$q\f$, so one can expect that \f$D(p, q) \neq D(q, p)\f$.

  \param top the \f$p\f$ in the above formula. (No default; must not be \c NULL)
  \param bottom the \f$q\f$ in the above formula. (No default; must not be \c NULL)
  \param draw_ct If I do the calculation via random draws, how many? (Default = 1e5)
  \param rng    A \c gsl_rng. If NULL, I'll take care of the RNG; see \ref autorng. (Default = \c NULL)

  This function can take empirical histogram-type models---\ref apop_pmf and \ref apop_histogram ---or continuous models like \ref apop_loess
  or \ref apop_normal.

 If there is an empirical model (I'll try \c top first, under the presumption that you are measuring the divergence of data from a `true' distribution), then I'll step
through it for the points in the summation.

\li If you have two empirical distributions, that they must be synced: if \f$p_i>0\f$
but \f$q_i=0\f$, then the function returns \c GSL_NEGINF. If <tt>apop_opts.verbose >=1</tt>
I print a message as well.

If neither distribution is empirical, then I'll take \c draw_ct random draws from \c bottom and evaluate at those points.

\li Set <tt>apop_opts.verbose = 3</tt> for observation-by-observation info.

This function uses the \ref designated syntax for inputs.
 */
APOP_VAR_HEAD double apop_kl_divergence(apop_model *top, apop_model *bottom, int draw_ct, gsl_rng *rng){
    apop_model * apop_varad_var(top, NULL);
    apop_model * apop_varad_var(bottom, NULL);
    Apop_assert(top, "The first model is NULL.");
    Apop_assert(bottom, "The second model is NULL.");
    double apop_varad_var(draw_ct, 1e5);
    static gsl_rng * spare_rng = NULL;
    gsl_rng * apop_varad_var(rng, NULL);
    if (!rng && !spare_rng) 
        spare_rng = apop_rng_alloc(++apop_opts.rng_seed);
    if (!rng)  rng = spare_rng;
APOP_VAR_ENDHEAD
    double div = 0;
    Apop_notify(3, "p(top)\tp(bot)\ttop*log(top/bot)\n");
    if (top->name && !strcmp(top->name, "PDF or sparse matrix")){
        apop_data *p = top->data;
        Get_vmsizes(p); //firstcol, vsize, msize1, msize2
        apop_data *a_row = apop_data_alloc(vsize, (msize1 ? 1 : 0), msize2);
        for (int i=0; i < (vsize ? vsize : msize1); i++){
            double pi = p->weights ? gsl_vector_get(p->weights, i) : 1./(vsize ? vsize : msize1);
            if (!pi){
                Apop_notify(3, "0\t--\t0");
                continue;
            } //else:
            get_one_row(p, a_row, i, firstcol, msize2);
            double qi = apop_p(a_row, bottom);
            Apop_assert_c(qi, GSL_NEGINF, 1, "The PMFs aren't synced: bottom has a value where "
                                                "top doesn't (which produces infinite divergence).");
            Apop_notify(3,"%g\t%g\t%g", pi, qi, pi ? pi * log(pi/qi):0);
            div += pi * log(pi/qi);
        }
        apop_data_free(a_row);
    } else { //the version with the RNG.
        apop_data *a_row = apop_data_alloc(1, top->dsize);
        for (int i=0; i < draw_ct; i++){
            apop_draw(a_row->matrix->data, rng, top);
            double pi = apop_p(a_row, top);
            double qi = apop_p(a_row, bottom);
            Apop_assert_c(qi, GSL_NEGINF, 1, "The PMFs aren't synced: bottom has a value where "
                                                "top doesn't (which produces infinite divergence).");
            Apop_notify(3,"%g\t%g\t%g", pi, qi, pi ? pi * log(pi/qi):0);
            if (pi) //else add zero.
                div += pi * log(pi/qi);
        }
        apop_data_free(a_row);
    }
    return div;
}


/** \defgroup tfchi t-, chi-squared, F-, Wishart distributions

Most of these distributions are typically used for testing purposes.  For such a situation, you don't need the models here.
Given a statistic of the right properties, you can find the odds that the statistic is above or below a cutoff on the t-, F, or chi-squared distribution using the \ref apop_test function. 

In that world, those three distributions are actually parameter free. The data is assumed to be normalized to be based on a mean zero, variance one process, you get the degrees of freedom from the size of the data, and the distribution is fixed.

For modeling purposes, more could be done. For example, the t-distribution is a favorite proxy for Normal-like situations where there are fat tails relative to the Normal (i.e., high kurtosis). Or, you may just prefer not to take the step of normalizing your data---one could easily rewrite the theorems underlying the t-distribution without the normalizations.

In such a case, the researcher would not want to fix the \f$df\f$, because \f$df\f$ indicates the fatness of the tails, which has some optimal value given the data. 
Thus, there are two modes of use for these distributions: 

\li Parameterized, testing style: the degrees of freedom are determined
from the data, and all necessary normalizations are assumed. Thus, this code---

\code
apop_data *t_for_testing = apop_estimate(data, apop_t)
\endcode

---will return exactly the type of \f$t\f$-distribution one would use for testing. 

\li By removing the \c estimate method---
\code
apop_model *spare_t = apop_model_copy(t);
spare_t.estimate = NULL;
apop_model *best_fitting_t = apop_estimate(your_data, spare_t);
\endcode
---I will find the best \f$df\f$ via maximum likelihood, which may be desirable for
to find the best-fitting model for descriptive purposes.

\c df works for all four distributions here; \c df2 makes sense only for the \f$F\f$, 

For the Wishart, the degrees of freedom and covariance matrix are always estimated via MLE.
*/

/** The multivariate generalization of the Gamma distribution.
\f[
\Gamma_p(a)=
\pi^{p(p-1)/4}\prod_{j=1}^p
\Gamma\left[ a+(1-j)/2\right]. \f]

Because \f$\Gamma(x)\f$ is undefined for \f$x\in\{0, -1, -2, ...\}\f$, this function returns \c GSL_NAN when \f$a+(1-j)/2\f$ takes on one of those values.

See also \ref apop_multivariate_lngamma, which is more numerically stable in most cases.
*/
double apop_multivariate_gamma(double a, double p){
    Apop_assert_c(!(-(a+(1-p)/2) == (int)-(a+(1-p)/2) && a+(1-p)/2 <=0), GSL_NAN, 1, "Undefined when a + (1-j)/2 = 0, -1, -2, ... [you sent a=%g, p=%g]", a, p);
    double out = pow(M_PI, p*(p-1.)/4.);
    double factor = 1;
    for (int i=1; i<=p; i++)
        factor *= gsl_sf_gamma(a+(1-i)/2.);
    return out * factor;
}

/** The log of the multivariate generalization of the Gamma; see also
 \ref apop_multivariate_gamma.
*/
double apop_multivariate_lngamma(double a, double p){
    Apop_assert_c(!(-(a+(1-p)/2) == (int)-(a+(1-p)/2) && a+(1-p)/2 <=0), GSL_NAN, 1, "Undefined when a + (1-j)/2 = 0, -1, -2, ... [you sent a=%g, p=%g]", a, p);
    double out = M_LNPI * p*(p-1.)/4.;
    double factor = 0;
    for (int i=1; i<=p; i++)
        factor += gsl_sf_lngamma(a+(1-i)/2.);
    return out + factor;
}

static void find_eigens(gsl_matrix **subject, gsl_vector *eigenvals, gsl_matrix *eigenvecs){
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc((*subject)->size1);
    gsl_eigen_symmv(*subject, eigenvals, eigenvecs, w);
    gsl_eigen_symmv_free (w);
    gsl_matrix_free(*subject); *subject  = NULL;
}

static void diagonal_copy(gsl_vector *v, gsl_matrix *m, char in_or_out){
    gsl_vector_view dv = gsl_matrix_diagonal(m);
    if (in_or_out == 'i') gsl_vector_memcpy(&(dv.vector), v);
    else                  gsl_vector_memcpy(v, &(dv.vector));
}

static double diagonal_size(gsl_matrix *m){
    gsl_vector_view dv = gsl_matrix_diagonal(m);
    return apop_sum(&dv.vector);
}

static double biggest_elmt(gsl_matrix *d){ 
    return  GSL_MAX(fabs(gsl_matrix_max(d)), fabs(gsl_matrix_min(d)));
}

/** Test whether the input matrix is positive semidefinite.

A covariance matrix will always be PSD, so this function can tell you whether your matrix is a valid covariance matrix.

Consider the 1x1 matrix in the upper left of the input, then the 2x2 matrix in the upper left, on up to the full matrix. If the matrix is PSD, then each of these has a positive determinant. This function thus calculates \f$N\f$ determinants for an \f$N\f$x\f$N\f$ matrix.

\param m The matrix to test. If \c NULL, I will return zero---not PSD.
\param semi If anything but 's', check for positive definite, not semidefinite. (default 's')

See also \ref apop_matrix_to_positive_semidefinite, which will change the input to something PSD.

This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD int apop_matrix_is_positive_semidefinite(gsl_matrix *m, char semi){
    gsl_matrix * apop_varad_var(m, NULL);
    Apop_assert_c(m, 0, 1, "You gave me a NULL matrix. I will take this as not positive semidefinite.");
    char apop_varad_var(semi, 's');
APOP_VAR_ENDHEAD
    for (int i=1; i<= m->size1; i++){
        gsl_matrix mv =gsl_matrix_submatrix (m, 0, 0, i, i).matrix;
        double det = apop_matrix_determinant(&mv);
        if ((semi == 'd' && det <0) || det <=0)
            return 0;
    }
    return 1;
}

void vfabs(double *x){*x = fabs(*x);}

/**  First, this function passes tests, but is under development.
  
It takes in a matrix and converts it to the `closest' positive semidefinite matrix.

\param m On input, any matrix; on output, a positive semidefinite matrix.
\return the distance between the original and new matrices.

\li See also the test function \ref apop_matrix_is_positive_semidefinite.
\li This function can be used as (the core of) a model constraint.

Adapted from the R Matrix package's nearPD, which is 
Copyright (2007) Jens Oehlschlägel [and is GPL].
*/
double apop_matrix_to_positive_semidefinite(gsl_matrix *m){
    if (apop_matrix_is_positive_semidefinite(m)) return 0; 
    double diffsize=0, dsize;
    apop_data *qdq; 
    gsl_matrix *d = apop_matrix_copy(m);
    gsl_matrix *original = apop_matrix_copy(m);
    double orig_diag_size = fabs(diagonal_size(d));
    int size = d->size1;
    gsl_vector *diag = gsl_vector_alloc(size);
    diagonal_copy(diag, d, 'o');
    apop_vector_apply(diag, vfabs);
    double origsize = biggest_elmt(d);
    do {
        //get eigenvals
        apop_data *eigenvecs = apop_data_alloc(size, size);
        gsl_vector *eigenvals = gsl_vector_calloc(size);
        gsl_matrix *junk_copy = apop_matrix_copy(d);
        find_eigens(&junk_copy, eigenvals, eigenvecs->matrix);//junk freed here.
        
        //prune positive only
        int j=0;
        int plussize = eigenvecs->matrix->size1;
        int *mask = calloc(eigenvals->size , sizeof(int));
        for (int i=0; i< eigenvals->size; i++)
            plussize -= 
            mask[i] = (gsl_vector_get(eigenvals, i) <= 0);
        
        //construct Q = pruned eigenvals
        apop_data_rm_columns(eigenvecs, mask);
        if (!eigenvecs->matrix) break;
        
        //construct D = positive eigen diagonal
        apop_data *eigendiag = apop_data_calloc(0, plussize, plussize);
        for (int i=0; i< eigenvals->size; i++)
            if (!mask[i]) {
                apop_data_set(eigendiag, j, j, eigenvals->data[i]);
                j++;
            }

        // Our candidate is QDQ', symmetrized, with the old diagonal subbed in.
        apop_data *qd = apop_dot(eigenvecs, eigendiag);
        qdq = apop_dot(qd, eigenvecs, .form2='t');
        for (int i=0; i< qdq->matrix->size1; i++)
            for (int j=i+1; j< qdq->matrix->size1; j++){
                double avg = (apop_data_get(qdq, i, j) +apop_data_get(qdq, j, i)) /2.;
                apop_data_set(qdq, i, j, avg);
                apop_data_set(qdq, j, i, avg);
            }
        diagonal_copy(diag, qdq->matrix, 'i');
        
        // Evaluate progress, clean up.
        dsize = biggest_elmt(d);
        gsl_matrix_sub(d, qdq->matrix);
        diffsize = biggest_elmt(d);
        apop_data_free(qd); gsl_matrix_free(d);
        apop_data_free(eigendiag); free(mask);
        apop_data_free(eigenvecs); gsl_vector_free(eigenvals);
        d = qdq->matrix;
        qdq->matrix=NULL; apop_data_free(qdq); qdq = NULL;
    } while (diffsize/dsize > 1e-3);

    apop_data *eigenvecs = apop_data_alloc(size, size);
    gsl_vector *eigenvals = gsl_vector_calloc(size);
    gsl_matrix *junk_copy = apop_matrix_copy(d);
    find_eigens(&junk_copy, eigenvals, eigenvecs->matrix);
    //make eigenvalues more positive
    double score =0;
    for (int i=0; i< eigenvals->size; i++){
        double v = gsl_vector_get(eigenvals, i);
        if (v < 1e-1){
            gsl_vector_set(eigenvals, i, 1e-1);
            score += 1e-1 - v;
        }
    }
    for (int i=0; i< size; i++)
        assert(eigenvals->data[i] >=0);
    //if (score){
        apop_data *eigendiag = apop_data_calloc(0, size, size);
        diagonal_copy(eigenvals, eigendiag->matrix, 'i');
        double new_diag_size = diagonal_size(eigendiag->matrix);
        gsl_matrix_scale(eigendiag->matrix, orig_diag_size/new_diag_size);
        apop_data *qd = apop_dot(eigenvecs, eigendiag);
        qdq = apop_dot(qd, eigenvecs, .form2='t');
        gsl_matrix_memcpy(m, qdq->matrix);
        apop_data_free(qd);
        apop_data_free(eigendiag);
    //}
    assert(apop_matrix_is_positive_semidefinite(m));
    apop_data_free(qdq); gsl_vector_free(diag);
    apop_data_free(eigenvecs); gsl_vector_free(eigenvals);
    gsl_matrix_sub(original, m);
    return biggest_elmt(original)/origsize;
}
