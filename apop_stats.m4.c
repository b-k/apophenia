/** \file apop_stats.c	Basic moments and some distributions. */
/* Copyright (c) 2006--2007, 2013 by Ben Klemens.  Licensed under the GPLv2; see COPYING.  */

#include "apop_internal.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_eigen.h>

#define Check_vw    \
    Apop_stopif(!v, return GSL_NAN, 0, "data vector is NULL. Returning NaN.\n");            \
    Apop_stopif(!v->size, return GSL_NAN, 0, "data vector has size 0. Returning NaN.\n");   \
    Apop_stopif(weights && weights->size != v->size, return GSL_NAN, 0, "data vector has size %zu; weighting vector has size %zu. Returning NaN.\n", v->size, weights->size);

/** Returns the sum of the data in the given vector.
*/
long double apop_vector_sum(const gsl_vector *in){
    Apop_stopif(!in, return 0, 1, "You just asked me to sum a NULL. Returning zero.");
    long double out = 0;
    for (size_t i=0; i< in->size; i++)
        out += gsl_vector_get(in, i);
	return out; 
}

/** \def apop_sum(in)
  An alias for \ref apop_vector_sum. Returns the sum of the data in the given vector.
*/

/** \def apop_mean(v)
 Returns the mean of the elements of the vector \c v.

\param v A \ref gsl_vector.
*/



/** \def apop_var(in)
An alias for \ref apop_vector_var.
Returns the variance of the data in the given vector.
*/

/** Returns an unbiased estimate of the sample skew of the data in the given vector.
*/
double apop_vector_skew(const gsl_vector *in){
	return apop_vector_skew_pop(in) * gsl_pow_2(in->size)/((in->size -1.)*(in->size -2.)); }

/** Returns the sample kurtosis of the data in the given
vector. Corrections are made to produce an unbiased result as per <a
href="http://modelingwithdata.org/pdfs/moments.pdf">Appendix M</a> (PDF) of <em>Modeling
with data</em>.

\li This is an estimate of the fourth central moment without normalization. The kurtosis
    of a \f${\cal N}(0,1)\f$ is \f$3 \sigma^4\f$, not three, one, or zero.
\see \ref apop_vector_kurtosis_pop
*/
double apop_vector_kurtosis(const gsl_vector *in){
    size_t n = in->size;
    long double coeff0= n*n/(gsl_pow_3(n-1)*(gsl_pow_2(n)-3*n+3));
    long double coeff1= n*gsl_pow_2(n-1)+ (6*n-9);
    long double coeff2= n*(6*n-9);
    return  coeff0 *(coeff1 * apop_vector_kurtosis_pop(in) - coeff2 * gsl_pow_2(apop_vector_var(in)*(n-1.)/n));
}

static double wskewkurt(const gsl_vector *v, const gsl_vector *w, const int exponent, const char *fn_name){
    long double wsum = 0, sumcu = 0, vv, ww, mu;
    //Using the E(x - \bar x)^3 form, which is lazy.
    mu  = apop_vector_mean(v, w);
    for (size_t i=0; i< w->size; i++){
        vv    = gsl_vector_get(v,i);
        ww    = gsl_vector_get(w,i);
        sumcu+= ww * gsl_pow_int(vv - mu, exponent); 
        wsum += ww; 
    }
    double len = wsum < 1.1 ? w->size : wsum;
    return sumcu/len;
}

/** Returns the population skew \f$(\sum_i (x_i - \mu)^3/n))\f$ of the data in the given vector. Observations may be weighted.

\param v       The data vector
\param weights The weight vector. Default: equal weights for all observations.
\return        The weighted skew.
 
\li  Some people like to normalize the skew by dividing by variance\f$^{3/2}\f$; that's not done here, so you'll have to do so separately if need be.

\li Apophenia tries to be smart about reading the weights. If weights
sum to one, then the system uses \c w->size as the number of elements,
and returns the usual sum over \f$n-1\f$. If weights > 1, then the
system uses the total weights as \f$n\f$. Thus, you can use the weights
as standard weightings or to represent elements that appear repeatedly.
*/
APOP_VAR_HEAD double apop_vector_skew_pop(gsl_vector const *v, gsl_vector const *weights){
    gsl_vector const * apop_varad_var(v, NULL);
    gsl_vector const * apop_varad_var(weights, NULL);
    Check_vw
APOP_VAR_ENDHEAD
    if (weights) return wskewkurt(v, weights, 3, "apop_vector_weighted_skew");

    //This code is cut/pasted/modified from the GSL. 
    //I reimplement the skew calculation here without the division by var^3/2 that the GSL does. 
    size_t n = v->size;
    long double avg = 0;
    long double mean = apop_vector_mean(v);
    for (size_t i = 0; i < n; i++) {
        const long double x = gsl_vector_get(v, i) - mean; 
        avg += (gsl_pow_3(x) - avg)/(i + 1);
    } 
    return avg;
}

/** Returns the population kurtosis [\f$\sum_i (x_i - \mu)^4/n)\f$] of the data in
the given vector, with an optional weighting.

\param v The data vector
\param weights The weight vector. If NULL, assume equal weights.
\return The weighted kurtosis.
 
  \li Some people like to normalize the fourth central moment by dividing by variance
squared, or by subtracting three; those things are not done here, so you'll have to
do them separately if need be.
  \li This function uses the \ref designated syntax for inputs.
\see \ref apop_vector_kurtosis for the unbiased sample version.
*/
APOP_VAR_HEAD double apop_vector_kurtosis_pop(gsl_vector const *v, gsl_vector const *weights){
    gsl_vector const * apop_varad_var(v, NULL);
    gsl_vector const * apop_varad_var(weights, NULL);
    Check_vw
APOP_VAR_ENDHEAD
    if (weights) return wskewkurt(v, weights, 4, "apop_vector_weighted_kurtosis");

    //This code is cut/pasted/modified from the GSL. 
    //I reimplement the kurtosis calculation here without the division by var^2 that the GSL does. 
    size_t n = v->size;
    long double avg  = 0;
    long double mean = apop_vector_mean(v);
    for (size_t i = 0; i < n; i++) {
        const long double x = gsl_vector_get(v, i) - mean; 
        avg += (gsl_pow_4(x) - avg)/(i + 1);
    } 
    return avg;
}

/** Returns the variance of the data in the given vector, given that you've already calculated the mean.
\param in	the vector in question
\param mean	the mean, which you've already calculated using \ref apop_vector_mean.
\see apop_vector_var
*/
double apop_vector_var_m(const gsl_vector *in, const double mean){
	return gsl_stats_variance_m(in->data,in->stride, in->size, mean); }

/** Returns the correlation coefficient of two vectors:
\f$ {\hbox{cov}(a,b)\over \sqrt(\hbox{var}(a)) \sqrt(\hbox{var}(b))}.\f$

An example
\code 
gsl_matrix *m = [fill matrix here];
printf("The correlation coefficient between rows two "
       "and three is %g.\n", apop_vector_correlation(Apop_mrv(m, 2), Apop_mrv(m, 3)));
\endcode 

  \li \c ina and \c inb are mandatory and must be non-NULL, but you can leave the weights
vector off in the typical case of uniform weights.
  \li This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD double apop_vector_correlation(const gsl_vector *ina, const gsl_vector *inb, const gsl_vector *weights){
    gsl_vector const * apop_varad_var(ina, NULL);
    gsl_vector const * apop_varad_var(inb, NULL);
    gsl_vector const * apop_varad_var(weights, NULL);
APOP_VAR_ENDHEAD
	return apop_vector_cov(ina, inb, weights) 
            / sqrt(apop_vector_var(ina, weights) * apop_vector_var(inb, weights)); }


/** Returns the distance between two vectors, where distance is defined
 based on the third (optional) parameter:

 - 'e' or 'E' (the default): scalar distance (standard Euclidean metric) between two vectors. \f$\sqrt{\sum_i{(a_i - b_i)^2}},\f$
where \f$i\f$ iterates over dimensions.
 - 'm' or 'M'  Returns the Manhattan metric distance  between two vectors: \f$\sum_i{|a_i - b_i|},\f$
where \f$i\f$ iterates over dimensions.
 - 'd' or 'D' The discrete norm: if \f$a = b\f$, return zero, else return one.
 - 's' or 'S' The sup norm: find the dimension where \f$|a_i - b_i|\f$ is largest, return the distance along that one dimension.
 - 'l' or 'L' The \f$L_p\f$ norm, \f$\left(\sum_i{|a_i - b_i|^2}\right)^{1/p}\f$. The value of \f$p\f$ is set by the fourth (optional) argument.

 \param ina First vector (No default, must not be \c NULL)
 \param inb Second vector (Default = zero)
 \param metric The type of metric, as above.
 \param norm  If you are using an \f$L_p\f$ norm, this is \f$p\f$. Must be strictly greater than zero. (default = 2)

\li  The defaults are such that
 \code
 apop_vector_distance(v);
 apop_vector_distance(v, .metric = 's');
 \endcode
 gives you the standard Euclidean length of \c v and its longest element.

\include test_distances.c

\li This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD double apop_vector_distance(const gsl_vector *ina, const gsl_vector *inb, const char metric, const double norm){
    static threadlocal gsl_vector *zero = NULL;
    const gsl_vector * apop_varad_var(ina, NULL);
    Apop_stopif(!ina, return NAN, 1, "The first vector is NULL. Returning NAN");
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
    Apop_stopif(ina->size != inb->size, return GSL_NAN, 0, "I need equal-sized vectors, but "
                "you sent a vector of size %zu and a vector of size %zu. Returning NaN.", ina->size, inb->size);
    double dist = 0;
    if (metric == 'e' || metric == 'E'){
        for (size_t i=0; i< ina->size; i++)
            dist += gsl_pow_2(gsl_vector_get(ina, i) - gsl_vector_get(inb, i));
        return sqrt(dist); 
    }
    if (metric == 'm' || metric == 'M'){ //redundant with vector_grid_distance, below.
        for (size_t i=0; i< ina->size; i++) 
            dist += fabs(gsl_vector_get(ina, i) - gsl_vector_get(inb, i));
        return dist; 
    }
    if (metric == 'd' || metric == 'D'){
        for (size_t i=0; i< ina->size; i++) 
            if (gsl_vector_get(ina, i) != gsl_vector_get(inb, i))
                return 1;
        return 0;
    }
    if (metric == 's' || metric == 'S'){
        for (size_t i=0; i< ina->size; i++) 
            dist = GSL_MAX(dist, fabs(gsl_vector_get(ina, i) - gsl_vector_get(inb, i)));
        return dist;
    }
    if (metric == 'l' || metric == 'L'){
        for (size_t i=0; i< ina->size; i++)
            dist += pow(fabs(gsl_vector_get(ina, i) - gsl_vector_get(inb, i)), norm);
        return pow(dist, 1./norm); 
    }
  Apop_stopif(1, return NAN, 1, "I couldn't find the metric type you gave, %c, in my list of supported types. Returning NaN", metric);
}

/** This function will normalize a vector, either such that it has mean
zero and variance one, or ranges between zero and one, or sums to one.

\param in	A \c gsl_vector with the un-normalized data. \c NULL
input gives \c NULL output. (No default)

\param out 	If normalizing in place, \c NULL.
If not, the address of a <tt>gsl_vector*</tt>. Do not allocate. (default = \c NULL.)

\param normalization_type 
\c 'p': normalized vector will sum to one. E.g., start with a set of observations in bins, end with the percentage of observations in each bin. (the default)<br>
\c 'r': normalized vector will range between zero and one. Replace each X with (X-min) / (max - min).<br>
\c 's': normalized vector will have mean zero and (sample) variance one. Replace
each X with \f$(X-\mu) / \sigma\f$, where \f$\sigma\f$ is the sample
standard deviation.<br>
\c 'm': normalize to mean zero: Replace each X with \f$(X-\mu)\f$<br>

\b Example 
\code
\endcode

\li This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD void apop_vector_normalize(gsl_vector *in, gsl_vector **out, const char normalization_type){
    gsl_vector * apop_varad_var(in, NULL);
    Apop_stopif(!in, return, 1, "Input vector is NULL. Doing nothing.");
    gsl_vector ** apop_varad_var(out, NULL);
    const char apop_varad_var(normalization_type, 'p');
APOP_VAR_END_HEAD
    double mu, min, max;
	if (!out) out = &in;
	else {
		*out = gsl_vector_alloc (in->size);
		gsl_vector_memcpy(*out, in);
	}
	if (normalization_type == 's'){
		mu = apop_vector_mean(in);
        Apop_stopif(!isfinite(mu), return, 0, "normalization failed: the mean of the vector is not finite.");
		gsl_vector_add_constant(*out, -mu);
        double scaling = 1./(sqrt(apop_vector_var_m(*out, 0)));
        Apop_stopif(!isfinite(scaling), return, 0, "normalization failed: 1/(std error)  of the vector is not finite.");
		gsl_vector_scale(*out, scaling);
	} 
	else if (normalization_type == 'r'){
        gsl_vector_minmax(in, &min, &max);
		gsl_vector_add_constant(*out, -min);
		gsl_vector_scale(*out, 1/(max-min));	

	}
	else if (normalization_type == 'p'){
		long double sum	= apop_sum(in);
        Apop_stopif(!sum, return, 0, "the vector sums to zero, so I can't normalize it to sum to one.");
		gsl_vector_scale(*out, 1/sum);	
	}
	else if (normalization_type == 'm'){
		mu = apop_vector_mean(in);
        Apop_stopif(!isfinite(mu), return, 0, "normalization failed: the mean of the vector is not finite.");
		gsl_vector_add_constant(*out, -mu);
	}
}

/** Returns the sum of the elements of a matrix. Occasionally convenient.
  \param m	the matrix to be summed. 
*/
long double apop_matrix_sum(const gsl_matrix *m){
    Apop_stopif(!m, return 0, 1, "You just asked me to sum a NULL. Returning zero.");
    long double	sum	= 0;
	for (size_t j=0; j< m->size1; j++)
		for (size_t i=0; i< m->size2; i++)
			sum += gsl_matrix_get(m, j, i);
	return sum;
}

/** Returns the mean of all elements of a matrix.

\param data	The matrix to be averaged. If \c NULL, return zero.
\return The mean of all cells of the matrix.
*/
double apop_matrix_mean(const gsl_matrix *data){
    if (!data) return 0;
    long double avg = 0;
    int cnt = 0;
    for(size_t i=0; i < data->size1; i++)
        for(size_t j=0; j < data->size2; j++){
            double x = gsl_matrix_get(data, i,j);
            long double ratio = cnt/(cnt+1.0);
            cnt++;
            avg*= ratio;
            avg+= x/cnt;
        }
	return avg;
}

/** Returns the mean and population variance of all elements of a matrix.
 
\li If \c NULL, return \f$\mu=0, \sigma^2=NaN\f$.
\li Gives the population variance (sum of squares divided by \f$N\f$).  
If you want sample variance, multiply the result by \f$N/(N-1)\f$:
\code
double mu, var;
apop_data *data= apop_query_to_data("select * from indata");
apop_matrix_mean_and_var(data->matrix, &mu, &var);
var *= (data->size1*data->size2)/(data->size1*data->size2-1.0);
\endcode

\param data	the matrix to be averaged. 
\param	mean	where to put the mean to be calculated.
\param	var	where to put the variance to be calculated.
*/
void apop_matrix_mean_and_var(const gsl_matrix *data, double *mean, double *var){
    if (!data) {*mean=0; *var=GSL_NAN; return;}
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

\param indata The table to be summarized. An \ref apop_data structure. May have a <tt>weights</tt> element.
\return     An \ref apop_data structure with one row for each column in the original
            table, and a column for each summary statistic.
\exception out->error='a'  Allocation error.

\li This function gives more columns than you probably want; use \ref apop_data_prune_columns to pick the ones you want to see.

\li See apop_data_prune_columns for an example.
*/
apop_data * apop_data_summarize(apop_data *indata){
    Apop_stopif(!indata, return NULL, 0, "You sent me a NULL apop_data set. Returning NULL.");
    Apop_stopif(!indata->matrix, return NULL, 0, "You sent me an apop_data set with a NULL matrix. Returning NULL.");
    apop_data *out = apop_data_alloc(indata->matrix->size2, 6);
    double mean, var;
    char rowname[10000]; //crashes on more than 10^9995 columns.
	apop_name_add(out->names, "mean", 'c');
	apop_name_add(out->names, "std dev", 'c');
	apop_name_add(out->names, "variance", 'c');
	apop_name_add(out->names, "min", 'c');
	apop_name_add(out->names, "median", 'c');
	apop_name_add(out->names, "max", 'c');
	if (indata->names !=NULL){
        apop_name_stack(out->names,indata->names, 'r', 'c');
        if (indata->names->title && strlen(indata->names->title)){
            char *title;
            Asprintf(&title, "summary for %s", indata->names->title);
            apop_name_add(out->names, title, 'h');
            free(title);
        }
    }
	else
		for (size_t i=0; i< indata->matrix->size2; i++){
			sprintf(rowname, "col %zu", i);
			apop_name_add(out->names, rowname, 'r');
		}
	for (size_t i=0; i< indata->matrix->size2; i++){
        gsl_vector *v = Apop_cv(indata, i);
        if (!indata->weights){
            mean = apop_vector_mean(v);
            var  = apop_vector_var_m(v, mean);
        } else {
            mean = apop_vector_mean(v, indata->weights);
            var  = apop_vector_var(v, indata->weights);
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

/** Find the mean, weighted or unweighted. 

\param v        The data vector
\param weights  The weight vector. Default: assume equal weights.
\return         The weighted mean 
\li This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD double apop_vector_mean(gsl_vector const *v, gsl_vector const *weights){
    gsl_vector const * apop_varad_var(v, NULL);
    gsl_vector const * apop_varad_var(weights, NULL);
    Check_vw
APOP_VAR_END_HEAD
    if (!weights) return gsl_stats_mean(v->data, v->stride, v->size);
    long double sum = 0, wsum = 0;
    for (size_t i=0; i< weights->size; i++){
        sum  += gsl_vector_get(weights, i) * gsl_vector_get(v,i); 
        wsum += gsl_vector_get(weights, i); 
    }
    return sum/wsum;
}

/** Find the sample variance of a vector, weighted or unweighted.

\param v       The data vector
\param weights The weight vector. If NULL (the default), assume equal weights.
\return        The weighted sample variance.  

  \li This uses (n-1) in the denominator of the sum; i.e., it corrects for the bias
introduced by using \f$\bar x\f$ instead of \f$\mu\f$.
  \li  At the moment, there is no var_pop function. Just multiply this by (n-1)/n if you need that.
  \li Apophenia tries to be smart about reading the weights. If weights
sum to one, then the system uses \c w->size as the number of elements,
and returns the usual sum over \f$n-1\f$. If weights > 1, then the
system uses the total weights as \f$n\f$. Thus, you can use the weights
as standard weightings or to represent elements that appear repeatedly.
  \li This function uses the \ref designated syntax for inputs.
\see apop_vector_var_m for the case where you already have the vector's mean.
*/
APOP_VAR_HEAD double apop_vector_var(gsl_vector const *v, gsl_vector const *weights){
    gsl_vector const * apop_varad_var(v, NULL);
    gsl_vector const * apop_varad_var(weights, NULL);
    Check_vw
APOP_VAR_END_HEAD
    if (!weights) return gsl_stats_variance(v->data, v->stride, v->size);
    //Using the E(x^2) - E^2(x) form.
    long double sum = 0, wsum = 0, sumsq = 0, vv, ww;
    for (size_t i=0; i< weights->size; i++){
        vv = gsl_vector_get(v, i);
        ww = gsl_vector_get(weights, i);
        sum   += ww * vv;
        sumsq += ww * gsl_pow_2(vv); 
        wsum  += ww; 
    }
    double len = (wsum < 1.1 ? weights->size : wsum);
    return (sumsq/len - gsl_pow_2(sum/len)) * len/(len -1.);
}

/** Find the sample covariance of a pair of vectors, with an optional weighting. This only
makes sense if the weightings are identical, so the function takes only one weighting vector for both.

\param  v1, v2  The data vectors (no default; must not be \c NULL)
\param  weights The weight vector. (default equal weights for all elements)
\return The sample covariance

\li This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD double apop_vector_cov(const gsl_vector *v1, const gsl_vector *v2, const gsl_vector *weights){
    gsl_vector const * apop_varad_var(v1, NULL);
    gsl_vector const * apop_varad_var(v2, NULL);
    gsl_vector const * apop_varad_var(weights, NULL);
    Apop_stopif(!v1, return GSL_NAN, 0, "first data vector is NULL. Returning NaN.");
    Apop_stopif(!v2, return GSL_NAN, 0, "second data vector is NULL. Returning NaN.");
    Apop_stopif(!v1->size, return GSL_NAN, 0, "first data vector has size 0. Returning NaN.");
    Apop_stopif(!v2->size, return GSL_NAN, 0, "second data vector has size 0. Returning NaN.");
    Apop_stopif(v1->size!= v2->size, return GSL_NAN, 0, "data vectors have sizes %zu and %zu. Returning NaN.", v1->size, v2->size);
    Apop_stopif(weights && ((weights->size != v1->size) || (weights->size != v2->size)), return GSL_NAN, 0, "data vectors have sizes %zu and %zu; weighting vector has size %zu. Returning NaN.", v1->size, v2->size, weights->size);

APOP_VAR_ENDHEAD
    if (!weights) return gsl_stats_covariance(v1->data, v1->stride, v2->data, v2->stride, v2->size);
    long double sum1 = 0, sum2 = 0, wsum = 0, sumsq = 0, vv1, vv2, ww;
    //Using the E(x^2) - E^2(x) form.
    for (size_t i=0; i< weights->size; i++){
        vv1   = gsl_vector_get(v1,i);
        vv2   = gsl_vector_get(v2,i);
        ww    = gsl_vector_get(weights,i);
        sum1 += ww * vv1;
        sum2 += ww * vv2;
        sumsq+= ww * vv1 * vv2;
        wsum += ww; 
    }
    double len = (wsum < 1.1 ? weights->size : wsum);
    return (sumsq/len  - sum1*sum2/gsl_pow_2(len)) *(len/(len-1));
}

/** Returns the sample variance/covariance matrix relating each column of the matrix to each other column.

\param in An \ref apop_data set. If the weights vector is set, I'll take it into account.

\li This is the sample covariance---dividing by \f$n-1\f$, not \f$n\f$. If you need the population variance, use 
\code
apop_data *popcov = apop_data_covariance(indata);
int size=indata->matrix->size1;
gsl_matrix_scale(popcov->matrix, size/(size-1.));
\endcode

\return Returns an \ref apop_data set the variance/covariance matrix.  
\exception out->error='a'  Allocation error.
*/
apop_data *apop_data_covariance(const apop_data *in){
    Apop_stopif(!in, return NULL, 1, "You sent me a NULL apop_data set. Returning NULL.");
    Apop_stopif(!in->matrix, return NULL, 1, "You sent me an apop_data set with a NULL matrix. Returning NULL.");
    apop_data *out = apop_data_alloc(in->matrix->size2, in->matrix->size2);
    Apop_stopif(out->error, return out, 0, "allocation error.");
    for (size_t i=0; i < in->matrix->size2; i++){
        for (size_t j=i; j < in->matrix->size2; j++){
            double var = apop_vector_cov(Apop_cv(in, i), Apop_cv(in, j), in->weights);
            gsl_matrix_set(out->matrix, i,j, var);
            if (i!=j) gsl_matrix_set(out->matrix, j,i, var);
        }
    }
    apop_name_stack(out->names, in->names, 'c');
    apop_name_stack(out->names, in->names, 'r', 'c');
    return out;
}

/** Returns the matrix of correlation coefficients \f$(\sigma^2_{xy}/(\sigma_x\sigma_y))\f$ relating each column with each other.

\param in 	A data matrix: rows are observations, columns are variables. If you give me a weights vector, I'll use it.

\return Returns the square variance/covariance matrix with dimensions equal to the number of input columns.
\exception out->error='a'  Allocation error.
*/
apop_data *apop_data_correlation(const apop_data *in){
    apop_data *out = apop_data_covariance(in);
    if (!out) return NULL;
    for(size_t i=0; i< in->matrix->size2; i++){
        double std_dev = sqrt(apop_vector_var(Apop_cv(in, i), in->weights));
        gsl_vector_scale(Apop_cv(out, i), 1.0/std_dev);
        gsl_vector_scale(Apop_rv(out, i), 1.0/std_dev);
    }
    return out;
}


/** Given a vector representing a probability distribution of observations, calculate the entropy, \f$\sum_i -\ln(v_i)v_i\f$.

\li The input vector need not be normalized to sum to one. You may input a vector giving frequencies.

\li The entropy of a data set depends only on the frequency with which elements are
observed, not the value of the elements themselves. The \ref apop_data_pmf_compress
function will reduce an input \ref apop_data set to one weighted line per observation, and
the weights would determine the entropy:

\code
apop_data *data = apop_text_to_data("indata");
apop_data_pmf_compress(data);
data_entropy = apop_vector_entropy(d->weights);
\endcode

\li The entropy is calculated using natural logs. To convert to base 2, divide by \f$\ln(2)\f$; see the example.

\li The entropy of an empty data set (\c NULL or a total weight of zero) is zero. Print a warning when given \c NULL
    input and <tt>apop_opts.verbose >=1</tt>.

\li If the input vector has negative elements, return \c NaN; print a warning when <tt>apop_opts.verbose >= 0</tt>.

Sample code:
\include entropy_vector.c
*/
long double apop_vector_entropy(gsl_vector *in){
    Apop_stopif(!in, return 0, 1, "Entropy of a NULL vector ≡ 0");
    Apop_stopif(!in->size, return 0, 1, "Entropy of a zero-length vector ≡ 0");//can't happen.

    //User may or may not have normalized in, so scale everything by the sum.
    long double sum = apop_vector_sum(in);
    Apop_stopif(sum<0, return NAN, 0, "Vector sums to a negative value (%Lg). Returning NaN.\n", sum);
    if (!sum) return 0;

    long double out=0;
    for (int i=0; i< in->size; i++){
        double val = gsl_vector_get(in, i)/sum;
        Apop_stopif(val<0, return NAN, 0, "negative value (%g) in vector position %i. Returning NaN.\n", val, i);
        if (!val) continue;
        out -= logl(val)*val;
    }
    return out;
}

static long double norment(apop_model *m){
    double sigma_sq = gsl_pow_2(apop_data_get(m->parameters, 1));
    return (log(2*M_PI*sigma_sq) +1)/2.;
}

double get_ll(apop_data *d, void *m){ return apop_log_likelihood(d, m); }


/** Calculate the entropy of a model: \f$\int -\ln(p(x))p(x)dx\f$, which is the expected
  value of \f$-\ln(p(x))\f$.

The default method is to make draws using \ref apop_model_draws, then
evaluate the log likelihood at those points using the model's \c log_likelihood method.

There are a number of routines for specific models, inlcuding the \ref apop_normal and \ref apop_pmf models.

\li  If you want the entropy of a data set, see \ref apop_vector_entropy.
\li The entropy is calculated using natural logs. If you prefer base-2 logs, just divide by \f$\ln(2)\f$: <tt>apop_model_entropy(my_model)/log(2)</tt>.

\param in A parameterized \ref apop_model. That is, you have already used \ref apop_estimate or \ref apop_model_set_parameters to estimate/set the model parameters.
\param draws If using the default method of making random draws, how many random draws to make (default=1,000)

Sample code:
\include entropy_model.c
*/
APOP_VAR_HEAD long double apop_model_entropy(apop_model *in, int draws){
    apop_model * apop_varad_var(in, NULL);
    Apop_stopif(!in, return NAN, 0, "NULL input model. Returning NaN.");
    int apop_varad_var(draws, 1000);
APOP_VAR_ENDHEAD
    static int setup=0; if (!(setup++)){
        apop_entropy_vtable_add(norment, apop_normal);
    }
    apop_entropy_type e_fn = apop_entropy_vtable_get(in);
    if (e_fn) return e_fn(in);

    apop_data *d = apop_model_draws(in, draws);
    apop_data *lls = apop_map(d, .fn_rp=get_ll, .param=in);

    long double out = -apop_vector_mean(lls->vector);
    apop_data_free(d);
    apop_data_free(lls);
    return out;
}

double a_div(gsl_vector *in){
    double pi = gsl_vector_get(in, 0);
    double qi = gsl_vector_get(in, 0);
    return pi ? pi * log(pi/qi):0;
}

/** Kullback-Leibler divergence.

This measure of the divergence of one distribution from another has the form \f$ D(p,q)
= \sum_i \ln(p_i/q_i) p_i \f$.  Notice that it is not a distance, because there is an
asymmetry between \f$p\f$ and \f$q\f$, so one can expect that \f$D(p, q) \neq D(q, p)\f$.

  \param from the \f$p\f$ in the above formula. (No default; must not be \c NULL)
  \param to the \f$q\f$ in the above formula. (No default; must not be \c NULL)
  \param draw_ct If I do the calculation via random draws, how many? (Default = 1e5)
  \param rng    A \c gsl_rng. If \c NULL or number of threads is greater than 1, I'll take care of the RNG; see \ref apop_rng_get_thread. (Default = \c NULL)

This function can take empirical histogram-type models (\ref apop_pmf) or continuous
models like \ref apop_loess or \ref apop_normal.

If there is a PMF (I'll try \c from first, under the presumption that you are measuring
the divergence of a fitted model from an observed data distribution), then I'll step
through it for the points in the summation.

\li If you have two empirical distributions in the form of \ref apop_pmf, they must
be synced: if \f$p_i>0\f$ but \f$q_i=0\f$, then the function returns \c GSL_NEGINF. If
<tt>apop_opts.verbose >=1</tt> I print a message as well.

If neither distribution is a PMF, then I'll take \c draw_ct random draws from \c from
and evaluate at those points.

\li Set <tt>apop_opts.verbose = 3</tt> for observation-by-observation info.

\li This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD long double apop_kl_divergence(apop_model *from, apop_model *to, int draw_ct, gsl_rng *rng){
    apop_model * apop_varad_var(from, NULL);
    apop_model * apop_varad_var(to, NULL);
    Apop_stopif(!from, return NAN, 0, "The first model is NULL; returning NaN.");
    Apop_stopif(!to, return NAN, 0, "The second model is NULL.");
    double apop_varad_var(draw_ct, 1e5);
    gsl_rng * apop_varad_var(rng, apop_rng_get_thread(-1));
APOP_VAR_ENDHEAD
    double div = 0;
    Apop_notify(3, "p(from)\tp(to)\tfrom*log(from/to)\n");
    if (from->name && !strcmp(from->name, "PDF or sparse matrix")){
        apop_data *p = from->data;
        apop_pmf_settings *settings = Apop_settings_get_group(from, apop_pmf);
        Get_vmsizes(p); //maxsize
        OMP_for_reduce (+:div,    int i=0; i < maxsize; i++){
            double pi = p->weights ? gsl_vector_get(p->weights, i)/settings->total_weight : 1./maxsize;
            if (!pi){
                Apop_notify(3, "0\t--\t0");
                continue;
            } //else:
            double qi = apop_p(Apop_r(p, i), to);
            Apop_notify(3,"%g\t%g\t%g", pi, qi, pi ? pi * log(pi/qi):0);
            Apop_stopif(!qi, div+=GSL_NEGINF; break, 1, "The PMFs aren't synced: from-distribution has a value where "
                                                "to-distribution doesn't (which produces infinite divergence).");
            div += pi * log(pi/qi);
        }
    } else { //the version with the RNG.
        Apop_stopif(!from->dsize, return GSL_NAN, 0, "I need to make random draws from the 'from' model, "
                                                     "but its dsize (draw size)==0. Returning NaN.");
        apop_data *draw_list = apop_data_alloc(draw_ct, 2);
        OMP_for_reduce(+:div,    int i=0; i < draw_ct; i++){
            double draw[from->dsize];
            apop_draw(draw, apop_rng_get_thread(-1), from);
            gsl_matrix_view dm = gsl_matrix_view_array(draw, 1, from->dsize);
            double pi = apop_p(&(apop_data){.matrix=&(dm.matrix)}, from);
            double qi = apop_p(&(apop_data){.matrix=&(dm.matrix)}, to);
            apop_data_set(draw_list, i, 0, pi);
            apop_data_set(draw_list, i, 1, qi);
            Apop_notify(3,"%g\t%g\t%g", pi, qi, pi ? pi * log(pi/qi):0);
            Apop_stopif(!qi, div+=GSL_NEGINF; break, 1, "From-distribution has a value where "
                                                "to-distribution doesn't (which produces infinite divergence).");
        }
        apop_vector_normalize(Apop_cv(draw_list, 0), NULL, 'p');
        apop_vector_normalize(Apop_cv(draw_list, 1), NULL, 'p');
        div = apop_map_sum(draw_list, .fn_v=a_div);
    }
    return div;
}

/** The multivariate generalization of the Gamma distribution.
\f[
\Gamma_p(a)=
\pi^{p(p-1)/4}\prod_{j=1}^p
\Gamma\left[ a+(1-j)/2\right]. \f]

Because \f$\Gamma(x)\f$ is undefined for \f$x\in\{0, -1, -2, ...\}\f$, this function returns \c NAN when \f$a+(1-j)/2\f$ takes on one of those values.

See also \ref apop_multivariate_lngamma, which is more numerically stable in most cases.
*/
long double apop_multivariate_gamma(double a, int p){
    Apop_stopif(-(a+(1-p)/2) == (int)-(a+(1-p)/2) && a+(1-p)/2 <=0, return NAN, 1, "Undefined when a + (1-p)/2 = 0, -1, -2, ... [you sent a=%g, p=%i]", a, p);
    long double out = pow(M_PI, p*(p-1.)/4.);
    long double factor = 1;
    for (int i=1; i<=p; i++)
        factor *= gsl_sf_gamma(a+(1-i)/2.);
    return out * factor;
}

/** The log of the multivariate generalization of the Gamma; see also
 \ref apop_multivariate_gamma.
*/
long double apop_multivariate_lngamma(double a, int p){
    Apop_stopif(-(a+(1-p)/2) == (int)-(a+(1-p)/2) && a+(1-p)/2 <=0, return NAN, 1, "Undefined when a + (1-p)/2 = 0, -1, -2, ... [you sent a=%g, p=%i]", a, p);
    long double out = M_LNPI * p*(p-1.)/4.;
    for (int i=1; i<=p; i++)
        out += gsl_sf_lngamma(a+(1-i)/2.);
    return out;
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

/** Test whether the input matrix is positive semidefinite (PSD).

A covariance matrix will always be PSD, so this function can tell you whether your matrix is a valid covariance matrix.

Consider the 1x1 matrix in the upper left of the input, then the 2x2 matrix in the
upper left, on up to the full matrix. If the matrix is PSD, then each of these has
a positive determinant. This function thus calculates \f$N\f$ determinants for an
\f$N\f$x\f$N\f$ matrix.

\param m The matrix to test. If \c NULL, I will return zero---not PSD.
\param semi If anything but \c 's', check for positive definite, not semidefinite. (default 's')

See also \ref apop_matrix_to_positive_semidefinite, which will change the input to something PSD.

\li This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD int apop_matrix_is_positive_semidefinite(gsl_matrix *m, char semi){
    gsl_matrix * apop_varad_var(m, NULL);
    Apop_stopif(!m, return 0, 1, "You gave me a NULL matrix. I will take this as not positive semidefinite; returning zero.");
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

/**  This function takes in a matrix and converts it in place to the `closest' positive semidefinite matrix.

\param m On input, any matrix; on output, a positive semidefinite matrix. If \c NULL, return \c NaN and print an error.
\return the distance between the original and new matrices.

\li See also the test function \ref apop_matrix_is_positive_semidefinite.
\li This function can be used as the core of a model constraint.
\li Adapted from the R Matrix package's nearPD, which is 
Copyright (2007) Jens Oehlschlägel [under the GPL].
*/
double apop_matrix_to_positive_semidefinite(gsl_matrix *m){
    Apop_stopif(!m, return NAN, 0, "Got a NULL matrix. Returning NaN.");
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
    find_eigens(&d, eigenvals, eigenvecs->matrix);//free d here.
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
