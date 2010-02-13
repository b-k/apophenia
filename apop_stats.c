/** \file apop_stats.c	Basic moments and some distributions. */
/* Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "db.h"     //just for apop_opts
#include "asst.h" //rng_alloc
#include "stats.h"
#include <gsl/gsl_rng.h>


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
  apop_assert(in, 0, 1,'c', "You just asked me to sum a NULL. Returning zero.")
  double  out = 0;
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
  //I reimplement the skew calculation here without the division by var^3/2
  //that the GSL does. 
  //Lawyers may want to know that this code is cut/pasted/modified from the GSL. 

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
 
  Some people like to normalize the skew by dividing by variance squared, or by subtracting three; those things are  not done here, so you'll have to do them separately if need be.
\ingroup vector_moments
*/
double apop_vector_kurtosis_pop(const gsl_vector *in){
  //I reimplement the kurtosis calculation here without the division by var^2
  //that the GSL does. 
  //Lawyers may want to know that this code is cut/pasted/modified from the GSL. 

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


  This does not normalize the output: the kurtosis of a \f${\cal N}(0,1)\f$ is three \f$\sigma^4\f$, not three, one, or zero.
\ingroup vector_moments
*/
double apop_vector_kurtosis(const gsl_vector *in){
  size_t n = in->size;
  double coeff1 = gsl_pow_3(n)/(n-1)/(gsl_pow_2(n)-3*n+3);
  double coeff2 = (6*n-9)/(gsl_pow_2(n)-3*n+3);
    return  coeff1 * apop_vector_kurtosis_pop(in) + coeff2 * gsl_pow_2(apop_vector_var(in)*(n-1.)/n);
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

\include eg/test_distances.c

This function uses the \ref designated syntax for inputs.
\ingroup convenience_fns
*/
APOP_VAR_HEAD double apop_vector_distance(const gsl_vector *ina, const gsl_vector *inb, const char metric, const double norm){
    static gsl_vector *zero = NULL;
    const gsl_vector * apop_varad_var(ina, NULL);
    apop_assert(ina, 0, 0, 's', "The first vector has to be non-NULL.");
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
    return apop_vector_distance_base(ina, inb, metric, norm);
APOP_VAR_ENDHEAD
  apop_assert(ina->size == inb->size, 0, 0,'s', 
                "I need equal-sized vectors, but "
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
        for (i=0; i< ina->size; i++){
            dist    += pow(gsl_vector_get(ina, i) - gsl_vector_get(inb, i), norm);
        }
        return pow(dist, 1./norm); 
    }
  apop_error(0,'s', "I couldn't find the metric type you gave, %c, in my list of supported types.", metric);
  return 0; //just to keep the compiler quiet.
}

/** Returns the scalar Manhattan metric distance  between two vectors. Simply \f$\sum_i{|a_i - b_i|},\f$
where \f$i\f$ iterates over dimensions.

Equivalent to \ref apop_vector_distance<tt>(ina, inb, .metric='M')</tt>.

\ingroup convenience_fns
*/
double apop_vector_grid_distance(const gsl_vector *ina, const gsl_vector *inb){
  apop_assert(ina, 0, 0, 'c', "first input vector is NULL. Returning 0.");
  apop_assert(inb, 0, 0, 'c', "second input vector is NULL. Returning 0.");
  apop_assert(ina->size == inb->size, 0, 0,'s', 
                "You sent a vector of size %zu and a vector of size %zu.", ina->size, inb->size);
  double  dist    = 0;
    for (size_t i=0; i< ina->size; i++)
        dist    += fabs(gsl_vector_get(ina, i) - gsl_vector_get(inb, i));
	return dist; 
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
      apop_assert_void(in, 1, 'c', "Input vector is NULL. Doing nothing.\n");
      gsl_vector ** apop_varad_var(out, NULL);
      const char  apop_varad_var(normalization_type, 'p');
      return apop_vector_normalize_base(in, out, normalization_type);
APOP_VAR_END_HEAD
  double		mu, min, max;
	if (!out) 	
		out	= &in;
	else {
		*out 	= gsl_vector_alloc (in->size);
		gsl_vector_memcpy(*out,in);
	}
        //the numbers are deprecated and will go away.
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
		mu	= apop_vector_mean(in);
		gsl_vector_scale(*out, 1/(mu * in->size));	
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
    apop_assert_void(data, 1, 'c', "input matrix is NULL. Doing nothing.");
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

/** Gives a random double between min and max [inclusive].

This function uses the \ref designated syntax for inputs. Notice that calling this function with no arguments, 

\code
apop_random_double()
\endcode
conveniently produces a number between zero and one. [To do this with less overhead, allocate your own RNG and use \c gsl_ran_uniform(r).]

\param min      Default = 0
\param  max 	Default = 1
\param rng    A \c gsl_rng. If NULL, I'll take care of the RNG; see \ref autorng. (Default = \c NULL)
*/
APOP_VAR_HEAD double apop_random_double(double min, double max, gsl_rng *rng){
    static gsl_rng * spare_rng = NULL;
    double apop_varad_var(min, 0);
    double apop_varad_var(max, 1);
    gsl_rng * apop_varad_var(rng, NULL);
    if (!rng && !spare_rng) 
        spare_rng = apop_rng_alloc(++apop_opts.rng_seed);
    if (!rng)  rng = spare_rng;
    return apop_random_double_base(min, max, rng);
APOP_VAR_ENDHEAD
  double		base = gsl_rng_uniform(rng);
	return base * (max - min) + min;
}

/** Gives a random integer between min and max [inclusive].

\param min  (default 0)
\param max 	(default 1)
\param rng    A \c gsl_rng. If NULL, I'll take care of the RNG; see \ref autorng. (Default = \c NULL)

Thus,
\code
x = apop_random_int()
\endcode
makes a binary zero-one draw, and
\code
data fivepoints[] = {1, 2, 3, 5, 7};
y = apop_random_int(0, 4)
x = apop_random_int(.max=4)
\endcode
gives two draws from a five-item vector. Notice that the max is the largest index, which is one minus the dimension.
*/
APOP_VAR_HEAD int apop_random_int(double min, double max, const gsl_rng *rng){
    static gsl_rng * spare_rng = NULL;
    double apop_varad_var(min, 0);
    double apop_varad_var(max, 1);
    const gsl_rng * apop_varad_var(rng, NULL);
    if (!rng && !spare_rng) 
        spare_rng = apop_rng_alloc(++apop_opts.rng_seed);
    if (!rng)  rng = spare_rng;
    return apop_random_int_base(min, max, rng);
APOP_VAR_ENDHEAD
  double		base = gsl_rng_uniform(rng);
	return (int) (base * (max - min + 1) + min);
}

/** Returns the sum of the elements of a matrix. Occasionally convenient.

  \param m	the matrix to be summed. 
\ingroup convenience_fns*/
long double apop_matrix_sum(const gsl_matrix *m){
  apop_assert(m, 0, 1,'c', "You just asked me to sum a NULL. Returning zero.")
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

/** Returns the variance of all elements of a matrix, given the
mean. If you want to calculate both the mean and the variance, use \ref
apop_matrix_mean_and_var.

\param data	the matrix to be averaged. 
\param mean	the pre-calculated mean
\ingroup convenience_fns*/
double apop_matrix_var_m(const gsl_matrix *data, double mean){
  double     avg2    = 0;
  int        cnt= 0;
  double     x, ratio;
    for(size_t i=0; i < data->size1; i++)
        for(size_t j=0; j < data->size2; j++){
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
  size_t         cnt= 0;
  long double x, ratio;
    for(size_t i=0; i < data->size1; i++)
        for(size_t j=0; j < data->size2; j++){
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

/** Put summary information about the columns of a table (mean, std dev, variance, min, median, max) in a table.

\param indata The table to be summarized. An \ref apop_data structure.
\return     An \ref apop_data structure with one row for each column in the original table, and a column for each summary statistic. May have a <tt>weights</tt> element.
\ingroup    output

\li \ref This function gives more columns than you probably want; use apop_data_prune_columns to pick the ones you want to see.
\todo We should probably let this summarize rows as well.  */
apop_data * apop_data_summarize(apop_data *indata){
  apop_assert(indata, NULL, 0, 'c', "You sent me a NULL apop_data set. Returning NULL.\n");
  apop_assert(indata->matrix, NULL, 0, 'c', "You sent me an apop_data set with a NULL matrix. Returning NULL.\n");
  apop_data	*out	= apop_data_alloc(0,indata->matrix->size2, 6);
  double		mean, stddev,var;
  char		rowname[10000]; //crashes on more than 10^9995 columns.
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
        double *pctiles = apop_vector_percentiles(v);
		gsl_matrix_set(out->matrix, i, 0, mean);
		gsl_matrix_set(out->matrix, i, 1, stddev);
		gsl_matrix_set(out->matrix, i, 2, var);
		gsl_matrix_set(out->matrix, i, 3, pctiles[0]);
		gsl_matrix_set(out->matrix, i, 4, pctiles[50]);
		gsl_matrix_set(out->matrix, i, 5, pctiles[101]);
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
    if (!w)
        return apop_vector_mean(v);
    apop_assert(v,  0, 0, 'c', "data vector is NULL. Returning zero.\n");
    apop_assert(v->size,  0, 1, 'c', "data vector has size 0. Returning zero.\n");
    apop_assert(w->size == v->size,  0, 0,'c', "data vector has size %zu; weighting vector has size %zu. Returning zero.\n", v->size, w->size);
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
  long double   sum = 0, wsum = 0, sumsq = 0, vv, ww;
    if (!w)
        return apop_vector_var(v);
    apop_assert(v,  0, 0, 'c', "data vector is NULL. Returning zero.\n");
    apop_assert(v->size,  0, 0,'c', "data vector has size 0. Returning zero.\n");
    apop_assert(w->size == v->size,  0, 0,'c', "data vector has size %zu; weighting vector has size %zu. Returning zero.\n", v->size, w->size);
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

static double skewkurt(const gsl_vector *v, const gsl_vector *w, const int exponent, const char *fn_name){
  int           i;
  long double   wsum = 0, sumcu = 0, vv, ww, mu;
    if (!w)
        return exponent ==3 ? apop_vector_skew(v) : apop_vector_kurtosis(v);
    apop_assert(v,  0, 0, 'c', "%s: data vector is NULL. Returning zero.\n", fn_name);
    apop_assert(v->size,  0, 1, 'c',"%s: data vector has size 0. Returning zero.\n", fn_name);
    apop_assert(w->size == v->size,  0, 1, 'c',"%s: data vector has size %zu; weighting vector has size %zu. Returning zero.\n", fn_name, v->size, w->size);
    //Using the E(x - \bar x)^3 form, which is lazy.
    mu  = apop_vector_weighted_mean(v, w);
    for (i=0; i< w->size; i++){
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
    return skewkurt(v,w,3, "apop_vector_weighted_skew");
}

/** Find the population kurtosis of a weighted vector.

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
  long double   sum1 = 0, sum2 = 0, wsum = 0, sumsq = 0, vv1, vv2, ww;
    if (!w)
        return apop_vector_cov(v1,v2);
    apop_assert(v1,  0, 0, 'c', "first data vector is NULL. Returning zero.\n");
    apop_assert(v2,  0, 0, 'c', "second data vector is NULL. Returning zero.\n");
    apop_assert(v1->size,  0, 1, 'c', "apop_vector_weighted_variance: data vector has size 0. Returning zero.\n");
    apop_assert((w->size == v1->size) && (w->size == v2->size),  0, 0, 'c', "apop_vector_weighted_variance: data vectors have sizes %zu and %zu; weighting vector has size %zu. Returning zero.\n", v1->size, v2->size, w->size);
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
    apop_assert(in,  NULL, 0, 'c', "Input matrix is NULL. Returning same.");
    const char apop_varad_var(normalize, 0)
    return apop_matrix_covariance_base(in, normalize);
APOP_VAR_ENDHEAD
  gsl_matrix	*out;
  double		means[in->size2];
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

This is the \c gsl_matrix  version of \ref apop_data_covariance; if you have column names, use that one.

\param in 	A data matrix: rows are observations, columns are variables. (No default, must not be \c NULL)
\param normalize
'n' or 'N' = subtract the mean from each column, thus changing the input data but speeding up the computation.<br>
anything else (like 0)= don't modify the input data (default = no modification)

\return Returns the variance/covariance matrix relating each column with each other. This function allocates the matrix for you.

This function uses the \ref designated syntax for inputs.
\ingroup matrix_moments */
APOP_VAR_HEAD gsl_matrix *apop_matrix_correlation(gsl_matrix *in, const char normalize){
    gsl_matrix *apop_varad_var(in, NULL)
    apop_assert(in,  NULL, 0, 'c', "Input matrix is NULL; returning NULL.");
    const char apop_varad_var(normalize, 0)
    return apop_matrix_correlation_base(in, normalize);
APOP_VAR_ENDHEAD
  gsl_matrix      *out    = apop_matrix_covariance(in, normalize);
  double          std_dev;
    for(size_t i=0; i< in->size2; i++){
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
  double      var;
    for (size_t i=0; i < in->matrix->size2; i++){
        for (size_t j=i; j < in->matrix->size2; j++){
            APOP_COL(in, i, v1);
            APOP_COL(in, j, v2);
            var = apop_vector_weighted_cov(v1, v2, in->weights);
            gsl_matrix_set(out->matrix, i,j, var);
            if (i!=j)
                gsl_matrix_set(out->matrix, j,i, var);
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
\ingroup matrix_moments */
apop_data *apop_data_correlation(const apop_data *in){
  apop_data *out = apop_data_covariance(in);
  double    std_dev;
    for(size_t i=0; i< in->matrix->size2; i++){
        APOP_COL(in, i, cvin);
        APOP_COL(out, i, cvout);
        APOP_ROW(out, i, rvout);
        std_dev     = sqrt(apop_vector_weighted_var(cvin,in->weights));
        gsl_vector_scale(cvout, 1.0/std_dev);
        gsl_vector_scale(rvout, 1.0/std_dev);
    }
    return out;
}
