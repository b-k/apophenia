#include "stats.h"

inline double mean(gsl_vector *in){
	return gsl_stats_mean(in->data,in->stride, in->size); }

inline double var(gsl_vector *in){
	return gsl_stats_variance(in->data,in->stride, in->size); }

inline double kurtosis(gsl_vector *in){
	return gsl_stats_kurtosis(in->data,in->stride, in->size); }

inline double kurt(gsl_vector *in){
	return gsl_stats_kurtosis(in->data,in->stride, in->size); }

	/*
inline double var_view(gsl_vector_view *in){
	return gsl_stats_variance(in->data,in->stride, in->size); }
	*/

inline double var_m(gsl_vector *in, double mean){
	return gsl_stats_variance_m(in->data,in->stride, in->size, mean); }

inline double covar(gsl_vector *ina, gsl_vector *inb){
	return gsl_stats_covariance(ina->data,ina->stride,inb->data,inb->stride,inb->size); }

inline double cov(gsl_vector *ina, gsl_vector *inb){return  covar(ina,inb);}

void normalize_vector(gsl_vector *in, gsl_vector **out){
double		mu,sigma;
int		i;
	*out 	=gsl_vector_alloc (in->size);
	mu	=mean(in);
	sigma	=sqrt(var_m(in, mu));
	for (i=in->size;i--;)
		gsl_vector_set(*out,i,(gsl_vector_get(in,i)-mu)/sigma);
	}

void normalize_data_matrix(gsl_matrix *data){
gsl_vector_view v;
double          m;
int             j,k;
        for (j = 0; j < data->size2; j++){
                v               = gsl_matrix_column(data, j);
                m               = mean(&(v.vector));
                for (k = 0; k < data->size1; k++)
                        gsl_matrix_set(data, k, j, gsl_matrix_get(data, k,j) - m);
        }
}


inline double test_chi_squared_var_not_zero(gsl_vector *in){
gsl_vector	*normed;
int		i;
double 		sum=0;
	normalize_vector(in,&normed);
	gsl_vector_mul(normed,normed);
	for(i=0;i< normed->size; 
			sum +=gsl_vector_get(normed,i++));
	return gsl_cdf_chisq_P(sum,in->size); }

inline double double_abs(double a) {if(a>0) return a; else return -a;}

void view_matrix(gsl_matrix *a){
int i,j;
	for(i=0;i< a->size1; i++){
		for(j=0;j< a->size2; j++)
			printf("%g ", gsl_matrix_get(a,i,j));
		printf("\n");
	}
}

/*
gsl_matrix *dot_prod(gsl_matrix *a, gsl_matrix *b){
int i,j, a1,a2, b1, b2;
	a1=a->size1; a2=a->size2;
	b1=b->size1; b2=b->size2;
	if(a2!= b1){ 
		printf("you're asking me to take the dot product of a %i X %i and a %i X %i matrix. With all due respect, no.\n", a1, a2,b1,b2);
		return a;}
	//Else:
	for (i=0;i<a1;i++)

}
	*/


double randombeta(double m, double v, gsl_rng *r) {
	/*Give me mean m and variance v, and I'll give you
	 * n draws from the appropriate beta dist.
	 * remember: 0<m<1, and v is tiny (<<1/12). You get NaNs if no
	 * appropriate distribution exists.*/
double 		k        = (m * (1- m)/ v) -1 ;
        return gsl_ran_beta(r, m* k ,  k*(1 - m) );
}




double multivariate_normal_prob(gsl_vector *x, gsl_vector* mu, gsl_matrix* sigma, int first_use){
	//Evaluate a multivariate normal(mu, sigma) at the point x.
//The equation:
//	exp(-1/2 (X-mu)' sigma^-1 (x-mu))
//	--------------------------
//	sqrt((2 Pi)^n det(sigma))
//
//The inverse and determinant are expensive, so keep them around where possible: on the first call, set 
//first_use to 1, then feed in as many new values of X as you want.

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
		determinant	= det_and_inv(sigma, inverse, 1,1);
	}
	if (determinant == 0) {printf("x"); return(GSL_NEGINF);} //tell minimizer to look elsewhere.
	numerator	= exp(- x_prime_sigma_x(x_minus_mu, inverse) / 2);
printf("(%g %g %g)", numerator, x_prime_sigma_x(x_minus_mu, inverse), (numerator / pow(2 * M_PI, (float)dimensions/2) * sqrt(determinant)));
	return(numerator / pow(2 * M_PI, (float)dimensions/2) * sqrt(determinant));
}
