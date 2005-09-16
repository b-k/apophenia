/** \file distributions.c

Here are some functions to describe distributions. Each typically includes
a likelihood function, giving P(data), the derivative of the likelihood
function (which you will need to do an MLE, along with the fdf function
which calls both at once), and maybe a random number generator. 

At the moment, most of the headers are in likelihoods.h. Maybe some day that will change.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/

/** \defgroup mle  Estimations requiring the maximization of a likelihood function
*/
/** \defgroup network_likelihoods  Likelihood fns associated with network analysis

These functions will estimate the parameter(s) in the named distribution for the given data.

[Wikipedia has notes on the <a
href="http://en.wikipedia.org/wiki/Zipf_distribution">Zipf
distribution</a>. I'm pretty sure their specification of the Yule is
wrong; will fix when I can check references.]

\param data 	Each row of the data matrix is one observation. The first column is the number of elements with one
link, the second is the number of elements with two links, et cetera.

\param starting_pt 	A vector of the appropriate size which indicates your
best initial guess for beta. if starting_pt=NULL, then (0,0,...0) will be assumed. [Gamma, Waring: two parameters; zipf, yule: one parameter.]

\param step_size
The scale of the initial steps the maximization algorithm
will take. Currently, it is a scalar, so every dimension will have the
same step_size.

\param uses
A pointer to an \ref apop_inventory with the info you would like. You always get the parameters and the log likelihood.

\param verbose 
Zero or one depending on whether you want to see the
maximizer's iterations.

\return
Returns an \ref apop_estimate with the appropriate info.

<b>example</b><br>
Here is a simple example; see also \ref apop_make_likelihood_vector for other examples.


\code
gsl_vector *    waring_parameters;
double          starting_pt[2] = {3, 0};
double          likelihood;
waring_parameters      = mle_waring(data, &likelihood, starting_pt, .01, 0);
printf("Your most likely waring parameters are %g and %g, with likelihood %g",
                        gsl_vector_get(waring_parameter, 0) gsl_vector_get(waring_parameter, 1), likelihood);
gsl_vector_free(waring_parameter);       //Don't forget to clean up when you're done.
\endcode
\ingroup mle */


#include <apophenia/headers.h>
void prep_inventory_mle(apop_inventory *in, apop_inventory *out); //in likelihoods.c

/** The Gamma distribution

\f$G(x, a, b) 	= 1/(\Gamma(a) b^a)  x^{a-1} e^{-x/b}\f$

\f$ln G(x, a, b)= -ln \Gamma(a) - a ln b + (a-1)ln(x) + -x/b\f$

\f$d ln G/ da	=  -\psi(a) - ln b + ln(x) \f$	(also, \f$d ln \gamma = \psi\f$)

\f$d ln G/ db	=  -a/b - x \f$
\ingroup network_likelihoods
*/
double apop_gamma_likelihood(const gsl_vector *beta, void *d){
float		a	= gsl_vector_get(beta, 0),
		b	= gsl_vector_get(beta, 1);
	if (a <= 0 || b <= 0 || gsl_isnan(a) || gsl_isnan(b)) return GSL_POSINF;	
						//a sign to the minimizer to look elsewhere.
int 		i, k;
gsl_matrix	*data		= d;
float 		llikelihood 	= 0,
		ln_ga		= gsl_sf_lngamma(a),
		ln_b		= log(b),
		x;
	for (i=0; i< data->size1; i++)
		for (k=0; k< data->size2; k++){
			x		 = gsl_matrix_get(data, i, k);
			if (x!=0)
				llikelihood	+= -ln_ga - a * ln_b + (a-1) * log(x) - x/b;
		}
	return -llikelihood;
}

/** The derivative of the Gamma distribution, for use in likelihood
 * minimization. You'll probably never need to call this directy.*/
void d_gamma_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
float		a	= gsl_vector_get(beta, 0),
		b	= gsl_vector_get(beta, 1);
	//if (a <= 0 || b <= 0 || gsl_isnan(a) || gsl_isnan(b)) return GSL_POSINF;	
						//a sign to the minimizer to look elsewhere.
int 		i, k;
gsl_matrix	*data	= d;
float 		d_a 	= 0,
		d_b	= 0,
		psi_a	= gsl_sf_psi(a),
		ln_b	= log(b),
		x;
	for (i=0; i< data->size1; i++)
		for (k=0; k< data->size2; k++){
			x		 = gsl_matrix_get(data, i, k);
			if (x!=0){
				d_a	+= -psi_a - ln_b + log(x);
				d_b	+= -a/b - x;
			}
		}
	gsl_vector_set(gradient,0, -d_a);
	gsl_vector_set(gradient,1, -d_b);
}

void gamma_fdf(const gsl_vector *beta, void *d, double *f, gsl_vector *df){
	*f	= apop_gamma_likelihood(beta, d);
	d_gamma_likelihood(beta, d, df);
}

/*
gsl_vector * apop_mle_gamma(gsl_matrix *data, double *likelihood, double *starting_pt, double step_size, int verbose){
gsl_vector	*beta;
double		ll;
	ll	= maximum_likelihood_w_d(data, &beta, 2, &apop_gamma_likelihood, d_gamma_likelihood, gamma_fdf, 
							starting_pt, step_size, verbose);
	if (likelihood != NULL)	*likelihood	= ll;
	return beta;
}
*/

apop_estimate * apop_mle_gamma(gsl_matrix *data, double *starting_pt, 
					double step_size, apop_name *n, apop_inventory *uses, int verbose){
apop_inventory	actual_uses;
	prep_inventory_mle(uses, &actual_uses);
apop_estimate	*out	= apop_estimate_alloc(data->size1, 2,n,  actual_uses);
	maximum_likelihood_w_d(data, out, &apop_gamma_likelihood, d_gamma_likelihood, 
			gamma_fdf, starting_pt, step_size, verbose);
	return out;
}



//////////////////
//The probit model
//////////////////

//This section includes some trickery to avoid recalculating beta dot x.
//
gsl_vector *beta_dot_x ;
int	beta_dot_x_is_current	= 0;

void	dot(const gsl_vector *beta, gsl_matrix *data){
gsl_matrix_view p 	= gsl_matrix_submatrix(data,0,1,data->size1,data->size2-1);
	beta_dot_x 	= gsl_vector_alloc(data->size1);			//global var
        gsl_blas_dgemv (CblasNoTrans, 1.0, &p.matrix, beta, 0.0, beta_dot_x);	//dot product
}

/**
find (data dot beta'), then find the integral of the \f$\cal{N}(0,1)\f$
up to that point. Multiply likelihood either by that or by 1-that, depending 
on the choice the data made.
\ingroup mle
*/
double apop_probit_likelihood(const gsl_vector *beta, void *d){
int		i;
long double	n, total_prob	= 0;
gsl_matrix 	*data 		= (gsl_matrix *) d;		//just type casting.

	dot(beta,data);
	for(i=0;i< data->size1; i++){
		n	=gsl_cdf_gaussian_P(gsl_vector_get(beta_dot_x,i),1);
		if (gsl_matrix_get(data, i, 0)==0) 	total_prob	+= log(n);
		else 					total_prob	+= log((1 - n));
	}
	return -total_prob;
}

/** The derivative of the probit distribution, for use in likelihood
 * minimization. You'll probably never need to call this directy.*/
void d_probit_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
	//derivative of the above. 
int		i, j;
long double	one_term, beta_term_sum;
gsl_matrix *data 		= d;		//just type casting.

if (beta_dot_x_is_current==0) 	dot(beta,data);

	for(j=0; j< beta->size; j++){
		beta_term_sum	= 0;
		for(i=0; i< data->size1; i++){
			one_term	 = gsl_matrix_get(data, i,j+1)
						* gsl_ran_gaussian_pdf(gsl_vector_get(beta_dot_x,i),1);
			if (gsl_matrix_get(data, i, 0)==0) 	
				one_term	/= gsl_cdf_gaussian_P(gsl_vector_get(beta_dot_x,i),1);
			else 	one_term	/= (gsl_cdf_gaussian_P(gsl_vector_get(beta_dot_x,i),1)-1);
			beta_term_sum	+= one_term;
		}
	gsl_vector_set(gradient,j,-beta_term_sum);
	}
	gsl_vector_free(beta_dot_x);
}


void probit_fdf( const gsl_vector *beta, void *d, double *f, gsl_vector *df){
	*f	= apop_probit_likelihood(beta, d);
	beta_dot_x_is_current	=1;
	d_probit_likelihood(beta, d, df);
	beta_dot_x_is_current	=0;
}


/*
gsl_vector * apop_mle_probit(gsl_matrix *data, double *likelihood, double *starting_pt, double step_size, int verbose){
gsl_vector	*beta;
double		ll;
	ll	= maximum_likelihood_w_d(data, &beta, data->size2 - 1, apop_probit_likelihood, 
			d_probit_likelihood, probit_fdf, starting_pt, step_size, verbose);
	if (likelihood != NULL)	*likelihood	= ll;
	return beta;
}
*/

/** Does a maximum likelihood estimate of the best parameters for a probit model.
\ingroup mle

Estimate the most likey parameter(s) which fit the dependent variable to the dependent variables via a probit model.

\param data 	The first column of the data matrix is the dependent
variable, and the remaining variables are the independent. This means
that the beta which will be output will be of size <tt>(data->size2 - 1)</tt>.
\param starting_pt 	A vector of the appropriate size which indicates your
best initial guess for beta. if starting_pt=NULL, then (0,0,...0) will
be assumed.
\param step_size	The scale of the initial steps the maximization algorithm
will take. Currently, it is a scalar, so every dimension will have the
same step_size.
\param uses 	An \ref apop_inventory describing what you'd like returned. You always get the parameters and the log likelihood.
\param n	The \ref apop_name structure, if any.
\param verbose	 Zero or one depending on whether you want to see the
maximizer's iterations.

\return 	Returns an ["apop_estimate"] with the appropriate data.

<b>example</b><br>
See \ref apop_make_likelihood_vector for an example of the use of the <tt>apop_mle_xxx</tt> functions.
 */
apop_estimate * apop_mle_probit(gsl_matrix *data, double *starting_pt, 
					double step_size, apop_name *n, apop_inventory *uses, int verbose){
apop_inventory	actual_uses;
	prep_inventory_mle(uses, &actual_uses);
apop_estimate	*out		= apop_estimate_alloc(data->size1, data->size2 - 1, n, actual_uses);
	maximum_likelihood_w_d(data, out, &apop_probit_likelihood, d_probit_likelihood, probit_fdf, 
			starting_pt, step_size, verbose);
	return out;
}


/////////////////////////
//The Waring distribution
/////////////////////////
/** The Waring distribution

\f$W(x, b,a) 	= (b-1) \gamma(b+a) \gamma(k+a) / [\gamma(a+1) \gamma(k+a+b)]\f$

\f$\ln W(x, b, a) = ln(b-1) + lng(b+a) + lng(k+a) - lng(a+1) - lng(k+a+b)\f$

\f$dlnW/db	= 1/(b-1)  + \psi(b+a) - \psi(k+a+b)\f$

\f$dlnW/da	= \psi(b+a) + \psi(k+a) - \psi(a+1) - \psi(k+a+b)\f$

\ingroup network_likelihoods
*/
double apop_waring_likelihood(const gsl_vector *beta, void *d){
float		bb	= gsl_vector_get(beta, 0),
		a	= gsl_vector_get(beta, 1);
	if (bb <=2 || a < 0) return GSL_POSINF;	//a sign to the minimizer to look elsewhere.
int 		i, k;
gsl_matrix*	data		= d;
double 		ln_a_k, ln_bb_a_k,
		likelihood 	= 0,
		ln_bb_a		= gsl_sf_lngamma(bb + a),
		ln_a_mas_1	= gsl_sf_lngamma(a + 1),
		ln_bb_less_1	= log(bb - 1);
	for (k=0; k< data->size2; k++){	//more efficient to go column-by-column
		ln_bb_a_k	 = gsl_sf_lngamma(k  + a + bb);
		ln_a_k		 = gsl_sf_lngamma(k  + a);
		for (i=0; i< data->size1; i++){
			likelihood	+= gsl_matrix_get(data, i, k) *(ln_bb_less_1 + ln_a_k + ln_bb_a - ln_a_mas_1 - ln_bb_a_k);
		}
	}
	return -likelihood;
}

/** The derivative of the Waring distribution, for use in likelihood
 * minimization. You'll probably never need to call this directy.*/
void d_waring_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
	//Psi is the derivative of the log gamma function.
float		bb		= gsl_vector_get(beta, 0),
		a		= gsl_vector_get(beta, 1);
int 		i, k;
gsl_matrix	*data		= d;
double		bb_minus_one_inv= 1/(bb-1),
		psi_a_bb	= gsl_sf_psi(bb + a),
		psi_a_mas_one	= gsl_sf_psi(a+1),
		psi_a_k,
		psi_bb_a_k,
		d_bb		= 0,
		d_a		= 0;
	for (k=0; k< data->size2; k++){	//more efficient to go column-by-column
		psi_bb_a_k	 = gsl_sf_psi(k  + a + bb);
		psi_a_k		 = gsl_sf_psi(k  + a);
		for (i=0; i< data->size1; i++){
			d_bb		+= gsl_matrix_get(data, i, k) *(bb_minus_one_inv + psi_a_bb - psi_bb_a_k);
			d_a		+= gsl_matrix_get(data, i, k) *(psi_a_bb + psi_a_k - psi_a_mas_one - psi_bb_a_k);
		}
	}
	gsl_vector_set(gradient, 0, -d_bb);
	gsl_vector_set(gradient, 1, -d_a);
}

void waring_fdf(const gsl_vector *beta, void *d, double *f, gsl_vector *df){
	*f	= apop_waring_likelihood(beta, d);
	d_waring_likelihood(beta, d, df);
}

apop_estimate * apop_mle_waring(gsl_matrix *data, double *starting_pt, 
					double step_size, apop_name *n, apop_inventory *uses, int verbose){
apop_inventory	actual_uses;
	prep_inventory_mle(uses, &actual_uses);
apop_estimate	*out		= apop_estimate_alloc(data->size1, 2, n, actual_uses);
	maximum_likelihood_w_d(data, out, &apop_waring_likelihood, d_waring_likelihood, waring_fdf, 
				starting_pt, step_size, verbose);
	return out;
}



///////////////////////
//The Yule distribution
///////////////////////
//Yule likelihood fn. The special case of Waring where alpha = 0.
//Y(x, b) 	= (b-1) g(b) g(k) / g(k+b), where g(x) =gamma(x)
//ln Y(x, b) 	= ln(b-1) + lng(b) + lng(k) - lng(k+b)
//dlnY/db	= 1/(b-1)  + psi(b) - psi(k+b)

double apop_yule_likelihood(const gsl_vector *beta, void *d){
float		bb		= gsl_vector_get(beta, 0);
	if (bb <=2) return GSL_POSINF;	//a sign to the minimizer to look elsewhere.
int 		i, k;
gsl_matrix 	*data	= d;
float 		ln_k, ln_bb_k,
	likelihood 	= 0,
	ln_bb		= gsl_sf_lngamma(bb),
	ln_bb_less_1	= log(bb-1);
	for (k=0; k< data->size2; k++)	{
		if (k>1) 	ln_k	= gsl_sf_lngamma(k);
		else		ln_k	= 0;
		ln_bb_k		= gsl_sf_lngamma(k+bb);
		for (i=0; i< data->size1; i++){
			likelihood	+= gsl_matrix_get(data,i,k) * (ln_bb_less_1 + ln_k + ln_bb - ln_bb_k);
		}
	}
	return -likelihood;
}

void d_yule_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
	//Psi is the derivative of the log gamma function.
float		bb		= gsl_vector_get(beta, 0);
int 		i, k;
gsl_matrix	*data		= d;
double		bb_minus_one_inv= 1/(bb-1),
		psi_bb		= gsl_sf_psi(bb),
		psi_bb_k,
		d_bb		= 0;
	for (k=0; k< data->size2; k++){
		psi_bb_k	 = gsl_sf_psi(k + bb);
		for (i=0; i< data->size1; i++){
			d_bb		+= gsl_matrix_get(data, i, k) *(bb_minus_one_inv + psi_bb - psi_bb_k);
		}
	}
	gsl_vector_set(gradient, 0, -d_bb);
}

void yule_fdf(const gsl_vector *beta, void *d, double *f, gsl_vector *df){
	*f	= apop_yule_likelihood(beta, d);
	d_yule_likelihood(beta, d, df);
}

apop_estimate * apop_mle_yule(gsl_matrix *data, double *starting_pt, 
					double step_size, apop_name *n, apop_inventory *uses, int verbose){
apop_inventory	actual_uses;
	prep_inventory_mle(uses, &actual_uses);
apop_estimate	*out		= apop_estimate_alloc(data->size1, 1, n, actual_uses);
	maximum_likelihood_w_d(data, out, &apop_yule_likelihood, d_yule_likelihood, yule_fdf, 
			starting_pt, step_size, verbose);
	return out;
}

/** This function is used to keep the minimizer away from bounds.

If you just return GSL_POSINF at the bounds, it's not necessarily smart enough to get it.
This helps it along.	*/
double keep_away(double value, double limit, char type, double base){
double		out;
	if (type == 't') 	out	= exp(limit - value);
	else 			out	= exp(value - limit);
	return out * base;
}

///////////////////////
//The Zipf distribution
///////////////////////
#include <gsl/gsl_sf_zeta.h>

/** The Zipf takes no data. It can't really be used to determine the likelihood
of data by itself; you need priors on \f$a\f$ to do that.
 */ 
double apop_zipf_likelihood(double a, int i){
	if (gsl_isnan(a))	return 0;
	return 1/(gsl_sf_zeta(a) * pow(i, a));
}

/** Find the mean squared error between the Zipf with the value beta[0] and the data's
distribution.*/
double apop_zipf_mse(const gsl_vector *beta, void *d){
double		bb	= gsl_vector_get(beta, 0);
static double	ka	= 0;
	if (bb <= 1) {
		if (ka ==0){
			gsl_vector *	b_ka	= gsl_vector_alloc(1);
			gsl_vector_set(b_ka,0, 1.01);
		 	ka	= apop_zipf_mse(b_ka , d);
			gsl_vector_free (b_ka);
		}
		return keep_away(bb, 1, 'b', ka);
	}
gsl_matrix	*data	= d;
int 		i, j;
long double	mse	= 0, z,
     		sub,
		denom	= 1.0/(data->size1 * data->size2);
	for(j=0; j< data->size2; j++){
		sub	= 0;
		for(i=0; i< data->size1; i++){
			z	 = apop_zipf_likelihood(bb, j+1);
			sub	+= gsl_pow_int(gsl_matrix_get(data,i,j) - z,2);
		}
		mse	+= sub *denom;
	}
printf("z: %g %Lg", bb, mse);
	return mse;
}	

/** MSE: \f$ \sum_i (x_i - {1\over \zeta(a) i^a})^2 = \sum_i (x_i - {2 x_i \over \zeta(a) i^a} + \left({1 \over \zeta(a) i^a}\right)^2)\f$
dMSE/da: \f$ 2x_i \zeta(a)^{-2} d\zeta(a) i^a + 2 x_i a \zeta(a)^-1 i ^{-(a+1)} - 2 \zeta(a)^{-3} (d\zeta(a))^{-2} i^{2a}
    - 2a \zeta(a)^{-2} i ^{-(2a+1)} \f$

where \f$d\zeta(a) = -a \zeta(a-1)\f$.
    */
void apop_zipf_d_mse(const gsl_vector *beta, void *d, gsl_vector *gradient){
double		a	= gsl_vector_get(beta, 0);
static double	ka	= 0;
	if (a <= 2) {
		if (ka ==0){
			gsl_vector 	*b_ka	= gsl_vector_alloc(1);
			gsl_vector_set(b_ka,0, 1.01);
		 	ka	= apop_zipf_mse(b_ka , d);
			gsl_vector_free (b_ka);
		}
		gsl_vector_set(gradient,0,keep_away(a, 2, 'b', ka));
		return;
	}
gsl_matrix	*data	= d;
int 		i, j;
double		dmse	= 0, 
		z, dz, x, i_to_a, i_to_a1;
	for(j=0; j< data->size2; j++)
		for(i=0; i< data->size1; i++){
			x	 = gsl_matrix_get(data,i,j);
			z	 = apop_zipf_likelihood(a, j+1);
			dz	 = -a * apop_zipf_likelihood(a-1, j+1);
			i_to_a	 = pow(j,a);
			i_to_a1	 = pow(j,a+1);
			dmse	+= 2* x/gsl_pow_int(z,2) * dz * i_to_a + 2 *x * a /z / i_to_a1  
				- 2 / gsl_pow_int(z,3) / gsl_pow_int(dz,2) * pow(j+1,2*a)
    				- 2*a /gsl_pow_int(z,2)* pow(j+1,-(2*a+1));
			if (gsl_isnan(dmse))
				printf("break here.");
		}
	dmse	/= (data->size1 * data->size2);
printf("dz: %g %g", a, dmse);
	gsl_vector_set(gradient,0,-dmse);
}	

void zipf_fdf_mse(const gsl_vector *beta, void *d, double *f, gsl_vector *df){
	*f	= apop_zipf_mse(beta, d);
	apop_zipf_d_mse(beta, d, df);
}


/** Find the value of the Zipf distribution's single parameter which
minimizes mean squared error. 

For the data I've got, it doesn't work right now.
\todo Fix this fn. */
apop_estimate * apop_zipf_min_mse(gsl_matrix *data, double *starting_pt, 
					double step_size, apop_name *n, apop_inventory *uses, int verbose){
apop_inventory	actual_uses;
	prep_inventory_mle(uses, &actual_uses);
apop_estimate	*out	= apop_estimate_alloc(data->size1, 1, n, actual_uses);
	maximum_likelihood_w_d(data, out, &apop_zipf_mse, &apop_zipf_d_mse, &zipf_fdf_mse, starting_pt, step_size, verbose);
	return out;
}
/*
double apop_zipf_min_mse(gsl_matrix *data, double *starting_pt, 
					double step_size, apop_name *n, apop_inventory *uses, int verbose){
apop_inventory	actual_uses;
	prep_inventory_mle(uses, &actual_uses);
gsl_vector	* betas;
	maximum_likelihood(data, &betas, 1, apop_zipf_mse, starting_pt, step_size, verbose);
	return gsl_vector_get(betas,0);
}
*/

/** The exponential distribution. A one-parameter likelihood fn.

Wikipedia calls something else the exponential distribution; one day
I'll have to get a reasonable sampling of the lit and resolve this.

 
\f$P(link ct==k) = C^{-k}   	\f$. The PDF is \f$ln(C) C^{-k}\f$, which sums to one.<br>
\f$Z(C,k) 	= ln(C) C^{-k} 	\f$ <br>
\f$ln Z(C,k) 	= -ln(ln(C)) * ln(C) * k	\f$ <br>
\f$dlnZ/dC	= -ln(ln(C)) * ln(C) * k	\f$ <br>
\f$d^2lnZ/dC^2	= k / C^2	\f$ <br>

\ingroup network_likelihoods
*/
double apop_exponential_likelihood(const gsl_vector *beta, void *d){
float		bb	= gsl_vector_get(beta, 0);
//double apop_zipf_likelihood(double bb, void *d){
	if (bb <=1) return GSL_POSINF;	//a sign to the minimizer to look elsewhere.
int 		i, k;
gsl_matrix	*data		= d;
double 		llikelihood 	= 0,
		ln_c		= log(bb),
		ln_ln_c		= log(ln_c);
	for (i=0; i< data->size1; i++)
		for (k=0; k< data->size2; k++)
//			llikelihood	+= gsl_matrix_get(data, i, k) * (ln_c * k);
			//llikelihood	+= -gsl_matrix_get(data, i, k) * ln_c ;
			llikelihood	+=  -ln_ln_c * ln_c * gsl_matrix_get(data, i, k);
	return -llikelihood;
printf("l: %g ", llikelihood);
	return -llikelihood;
}


void d_exponential_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
double		bb		= gsl_vector_get(beta, 0);
int 		i, k;
gsl_matrix	*data		= d;
double 		d_likelihood 	= 0,
		ln_c		= log(bb),
		ln_ln_c		= log(ln_c);
	for (i=0; i< data->size1; i++)
		for (k=0; k< data->size2; k++){
			//d_likelihood	+= gsl_matrix_get(data, i, k)  *  -k/bb;
			//d_likelihood	+= -gsl_matrix_get(data, i, k)  /bb;
			d_likelihood	+= -(ln_ln_c +1) * gsl_matrix_get(data, i, k)/bb;
		}
printf("d: %g ", -d_likelihood);
	gsl_vector_set(gradient,0, -d_likelihood);
}

void dd_exponential_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
double		bb		= gsl_vector_get(beta, 0);
int 		i, k;
gsl_matrix	*data		= d;
double 		dd_likelihood 	= 0;
	for (i=0; i< data->size1; i++)
		for (k=0; k< data->size2; k++){
			dd_likelihood	+= -gsl_matrix_get(data, i, k) / gsl_pow_int(bb,2);
			//dd_likelihood	+= gsl_matrix_get(data, i, k)  *  k/gsl_pow_int(bb,2);
		}
	gsl_vector_set(gradient,0, -dd_likelihood);
}

void exponential_fdf(const gsl_vector *beta, void *d, double *f, gsl_vector *df){
	*f	= apop_exponential_likelihood(beta, d);
	d_exponential_likelihood(beta, d, df);
}

apop_estimate * apop_mle_exponential(gsl_matrix *data, double *starting_pt, 
					double step_size, apop_name *n, apop_inventory *uses, int verbose){
apop_inventory	actual_uses;
	prep_inventory_mle(uses, &actual_uses);
apop_estimate	*out	= apop_estimate_alloc(data->size1, 1, n, actual_uses);
	maximum_likelihood_w_d(data, out, &apop_exponential_likelihood, d_exponential_likelihood, exponential_fdf, 
			starting_pt, step_size, verbose);
	return out;
}

/** Draw from a Zipf distribution with parameter a

For example:
\code
gsl_rng *       r;
gsl_rng_env_setup();
r=gsl_rng_alloc(gsl_rng_taus);	//for example. 
apop_rng_zipf(r, 1.4);
\endcode

Cribbed from <a href="http://cgm.cs.mcgill.ca/~luc/mbookindex.html>Devroye (1986)</a>, p 551.  */
double apop_rng_zipf(gsl_rng* r, double a){
if (a  <= 1)	
	{printf("apop_rng_zipf: Zipf needs a parameter >=1. Returning 0.\n"); return 0;};
double		x;
double		u, v, t, 
		b 	= pow(2, a-1), 
		ainv	= -(1.0/(a-1));
	do {
		u	= gsl_rng_uniform(r);
		v	= gsl_rng_uniform(r);
		x	=  pow(u, ainv);
		t	= pow((1.0 + 1.0/x), (a-1));
	} while (v * x * (t-1.0)/(b-1) > t/b || x < 0);
	return x;
}
