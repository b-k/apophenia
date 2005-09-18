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


Oh, and guess what: this documentation is out of date right now, and refers at points to how the MLEs were done before version 0.12.
Will be fixed shortly.

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
\ingroup mle 
\todo Update this documentation.*/


#include <apophenia/headers.h>
void prep_inventory_mle(apop_inventory *in, apop_inventory *out); //in likelihoods.c

/** The Gamma distribution

\f$G(x, a, b) 	= 1/(\Gamma(a) b^a)  x^{a-1} e^{-x/b}\f$

\f$ln G(x, a, b)= -ln \Gamma(a) - a ln b + (a-1)ln(x) + -x/b\f$

\f$d ln G/ da	=  -\psi(a) - ln b + ln(x) \f$	(also, \f$d ln \gamma = \psi\f$)

\f$d ln G/ db	=  -a/b - x \f$
\ingroup network_likelihoods
*/
double apop_gamma_log_likelihood(const gsl_vector *beta, void *d){
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
void apop_gamma_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
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
double apop_probit_log_likelihood(const gsl_vector *beta, void *d){
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
void apop_probit_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
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


/** Saves some time in calculating both log likelihood and dlog
likelihood for probit.	*/
void apop_probit_fdf( const gsl_vector *beta, void *d, double *f, gsl_vector *df){
	*f	= apop_probit_log_likelihood(beta, d);
	beta_dot_x_is_current	=1;
	apop_probit_dlog_likelihood(beta, d, df);
	beta_dot_x_is_current	=0;
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
double apop_waring_log_likelihood(const gsl_vector *beta, void *d){
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
 minimization. You'll probably never need to call this directy.*/
void apop_waring_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
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


double apop_yule_log_likelihood(const gsl_vector *beta, void *d){
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

void apop_yule_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
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


/** Draw from a Yule distribution with parameter a

Call this fn using \ref apop_yule.rng().

\param	a	The parameter.
\param	r	A gsl_rng that you've already set up.

For example:
\code
gsl_rng *       r;
gsl_rng_env_setup();
r=gsl_rng_alloc(gsl_rng_taus);	//for example. 
apop_yule_rng(r, 1.4);
\endcode

Cribbed from <a href="http://cgm.cs.mcgill.ca/~luc/mbookindex.html>Devroye (1986)</a>, p 553.  */
double apop_yule_rng(gsl_rng * r, double a){
double 		e1, e2;
int		x;
	e1	= gsl_ran_exponential(r, 1);
	e2	= gsl_ran_exponential(r, 1);
	x	= - e1  / log(1 - exp(-e2 / (a -1)));
	return  x + 1;	//we rounded down to floor, but want ceil.
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

/** The Zipf distribution.

\f$Z(a)		= {1\over \zeta(a) * i^a}		\f$<br>

 \todo link this fn in with the object */ 
double apop_zipf_likelihood(double a, int i){
double		z	= 1/(gsl_sf_zeta(a) * pow(i, a));
	return z;
}

double apop_zipf_log_likelihood(const gsl_vector *beta, void *d){
double		bb	= gsl_vector_get(beta, 0);
static double	ka	= 0;
gsl_matrix	*data	= d;
int 		i, j;
double		like	= 0, z;
	if (bb <= 1) {		//run away
		if (ka ==0){
			gsl_vector *	b_ka	= gsl_vector_alloc(1);
			gsl_vector_set(b_ka,0, 1.01);
		 	ka	= apop_zipf_log_likelihood(b_ka , d);
			gsl_vector_free (b_ka);
		}
		return keep_away(bb, 1, 'b', ka);
	}			//else:
	for(j=0; j< data->size2; j++){
		for(i=0; i< data->size1; i++){
			z	 = -log(gsl_sf_zeta(bb)) - bb * log(j+1);
			like	+= gsl_matrix_get(data,i,j) * z;
		}
	}
//printf("z: %g %g", bb, like);
	return -like;
}	

void apop_zipf_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
double		a	= gsl_vector_get(beta, 0);
static double	ka	= 0;
gsl_matrix	*data	= d;
int 		i, j;
double		dlike	= 0, 
		colsum, dz;
	if (a <= 1) {		//keep away
		if (ka ==0){
			gsl_vector 	*b_ka	= gsl_vector_alloc(1);
			gsl_vector 	*b_kg	= gsl_vector_alloc(1);
			gsl_vector_set(b_ka,0, 1.01);
		 	apop_zipf_dlog_likelihood(b_ka , d, b_kg);
			ka	= gsl_vector_get(b_kg, 0);
			gsl_vector_free (b_ka);
			gsl_vector_free (b_kg);
		}
		gsl_vector_set(gradient,0,-keep_away(a, 1, 'b', ka));
		return;
	}			//else:
	for(j=0; j< data->size2; j++){
		colsum	= 0;
		dz	 = a*gsl_sf_zeta(a+1)/gsl_sf_zeta(a) - log(j+1);
		for(i=0; i< data->size1; i++){
			colsum	+= gsl_matrix_get(data,i,j);
		}
		dlike	+= colsum * dz;
	}
printf("dz: %g %g", a, dlike);
	gsl_vector_set(gradient,0,-dlike);
}	


/** Draw from a Zipf distribution with parameter \f$ a \f$

Call this fn using \ref apop_zipf.rng().
  
For example:
\code
gsl_rng *       r;
gsl_rng_env_setup();
r=gsl_rng_alloc(gsl_rng_taus);	//for example. 
apop_zipf.rng(r, 1.4);
\endcode

Cribbed from <a href="http://cgm.cs.mcgill.ca/~luc/mbookindex.html>Devroye (1986)</a>, p 551.  */
double apop_zipf_rng(gsl_rng* r, double a){
if (a  <= 1)	
	{printf("apop_zipf.rng: Zipf needs a parameter >=1. Returning 0.\n"); return 0;};
int		x;
double		u, v, t, 
		b 	= pow(2, a-1), 
		ainv	= -(1.0/(a-1));
	do {
		u	= gsl_rng_uniform(r);
		v	= gsl_rng_uniform(r);
		x	= pow(u, ainv);
		t	= pow((1.0 + 1.0/x), (a-1));
	} while (v * x * (t-1.0)/(b-1) > t/b);
	return x;
}

/** The exponential distribution. A one-parameter likelihood fn.

This one is just wrong. Ignore it for now.
 
\f$P(link ct==k) = (1/C) C^{-k}   		\f$. <br>
\f$Z(C,k) 	= ln(C) C^{-k} 			\f$ <br>
\f$ln Z(C,k) 	= -ln(ln(C)) * ln(C) * k	\f$ <br>
\f$dlnZ/dC	= -ln(ln(C)) * ln(C) * k	\f$ <br>
\f$d^2lnZ/dC^2	= k / C^2	\f$ <br>

\ingroup network_likelihoods
\todo Set up an exponential object which makes use of the GSL.
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

/** The Zipf distribution.

\f$Z(a)		= {1\over \zeta(a) * i^a}		\f$<br>
\f$lnZ(a)	= -(\log(\zeta(a)) + a \log(i))	\f$<br>
\f$dlnZ(a)/da	= -\zeta(a)\over a \log(\zeta(a-1))} -  \log(i)		\f$<br>
*/
apop_distribution apop_zipf = {1, apop_zipf_log_likelihood, apop_zipf_dlog_likelihood, NULL, apop_zipf_rng};

/** The Yule distribution

Yule likelihood fn. The special case of Waring where \f$ \alpha = 0.	\f$<br>

\f$ Y(x, b) 	= (b-1) \gamma(b) \gamma(k) / \gamma(k+b)			\f$<br>
\f$ \ln Y(x, b)	= \ln(b-1) + ln\gamma(b) + \ln\gamma(k) - \ln\gamma(k+b)	\f$<br>
\f$ d\ln Y/db	= 1/(b-1)  + \psi(b) - \psi(k+b)				\f$<br>
*/
apop_distribution apop_yule = {1, apop_yule_log_likelihood, apop_yule_dlog_likelihood, NULL, apop_yule_rng};

/** The Waring distribution

\f$W(x, b,a) 	= (b-1) \gamma(b+a) \gamma(k+a) / [\gamma(a+1) \gamma(k+a+b)]\f$

\f$\ln W(x, b, a) = ln(b-1) + lng(b+a) + lng(k+a) - lng(a+1) - lng(k+a+b)\f$

\f$dlnW/db	= 1/(b-1)  + \psi(b+a) - \psi(k+a+b)\f$

\f$dlnW/da	= \psi(b+a) + \psi(k+a) - \psi(a+1) - \psi(k+a+b)\f$
*/
apop_distribution apop_waring = {2, apop_waring_log_likelihood, apop_waring_dlog_likelihood, NULL, NULL};

/** The Gamma distribution

\f$G(x, a, b) 	= 1/(\Gamma(a) b^a)  x^{a-1} e^{-x/b}\f$

\f$ln G(x, a, b)= -ln \Gamma(a) - a ln b + (a-1)ln(x) + -x/b\f$

\f$d ln G/ da	=  -\psi(a) - ln b + ln(x) \f$	(also, \f$d ln \gamma = \psi\f$)

\f$d ln G/ db	=  -a/b - x \f$
\ingroup network_likelihoods
*/
apop_distribution apop_gamma = {2, apop_gamma_log_likelihood, apop_gamma_dlog_likelihood, NULL, NULL};

apop_distribution apop_probit = {2, apop_probit_log_likelihood, apop_probit_dlog_likelihood, apop_probit_fdf, NULL};
