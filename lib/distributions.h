#ifndef __apop_distributions_h__
#define __apop_distributions_h__

/** This is an object to describe a distribution. It would primarily be
used for maximum likelihood estimation, but is intended to have anything
else you would want a distribution to have too, like a random number
generator. */
typedef struct distribution{
	char	name[100];
	int	parameter_ct;
	/** the likelihood fn given data*/
	double (*log_likelihood)(const gsl_vector *beta, void *d);
	/** the derivative of the likelihood fn */
	void (*dlog_likelihood)(const gsl_vector *beta, void *d, gsl_vector *gradient);
	/** Do both of the above at once. Can be NULL if it'd just call them separately. */
	void (*fdf)( const gsl_vector *beta, void *d, double *f, gsl_vector *df);
	/** a random number generator. */
	double (*rng)(gsl_rng* r, double *a);
	//to add: 
	//the 2nd derivative of the likelihood fn
	//the MLE
} apop_distribution;


apop_distribution apop_exponential;
apop_distribution apop_gamma;
apop_distribution apop_probit;
apop_distribution apop_waring;
apop_distribution apop_yule;
apop_distribution apop_zipf;



//This is in asst.c.
double apop_generalized_harmonic(int N, double s);

#endif
