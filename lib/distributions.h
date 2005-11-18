#ifndef __apop_likelihoods_h__
#define __apop_likelihoods_h__

/** This is an object to describe a likelihood function. It would primarily be
used for maximum likelihood estimation, but is intended to have anything
else you would want a probability distribution to have too, like a random number
generator.  */
typedef struct apop_likelihood{
	char	name[100];
	int	parameter_ct;
	
	/** the likelihood fn given data*/
	double (*log_likelihood)(const gsl_vector *beta, void *d);

	/** the derivative of the likelihood fn */
	void (*dlog_likelihood)(const gsl_vector *beta, void *d, gsl_vector *gradient);

	/** Do both of the above at once. Can be NULL if it'd just call them separately. */
	void (*fdf)( const gsl_vector *beta, void *d, double *f, gsl_vector *df);

	/** The constraint count */
	int constraint_ct;

	/** The vector of constraints */
	double (* constraint) (gsl_vector *beta, void *d, gsl_vector *inside_constraint);
	
	/** a random number generator. */
	double (*rng)(gsl_rng* r, double *a);
	//to add: 
	//the 2nd derivative of the likelihood fn
	//the MLE
} apop_likelihood;


extern apop_likelihood apop_exponential;
extern apop_likelihood apop_gamma;
extern apop_likelihood apop_gamma_rank;
extern apop_likelihood apop_gaussian;//synonym for apop_normal
extern apop_likelihood apop_normal;
extern apop_likelihood apop_probit;
extern apop_likelihood apop_waring;
extern apop_likelihood apop_yule;
extern apop_likelihood apop_zipf;



//This is in asst.c.
double apop_generalized_harmonic(int N, double s);

#endif
