#ifndef __apop_distributions_h__
#define __apop_distributions_h__

typedef struct distribution{
	int	parameter_ct;
	double (*log_likelihood)(const gsl_vector *beta, void *d);	//the likelihood fn given data
	void (*dlog_likelihood)(const gsl_vector *beta, void *d, gsl_vector *gradient);//the derivative of the likelihood fn
	void (*fdf)( const gsl_vector *beta, void *d, double *f, gsl_vector *df);// do both of the above at once.
	//the 2nd derivative of the likelihood fn
	//the MLE
	double (*rng)(gsl_rng* r, double a);//a random number generator.
} apop_distribution;


apop_distribution apop_zipf;
apop_distribution apop_yule;
apop_distribution apop_waring;
apop_distribution apop_gamma;



//This is in asst.c.
double apop_generalized_harmonic(int N, double s);

#endif
