#include <likelihoods.h>

/*
double zipf_likelihood(const gsl_vector *beta, void *d){
// P(link ct==k) = C^{-k}. So the PDF is [ln(C) * C^{-k}], which integrates to one;
// The log likelihood that a draw has degree k is ln(ln(C)) - [ln(C) * k],
// I don't need the first part in the search for the best C.
// A one-parameter likelihood fn.
float		bb	= gsl_vector_get(beta, 0),
		actual, total_error=0;
	if (bb <=1) return GSL_POSINF;	//a sign to the minimizer to look elsewhere.
int 		i, k;
int* 		classroom_matrix= d;
float 		likelihood 	= 0,
		l,
		ln_c		= log(bb),
		ln_ln_c		= log(ln_c);
	for (i=1; i< class_ct; i++)
		for (k=0; k< max_popularity; k++){
			likelihood	+= 
			l		 = ((float) CM(k, i)) * (ln_ln_c - ln_c * k);
			if (find_errors){		//a global var.
				actual		 = CM(k, i)/classroom_totals[i];
				total_error	+= gsl_pow_2(actual -  l);
			}
		}
	if (!find_errors)	return -likelihood;
	else			return total_error;
}


void d_zipf_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
float		bb		= gsl_vector_get(beta, 0);
int 		i, k;
int* 		classroom_matrix= d;
float 		d_likelihood 	= 0,
	ln_c		= log(bb);
	for (i=1; i< class_ct; i++)
		for (k=0; k< max_popularity; k++){
			d_likelihood	+= ((float) CM(k, i)) * (1/ln_c - k)/bb;
		}
	gsl_vector_set(gradient,0, -d_likelihood);
}

void zipf_fdf(const gsl_vector *beta, void *d, double *f, gsl_vector *df){
	*f	= zipf_likelihood(beta, d);
	d_zipf_likelihood(beta, d, df);
}



//Yule likelihood fn

double yule_likelihood(const gsl_vector *beta, void *d){
float		bb		= gsl_vector_get(beta, 0),
		actual, total_error=0;
	if (bb <=2) return GSL_POSINF;	//a sign to the minimizer to look elsewhere.
int 		i, k;
int* 		classroom_matrix= d;
float 		ln_k, ln_bb_k, el,
	likelihood 	= 0,
	ln_bb		= gsl_sf_lngamma(bb),
	ln_bb_less_1	= log(bb-1);
	for (k=1; k< max_popularity; k++)
		for (i=0; i< class_ct; i++){
			//if (k>1) 	ln_k	= gsl_sf_lngamma(k+1);
			if (k>1) 	ln_k	= gsl_sf_lngamma(k);
			else		ln_k	= 0;
			//ln_bb_k	= gsl_sf_lngamma(k+1+bb);
			ln_bb_k	= gsl_sf_lngamma(k+bb);
			likelihood	+= 
			el		 = ((float)CM(k, i)) *(ln_bb_less_1 + ln_k + ln_bb - ln_bb_k);
			if (find_errors){		//a global var.
				actual		 = ((float)CM(k, i))/classroom_totals[i];
				total_error	+= gsl_pow_2(actual - el);
			}
		}
	if (!find_errors)	return -likelihood;
	else			return total_error;
}

*/


//The probit section. It includes some trickery to avoid recalculating beta dot x.
//
//
gsl_vector *beta_dot_x ;
int	beta_dot_x_is_current	= 0;

void	dot(const gsl_vector *beta, gsl_matrix *data){
gsl_matrix_view p      		= gsl_matrix_submatrix(data,0,1,data->size1,data->size2-1);
		beta_dot_x 	= gsl_vector_alloc(data->size1);

        gsl_blas_dgemv (CblasNoTrans, 1.0, &p.matrix, beta, 0.0, beta_dot_x);//dot product
}

double probit_likelihood(const gsl_vector *beta, void *d){
	//find (data dot beta'), then find the integral of the Normal (0,1)
	//up to that point. Multiply likelihood either by that or by 1-that, depending 
	//on the choice the data made.
int		i;
long double	n, total_prob	= 0;
gsl_matrix *data 		= (gsl_matrix *) d;		//just type casting.

	dot(beta,data);
	for(i=0;i< data->size1; i++){
		n	=gsl_cdf_gaussian_P(gsl_vector_get(beta_dot_x,i),1);
		if (gsl_matrix_get(data, i, 0)==0) 	total_prob	+= log(n);
		else 					total_prob	+= log((1 - n));
	}
	return -total_prob;
}

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
	*f=probit_likelihood(beta, d);
	beta_dot_x_is_current	=1;
	d_probit_likelihood(beta, d, df);
	beta_dot_x_is_current	=0;
}




//The max likelihood functions themselves. Mostly straight out of the GSL manual.



double	maximum_likelihood(void * data, gsl_vector **betas, int betasize,
					double (* likelihood)(const gsl_vector *beta, void *d), int verbose){
int			iter =0, status;
gsl_multimin_function 	minme;
gsl_multimin_fminimizer *s;
gsl_vector 		*x, *ss;
double			size;
	s	= gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, betasize);
	x	= gsl_vector_alloc(betasize);
	*betas	= gsl_vector_alloc(betasize);
	ss	= gsl_vector_alloc(betasize);
  //	gsl_vector_set_all (x, 3);
  	gsl_vector_set_all (ss, .01);

	minme.f		= likelihood;
	minme.n		= betasize;
	minme.params	= data;
	gsl_multimin_fminimizer_set (s, &minme, x,  ss);

      	do { 	iter++;
		status 	= gsl_multimin_fminimizer_iterate(s);
		if (status) 	break; 
		size	= gsl_multimin_fminimizer_size(s);
        	status 	= gsl_multimin_test_size (size, 1e-3); 
		if(verbose){
		int i;
			printf ("%5d ", iter);
			for (i = 0; i < betasize; i++) {
				printf ("%8.3e ", gsl_vector_get (s->x, i)); } 
			printf ("f()=%7.3f size=%.3f\n", s->fval, size);
       			if (status == GSL_SUCCESS) {
	   			printf ("Minimum found at:\n");
				printf ("%5d ", iter);
				for (i = 0; i < betasize; i++) {
					printf ("%8.3e ", gsl_vector_get (s->x, i)); } 
				printf ("f()=%7.3f size=%.3f\n", s->fval, size);
			}
		}
      	} while (status == GSL_CONTINUE && iter < MAX_ITERATIONS);
	if (iter == MAX_ITERATIONS)
		printf("Minimization reached maximum number of iterations.");

	gsl_vector_memcpy(*betas, s->x);
	return 0; 
}



void	maximum_likelihood_w_d(void * data, gsl_vector **betas, int betasize,
					double (* likelihood)(const gsl_vector *beta, void *d),
					void (* d_likelihood)(const gsl_vector *beta, void *d, gsl_vector *df), 
					void (* fdf)(const gsl_vector *beta, void *d, double *f, gsl_vector *df), int verbose){


int			iter =0, status;
gsl_multimin_function_fdf 	minme;
gsl_multimin_fdfminimizer 	*s;
gsl_vector 		*x, *ss;
	s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, betasize);
	//s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, betasize);
	x	= gsl_vector_alloc(betasize);
	*betas	= gsl_vector_alloc(betasize);
	ss	= gsl_vector_alloc(betasize);
  	gsl_vector_set_all (x,  1.2);
  	gsl_vector_set_all (ss, 0.01);
	minme.f		= likelihood;
	minme.df	= d_likelihood;
	minme.fdf	= fdf;
	minme.n		= betasize;
	minme.params	= (void *)data;
	gsl_multimin_fdfminimizer_set (s, &minme, x,  .001, .0001);

      	do { 	iter++;
		status 	= gsl_multimin_fdfminimizer_iterate(s);
		if (status) 	break; 
		status = gsl_multimin_test_gradient(s->gradient, .0001);
		//size	= gsl_multimin_fminimizer_size(s);
        	//status 	= gsl_multimin_test_size (size, 1e-5); 
        	if (verbose){
        		if (status == GSL_SUCCESS) 
		   		printf ("Minimum found.\n");
	        	printf ("%5i %.5f  %10.5f\n", iter,
	        		gsl_vector_get (s->x, 0),  s->f);
		}
       	 }
       	while (status == GSL_CONTINUE && iter < MAX_ITERATIONS_w_d);
	if(iter==MAX_ITERATIONS_w_d) printf("No min!!\n");

	gsl_vector_memcpy(*betas, s->x);
	gsl_multimin_fdfminimizer_free(s);
}
