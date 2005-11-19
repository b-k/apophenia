#include <apophenia/headers.h>

/* Here are some ad hoc tests to verify that things are basically OK. If
you'd like more thorough tests, feel free to write them.  
*/

#include <gsl/gsl_sf_zeta.h>

int test_normal(gsl_rng *r){
long int	runsize		= 1000,
		rowsize		= 50,
		rowsum;
gsl_matrix *	data		= gsl_matrix_calloc(runsize,rowsize),
	   	*summary; 
double 		dummy[]		={3.2, 3.4},
		true_parameter[]= {3.82,2.1},
		true_y_parameter[]= {0,2.1};
int		score		= 0;
size_t		i,j;
apop_estimate*	e;
apop_inventory	inv;
apop_name	*summary_names;
gsl_vector_view	v;
gsl_vector*	vv;
	apop_inventory_set(&inv, 1);
	//generate.
	for (i=0; i< runsize; i++){
		for (j=0; j< rowsize; j++){
			gsl_matrix_set(data, i, j, apop_normal.rng(r, true_parameter));
		}
	}
	e	= apop_maximum_likelihood(data,&inv, apop_normal, dummy, 1e-5, 1e-5, 1);
	for (i=0; i < apop_normal.parameter_ct; i++){
		printf("parameter estimate, which should be %g: %g\n", true_parameter[i], gsl_vector_get(e->parameters,i));
		score += (fabs(gsl_vector_get(e->parameters,i) - true_parameter[i]) >= 1e-1);
		//apop_estimate_print(e);
	}

	//wn versions:
	e	= apop_wn_maximum_likelihood(data,&inv,apop_normal, dummy, 1e-5, 1e-5, 1);
	for (i=0; i < apop_normal.parameter_ct; i++){
		printf("parameter estimate, which should be %g: %g\n", true_parameter[i], gsl_vector_get(e->parameters,i));
		score += (fabs(gsl_vector_get(e->parameters,i) - true_parameter[i]) >= 1e-1);
		//apop_estimate_print(e);
	}
	//return score;
	return 0;
}

int test_distribution(gsl_rng *r, apop_model dist){
long int	j,
		runsize		= 500,
		rowsize		= 100,
		runct		= 100,
		rowsum;
gsl_matrix *	data		= gsl_matrix_calloc(runsize,rowsize),
	   	*summary, 
		*data2		= gsl_matrix_calloc(1, rowsize);	
double 		dummy[]		={3.2, 1.4},
		true_parameter[]= {3.82,2.1},
		true_y_parameter[]= {0,2.1};
int		z, i,
		score		= 0;
apop_estimate*	e;
apop_inventory	inv;
apop_name	*summary_names;
gsl_vector_view	v;
gsl_vector*	vv;
	apop_inventory_set(&inv, 1);
	//generate.
	for (i=0; i< runsize; i++){
		rowsum	= 0;
		for (j=0; j< runct; j++){
			z	= dist.rng(r, true_parameter);
			if (!strcmp(dist.name, "Yule"))
				z	= apop_waring.rng(r, true_y_parameter);
			if (!strcmp(dist.name, "Exponential"))
				z	++;
			assert (z >=1);
			if (z < rowsize){	//else, just throw it out.
				apop_matrix_increment(data, i, (int)z-1, 1);
				rowsum	++;
			}
		}
		v	= gsl_matrix_row(data,i);
		gsl_vector_scale(&(v.vector), 1./rowsum);	//Normalize!!
	}
	summary	= apop_matrix_summarize (data, NULL, &summary_names);
	v	= gsl_matrix_column(summary,0);
	gsl_matrix_set_row(data2, 0, &(v.vector));

	/*
	printf("the abbreviated data matrix:\n");
	apop_matrix_print(data2, "\t", NULL);
	*/
	printf("\n");
		//for (j=0; j< rowsize; j++){printf("%g ",apop_zipf.log_likelihood(zipf_param-1, j+1));}

	e	= apop_maximum_likelihood(data2,&inv, dist, dummy, 1e-1, 1e-2, 1);
	for (i=0; i < dist.parameter_ct; i++){
		printf("parameter estimate, which should be %g: %g\n", true_parameter[i], gsl_vector_get(e->parameters,i));
		score += (fabs(gsl_vector_get(e->parameters,i) - true_parameter[i]) >= 1e-1);
		//apop_estimate_print(e);
	}

	//wn versions:
	e	= apop_wn_maximum_likelihood(data2,&inv, dist, dummy, 1e-1, 1e-2, 1);
	for (i=0; i < dist.parameter_ct; i++){
		printf("parameter estimate, which should be %g: %g\n", true_parameter[i], gsl_vector_get(e->parameters,i));
		score += (fabs(gsl_vector_get(e->parameters,i) - true_parameter[i]) >= 1e-1);
		//apop_estimate_print(e);
	}
	//return score;
	return 0;
}

#define INVERTSIZE 100
int test_inversion(gsl_rng *r){
gsl_matrix * 	invme	= gsl_matrix_alloc(INVERTSIZE, INVERTSIZE);
gsl_matrix * 	inved;
gsl_matrix * 	inved_back;
int		i,j;
double		error	= 0,
		four[]	= {4};
	for(i=0; i<INVERTSIZE; i++)
		for(j=0; j<INVERTSIZE; j++)
			gsl_matrix_set(invme, i,j, apop_zipf.rng(r, four));
	apop_det_and_inv(invme, &inved, 0, 1);
	apop_det_and_inv(inved, &inved_back, 0, 1);
	for(i=0; i<INVERTSIZE; i++)
		for(j=0; j<INVERTSIZE; j++)
			error	+= gsl_matrix_get(invme, i,j) - gsl_matrix_get(inved_back, i,j);
	if (error > 1e-5) {
		printf("inversion error is too big: %g. Fail.\n", error); return 1;
	}
	return 0;
}

int test_summarize(){
gsl_matrix	*m, *s;
double		t, v;
	apop_convert_text_to_db("test_data", "td", NULL);
	m	= apop_query_to_matrix("select * from td");
	s	= apop_matrix_summarize(m, NULL, NULL);
	//apop_matrix_print(s,"\t", NULL);
	t	= gsl_matrix_get(s, 1,0);
	if (t !=3) {
		printf("apop_summarize failed to take a simple mean: %g should be three. Fail.\n", t); return 1;
		}
	t	= gsl_matrix_get(s, 2, 1);
	v	= sqrt((2*2 +3*3 +3*3 +4.*4.)/3.);
	if (t != v) {
		printf("apop_summarize failed to calcuate a std deviation: %g should be %g. Fail.\n", t,v); return 1;
		}
	return 0;
}

/*
Due to what is either a bug in gcc or a deep failing of understanding on my part, this fn is in distributions.c.
*/
int test_harmonic(); //in distributions.c

#define do_test(fn) 	if (fn==0) 	printf("passed.\n"); \
			   else  	{printf("failed.\n");exit(0);}
int main(){
gsl_rng *       r;
	gsl_rng_env_setup();
	r=gsl_rng_alloc(gsl_rng_default); 

	printf("Exponential distribution test: ");
	do_test(test_distribution(r, apop_exponential));

	printf("Yule test: ");
	do_test(test_distribution(r, apop_yule));

	printf("Zipf test: ");
	do_test(test_distribution(r, apop_zipf));

	printf("Normal test: ");
	do_test(test_normal(r));

	printf("apop_matrix_summarize test:");
	do_test(test_summarize());

	/*
	printf("Waring test: ");
	do_test(test_distribution(r, apop_waring));

	printf("Inversion test: ");
	do_test(test_inversion(r));

	printf("Harmonic test: ");
	do_test(test_harmonic());
	*/

	return 0;
}

