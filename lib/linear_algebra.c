//linear_algebra.c		  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
#include <gsl/gsl_blas.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>	//popen, I think.
#include "linear_algebra.h" 
#include "stats.h" 
#include "math.h" //pow!
#include "gnulib/vasprintf.h"

double apop_det_and_inv(gsl_matrix *in, gsl_matrix *out, int calc_det, int calc_inv) {
int 		sign;
double 		the_determinant = 0;
	gsl_matrix *invert_me = gsl_matrix_alloc(in->size1, in->size1);
	gsl_permutation * perm = gsl_permutation_alloc(in->size1);
	invert_me = gsl_matrix_alloc(in->size1, in->size1);
	gsl_matrix_memcpy (invert_me, in);
	gsl_linalg_LU_decomp(invert_me, perm, &sign);
	if (calc_inv)
		gsl_linalg_LU_invert(invert_me, perm, out);
	if (calc_det)
		the_determinant	= gsl_linalg_LU_det(invert_me, sign);
	gsl_matrix_free(invert_me);
	gsl_permutation_free(perm);
	return(the_determinant);
}

double apop_x_prime_sigma_x(gsl_vector *x, gsl_matrix *sigma){
//This comes up often enough that it deserves its own convenience function.
gsl_vector * 	sigma_dot_x	= gsl_vector_calloc(x->size);
double		the_result;
	gsl_blas_dsymv(CblasUpper, 1, sigma, x, 0, sigma_dot_x); //sigma should be symmetric
	gsl_blas_ddot(x, sigma_dot_x, &the_result);
	gsl_vector_free(sigma_dot_x);
	return(the_result);
}

void apop_normalize_for_svd(gsl_matrix *in){
//Greene (2nd ed, p 271) recommends pre- and post-multiplying by sqrt(diag(X'X)) so that X'X = I.
gsl_vector_view	v;
gsl_vector	*diagonal = gsl_vector_alloc(in->size1);
int 		i;
	//Get the diagonal, take the square root
	v	= gsl_matrix_diagonal(in);
	gsl_vector_memcpy(diagonal, &(v.vector));
	for (i=0; i<diagonal->size; i++)
		gsl_vector_set(diagonal, i, pow(gsl_vector_get(diagonal,i), .5));
	//mulitply each row and column by the diagonal vector.
	for (i=0; i<diagonal->size; i++){
		v	= gsl_matrix_column(in, i);
		gsl_vector_mul(&(v.vector), diagonal);
		v	= gsl_matrix_row(in, i);
		gsl_vector_mul(&(v.vector), diagonal);
	}
	gsl_vector_free(diagonal);
}

void apop_sv_decomposition(gsl_matrix *data, int dimensions_we_want, gsl_matrix ** pc_space, gsl_vector **total_explained) {
//Get X'X
gsl_matrix * 	eigenvectors 	= gsl_matrix_alloc(data->size2, data->size2);
gsl_vector * 	dummy_v 	= gsl_vector_alloc(data->size2);
gsl_vector * 	all_evalues 	= gsl_vector_alloc(data->size2);
gsl_matrix * 	square  	= gsl_matrix_calloc(data->size2, data->size2);
gsl_vector_view v;
int 		i;
double		eigentotals	= 0;
	*pc_space	= gsl_matrix_alloc(data->size2, dimensions_we_want);
	*total_explained= gsl_vector_alloc(dimensions_we_want);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, data, data, 0, square);
	apop_normalize_for_svd(square);	
	gsl_linalg_SV_decomp(square, eigenvectors, all_evalues, dummy_v);
	for (i=0; i< all_evalues->size; i++)
		eigentotals	+= gsl_vector_get(all_evalues, i);
	for (i=0; i<dimensions_we_want; i++){
		v	= gsl_matrix_column(eigenvectors, i);
		gsl_matrix_set_col(*pc_space, i, &(v.vector));
		gsl_vector_set(*total_explained, i, gsl_vector_get(all_evalues, i)/eigentotals);
	}
	gsl_vector_free(dummy_v); 	gsl_vector_free(all_evalues);
	gsl_matrix_free(square); 	gsl_matrix_free(eigenvectors);
}



inline void apop_vector_increment(gsl_vector * v, int i, double amt){
	v->data[i * v->stride]	+= amt;
}

inline void apop_matrix_increment(gsl_matrix * m, int i, int j, double amt){
	m->data[i * m->tda +j]	+= amt;
}


////////////////////////////
/////The printing functions.
////////////////////////////

void print_core_v(gsl_vector *data, char *separator, char *filename, 
			void (* p_fn)(FILE * f, double number)){
int 		i;
FILE * 		f;
	if (filename == NULL)
		f	= stdout;
	else	f	= fopen(filename, "a");
	for (i=0; i<data->size; i++){
		p_fn(f, gsl_vector_get(data, i));
		if (i< data->size -1)	fprintf(f, "%s", separator);
	}
	fprintf(f,"\n");
	if (filename !=NULL)	fclose(f);
}

void print_core_m(gsl_matrix *data, char *separator, char *filename, 
			void (* p_fn)(FILE * f, double number)){
FILE * 		f;
int 		i,j;
	if (filename == NULL)
		f	= stdout;
	else	f	= fopen(filename, "a");
	for (i=0; i<data->size1; i++){
		for (j=0; j<data->size2; j++){
			p_fn(f, gsl_matrix_get(data, i,j));
			if (j< data->size2 -1)	fprintf(f, "%s", separator);
		}
		fprintf(f,"\n");
	}
	if (filename !=NULL)	fclose(f);
}

void dumb_little_pf_f(FILE * f, double data){
	fprintf(f, "% 5f", data); }

void dumb_little_pf_i(FILE * f, double data){
	fprintf(f, "% 5i", (int) data); }

void apop_print_vector(gsl_vector *data, char *separator, char *filename){
	print_core_v(data, separator, filename, dumb_little_pf_f); }

void apop_print_vector_int(gsl_vector *data, char *separator, char *filename){
	print_core_v(data, separator, filename, dumb_little_pf_i); }

void apop_print_matrix(gsl_matrix *data, char *separator, char *filename){
	print_core_m(data, separator, filename, dumb_little_pf_f); }

void apop_print_matrix_int(gsl_matrix *data, char *separator, char *filename){
	print_core_m(data, separator, filename, dumb_little_pf_i); }


void apop_plot(gsl_matrix *data, char plot_type, int delay){
/* This is a dumb little function to call gnuplot for you,
   in case you're so exceptionally lazy that you can't call
   apop_print_matrix(data, "\t", "outfile") yourself.
   It's so silly, I don't even document it.
   plot_type: 's'=surface plot; anything else = 2D x-y plot
   delay: the amount of time before gnuplot closes itself.
*/
FILE 		*output;
int		i,j;
	output = popen ("gnuplot", "w");
	if (!output) {
		fprintf (stderr, "Can't find gnuplot.\n");
		return;
	}
  	if (plot_type == 's')
		fprintf(output, "splot \"-\"\n");
  	if (plot_type != 's')
		fprintf(output, "plot \"-\" using 1:2\n");
	for (i=0; i<data->size1; i++){
		for (j=0; j<data->size2; j++){
			fprintf(output, "%g", gsl_matrix_get(data, i,j));
			if (j< data->size2 -1)	fprintf(output, "\t");
		}
		fprintf(output,"\n");
	}
	fprintf(output,"e\n pause %i\n", delay);
	pclose (output);
}

void apop_matrix_summarize(gsl_matrix *data, apop_name *names){
int		i;
gsl_vector_view	v;
	if (names !=NULL)
		printf("names");
	printf("\tmean:\tstd dev:\n");
	for (i=0; i< data->size2; i++){
                v       = gsl_matrix_column(data, i);
		if (names !=NULL)
			printf("%s\t%5f\t%5f\n",names->colnames[i],apop_mean(&(v.vector)),sqrt(apop_var(&(v.vector))));
		else
			printf("col %i\t%5f\t%5f\n",i,apop_mean(&(v.vector)),sqrt(apop_var(&(v.vector))));
	}	
}
