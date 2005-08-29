/** \file linear_algebra.c	Assorted things to do with matrices,
such as take determinants or do singular value decompositions.



	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
*/

/** \defgroup linear_algebra 	Singular value decompositions, determinants, et cetera.  */

/** \defgroup convenience_fns 	Things to make life easier with the GSL.
 */

/** \defgroup apop_print 	Asst printing functions		

Many have multiple aliases, because I could never remember which way to write them.
*/

#include <gsl/gsl_blas.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>	//popen, I think.
#include <apophenia/linear_algebra.h> 
#include <apophenia/stats.h>
#include "math.h" //pow!
#include <apophenia/vasprintf.h>

gsl_matrix *apop_covariance_matrix(gsl_matrix *in, int normalize){
gsl_matrix	*out;
int		i,j,k;
double		means[in->size2];
gsl_vector_view	v;
	if (normalize){
		out	= gsl_matrix_alloc(in->size2, in->size2);
		apop_normalize_matrix(in);
		gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1, in, in, 0, out);
	}
	else{
		out	= gsl_matrix_calloc(in->size2, in->size2);
		for(i=0; i< in->size2; i++){
			v		= gsl_matrix_column(in, i);
			means[i]	= apop_mean(&(v.vector));
		}
		for(i=0; i< in->size2; i++)
			for(j=i; j< in->size2; j++){
				for(k=0; k< in->size1; k++)
					apop_matrix_increment(out, i, j, 
					   (gsl_matrix_get(in,k,i)-means[i])* (gsl_matrix_get(in,k,j)-means[j]));
				if (i != j)	//set the symmetric element.
					gsl_matrix_set(out, j, i, gsl_matrix_get(out, i, j));
			}
	}
	gsl_matrix_scale(out, 1.0/in->size2);
	return out;
}

/**
Calculate the determinant of a matrix, its inverse, or both. The \c in matrix is not destroyed in the process.

\param in
The matrix to be inverted/determined.

\param out 
If you want an inverse, this is where to place the matrix to be filled with the inverse. If <tt>calc_inv == 0</tt>, then use <tt>NULL</tt>. 

\param calc_det 
0: Do not calculate the determinant.

1: Do.

\param calc_inv
0: Do not calculate the inverse.

1: Do.

\return
If <tt>calc_det == 1</tt>, then return the determinant. Otherwise, just returns zero.
\ingroup linear_algebra
*/
double apop_det_and_inv(gsl_matrix *in, gsl_matrix **out, int calc_det, int calc_inv) {
int 		sign;
double 		the_determinant = 0;
	gsl_matrix *invert_me = gsl_matrix_alloc(in->size1, in->size1);
	gsl_permutation * perm = gsl_permutation_alloc(in->size1);
	invert_me = gsl_matrix_alloc(in->size1, in->size1);
	gsl_matrix_memcpy (invert_me, in);
	gsl_linalg_LU_decomp(invert_me, perm, &sign);
	if (calc_inv){
		*out	= gsl_matrix_alloc(in->size1, in->size1); //square.
		gsl_linalg_LU_invert(invert_me, perm, *out);
		}
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

/**
Singular value decomposition, aka principal component analysis, aka factor analysis.

\param data 
The input matrix.

\param dimensions_we_want 
The singular value decomposition will return this many of the eigenvectors with the largest eigenvalues.

\param pc_space 
This will be the principal component space. Each column of the returned matrix will be another eigenvector; the columns will be ordered by the eigenvalues. Input the address of an un-allocated {{{gsl_matrix}}}.

\param total_explained
This will return the largest eigenvalues, scaled by the total of all eigenvalues (including those that were thrown out). The sum of these returned values will give you the percentage of variance explained by the factor analysis.
\ingroup linear_algebra */
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


/** Just add <tt>amt</tt> to a \c gsl_vector element. Equivalent to <tt>gsl_vector_set(gsl_vector_get(v, i) + amt, i)</tt>, but more readable (and potentially faster).

\param v The \c gsl_vector in question
\param i The location in the vector to be incremented.
\param amt The amount by which to increment. Of course, one can decrement by specifying a negative amt.
\ingroup convenience_fns
 */
inline void apop_vector_increment(gsl_vector * v, int i, double amt){
	v->data[i * v->stride]	+= amt;
}

/** Just add <tt>amt</tt> to a \c gsl_matrix element. Equivalent to <tt>gsl_matrix_set(gsl_matrix_get(m, i, j) + amt, i, j)</tt>, but more readable (and potentially faster).

\param m The \c gsl_matrix in question
\param i The row of the element to be incremented.
\param j The column of the element to be incremented.
\param amt The amount by which to increment. Of course, one can decrement by specifying a negative amt.
\ingroup convenience_fns
 */
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

/** Print a vector in real format.

\ingroup apop_print */
void apop_print_vector(gsl_vector *data, char *separator, char *filename){
	print_core_v(data, separator, filename, dumb_little_pf_f); }

/** Print a vector in int format.

\ingroup apop_print */
void apop_print_vector_int(gsl_vector *data, char *separator, char *filename){
	print_core_v(data, separator, filename, dumb_little_pf_i); }

/** Print a matrix in real format.

\ingroup apop_print */
void apop_print_matrix(gsl_matrix *data, char *separator, char *filename){
	print_core_m(data, separator, filename, dumb_little_pf_f); }

/** Print a matrix in int format.

\ingroup apop_print */
void apop_print_matrix_int(gsl_matrix *data, char *separator, char *filename){
	print_core_m(data, separator, filename, dumb_little_pf_i); }

/** Print a vector in float format.

\ingroup apop_print */
void apop_vector_print(gsl_vector *data, char *separator, char *filename){
	print_core_v(data, separator, filename, dumb_little_pf_f); }

/** Print a vector in int format.

\ingroup apop_print */
void apop_vector_print_int(gsl_vector *data, char *separator, char *filename){
	print_core_v(data, separator, filename, dumb_little_pf_i); }

/** Print a matrix in float format.

\ingroup apop_print */
void apop_matrix_print(gsl_matrix *data, char *separator, char *filename){
	print_core_m(data, separator, filename, dumb_little_pf_f); }

/** Print a matrix in int format.

\ingroup apop_print */
void apop_matrix_print_int(gsl_matrix *data, char *separator, char *filename){
	print_core_m(data, separator, filename, dumb_little_pf_i); }


/** This is a dumb little function to call gnuplot for you,
   in case you're so exceptionally lazy that you can't call
   <tt>apop_print_matrix(data, "\t", "outfile")</tt> yourself.

\param data the data to be plotted.
\param plot_type 's'=surface plot; anything else = 2D x-y plot
\param delay the amount of time before gnuplot closes itself.
*/
void apop_plot(gsl_matrix *data, char plot_type, int delay){
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

/** Print a summary of each column of a table to the screen (i.e., STDOUT). 

\todo At the moment, only gives the mean and the standard deviation
of the data in each column; should give more in the near future.

\param data
The table to be summarized.

\param names
The \c apop_names structure associated with the table. If there is no such structure, use <tt>NULL</tt>.
*/
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
