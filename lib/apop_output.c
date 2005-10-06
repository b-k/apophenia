/** \file apop_output.c	 Some printing and gnuplot interface functions. 

No, gnuplot is not associated with the GNU, and it is not as free as GNU
software (they don't want forking), but it's free enough for virtually
any purposes. It interfaces with text files which are often painfully
ornery, so autogeneration of these files is very desirable, and that's
where this file comes in. It includes a few simple functions to produce
files which gnuplot can plot directly.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL.
*/

#include <apophenia/output.h>
/** Prep for gnuplot one of those cute scatterplots with a regression line through it.

Currently, you only get two dimensions.

\param	data	This is a copy of what you'd sent to the regression fn. That is, the first column is the dependent variable and the second is the independent. That is, what will be on Y axis is the <i>first</i> column, and what is on the X axis is the second. Custom for regressions and custom for graphs just clash on this one.
\param	est	The \ref apop_estimate structure your regression function gave you.
\param	n	The \ref apop_name structure, if any. If none, send NULL
\param outfile	The file to write to. It appends instead of overwriting, so you can prep the file if you want; see sample code.

The sample program below will pull data from a database (you'll need to
modify it to produce your own two-column table), then runs OLS, produces
a gnuplot file to write to a file named "scatter.eps", and finally calls
gnuplot for you. You'll still have the gnuplot file ("auto") if you want
to make further modifications.

\code
int main(){
gsl_matrix 	*data, *data_copy;
apop_estimate   *est;
apop_name       *n;
FILE		*f;
char		outfile[]	= "auto",
		do_me[10000];

	apop_open_db("cpia.db");
	data      =apop_query_to_matrix("select gnppercap, cpia from cpiagnp;");
	n          = apop_db_get_names();
	apop_close_db(0);

	//The regression destroys your data, so copy it first.
	data_copy	= gsl_matrix_alloc(data->size1, data->size2);
	gsl_matrix_memcpy(data_copy, data);

	//Run OLS, display results on terminal
	est  = apop_OLS(data, n, NULL);
	apop_estimate_print(est);

	//Prep the file with a header, then call the function.
	f    = fopen(outfile, "w");
	fprintf(f,"set term postscript;\n set output \"scatter.eps\"\n set yrange [0:*]\n");
	fclose(f);
	apop_plot_line_and_scatter(data_copy,est, n, outfile);

	//Have the machine run gnuplot for you. 
	sprintf(do_me, "gnuplot -persist %s", outfile);
	system (do_me);
	return 0;
}
\endcode
\todo The sample data here should correspond to that which apophenia ships with.
\ingroup output
*/
void apop_plot_line_and_scatter(gsl_matrix *data, apop_estimate *est, apop_name *n, char *outfile){
FILE *          f;
        if (outfile == NULL) 	f       = stdout;
        else    		f       = fopen(outfile, "a");
	fprintf(f, "f(x) = %g  + %g * x\n", gsl_vector_get(est->parameters,0), gsl_vector_get(est->parameters,1));
	if (n){
		fprintf(f, "set xlabel \"%s\"\n", n->colnames[1]);
		fprintf(f, "set ylabel \"%s\"\n", n->depnames[0]);
	}
	fprintf(f, "set key off\n");
	fprintf(f, "plot \"-\" using 2:1 , f(x) with lines;\n");
	if (outfile !=NULL)    fclose(f);
	apop_print_matrix(data, ", ", outfile);
}


/** This is a dumb little function to call gnuplot for you,
   in case you're so exceptionally lazy that you can't call
   <tt>apop_print_matrix(data, "\t", "outfile")</tt> yourself.

   I have such disdain for this function that it will probably be replaced shortly with something which works more like 
   \ref apop_plot_line_and_scatter, producing a text file which you then get to plot yourself.

\param data the data to be plotted.
\param plot_type 's'=surface plot; anything else = 2D x-y plot
\param delay the amount of time before gnuplot closes itself.
\ingroup output
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



/** Print a summary of each column of a table to the screen (i.e., STDOUT). 

\todo At the moment, only gives the mean and the standard deviation
of the data in each column; should give more in the near future.

\param data
The table to be summarized.

\param names
The \ref apop_name structure associated with the table. If there is no such structure, use <tt>NULL</tt>.
\ingroup output
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
*/

/** Put summary information about the columns of a table (mean, var) in a table.


\param data
The table to be summarized.

\param names_in The \ref apop_name structure associated with the table. If there is no such structure, use <tt>NULL</tt>.
\param names_out The \ref apop_name structure which will be associated with the output table.
\ingroup output
\todo At the moment, only gives the mean and the standard deviation
of the data in each column; should give more in the near future.
\todo We should probably let this summarize rows as well.
*/
gsl_matrix * apop_matrix_summarize(gsl_matrix *data, apop_name *names_in, apop_name **names_out){
int		i;
gsl_vector_view	v;
gsl_matrix	*out	= gsl_matrix_alloc(data->size2, 2);
double		mean, stddev;
char		rowname[10000]; //crashes on more than 10^9995 columns.
	if (names_out !=NULL){
		*names_out	= apop_name_alloc();
		apop_name_add(*names_out, "mean", 'c');
		apop_name_add(*names_out, "std dev", 'c');
		if (names_in !=NULL)
			for (i=0; i< data->size2; i++)
				apop_name_add(*names_out, names_in->colnames[i], 'r');
		else
			for (i=0; i< data->size2; i++){
				sprintf(rowname, "col %i", i);
				apop_name_add(*names_out, rowname, 'r');
			}
	}
	for (i=0; i< data->size2; i++){
                v       = gsl_matrix_column(data, i);
		mean	= apop_mean(&(v.vector));
		stddev	= sqrt(apop_var_m(&(v.vector),mean));
		gsl_matrix_set(out, i, 0, mean);
		gsl_matrix_set(out, i, 1, stddev);
	}	
	return out;
}
