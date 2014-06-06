
#define Apop_submatrix APOP_SUBMATRIX
/** \file 
  Some printing and output interface functions. */
/* Copyright (c) 2006--2007, 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

//The reader will find a few function headers for this file in asst.h
#include "apop_internal.h"

/** \defgroup output		Printing to the screen or a text file

Most functions print only to the screen, but the 
\ref apop_print "matrix and vector printing functions" will let you print to a text file as
well. The presumption is that statistic estimates are for your own
consumption, while you are printing a matrix for import into another program.

*/
/** \defgroup apop_print 	Assorted printing functions		

The <tt>apop_*_print</tt> functions will print to screen, text file,
or database, depending on how you set \c .output_type.
The <tt>apop_*_show</tt> functions print only to screen, and are basically
just a convenience shell to the corresponding <tt>apop_*_print</tt>
function.

\ingroup output
*/
 
#define Output_vars output_name, output_pipe, output_type, output_append

#define Output_declares char const * output_name, FILE * output_pipe, char output_type, char output_append

/** If you're reading this, it is probably because you were referred by another function
  that uses this internally. You should never call this function directly, but do read
  this documentation.

  There are four settings that affect how output happens, which can be set when you call the
  function that sent you to this documentation, e.g:

  \code
  apop_data_print(your_data, .output_type ='f', .output_append = 'w');
  \endcode

  \param output_name The name of the output file, if any. For a database, the table to write.
  \param output_pipe If you have already opened a file and have a \c FILE* on hand, use
  this instead of giving the file name.
  \param output_type \c 'p' = pipe, \c 'f'= file, \c 'd' = database, \c 's' = stdout
  \param output_append \c 'a' = append (default), \c 'w' = write over.

At the end, \c output_name, \c output_pipe, and \c output_type are all set.
Notably, the local \c output_pipe will have the correct location for the calling function to \c fprintf to.
*/
int apop_prep_output(char const *output_name, FILE ** output_pipe, char *output_type, char *output_append){
    *output_append = *output_append ? *output_append : 'w';

    if (!output_name && !*output_pipe && !*output_type)     *output_type = 's';              
    else if (output_name && !*output_pipe && !*output_type) *output_type = 'f'; 
    else if (!output_name && *output_pipe && !*output_type) *output_type = 'p';     

    if (*output_type =='p')      *output_pipe = *output_pipe ? *output_pipe: stdout;      
    else if (*output_type =='s') *output_pipe = stdout; 
    else if (*output_type =='d') *output_pipe = stdout;  //won't be used.
    else *output_pipe = output_name
                        ? fopen(output_name, *output_append == 'a' ? "a" : "w")
                        : stdout;
    Apop_stopif(!output_pipe && output_name, return -1, 0, "Trouble opening file %s.", output_name);
    return 0;
}

#define Dispatch_output                        \
    char const *apop_varad_var(output_name, NULL);  \
    FILE * apop_varad_var(output_pipe, NULL);  \
    char apop_varad_var(output_type, 0);       \
    char apop_varad_var(output_append, 0);     \
    Apop_stopif(apop_prep_output(output_name, &output_pipe, &output_type, &output_append), \
            return, 0, "Trouble preparing to write output.");

/** Prep for Gnuplot one of those cute scatterplots with a regression line through it.

Currently, you only get two dimensions.

Set the global \ref apop_opts_type "apop_opts.output_name" to the filename you want before running this.
It appends instead of overwriting, so you can prep the file if you want; see sample code. [to overwrite a file, just remove it first with the standard C function <tt>remove("filename");</tt>]

\param	data	This is a copy of what you'd sent to the regression fn. That is, the first column is the dependent variable and the second is the independent. That is, what will be on Y axis is the <i>first</i> column, and what is on the X axis is the second. Custom for regressions and custom for graphs just clash on this one.
\param	est	The \ref apop_model structure your regression function gave
you. (if \c NULL, I'll estimate an OLS model for you).

The sample program below will pull data from a database (ridership at
the Silver Spring, MD Metro station; get the database in the {\em Modeling
with Data} sample code, at http://modelingwithdata.org/appendices.html), then runs OLS, and produce
a Gnuplot file to write to a file named "scatter.eps". You can run the
result through Gnuplot via <tt> gnuplot scatter.gplot</tt>, and if you
don't like the results,  you have the Gnuplot file ("scatter.gplot") on
hand for modifications.

\include scatter.c

\li See \ref apop_prep_output for more on how settings are set.
\li This function uses the \ref designated syntax for inputs.
\ingroup output
*/
#ifdef APOP_NO_VARIADIC
void apop_plot_line_and_scatter(apop_data *data, apop_model *est, Output_declares){
#else
apop_varad_head(void, apop_plot_line_and_scatter){
    int newmodel=0;
    apop_data * apop_varad_var(data, NULL);
    apop_model * apop_varad_var(est, NULL);
    if (!est) {
        newmodel++;
        est = apop_estimate(data, apop_ols);
    }
    Dispatch_output
    apop_plot_line_and_scatter_base(data, est, Output_vars);
    apop_model_free(est);
    return;
     apop_plot_line_and_scatter_base(data, est, Output_vars);
}

 void apop_plot_line_and_scatter_base(apop_data *data, apop_model *est, Output_declares){
#endif
    char  exdelimiter[100];
    FILE *f = output_pipe;
	fprintf(f, "f(x) = %g  + %g * x\n", gsl_vector_get(est->parameters->vector,0), gsl_vector_get(est->parameters->vector,1));
	if (data->names){
		fprintf(f, "set xlabel \"%s\"\n", data->names->col[1]);
		fprintf(f, "set xlabel \"%s\"\n", data->names->vector);
	}
	fprintf(f, "plot \"-\" using 2:1 , f(x) with lines;\n");
    if (output_type == 'f') 	fclose(f);

    //force the delimiter to be a comma space; don't tell the user.
    strcpy(exdelimiter, apop_opts.output_delimiter);
    strcpy(apop_opts.output_delimiter, ", ");
	apop_matrix_print(data->matrix, output_name, output_pipe, output_type, 'a');
    strcpy(apop_opts.output_delimiter, exdelimiter);
}

/** This convenience function will take in a \c gsl_vector of data and put out a histogram, ready to pipe to Gnuplot.

The function respects the <tt>output_type</tt> option, so code like:
\code
FILE* f = popen("/usr/bin/gnuplot", "w");
apop_plot_histogram(data, .bin_count=100, .output_pipe = f);
\endcode
will print directly to Gnuplot.

\param data A \c gsl_vector holding the data. Do not pre-sort or bin; this function does that for you. (no default, must not be \c NULL)
\param bin_count   The number of bins in the output histogram (default = \f$\sqrt(N)\f$, where \f$N\f$ is the length of the vector.)
\param with The method for Gnuplot's plotting routine. Default is \c "boxes", so the gnuplot call will read <tt>plot '-' with boxes</tt>. The \c "lines" option is also popular, and you can add extra terms if desired, like <tt> "boxes linetype 3"</tt>.

\li See \ref apop_prep_output for more on how printing settings are set.
\li See also the legible output section of the \ref outline for more details and examples.
\li This function uses the \ref designated syntax for inputs.
  \ingroup output
*/
#ifdef APOP_NO_VARIADIC
void apop_plot_histogram(gsl_vector *data, size_t bin_count, char *with, Output_declares){
#else
apop_varad_head(void, apop_plot_histogram){
    gsl_vector * apop_varad_var(data, NULL);
    Apop_assert_n(data, "Input vector is NULL.");
    size_t apop_varad_var(bin_count, 0);
    char * apop_varad_var(with, "impulses");
    Dispatch_output
     apop_plot_histogram_base(data, bin_count, with, Output_vars);
}

 void apop_plot_histogram_base(gsl_vector *data, size_t bin_count, char *with, Output_declares){
#endif
    apop_data vector_as_data = (apop_data){.vector=data};
    apop_data *histodata = apop_data_to_bins(&vector_as_data, .bin_count=bin_count, .close_top_bin='y');
    apop_data_sort(histodata);
    apop_data_free(histodata->more); //the binspec.

    fprintf(output_pipe, "set key off	;\n"
               "plot '-' with %s\n", with);
    apop_data_print(histodata, .output_pipe=output_pipe);
    fprintf(output_pipe, "e\n");

    if (output_type == 'p') fflush(output_pipe);
    else if (output_name)   fclose(output_pipe);
    apop_data_free(histodata);
}

/////The printing functions.

static void white_pad(int ct){
    for(size_t i=0; i < ct; i ++)
        printf(" ");
}

/** This function prettyprints the \c apop_data set to a screen.

This takes a lot of machinery. I write every last element to a text array, then measure column widths, then print to screen with padding to guarantee that everything lines up.  There's no way to have the first element of a column line up with the last unless you interrogate the width of every element in the column, so printing columns really can't be a one-pass process.

So, I produce an \ref apop_data set with no numeric elements and a text element to be filled with the input data set, and then print that. That means that I'll be using (more than) twice the memory to print this. If this is a problem, you can use \ref apop_print to dump your data to a text file, and view the text file, or print subsets.

For more machine-readable printing, see \ref apop_data_print.

\ingroup output
*/
void apop_data_show(const apop_data *in){
    if (!in) {printf("NULL\n"); return;}
    Get_vmsizes(in) //vsize, msize1, msize2, tsize
//Take inventory and get sizes
    size_t hasrownames = (in->names && in->names->rowct) ? 1 : 0;
    size_t hascolnames = in->names && 
                    (in->names->vector || in->names->colct || in->names->textct);
    size_t hasweights = (in->weights != NULL);

    size_t outsize_r = GSL_MAX(in->matrix ? in->matrix->size1 : 0, in->vector ? in->vector->size: 0);
    outsize_r = GSL_MAX(outsize_r, in->textsize[0]);
    outsize_r = GSL_MAX(outsize_r, wsize);
    if (in->names) outsize_r = GSL_MAX(outsize_r, in->names->rowct);
    outsize_r += hascolnames;

    size_t outsize_c = msize2;
    outsize_c += in->textsize[1];
    outsize_c += (vsize>0);
    outsize_c += (wsize>0);
    outsize_c += hasrownames + hasweights;

//Write to the printout data set.
    apop_data *printout = apop_text_alloc(NULL , outsize_r, outsize_c);
    if (hasrownames)
        for (size_t i=0; i < in->names->rowct; i ++)
            apop_text_add(printout, i + hascolnames, 0, "%s", in->names->row[i]);
    for (size_t i=0; i < vsize; i ++) //vsize may be zero.
        apop_text_add(printout, i + hascolnames, hasrownames, "%g", gsl_vector_get(in->vector, i));
    for (size_t i=0; i < msize1; i ++) //msize1 may be zero.
        for (size_t j=0; j < msize2; j ++)
            apop_text_add(printout, i + hascolnames, hasrownames + (vsize >0)+ j, "%g", gsl_matrix_get(in->matrix, i, j));
    if (in->textsize[0])
        for (size_t i=0; i < in->textsize[0]; i ++)
            for (size_t j=0; j < in->textsize[1]; j ++)
                apop_text_add(printout, i + hascolnames, hasrownames + (vsize>0)+ msize2 + j, "%s", in->text[i][j]);
    if (hasweights)
        for (size_t i=0; i < in->weights->size; i ++)
            apop_text_add(printout, i + hascolnames, outsize_c-1, "%g", gsl_vector_get(in->weights, i));

//column names
    if (hascolnames){
        if (vsize && in->names->vector)
            apop_text_add(printout, 0 , hasrownames, "%s", in->names->vector);
        if (msize2 && in->names)
            for (size_t i=0; i < in->names->colct; i ++)
                apop_text_add(printout, 0 , hasrownames + (vsize>0) + i, "%s", in->names->col[i]);
        if (in->textsize[1] && in->names)
            for (size_t i=0; i < in->names->textct; i ++)
                apop_text_add(printout, 0 , hasrownames + (vsize>0) + msize2 + i, "%s", in->names->text[i]);
        if (hasweights)
            apop_text_add(printout, 0 , outsize_c-1, "Weights");
    }

//get column sizes
    int colsizes[outsize_c];
    for (size_t i=0; i < outsize_c; i ++){
        colsizes[i] = strlen(printout->text[0][i]);
        for (size_t j=1; j < outsize_r; j ++)
            colsizes[i] = GSL_MAX(colsizes[i], strlen(printout->text[j][i]));
    }

//Finally, print
    if (in->names && in->names->title && strlen(in->names->title))
        printf("\t%s\n\n", in->names->title);
    for (size_t j=0; j < outsize_r; j ++){
        for (size_t i=0; i < outsize_c; i ++){
            white_pad(colsizes[i] - strlen(printout->text[j][i]) + 1);//one spare space.
            printf("%s", printout->text[j][i]);
            if (i > 0 && i< outsize_c-1) 
                printf(" %s ", apop_opts.output_delimiter);
        }
        printf("\n");
    }

    if (in->more) {
        printf("\n");
        apop_data_show(in->more);
    }
    apop_data_free(printout);
}

void p_fn(FILE * f, double data){
    if (data == (int) data) fprintf(f, "% 5i", (int) data); 
    else                    fprintf(f, "% 5f", data);
}

static void print_core_v(const gsl_vector *data, char *separator, Output_declares){
    FILE *f = output_pipe;
    if (!data) fprintf(f, "NULL\n");
    else {
	    for (size_t i=0; i<data->size; i++){
		    p_fn(f, gsl_vector_get(data, i));
		    if (i< data->size -1) fprintf(f, "%s", separator);
	    }
	    fprintf(f,"\n");
    }
	if (output_name) fclose(f);
}

/** Print a vector in float format.
You may want to set \ref apop_opts_type "apop_opts.output_delimiter"; the default is a tab, which puts the vector on one line, but a newline would print the vector vertically.

\li See \ref apop_prep_output for more on how printing settings are set.
\li See also the legible output section of the \ref outline for more details and examples.
\li This function uses the \ref designated syntax for inputs.
\ingroup apop_print */
#ifdef APOP_NO_VARIADIC
void apop_vector_print(gsl_vector *data, Output_declares){
#else
apop_varad_head(void, apop_vector_print){
    gsl_vector *apop_varad_var(data, NULL);
    Dispatch_output
     apop_vector_print_base(data, Output_vars);
}

 void apop_vector_print_base(gsl_vector *data, Output_declares){
#endif
	print_core_v(data, apop_opts.output_delimiter, Output_vars);
 }

/** Dump a <tt>gsl_vector</tt> to the screen. 
    You may want to set \ref apop_opts_type "apop_opts.output_delimiter".

\li See \ref apop_prep_output for more on how printing settings are set.
\li See also the legible output section of the \ref outline for more details and examples.
\li This function uses the \ref designated syntax for inputs.
\ingroup apop_print */
void apop_vector_show(const gsl_vector *data){
	print_core_v(data, apop_opts.output_delimiter, NULL, stdout, 's', 0); 
}

static int get_max_strlen(char **names, size_t len){
    int max  = 0;
    for (int i=0; i< len; i++)
        max = GSL_MAX(max, strlen(names[i]));
    return max;
}

//On screen, display a pipe, else use the usual output delimiter.
static void a_pipe(FILE *f, char displaytype){
    if (displaytype == 's') fprintf(f, " | ");
    else                    fprintf(f, "%s", apop_opts.output_delimiter);
}

static void apop_data_print_core(const apop_data *data, FILE *f, char displaytype){
    if (!data){
        fprintf(f, "NULL\n");
        return;
    }
    int i, j, L = 0, 
        start   = (data->vector)? -1 : 0,
        end     = (data->matrix)? data->matrix->size2 : 0,
        rowend  = (data->matrix)? data->matrix->size1 : (data->vector) ? data->vector->size : data->text ? data->textsize[0] : -1;
    if (data->names->title && strlen(data->names->title))
        fprintf(f, "\t%s\n\n", data->names->title);
    if (data->names->rowct)
        L   = get_max_strlen(data->names->row, data->names->rowct);
    if (data->names->rowct && (data->names->vector || data->names->colct || data->names->textct))
        fprintf(f, "%*s  ", L+2, " ");
    if (data->vector && data->names->vector){
        fprintf(f, "%s", data->names->vector);
    }
    if (data->matrix){
        if (data->vector && data->names->colct){
            fprintf(f, "%c ", data->names->vector ? ' ' : '\t' );
            a_pipe(f, displaytype);
        }
        for(i=0; i< data->names->colct; i++){
            if (i < data->names->colct -1)
                fprintf(f, "%s%s", data->names->col[i], apop_opts.output_delimiter);
            else
                fprintf(f, "%s", data->names->col[i]);
        }
    }
    if (data->textsize[1] && data->names->textct){
        if ((data->vector && data->names->vector) || (data->matrix && data->names->colct))
            a_pipe(f, displaytype);
        for(i=0; i< data->names->textct; i++){
            if (i < data->names->textct -1)
                fprintf(f, "%s%s", data->names->text[i], apop_opts.output_delimiter);
            else
                fprintf(f, "%s", data->names->text[i]);
        }
    }
    if(data->names->vector || data->names->colct || data->names->textct)
        fprintf(f, "\n");
    for(j=0; j< rowend; j++){
        if (data->names->rowct > j)
            fprintf(f, "%*s%s", L+2, data->names->row[j], apop_opts.output_delimiter);
        for(i=start; i< end; i++){
            if ((i < 0 && j < data->vector->size) || (i>= 0 && j < data->matrix->size1 && i < data->matrix->size2))
                p_fn(f,  apop_data_get(data, j, i));
            else
                fprintf(f, " ");
            if (i==-1 && data->matrix) 
                a_pipe(f, displaytype);
            if (i < end-1)
                fprintf(f, "%s", apop_opts.output_delimiter);
        }
        if (data->text){
            if (data->vector || data->matrix)
                a_pipe(f, displaytype);
            if (j < data->textsize[0])
                for(i=0; i< data->textsize[1]; i++){
                    fprintf(f, "%s", data->text[j][i]);
                    if (i < data->textsize[1]-1) fprintf(f, "%s", apop_opts.output_delimiter);
                }
        }
        if (data->weights && j < data->weights->size){
            a_pipe(f, displaytype);
            p_fn(f, data->weights->data[j]);
        }
        fprintf(f, "\n");
    }
}

/** Print an \ref apop_data set to a file, the database, or the screen,
  as determined by the \c .output_type.

\li See \ref apop_prep_output for more on how printing settings are set.
\li See also the legible output section of the \ref outline for more details and examples.
\li This function uses the \ref designated syntax for inputs.
\ingroup apop_print */
#ifdef APOP_NO_VARIADIC
void apop_data_print(const apop_data *data, Output_declares){
#else
apop_varad_head(void, apop_data_print){
    const apop_data * apop_varad_var(data, NULL);
    Dispatch_output
     apop_data_print_base(data, Output_vars);
}

 void apop_data_print_base(const apop_data *data, Output_declares){
#endif 
    if (output_type  == 'd'){
        if (output_append == 'w') apop_table_exists(output_name, 'd');
        apop_data_to_db(data, output_name, output_append);
        return;
    }
    apop_data_print_core(data, output_pipe, output_type);
    if (data && data->more) {
        output_append='a';
        apop_data_print(data->more, Output_vars);
    }
    if (output_name)
        fclose(output_pipe);
}

/** Print a matrix in float format.
    You may want to set \ref apop_opts_type "apop_opts.output_delimiter".

\li See \ref apop_prep_output for more on how printing settings are set.
\li See also the legible output section of the \ref outline for more details and examples.
\li This function uses the \ref designated syntax for inputs.
\ingroup apop_print */
#ifdef APOP_NO_VARIADIC
void apop_matrix_print(const gsl_matrix *data, Output_declares){
#else
apop_varad_head(void, apop_matrix_print){
    const gsl_matrix *apop_varad_var(data, NULL);
    Dispatch_output
     apop_matrix_print_base(data, Output_vars);
}

 void apop_matrix_print_base(const gsl_matrix *data, Output_declares){
#endif
    if (output_type == 'd'){
        Apop_assert_c(data, , 1, "You sent me a NULL matrix. No database table will be created.");
    } else if (!data){
        fprintf(output_pipe, "NULL\n");
        return;
    }
    apop_data *d = apop_data_alloc();
    d->matrix=(gsl_matrix *) data; //cheating on the const qualifier
    apop_data_print(d, Output_vars);
    d->matrix=NULL;
    apop_data_free(d);
}

/** Dump a <tt>gsl_matrix</tt> to the screen.
    You may want to set \ref apop_opts_type "apop_opts.output_delimiter".
\li This function uses the \ref designated syntax for inputs.
\ingroup apop_print */
void apop_matrix_show(const gsl_matrix *data){
    apop_data *dtmp = apop_matrix_to_data((gsl_matrix*) data);
    apop_data_print_core(dtmp,  stdout, 's');
    dtmp->matrix = NULL;
    apop_data_free(dtmp);
}

/* the next function plots a single graph for the \ref apop_plot_lattice  fn */
static void printone(FILE *f, double width, double height, double margin, int xposn, int yposn, const apop_data *d){
    //pull two columns
    size_t count    = d->matrix->size2;
    double nudge    = 0.08;
    gsl_vector  v1  = gsl_matrix_column(d->matrix, xposn).vector;
    gsl_vector  v2  = gsl_matrix_column(d->matrix, yposn).vector;
    gsl_matrix  *m  = gsl_matrix_alloc(d->matrix->size1, 2);
    gsl_matrix_set_col(m, 0, &v1);
    gsl_matrix_set_col(m, 1, &v2);
    double sizex    = (double)(width - margin * (count -1))/count;
    double sizey    = (double)(height - margin * (count -1))/count;
    double offx     = width - (sizex +margin-nudge)* (count - xposn) - nudge*count;
    double offy     = height - (sizey +margin-nudge)* (1 + yposn) - nudge*count;
        fprintf(f, "unset y2tics; unset y2label; unset ytics; unset ylabel\n");
        fprintf(f, "unset x2tics; unset x2label; unset xtics; unset xlabel\n");
    fprintf(f, "set size   %g, %g\n"
               "set origin %g, %g\n"
                ,  sizex, sizey, offx, offy);
    //Gnuplot has no way to just set a label with no plot, so 
    //labels have to be drawn with a long offset from one of the
    //plots we do display.
    if (xposn  == yposn+1)
        fprintf(f, "set label %i '%s' center at graph %g, %g\n",yposn+1, (d->names->colct >yposn)? d->names->col[yposn]: "",  -0.5, 0.5);
    if ((yposn == count-1) && (xposn  == count-2))
        fprintf(f, "set label %zu '%s' center at graph %g, %g\n",count, (d->names->colct >count -1)? d->names->col[count -1]: "",  1.5, 0.5);
    if (xposn != yposn){
        fprintf(f, "plot '-'\n");
        fflush(f);
        apop_matrix_print(m, NULL, f, 'p');
        fprintf(f,"e\n");
        gsl_matrix_free(m);
    } 
    if (xposn  == yposn+1)
        fprintf(f, "unset label %i \n",yposn+1);
    if ((yposn == count-1) && (xposn  == count-2))
        fprintf(f, "unset label %zu\n",count);
}

/** This produces a Gnuplot file that will produce an array of 2-D
 plots, one for each pair of columns in the data set. Along the diagonal
 is a plot of the variable against itself---a density plot of the variable.

 \param d       The data set whose (matrix) columns will be compared. (No default, must not be \c NULL.)

\image latex "lattice.png" "A lattice showing three variables graphed against each other."
\image html "lattice.png" "A lattice showing three variables graphed against each other."

\li See \ref apop_prep_output for more on how printing settings are set.
\li See also the legible output section of the \ref outline for more details and examples.
\li This function uses the \ref designated syntax for inputs.
\ingroup output
*/
#ifdef APOP_NO_VARIADIC
void apop_plot_lattice(const apop_data *d, Output_declares){
#else
apop_varad_head(void, apop_plot_lattice){
    const apop_data * apop_varad_var(d, NULL);
    Apop_assert_n(d, "Input data set is NULL.\n");
    Dispatch_output
     apop_plot_lattice_base(d, Output_vars);
}

 void apop_plot_lattice_base(const apop_data *d, Output_declares){
#endif
    double  width  = 1.2,//these used to be options, but who's ever gonna set them to something else.
            height = 1.2;
    double  margin = 0;
    FILE *f = output_pipe;
    fprintf(f, "set size %g, %g\n"
                "set rmargin 5\n"
                "set lmargin -1\n"
                "set tmargin 2.4\n"
                "set bmargin -2\n"
                "set origin %g, %g       \n"
                "set multiplot   #layout %zu, %zu downwards        \n"
                "unset xtics; unset xlabel; unset ytics; unset ylabel\n"
                "set nokey           \n"
        , width, height, margin,margin, d->matrix->size2, d->matrix->size2);
    for (size_t i = 0; i< d->matrix->size2; i++)
        for (size_t j = 0; j< d->matrix->size2; j++)
            printone(f, width, height, margin, i, j, d);
    fprintf(f, "unset multiplot\n"); 
    if (output_type == 'f' && output_name)
        fclose(f);
}

/** Plot the percentiles of a data set against the percentiles of a distribution.

The distribution percentiles will be on the $x$-axis, your data percentiles on the $y$-.

\param v    The data (No default, must not be \c NULL.)
\param m    The distribution, such as apop_normal. I'll be using the \c draw method. (Default = best-fitting Normal)
\param bins The number of bins in the histogram. The number of points on the plot will always be 101 (i.e. percentiles). (default = MAX(10, data->size/10); denominator subject to future adjustment)
\param r    A \c gsl_rng. If NULL, I'll get an RNG via \ref apop_rng_get_thread. (Default = \c NULL)

\li See \ref apop_prep_output for more on how printing settings are set.
\li See also the legible output section of the \ref outline for more details and examples.
\li This function uses the \ref designated syntax for inputs.
*/
#ifdef APOP_NO_VARIADIC
void apop_plot_qq(gsl_vector *v, apop_model *m, Output_declares, size_t bins, gsl_rng *r){
#else
apop_varad_head(void, apop_plot_qq){
    int free_m = 0;
    gsl_vector * apop_varad_var(v, NULL);
    Apop_assert_n(v, "Input vector is NULL.");
    apop_model  *apop_varad_var(m, NULL);
    if (!m){
        free_m++;
        apop_data *d = apop_vector_to_data(v);
        m = apop_estimate(d, apop_normal);
        d->vector = NULL;
        apop_data_free(d);
    }
    Dispatch_output
    size_t apop_varad_var(bins, GSL_MAX(10, v->size/10));
    gsl_rng *apop_varad_var(r, apop_rng_get_thread())

    apop_plot_qq_base(v, m, Output_vars, bins, r);
    if (free_m) apop_model_free(m);
    return;
     apop_plot_qq_base(v, m, Output_vars, bins, r);
}

 void apop_plot_qq_base(gsl_vector *v, apop_model *m, Output_declares, size_t bins, gsl_rng *r){
#endif
    double *pctdata = apop_vector_percentiles(v, 'a');

    //produce percentiles from the model via RNG.
    gsl_vector  *vd  = gsl_vector_alloc(bins);
    for(int i=0; i< bins; i++)
        m->draw(gsl_vector_ptr(vd, i), r, m);
    double *pctdist = apop_vector_percentiles(vd, 'a');

    fprintf(output_pipe, "set key off; set size square;\n"
               "plot x;\n"
               "replot '-' with points\n");
    for (int i=0; i < 101; i++)
        fprintf(output_pipe, "%g\t %g\n", pctdist[i], pctdata[i]);
    fprintf(output_pipe, "e\n");
    if (output_type == 'p')
        fflush(output_pipe);
    else if (output_name)
        fclose(output_pipe);
    gsl_vector_free(vd);
}

/** This produces a nifty triangle plot from an input matrix with three
 columns. Each row is plotted inside an equilateral triangle such that
 (1, 0, 0) is the lower left corner, (0, 1, 0) the lower right corner,
 and (0, 0, 1) the middle, upper corner. 

\image html "triangle.png" "A triangle plot with a smattering of data points."

\li Gnuplot will rescale for you, so you don't need to worry about whether the row sums to one. 
\li See \ref apop_prep_output for more on how printing settings are set.
\li See also the legible output section of the \ref outline for more details and examples.
\li This function uses the \ref designated syntax for inputs.
*/
#ifdef APOP_NO_VARIADIC
void apop_plot_triangle(apop_data *in, Output_declares){
#else
apop_varad_head(void, apop_plot_triangle){
    apop_data *apop_varad_var(in, NULL);
    Apop_assert_n(in, "You sent me a NULL data set.");
    Dispatch_output
     apop_plot_triangle_base(in, Output_vars);
}

 void apop_plot_triangle_base(apop_data *in, Output_declares){
#endif 
    FILE *f=output_pipe;
    Apop_assert_n(f, "Error opening file %s for writing.", output_name);
    if (in->names && in->names->colct>=3){
        fprintf(f, "set label '%s' at -0.03, 0 right; \n", in->names->col[0]);
        fprintf(f, "set label '%s' at 1.03, 0 left; \n", in->names->col[1]);
        fprintf(f, "set label '%s' at 0.5, 1/sqrt(2)+0.05 center; \n", in->names->col[2]);
    }
    fprintf(f, 
        " set size square;      \n"
        " set noborder; set nogrid;         \n"
        " set noxtics; set noytics;         \n"
        " set nokey;        \n"
        " set lmargin 10;       \n"
        " set bmargin 3;        \n"
        " set tmargin 5;        \n"
        " set arrow 1 from 0,0 to 1,0 nohead        \n"
        " set arrow 2 from 0.5,1/sqrt(2) to 1,0 nohead      \n"
        " set arrow 3 from 0.5,1/sqrt(2) to 0,0 nohead      \n"
        " unset title       \n"
        " set xrange [-.1:1.1]      \n"
        " set yrange [-.1:0.81]     \n"
        " plot '-' using (($2 + $3/2)/($1+$2+$3)):($3/($1+$2+$3)/sqrt(2)) with points;         \n"
    );
    Apop_submatrix(in->matrix, 0,0, in->matrix->size1, 3, triplets);
    apop_matrix_print(triplets, .output_pipe=f, .output_type='p');
    if (output_name) fclose(f);
}
