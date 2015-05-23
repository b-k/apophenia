/** \file 
  Some printing and output interface functions. */
/* Copyright (c) 2006--2007, 2009 by Ben Klemens.  Licensed under the GPLv2; see COPYING.  */

//The reader will find a few function headers for this file in asst.h
#include "apop_internal.h"
 
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

\li Tip: if writing to the database, you can get a major speed boost by wrapping the call in a begin/commit wrapper:

\code
apop_query("begin;");
apop_data_print(your_data, .output_name="dbtab", .output_type='d');
apop_query("commit;");
\endcode
\ingroup all_public
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


/////The printing functions.

static void white_pad(int ct){
    for(size_t i=0; i < ct; i ++)
        printf(" ");
}

/** This function prettyprints the \c apop_data set to a screen.

This takes a lot of machinery. I write every last element to a text array, then measure column widths, then print to screen with padding to guarantee that everything lines up.  There's no way to have the first element of a column line up with the last unless you interrogate the width of every element in the column, so printing columns really can't be a one-pass process.

So, I produce an \ref apop_data set with no numeric elements and a text element to be filled with the input data set, and then print that. That means that I'll be using (more than) twice the memory to print this. If this is a problem, you can use \ref apop_print to dump your data to a text file, and view the text file, or print subsets.

For more machine-readable printing, see \ref apop_data_print.
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
\li See also \ref Legi for more details and examples.
\li This function uses the \ref designated syntax for inputs.
\ingroup all_public
*/
APOP_VAR_HEAD void apop_vector_print(gsl_vector *data, Output_declares){
    gsl_vector *apop_varad_var(data, NULL);
    Dispatch_output
APOP_VAR_ENDHEAD
	print_core_v(data, apop_opts.output_delimiter, Output_vars);
 }

/** Dump a <tt>gsl_vector</tt> to the screen. 
    You may want to set \ref apop_opts_type "apop_opts.output_delimiter".

\li See \ref apop_prep_output for more on how printing settings are set.
\li See also \ref Legi for more details and examples.
\li This function uses the \ref designated syntax for inputs.
*/
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
    if (data->names && data->names->title && strlen(data->names->title))
        fprintf(f, "\t%s\n\n", data->names->title);
    if (data->names && data->names->rowct)
        L   = get_max_strlen(data->names->row, data->names->rowct);
    if (data->names && data->names->rowct && (data->names->vector || data->names->colct || data->names->textct)){
        if (*apop_opts.db_name_column=='\0' || 
                !strcmp(apop_opts.db_name_column, "row_names"))
            fprintf(f, "%*s  ", L+2, " ");
        else { fprintf(f, "%s", apop_opts.db_name_column); a_pipe(f, displaytype); }
    }
    if (data->vector && data->names && data->names->vector){
        fprintf(f, "%s", data->names->vector);
    }
    if (data->matrix){
        if (data->vector && data->names && data->names->colct){
            fprintf(f, "%c ", data->names->vector ? ' ' : '\t' );
            a_pipe(f, displaytype);
        }
        if (data->names) 
          for(i=0; i< data->names->colct; i++){
            if (i < data->names->colct -1)
                fprintf(f, "%s%s", data->names->col[i], apop_opts.output_delimiter);
            else
                fprintf(f, "%s", data->names->col[i]);
        }
    }
    if (data->textsize[1] && data->names && data->names->textct){
        if ((data->vector && data->names && data->names->vector) || (data->matrix && data->names->colct))
            a_pipe(f, displaytype);
        if (data->names)
          for(i=0; i< data->names->textct; i++){
            if (i < data->names->textct -1)
                fprintf(f, "%s%s", data->names->text[i], apop_opts.output_delimiter);
            else
                fprintf(f, "%s", data->names->text[i]);
        }
    }
    if(data->names && (data->names->vector || data->names->colct || data->names->textct))
        fprintf(f, "\n");
    for(j=0; j< rowend; j++){
        if (data->names && data->names->rowct > j)
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
\li See also \ref Legi for more details and examples.
\li This function uses the \ref designated syntax for inputs.
\ingroup all_public
*/
APOP_VAR_HEAD void apop_data_print(const apop_data *data, Output_declares){
    const apop_data * apop_varad_var(data, NULL);
    Dispatch_output
APOP_VAR_ENDHEAD 
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
\li See also \ref Legi for more details and examples.
\li This function uses the \ref designated syntax for inputs.
\ingroup all_public
*/
APOP_VAR_HEAD void apop_matrix_print(const gsl_matrix *data, Output_declares){
    const gsl_matrix *apop_varad_var(data, NULL);
    Dispatch_output
APOP_VAR_ENDHEAD
    if (output_type == 'd'){
        Apop_assert_c(data, , 1, "You sent me a NULL matrix. No database table will be created.");
    } else if (!data){
        fprintf(output_pipe, "NULL\n");
        return;
    }
    apop_data_print(&(apop_data){.matrix=(gsl_matrix*)data}, Output_vars); //cheating on the const qualifier
}

/** Convenience function to dump a <tt>gsl_matrix</tt> to the screen.
*/
void apop_matrix_show(const gsl_matrix *data){
    apop_data_print_core(&(apop_data){.matrix=(gsl_matrix*)data},  stdout, 's');
}
