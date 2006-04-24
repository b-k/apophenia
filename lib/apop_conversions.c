/** \file apop_conversions.c	The various functions to convert from one format to another.

Copyright (c) 2006 by Ben Klemens. Licensed under the GNU GPL v2.
 */
#include "conversions.h"
#include "assert.h"

#define Text_Line_Limit 100000

/** \defgroup conversions Conversion functions

The functions to shunt data between text files, database tables, GSL matrices, and plain old arrays.*/
/** \defgroup convertfromarray  Functions to convert from an array (i.e., double **).
\ingroup conversions*/
/** \defgroup convertfromtext Functions to convert from a text file
\ingroup conversions*/
/** \defgroup convertfromdb  Functions to convert from a database
\ingroup conversions*/
/** \defgroup convertfrommatrix  Functions to convert from a gsl_matrix.
See also \ref output for funcions to write a matrix to a text file.
\ingroup conversions*/
/** \defgroup convertfromvector  Functions to convert from a gsl_vector.
See also \ref output for funcions to write a vector to a text file.
\ingroup conversions*/

/** \page gsl_views 	Using GSL Views
The GSL includes a convenient structure for pulling a vector from a matrix. Here's how to get the fifth row of <tt>a_matrix</tt> into a vector view:

\code
gsl_vector_view v;
v = gsl_matrix_col(a_matrix, 4);
\endcode

For rows, use <tt>gsl_matrix_row(a_matrix, n)</tt>. The vector view is a
data structure which includes an element of type <tt>gsl_vector</tt> named
<tt>vector</tt>; this is the only element you will be interested in. The
expression <tt>&(v.vector)</tt> is of type <tt>gsl_vector *</tt>, and therefore
can be used as you would any other pointer to a <tt>gsl_vector</tt>. For
example, try \ref apop_mean<tt>(&(v.vector));</tt>.

The view is intended to be a common variable, not a pointer. If you want
to retain the data after the function exits, copy it to another vector: 
\code
gsl_vector_view v;
gsl_vector *a_new_vector = gsl_vector_alloc(a_matrix->size1);
v = gsl_matrix_col(a_matrix, 4);
gsl_vector_memcpy(a_new_vector, &(v.vector));
\endcode
\ingroup convertfrommatrix
*/


/** Um, converts a GSL vector to an array.

\param in
A GSL vector

\param out
A pointer to a <tt>double*</tt>, which will be <tt>malloc</tt>ed inside the function.

\return Returns the size of the vector, i.e., <tt>in->size</tt>.

\note Does not use memcpy, because we don't know the stride of the vector.
\ingroup convertfromvector 
*/
int apop_vector_to_array(gsl_vector *in, double **out){
int		i;	
	*out	= malloc(sizeof(double) * in->size);
	for (i=0; i < in->size; i++)
		(*out)[i]	= gsl_vector_get(in, i);
	return in->size;
}

/** Just copies a one-dimensional array to a <tt>gsl_vector</tt>. The input array is undisturbed.

\param in 	    A vector.
\param size 	You will have to tell the function how long <tt>in</tt> is.
\return         A <tt>gsl_vector</tt>. Declare but do not allocate.
\ingroup convertfromarray 
*/ 
gsl_vector * apop_array_to_vector(double *in, int size){
int		    i;
gsl_vector  *out;
    out	= gsl_vector_alloc(size);
	for(i=0; i < size; i++)
		gsl_vector_set(out, i, in[i]);
    return out;
}
/*
void apop_array_to_vector(double *in, gsl_vector **out, int size){
int		i;
	*out	= gsl_vector_alloc(size);
	for(i=0; i < size; i++)
		gsl_vector_set(*out, i, in[i]);
}
*/

static void convert_array_to_line(double **in, double **out, int rows, int cols){
	//go from in[i][j] form to the GSL's preferred out[i*cols + j] form
int		i, j;
	*out	= malloc(sizeof(double) * rows * cols);
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			(*out)[i * cols + j]	= in[i][j];
}

/** convert a <tt>double **</tt> array to a <tt>gsl_matrix</tt>

\param in	the array to read in
\param rows, cols	the size of the array.
\return the <tt>gsl_matrix</tt>, allocated for you and ready to use.

usage: \code gsl_matrix *m = apop_array_to_matrix(indata, 34, 4); \endcode

If you want to intialize on the allocation line, this isn't what you want. See \ref apop_line_to_matrix.
\ingroup convertfromarray 
*/
gsl_matrix * apop_array_to_matrix(double **in, int rows, int cols){
gsl_matrix_view	m;
gsl_matrix      *out;
double		    *line;
	out	= gsl_matrix_alloc(rows, cols);
	convert_array_to_line(in, &line, rows, cols);
	m	= gsl_matrix_view_array(line, rows,cols);
	gsl_matrix_memcpy(out,&(m.matrix));
	free(line);
    return out;
}

/** convert a <tt>double **</tt> array to an \ref apop_data set. It will
have no names.

\param in	the array to read in
\param rows, cols	the size of the array.
\return the \ref apop_data set, allocated for you and ready to use.

usage: \code apop_data *d = apop_array_to_data(indata, 34, 4); \endcode
\ingroup convertfromarray 

If you want to intialize on the allocation line, this isn't what you want. See \ref apop_line_to_data.
*/
apop_data * apop_array_to_data(double **in, int rows, int cols){
    return apop_matrix_to_data(apop_array_to_matrix(in, rows, cols));
}

/** convert a <tt>double **</tt> array to a <tt>gsl_matrix</tt>

\param line	the array to read in
\param rows, cols	the size of the array.
\return the <tt>gsl_matrix</tt>, allocated for you and ready to use.

usage: \code gsl_matrix *m = apop_array_to_matrix(indata, 34, 4); \endcode
\ingroup convertfromarray 
*/
gsl_matrix * apop_line_to_matrix(double *line, int rows, int cols){
gsl_matrix_view	m;
gsl_matrix      *out;
	out	= gsl_matrix_alloc(rows, cols);
	m	= gsl_matrix_view_array(line, rows,cols);
	gsl_matrix_memcpy(out,&(m.matrix));
    return out;
}

/** convert a <tt>double **</tt> array to an \ref apop_data set. It will
have no names.

\param in	the array to read in
\param rows, cols	the size of the array.
\return the \ref apop_data set, allocated for you and ready to use.

usage: \code apop_data *d = apop_array_to_data(indata, 34, 4); \endcode
\ingroup convertfromarray 
*/
apop_data * apop_line_to_data(double *in, int rows, int cols){
    return apop_matrix_to_data(apop_line_to_matrix(in, rows, cols));
}

static int find_cat_index(char **d, char * r, int start_from, int size){
//used for apop_db_to_crosstab.
int	i	= start_from;	//i is probably the same or i+1.
	do {
		if(!strcmp(d[i], r))
			return i;
		i	++;
		i	%= size;	//loop around as necessary.
	} while(i!=start_from); 
	printf(" apop %s, %i: something went wrong in the crosstabbing; couldn't find %s.\n", __FILE__, __LINE__, r);
	return 0;
}

/**Give the name of a table in the database, and names of three of its
columns: the x-dimension, the y-dimension, and the data.
the output is a 2D matrix with rows indexed by r1 and cols by
r2. if !=NULL, d1 and d2 will list the labels on the dimensions.

\ingroup db
*/
apop_data  *apop_db_to_crosstab(char *tabname, char *r1, char *r2, char *datacol){
gsl_matrix	*out;
int		    i	= 0,
		    j	= 0,
		    k, ct_r, ct_c, datasize; 
char        ***pre_d1, ***pre_d2, ***datachars;
apop_data   *outdata    = apop_data_alloc(1,1);
	datachars	= apop_query_to_chars("select %s, %s, %s from %s", r1, r2, datacol, tabname);
    datasize    = apop_db_get_rows();   

    //A bit inefficient, but well-encapsulated.
    //Pull the distinct (sorted) list of headers, copy into outdata->names.
	pre_d1	    = apop_query_to_chars("select distinct %s, 1 from %s order by %s", r1, tabname, r1);
    ct_r        = apop_db_get_rows();
	if (pre_d1 == NULL) 
		printf (" apop %s, %i: selecting %s from %s returned an empty table.\n", __FILE__, __LINE__, r1, tabname);
    for (i=0; i < ct_r; i++)
        apop_name_add(outdata->names, pre_d1[i][0], 'r');
    apop_cats_free(pre_d1, ct_r, 1);

	pre_d2	= apop_query_to_chars("select distinct %s from %s order by %s", r2, tabname, r2);
    ct_c    = apop_db_get_rows();
	if (pre_d2 == NULL) 
		printf (" apop %s, %i: selecting %s from %s returned an empty table.\n", __FILE__, __LINE__, r2, tabname);
    for (i=0; i < ct_c; i++)
        apop_name_add(outdata->names, pre_d2[i][0], 'c');
    apop_cats_free(pre_d2, ct_c, 1);

	out	= gsl_matrix_calloc(ct_r, ct_c);
	for (k =0; k< datasize; k++){
		i	= find_cat_index(outdata->names->rownames, datachars[k][0], i, ct_r);
		j	= find_cat_index(outdata->names->colnames, datachars[k][1], j, ct_c);
		gsl_matrix_set(out, i, j, atof(datachars[k][2]));
	}
    apop_cats_free(datachars, datasize, 3);
    gsl_matrix_free(outdata->matrix);
    outdata->matrix   = out;
	return outdata;
}

/*
Much of the magic below is due to the following regular expression.

Without backslashes and spaced out in perl's /x style, it would look like this:
("[^"][^"]*"        (starts with a ", has no "" in between, ends with a ".
|               or
[^%s"][^%s"]*)  anything but a "" or the user-specified delimiters. At least 1 char long.)
([%s\n]|$)      and ends with a delimiter or the end of line.
*/
const char      divider[]="\\(\"[^\"][^\"]*\"\\|[^\"%s][^\"%s]*\\)[%s\n]";

//in: the line being read, the allocated outstring, the result from the regexp search, the offset
//out: the outstring is filled with a bit of match, last_match is updated.
static void pull_string(char *line, char * outstr, regmatch_t *result, size_t * last_match){
int     length_of_match = result[1].rm_eo - result[1].rm_so;
    memcpy(outstr, line + (*last_match)+result[1].rm_so, length_of_match);
    outstr[length_of_match]       = '\0';
    (*last_match)                += result[1].rm_eo+1;
}

/** Open file, find the first non-comment row, count columns, close file.
 */
int apop_count_cols_in_text(char *text_file){
FILE * 		infile;
char		instr[Text_Line_Limit], outstr[Text_Line_Limit],
            full_divider[1000];
int		    ct	                = 0,
            length_of_string    = 0;
size_t      last_match          = 0;
regex_t     *regex              = malloc(sizeof(regex_t));
regmatch_t  result[3];
    sprintf(full_divider, divider, apop_opts.input_delimiters, apop_opts.input_delimiters, apop_opts.input_delimiters);
    regcomp(regex, full_divider, 0);
	infile	= fopen(text_file,"r");
    if (infile == NULL){
        printf("Error opening file %s. apop_count_cols_in_text returning 0.", text_file);
        return 0;
    }
	if (infile==NULL) {printf("Error opening %s", text_file); return 1;}
	fgets(instr, Text_Line_Limit, infile);
	while(instr[0]=='#')	//burn off comment lines
		fgets(instr, Text_Line_Limit, infile);
    length_of_string= strlen(instr);
    while (last_match < length_of_string && !regexec(regex, (instr+last_match), 2, result, 0)){
        pull_string(instr,  outstr, result,  &last_match);
		ct++;
	}
	fclose(infile);
    regfree(regex);
	return ct;
}

/** Open file, count lines that don't start with #, close file.
 */
int apop_count_rows_in_text(char *text_file){
FILE * 		infile;
char		instr[Text_Line_Limit];
int		ct	= 0;
	infile	= fopen(text_file,"r");
	if (infile==NULL) {printf("Error opening %s", text_file); return 1;}
	while(fgets(instr,Text_Line_Limit,infile)!=NULL)
	    if(instr[0]!='#')	//burn off comment lines
		    ct  ++;
	fclose(infile);
	return ct-1;
}

//OK, OK. C sucks. This fn just strips leading and trailing blanks.
static char * strip(char *in){
regex_t     *regex      = malloc(sizeof(regex_t));
size_t      dummy       = 0;
char 		*out 	    = malloc(sizeof(char) * (1+strlen(in)));
char        stripregex[]= "[ \n\t]*\\([^ \n\t]*.*[^ \n\t]*\\)[ \n\t]*";
regmatch_t  result[3];
    regcomp(regex, stripregex, 0);
    if(!regexec(regex, in, 2, result, 0))
        pull_string(in, out, result,  &dummy);
    regfree(regex);
	return out;
}

/** \page text_format Notes on input text file formatting

If you are reading into an array or <tt>gsl_matrix</tt> or \ref
apop_data set, all text fields are taken as zeros. You will be warned
of such substitutions unless you set \ref apop_opts.verbose<tt>==0</tt>
beforehand.

You will also be interested in \ref apop_opts.input_delimiters. By
default, it is set to "| ,\t", meaning that a pipe,
comma, space, or tab will delimit separate entries. Try \code
strcpy(apop_opts.input_delimiters, ";")\endcode to set the delimiter to
a semicolon, for example.

There are often two delimiters in a row, e.g., "23, 32,, 12". When
it's two commas like this, the user typically means that there is a
missing value; when it is two tabs in a row, this is typically just
a formatting glitch. Apophenia is not smart enough to work out the difference:
a sequence of delimiters is always taken to be a single delimiter. If you do
have double-commas like this, you will have to explicitly insert a note
that there is a missing data point. Try:
\code
		perl -pi.bak -e 's/,,/,NaN,/g' data_file
\endcode

Text is always delimited by quotes. Delimiters inside quotes are perfectly
OK, e.g., "Males, 30-40", is an OK column name.

Lines beginning with # are taken to be comments and ignored. 

If there are row names and column names, then the input will not be
perfectly square: there should be no first entry in the row with column
names like 'row names'. That is, for a 100x100 data set with row and
column names, there are 100 names in the top row, and 101 entries in
each subsequent row (name plus 100 data points).

The maximum line length is 100,000 characters. If you have a line
longer than this, you will need to open up apop_conversions.c, modify
<tt>Text_Line_Limit</tt>, and recompile.

*/

/** Read a delimited text file into an array. 

  See \ref text_format.

\param text_file    The name of the text file to be read in.
\param has_row_names Does the lines of data have row names?
\param has_col_names Is the top line a list of column names? If there are row names, then there should be no first entry in this line like 'row names'. That is, for a 100x100 data set with row and column names, there are 100 names in the top row, and 101 entries in each subsequent row (name plus 100 data points).
\return 	Returns an apop_data set.

<b>example:</b> See \ref apop_OLS.
\ingroup convertfromtext	*/
apop_data * apop_text_to_data(char *text_file, int has_row_names, int has_col_names){
apop_data   *set;
FILE * 		infile;
char		instr[Text_Line_Limit], 
            *str, *stripped, 
            outstr[Text_Line_Limit],
            full_divider[1000];
int 		i	        = 0,
            line_no     = 0,
            length_of_string,
		    ct, colno, rowct;
size_t      last_match;
regex_t     *regex  = malloc(sizeof(regex_t));
regmatch_t  result[2];
	ct	    = apop_count_cols_in_text(text_file);
	rowct	= apop_count_rows_in_text(text_file);
    set     = apop_data_alloc(rowct+1-has_col_names,ct);
	infile	= fopen(text_file,"r");
    if (infile == NULL){
        printf("Error opening file %s. apop_text_to_data returning NULL.", text_file);
        return NULL;
    }
    sprintf(full_divider, divider, apop_opts.input_delimiters, apop_opts.input_delimiters, apop_opts.input_delimiters);
    regcomp(regex, full_divider, 0);

    //First, handle the top line, which is assumed to be column names.
    if (has_col_names){
	    fgets(instr, Text_Line_Limit, infile);
        line_no   ++;
	    while(instr[0]=='#')	//burn off comment lines
		    fgets(instr, Text_Line_Limit, infile);
	    set->names->colnamect= 0;
	    set->names->colnames	= malloc(sizeof(char*));
        last_match      = 0;
        length_of_string= strlen(instr);
        while (last_match < length_of_string && !regexec(regex, (instr+last_match), 2, result, 0)){
            pull_string(instr,  outstr, result,  &last_match);
	        stripped	= strip(outstr);
            apop_name_add(set->names, stripped, 'c');
		    free(stripped);
	    }
    }

    //Now do the body. First elmt may be a row name.
	while(fgets(instr,Text_Line_Limit,infile)!=NULL){
		colno	= 0;
		if(instr[0]!='#') {
			i	            ++;
            last_match      = 0;
            length_of_string= strlen(instr);
            regexec(regex, (instr+last_match), 2, result, 0);   //one for the headers.
			if (has_row_names){
                pull_string(instr,  outstr, result,  &last_match);
	            stripped	= strip(outstr);
                apop_name_add(set->names, stripped, 'r');
			    free(stripped);
			}
            while (last_match < length_of_string && !regexec(regex, (instr+last_match), 2, result, 0)){
                pull_string(instr,  outstr, result,  &last_match);
				colno++;
				gsl_matrix_set(set->matrix, i-1, colno-1,	 strtod(outstr, &str));
				if (apop_opts.verbose && !strcmp(outstr, str))
				    printf("trouble converting item %i on line %i; using zero.\n", colno, i);
			}
		}
	}
	fclose(infile);
    regfree(regex);
	return set;
}

/**If the string isn't a number, it needs quotes, and sqlite wants 0.1,
  not just .1. 
  \todo This could be easier with regexes.
 */
static char * prep_string_for_sqlite(char *astring){
char		*tmpstring, 
		*out	= NULL,
		*str	= NULL;
	strtod(astring, &str);
	if (!strcmp(astring, str)){	//then it's not a number.
		tmpstring=strip(astring);
		if (strlen (tmpstring)==0){
			out	= malloc(sizeof(char) * 2);
			sprintf(out, " ");
			free(tmpstring);
			return out;
		}
		if (tmpstring[0]!='"'){
			out	= malloc(sizeof(char) * (strlen(tmpstring)+3));
			sprintf(out, "\"%s\"",tmpstring);
			free(tmpstring);
			return out;
		} else return tmpstring;
	} else {			//sqlite wants 0.1, not .1
		tmpstring=strip(astring);
		assert(strlen (tmpstring)!=0);
		/*if (strlen (tmpstring)==0){
			out	= malloc(sizeof(char) * 2);
			free(tmpstring);
			sprintf(out, " ");
			return out;
		}*/
		if (tmpstring[0]=='.'){
			out	= malloc(sizeof(char) * (strlen(tmpstring)+2));
			sprintf(out, "0%s",tmpstring);
			free(tmpstring);
			return out;
		} else return tmpstring;
	} //If you're here, then it's a number which needs no fixing.
	out	= malloc(sizeof(char) * (strlen(astring)+1));
	strcpy(out, astring);
	return out;
}

/** This function will print a string to another string, allocating the
appropriate amount of space along the way.
 
That is, it will (1) reallocate base to exactly the needed length,
and then (2) write addme. 

\param base     The pointer to be written to. May be NULL. If base is
automatically allocated (i.e., you declared it with char base[]), then
this will crash.
\param addme    A string.
\return a pointer to base. 

\ingroup convenience_fns
*/
char *apop_strcpy(char **base, char *addme){
int     addlen  = (addme) ? strlen(addme): 0;
    *base    = realloc(*base, sizeof(char)*(addlen+1));
    if (!*base)
        printf("Ran out of memory in apop_strcpy. Returning NULL.\n");
    strcpy(*base, addme);
    return *base;
}

/** In the proud tradition of every library providing its own haphazard string handling functions, this function will safely append one string on to another.
 
That is, it will (1) reallocate base to exactly the needed length,
and then (2) append addme. 

\param base     The pointer to be extended. May be NULL. If base is
automatically allocated (i.e., you declared it with char base[]), then
this will crash.
\param addme    A string.
\return a pointer to base. 

\ingroup convenience_fns
*/
char *apop_strcat(char **base, char *addme){
    if (!*base){
        *base    = malloc(sizeof(char));
        (*base)[0] = '\0';
    }
int     baselen = strlen(*base),
        addlen  = (addme) ? strlen(addme): 0;
    *base    = realloc(*base, sizeof(char)*(baselen+addlen+1));
    if (!*base)
        printf("Ran out of memory in apop_strcat. Returning NULL.\n");
    strcat(*base, addme);
    return *base;
}

/** Read a text file into a database table.

  See \ref text_format.

\param text_file    The name of the text file to be read in.
\param tabname      The name to give the table in the database
\param has_row_names Does the lines of data have row names?
\param has_col_names Is the top line a list of column names? All dots in the column names are converted to underscores, by the way.
\param field_names The list of field names, which will be the columns for the table. If <tt>has_col_names==1</tt>, read the names from the file (and just set this to <tt>NULL</tt>). If has_col_names == 1 && field_names !=NULL, I'll use the field names. 

\return Returns the number of rows.

Using the data set from the example on the \ref apop_OLS "apop_OLS" page, here's another way to do the regression:

\code
#include <apophenia/headers.h>

int main(void){ 
apop_data       *data; 
apop_estimate   *est;
    apop_db_open(NULL);
    apop_text_to_db("data", "d", 0,1,NULL);
    data       = apop_query_to_data("select * from d");
    estimate   = apop_OLS.estimate(data, NULL, NULL);
    printf("The OLS coefficients:\n");
    apop_estimate_print(est);
    return 0;
} 
\endcode
\ingroup convertfromtext
*/
int apop_text_to_db(char *text_file, char *tabname, int has_row_names, int has_col_names, char **field_names){
FILE * 		infile;
char		*q  = NULL, instr[Text_Line_Limit], **fn, *prepped;
char		*stripped, *stripme, outstr[Text_Line_Limit],
            full_divider[1000];
int 		ct, one_in,
            length_of_string,
		    i			        = 0, 
		    use_names_in_file   = 0,
		    rows			    = 0;
size_t      last_match;
regex_t     *regex  = malloc(sizeof(regex_t));
regmatch_t  result[2];
	ct	= apop_count_cols_in_text(text_file);
	if (apop_table_exists(tabname,0)){
	       	printf("apop: %s table exists; not recreating it.\n", tabname);
		return 0; //to do: return the length of the table.
	} else{
        sprintf(full_divider, divider, apop_opts.input_delimiters, apop_opts.input_delimiters, apop_opts.input_delimiters);
        regcomp(regex, full_divider, 0);
		infile	= fopen(text_file,"r");
	       	if (infile==NULL) {
			printf("Trouble opening %s. apop_text_to_db bailing.\n", text_file);
			return 0;
		}
		if (has_col_names && field_names == NULL){
			use_names_in_file++;
			fgets(instr, Text_Line_Limit, infile);
			while(instr[0]=='#')	//burn off comment lines
				fgets(instr, Text_Line_Limit, infile);
			fn	= malloc(ct * sizeof(char*));

            last_match      = 0;
            length_of_string= strlen(instr);
            while (last_match < length_of_string && !regexec(regex, (instr+last_match), 2, result, 0)){
                pull_string(instr,  outstr, result,  &last_match);
	            stripme	    = strip(outstr);
                stripped    = apop_strip_dots(stripme,'d');
			    fn[i]	= NULL;
                /*
			    fn[i]	= malloc(1000 * sizeof(char));
			    strcpy(fn[i], stripped);
                */
                apop_strcpy(&fn[i], stripped);
                free(stripme);
		        free(stripped);
			    i++;
	        }
		} else	{
            if (field_names)
                fn	= field_names;
            else{
			    fn	= malloc(ct * sizeof(char*));
                for (i =0; i < ct; i++){
			        fn[i]	= malloc(1000 * sizeof(char));
                    sprintf(fn[i], "col_%i", i);
                }
            }
        }
		apop_strcpy(&q, "begin; CREATE TABLE ");
		apop_strcat(&q, tabname);
		for (i=0; i<ct; i++){
			if (i==0) 	{
                if (has_row_names)
                    apop_strcat(&q, " (row_names, ");
                else
                    apop_strcat(&q, " (");
            } else		apop_strcat(&q, " , ");
            apop_strcat(&q, " ");
            apop_strcat(&q, fn[i]);
		}
		apop_query_db("%s ); commit; begin;", q);
		if (use_names_in_file){
		    for (i=0; i<ct; i++)
			    free(fn[i]);
			free(fn);
		}

        //convert a data line into SQL: insert into TAB values (0.3, 7, "et cetera");
		while(fgets(instr,Text_Line_Limit,infile)!=NULL){
			if((instr[0]!='#') && (instr[0]!='\n')) {	//comments and blank lines.
				rows	        ++;
				one_in          = 
                last_match      = 0;
                length_of_string= strlen(instr);
				sprintf(q, "INSERT INTO %s VALUES (", tabname);
                while (last_match < length_of_string && !regexec(regex, (instr+last_match), 2, result, 0)){
					if(one_in++) 	apop_strcat(&q, ", ");
                    pull_string(instr,  outstr, result,  &last_match);
					prepped	=prep_string_for_sqlite(outstr);
                    if (strlen(prepped) > 0)
					    apop_strcat(&q, prepped);
					free(prepped);
				}
				apop_query_db("%s);",q);
			}
		}
		apop_query_db("commit;");
		fclose(infile);
        free(q);
		return rows;
	}
}

/** See \ref apop_db_to_crosstab for the storyline; this is the complement.
 \ingroup db
 */
int apop_crosstab_to_db(apop_data *in,  char *tabname, char *row_col_name, 
						char *col_col_name, char *data_col_name){
int		    i,j;
apop_name   *n = in->names;
	apop_query_db("CREATE TABLE %s (%s , %s , %s);", tabname, 
            apop_strip_dots(row_col_name, 'd'), 
            apop_strip_dots(col_col_name, 'd'), 
            apop_strip_dots(data_col_name, 'd'));
	for (i=0; i< n->colnamect; i++){
		apop_query_db("begin;");
		for (j=0; j< n->rownamect; j++)
			apop_query_db("INSERT INTO %s VALUES (%s, %s,%g);", tabname, 
					n->rownames[j], n->colnames[i], 
                                    gsl_matrix_get(in->matrix, j, i));
		apop_query_db("commit;");
	}
	return 0;
}

/** \page dbtomatrix converting from database table to <tt>gsl_matrix</tt> or \ref apop_data

Use <tt>fill_me = apop_query_to_matrix("select * from table_name;");</tt>
or <tt>fill_me = apop_query_to_data("select * from table_name;");</tt>. [See \ref apop_query_to_matrix; \ref apop_query_to_data.]
\ingroup convertfromdb
*/


/** Copy one  <tt>gsl_vector</tt> to another. That is, all data is duplicated.
 Unlike <tt>gsl_vector_memcpy</tt>, this function allocates and returns the destination, so you can use it like this:

 \code
 gsl_vector *a_copy = apop_vector_copy(original);
 \endcode

  \param in    the input data
  \return       a structure that this function will allocate and fill
\ingroup convenience_fns
  */
gsl_vector *apop_vector_copy(gsl_vector *in){
gsl_vector *out = gsl_vector_alloc(in->size);
    gsl_vector_memcpy(out, in);
    return out;
}

/** Copy one  <tt>gsl_matrix</tt> to another. That is, all data is duplicated.
 Unlike <tt>gsl_matrix_memcpy</tt>, this function allocates and returns the destination, so you can use it like this:

 \code
 gsl_matrix *a_copy = apop_matrix_copy(original);
 \endcode

  \param in    the input data
  \return       a structure that this function will allocate and fill
\ingroup convenience_fns
  */
gsl_matrix *apop_matrix_copy(gsl_matrix *in){
gsl_matrix *out = gsl_matrix_alloc(in->size1, in->size2);
    gsl_matrix_memcpy(out, in);
    return out;
}
