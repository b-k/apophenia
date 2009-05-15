/** \file apop_conversions.c	The various functions to convert from one format to another.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#include "conversions.h"
#include <gsl/gsl_math.h> //GSL_NAN
#include <assert.h>


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


/** Converts a GSL vector to an array.

\param in A GSL vector

\return A pointer to a <tt>double*</tt>, which will be <tt>malloc</tt>ed inside the function.

\note Does not use memcpy, because we don't know the stride of the vector.
\ingroup convertfromvector 
*/
double * apop_vector_to_array(const gsl_vector *in){
  apop_assert(in, NULL, 1, 'c', "You sent me a NULL vector; returning NULL");
  int		i;	
  double	*out	= malloc(sizeof(double) * in->size);
	for (i=0; i < in->size; i++)
		out[i]	= gsl_vector_get(in, i);
	return out;
}

/** Just copies a one-dimensional array to a <tt>gsl_vector</tt>. The input array is undisturbed.

\param line     A vector.
\param vsize 	You will have to tell the function how long \c line is.
\return         A <tt>gsl_vector</tt>. Declare but do not allocate.
\ingroup convertfromarray 
*/ 
gsl_vector * apop_array_to_vector(const double *line, const int vsize){
  apop_assert(line, NULL, 1, 'c', "You sent me NULL data; returning NULL.");
  apop_assert(vsize, NULL, 1, 'c', "You sent me vsize==0, so I'm returning NULL.");
  gsl_vector        *out    = gsl_vector_alloc(vsize);
  gsl_vector_view	v	    = gsl_vector_view_array((double*)line, vsize);
	gsl_vector_memcpy(out,&(v.vector));
    return out;
}

/** Mathematically, a vector of size \f$N\f$ and a matrix of size \f$N
 \times 1 \f$ are equivalent, but they're two different types in C. This function copies the data in a vector to a new one-column matrix and returns the newly-allocated and filled matrix.

 \param in a <tt>gsl_vector</tt>
 \return a newly-allocated <tt>gsl_matrix</tt> with one column.

 \ingroup convenience_fns
 */
gsl_matrix * apop_vector_to_matrix(const gsl_vector *in){
    apop_assert(in,  NULL, 1,'c', "Converting NULL vector to NULL matrix.");
    gsl_matrix *out = gsl_matrix_alloc(in->size, 1);
    apop_assert(out,  NULL, 0,'s', "gsl_matrix_alloc failed; probably out of memory.");
    gsl_matrix_set_col(out, 0, in);
    return out;
}

static void convert_array_to_line(const double **in, double **out, const int rows, const int cols){
	//go from in[i][j] form to the GSL's preferred out[i*cols + j] form
  int		i, j;
	*out	= malloc(sizeof(double) * rows * cols);
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			(*out)[i * cols + j]	= in[i][j];
}

/** Convert a <tt>double **</tt> array to a <tt>gsl_matrix</tt>

\param in	the array to read in
\param rows, cols	the size of the array.
\return the <tt>gsl_matrix</tt>, allocated for you and ready to use.

usage: \code gsl_matrix *m = apop_array_to_matrix(indata, 34, 4); \endcode

If you want to intialize on the allocation line, this isn't what you want. See \ref apop_line_to_matrix.
\ingroup convertfromarray 
*/
gsl_matrix * apop_array_to_matrix(const double **in, const int rows, const int cols){
  gsl_matrix_view   m;
  gsl_matrix        *out;
  double		    *line;
	out	= gsl_matrix_alloc(rows, cols);
	convert_array_to_line(in, &line, rows, cols);
	m	= gsl_matrix_view_array(line, rows,cols);
	gsl_matrix_memcpy(out,&(m.matrix));
	free(line);
    return out;
}

/** Convert a <tt>double **</tt> array to an \ref apop_data set. It will
have no names.

\param in	the array to read in
\param rows, cols	the size of the array.
\return the \ref apop_data set, allocated for you and ready to use.

usage: \code apop_data *d = apop_array_to_data(indata, 34, 4); \endcode
\ingroup convertfromarray 

If you want to intialize on the allocation line, this isn't what you want. See \ref apop_line_to_data.
*/
apop_data * apop_array_to_data(const double **in, const int rows, const int cols){
    return apop_matrix_to_data(apop_array_to_matrix(in, rows, cols));
}


/** Convert a <tt>double *</tt> array to a <tt>gsl_matrix</tt>

\param line	the array to read in
\param rows, cols	the size of the array.
\return the <tt>gsl_matrix</tt>, allocated for you and ready to use.

usage: \code gsl_matrix *m = apop_array_to_matrix(indata, 34, 4); \endcode
\ingroup convertfromarray 
*/
gsl_matrix * apop_line_to_matrix(double *line, int rows, int cols){
  gsl_matrix_view	m;
  gsl_matrix        *out;
	out	= gsl_matrix_alloc(rows, cols);
	m	= gsl_matrix_view_array(line, rows,cols);
	gsl_matrix_memcpy(out,&(m.matrix));
    return out;
}

/** Convert a <tt>double **</tt> array to an \ref apop_data set. It will
have no names.

\param in	The array to read in
\param vsize    The vector size
\param rows, cols	the size of the array.
\return the \ref apop_data set, allocated for you and ready to use.

\ingroup convertfromarray 
*/
apop_data * apop_line_to_data(double *in, int vsize, int rows, int cols){
    if (vsize==0 && (rows>0 && cols>0))
      return apop_matrix_to_data(apop_line_to_matrix(in, rows, cols));
    if ((rows==0 || cols==0) && vsize>0)
      return apop_vector_to_data(apop_array_to_vector(in, vsize));
    apop_assert(vsize==rows,  NULL, 0,'c',"apop_line_to_data expects either only a matrix, only a vector, or that matrix row count and vector size are equal. You gave me a row size of %i and a vector size of %i. Returning NULL.\n", rows, vsize);
  int i, j, ctr = 0;
  apop_data *out  = apop_data_alloc(vsize, rows, cols);
    for (i=0; i< rows; i++)
        for (j=-1; j< cols; j++)
            apop_data_set(out, i, j, in[ctr++]);
    return out;
}

static int find_cat_index(char **d, char * r, int start_from, int size){
//used for apop_db_to_crosstab.
  int	i	= start_from % size;	//i is probably the same or i+1.
	do {
		if(!strcmp(d[i], r))
			return i;
		i	++;
		i	%= size;	//loop around as necessary.
	} while(i!=start_from); 
    Apop_assert(0, 0, 0, 'c', "Something went wrong in the crosstabbing; couldn't find %s.", r);
}

/**Give the name of a table in the database, and names of three of its
columns: the x-dimension, the y-dimension, and the data.
the output is a 2D matrix with rows indexed by r1 and cols by
r2. if !=NULL, d1 and d2 will list the labels on the dimensions.

\ingroup db
*/
apop_data  *apop_db_to_crosstab(char *tabname, char *r1, char *r2, char *datacol){
  gsl_matrix	*out;
  int		i = 0,
		j = 0,
		k; 
  apop_data     *pre_d1, *pre_d2, *datachars;
  apop_data     *outdata    = apop_data_alloc(0,0,0);

    char p = apop_opts.db_name_column[0];
    apop_opts.db_name_column[0]= '\0';//we put this back at the end.
    datachars	= apop_query_to_text("select %s, %s, %s from %s", r1, r2, datacol, tabname);
    apop_assert(datachars,  NULL, 0, 's', "selecting %s, %s, %s from %s returned an empty table.\n",  r1, r2, datacol, tabname);

    //A bit inefficient, but well-encapsulated.
    //Pull the distinct (sorted) list of headers, copy into outdata->names.
    pre_d1	    = apop_query_to_text("select distinct %s, 1 from %s order by %s", r1, tabname, r1);
    apop_assert(pre_d1,  NULL, 0, 's', "selecting %s from %s returned an empty table.", r1, tabname);
    for (i=0; i < pre_d1->textsize[0]; i++)
        apop_name_add(outdata->names, pre_d1->text[i][0], 'r');

	pre_d2	= apop_query_to_text("select distinct %s from %s order by %s", r2, tabname, r2);
	apop_assert(pre_d2,  NULL, 0, 's', "selecting %s from %s returned an empty table.", r2, tabname);
    for (i=0; i < pre_d2->textsize[0]; i++)
        apop_name_add(outdata->names, pre_d2->text[i][0], 'c');

	out	= gsl_matrix_calloc(pre_d1->textsize[0], pre_d2->textsize[0]);
	for (k =0; k< datachars->textsize[0]; k++){
		i	= find_cat_index(outdata->names->row, datachars->text[k][0], i, pre_d1->textsize[0]);
		j	= find_cat_index(outdata->names->column, datachars->text[k][1], j, pre_d2->textsize[0]);
		gsl_matrix_set(out, i, j, atof(datachars->text[k][2]));
	}
    apop_data_free(pre_d1);
    apop_data_free(pre_d2);
    apop_data_free(datachars);
    outdata->matrix   = out;
    apop_opts.db_name_column[0]= p;
	return outdata;
}


/*
Empirically, this divider seems to work better for the header line.

Without backslashes and spaced out in perl's /x style, it would look like this:
("[^"][^"]*"        (starts with a ", has no "" in between, ends with a ".
|               or
[^%s"][^%s"]*)  anything but a "" or the user-specified delimiters. At least 1 char long.)
([%s\n]|$)      and ends with a delimiter or the end of line.
*/
static const char      headerdivider[]="(\"[^\"][^\"]*\"|[^\"%s][^\"%s]*)[%s\n]";

/*
Much of the magic below is due to the following regular expression.

Without backslashes and spaced out in perl's /x style, it would look like this:
([^%s][^%s]*)  anything but the user-specified delimiters. At least 1 char long.
[%s\n]      and ends with a delimiter or the end of line.

There are special bits to handle the case of quoted strings, outside of this basic regex.
*/
//static const char      divider[]="([^%s][^%s]*)[%s\n]";
static const char      divider[]="([^%s]*)[%s\n]";
static const char      divider2[]="([^%s][^%s]*)[%s\n][%s\n]*";

//in: the line being read, the allocated outstring, the result from the regexp search, the offset
//out: the outstring is filled with a bit of match, last_match is updated.
static void pull_string(char *line, char * outstr, regmatch_t *result, size_t * last_match){
  int     length_of_match = result[1].rm_eo - result[1].rm_so;
    memcpy(outstr, line + (*last_match)+result[1].rm_so, length_of_match);
    /*if (outstr[length_of_match -1] == '"')
        length_of_match     --;*/
    outstr[length_of_match]       = '\0';
    (*last_match)                += result[1].rm_eo+1;
}

/* Open file, find the first non-comment row, count columns, close file.  */
static int apop_count_cols_in_text(char *text_file){
  FILE * 		infile;
  char		    instr[Text_Line_Limit], outstr[Text_Line_Limit],
                full_divider[1000];
  int		    ct	                = 0,
                length_of_string    = 0;
  size_t        last_match          = 0;
  regex_t       *regex              = malloc(sizeof(regex_t));
  regmatch_t    result[3];
    sprintf(full_divider, headerdivider, apop_opts.input_delimiters, apop_opts.input_delimiters, apop_opts.input_delimiters);
    regcomp(regex, full_divider, 1);
	infile	= fopen(text_file,"r");
    apop_assert(infile, 0, 0, 'c', "Error opening file %s. Returning 0.", text_file);
	fgets(instr, Text_Line_Limit, infile);
	while(instr[0]=='#')	//burn off comment lines
		fgets(instr, Text_Line_Limit, infile);
    length_of_string= strlen(instr);
    while (last_match < length_of_string && !regexec(regex, (instr+last_match), 2, result, 0)){
        pull_string(instr,  outstr, result,  &last_match);
        if (strlen(outstr))
            ct++;
	}
	fclose(infile);
    regfree(regex);
    free(regex);
	return ct;
}

/* Open file, count lines that don't start with #, close file.  */
static int apop_count_rows_in_text(char *text_file){
  FILE * 	infile;
  char		instr[Text_Line_Limit];
  int		ct	= 0;
	infile	= fopen(text_file,"r");
    apop_assert(infile, 0, 0, 'c', "Error opening file %s. Returning 0.", text_file);
	while(fgets(instr,Text_Line_Limit,infile)!=NULL)
	    if(instr[0]!='#')	//burn off comment lines
		    ct  ++;
	fclose(infile);
	return ct-1;
}

regex_t     *strip_regex      = NULL;

static void strip_regex_alloc(){
  char      w[]           = "[:space:]\"";
  char      stripregex[1000];
    if(!strip_regex){
        sprintf(stripregex, "[%s]*\\([^%s]*.*[^%s][^%s]*\\)[%s]*", w,w,w,w,w);
        strip_regex     = malloc(sizeof(regex_t));
        regcomp(strip_regex, stripregex, 0);
    }
}

//OK, OK. C sucks. This fn just strips leading and trailing blanks.
static char * strip(char *in){
  size_t      dummy       = 0;
  char 		*out 	    = malloc(1+strlen(in));
  regmatch_t  result[3];
    if(!regexec(strip_regex, in, 2, result, 0))
        pull_string(in, out, result,  &dummy);
	return out;
}

/** \page text_format Notes on input text file formatting

If you are reading into an array or <tt>gsl_matrix</tt> or \ref
apop_data set, all text fields are taken as zeros. You will be warned
of such substitutions unless you set \code apop_opts.verbose==0\endcode
beforehand.

You will also be interested in \c apop_opts.input_delimiters. By
default, it is set to "| ,\t", meaning that a pipe,
comma, space, or tab will delimit separate entries. Try \code
strcpy(apop_opts.input_delimiters, ";")\endcode to set the delimiter to
a semicolon, for example.

There are often two delimiters in a row, e.g., "23, 32,, 12". When
it's two commas like this, the user typically means that there is a
missing value and the system should insert an NAN; when it is two tabs
in a row, this is typically just a formatting glitch. Thus, if there
are multiple delimiters in a row, Apophenia checks whether the second
(and subsequent) is a space or a tab; if it is, then it is ignored,
and if it is any other delimiter (including the end of the line) then
an NAN is inserted.

If this rule doesn't work for your situation, you can 
explicitly insert a note that there is a missing data
point. E.g., try: \code
		perl -pi.bak -e 's/,,/,NaN,/g' data_file
\endcode

If you have missing data delimiters, you will need to set \ref
apop_opts_type "apop_opts.db_nan" to a regular expression that matches
the given format. Some examples:

\code
//Apophenia's default NaN string, matching NaN, nan, or NAN:
strcpy(apop_opts.db_nan, "\\(NaN\\|nan\\|NAN\\)");
//Literal text:
strcpy(apop_opts.db_nan, "Missing");
//Literal text, capitalized or not:
strcpy(apop_opts.db_nan, "[mM]issing");
//Matches two periods. Periods are special in regexes, so they need backslashes.
strcpy(apop_opts.db_nan, "\\.\\.");
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

/** Read a delimited text file into the matrix element of an \ref
 apop_data set.

  See \ref text_format.

\param text_file  = "-"  The name of the text file to be read in. If "-" (the default), use stdin.
\param has_row_names = 0. Does the lines of data have row names?
\param has_col_names = 1. Is the top line a list of column names? If there are row names, then there should be no first entry in this line like 'row names'. That is, for a 100x100 data set with row and column names, there are 100 names in the top row, and 101 entries in each subsequent row (name plus 100 data points).
\return 	Returns an apop_data set.

<b>example:</b> See \ref apop_ols.

This function uses the \ref designated syntax for inputs.

\ingroup convertfromtext	*/
#ifdef APOP_NO_VARIADIC
apop_data * apop_text_to_data(char *text_file, int has_row_names, int has_col_names){
#else
apop_data * apop_text_to_data_base(char *text_file, int has_row_names, int has_col_names){
#endif
  apop_data     *set;
  FILE * 		infile;
  char		    instr[Text_Line_Limit], 
                *str, *stripped,
                outstr[Text_Line_Limit],
                full_divider[1000],
                header_divider[1000];
  int 		    i	        = 0,
                line_no     = 0,
                length_of_string,
		        ct, colno, rowct;
  size_t        last_match;
  regex_t       *regex      = malloc(sizeof(regex_t));
  regex_t       *headregex  = malloc(sizeof(regex_t));
  regmatch_t    result[2];
	ct	    = apop_count_cols_in_text(text_file);
	rowct	= apop_count_rows_in_text(text_file);
    set     = apop_data_alloc(0,rowct+1-has_col_names,ct);
    if (strcmp(text_file,"-"))
        infile	= fopen(text_file,"r");
    else
        infile  = stdin;
    apop_assert(infile, NULL, 0, 'c', "Error opening file %s. Returning NULL.", text_file);
    sprintf(full_divider, divider, apop_opts.input_delimiters, apop_opts.input_delimiters);
    sprintf(header_divider, headerdivider, apop_opts.input_delimiters, apop_opts.input_delimiters, apop_opts.input_delimiters);
    regcomp(regex, full_divider, 1);
    regcomp(headregex, header_divider, 1);
    strip_regex_alloc();

    //First, handle the top line, if we're told that it has column names.
    if (has_col_names){
	    fgets(instr, Text_Line_Limit, infile);
        line_no   ++;
	    while(instr[0]=='#')	//burn off comment lines
		    fgets(instr, Text_Line_Limit, infile);
	    set->names->colct= 0;
	    set->names->column	= malloc(sizeof(char*));
        last_match      = 0;
        length_of_string= strlen(instr);
        while (last_match < length_of_string && !regexec(headregex, (instr+last_match), 2, result, 0)){
            pull_string(instr,  outstr, result,  &last_match);
	        stripped	    = strip(outstr);
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
                if (strlen(outstr)){
                    colno++;
                    gsl_matrix_set(set->matrix, i-1, colno-1,	 strtod(outstr, &str));
                    if (apop_opts.verbose && !strcmp(outstr, str))
                        printf("trouble converting item %i on line %i; using zero.\n", colno, i);
                } else{
                    char d = instr[last_match-1];
                    if (d!='\t' && d!=' '){
                        gsl_matrix_set(set->matrix, i-1, colno-1, GSL_NAN);
                        colno++;
                    }
                }
            }
		}
	}
	fclose(infile);
    regfree(regex);     free(regex);
    regfree(headregex); free(headregex);
	return set;
}

#ifndef APOP_NO_VARIADIC
apop_varad_head(apop_data *, apop_text_to_data){
    char *apop_varad_var(text_file, "-")
    int apop_varad_var(has_row_names, 0)
    int apop_varad_var(has_col_names, 1)
    return apop_text_to_data_base(text_file,has_row_names, has_col_names);
}
#endif

/** See \ref apop_db_to_crosstab for the storyline; this is the complement.
 \ingroup db
 */
int apop_crosstab_to_db(apop_data *in,  char *tabname, char *row_col_name, 
						char *col_col_name, char *data_col_name){
  int		    i,j;
  apop_name   *n = in->names;
	apop_query("CREATE TABLE %s (%s , %s , %s);", tabname, 
            apop_strip_dots(row_col_name, 'd'), 
            apop_strip_dots(col_col_name, 'd'), 
            apop_strip_dots(data_col_name, 'd'));
	apop_query("begin;");
	if (in->matrix)
		for (i=0; i< n->colct; i++)
			for (j=0; j< n->rowct; j++)
				apop_query("INSERT INTO %s VALUES ('%s', '%s',%g);", tabname, 
					n->row[j], n->column[i], gsl_matrix_get(in->matrix, j, i));
	if (in->text)
		for (i=0; i< n->textct; i++)
			for (j=0; j< n->rowct; j++)
				apop_query("INSERT INTO %s VALUES ('%s', '%s','%s');", tabname, 
					n->row[j], n->text[i], in->text[j][i]);
	apop_query("commit;");
	return 0;
}

/** \page dbtomatrix Converting from database table to <tt>gsl_matrix</tt> or \ref apop_data

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
gsl_vector *apop_vector_copy(const gsl_vector *in){
    if (!in) return NULL;
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
gsl_matrix *apop_matrix_copy(const gsl_matrix *in){
  gsl_matrix *out = gsl_matrix_alloc(in->size1, in->size2);
    gsl_matrix_memcpy(out, in);
    return out;
}





/////////////////////////////
//The text processing section:
/////////////////////////////

static regex_t   *regex;
static regex_t   *headregex;
static regex_t   *nan_regex = NULL;
static int  use_names_in_file;
static char *add_this_line  = NULL;
static char **fn            = NULL;
static regex_t   *one_quote = NULL, *end_quote = NULL, *full_quote = NULL;

/**
--If the string has zero length, then it's probably a missing value.
 --If the string isn't a number, it needs quotes, and sqlite wants 0.1,
  not just .1. 
  --It may be text with no "delimiters"


 */
static char * prep_string_for_sqlite(char *astring){
  regmatch_t  result[2];
  char		*tmpstring, 
		*out	= NULL,
		*tail	= NULL;
	strtod(astring, &tail);
    tmpstring=strip(astring); 
    if (!strlen(astring)){
        asprintf(&out, "");
        return out;
    }
    if (!regexec(nan_regex, tmpstring, 1, result, 0))
        asprintf(&out, "NULL");
    else if (tail){	//then it's not a number.
		if (strlen (tmpstring)==0)
			asprintf(&out, " ");
		else        //it's text but nothing special
            asprintf(&out, "'%s'",tmpstring);
	} else {	    //sqlite wants 0.1, not .1
		assert(strlen (tmpstring)!=0);
		if (tmpstring[0]=='.')
			asprintf(&out, "0%s",tmpstring);
		else
            asprintf(&out, tmpstring);
	}
    free(tmpstring);
    return out;
}

/** Open file, find the first non-comment row, count columns, close file.  */
static int count_cols_in_row(char *instr){
  int       length_of_string    = strlen(instr);
  int       ct                  = 0;
  size_t    last_match          = 0;
  char	    outstr[Text_Line_Limit];
  regmatch_t result[3];
    while (last_match < length_of_string && !regexec(headregex, (instr+last_match), 2, result, 0)){
        pull_string(instr,  outstr, result,  &last_match);
        if (strlen(outstr))
            ct++;
	}
	return ct;
}

static int get_field_names(int has_col_names, char **field_names, FILE *infile){
  char		instr[Text_Line_Limit], *stripped, *stripme, outstr[Text_Line_Limit];
  int       i = 0, ct, length_of_string;
  size_t    last_match;
  regmatch_t  result[2];

  char full_divider2[1000];
  regex_t       *regex_no_blank_fields = malloc(sizeof(regex_t));
    //sprintf(full_divider2, headerdivider, apop_opts.input_delimiters, apop_opts.input_delimiters, apop_opts.input_delimiters, apop_opts.input_delimiters);
    sprintf(full_divider2, headerdivider, apop_opts.input_delimiters, apop_opts.input_delimiters, apop_opts.input_delimiters);
    regcomp(regex_no_blank_fields, full_divider2, 1);

    fgets(instr, Text_Line_Limit, infile);
    while(instr[0]=='#' || instr[0]=='\n')	//burn off comment lines
        fgets(instr, Text_Line_Limit, infile);
    ct              = count_cols_in_row(instr);

    if (has_col_names && field_names == NULL){
        use_names_in_file++;
        if (!has_col_names) //then you have a data line, which you should save
            asprintf(&add_this_line, instr);
        fn	            = malloc(ct * sizeof(char*));
        last_match      = 0;
        length_of_string= strlen(instr);
        while (last_match < length_of_string && !regexec(regex_no_blank_fields, (instr+last_match), 2, result, 0)){
            pull_string(instr,  outstr, result,  &last_match);
                stripme	    = strip(outstr);
                if (strlen(outstr)){
                    stripped    = apop_strip_dots(stripme,'d');
                    fn[i]	    = NULL;
                    asprintf(&fn[i], stripped);
                    free(stripped);
                    i++;
                }
            free(stripme);
        }
    } else	{
        if (field_names)
            fn	= field_names;
        else{
            asprintf(&add_this_line, instr); //save this line for later.
            fn	= malloc(ct * sizeof(char*));
            for (i =0; i < ct; i++){
                fn[i]	= malloc(1000);
                sprintf(fn[i], "col_%i", i);
            }
        }
    }
    return ct;
}

static void tab_create_mysql(char *tabname, int ct, int has_row_names){
  char  *q = NULL;
  int   i;
    asprintf(&q, "CREATE TABLE %s", tabname);
    for (i=0; i<ct; i++){
        if (i==0) 	{
            if (has_row_names)
                asprintf(&q, "%s (row_names varchar(100), ", q);
            else
                asprintf(&q, "%s (", q);
        } else		asprintf(&q, "%s varchar(100) , ", q);
        asprintf(&q, "%s ", q);
        asprintf(&q, "%s%s", q, fn[i]);
    }
    asprintf(&q, "%s varchar(100) );", q);
    apop_query(q);
    apop_assert_void(apop_table_exists(tabname, 0), 0, 's', "query \"%s\" failed.\n", q);
    if (use_names_in_file){
        for (i=0; i<ct; i++)
            free(fn[i]);
        free(fn);
        fn  = NULL;
    }
}

static void tab_create(char *tabname, int ct, int has_row_names){
  char  *r, *q = NULL;
  int   i;
    asprintf(&q, "create table %s", tabname);
    for (i=0; i<ct; i++){
        r = q;
        if (i==0) 	{
            if (has_row_names)
                asprintf(&q, "%s (row_names, ", q);
            else
                asprintf(&q, "%s (", q);
        } else		asprintf(&q, "%s' , ", q);
        free(r); r=q;
        asprintf(&q, "%s '%s", q, fn[i]);
        free(r);
    }
    asprintf(&q, "%s' ); begin;", q);
    apop_query(q);
    apop_assert_void(apop_table_exists(tabname, 0), 0, 's', "query \"%s\" failed.\n", q);
    if (use_names_in_file){
        for (i=0; i<ct; i++)
            free(fn[i]);
        free(fn);
        fn  = NULL;
    }
}

static void line_to_insert(char instr[], char *tabname){
  int       pc   = 0, //pending comma.
            length_of_string= strlen(instr);
  size_t    last_match      = 0;
  char	    outstr[Text_Line_Limit], *prepped;
  char      *r, *q  = malloc(100+ strlen(tabname));
  regmatch_t  result[2], other_result[2];
    sprintf(q, "INSERT INTO %s VALUES (", tabname);
    while (last_match < length_of_string && !regexec(regex, instr+last_match, 2, result, 0)){
        pull_string(instr,  outstr, result,  &last_match);
        //We need to check for , "fields with, uh, delimiters",
        if (!regexec(one_quote,outstr, 1, result, 0) 
                && regexec(full_quote,outstr, 1, result, 0)){
            r = q;
            asprintf(&q, "%s%s%s%c", r, (pc ? ", ": ""), outstr, instr[last_match-1]);
            free(r);
            regexec(end_quote, instr+last_match, 2, other_result, 0);
            pull_string(instr,  outstr, other_result,  &last_match);
            r = q;
            asprintf(&q, "%s%s", r, outstr);
            free(r);
            pc =1;
        } else {
            prepped	= prep_string_for_sqlite(outstr);
            if (strlen(prepped) > 0 && !(strlen(outstr) < 2 && (outstr[0]=='\n' || outstr[0]=='\r'))){
                r = q;
                asprintf(&q, "%s%s%s", r, (pc ? ", ": ""),  prepped);
                free(r);
                pc = 1;
            } else {
                char d = instr[last_match-1];
                if (d!='\t' && d!='\r' && d!='\n' && d!=' '){
                    r = q;
                    asprintf(&q, "%s%sNULL", (pc ? ", ": ""), r);
                    free(r);
                    pc = 1;
                }
            }
            free(prepped);
        }
    }
    apop_query("%s);",q); 
    free (q);
}

/** Read a text file into a database table.

  See \ref text_format.

Using the data set from the example on the \ref apop_ols page, here's another way to do the regression:

\code
#include <apop.h>

int main(void){ 
  apop_data       *data; 
  apop_model   *est;
    apop_db_open(NULL);
    apop_text_to_db("data", "d", 0,1,NULL);
    data       = apop_query_to_data("select * from d");
    est        = apop_estimate(data, apop_ols);
    printf("The OLS coefficients:\n");
    apop_model_print(est);
} 
\endcode

By the way, there is a begin/commit wrapper that bundles the process into bundles of 2000 inserts per transaction. if you want to change this to more or less frequent commits, you'll need to modify and recompile the code.

\param text_file    The name of the text file to be read in. If \c "-", then read from \c STDIN. (default = "-")
\param tabname      The name to give the table in the database (default = <tt> apop_strip_dots (text_file, 'd')</tt>)
\param has_row_names Does the lines of data have row names? (default = 0)
\param has_col_names Is the top line a list of column names? All dots in the column names are converted to underscores, by the way. (default = 1)
\param field_names The list of field names, which will be the columns for the table. If <tt>has_col_names==1</tt>, read the names from the file (and just set this to <tt>NULL</tt>). If has_col_names == 1 && field_names !=NULL, I'll use the field names.  (default = NULL)

\return Returns the number of rows.
This function uses the \ref designated syntax for inputs.
\ingroup convertfromtext
*/
#ifdef APOP_NO_VARIADIC
int apop_text_to_db(char *text_file, char *tabname, int has_row_names, int has_col_names, char **field_names){
#else
apop_varad_head(int, apop_text_to_db){
    char *apop_varad_var(text_file, "-")
    char *apop_varad_var(tabname, apop_strip_dots(text_file, 'd'))
    int apop_varad_var(has_row_names, 0)
    int apop_varad_var(has_col_names, 1)
    char ** apop_varad_var(field_names, NULL)
    return apop_text_to_db_base(text_file,tabname, has_row_names, has_col_names, field_names);
}

int apop_text_to_db_base(char *text_file, char *tabname, int has_row_names, int has_col_names, char **field_names){
#endif
  int       batch_size  = 2000,
      		ct, 
		    rows    = 0;
  FILE * 	infile;
  char		*q  = NULL, instr[Text_Line_Limit];
  char		full_divider[1000], head_divider[1000], nan_string[500];

    if (!one_quote){ //initialize a few global regexes.
        one_quote   = malloc(sizeof(regex_t));
        end_quote   = malloc(sizeof(regex_t));
        full_quote  = malloc(sizeof(regex_t));
        regcomp(one_quote, "\"[^\"]*", 0);
        regcomp(end_quote, "\\([^\\\"]*\"\\)", 0);
        regcomp(full_quote, "\"\\([^\"]*\\)\"", 0);
        //regcomp(one_quote, "(\\.|[^\\\"])*\"(\\.|[^\\\"])*", 0);
        //regcomp(end_quote, "(\\.|[^\\\"])*\"", 0);
    }

    strip_regex_alloc();
	use_names_in_file   = 0;    //file-global.
    if (!regex){
        regex               = malloc(sizeof(regex_t));//file-global, above.
        sprintf(full_divider, divider, apop_opts.input_delimiters, apop_opts.input_delimiters);
        regcomp(regex, full_divider, 1);
        headregex           = malloc(sizeof(regex_t));//file-global, above.
        sprintf(head_divider, headerdivider, apop_opts.input_delimiters, apop_opts.input_delimiters, apop_opts.input_delimiters);
        regcomp(headregex, head_divider, 1);
    }

	if (apop_table_exists(tabname,0)){
	       	printf("apop: %s table exists; not recreating it.\n", tabname);
		return 0; //to do: return the length of the table.
	} else{
        //divider regex:
        //NaN regex:
        if (strlen(apop_opts.db_nan)){
            //sprintf(nan_string, "\\\"*%s\\\"*", apop_opts.db_nan);
            //sprintf(nan_string, "^%s$", apop_opts.db_nan);
            sprintf(nan_string, "^%s$", apop_opts.db_nan);
            nan_regex   = malloc(sizeof(regex_t));
            regcomp(nan_regex, nan_string, 0);
        }
        if (strcmp(text_file,"-"))
		    infile	= fopen(text_file,"r");
        else
            infile  = stdin;
        apop_assert(infile, 0,  0, 'c', "Trouble opening %s. Bailing.\n", text_file);
        ct  = get_field_names(has_col_names, field_names, infile);
        if (apop_opts.db_engine=='m')
            tab_create_mysql(tabname, ct, has_row_names);
        else
            tab_create(tabname, ct, has_row_names);
        //convert a data line into SQL: insert into TAB values (0.3, 7, "et cetera");
        ct  = 0;
        if (add_this_line){
            rows    ++;
            line_to_insert(add_this_line, tabname);
        }
		while(fgets(instr,Text_Line_Limit,infile)!=NULL){
			if((instr[0]!='#') && (instr[0]!='\n')) {	//comments and blank lines.
				rows	        ++;
                line_to_insert(instr, tabname);
                if (apop_opts.db_engine != 'm' && !(ct++ % batch_size)) apop_query("commit; begin;");
			}
		}
		if (apop_opts.db_engine != 'm') apop_query("commit;");
        if (strcmp(text_file,"-"))
		    fclose(infile);
        if (nan_regex){
            regfree(nan_regex);
            free(nan_regex);
            nan_regex   = NULL;
        }
        free(q);
		return rows;
	}
}



/** This is the complement to \c apop_data_pack, qv. It converts the \c gsl_vector produced by that function back
    to an \c apop_data set with the given dimensions. It overwrites the data in the vector and matrix elements (and that's it).

 \param in a \c gsl_vector of the form produced by \c apop_data_pack.
\param d   that data set to be filled. Must be allocated to the correct size.
\ingroup conversions
*/
void apop_data_unpack(const gsl_vector *in, apop_data *d){
  int           i, offset   = 0;
  gsl_vector    vin, vout;
    if(d->vector){
        vin = gsl_vector_subvector((gsl_vector *)in, 0, d->vector->size).vector;
        gsl_vector_memcpy(d->vector, &vin);
        offset  += d->vector->size;
    }
    if(d->matrix)
        for (i=0; i< d->matrix->size1; i++){
            vin     = gsl_vector_subvector((gsl_vector *)in, offset, d->matrix->size2).vector;
            vout    = gsl_matrix_row(d->matrix, i).vector;
            gsl_vector_memcpy(&vout, &vin);
            offset  += d->matrix->size2;
        }
}

/** Sometimes, you need to turn an \c apop_data set into a column of
 numbers. Thus, this function, that takes in an apop_data set and outputs a \c gsl_vector.
 It is valid to use the \c out_vector->data element as an array of \c doubles of size \c out_vector->data->size.

 The complement is \c apop_data_unpack. I.e., 
\code
apop_data_unpack(apop_data_pack(in_data, vsize, m1size, m2size)) 
\endcode
will return the original data set (stripped of text and names).

 \param in an \c apop_data set.
 \return A \c gsl_vector with the vector data (if any), then each row of data (if any).
\ingroup conversions
 */
gsl_vector * apop_data_pack(const apop_data *in){
  int total_size    = (in->vector ? in->vector->size : 0)
                       + (in->matrix ? in->matrix->size1 * in->matrix->size2 : 0);
    if (!total_size)
        return NULL;
  int i, offset        = 0;
  gsl_vector *out   = gsl_vector_alloc(total_size);
  gsl_vector vout, vin;
    if (in->vector){
        vout     = gsl_vector_subvector((gsl_vector *)out, 0, in->vector->size).vector;
        gsl_vector_memcpy(&vout, in->vector);
        offset  += in->vector->size;
    }
    if (in->matrix)
        for (i=0; i< in->matrix->size1; i++){
            vin = gsl_matrix_row(in->matrix, i).vector;
            vout= gsl_vector_subvector((gsl_vector *)out, offset, in->matrix->size2).vector;
            gsl_vector_memcpy(&vout, &vin);
            offset  += in->matrix->size2;
        }
    return out;
}


#include <stdarg.h>

/** Fill a pre-allocated data set with values.

  For example:

\code
#include <apop.h>

int main(){
  apop_data *a =apop_data_alloc(2,2,2);
  double    eight   = 8.0;
    apop_data_fill(a, 8.,    2.0, eight/2,
                      0.,    6.0, eight);
    apop_data_show(a);
}
\endcode

This function has two important and subtle caveats, which are just inevitable facts
of C's handling of variadic functions.

* You must have exactly as many values as spaces in the data set. There
is no partial filling, though see \ref apop_matrix_fill and \ref apop_vector_fill.
Too many values will be ignored; too few will segfault.

* Every value must be floating point. Int values will cauase a segfault
or erratic results. So use \c 0. (a double) not \c 0 (an int). 

\param in   An \c apop_data set (that you have already allocated).
\param ...  A series of exactly as many floating-point values as there are blanks in the data set.
\return     A pointer to the same matrix that was input.
*/
apop_data *apop_data_fill(apop_data *in, ...){
    if (!in) 
        return NULL;
  va_list  ap;
  int       i, j, start=0, fin=0, height=0;
    if (in->vector){
        start   = -1;
        height  = in->vector->size;
    }
    if (in->matrix){
        fin   = in->matrix->size2;
        height  = in->matrix->size1;
    }
    va_start(ap, in);
    for (i=0; i< height; i++)
        for (j=start; j< fin; j++)
            apop_data_set(in, i, j, va_arg(ap, double));
    return in;
}

/** The version of \c apop_vector_fill that takes in a <tt>va_list</tt>.
  If this doesn't make sense to you, then you should be using \c apop_vector_fill.
  */
gsl_vector *apop_vector_vfill(gsl_vector *in, va_list ap){
  int i;
    for (i=0; i< in->size; i++)
        gsl_vector_set(in, i, va_arg(ap, double));
    return in;
}

/** Fill a pre-allocated \c gsl_vector with values.

  See \c apop_data_alloc for a relevant example. See also \c apop_matrix_alloc.

This function has two important caveats, which are just inevitable facts
of C's handling of variadic functions.

* You must have exactly as many values as spaces in the data set.
 Too many values will be ignored; too few will segfault.

* Every value must be floating point. Int values will cauase a segfault
or erratic results.

\param in   A \c gsl_vector (that you have already allocated).
\param ...  A series of exactly as many floating-point values as there are blanks in the vector.
\return     A pointer to the same vector that was input.
*/
gsl_vector *apop_vector_fill(gsl_vector *in, ...){
    if (!in) 
        return NULL;
  va_list  ap;
  int       i;
    va_start(ap, in);
    for (i=0; i< in->size; i++)
        gsl_vector_set(in, i, va_arg(ap, double));
    va_end(ap);
    return in;
}

/** Fill a pre-allocated \c gsl_matrix with values.

  See \c apop_data_alloc for a relevant example. See also \c apop_vector_alloc.

This function has two important caveats, which are just inevitable facts
of C's handling of variadic functions.

* You must have exactly as many values as spaces in the data set. 
Too many values will be ignored; too few will segfault.

* Every value must be floating point. Int values will cauase a segfault
or erratic results.

\param in   A \c gsl_matrix (that you have already allocated).
\param ...  A series of exactly as many floating-point values as there are blanks in the matrix.
\return     A pointer to the same matrix that was input.
*/
gsl_matrix *apop_matrix_fill(gsl_matrix *in, ...){
    if (!in) 
        return NULL;
  va_list  ap;
  int       i, j;
    va_start(ap, in);
    for (i=0; i< in->size1; i++)
        for (j=0; j< in->size2; j++)
            gsl_matrix_set(in, i, j, va_arg(ap, double));
    return in;
}
