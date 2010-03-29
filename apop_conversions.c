/** \file apop_conversions.c	The various functions to convert from one format to another. */
/* Copyright (c) 2006--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#include "internal.h"
#include "conversions.h"
#include <gsl/gsl_math.h> //GSL_NAN
#include <assert.h>

/** \defgroup conversions Conversion functions
The functions to shunt data between text files, database tables, GSL matrices, and plain old arrays.*/

/** Converts a \c gsl_vector to an array.

\param in A \c gsl_vector

\return A \c double*, which will be <tt>malloc</tt>ed inside the function.

\ingroup conversions
*/
double * apop_vector_to_array(const gsl_vector *in){
//Does not use memcpy, because we don't know the stride of the vector.
  apop_assert(in, NULL, 1, 'c', "You sent me a NULL vector; returning NULL");
  double	*out	= malloc(sizeof(double) * in->size);
	for (size_t i=0; i < in->size; i++)
		out[i]	= gsl_vector_get(in, i);
	return out;
}

/** Just copies a one-dimensional array to a <tt>gsl_vector</tt>. The input array is undisturbed.

\param in     An array of <tt>double</tt>s. (No default. Must not be \c NULL);
\param size 	How long \c line is. If this is zero or omitted, I'll
guess using the <tt>sizeof(line)/sizeof(line[0])</tt> trick, which will
work for most arrays allocated using <tt>double []</tt> and won't work
for those allocated using <tt>double *</tt>. (default = auto-guess)
\return         A <tt>gsl_vector</tt> (which I will allocate for you).
\ingroup conversions
This function uses the \ref designated syntax for inputs.

*/ 
APOP_VAR_HEAD gsl_vector * apop_array_to_vector(double *in, int size){
    double * apop_varad_var(in, NULL);
    apop_assert(in, NULL, 1, 'c', "You sent me NULL data; returning NULL.");
    int apop_varad_var(size, sizeof(in)/sizeof(in[0]));
    return apop_array_to_vector_base(in, size);
APOP_VAR_END_HEAD
    gsl_vector      *out = gsl_vector_alloc(size);
    gsl_vector_view	v	 = gsl_vector_view_array((double*)in, size);
	gsl_vector_memcpy(out,&(v.vector));
    return out;
}

/** Mathematically, a vector of size \f$N\f$ and a matrix of size \f$N \times 1 \f$ are equivalent, but they're two different types in C. This function copies the data in a vector to a new one-column (or one-row) matrix and returns the newly-allocated and filled matrix.

  For the reverse, try \ref apop_data_pack.

\param in a \c gsl_vector (No default. If \c NULL, I return \c NULL, with a warning if <tt>apop_opts.verbose >=1 </tt>)
\param row_col If \c 'r', then this will be a row (1 x N) instead of the default, a column (N x 1). (default: \c 'c')
\return a newly-allocated <tt>gsl_matrix</tt> with one column.

This function uses the \ref designated syntax for inputs.
 \ingroup conversions
 */
APOP_VAR_HEAD gsl_matrix * apop_vector_to_matrix(const gsl_vector *in, char row_col){
    const gsl_vector * apop_varad_var(in, NULL);
    apop_assert(in,  NULL, 1,'c', "Converting NULL vector to NULL matrix.");
    char apop_varad_var(row_col, 'c');
    return apop_vector_to_matrix_base(in, row_col);
APOP_VAR_ENDHEAD
    gsl_matrix *out = 
        (row_col == 'r' || row_col == 'R') 
           ? gsl_matrix_alloc(1, in->size)
           : gsl_matrix_alloc(in->size, 1);
    apop_assert(out,  NULL, 0,'s', "gsl_matrix_alloc failed; probably out of memory.");
    if (row_col == 'r' || row_col == 'R') 
        gsl_matrix_set_row(out, 0, in);
    else
        gsl_matrix_set_col(out, 0, in);
    return out;
}

static void convert_array_to_line(const double **in, double **out, const int rows, const int cols){
	//go from in[i][j] form to the GSL's preferred out[i*cols + j] form
	*out	= malloc(sizeof(double) * rows * cols);
	for (size_t i=0; i<rows; i++)
		for (size_t j=0; j<cols; j++)
			(*out)[i * cols + j]	= in[i][j];
}

/** Convert a <tt>double **</tt> array to a <tt>gsl_matrix</tt>

\param in	the array to read in
\param rows, cols	the size of the array.
\return the <tt>gsl_matrix</tt>, allocated for you and ready to use.

usage: \code gsl_matrix *m = apop_array_to_matrix(indata, 34, 4); \endcode

If you want to initialize on the allocation line, this isn't what you want. See \ref apop_line_to_matrix.
\ingroup conversions
*/
gsl_matrix * apop_array_to_matrix(const double **in, const int rows, const int cols){
  double		    *line;
    gsl_matrix  *out = gsl_matrix_alloc(rows, cols);
	convert_array_to_line(in, &line, rows, cols);
    gsl_matrix_view   m = gsl_matrix_view_array(line, rows,cols);
	gsl_matrix_memcpy(out,&(m.matrix));
	free(line);
    return out;
}

/** Convert a <tt>double **</tt> array to an \ref apop_data set. It will
have no names. Input data is copied.

\param in	the array to read in
\param rows, cols	the size of the array.
\return the \ref apop_data set, allocated for you and ready to use.

usage: \code apop_data *d = apop_array_to_data(indata, 34, 4); \endcode
\ingroup conversions

If you want to initialize on the allocation line, this isn't what you want. See \ref apop_line_to_data.
*/
apop_data * apop_array_to_data(const double **in, const int rows, const int cols){
    return apop_matrix_to_data(apop_array_to_matrix(in, rows, cols));
}

/** Convert a <tt>double *</tt> array to a <tt>gsl_matrix</tt>. Input data is copied.

\param line	the array to read in
\param rows, cols	the size of the array.
\return the <tt>gsl_matrix</tt>, allocated for you and ready to use.

usage: \code gsl_matrix *m = apop_line_to_matrix(indata, 34, 4); \endcode
\see apop_arrary_to_matrix
\ingroup conversions
*/
gsl_matrix * apop_line_to_matrix(double *line, int rows, int cols){
    gsl_matrix    *out= gsl_matrix_alloc(rows, cols);
    gsl_matrix_view	m = gsl_matrix_view_array(line, rows,cols);
	gsl_matrix_memcpy(out,&(m.matrix));
    return out;
}

/** A convenience function to convert a <tt>double *</tt> array to an \ref apop_data set. It will
have no names. The input data is copied, not pointed to.

See also \ref apop_line_to_matrix or \ref apop_array_to_vector; this function will use these and then wrap an \ref apop_data struct around the output(s).

\param in	The array to read in. If there were appropriately placed line breaks, then this would look like the eventual data set. For example,
\code
double params[] = {0, 1, 2
                   3, 4, 5};
apop_data *out = apop_line_to_data(params, 2, 2, 2);
\endcode
will produce an \ref apop_data set with a vector \f$\left[\matrix{0 \cr 3}\right]\f$ and a matrix \f$\left[\matrix{1 & 2 \cr 4 & 5}\right]\f$.
\param vsize    The vector size. If there are also rows/cols, I expect this to equal the number or rows.
\param rows, cols	the size of the array.
\return the \ref apop_data set, allocated for you and ready to use.

\ingroup conversions
*/
apop_data * apop_line_to_data(double *in, int vsize, int rows, int cols){
    if (vsize==0 && (rows>0 && cols>0))
      return apop_matrix_to_data(apop_line_to_matrix(in, rows, cols));
    if ((rows==0 || cols==0) && vsize>0)
      return apop_vector_to_data(apop_array_to_vector(in, vsize));
    apop_assert(vsize==rows,  NULL, 0,'c',"apop_line_to_data expects either only a matrix, only a vector, or that matrix row count and vector size are equal. You gave me a row size of %i and a vector size of %i. Returning NULL.\n", rows, vsize);
  int ctr = 0;
  apop_data *out  = apop_data_alloc(vsize, rows, cols);
    for (size_t i=0; i< rows; i++)
        for (int j=-1; j< cols; j++)
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
r2.

\param tabname The database table I'm querying. Anything that will work inside a \c from clause is OK, such as a subquery in parens.
\param r1 The column of the data set that will indicate the rows of the output crosstab
\param r2 The column of the data set that will indicate the columns of the output crosstab
\param datacol The column of the data set holding the data for the cells of the crosstab

\ingroup db
*/
apop_data  *apop_db_to_crosstab(char *tabname, char *r1, char *r2, char *datacol){
  gsl_matrix	*out;
  int		    i, j=0;
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
	for (size_t k =0; k< datachars->textsize[0]; k++){
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

/** See \ref apop_db_to_crosstab for the storyline; this is the complement.
 \ingroup db
 */
void apop_crosstab_to_db(apop_data *in,  char *tabname, char *row_col_name, 
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
			for (j=0; j< n->rowct; j++){
                double x = gsl_matrix_get(in->matrix, j, i); 
                if (!isnan(x))
                    apop_query("INSERT INTO %s VALUES ('%s', '%s',%g);", tabname, 
                        n->row[j], n->column[i], x);
            }
	if (in->text)
		for (i=0; i< n->textct; i++)
			for (j=0; j< n->rowct; j++)
				apop_query("INSERT INTO %s VALUES ('%s', '%s','%s');", tabname, 
					n->row[j], n->text[i], in->text[j][i]);
	apop_query("commit;");
}

/** \page dbtomatrix Converting from database table to <tt>gsl_matrix</tt> or \ref apop_data

Use <tt>fill_me = apop_query_to_matrix("select * from table_name;");</tt>
or <tt>fill_me = apop_query_to_data("select * from table_name;");</tt>. [See \ref apop_query_to_matrix; \ref apop_query_to_data.]
\ingroup conversions
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
    apop_assert(out, NULL, 0, 's', "failed to allocate a gsl_vector of size %zu. Out of memory?", in->size);
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
    if (!in) return NULL;
  gsl_matrix *out = gsl_matrix_alloc(in->size1, in->size2);
    apop_assert(out, NULL, 0, 's', "failed to allocate a gsl_matrix of size %zu x %zu. Out of memory?", in->size1, in->size2);
    gsl_matrix_memcpy(out, in);
    return out;
}


///////////////The text processing section
static int  use_names_in_file;
static char **fn            = NULL;

/*
Much of the magic below is due to the following regular expression, which breaks a line into fields.

Without backslashes and spaced out in Perl's /x style, it would look like this:
[[:space:]]*    all the spaces you can eat.
("([^"]|[\]")+" (starts with a ", has no "" in between (may have a \"), ends with a "., at least one char.
|               or
[^%s"]+)        anything but a "" or the user-specified delimiters. At least 1 char long.)
[[:space:]]*    has all the spaces you can eat,
[%s\n]          and ends with a delimiter or the end of line.
*/
static const char      divider[]="[[:space:]]*(\"([^\"]|[\\]\")+\"|[^\"%s]+)[[:space:]]*[%s\n]";

//in: the line being read, the allocated outstring, the result from the regexp search, the offset
//out: the outstring is filled with a bit of match, last_match is updated.
static void pull_string(char *line, char * outstr, regmatch_t *result,
        size_t * last_match, int prev_end, int last_end){
  int     length_of_match, match_start;
    if (last_end){
        length_of_match = last_end - prev_end;
        match_start = prev_end;
    } else {
        length_of_match = result[1].rm_eo - result[1].rm_so;
        match_start     = (*last_match)+result[1].rm_so;
        (*last_match)  += result[1].rm_eo+1;
    }
    memcpy(outstr, line + match_start, length_of_match);
    outstr[length_of_match]   = '\0';
}

//OK, OK. C sucks. This fn just strips leading and trailing blanks.
static char * strip(char *in){
  size_t      dummy       = 0;
  char 		*out 	    = malloc(1+strlen(in));
  out[0] = '\0';
  static regex_t *strip_regex     = NULL;
  if (!strip_regex){
      char    w[]     = "[:space:]\"\n";
      strip_regex     = malloc(sizeof(regex_t));
      /* pattern is: maybe blanks, (either one non-blank, or non-blank/anything/non-blank), 
                maybe more blanks.  Then retain only the part in parens. */
      char    stripregex[1000];
      sprintf(stripregex, "[%s]*([^%s]|[^%s]+.*[^%s]+)[%s]*", w,w,w,w,w);
      regcomp(strip_regex, stripregex, REG_EXTENDED);
  }
  regmatch_t  result[3];
    if(!regexec(strip_regex, in, 2, result, 0))
        pull_string(in, out, result,  &dummy, 0, 0);
	return out;
}

static char* read_a_line(FILE *infile, char *filename){
    int chunk = 1000;
    char *line = NULL, instring[chunk];
    int len_so_far = 0;
    int total_so_far, LL;
    do {
        LL = line ? strlen(line) : 0;
        if (fgets(instring, chunk, infile)){
            total_so_far = LL +strlen(instring);
            if (len_so_far < total_so_far){
                len_so_far = total_so_far+1;
                line = realloc(line, len_so_far);
            }
            memcpy(line+LL, instring, strlen(instring)+1);
        }
        else return line;
    } while (!memchr(instring, '\n', strlen(instring)));//while no \n in string.
    /*if (!line)
            apop_error(0, 's', "Trouble reading %s, such as no data lines "\
                    " or garbage in the file. Exiting.", filename);*/
    return line;
}

/** \page text_format Notes on input text file formatting

If you are reading into an array or <tt>gsl_matrix</tt> or \ref apop_data set, all text fields are taken as zeros. You will be warned of such substitutions unless you set \code apop_opts.verbose==0\endcode beforehand.

You will also be interested in \c apop_opts.input_delimiters. By default, it is set to "| ,\t", meaning that a pipe, comma, space, or tab will delimit separate entries. Try \code strcpy(apop_opts.input_delimiters, ";")\endcode to set the delimiter to a semicolon, for example.

There are often two delimiters in a row, e.g., "23, 32,, 12". When it's two commas like this, the user typically means that there is a missing value and the system should insert an NAN; when it is two tabs in a row, this is typically just a formatting glitch. Thus, if there are multiple delimiters in a row, Apophenia checks whether the second (and subsequent) is a space or a tab; if it is, then it is ignored, and if it is any other delimiter (including the end of the line) then an NAN is inserted.

If this rule doesn't work for your situation, you can explicitly insert a note that there is a missing data
point. E.g., try: \code
		perl -pi.bak -e 's/,,/,NaN,/g' data_file
\endcode

If you have missing data delimiters, you will need to set \ref apop_opts_type "apop_opts.db_nan" to a regular expression that matches the given format. Some examples:

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

The system also uses the standards for C's atof() function for
floating-point numbers: INFINITY, -INFINITY, and NaN work as expected.
I use some tricks to get SQLite to accept these values, but they work.

Text is always delimited by quotes. Delimiters inside quotes are
perfectly OK, e.g., "Males, 30-40", is an OK column name, as is "Males named \\"Joe\\"".

Lines beginning with # (i.e. in the first column) are taken to be comments and ignored.  Blank lines are also ignored.

If there are row names and column names, then the input will not be perfectly square: there should be no first entry in the row with column names like 'row names'. That is, for a 100x100 data set with row and column names, there are 100 names in the top row, and 101 entries in each subsequent row (name plus 100 data points).
*/

static int prep_text_reading(char *text_file, FILE **infile, regex_t *regex, regex_t *nan_regex){
  char		full_divider[1000], nan_string[500];
    if (strcmp(text_file,"-"))
	    *infile	= fopen(text_file,"r");
    else
        *infile  = stdin;
    apop_assert(infile, 0,  0, 'c', "Trouble opening %s. Returning NULL.\n", text_file);

    sprintf(full_divider, divider, apop_opts.input_delimiters, apop_opts.input_delimiters);
    regcomp(regex, full_divider, REG_EXTENDED);

    if (strlen(apop_opts.db_nan)){
        sprintf(nan_string, "^%s$", apop_opts.db_nan);
        regcomp(nan_regex, nan_string, 0);
    }
    return 1;
}

/**
--If the string has zero length, then it's probably a missing value.
 --If the string isn't a number, it needs quotes, and SQLite wants 0.1,
  not just .1. 
  --It may be text with no "delimiters"
 */
static char * prep_string_for_sqlite(char *astring, regex_t *nan_regex){
  regmatch_t  result[2];
  char  *out	    = NULL,
		*tail	    = NULL,
		*stripped	= strip(astring);
	if(strtod(stripped, &tail)) /*do nothing.*/;

    if (!strlen(astring)){ //it's empty
        out = malloc(1); 
        out[0] ='\0';
    } else if (!regexec(nan_regex, stripped, 1, result, 0) 
                    || !strlen (stripped)) //nan_regex match or blank field = NaN.
        asprintf(&out, "NULL");
    else if (strlen(tail)){	//then it's not a number.
#ifdef HAVE_LIBSQLITE3 
        char *out = sqlite3_mprintf("%Q", stripped);//extra checks for odd chars.
#else
        asprintf(&out, "\"%s\"",stripped);
#endif
	} else {	    //number, but sqlite wants 0.1, not .1
		assert(strlen (stripped)!=0);
        if (isinf(atof(stripped))==1)
			asprintf(&out, "9e9999999");
        else if (isinf(atof(stripped))==-1)
			asprintf(&out, "-9e9999999");
        else if (gsl_isnan(atof(stripped)))
			asprintf(&out, "0.0/0.0");
        else if (stripped[0]=='.')
			asprintf(&out, "0%s",stripped);
		else
            asprintf(&out, "%s", stripped);
	}
    free(stripped);
    return out;
}

static int count_cols_in_row(char *instr, regex_t *regex, int *field_ends){
  int       length_of_string    = strlen(instr);
  int       ct                  = 0;
  size_t    last_match          = 0;
  char	    outstr[strlen(instr+1)];
  regmatch_t result[3];
    if (field_ends){
        int ct = 0;
        int len = strlen(instr);
        while (field_ends[ct++] < len-1)//minus the newline
            ;
        return ct;
    } else while (last_match < length_of_string && !regexec(regex, (instr+last_match), 2, result, 0)){
        pull_string(instr,  outstr, result,  &last_match, 0, 0);
        if (strlen(outstr))
            ct++;
	}
	return ct;
}

static int get_field_names(int has_col_names, char **field_names, FILE
        *infile, char *filename, regex_t *regex, char **add_this_line, int *field_ends){
  char		*instr;
  int       i = 0, ct, length_of_string, prev_end, last_end = 0;
  size_t    last_match;
  regmatch_t  result[2];
    instr =  read_a_line(infile, filename);
    while(instr[0]=='#' || instr[0]=='\n')	//burn off comment lines
        instr =  read_a_line(infile, filename);
    ct   = count_cols_in_row(instr, regex, field_ends);

    if (has_col_names && field_names == NULL){
        use_names_in_file++;
        last_end = field_names ? field_ends[0] : 0;
        if (!has_col_names) //then you have a data line, which you should save
            asprintf(add_this_line, "%s", instr);
        fn	            = malloc(ct * sizeof(char*));
        last_match      = 0;
        length_of_string= strlen(instr);
        char outstr[length_of_string+1];
        while (last_match < length_of_string 
                && last_end < length_of_string -1
                && !regexec(regex, (instr+last_match), 2, result, 0)){
            prev_end = last_end;
            last_end = field_ends ? field_ends[i] : 0;
            pull_string(instr,  outstr, result,  &last_match, prev_end, last_end);
            if (strlen(outstr)){
                char *tmpstring=strip(outstr); //remove "exraneous quotes".
                fn[i]  = apop_strip_dots(tmpstring,'d');
                free(tmpstring);
                i++;
            }
        }
    } else	{
        if (field_names)
            fn	= field_names;
        else{
            asprintf(add_this_line, "%s", instr); //save this line for later.
            fn	= malloc(ct * sizeof(char*));
            for (i =0; i < ct; i++){
                fn[i]	= malloc(1000);
                sprintf(fn[i], "col_%i", i);
            }
        }
    }
    return ct;
}

/** Read a delimited text file into the matrix element of an \ref
 apop_data set.

  See \ref text_format.

\param text_file  = "-"  The name of the text file to be read in. If "-" (the default), use stdin.
\param has_row_names = 'n'. Does the lines of data have row names?
\param has_col_names = 'y'. Is the top line a list of column names? If there are row names, then there should be no first entry in this line like 'row names'. That is, for a 100x100 data set with row and column names, there are 100 names in the top row, and 101 entries in each subsequent row (name plus 100 data points).
\param field_ends If fields have a fixed size, give the end of each field, e.g. {3, 8 11}.
\return 	Returns an apop_data set.

<b>example:</b> See \ref apop_ols.

This function uses the \ref designated syntax for inputs.

\ingroup conversions	*/
APOP_VAR_HEAD apop_data * apop_text_to_data(char *text_file, int has_row_names, int has_col_names, int *field_ends){
    char *apop_varad_var(text_file, "-")
    int apop_varad_var(has_row_names, 'n')
    int apop_varad_var(has_col_names, 'y')
    if (has_row_names==1) has_row_names ='y';
    if (has_col_names==1) has_col_names ='y';
    int * apop_varad_var(field_ends, NULL);
    return apop_text_to_data_base(text_file,has_row_names, has_col_names, field_ends);
APOP_VAR_END_HEAD
  apop_data     *set = NULL;
  FILE * 		infile;
  char		    *instr, *str, *add_this_line= NULL;
  int 		    i	        = 0,
                length_of_string, colno,
                prev_end, last_end=0, 
                hasrows = (has_row_names == 'y');
  size_t        last_match;
  regmatch_t    result[2];
  regex_t       regex, nan_regex;
    if (!prep_text_reading(text_file, &infile, &regex, &nan_regex)) 
        return NULL;

    //First, handle the top line, if we're told that it has column names.
    if (has_col_names=='y'){
        int col_ct  = get_field_names(1, NULL, infile, text_file, &regex, &add_this_line, field_ends);
        set = apop_data_alloc(0,1, col_ct);
	    set->names->colct   = 0;
	    set->names->column	= malloc(sizeof(char*));
        for (int j=0; j< col_ct; j++)
            apop_name_add(set->names, fn[j], 'c');
    }

    //Now do the body. First elmt may be a row name.
	while((instr= add_this_line) || (instr= read_a_line(infile, text_file))!=NULL){
		colno	= 0;
		if(instr[0]!='#') {
            if (!set)
                set = apop_data_alloc(0,1, count_cols_in_row(instr, &regex, field_ends));
			i ++;
            set->matrix = apop_matrix_realloc(set->matrix, i, set->matrix->size2);
            last_match      = 0;
            length_of_string= strlen(instr);
            char outstr[length_of_string+1];
            regexec(&regex, (instr+last_match), 2, result, 0);   //one for the headers.
            prev_end = last_end;
            last_end = field_ends ? field_ends[colno] : 0;
			if (hasrows){
                pull_string(instr,  outstr, result,  &last_match, prev_end, last_end);
                apop_name_add(set->names, outstr, 'r');
			}
            while (last_match < length_of_string 
                    && last_end < length_of_string -1
                    && !regexec(&regex, (instr+last_match), 2, result, 0)){
                prev_end = last_end;
                last_end = field_ends ? field_ends[hasrows+colno] : 0;
                pull_string(instr,  outstr, result,  &last_match, prev_end, last_end);
                if (strlen(outstr)){
                    colno++;
                    gsl_matrix_set(set->matrix, i-1, colno-1,	 strtod(outstr, &str));
                    if (apop_opts.verbose && !strcmp(outstr, str))
                        printf("trouble converting item %i on line %i; using zero.\n", colno, i);
                } else{
                    char d = instr[last_match-1];
                    if (d!='\t' && d!=' '){
                        colno++;
                        gsl_matrix_set(set->matrix, i-1, colno-1, GSL_NAN);
                    }
                }
            }
		}
        add_this_line = NULL;
	}
    if (strcmp(text_file,"-"))
        fclose(infile);
    regfree(&regex); regfree(&nan_regex);
	return set;
}

static void tab_create_mysql(char *tabname, int ct, int has_row_names){
  char  *q = NULL;
    asprintf(&q, "CREATE TABLE %s", tabname);
    for (int i=0; i<ct; i++){
        if (i==0) 	{
            if (has_row_names)
                asprintf(&q, "%s (row_names varchar(100), ", q);
            else
                asprintf(&q, "%s (", q);
        } else		asprintf(&q, "%s varchar(100) , ", q);
        asprintf(&q, "%s %s", q, fn[i]);
    }
    asprintf(&q, "%s varchar(100) );", q);
    apop_query("%s", q);
    apop_assert_void(apop_table_exists(tabname, 0), 0, 's', "query \"%s\" failed.\n", q);
    if (use_names_in_file){
        for (int i=0; i<ct; i++)
            free(fn[i]);
        free(fn);
        fn  = NULL;
    }
}

static void tab_create(char *tabname, int ct, int has_row_names){
  char  *r, *q = NULL;
    asprintf(&q, "create table %s", tabname);
    for (int i=0; i<ct; i++){
        r = q;
        if (i==0) 	{
            if (has_row_names)
                asprintf(&q, "%s (row_names, ", q);
            else
                asprintf(&q, "%s (", q);
        } else		asprintf(&q, "%s' numeric, ", q);
        free(r); r=q;
        asprintf(&q, "%s '%s", q, fn[i]);
        free(r);
    }
    apop_query("%s' numeric); begin;", q);
    apop_assert_void(apop_table_exists(tabname, 0), 0, 's', "query \"%s' ); begin;\" failed.\n", q);
    if (use_names_in_file){
        for (int i=0; i<ct; i++)
            free(fn[i]);
        free(fn);
        fn  = NULL;
    }
}

static void line_to_insert(char instr[], char *tabname, regex_t *regex, regex_t *nan_regex, int * field_ends){
  int       length_of_string= strlen(instr),
            prev_end    = 0, last_end    = 0, ctr = 0;
  size_t    last_match  = 0;
  char	    *prepped, comma = ' ';
  char      *r, *q  = malloc(100+ strlen(tabname));
  regmatch_t  result[2];
    sprintf(q, "INSERT INTO %s VALUES (", tabname);
    char outstr[length_of_string+1];
    while (last_match < length_of_string 
            && last_end < length_of_string-1
            && !regexec(regex, instr+last_match, 2, result, 0)){
        prev_end = last_end;
        last_end = field_ends ? field_ends[ctr++] : 0;
        pull_string(instr,  outstr, result,  &last_match, prev_end, last_end);
        prepped	= prep_string_for_sqlite(outstr, nan_regex);
        if (strlen(prepped) > 0 && !(strlen(outstr) < 2 && (outstr[0]=='\n' || outstr[0]=='\r'))){
            r = q;
            asprintf(&q, "%s%c %s", r, comma,  prepped);
            free(r);
            comma = ',';
        } else {
            char d = instr[last_match-1];
            if (d!='\t' && d!='\r' && d!='\n' && d!=' '){
                r = q;
                asprintf(&q, "%s%cNULL", r, comma);
                free(r);
                comma = ',';
            }
        }
        free(prepped);
    }
    apop_query("%s);",q); 
    free (q);
}

/** Read a text file into a database table.

  See \ref text_format.

Using the data set from the example on the \ref apop_ols page, here's another way to do the regression:

\include ols.c

By the way, there is a begin/commit wrapper that bundles the process into bundles of 2000 inserts per transaction. if you want to change this to more or less frequent commits, you'll need to modify and recompile the code.

\param text_file    The name of the text file to be read in. If \c "-", then read from \c STDIN. (default = "-")
\param tabname      The name to give the table in the database (default
= <tt> apop_strip_dots (text_file, 'd')</tt>; default in Python/R interfaces="t")
\param has_row_names Does the lines of data have row names? (default = 0)
\param has_col_names Is the top line a list of column names? All dots in the column names are converted to underscores, by the way. (default = 1)
\param field_names The list of field names, which will be the columns for the table. If <tt>has_col_names==1</tt>, read the names from the file (and just set this to <tt>NULL</tt>). If has_col_names == 1 && field_names !=NULL, I'll use the field names.  (default = NULL)
\param field_ends If fields have a fixed size, give the end of each field, e.g. {3, 8 11}.

\return Returns the number of rows.
This function uses the \ref designated syntax for inputs.
\ingroup conversions
*/
APOP_VAR_HEAD int apop_text_to_db(char *text_file, char *tabname, int has_row_names, int has_col_names, char **field_names, int *field_ends){
    char *apop_varad_var(text_file, "-")
    char *apop_varad_var(tabname, apop_strip_dots(text_file, 'd'))
    int apop_varad_var(has_row_names, 0)
    int apop_varad_var(has_col_names, 1)
    int *apop_varad_var(field_ends, NULL)
    char ** apop_varad_var(field_names, NULL)
    return apop_text_to_db_base(text_file,tabname, has_row_names, has_col_names, field_names, field_ends);
APOP_VAR_END_HEAD
  int       batch_size  = 2000,
      		ct, rows    = 0;
  FILE * 	infile;
  char		*q  = NULL, *instr, *add_this_line=NULL;
  regex_t   regex, nan_regex;

    if(!prep_text_reading(text_file, &infile, &regex, &nan_regex))
        return 0;
	use_names_in_file   = 0;    //file-global.

	apop_assert(!apop_table_exists(tabname,0), 0, 0, 'c', "table %s exists; not recreating it.", tabname);
    ct  = get_field_names(has_col_names, field_names, infile, text_file, &regex, &add_this_line, field_ends);
    if (apop_opts.db_engine=='m')
        tab_create_mysql(tabname, ct, has_row_names);
    else
        tab_create(tabname, ct, has_row_names);
    //convert a data line into SQL: insert into TAB values (0.3, 7, "et cetera");
    ct  = 0;
	while((instr=add_this_line) || (instr= read_a_line(infile, text_file))!=NULL){
		if((instr[0]!='#') && (instr[0]!='\n')) {	//comments and blank lines.
			rows ++;
            line_to_insert(instr, tabname, &regex, &nan_regex, field_ends);
            if (!(ct++ % batch_size)){
                if (apop_opts.db_engine != 'm') apop_query("commit; begin;");
                if (apop_opts.verbose >= 0) {printf(".");fflush(NULL);}
            }
		}
        add_this_line = NULL;
	}
	if (apop_opts.db_engine != 'm') apop_query("commit;");
    if (strcmp(text_file,"-"))
	    fclose(infile);
    free(q);
    regfree(&regex); regfree(&nan_regex);
	return rows;
}


/** This is the complement to \c apop_data_pack, qv. It writes the \c gsl_vector produced by that function back
    to the \c apop_data set you provide. It overwrites the data in the vector and matrix elements and, if present, the \c weights (and that's it, so names or text are as before).

 \param in a \c gsl_vector of the form produced by \c apop_data_pack. No default; must not be \c NULL.
\param d   that data set to be filled. Must be allocated to the correct size. No default; must not be \c NULL.
\param all_pages If \c 'y', then follow the <tt> ->more</tt> pointer to fill subsequent
pages; else fill only the first page. Default = \c 'n'.

\li If I get to the end of the first page and have more vector to unpack, and the data to
fill has a \c more element, then I will continue into subsequent pages.

This function uses the \ref designated syntax for inputs.
\ingroup conversions
*/
APOP_VAR_HEAD void apop_data_unpack(const gsl_vector *in, apop_data *d){
    const gsl_vector * apop_varad_var(in, NULL);
    apop_data* apop_varad_var(d, NULL);
    //char apop_varad_var(all_pages, 'n');
    //apop_data_unpack_base(in, d, all_pages);
    apop_data_unpack_base(in, d);
APOP_VAR_ENDHEAD
  int           offset   = 0;
  gsl_vector    vin, vout;
    if(d->vector){
        vin = gsl_vector_subvector((gsl_vector *)in, 0, d->vector->size).vector;
        gsl_vector_memcpy(d->vector, &vin);
        offset  += d->vector->size;
    }
    if(d->matrix)
        for (size_t i=0; i< d->matrix->size1; i++){
            vin     = gsl_vector_subvector((gsl_vector *)in, offset, d->matrix->size2).vector;
            vout    = gsl_matrix_row(d->matrix, i).vector;
            gsl_vector_memcpy(&vout, &vin);
            offset  += d->matrix->size2;
        }
    if(d->weights){
        vin = gsl_vector_subvector((gsl_vector *)in, offset, d->weights->size).vector;
        gsl_vector_memcpy(d->weights, &vin);
        offset  += d->weights->size;
    }
    if (offset != in->size && d->more){
        vin = gsl_vector_subvector((gsl_vector *)in, offset, in->size - offset).vector;
        apop_data_unpack(&vin, d->more);
    }
}

static size_t sizecount(const apop_data *in, const int all_pp){ 
    if (!in)
        return 0;
    return (in->vector ? in->vector->size : 0)
             + (in->matrix ? in->matrix->size1 * in->matrix->size2 : 0)
             + (in->weights ? in->weights->size : 0)
             + (all_pp ? sizecount(in->more, all_pp) : 0);
}

/** Sometimes, you need to turn an \c apop_data set into a column of
 numbers. Thus, this function, that takes in an apop_data set and outputs a \c gsl_vector.
 It is valid to use the \c out_vector->data element as an array of \c doubles of size \c out_vector->data->size (i.e. its <tt>stride==1</tt>).

 The complement is \c apop_data_unpack. I.e., 
\code
apop_data_unpack(apop_data_pack(in_data), data_copy) 
\endcode
will return the original data set (stripped of text and names).

 \param in an \c apop_data set. No default; if \c NULL, return \c NULL.
 \param out If this is not \c NULL, then put the output here. The dimensions must match exactly. If \c NULL, then allocate a new data set. Default = \c NULL. 
  \param all_pages If \c 'y', then follow the <tt> ->more</tt> pointer to fill subsequent
pages; else fill only the first page. Default = \c 'n'.

 \return A \c gsl_vector with the vector data (if any), then each row of data (if any), then the weights (if any), then the same for subsequent pages (if any <tt>&& .all_pages=='y'</tt>). If \c out is not \c NULL, then this is \c out.

This function uses the \ref designated syntax for inputs.
\ingroup conversions
 */
APOP_VAR_HEAD gsl_vector * apop_data_pack(const apop_data *in, gsl_vector *out, char all_pages){
    const apop_data * apop_varad_var(in, NULL);
    if (!in) return NULL;
    gsl_vector * apop_varad_var(out, NULL);
    char apop_varad_var(all_pages, 'n');
    if (out) {
        int total_size    = sizecount(in, (all_pages == 'y' || all_pages == 'Y'));
        apop_assert(out->size == total_size, NULL, 0, 's', "The input data set has %i elements,"
               " but the output vector you want to fill has size %zu. Please make these sizes equal."
               , total_size, out->size);
    }
    return apop_data_pack_base(in, out, all_pages);
APOP_VAR_ENDHEAD
        int total_size    = sizecount(in, (all_pages == 'y' || all_pages == 'Y'));
    if (!total_size)
        return NULL;
  int offset        = 0;
  if (!out)
        out   = gsl_vector_alloc(total_size);
  gsl_vector vout, vin;
    if (in->vector){
        vout     = gsl_vector_subvector((gsl_vector *)out, 0, in->vector->size).vector;
        gsl_vector_memcpy(&vout, in->vector);
        offset  += in->vector->size;
    }
    if (in->matrix)
        for (size_t i=0; i< in->matrix->size1; i++){
            vin = gsl_matrix_row(in->matrix, i).vector;
            vout= gsl_vector_subvector((gsl_vector *)out, offset, in->matrix->size2).vector;
            gsl_vector_memcpy(&vout, &vin);
            offset  += in->matrix->size2;
        }
    if (in->weights){
        vout     = gsl_vector_subvector((gsl_vector *)out, offset, in->weights->size).vector;
        gsl_vector_memcpy(&vout, in->weights);
        offset  += in->weights->size;
    }
    if ((all_pages == 'y' ||all_pages =='Y') && in->more){
        vout = gsl_vector_subvector((gsl_vector *)out, offset, out->size - offset).vector;
        apop_data_pack(in->more, &vout, .all_pages='y');
    }
    return out;
}

/* This function might give users a way to look up the location of an element in a packed
   list. It seemes theoretically useful, but am not yet convinced. So, here's a
   placeholder.
APOP_VAR_HEAD apop_data_pack_order(apop_data *ref, size_t row, size_t col, char *page){

APOP_VAR_ENDHEAD

} */

/** \def apop_data_fill (in, ap)
Fill a pre-allocated data set with values.

  For example:

\code
#include <apop.h>

int main(){
  apop_data *a =apop_data_alloc(2,2,2);
  double    eight   = 8.0;
    apop_data_fill(a, 8, 2.2, eight/2,
                      0, 6.0, eight);
    apop_data_show(a);
}
\endcode

Warning: I need as many arguments as the size of the data set, and can't count them for you. Too many will be ignored; too few will produce unpredictable results, which may include padding your matrix with garbage or a simple segfault.

I assume that <tt>vector->size==matrix->size1</tt>; otherwise I just use
\c matrix->size1.

\param in   An \c apop_data set (that you have already allocated).
\param ...  A series of at least as many floating-point values as there are blanks in the data set.
\return     A pointer to the same data set that was input.
*/

apop_data *apop_data_fill_base(apop_data *in, double ap[]){
/* In conversions.h, you'll find this header, which turns all but the first input into an array of doubles of indeterminate length:
#define apop_data_fill(in, ...) apop_data_fill_base((in), (double []) {__VA_ARGS__})
*/
    if (!in) 
        return NULL;
  int  k=0, start=0, fin=0, height=0;
    if (in->vector){
        start   = -1;
        height  = in->vector->size;
    }
    if (in->matrix){
        fin   = in->matrix->size2;
        height  = in->matrix->size1;
    }
    for (int i=0; i< height; i++)
        for (int j=start; j< fin; j++)
            apop_data_set(in, i, j, ap[k++]);
    return in;
}

/** \def apop_vector_fill(in, ap)
 Fill a pre-allocated \c gsl_vector with values.

  See \c apop_data_alloc for a relevant example. See also \c apop_matrix_alloc.

Warning: I need as many arguments as the size of the vector, and can't count them for you. Too many will be ignored; too few will produce unpredictable results, which may include padding your vector with garbage or a simple segfault.


\param in   A \c gsl_vector (that you have already allocated).
\param ...  A series of exactly as many values as there are spaces in the vector.
\return     A pointer to the same vector that was input.
*/
gsl_vector *apop_vector_fill_base(gsl_vector *in, double ap[]){
    if (!in) 
        return NULL;
    for (int i=0; i< in->size; i++)
        gsl_vector_set(in, i, ap[i]);
    return in;
}

/** \def apop_matrix_fill(in, ap)
 Fill a pre-allocated \c gsl_matrix with values. 
 
The values should be in row-major order (i.e., list the entire first row, followed by the second row, et cetera). If your data is column-major, then try calling \c gsl_matrix_transpose after this function.

  See \c apop_data_alloc for a relevant example. See also \c apop_vector_alloc.

  I need as many arguments as the size of the matrix. Too many will be ignored; too few will produce unpredictable results, which may include padding your matrix with garbage or a simple segfault.

\param in   A \c gsl_matrix (that you have already allocated).
\param ...  A series of exactly as many floating-point values as there are blanks in the matrix.
\return     A pointer to the same matrix that was input.
*/
gsl_matrix *apop_matrix_fill_base(gsl_matrix *in, double ap[]){
    if (!in) 
        return NULL;
  int       k = 0;
    for (int i=0; i< in->size1; i++)
        for (int j=0; j< in->size2; j++)
            gsl_matrix_set(in, i, j, ap[k++]);
    return in;
}

apop_data_row apop_data_get_row(apop_data *d, int row_number){
    apop_data_row out = {.vector_pt = d->vector ? gsl_vector_ptr(d->vector, row_number) : NULL,
                    .text_row = d->text ? d->text[row_number] : NULL,
                    .column_names = d->names->column,
                    .textsize = d->textsize[1],
                    .weight = d->weights ? gsl_vector_ptr(d->weights, row_number) : NULL,    
                    .index = row_number};
    if (d->matrix){
        out.mrv = gsl_matrix_row(d->matrix, row_number);
        out.matrix_row = &(out.mrv.vector);
    }
    return out;
}
