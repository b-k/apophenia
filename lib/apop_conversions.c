/** \file apop_conversions.c	The various functions to convert from one format to another.

 Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
 */
#include "conversions.h"
#include "assert.h"

#define Text_Size_Limit 1000000

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

\param in 	A vector.
\param out 	A <tt>gsl_vector</tt>. Declare but do not allocate.
\param size 	You will have to tell the function how long <tt>in</tt> is.
\ingroup convertfromarray 
*/ 
void apop_array_to_vector(double *in, gsl_vector **out, int size){
int		i;
	*out	= gsl_vector_alloc(size);
	for(i=0; i < size; i++)
		gsl_vector_set(*out, i, in[i]);
}

void convert_array_to_line(double **in, double **out, int rows, int cols){
	//go from in[i][j] form to the GSL's preferred out[i*cols + j] form
int		i, j;
	*out	= malloc(sizeof(double) * rows * cols);
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			(*out)[i * cols + j]	= in[i][j];
}

/** convert a <tt>double **</tt> array to a <tt>gsl_matrix</tt>

\param in	the array to read in
\param out	the <tt>gsl_matrix</tt> out. Declare it but don't allocate it.
\param rows, cols	the size of the array.
\ingroup convertfromarray 
*/
void apop_array_to_matrix(double **in, gsl_matrix **out, int rows, int cols){
gsl_matrix_view	m;
double		*line;
	*out	= gsl_matrix_alloc(rows, cols);
	convert_array_to_line(in, &line, rows, cols);
	m	= gsl_matrix_view_array(line, rows,cols);
	gsl_matrix_memcpy(*out,&(m.matrix));
	free(line);
}

static int apop_find_index(gsl_matrix *d, double r, int start_from){
//used for apop_db_to_crosstab.
int	i	= start_from;	//i is probably the same or i+1.
	do {
		if(gsl_matrix_get(d, i,0) == r) 
			return i;
		i	++;
		i	%= d->size1;	//loop around as necessary.
		//if (i == d->size1)
			//i = 0;
	} while(i!=start_from); 
	printf(" apop %s, %i: something went wrong in the crosstabbing; couldn't find %g.\n", __FILE__, __LINE__, r);
	return 0;
}

/**Give the name of a table in the database, and names of three of its
columns: the x-dimension, the y-dimension, and the data.
the output is a 2D matrix with rows indexed by r1 and cols by
r2. if !=NULL, d1 and d2 will list the labels on the dimensions.

\ingroup db
*/
gsl_matrix * apop_db_to_crosstab(char *tabname, char *r1, char *r2, char *datacol, gsl_vector **d1, gsl_vector **d2){

gsl_matrix	*pre_d1	= NULL, 
		*pre_d2	= NULL, *datatab, *out;
int		i	= 0,
		j	= 0,
		k; 
double		datum, r, c;
gsl_vector_view	v;
	pre_d1	= apop_query_to_matrix("select distinct %s, 1 from %s order by %s", r1, tabname, r1);
	if (pre_d1 == NULL) 
		printf (" apop %s, %i: selecting %s from %s returned an empty table.\n", __FILE__, __LINE__, r1, tabname);
	pre_d2	= apop_query_to_matrix("select distinct %s from %s order by %s", r2, tabname, r2);
	if (pre_d2 == NULL) 
		printf (" apop %s, %i: selecting %s from %s returned an empty table.\n", __FILE__, __LINE__, r2, tabname);
	datatab	= apop_query_to_matrix("select %s, %s, %s from %s", r1, r2, datacol, tabname);
	out	= gsl_matrix_calloc(pre_d1->size1, pre_d2->size1);
	for (k =0; k< datatab->size1; k++){
		r	= gsl_matrix_get(datatab, k, 0);
		c	= gsl_matrix_get(datatab, k, 1);
		datum	= gsl_matrix_get(datatab, k, 2);
		i	= apop_find_index(pre_d1, r, i);
		j	= apop_find_index(pre_d2, c, j);
		gsl_matrix_set(out, i, j, datum);
	}
	if(d1!=NULL && d2!= NULL){
		*d1	= gsl_vector_alloc(pre_d1->size1);
		*d2	= gsl_vector_alloc(pre_d2->size1);
		v	= gsl_matrix_column(pre_d1, 0);
		gsl_vector_memcpy(*d1, &(v.vector));
		v	= gsl_matrix_column(pre_d2, 0);
		gsl_vector_memcpy(*d2, &(v.vector));
	}
	gsl_matrix_free(pre_d1); gsl_matrix_free(pre_d2); gsl_matrix_free(datatab);
	return out;
}

/*
Much of the magic below is due to the following regular expression.

Without backslashes and spaced out in perl's /x style, it would look like this:
("[^"]*"        (starts with a ", has no "" in between, ends with a ".
|               or
[^%s"][^%s"]*)  anything but a "" or the user-specified delimiters. At least 1 char long.)
([%s\n]|$)      and ends with a delimiter or the end of line.
*/
const char      divider[]="\\(\"[^\"]*\"\\|[^\"%s][^\"%s]*\\)[%s\n]";

//in: the line being read, the allocated outstring, the result from the regexp search, the offset
//out: the outstring is filled with a bit of match, last_match is updated.
static void pull_string(char *line, char * outstr, regmatch_t *result, int * last_match){
int     length_of_match = result[1].rm_eo - result[1].rm_so;
    memcpy(outstr, line + (*last_match)+result[1].rm_so, length_of_match);
    outstr[length_of_match]       = '\0';
    (*last_match)                += result[1].rm_eo+1;
}

/** Open file, find the first non-comment row, count columns, close file.
 */
int apop_count_cols_in_text(char *text_file){
FILE * 		infile;
char		instr[Text_Size_Limit], outstr[Text_Size_Limit],
            full_divider[1000];
int		    ct	                = 0,
            length_of_string    = 0,
            last_match          = 0;
regex_t     *regex              = malloc(sizeof(regex_t));
regmatch_t  result[3];
    sprintf(full_divider, divider, apop_opts.input_delimiters, apop_opts.input_delimiters, apop_opts.input_delimiters);
    regcomp(regex, full_divider, 0);
	infile	= fopen(text_file,"r");
	if (infile==NULL) {printf("Error opening %s", text_file); return 1;}
	fgets(instr, Text_Size_Limit, infile);
	while(instr[0]=='#')	//burn off comment lines
		fgets(instr, Text_Size_Limit, infile);
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
char		instr[Text_Size_Limit];
int		ct	= 0;
	infile	= fopen(text_file,"r");
	if (infile==NULL) {printf("Error opening %s", text_file); return 1;}
	while(fgets(instr,Text_Size_Limit,infile)!=NULL)
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

Since you're reading into an array, all text fields are taken
as zeros. You will be warned of this unless you set \ref apop_opts.verbose<tt>==0</tt> beforehand.

You will also be interested in \ref apop_opts.input_delimiters. By
default, it is set to "| ,\t", meaning that a pipe, comma, space, or tab will
delimit separate entries.

There are often two delimiters in a row, e.g., "23, 32,, 12". When
it's two commas like this, the user typically means that there is a
missing value; when it is two tabs in a row, this is typically just
formatting. Apophenia is not smart enough to work out the difference:
a sequence of delimiters is taken to be a single delimiter. If you do
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
char		instr[Text_Size_Limit], 
            *str, *stripped, 
            outstr[Text_Size_Limit],
            full_divider[1000];
int 		i	        = 0,
            line_no     = 0,
            last_match,
            length_of_string,
		    ct, colno, rowct;
regex_t     *regex  = malloc(sizeof(regex_t));
regmatch_t  result[2];
	ct	    = apop_count_cols_in_text(text_file);
	rowct	= apop_count_rows_in_text(text_file);
    set     = apop_data_alloc(rowct+1-has_col_names,ct);
	infile	= fopen(text_file,"r");
    sprintf(full_divider, divider, apop_opts.input_delimiters, apop_opts.input_delimiters, apop_opts.input_delimiters);
    regcomp(regex, full_divider, 0);

    //First, handle the top line, which is assumed to be column names.
    if (has_col_names){
	    fgets(instr, Text_Size_Limit, infile);
        line_no   ++;
	    while(instr[0]=='#')	//burn off comment lines
		    fgets(instr, Text_Size_Limit, infile);
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
	while(fgets(instr,Text_Size_Limit,infile)!=NULL){
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
				gsl_matrix_set(set->data, i-1, colno-1,	 strtod(outstr, &str));
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

/** Read a text file into a database table.

  See \ref text_format.

\param text_file    The name of the text file to be read in.
\param tabname      The name to give the table in the database
\param has_row_names Does the lines of data have row names?
\param has_col_names Is the top line a list of column names? 
\param field_names The list of field names, which will be the columns for the table. If <tt>has_col_names==1</tt>, read the names from the file (and just set this to <tt>NULL</tt>).

\return Returns the number of rows.

Using the data set from the example on the \ref apop_OLS "apop_OLS" page, here's another way to do the regression:

\code
#include <apophenia/headers.h>

int main(void){ 
apop_data       *data; 
apop_estimate   *est;
    apop_open_db(NULL);
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
char		q[20000], instr[Text_Size_Limit], **fn, *prepped;
char		*stripped, outstr[Text_Size_Limit],
            full_divider[1000];
int 		ct, one_in,
            last_match,
            length_of_string,
		    i			    = 0, 
		    use_names_in_file = 0,
		    rows			= 0;
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
			printf("apop, %s: %i. Trouble opening %s.\n", __FILE__, __LINE__, text_file);
			return 0;
		}
		if (field_names == NULL){
			use_names_in_file++;
			fgets(instr, Text_Size_Limit, infile);
			while(instr[0]=='#')	//burn off comment lines
				fgets(instr, Text_Size_Limit, infile);
			fn	= malloc(ct * sizeof(char*));

            last_match      = 0;
            length_of_string= strlen(instr);
            while (last_match < length_of_string && !regexec(regex, (instr+last_match), 2, result, 0)){
                pull_string(instr,  outstr, result,  &last_match);
	            stripped	= strip(outstr);
			    fn[i]	= malloc(1000 * sizeof(char));
			    strcpy(fn[i], stripped);
		        free(stripped);
			    i++;
	        }
		} else	fn	= field_names;
		strcpy(q, "begin; CREATE TABLE ");
		strcat(q, tabname);
		for (i=0; i<ct; i++){
			if (i==0) 	{
                if (has_row_names)
                    strcat(q, " (row_names, ");
                else
                    strcat(q, " (");
            } else		strcat(q, " , ");
			sprintf(q, "%s %s", q, fn[i]);
		}
		strcat(q, "); commit; begin;");
		apop_query_db(q);



		while(fgets(instr,Text_Size_Limit,infile)!=NULL){
			if((instr[0]!='#') && (instr[0]!='\n')) {	//comments and blank lines.
				rows	++;
				//i       ++;
				one_in  =0;
				sprintf(q, "INSERT INTO %s VALUES (", tabname);
                last_match      = 0;
                length_of_string= strlen(instr);
                while (last_match < length_of_string && !regexec(regex, (instr+last_match), 2, result, 0)){
					if(one_in++) 	strcat(q, ", ");
                    pull_string(instr,  outstr, result,  &last_match);
					prepped	=prep_string_for_sqlite(outstr);
                    if (strlen(prepped) > 0)
					    strcat(q, prepped);
					free(prepped);
				}
				apop_query_db("%s);",q);
			}
		}
		apop_query_db("commit;");
		fclose(infile);
		if (use_names_in_file){
			free(fn);
		}
		return rows;
	}
}

/** See \ref apop_db_to_crosstab for the storyline; this is the complement.
 \ingroup db
 */
int apop_crosstab_to_db(gsl_matrix *in, apop_name n, char *tabname, char *row_col_name, 
						char *col_col_name, char *data_col_name){
int		i,j;
	apop_query_db("CREATE TABLE %s (%s , %s , %s);", tabname, row_col_name, col_col_name, data_col_name);
	for (i=0; i< n.colnamect; i++){
		apop_query_db("begin;");
		for (j=0; j< n.rownamect; j++)
			apop_query_db("INSERT INTO %s VALUES (%s, %s,%g);", tabname, 
					n.rownames[j], n.colnames[i], gsl_matrix_get(in, j, i));
		apop_query_db("commit;");
	}
	return 0;
}

/** \page dbtomatrix converting from database table to <tt>gsl_matrix</tt> or \ref apop_data

Use <tt>fill_me = apop_query_to_matrix("select * from table_name;");</tt>
or <tt>fill_me = apop_query_to_data("select * from table_name;");</tt>. [See \ref apop_query_to_matrix; \ref apop_query_to_data.]
\ingroup convertfromdb
*/

