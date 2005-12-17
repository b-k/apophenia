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
int apop_convert_vector_to_array(gsl_vector *in, double **out){
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
void apop_convert_array_to_vector(double *in, gsl_vector **out, int size){
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
void apop_convert_array_to_matrix(double **in, gsl_matrix **out, int rows, int cols){
gsl_matrix_view	m;
double		*line;
	*out	= gsl_matrix_alloc(rows, cols);
	convert_array_to_line(in, &line, rows, cols);
	m	= gsl_matrix_view_array(line, rows,cols);
	gsl_matrix_memcpy(*out,&(m.matrix));
	free(line);
}

int apop_find_index(gsl_matrix *d, double r, int start_from){
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



int apop_count_cols_in_text(char *text_file){
//Open file, find the first valid row, count columns, close file.
FILE * 		infile;
char		instr[Text_Size_Limit], *astring;
int		ct	= 0;
	infile	= fopen(text_file,"r");
	if (infile==NULL) {printf("Error opening %s", text_file); return 1;}
	fgets(instr, Text_Size_Limit, infile);
	while(instr[0]=='#')	//burn off comment lines
		fgets(instr, Text_Size_Limit, infile);
	astring	= strtok(instr,",");
	while (astring !=NULL){
		ct++;
		astring	= strtok(NULL,",");
	}
	fclose(infile);
	return ct;
}

char * strip(char *in){
//OK, OK. C sucks.
char 		*out 	= malloc(sizeof(char) * (1+strlen(in)));
int		i	= 0, 
		done	= 0,
		first_ok= 0,
		last_ok	= 0;
	while (!done)
		if (in[first_ok]==' ' || in[first_ok]=='\n' || in[first_ok]=='\t')
			first_ok	++;
		else	done		++;
	for (i=first_ok; i< strlen(in); i++){
		out[i-first_ok]	= in[i];
		if (in[i]!=' ' && in[i]!='\n' && in[i]!='\t')
			last_ok	= i-first_ok;
	}
	if (last_ok <= first_ok-1) //blank string
		sprintf(out, " ");
	out[last_ok+1]	= '\0';	//redundant for blanks, but whatever.
	return out;
}

void add_a_name(int *ct, char*** list, char *addme){
char		*stripped;
	(*ct)	++;
	*list	= realloc(*list, (*ct) * sizeof(char*));
	(*list)[(*ct)-1]	=malloc(sizeof(char) * (strlen(addme)+1));
	stripped	= strip(addme);
	strcpy((*list)[(*ct)-1],stripped);
	free(stripped);
}

/** Read a delimited text file into an array. 
\param text_file	The input file. At the moment, it needs to be
comma delimited. Lines with a # at the head are taken to be comments
and ignored. If field_names is NULL, then the first non-comment line
of the file is taken to be strings giving the (comma-delimited) field
names. 

Since you're reading into an array, all text fields are taken
as zeros. You will be warned of this unless you set \ref
apop_verbose<tt>==0</tt> beforehand.

\param delimiters A list of delimiters. "," is typical, as is ",|", for example.
\param tab 	A table, to be allocated and filled with data.
\param names 	An apop_name structure. If the data has column names as the first (noncomment) row, then set
<tt>name.colname= 1;</tt>
before running; if the first element of each row is a row name, set
<tt>name.rowname= 1;</tt>. Else, set these to zero.
\return 	Returns the number of rows.

<b>example:</b> See \ref apop_OLS.
\bug I suspect apop_convert_text_to_array doesn't work very well with delimiters besides the standard ones.  
\ingroup convertfromtext	*/
int apop_convert_text_to_array(char *text_file, char *delimiters, double ***tab, apop_name *names){
FILE * 		infile;
char		instr[Text_Size_Limit], *astring, *str;
int 		i	= 0,
		ct, colno;
	ct	= apop_count_cols_in_text(text_file)+1;
	*tab	= malloc(sizeof(double));
	infile	= fopen(text_file,"r");
	if (names != NULL){
		fgets(instr, Text_Size_Limit, infile);
		while(instr[0]=='#')	//burn off comment lines
			fgets(instr, Text_Size_Limit, infile);
		if (names->colnamect !=0){
			astring	= strtok(instr,delimiters);
			names->colnamect= 0;
			names->colnames	= malloc(sizeof(char*));
			while (astring !=NULL){
				add_a_name(&(names->colnamect), &(names->colnames), astring);
				astring	= strtok(NULL,delimiters);
			}
		}
		if (names->rownamect !=0){
			names->rownamect= 0;
			names->rownames	= malloc(sizeof(char*));
		} else 	names->rownames	= NULL;
	}
	while(fgets(instr,Text_Size_Limit,infile)!=NULL){
		colno	= 0;
		if(instr[0]!='#') {
			i	++;
			*tab	= realloc(*tab, sizeof(double*) * i);
			(*tab)[i-1]= malloc(sizeof(double)*ct);
			astring	= strtok(instr,delimiters);
				if (names !=NULL && names->rownames !=NULL){
					add_a_name(&(names->rownamect), &(names->rownames), astring);
					astring	= strtok(NULL,delimiters);
					if (astring==NULL){
						printf("row name with no data on line %i.\n", i);
						return 1;
					}
				}
			while (astring !=NULL){
				colno++;
				(*tab)[i-1][colno-1]	= strtod(astring, &str);
				if (apop_verbose && !strcmp(astring, str))
					printf("trouble converting item %i on line %i; using zero.\n", colno, i);
				astring	= strtok(NULL,delimiters);
			}
		}
	}
	fclose(infile);
	return i;
}

char * prep_string_for_sqlite(char *astring){
//If the string isn't a number, it needs quotes.
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

/** Read a textfile into a database table.

\param text_file
The input file. At the moment, it needs to be comma delimited. Lines with
a # at the head are taken to be comments and ignored. If field_names is
<tt>NULL</tt>, then the first non-comment line of the file is taken to be strings
giving the (comma-delimited) field names.

\param tabname 
The name to give the table in the database

\param field_names 
The list of field names, which will be the columns for the table. If <tt>NULL</tt>, read the names from the file.

\return Returns the number of rows.


Using the data set from the example on the \ref apop_OLS "apop_OLS" page, here's another way to do the regression:

\verbatim
#include <gsl/gsl_matrix.h>
#include <apophenia/db.h>
#include <apophenia/conversions.h>
#include <apophenia/regression.h>
#include <apophenia/linear_algebra.h> //Print_vector

int main(void){
gsl_vector      *beta;
gsl_matrix      *data;
     apop_open_db(NULL);
     apop_convert_text_to_db("data", "d", NULL);
     apop_query_to_matrix(&data, "select * from d");
     apop_OLS(data, &beta);
     printf("The OLS coefficients:\n");
     apop_print_vector(beta);
return 0;
}
\endverbatim 
\ingroup convertfromtext
*/
int apop_convert_text_to_db(char *text_file, char *tabname, char **field_names){
FILE * 		infile;
char		q[20000], instr[Text_Size_Limit], **fn, *astring, *prepped;
char		delimiters[]	=",";
int 		ct, one_in,
		i			= 0, 
		use_names_in_file	= 0,
		rows			= 0;
	ct	= apop_count_cols_in_text(text_file);
	if (apop_table_exists(tabname,0)){
	       	printf("apop: %s table exists; not recreating it.\n", tabname);
		return 0; //to do: return the length of the table.
	} else{
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
			astring	= strtok(instr,delimiters);
			while(astring !=NULL){
				fn[i]	= malloc(1000 * sizeof(char));
				strcpy(fn[i], astring);
				astring	= strtok(NULL,delimiters);
				i++;
			}
		} else	fn	= field_names;
		strcpy(q, "begin; CREATE TABLE ");
		strcat(q, tabname);
		for (i=0; i<ct; i++){
			if (i==0) 	strcat(q, " (");
			else		strcat(q, " , ");
			sprintf(q, "%s %s", q, fn[i]);
		}
		strcat(q, "); commit; begin;");
		apop_query_db(q);
		while(fgets(instr,Text_Size_Limit,infile)!=NULL){
			if((instr[0]!='#') && (instr[0]!='\n')) {	//comments and blank lines.
				rows	++;
				i++;
				one_in=0;
				sprintf(q, "INSERT INTO %s VALUES (", tabname);
				astring	= strtok(instr,delimiters);
				while(astring !=NULL){
					if(one_in++) 	strcat(q, ", ");
					prepped	=prep_string_for_sqlite(astring);
					strcat(q, prepped);
					free(prepped);
					astring	= strtok(NULL,delimiters);
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

/** an alias for \ref apop_convert_text_to_db  */
int apop_text_to_db(char *text_file, char *tabname, char **field_names){
	return apop_convert_text_to_db(text_file, tabname, field_names);
}

/** See \ref apop_db_to_crosstab for the storyline; this is the
 * complement.

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

/** \page dbtomatrix converting from database table to gsl_matrix

Use <tt>fill_me = apop_query_to_matrix("select * from table_name;");</tt>. [See \ref apop_query_to_matrix.]
\ingroup convertfromdb
*/

