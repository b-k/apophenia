/*Conversions.h

  Convenience functions to convert among vectors (gsl_vector), matrices (gsl_matrix), 
  arrays (double *), and (database) tables

  To the extent that it makes sense, all functions are of the form convert_in_to_out(in, out),
  and allocate the output object for you.
*/


#include "db.h"
#include <string.h>
#include <gsl/gsl_matrix.h>

// vector ==>array
int convert_vector_to_array(gsl_vector *in, double **out);
//Returns the length of the array (i.e., in->size);

//array ==> vector
void convert_array_to_vector(double *in, gsl_vector **out, int size);




//text ==> db
void convert_text_to_db(char *text_file, char *db, char *tabname, int ct, 
				char **field_names, char *field_types);

/* text_file: the input file. At the moment, it needs to be comma delimited.
	Lines with a # at the head are taken to be comments and ignored.
	If field_names is NULL, then the first non-comment line of
	the file is taken to be strings giving the (comma-delimited)
	field names.
   db: name of the database file on the hard drive. If NULL, then keep the 
   	database in memory.
   tabname: the name to give the table in the database
   ct: the number of fields the table will have
   field_names: the list of field names, which will be the columns for
   	the table. If NULL, read the names from the file.
   field_types: A string of ct Is, like "IIIII". Come back later and
   	I'll take other types as well.
*/
