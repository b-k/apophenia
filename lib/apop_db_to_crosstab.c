/** \file apop_db_to_crosstab.c	Command line utility to convert a three-column table to a crosstab.

 Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
 */
#include "apophenia/db.h"
#include "apophenia/linear_algebra.h"
#include "apophenia/conversions.h"
#include "apophenia/output.h"
#include <unistd.h>


int main(int argc, char **argv){
char		c, 
		*delimiter,
		msg[1000],
		*outfile	=NULL;
gsl_matrix	*m;

	sprintf(msg, "%s [opts] dbname table_name rows columns data\n\n\
-d\tdelimiter\t\tdefault= \",\"\n\
-f\tfile to dump to\t\tdefault=STDOUT\n", argv[0]); 

	if(argc<5){
		printf(msg);
		return 0;
	}
	delimiter	= malloc(sizeof(char) * 5);
	strcpy(delimiter, ",");
	while ((c = getopt (argc, argv, "d:f:h")) != -1){
		switch (c){
		  case 'd':
			  delimiter	= optarg;
			  break;
		  case 'f':
			  outfile	= optarg;
			  break;
		  case 'h':
			printf(msg);
			return 0;
		}
	}
	apop_open_db(argv[optind]);
	m	= apop_db_to_crosstab(argv[optind +1], argv[optind+2], argv[optind+3], argv[optind+4], NULL, NULL);
	apop_matrix_print(m,delimiter, outfile);
	return 0;
}
