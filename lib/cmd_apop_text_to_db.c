/** \file cmd_apop_text_to_db.c	A command line script to read a text file into a database.

 Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
 */
#include "apophenia/db.h"
#include "apophenia/conversions.h"
#include <unistd.h>


int main(int argc, char **argv){
char		c, 
		*delimiters	= NULL, 
		msg[1000];
int     colnames    = 0,
        rownames    = 0;

	sprintf(msg, "%s [-d delimiters] text_file table_name dbname\n\
e.g.: %s -d\",|\" infile.txt a_table info.db\n\
-c\t\tData includes column names\n\
-r\t\tData includes row names\n\
-v\t\tVerbose\n\
-h\t\tPrint this help\n\
\n", argv[0], argv[0]); 

	if(argc<3){
		printf(msg);
		return 0;
	}
	while ((c = getopt (argc, argv, "cd:hrv")) != -1){
		switch (c){
		  case 'c':
			colnames    ++;
			break;
		  case 'd':
			delimiters	= optarg;
			break;
		  case 'h':
			printf(msg);
			return 0;
		  case 'r':
			rownames    ++;
			break;
		  case 'v':
			apop_opts.verbose ++;
			break;
		}
	}
	apop_db_open(argv[optind + 2]);
	apop_text_to_db(argv[optind], argv[optind+1], rownames,colnames, NULL);
	return 0;
}
