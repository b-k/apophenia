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

	sprintf(msg, "%s [-d delimiters] text_file table_name dbname\n\
e.g.: %s -d\",|\" infile.txt a_table info.db\n\
I ignore delimiters right now; sorry.\n", argv[0], argv[0]); 

	if(argc<3){
		printf(msg);
		return 0;
	}
	while ((c = getopt (argc, argv, "d:h")) != -1){
		switch (c){
		  case 'd':
			delimiters	= optarg;
			break;
		  case 'h':
			printf(msg);
			return 0;
		}
	}
	apop_open_db(argv[optind + 2]);
	apop_text_to_db(argv[optind], argv[optind+1], 0,0, NULL);
	return 0;
}
