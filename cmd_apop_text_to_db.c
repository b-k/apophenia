/** \file cmd_apop_text_to_db.c	A command line script to read a text file into a database.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "db.h"
#include "conversions.h"
#include <unistd.h>


int main(int argc, char **argv){
char		c, 
		msg[1000];
int     colnames            = 1,
        rownames            = 0,
        tab_exists_check    = 0;

	sprintf(msg, "%s [-d delimiters] text_file table_name dbname\n\
e.g.: %s -d\",|\" infile.txt a_table info.db\n\
If the input text file name is a single dash, -, then read from STDIN.\n\
-nc\t\tData does not include column names\n\
-n regex\t\tCase-insensitive regular expression indicating Null values. Default: NaN \n\
-m\t\tUse a mysql database (default: SQLite)\n\
-u\t\tmysql username\n\
-p\t\tmysql password\n\
-r\t\tData includes row names\n\
-v\t\tVerbose\n\
-O\t\tIf table exists, erase it and write from scratch (i.e., Overwrite)\n\
-h\t\tPrint this help\n\
\n", argv[0], argv[0]); 

	if(argc<3){
		printf("%s", msg);
		return 0;
	}
	while ((c = getopt (argc, argv, "n:d:hmp:ru:vO")) != -1){
		switch (c){
		  case 'n':
              if (optarg[0]=='c')
			    colnames    --;
              else
                strcpy(apop_opts.db_nan, optarg);
			break;
		  case 'd':
			strcpy(apop_opts.input_delimiters, optarg);
			break;
		  case 'h':
			printf("%s", msg);
			return 0;
		  case 'm':
			apop_opts.db_engine = 'm';
            break;
		  case 'u':
			strcpy(apop_opts.db_user, optarg);
			break;
		  case 'p':
			strcpy(apop_opts.db_pass, optarg);
			break;
		  case 'r':
			rownames    ++;
			break;
		  case 'v':
			apop_opts.verbose ++;
			break;
		  case 'O':
            tab_exists_check    ++;
			break;
		}
	}
	apop_db_open(argv[optind + 2]);
    if (tab_exists_check)
        apop_table_exists(argv[optind+1],1);
	apop_text_to_db(argv[optind], argv[optind+1], rownames,colnames, NULL);
}
