/** \file 
 A command line script to read a text file into a database.

Copyright (c) 2006--2007, 2013 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
#include <unistd.h>

int *break_down(char *in){
    int *out = NULL;
    int ctr = 0;
    char *cp = strtok (in, ",");
    while (cp != NULL) {
      out = realloc(out, sizeof(int)*(ctr+1));
      out[ctr++] = atoi(cp);
      cp = strtok (NULL, ",");
    }
    return out;
}

int main(int argc, char **argv){
    char c, *msg;
    int colnames = 'y',
        rownames = 0,
        tab_exists_check = 0;
    char **field_names = NULL;

	Asprintf(&msg, "%s [-d delimiters] text_file table_name dbname\n"
                "e.g.: %s -d\",|\" infile.txt a_table info.db\n"
"If the input text file name is a single dash, -, then read from STDIN.\n"
"Input must be plain ASCII or UTF-8.\n"
"-d\t\tThe single-character delimiters to use, e.g., -d \" ,\" or -d \"\\t\" (which you \n"
"\t\t\twill almost certainly have to write as -d \"\\\\t\"). Default: \"|,\\t\", meaning \n"
"\t\t\tthat any of a pipe, comma, or tab will delimit separate entries\n"
"-nc\t\tData does not include column names\n"
"-n regex\t\tCase-insensitive regular expression indicating Null values. Default: NaN \n"
"-m\t\tUse a mysql database (default: SQLite)\n"
"-f\t\tfixed width field ends: -f\"3,8,12,17\" (first char is one, not zero)\n"
"-u\t\tmysql username\n"
"-p\t\tmysql password\n"
"-r\t\tData includes row names\n"
"-v\t\tVerbose\n"
"-N\t\tA comma-separated list of column names: -N\"apple,banana,carrot,durian\"\n"
"-en\t\tIf table exists, do nothing; exit.\n"
"-ed\t\tIf table exists, retain the table, delete all data, refill with the new data (i.e., call 'delete * from your_table').\n"
"-eo\t\tIf table exists, overwrite the table from scratch; deleting the previous table entirely.\n"
"-ea\t\tIf table exists, append new data to the existing table.\n"
"-h\t\tPrint this help\n\n"
, argv[0], argv[0]); 
    int * field_list = NULL;
    char if_exists = 'n';

	if(argc<3){
		printf("%s", msg);
		return 0;
	}
	while ((c = getopt (argc, argv, "n:d:e:f:hmp:ru:vN:O")) != -1)
        if (c=='n') {
              if (optarg[0]=='c') colnames='n';
              else                apop_opts.nan_string = optarg;
        }
		else if (c=='N') {
            apop_data *field_name_data;
            apop_regex(optarg, " *([^,]*[^ ]) *(,|$) *", &field_name_data);
            Apop_stopif(!field_name_data, return 1, 0, "'%s' should be a "
                    "comma-delimited list of field names, but I had trouble "
                    "parsing it as such.", optarg);
            apop_data_transpose(field_name_data);
            field_names = field_name_data->text[0];
        }
        else if (c=='d') strcpy(apop_opts.input_delimiters, optarg);
		else if (c=='f') field_list = break_down(optarg);
		else if (c=='h') {printf("%s", msg); return 0;}
		else if (c=='m') apop_opts.db_engine = 'm';
		else if (c=='u') strcpy(apop_opts.db_user, optarg);
		else if (c=='p') strcpy(apop_opts.db_pass, optarg);
		else if (c=='r') rownames++;
		else if (c=='v') apop_opts.verbose=2;
		else if (c=='O') tab_exists_check++; //deprecated as of December 2013.
		else if (c=='e') {
            if (optarg[0]=='n')       if_exists='n'; //the default anyway.
            else if (optarg[0]=='d')  if_exists='d';
            else if (optarg[0]=='a')  if_exists='a';
            else if (optarg[0]=='o') {if_exists='o';
                                      tab_exists_check++;
                                     }

        }
	apop_db_open(argv[optind + 2]);
    if (tab_exists_check) apop_table_exists(argv[optind+1],1);
    apop_query("begin");
	apop_text_to_db(argv[optind], argv[optind+1], rownames, colnames, field_names, .field_ends=field_list, .if_table_exists=if_exists);
    apop_query("commit");
}
