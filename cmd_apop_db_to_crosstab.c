/** \file cmd_apop_db_to_crosstab.c	Command line utility to convert a three-column table to a crosstab.*/

/*Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
#include <unistd.h>

int main(int argc, char **argv){
char		c,  verbose=0,
		    *delimiter,
            *outfile    = NULL,
		    msg[1000];
apop_data	*m;

	sprintf(msg, "%s [opts] dbname table_name rows columns data\n\n"
            "-d\tdelimiter\t\tdefault= \"|,<space><tab>\"\n"
            "-a\tappend\t\t\tdefault= append\n"
            "-o\toverwrite\t\tdefault= append\n"
            "-v\tverbose: prints status info on stderr and raises apop_opts.verbose by one for each use (so use -v -v for extra-verbose)\n"
            "-f\tfile to dump to\t\tdefault=STDOUT\n", argv[0]); 

	if(argc<5){
		printf("%s", msg);
		return 0;
	}
	delimiter	= malloc(5);
	strcpy(delimiter, ",");
	while ((c = getopt (argc, argv, "ad:f:ho")) != -1){
		switch (c){
		  case 'a':
              apop_opts.output_append = 1;
			  break;
		  case 'd':
			  strcpy(apop_opts.output_delimiter,optarg);
			  break;
		  case 'o':
              apop_opts.output_append = 0;
			  break;
		  case 'f':
              outfile   = malloc(1000);
			  sprintf(outfile, "%s", optarg);
			  apop_opts.output_type	= 'f';
			  break;
		  case 'h':
			printf("%s", msg);
			return 0;
		  case 'v':
            verbose++;
            apop_opts.verbose++;
            break;
		}
	}
    Apop_assert(optind+4 <= argc, "I need five arguments past the options: database, table, row col, column col, data col");
    if (verbose){
        fprintf(stderr, "database:%s\ntable: %s\nrow col: %s\ncol col:%s\ndata col:%s\n",
            argv[optind], argv[optind +1], argv[optind+2], argv[optind+3], argv[optind+4]);
        if (outfile) fprintf(stderr, "outfile: %s\n", outfile);
        else  fprintf(stderr, "output to stdout\n");
        if (apop_opts.output_append) fprintf(stderr, "appending to output\n");
        else fprintf(stderr, "overwriting output\n");
    }
	apop_db_open(argv[optind]);
	m	= apop_db_to_crosstab(argv[optind +1], argv[optind+2], argv[optind+3], argv[optind+4]);
	apop_data_print(m, outfile);
}
