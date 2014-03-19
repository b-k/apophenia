/** \file 
Command line utility to convert a three-column table to a crosstab.*/

/*Copyright (c) 2005--2007, 2013 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
#include <unistd.h>

int main(int argc, char **argv){
    char c, verbose=0;
    char const *msg="%s [opts] dbname table_name rows columns data\n\n"
        "A command-line wrapper for the apop_db_to_crosstab function.\n"
        "See Apophenia's online documentation for that function for details and tricks.\n"
          "-d\tdelimiter\t\tdefault=<tab>\n"
          "-v\tverbose: prints status info on stderr\n"
          "-v -v\textra verbose: also print queries executed on stderr\n";

	Apop_stopif(argc<5, return 1, 0, msg, argv[0]);

    apop_opts.verbose=0;  //so don't print queries until -v -v.

	while ((c = getopt (argc, argv, "d:f:hv-")) != -1)
		if      (c=='d') strcpy(apop_opts.output_delimiter,optarg);
        else if (c=='h'||c=='-') {printf(msg, argv[0]); exit(0);}
        else if (c=='v') {
            verbose++;
            apop_opts.verbose++;
        }

    Apop_stopif(optind+4 > argc, return 1, 0, "I need five arguments past the options: database, table, row col, column col, data col");
    if (verbose){
        fprintf(stderr, "database:%s\ntable: %s\nrow col: %s\ncol col:%s\ndata col:%s\n---------\n",
            argv[optind], argv[optind +1], argv[optind+2], argv[optind+3], argv[optind+4]);
    }
	apop_db_open(argv[optind]);
    apop_data *m = apop_db_to_crosstab(argv[optind +1], argv[optind+2], argv[optind+3], argv[optind+4]);
	apop_data_print(m);
}
