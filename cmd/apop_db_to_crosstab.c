/** \file 
Command line utility to convert a three-column table to a crosstab.*/

/*Copyright (c) 2005--2007, 2013 by Ben Klemens.  Licensed under the GPLv2; see COPYING.  */

#include "apop_internal.h"
#include <unistd.h>

int main(int argc, char **argv){
    int c;
    char verbose=0;
    char const *msg="Usage: %s [opts] dbname table_name rows columns data\n"
"\n"
"A command-line wrapper for the apop_db_to_crosstab function.\n"
"See Apophenia's online documentation for that function for details and tricks.\n"
"The default for the data column is a count [count(*)]\n"
"The column is optional; leave it out if you want a single-dimensional crosstab.\n"
"If you need a non-default data column but want a 1-D crosstab, use 1 as your column.\n"
"\n"
" -d\tdelimiter (default: <tab>)\n"
" -v\tverbose: prints status info on stderr\n"
" -v -v\tvery verbose: also print queries executed on stderr\n"
" -h\tdisplay this help and exit\n"
"\n";

	apop_opts.verbose=0;  //so don't print queries until -v -v.

	while ((c = getopt (argc, argv, "d:f:hv-")) != -1)
		if      (c=='d') strcpy(apop_opts.output_delimiter,optarg);
        else if (c=='h'||c=='-') {printf(msg, argv[0]); exit(0);}
        else if (c=='v') {
            verbose++;
            apop_opts.verbose++;
        }

    Apop_stopif(optind+2 > argc, return 1, 0, "I need at least two arguments past the options: database table [optional rowcol] [optional columncol] [optional datacol]");
    _Bool no_rowcol = optind+2 > argc;
    _Bool no_columncol = optind+3 > argc;
    _Bool no_datacol = optind+4 > argc;
    char *rowcol = no_rowcol    ? "1" : argv[optind+2];
    char *colcol = no_columncol ? "1" : argv[optind+3];
    char *datacol = no_datacol  ? NULL: argv[optind+4];
    if (verbose){
        fprintf(stderr, "database:%s\ntable: %s\nrow col: %s\ncol col:%s%s%s\n---------\n",
            argv[optind], argv[optind +1], rowcol, colcol,
            no_datacol ?"":"\ndata col:", datacol);
    }
	apop_db_open(argv[optind]);
	apop_data *m = apop_db_to_crosstab(argv[optind +1], rowcol, colcol, datacol);
	apop_data_print(m);
}
