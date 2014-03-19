/** \file 
 Command line utility to take in a query and put out a Gnuplottable file.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
#include <unistd.h>

char *plot_type = NULL;
int histobins = 0;
int histoplotting = 0;

FILE *open_output(char *outfile, int sf){
    FILE  *f;
    if (sf && !strcmp (outfile, "-"))
        return stdout;
    if (sf && outfile){
        f = fopen(outfile, "w");
        Apop_stopif(!f, exit(0), 0, "Trouble opening %s.", outfile);
        return f;
    }
    f = popen("`which gnuplot` -persist", "w");
    Apop_stopif(!f, exit(0), 0, "Trouble opening %s.", "gnuplot");
    return f;
}

char *read_query(char *infile){
    char in[1000];
    char *q = malloc(10);
    q[0] = '\0';
    FILE *inf = fopen(infile, "r");
    Apop_stopif(!inf, exit(0), 0, "Trouble opening %s. Look into that.\n", infile);
    while(fgets(in, 1000, inf)){
        q = realloc(q, strlen(q) + strlen(in) + 4);
        sprintf(q, "%s%s", q, in);
    }
    sprintf(q, "%s;\n", q);
    fclose(inf);
    return q;
}

gsl_matrix *query(char *d, char *q, int no_plot){
	apop_db_open(d);
    apop_data *result = apop_query_to_data("%s", q);
	apop_db_close(0);
    Apop_stopif(!result, exit(2), 0, "Your query returned a blank table. Quitting.");
    Apop_stopif(result->error, exit(2), 0, "Error running your query. Quitting.");
    if (no_plot){
        apop_data_show(result);
        exit(0);
    }
    return result->matrix;
}

void print_out(FILE *f, char *outfile, gsl_matrix *m){
    if (!histoplotting){
        fprintf(f,"plot '-' with %s\n", plot_type);
	    apop_matrix_print(m, NULL, .output_type='p', .output_pipe=f);
    } else {
        Apop_col_v(&(apop_data){.matrix=m}, 0, v);
        apop_plot_histogram(v, histobins, .output_type='p', .output_pipe=f);
    }
    if (outfile) fclose(f);
}

int main(int argc, char **argv){
    char c, *q = NULL, 
         *d = NULL,
         *outfile = NULL;
    int sf = 0,
        no_plot = 0;

    const char* msg= "%s [opts] dbname query\n\n"
"Runs a query, and pipes the output directly to gnuplot. Use -f to dump to stdout or a file.\n"
"-d\tdatabase to use\t\t\t\t\tmandatory \n"
"-q\tquery to run\t\t\t\t\tmandatory (or use -Q)\n"
"-Q\tfile from which to read the query\t\t\n"
"-n\tno plot: just run the query and display results to stdout\t\t\n"
"-t\tplot type (points, bars, ...)\t\t\tdefault=\"lines\"\n"
"-H\tplot histogram with this many bins (e.g., -H100). To let the system auto-select bin sizes, use -H0 .\n"
"-f\tfile to dump to. If -f- then use stdout.\tdefault=pipe to Gnuplot\n";

	Apop_stopif(argc<2, return 1, 0, msg, argv[0]);
	while ((c = getopt (argc, argv, "ad:f:hH:nQ:q:st:-")) != -1)
	    if (c=='f'){
              outfile = strdup(optarg);
              sf++;
        } else if (c=='H'){
              histoplotting = 1;
              histobins = atoi(optarg);
        } 
        else if (c=='h'||c=='-') {
            printf(msg, argv[0]);
			return 0;
		}
        else if (c=='d') d = strdup(optarg);
        else if (c=='n') no_plot ++;
        else if (c=='Q') q = read_query(optarg);
        else if (c=='q') q = strdup(optarg);
        else if (c=='t') plot_type = strdup(optarg);

    if (optind == argc -2){
        d = argv[optind];
        q = argv[optind+1];
    } else if (optind == argc-1)
        q = argv[optind];
    
    Apop_stopif(!q, return 1, 0, "I need a query specified with -q.\n");

    if (!plot_type) plot_type = strdup("lines");

    FILE *f = open_output(outfile, sf);
    gsl_matrix *m = query(d, q, no_plot);
    print_out(f, outfile, m);
}
