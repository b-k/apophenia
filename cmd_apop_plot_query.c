/** \file cmd_apop_plot_query.c	Command line utility to take in a query and put out a Gnuplottable file.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "db.h"
#include "model.h"
#include "output.h"
#include "conversions.h"
#include "linear_algebra.h"
#include <unistd.h>

char *plot_type = NULL;
int histobins   = 0;

FILE *open_output(char *outfile, int sf){
  FILE  *f;
    if (sf && !strcmp (outfile, "-"))
        return stdout;
    if (sf && outfile){
        f   = fopen(outfile, "w");
        if (!f){
            fprintf(stderr, "Trouble opening %s. Look into that.\n", outfile);
            exit(0);
        }
        return f;
    } 
    f   = popen("`which gnuplot` -persist", "w");
    if (!f){
        fprintf(stderr, "Trouble opening %s. Look into that.\n", "gnuplot");
        exit(0);
    }
    return f;
}

char *read_query(char *infile){
  char in[1000];
  char *q       = malloc(10);
    q[0]        = '\0';
  FILE  *inf    = fopen(infile, "r");
    if (!inf){
        fprintf(stderr, "Trouble opening %s. Look into that.\n", infile);
        exit(0);
    }
    while(fgets(in, 1000, inf)){
        q   = realloc(q, strlen(q) + strlen(in) + 4);
        sprintf(q, "%s%s", q, in);
    }
    sprintf(q, "%s;\n", q);
    fclose(inf);
    return q;
}

gsl_matrix *query(char *d, char *q, int no_plot){
	apop_db_open(d);
apop_data *result 	= apop_query_to_data("%s", q);
	apop_db_close(0);
    Apop_assert(result, "Your query returned a blank table. Quitting.");
    if (no_plot){
        apop_data_show(result);
        exit(0);
    }
    return result->matrix;
}

void print_out(FILE *f, char *outfile, gsl_matrix *m){
    apop_opts.output_type = 'p';
    apop_opts.output_pipe = f;
    if (!histobins){
        fprintf(f,"plot '-' with %s\n", plot_type);
	    apop_matrix_print(m, NULL);
    }
    else {
        APOP_MATRIX_COL(m, 0, v);
        apop_plot_histogram(v, histobins, NULL);
    }
    if (outfile) fclose(f);
}

int main(int argc, char **argv){
  char		    c, *q       = NULL, 
                *d          = NULL,
                *outfile    = NULL,
		        msg[2000];
  gsl_matrix	*m;
  FILE          *f;
  int           sf          = 0,
                no_plot     = 0;

	sprintf(msg, "%s [opts] dbname query\n\n\
Runs a query, and pipes the output directly to gnuplot. Use -f to dump to STDOUT or a file.\n\
-d\tdatabase to use\t\t\t\t\tmandatory \n\
-q\tquery to run\t\t\t\t\tmandatory (or use -Q)\n\
-Q\tfile from which to read the query\t\t\n\
-n\tno plot: just run the query and display results to STDOUT\t\t\n\
-t\tplot type (points, bars, ...)\t\t\tdefault=\"lines\"\n\
-H\tplot histogram with this many bins (e.g., -H100)\n\
-f\tfile to dump to. If -f- then use STDOUT.\tdefault=pipe to Gnuplot\n", argv[0]); 

	if(argc<2){
		printf("%s", msg);
		return 0;
	}
	while ((c = getopt (argc, argv, "ad:f:hH:nQ:q:st:")) != -1){
		switch (c){
		  case 'd':
              d   = malloc(2+strlen(optarg));
			  sprintf(d, "%s", optarg);
			  break;
		  case 'f':
              outfile   = malloc(1000);
			  sprintf(outfile, "%s", optarg);
			  apop_opts.output_type	= 'f';
              sf  ++;
			  break;
          case 'H':
              histobins = atoi(optarg);
              break;
		  case 'n':
              no_plot ++;
			  break;
		  case 'Q':
              q   = read_query(optarg);
			  break;
		  case 'q':
              q   = malloc(2+strlen(optarg));
			  sprintf(q, "%s", optarg);
			  break;
		  case 't':
              plot_type   = malloc(2+strlen(optarg));
			  sprintf(plot_type, "%s", optarg);
			  break;
		  case 'h':
			printf("%s", msg);
			return 0;
		}
	}
    if (!plot_type){
              plot_type   = malloc(20);
			  sprintf(plot_type, "lines");
    }
    if (optind == argc -2){
        d = argv[optind];
        q = argv[optind+1];
    } else if (optind == argc-1)
        q = argv[optind];
    
    if (!q){
        fprintf(stderr, "I need a query specified with -q.\n");
        return 0;
    }

    f   = open_output(outfile, sf);
    m   = query(d, q, no_plot);
    print_out(f, outfile, m);
}
