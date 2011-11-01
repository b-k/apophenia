/** \file cmd_apop_lookup.c	Command line utility to look things up in distribution tables. */

/*  As with many pleasant utility programs, it is mostly parsing of the input text, and then about three lines at the end doing actual work.

Copyright (c) 2008 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
#include "apop_internal.h"
#include <unistd.h>

char *plot_type = NULL;
int histobins   = 0;

typedef enum {Beta, Binomial, F, Negbinom, Normal, Poisson, T} distlist ;

int main(int argc, char **argv){
  distlist      distribution = Normal;
  char	        msg[10000], c;
  int           pval = 0, qval = 0;
  double        param1 = GSL_NAN, param2 =GSL_NAN, findme = GSL_NAN;
  char          number[1000];
	sprintf(msg, "%s [opts] number_to_lookup\n\n"
    "Look up a probability or p-value for a given standard distribution.\n"
    "[This is still loosely written and counts as beta. Notably, negative numbers are hard to parse.]\n"
    "E.g.:\n"
    "%s -dbin 100 .5 34\n"
    "sets the distribution to a Binomial(100, .5), and find the odds of 34 appearing.\n"
    "%s -p 2     \n"
    "find the area of the Normal(0,1) between -infty and 2.  \n"
    "\n"
    "-pval Find the p-value: integral from -infinity to your value\n"
    "-qval Find the q-value: integral from your value to infinity\n"
    "\n"
    "After giving an optional -p or -q, specify the distribution. \n"
    "Default is Normal(0, 1). Other options:\n"
    "\t\t-binom Binomial(n, p)\n"
    "\t\t-beta Beta(a, b)\n"
    "\t\t-f F distribution(df1, df2)\n"
    "\t\t-norm Normal(mu, sigma)\n"
    "\t\t-negative bin Negative binomial(n, p)\n"
    "\t\t-poisson Poisson(L)\n"
    "\t\t-t t distribution(df)\n"
    "I just need enough letters to distinctly identify a distribution.\n"
, argv[0], argv[0], argv[0]); 

    opterr=0;
	if(argc==1){
		printf("%s", msg);
		return 0;
	}
	while ((c = getopt (argc, argv, "B:b:F:f:N:n:pqT:t:")) != -1){
		switch (c){
		  case 'B':
		  case 'b':
              if (optarg[0]=='i')
                  distribution = Binomial;
              else if (optarg[0]=='e')
                  distribution = Beta;
            else {
                printf("I can't parse the option -b%s\n", optarg);
                exit(0);
            }
              param1 = atof(argv[optind]);
              param2 = atof(argv[optind+1]);
              findme =  atof(argv[optind+2]);
			  break;
          case 'F':
          case 'f':
            distribution = F;
            param1 = atof(argv[optind]);
            findme =  atof(argv[optind+1]);
            break;
          case 'H':
		  case 'h':
			printf("%s", msg);
			return 0;
          case 'n':
          case 'N':
            if (optarg[0]=='o'){ //normal
                  param1 = atof(argv[optind]);
                  param2 = atof(argv[optind+1]);
                  findme =  atof(argv[optind+2]);
            } else if (optarg[0]=='e'){
                  distribution = Negbinom;
                  param1 = atof(argv[optind]);
                  param2 = atof(argv[optind+1]);
                  findme =  atof(argv[optind+2]);
            } else {
                printf("I can't parse the option -n%s\n", optarg);
                exit(0);
            }
			  break;
          case 'p':
            if (!optarg || optarg[0] == 'v')
                pval++;
            else if (optarg[0] == 'o'){
                distribution = Poisson;
                param1 = atof(argv[optind]);
                findme =  atof(argv[optind+1]);
            } else {
                printf("I can't parse the option -p%s\n", optarg);
                exit(0);
            }
            break;
          case 'q':
            qval++;
            break;
          case 'T':
          case 't':
            distribution = T;
            param1 = atof(argv[optind]);
            findme =  atof(argv[optind+1]);
            break;
          case '?'://probably a negative number
            if (optarg)
                 snprintf(number, 1000, "%c%s", optopt, optarg);
            else snprintf(number, 1000, "%c", optopt);
            if (gsl_isnan(param1)) param1 = -atof(number);
            else if (gsl_isnan(param2)) param2 = -atof(number);
            else if (gsl_isnan(findme)) findme = -atof(number);
		}
	}
    if (gsl_isnan(findme)) findme =  atof(argv[optind]);
    //defaults, as promised
    if (gsl_isnan(param1)) param1 = 0;
    if (gsl_isnan(param2)) param2 = 1;
    if (!pval && !qval){
        double val =
        distribution == Beta ? gsl_ran_beta_pdf(findme, param1, param2)
        : distribution == Binomial ? gsl_ran_binomial_pdf(findme, param2, param1)
        : distribution == F ? gsl_ran_fdist_pdf(findme, param1, param2) 
        : distribution == Negbinom ? gsl_ran_negative_binomial_pdf(findme, param2, param1)
        : distribution == Normal ? gsl_ran_gaussian_pdf(findme, param2)+param1
        : distribution == Poisson ? gsl_ran_poisson_pdf(findme, param1) 
        : distribution == T ? gsl_ran_tdist_pdf(findme, param1) : GSL_NAN;
        printf("%g\n", val); 
        return 0;
    }
    if (distribution == Binomial){
        printf("Sorry, the GSL doesn't have a Binomial CDF.\n");
        return 0; }
    if (distribution == Negbinom){
        printf("Sorry, the GSL doesn't have a Negative Binomial CDF.\n");
        return 0; }
    if (distribution == Poisson){
        printf("Sorry, the GSL doesn't have a Poisson CDF.\n");
        return 0; }
    if (pval){
        double val =
        distribution == Beta ? gsl_cdf_beta_P(findme, param1, param2)
        : distribution == F ? gsl_cdf_fdist_P(findme, param1, param2) 
        : distribution == Normal ? gsl_cdf_gaussian_P(findme-param1, param2)
        : distribution == T ? gsl_cdf_tdist_P(findme, param1) : GSL_NAN;
        printf("%g\n", val); 
        return 0;
    }
    if (qval){
        double val =
        distribution == Beta ? gsl_cdf_beta_Q(findme, param1, param2)
        : distribution == F ? gsl_cdf_fdist_Q(findme, param1, param2) 
        : distribution == Normal ? gsl_cdf_gaussian_Q(findme-param1, param2)
        : distribution == T ? gsl_cdf_tdist_Q(findme, param1) : GSL_NAN;
        printf("%g\n", val); 
    }
}
