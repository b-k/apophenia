#include <unistd.h>
#include <stdio.h>
#include <stdlib.h> //atof
#include <math.h>   //powf

int main(int argc, char ** argv){
    char opts[]= "M:m:i:h";
    char help[]= "\
A program to take powers of a function. Usage\n\
\t\tpowers [options] [a number]\n\
-h\t This help\n\
-m\t The minimum exponent at which to start.\n\
-M\t The maximum exponent at which to finish.\n\
-i\t Increment by this.\n";
    char c;
    double  i;
    double  min    = 0.;
    double  max    = 10.;
    double  incr    = 1.;
    double base = 2.;

    if (argc==1) {
        printf(help);
        return 0;
    }
    while ( (c=getopt(argc, argv, opts)) != -1){
        if (c=='h'){
            printf(help);
            return 0;
        } else if (c=='m'){
            min = atof(optarg);
        } else if (c=='M'){
            max = atof(optarg);
        } else if (c=='i'){
            incr = atof(optarg);
        }
    }
    if (optind == argc-1)
        base    = atof(argv[optind]);
    for (i=min; i<=max; i+= incr)
        printf("%g^%g: %g\n", base, i, powf(base, i));
    return 0;
}
