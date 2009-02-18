/** \file apop_asst.c  The odds and ends bin. 
Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "types.h"
#include "vasprintf/vasprintf.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

/** Calculate \f$\sum_{n=1}^N {1\over n^s}\f$

There are no doubt efficient shortcuts do doing this, but I use brute
force. To speed things along, I save the results so that they can later
just be looked up. Each row in the saved structure is an \f$s\f$, and each
column is \f$1\dots n\f$, up to the largest \f$n\f$ calculated to date.

When reading the code, remember that the zeroth element holds the value
for N=1, and so on.

\todo Look up the tricks for calculating this.
*/
double apop_generalized_harmonic(int N, double s){
static double * 	eses	= NULL;
static int * 		lengths	= NULL;
static int		    count	= 0;
static double **	precalced=NULL;
int			        j, old_len, i;
	for (i=0; i< count; i++)
		if (eses == NULL || eses[i] == s) 	
            break;
	if (i == count){	//you need to build the vector from scratch.
		count			++;
        i               = count - 1;
		precalced 		= realloc(precalced, sizeof (double*) * count);
		lengths 		= realloc(lengths, sizeof (int*) * count);
		eses 			= realloc(eses, sizeof (double) * count);
		precalced[i]	= malloc(sizeof(double) * N);
		lengths[i]	    = N;
		eses[i]		    = s;
		precalced[i][0]	= 1;
		old_len			= 1;
	}
	else {	//then you found it.
		old_len		= lengths[i];
	}
	if (N-1 >= old_len){	//It's there, but you need to extend what you have.
		precalced[i]	= realloc(precalced[i],sizeof(double) * N);
		for (j=old_len; j<N; j++)
			precalced[i][j] = precalced[i][j-1] + 1/pow((j+1),s);
	}
	return 	precalced[i][N-1];
}

/** RNG from a Generalized Hypergeometric type B3.

 Devroye uses this as the base for many of his
 distribution-generators, e.g., \ref apop_waring "apop_waring.rng". 
*/  //Header in stats.h
double apop_rng_GHgB3(gsl_rng * r, double* a){
if ((a[0]<=0) || (a[1] <= 0) || (a[2] <=0)){
	printf("apop_GHgB3_rng took a zero parameter; bad.\n");
	return 0;
	}
double		aa	= gsl_ran_gamma(r, a[0], 1),
		b	= gsl_ran_gamma(r, a[1], 1),
		c	= gsl_ran_gamma(r, a[2], 1);
int		p;
	p	= gsl_ran_poisson(r, aa*b/c);
	return p;
}


/** Strip dots from a name.

\param  in          A string
\param  strip_type  'd': replace all '.' with '_'.<br>
                    'b': return only the string before the '.', so 'table.col' becomes 'col'. If there are multiple dots, cuts off at the first dot.
                    'a': return only the string after the '.', so 'table.col' becomes 'col'. If there are multiple dots, cuts off at the last dot.
\ingroup convenience_fns
 */
char * apop_strip_dots(char *in, char strip_type){
int     i;
char    *out    = NULL;
    if ((strip_type ==0) || (strip_type == 'd')){
        out    = malloc(strlen(in)+1);
        for (i=0; i< strlen(in)+1; i++){//will copy over the '/0' too.
            if (in[i] == '.')
                out[i] = '_';
            else
                out[i] = in[i];
        }
    }
    else if ((strip_type ==1) || (strip_type == 'b')){
        out    = malloc(strlen(in)+1);
        strcpy(out, in);
        for (i=strlen(in)+1; i--; )
            if (in[i] == '.'){
                out[i] = '\0';
                break;
            }
    }
    else if ((strip_type ==2) || (strip_type == 'a')){
        for (i=0; i< strlen(in)+1; i++)
            if (in[i] == '.')
                break;
        out    = malloc(strlen(in)-i);
        strcpy(out, (in+i+1));
    }
    return out;
}

/** Inform the user of a faux pas.

  Notice that the message is the last parameter, since it is probably long.

 \param level   At what verbosity level should the user be warned? E.g., if level==2, then print iff apop_opts.verbosity >= 2. You can set apop_opts.verbose==-1 to turn off virtually all messages, but this is probably ill-advised.
 \param stop   Either 's' or 'c', indicating whether the program should stop or continue. If stopping, uses \c assert(0) for easy debugging. You can use 'h' (halt) as a synonym for 's'.
 \param msg The message to write to STDERR (presuming the verbosity level is high enough). This can be a printf-style format with following arguments. You can produce much more informative error messages this way, e.g., \c apop_error(0, 's', "Beta is %g but should be greater than zero.", beta);.
*/
void apop_error(int level, char stop, char *msg, ...){
  va_list   argp;
  char      *message;
    va_start(argp, msg);
    vasprintf(&message, msg, argp);
    va_end(argp);

    if (apop_opts.verbose >= level)
        fprintf(stderr, message);
    free(message);
    if (stop == 's' || stop == 'h')
        assert(0);
}


/** Call \c system(), but with <tt>printf</tt>-style arguments. E.g.,
  
 \code
char filenames[] = "apop_asst.c apop_asst.o"
apop_system("ls -l %s", filenames);
\endcode

\return The return value of the \c system() call.

 */
int apop_system(const char *fmt, ...){
  char 		*q;
  va_list   argp;
	va_start(argp, fmt);
	vasprintf(&q, fmt, argp);
	va_end(argp);
    int out = system(q);
    free(q);
    return out;
}
