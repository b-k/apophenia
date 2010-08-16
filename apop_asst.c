/** \file apop_asst.c  The odds and ends bin. 
Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "asst.h"
#include "types.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

/** Calculate \f$\sum_{n=1}^N {1\over n^s}\f$

\li There are no doubt efficient shortcuts do doing this, but I use brute force. [Though Knuth's Art of Programming v1 doesn't offer anything, which is strong indication of nonexistence.] To speed things along, I save the results so that they can just be looked up should you request the same calculation. 

For example: \include test_harmonic.c
*/
double apop_generalized_harmonic(int N, double s){
/* 
Each row in the saved-results structure is an \f$s\f$, and each column is \f$1\dots n\f$, up to the largest \f$n\f$ calculated to date.

When reading the code, remember that the zeroth element holds the value for N=1, and so on.
*/
  static double * 	eses	= NULL;
  static int * 		lengths	= NULL;
  static int		  count	= 0;
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

/** Are two strings equal?
 
<tt>strcmp()</tt> is a distance function, but its most common use, by
far, is just to ask whether two strings are equal. This function directly addresses
that question. It differs from <tt>strcmp()</tt> in three ways, one
especially important.

Tests for equality: I answer the question {\em Is string one == string
two?}, so if two strings are identical, I return 1; if they differ, I
return 0. {\em This is the opposite of strcmp}, which returns zero when
two strings are identical (and are thus zero distance apart). It is common for
people to think strcmp answers the equality question, and then wind up writing the wrong
thing. If that's you, consider using this function (or one like it) and
never touching strcmp. If you think the fact that it's confusing that \c apop_strcmp and \c strcmp return opposite
outputs, feel free to never use this function.

Faster: as soon as I determine a difference, I leave. That is, I'll step
through both strings until I find a character that differs, and will
leave at that point. <tt>strcmp</tt> will always step through the whole of both
strings.

Accepts <tt>NULL</tt>s: 
If one string is <tt>NULL</tt> and the other isn't I return 0 (because a
<tt>NULL</tt> string surely differs from a non-<tt>NULL</tt>. If both strings
are <tt>NULL</tt>, that's equality, and I return 1. If you're looking for
different behavior on <tt>NULL</tt>s, you're best off just testing your
inputs.

\param one The first string to compare
\param two The other string to compare

For example: \include test_strcmp.c
*/
int apop_strcmp(char *one, char *two){
    if (!one && !two)
        return 1;
    if ((!one && two) || (one && !two))
        return 0;
    for ( ; *one && *two; one++, two++)
        if(*one != *two)
            return 0;
    if (*one || *two) //different length strings.
        return 0;
    return 1;
}

/** Strip dots from a name.

\param  in          A string
\param  strip_type  'd': replace all '.' with '_'.<br>
                    'b': return only the string before the '.', so 'table.col' becomes 'table'. If there are multiple dots, cuts off at the first dot.
                    'a': return only the string after the '.', so 'table.col' becomes 'col'. If there are multiple dots, cuts off at the last dot.

For example: \include test_strip_dots.c
\ingroup convenience_fns
 */
char * apop_strip_dots(char *in, char strip_type){
int     i;
char    *out    = NULL;
    if ((strip_type ==0) || (strip_type == 'd') || (strip_type == 'D')){
        out    = malloc(strlen(in)+1);
        for (i=0; i< strlen(in)+1; i++)   //will copy over the '/0' too.
            out[i] = (in[i] == '.') ? '_' : in[i];
    }
    else if ((strip_type ==1) || (strip_type == 'b') || (strip_type == 'B')){
        out    = malloc(strlen(in)+1);
        strcpy(out, in);
        for (i=strlen(in)+1; i--; )
            if (in[i] == '.'){
                out[i] = '\0';
                break;
            }
    }
    else if ((strip_type ==2) || (strip_type == 'a') || (strip_type == 'A')){
        for (i=0; i< strlen(in)+1; i++)
            if (in[i] == '.')
                break;
        out    = malloc(strlen(in)-i);
        strcpy(out, (in+i+1));
    }
    return out;
}

/** Inform the user of a faux pas. See also \ref Apop_assert, which allows the function to return a value.

  Deprecated. Use \ref Apop_assert, \ref Apop_assert_c, or \ref apop_notify.

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
        fprintf(stderr, "%s", message);
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


/* \defgroup sorting Sorting functions

 A few functions to sort data. One sorts an \c apop_data set in place, and one returns percentiles for a sorted vector.
  \{ */

#include "variadic.h"
#include "likelihoods.h"
#include <gsl/gsl_sort_vector.h>

static int find_min_unsorted(size_t *sorted, size_t height, size_t min){
    while (min<height)
        if (!sorted[min])   return min;
        else                min++;
    return -1;
}

/** This function sorts the whole of a \c apop_data set based on one column. Sorts in place, with little additional memory used.

 Uses the \c gsl_sort_vector_index function internally, and that function just ignores NaNs; therefore this function just leaves NaNs exactly where they lay.

 \param data    The input set to be modified. (No default, must not be \c NULL.)
 \param sortby  The column of data by which the sorting will take place. As usual, -1 indicates the vector element. (default: column zero of the matrix)
 \param asc   If 'd' or 'D', sort in descending order; else sort in ascending order. (Default: ascending)
 \return A pointer to the data set, so you can do things like \c apop_data_show(apop_data_sort(d, -1)).

This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD apop_data * apop_data_sort(apop_data *data, int sortby, char asc){
    apop_data * apop_varad_var(data, NULL);
    apop_assert_s(data, "You gave me NULL data to sort.");
    int apop_varad_var(sortby, 0);
    char apop_varad_var(asc, 0);
APOP_VAR_ENDHEAD
  size_t            height  = (sortby==-1) ? data->vector->size: data->matrix->size1;
  size_t            sorted[height];
  size_t            i, *perm, start=0;
  gsl_permutation   *p  = gsl_permutation_alloc(height);
    memset(sorted, 0, sizeof(size_t)*height);
    if (sortby == -1)
        gsl_sort_vector_index (p, data->vector);
    else {
        APOP_COL(data, sortby, v);
        gsl_sort_vector_index (p, v);
    }
    perm    = p->data;
    if (asc=='d' || asc=='D') //reverse the perm matrix.
        for (size_t j=0; j< height/2; j++){
            double t         = perm[j];
            perm[j]          = perm[height-1-j];
            perm[height-1-j] = t;
        }
    while (1){
        i           =
        start       = find_min_unsorted(sorted, height, start);
        if (i==-1) 
            break;
        Apop_data_row(data, start, firstrow);
        apop_data *first_row_storage = apop_data_copy(firstrow);
        sorted[start]++;
        while (perm[i]!=start){
            //copy from perm[i] to i
            Apop_data_row(data, perm[i], onerow);
            apop_data_set_row(data, onerow, i);
            sorted[perm[i]]++;
            i   = perm[i];
        }
        apop_data_set_row(data, first_row_storage, i);
        apop_data_free(first_row_storage);
    }
    gsl_permutation_free(p);
    return data;
}


/** Returns an array of size 101, where returned_vector[95] gives the value of the 95th percentile, for example. Returned_vector[100] is always the maximum value, and returned_vector[0] is always the min (regardless of rounding rule).

\param data	a gsl_vector of data. (No default, must not be \c NULL.)
\param rounding This will either be 'u', 'd', or 'a'. Unless your data is exactly a multiple of 101, some percentiles will be ambiguous. If 'u', then round up (use the next highest value); if 'd' (or anything else), round down to the next lowest value; if 'a', take the mean of the two nearest points. If 'u' or 'a', then you can say "5% or more  of the sample is below returned_vector[5]"; if 'd' or 'a', then you can say "5% or more of the sample is above returned_vector[5]".   (Default = 'd'.)

\li You may eventually want to free() the array returned by this function.
This function uses the \ref designated syntax for inputs.
*/ 
APOP_VAR_HEAD double * apop_vector_percentiles(gsl_vector *data, char rounding){
    gsl_vector *apop_varad_var(data, NULL);
    apop_assert_s(data, "You gave me NULL data.");
    char apop_varad_var(rounding, 'd');
APOP_VAR_ENDHEAD
  gsl_vector	*sorted	= gsl_vector_alloc(data->size);
  double		*pctiles= malloc(sizeof(double) * 101);
	gsl_vector_memcpy(sorted,data);
	gsl_sort_vector(sorted);
	for(int i=0; i<101; i++){
		int index = i*(data->size-1)/100.0;
		if (rounding == 'u' && index != i*(data->size-1)/100.0)
			index ++; //index was rounded down, but should be rounded up.
		if (rounding == 'a' && index != i*(data->size-1)/100.0)
            pctiles[i]	= (gsl_vector_get(sorted, index)+gsl_vector_get(sorted, index+1))/2.;
        else pctiles[i]	= gsl_vector_get(sorted, index);
	}
	gsl_vector_free(sorted);
	return pctiles;
}

/** \} */

static int count_parens(const char *string){
    int out = 0;
    int last_was_backslash = 0;
    for(const char *step =string; *step !='\0'; step++){
        if (*step == '\\' && !last_was_backslash){
            last_was_backslash = 1;
            continue;
        }
        if (*step == ')' && !last_was_backslash)
            out++;
        last_was_backslash = 0;
    }
    return out;
}

/** A convenience function for regular expression searching 

\li There are three common flavors of regular expression: Basic, Extended,
and Perl-compatible (BRE, ERE, PCRE). I use EREs, as per the specs of
your C library, which should match POSIX's ERE specification. 

For example, "p.val" will match "P value", "p.value", "p values" (and even "tempeval", so be
careful).

If you give a non-\c NULL address in which to place a table of paren-delimited substrings, I'll return them as a row in the text element of the returned \ref apop_data set. I'll return <em>all</em> the matches, filling the first row with substrings from the first application of your regex, then filling the next row with another set of matches (if any), and so on to the end of the string. Useful when parsing a list of items, for example.


\param string        The string to search (no default; if \c NULL, I return 0---no match)
\param regex       The regular expression (no default)
\param substrings   Parens in the regex indicate that I should return matching substrings. Give me the _address_ of an \ref apop_data* set, and I will allocate and fill the text portion with matches. Default= \c NULL, meaning do not return substrings (even if parens exist in the regex).
\param use_case         Should I be case sensitive, \c 'y' or \c 'n'? (default = \c 'n', which is not the POSIX default.)

\return         1 == match; 0 == no match. \c substrings may be allocated and filled if needed.
\ingroup names

\li Here is the test function. Notice that the substring-pulling
function call passes \c &subs, not plain \c subs. Also, the non-match
has a zero-length blank in <tt>subs->text[0][1]</tt>.
\include test_regex.c
*/
APOP_VAR_HEAD int  apop_regex(const char *string, const char* regex, apop_data **substrings, const char use_case){
    const char * apop_varad_var(string, NULL);
    if (!string) return 0;
    const char * apop_varad_var(regex, NULL);
    apop_assert_s(regex, "You gave me a NULL regex.");
    apop_data **apop_varad_var(substrings, NULL);
    const char apop_varad_var(use_case, 'y');
APOP_VAR_ENDHEAD
  regex_t    re;
  int        matchcount=count_parens(regex);
  int        found;
  regmatch_t result[matchcount];
    int compiled_ok = !regcomp(&re, regex, REG_EXTENDED 
                                            + (use_case=='y' ? REG_ICASE : 0)
                                            + (substrings ? 0 : REG_NOSUB) );
    apop_assert_s(compiled_ok, "This regular expression didn't compile: \"%s\"", regex)

    int matchrow = 0;
    do {
        found = !regexec(&re, string, matchcount+1, result, matchrow ? REG_NOTBOL : 0);
        if (substrings && found){
            *substrings = apop_text_alloc(matchrow ? *substrings : NULL, matchrow+1, matchcount);
            //match zero is the whole string; ignore.
            for (int i=0; i< matchcount; i++){
                if (result[i+1].rm_eo > 0){//GNU peculiarity: match-to-empty marked with -1.
                    int length_of_match = result[i+1].rm_eo - result[i+1].rm_so;
                    (*substrings)->text[matchrow][i] = malloc(strlen(string)+1);
                    memcpy((*substrings)->text[matchrow][i], string + result[i+1].rm_so, length_of_match);
                    (*substrings)->text[matchrow][i][length_of_match] = '\0';
                } else{ //matches nothing
                    (*substrings)->text[matchrow][i] = malloc(1);
                    (*substrings)->text[matchrow][i][0]='\0';
                }
            }
            string += result[0].rm_eo; //end of whole match;
            matchrow++;
        }
    } while (substrings && found && string[0]!='\0');
    regfree(&re);
    return found;
}
