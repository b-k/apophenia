/** \file asst.c  The odds and ends bin.  */

/** Calculate \f$\sum_{n=1}^N {1\over n^s}\f$

There are no doubt efficient shortcuts do doing this, but I use brute
force. To speed things along, I save the results so that they can later
just be looked up.

When reading the code, remember that the zeroth element holds the value
for N=1, and so on.

\todo Look up the tricks for calculating this.
*/

#include <apophenia/headers.h>
double apop_generalized_harmonic(int N, double s){
static double * 	eses	= NULL;
static int * 		lengths	= NULL;
static	int		count	= 0;
static	double **	precalced=NULL;
int			j, old_len,
			i	= 0;
	if (count == 0){
		precalced 	= malloc(sizeof (double*));
		eses 		= malloc(sizeof (double));
		lengths		= malloc(sizeof (int));
	}
	while (i< count)
		if (eses[i] == s) 	break;
		else 			i++;
	if (i == count){	//you need to build the vector from scratch.
		count			++;
		precalced 		= realloc(precalced, sizeof (double*) * count);
		lengths 		= realloc(lengths, sizeof (int*) * count);
		eses 			= realloc(eses, sizeof (double) * count);
		precalced[count-1]	= malloc(sizeof(double) * N);
		lengths[count-1]	= N;
		eses[count-1]		= s;
		precalced[count-1][0]	= 1;
		old_len			= 1;
		i			= count -1;
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

/** test the generalized harmonic summing thing.

\bug If this is called from outside the library (e.g., the test program), it sends back wrong answers. Potentially a gcc bug!  */
int test_harmonic(){
double		out;
int		count = 0;
	out	= apop_generalized_harmonic(270, 0.0);
	if(out !=270){
		printf("Generalized harmonic(270,0) should be 270, but it's %g. Fail.\n", out);
		count++;
	}
	out	= apop_generalized_harmonic(370, -1.0);
	if(out !=370*371/2){
		printf("Generalized harmonic(370,-1) should be 370*371/2, but it's %g. Fail.\n", out);
		count++;
	}
	out	= apop_generalized_harmonic(12, -1.0);
	if(out !=12*13/2){
		printf("Generalized harmonic(12,-1) should be 12*13/2, but it's %g. Fail.\n", out);
		count++;
	}
	return count;
}

