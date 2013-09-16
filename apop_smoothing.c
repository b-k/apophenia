/** \file apop_smoothing.c	A few smoothing-type functions, like moving averages. */

/* Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"

/** Return a new vector that is the moving average of the input vector.
 \param v The input vector, unsmoothed
 \param bandwidth The number of elements to be smoothed. 
 */
gsl_vector *apop_vector_moving_average(gsl_vector *v, size_t bandwidth){
    Apop_stopif(!v, return NULL, 0, "You asked me to smooth a NULL vector; returning NULL.\n");
    Apop_stopif(!bandwidth, return apop_vector_copy(v), 0, "Bandwidth must be >=1. Returning a copy of original vector with no smoothing.");
    int halfspan = bandwidth/2;
    gsl_vector *vout = gsl_vector_calloc(v->size - halfspan*2);
    for(size_t i=0; i < vout->size; i ++){
        double *item = gsl_vector_ptr(vout, i);
        for (int j=-halfspan; j < halfspan+1; j ++)
            *item += gsl_vector_get(v, j+ i+ halfspan);
        *item /= halfspan*2 +1;
    }
    return vout;
}
