/* The Zipf distribution.

Copyright (c) 2005--2009, 2011 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.

\amodel apop_zipf
Wikipedia has notes on the <a href="http://en.wikipedia.org/wiki/Zipf_distribution">Zipf distribution</a>. 
\f$Z(a)   = {1\over \zeta(a) * i^a}        \f$

\f$lnZ(a) = -(\log(\zeta(a)) + a \log(i))    \f$

apop_zipf.estimate() is an MLE, so feed it appropriate \ref apop_mle_settings.

\adoc    Input_format    Ignores the matrix structure of the input data, so send in a 1 x N, an N x 1, or an N x M.

See also \ref apop_data_rank_compress for means of dealing with one more input data format.

\adoc    Parameter_format One item in the parameter set's vector.    
\adoc    settings  \ref apop_mle_settings, \ref apop_parts_wanted_settings    
*/

#include "apop_internal.h"
#include <gsl/gsl_sf_zeta.h>

static long double zipf_constraint(apop_data *returned_beta, apop_model *m){
    //constraint is 1 < beta_1
    Nullcheck_mp(m, GSL_NAN);
    Staticdef(apop_data *, constraint, apop_data_falloc((1,1,1), 1, 1));
    return apop_linear_constraint(m->parameters->vector, constraint, 1e-4);
}

static long double zipf_log_likelihood(apop_data *d, apop_model *m){
    Nullcheck_mpd(d, m, GSL_NAN);
    Get_vmsizes(d) //tsize
    long double bb = apop_data_get(m->parameters, 0, -1);
    Apop_stopif(isnan(bb) || bb < 1, return GSL_NAN, 0, "Zipf needs a parameter >=1; "
                                              "got %Lg. Returning NaN.", bb); 
    double like = -apop_map_sum(d, log) * bb;
    like -= log(gsl_sf_zeta(bb)) * tsize;
    return like;
}    

/*  \adoc RNG Returns a ranking: If the population were Zipf distributed, you're most
likely to get the 1st most common item, so this produces a lot of ones,
a great deal of twos, and so on.

Cribbed from <a href="http://cgm.cs.mcgill.ca/~luc/mbookindex.html>Devroye (1986)</a>, Chapter 10, p 551.  */
static int zipf_rng(double *out, gsl_rng* r, apop_model *param){
    Nullcheck_mp(param, 1);
    double a = apop_data_get(param->parameters, 0, -1);
    Apop_stopif(isnan(a) || a < 1, *out=GSL_NAN; return 1, 
            0, "Zipf needs a parameter >=1; got %g. Setting *out to NAN.", a); 
    int x;
    long double u, v, t, 
            b    = powl(2, a-1), 
            ainv = -(1.0/(a-1));
    do {
        u = gsl_rng_uniform(r);
        v = gsl_rng_uniform(r);
        x = powl(u, ainv);
        t = powl((1.0 + 1.0/x), (a-1));
    } while (v * x * (t-1.0)/(b-1) > t/b);
    *out = x;
    return 0;
}

apop_model *apop_zipf = &(apop_model){"Zipf distribution", 1,0,0, .dsize=1,
     .log_likelihood = zipf_log_likelihood, .constraint = zipf_constraint, .draw = zipf_rng};
