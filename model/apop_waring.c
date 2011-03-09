/* The Waring distribution.  
Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  

\amodel apop_waring  The Waring distribution

\f$W(x,k, b,a) 	= (b-1) \gamma(b+a) \gamma(k+a) / [\gamma(a+1) \gamma(k+a+b)]\f$

\f$\ln W(x, b, a) = \ln(b-1) + \ln\gamma(b+a) + \ln\gamma(k+a) - \ln\gamma(a+1) - \ln\gamma(k+a+b)\f$

\f$dlnW/db	= 1/(b-1)  + \psi(b+a) - \psi(k+a+b)\f$

\f$dlnW/da	= \psi(b+a) + \psi(k+a) - \psi(a+1) - \psi(k+a+b)\f$

\adoc    Input_format     
Ignores the matrix structure of the input data, so send in a 1 x N, an N x 1, or an N x M.

See also \ref apop_data_rank_compress for means of dealing with one more input data format.
     
\adoc    Parameter_format  Two elements in the parameter set's vector.    
\adoc    settings   MLE-type: \ref apop_mle_settings, \ref apop_parts_wanted_settings    
*/

#include "apop_internal.h"

typedef struct {
    double bb, a;
} ab_type;

static double beta_zero_and_one_greater_than_x_constraint(apop_data *returned_beta, apop_model *m){
    //constraint is 1 < beta_1 and  0 < beta_2
    Staticdef(apop_data *, constraint, apop_data_fill(apop_data_calloc(2,2,2),
                             1., 1., 0.,
                             0., 0., 1.));
    return apop_linear_constraint(m->parameters->vector, constraint, 1e-4);
}

static double apply_me(double val, void *in){
    ab_type *ab = in;
        double ln_a_k		 = gsl_sf_lngamma(val + ab->a);
        double ln_bb_a_k	 = gsl_sf_lngamma(val + ab->a + ab->bb);
        return  ln_a_k - ln_bb_a_k;
}

static double waring_log_likelihood(apop_data *d, apop_model *m){
  Nullcheck_mpd(d, m);
  Get_vmsizes(d) //tsize
  ab_type abstruct = {.bb = gsl_vector_get(m->parameters->vector, 0),
                      .a  = gsl_vector_get(m->parameters->vector, 1)};
  double ln_bb_a = gsl_sf_lngamma(abstruct.bb + abstruct.a),
         ln_a_mas_1	= gsl_sf_lngamma(abstruct.a + 1),
         ln_bb_less_1= log(abstruct.bb - 1);
    double likelihood  = apop_map_sum(d, .fn_dp = apply_me, .param=&abstruct);
	likelihood += (ln_bb_less_1 + ln_bb_a - ln_a_mas_1)* tsize;
    printf("..%g\t %g\n", abstruct.bb, abstruct.a);
	return likelihood;
}

static double dapply_a(double val, void *in){
    ab_type *ab = in;
        long double psi_a_k		 = gsl_sf_psi(val + ab->a);
        long double psi_bb_a_k	 = gsl_sf_psi(val + ab->a + ab->bb);
        return psi_a_k - psi_bb_a_k;
}

static double dapply_b(double val, void *in){
    ab_type *ab = in;
    long double psi_bb_a_k = gsl_sf_psi(val + ab->a + ab->bb);
    return -psi_bb_a_k;
}

/*static*/ void waring_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
	//Psi is the derivative of the log gamma function.
    Nullcheck_mpd(d, m);
    Get_vmsizes(d) //tsize
    ab_type ab = {.bb = gsl_vector_get(m->parameters->vector, 0),
                  .a  = gsl_vector_get(m->parameters->vector, 1)};
    apop_data_show(m->parameters);
    double d_bb = apop_map_sum(d, .fn_dp=dapply_b, .param=&ab);
    double d_a  = apop_map_sum(d, .fn_dp=dapply_a, .param=&ab);
    double	bb_minus_one_inv= 1./(ab.bb-1),
      	    psi_a_bb	    = gsl_sf_psi(ab.bb + ab.a),
      	    psi_a_mas_one	= gsl_sf_psi(ab.a+1);
	d_bb += (bb_minus_one_inv + psi_a_bb) *tsize;
	d_a	 += (psi_a_bb- psi_a_mas_one) * tsize;
	gsl_vector_set(gradient, 0, d_bb);
	gsl_vector_set(gradient, 1, d_a);
    apop_vector_show(gradient);
}

/* \adoc RNG Give me parameters, and I'll draw a ranking from the appropriate
Waring distribution. [I.e., if I randomly draw from a Waring-distributed
population, return the ranking of the item I just drew.]

Source: page seven of
L. Devroye, <a href="http://cgm.cs.mcgill.ca/~luc/digammapaper.ps">Random
variate generation for the digamma and trigamma distributions</a>, Journal
of Statistical Computation and Simulation, vol. 43, pp. 197-216, 1992.
*/
static void waring_rng(double *out, gsl_rng *r, apop_model *eps){
//The key to convert from Devroye's GHgB3 notation to what I
//consider to be the standard Waring notation in \ref apop_waring:
// a = a + 1
// b = 1 
// c = b - 1 
// n = k - 1 , so if it returns 0, that's first rank.
// OK, I hope that clears everything up.
    double		x, u,
                b   = gsl_vector_get(eps->parameters->vector, 0),
                a   = gsl_vector_get(eps->parameters->vector, 1);
    double		params[]	={a+1, 1, b-1};
	do{
		x	= apop_rng_GHgB3(r, params)+1;
		u	= gsl_rng_uniform(r);
	} while (u >= (x + a)/(GSL_MAX(a+1,1)*x));
    *out = x;
}

apop_model apop_waring = {"Waring distribution", 2,0,0, .dsize=1,
	 .log_likelihood =  waring_log_likelihood,
     .constraint =  beta_zero_and_one_greater_than_x_constraint,  .draw = waring_rng};
 /*.score = waring_dlog_likelihood,  //seems numerically unstable*/
