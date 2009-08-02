/** \file apop_yule.c

  The Yule distribution. A special case of the Waring.*/ 
/*Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "types.h"
#include "mapply.h"
#include "likelihoods.h"

static double yule_log_likelihood_rank(apop_data *d, apop_model *m){
  float         bb	            = gsl_vector_get(m->parameters->vector, 0);
  int 		    k;
  float 	    ln_k, ln_bb_k,
		        likelihood 	    = 0,
		        ln_bb		    = gsl_sf_lngamma(bb),
		        ln_bb_less_1    = log(bb-1);
	for (k=0; k< d->matrix->size2; k++)	{
		if (k>=1) 	ln_k	= gsl_sf_lngamma(k+1);
		else		ln_k	= 0;
		ln_bb_k		  = gsl_sf_lngamma(k+1+bb);
        APOP_COL(d, k, v);
		likelihood   += apop_sum(v) *  (ln_k - ln_bb_k);
	}
	likelihood   += (ln_bb_less_1 + ln_bb) * d->matrix->size1 * d->matrix->size2;
	return likelihood;
}

static void yule_dlog_likelihood_rank(apop_data *d, gsl_vector *gradient, apop_model *p){
  float         bb		        = gsl_vector_get(p->parameters->vector, 0);
  gsl_matrix    *data	        = d->matrix;
  int 		    k;
  double	    bb_minus_one_inv= 1/(bb-1),
		        psi_bb	        = gsl_sf_psi(bb),
		        psi_bb_k,
		        d_bb		= 0;
	for (k=0; k< data->size2; k++){
		psi_bb_k= gsl_sf_psi(k +1 + bb);
        APOP_COL(d, k, v);
		d_bb		-= apop_sum(v) * psi_bb_k;
	}
	d_bb   += (bb_minus_one_inv + psi_bb) * data->size1 * data->size2;
	gsl_vector_set(gradient, 0, d_bb);
}



static double beta_greater_than_x_constraint(apop_data *returned_beta, apop_model *m){
  apop_assert(m->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
    //constraint is 1 < beta_1
  static apop_data *constraint = NULL;
    if (!constraint){
        constraint= apop_data_calloc(1,1,1);
        apop_data_set(constraint, 0, 0, 1);
        apop_data_set(constraint, 0, -1, 1);
        }
    return apop_linear_constraint(m->parameters->vector, constraint, 1e-3);
}

static long double  bb;

static double  apply_me(gsl_vector *v){
  long double   likelihood =0, 
                ln_k, ln_bb_k;
  int 	        k;
  double        pt;
    for (k=0; k< v->size; k++)	{
        pt            = gsl_vector_get(v, k);
		if (pt>=1) 	ln_k	= gsl_sf_lngamma(pt);
		else		ln_k	= 0;
        ln_bb_k		  = gsl_sf_lngamma(pt+bb);
        likelihood   +=  ln_k - ln_bb_k;
    }
    return likelihood;
}

static double  dapply_me(gsl_vector *v){
  int           k;
  double        pt; 
  long double   d   = 0;
    for (k=0; k< v->size; k++)	{
        pt   = gsl_vector_get(v, k);
        d	-= gsl_sf_psi(pt + bb);
    }
    return d;
}

static double yule_log_likelihood(apop_data *d, apop_model *m){
  apop_assert(m->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
  if (apop_settings_get_group(m, "apop_rank"))
      return yule_log_likelihood_rank(d, m);
    bb	            = gsl_vector_get(m->parameters->vector, 0);
    long double   ln_bb		    = gsl_sf_lngamma(bb),
                  ln_bb_less_1    = log(bb-1);
    gsl_vector *  v               = apop_matrix_map(d->matrix, apply_me);
    double        likelihood      = apop_vector_sum(v);
    gsl_vector_free(v);
	return likelihood + (ln_bb_less_1 + ln_bb) * d->matrix->size1 * d->matrix->size2;
}

static void yule_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
  apop_assert_void(m->parameters, 0,'s', "You asked me to evaluate an un-parametrized model.");
    if (apop_settings_get_group(m, "apop_rank"))
      return yule_dlog_likelihood_rank(d, gradient, m);
	//Psi is the derivative of the log gamma function.
    bb		= gsl_vector_get(m->parameters->vector, 0);
double		    bb_minus_one_inv= 1/(bb-1),
		        psi_bb	        = gsl_sf_psi(bb);
  gsl_vector *  v               = apop_matrix_map(d->matrix, dapply_me);
  double        d_bb            = apop_vector_sum(v);
    gsl_vector_free(v);
    d_bb    += (bb_minus_one_inv + psi_bb) * d->matrix->size1 * d->matrix->size2;
	gsl_vector_set(gradient, 0, d_bb);
}


/** Draw from a Yule distribution with parameter a

Call this fn using <tt> apop_draw(*out, r, apop_yule)</tt>.

\param	a	The parameter.
\param	r	A gsl_rng that you've already set up.

Cribbed from <a href="http://cgm.cs.mcgill.ca/~luc/mbookindex.html>Devroye (1986)</a>, p 553.  */
static void yule_rng( double *out, gsl_rng * r, apop_model *a){
double 	e1, e2;
int		x;
	e1	= gsl_ran_exponential(r, 1);
	e2	= gsl_ran_exponential(r, 1);
	x	= GSL_MAX((int) (- e1  / log(1 - exp(-e2 / (*a->parameters->vector->data -1)))), 0);
	*out =  x + 1;	//we rounded down to floor, but want ceil.
}


/** The Yule distribution

The special case of Waring where \f$ \alpha = 0.	\f$<br>

apop_yule.estimate() is an MLE, so feed it appropriate \ref apop_mle_settings.

\f$ Y(x, b) 	= (b-1) \gamma(b) \gamma(k) / \gamma(k+b)			\f$

\f$ \ln Y(x, b)	= \ln(b-1) + ln\gamma(b) + \ln\gamma(k) - \ln\gamma(k+b)	\f$

\f$ d\ln Y/db	= 1/(b-1)  + \psi(b) - \psi(k+b)				\f$

To specify that you have frequency or ranking data, use 
\code
Apop_settings_add_group(your_model, apop_rank, NULL);
\endcode

\ingroup models
\todo I'm pretty sure Wikipedia's specification of the Yule is wrong; I should check and fix when I have references on hand.
*/
apop_model apop_yule = {"Yule", 1,0,0, .log_likelihood = yule_log_likelihood, 
    .score = yule_dlog_likelihood, .constraint = beta_greater_than_x_constraint, 
    .draw = yule_rng};
//estimate via the default MLE method
