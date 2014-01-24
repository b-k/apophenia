/* The Bernoulli distribution as an \ref apop_model.
Copyright (c) 2007--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

/* \amodel apop_bernoulli The Bernoulli model: A single random draw with probability \f$p\f$.

\adoc   Input_format
  The matrix or vector can have any size, and I just count up zeros
  and non-zeros. The Bernoulli parameter \f$p\f$ is the percentage of non-zero
  values in the matrix. Its variance is \f$p(1-p)\f$.

\adoc    Parameter_format A vector of length one */

#include "apop_internal.h"

static double bernie_ll(double x, void * pin){ 
    double *p = pin; 
    return x ? log(*p) : log(1-*p); 
}

static long double bernoulli_log_likelihood(apop_data *d, apop_model *params){
    Nullcheck_mpd(d, params, GSL_NAN);
    double p = apop_data_get(params->parameters, 0, -1);
	return apop_map_sum(d, .fn_dp = bernie_ll, .param=&p);
}

static double nonzero (double in) { return in !=0; }

/* \adoc estimated_parameters \f$p\f$ is the only element in the vector. A
<tt>\<Covariance\></tt> page has the variance of \f$p\f$ in the (0,0)th element of the matrix.
\adoc estimated_info   Reports <tt>log likelihood</tt>.
*/
static void bernoulli_estimate(apop_data * data,  apop_model *est){
    Nullcheck_mpd(data, est, ); Get_vmsizes(data); //tsize;
    double n = tsize;
    double p = apop_map_sum(data, nonzero)/n;
    apop_name_add(est->parameters->names, "p", 'r');
	gsl_vector_set(est->parameters->vector, 0, p);
    apop_data_add_named_elmt(est->info, "log likelihood", bernoulli_log_likelihood(data, est));
    apop_data *cov = apop_data_add_page(est->parameters, apop_data_alloc(1,1), "<Covariance>");
    apop_data_set(cov, 0,0, p*(1-p));
}

static long double bernoulli_constraint(apop_data *data, apop_model *inmodel){
    //constraint is 0 < b and  1 > b
    Staticdef(apop_data *, constraint, apop_data_falloc((2,2,1), 0., 1.,
                                                                -1., -1.));
    return apop_linear_constraint(inmodel->parameters->vector, constraint, 1e-3);
}

/* \adoc    RNG Returns a single zero or one. */
static int bernoulli_rng(double *out, gsl_rng *r, apop_model* eps){
    *out = gsl_rng_uniform (r) < eps->parameters->vector->data[0]; 
    return 0;
}

static long double bernoulli_cdf(apop_data *d, apop_model *params){
//One of those functions that just fills out the form.
//CDF to zero = 1-p
//CDF to one = 1
    Nullcheck_mpd(d, params, GSL_NAN); Get_vmsizes(d)  //firstcol
    double val = apop_data_get(d, .col=firstcol);
    double p = *params->parameters->vector->data;
    return isnan(val) ? GSL_NAN
            : val < 0 ? 0
            : val < 1 ? 1 - p
                      : 1;
}

static void bernie_print(apop_model *m, FILE *out){
    fprintf(out, "Bernoulli distribution with p = %g.\n", apop_data_get(m->parameters,0,-1));
}

static void bernie_prep(apop_data *data, apop_model *params){
    apop_model_print_vtable_add(bernie_print, apop_bernoulli);
    apop_model_clear(data, params);
}

/* \adoc Settings None. */
apop_model *apop_bernoulli = &(apop_model){"Bernoulli distribution", 1, .dsize=1, .prep=bernie_prep,
	.estimate = bernoulli_estimate, .log_likelihood = bernoulli_log_likelihood, 
   .constraint = bernoulli_constraint, .cdf = bernoulli_cdf, .draw = bernoulli_rng};
