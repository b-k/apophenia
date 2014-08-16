/** \file 
  The \ref apop_update function.  */ 
/* Copyright (c) 2006--2009, 2014 by Ben Klemens. Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
#include <stdbool.h>

/* This file in four parts:
   --an apop_model named product, purpose-built for apop_update to send to apop_model_metropolis
   --apop_mcmc settings and their defaults.
   --apop_update and its equipment, which has three cases:
        --conjugates, in which case see the functions
        --call Metropolis
*/


/* This will be used by apop_update to send to apop_mcmc below.

   To set it up, add a more pointer to an array of two models, the prior and likelihood. 
   The total likelihood of a data point is (likelihood these parameters are drawn from
   prior)*(likelihood of these parameters and the data set using the likelihood fn)
*/
static long double product_ll(apop_data *d, apop_model *m){
    apop_model **pl = m->more;
    gsl_vector *v = apop_data_pack(m->parameters);
    apop_data_unpack(v, pl[1]->parameters);
    gsl_vector_free(v);
    return apop_log_likelihood(m->parameters, pl[0]) + apop_log_likelihood(d, pl[1]);
}

static long double product_constraint(apop_data *data, apop_model *m){
    apop_model **pl = m->more;
    gsl_vector *v = apop_data_pack(m->parameters);
    apop_data_unpack(v, pl[1]->parameters);
    gsl_vector_free(v);
    return pl[1]->constraint(data, pl[1]);
}

apop_model *product = &(apop_model){"product of two models", 
    .log_likelihood=product_ll, .constraint=product_constraint};


///////////the conjugate table

static apop_model *betabinom(apop_data *data, apop_model *prior, apop_model *likelihood){
    apop_model *outp = apop_model_copy(prior);
    if (!data && likelihood->parameters){
        double n = likelihood->parameters->vector->data[0];
        double p = likelihood->parameters->vector->data[1];
        *gsl_vector_ptr(outp->parameters->vector, 0) += n*p;
        *gsl_vector_ptr(outp->parameters->vector, 1) += n*(1-p);
    } else {
        gsl_vector *hits = Apop_cv(data, 1);
        gsl_vector *misses = Apop_cv(data, 0);
        *gsl_vector_ptr(outp->parameters->vector, 0) += apop_sum(hits);
        *gsl_vector_ptr(outp->parameters->vector, 1) += apop_sum(misses);
    }
    return outp;
}

double countup(double in){return in!=0;}

static apop_model *betabernie(apop_data *data, apop_model *prior, apop_model *likelihood){
    apop_model *outp = apop_model_copy(prior);
    Get_vmsizes(data);//tsize
    double sum = apop_map_sum(data, .fn_d=countup, .part='a');
    *gsl_vector_ptr(outp->parameters->vector, 0) += sum;
    *gsl_vector_ptr(outp->parameters->vector, 1) += tsize - sum;
    return outp;
}

static apop_model *gammaexpo(apop_data *data, apop_model *prior, apop_model *likelihood){
    apop_model *outp = apop_model_copy(prior);
    Get_vmsizes(data); //maxsize
    *gsl_vector_ptr(outp->parameters->vector, 0) += maxsize;
    apop_data_set(outp->parameters, 1, .val=1./
                          (1./apop_data_get(outp->parameters, 1) 
                        + (data->matrix ? apop_matrix_sum(data->matrix) : 0)
                        + (data->vector ? apop_sum(data->vector) : 0)));
    return outp;
}

static apop_model *gammapoisson(apop_data *data, apop_model *prior, apop_model *likelihood){
    /* Posterior alpha = alpha_0 + sum x; posterior beta = beta_0/(beta_0*n + 1) */
    apop_model *outp = apop_model_copy(prior);
    Get_vmsizes(data); //vsize, msize1,maxsize
    *gsl_vector_ptr(outp->parameters->vector, 0) +=
                         (vsize  ? apop_sum(data->vector): 0) +
                         (msize1 ? apop_matrix_sum(data->matrix): 0);

    double *beta = gsl_vector_ptr(outp->parameters->vector, 1);
    *beta = *beta/(*beta * maxsize + 1);
    return outp;
}

static apop_model *normnorm(apop_data *data, apop_model *prior, apop_model *likelihood){
/*
output \f$(\mu, \sigma) = (\frac{\mu_0}{\sigma_0^2} + \frac{\sum_{i=1}^n x_i}{\sigma^2})/(\frac{1}{\sigma_0^2} + \frac{n}{\sigma^2}), (\frac{1}{\sigma_0^2} + \frac{n}{\sigma^2})^{-1}\f$

That is, the output is weighted by the number of data points for the
likelihood. If you give me a parametrized normal, with no data, then I'll take the weight to be \f$n=1\f$. 
*/
    double mu_like, var_like;
    long int n;
    apop_model *outp = apop_model_copy(prior);
    apop_prep(data, outp);
    long double  mu_pri = prior->parameters->vector->data[0];
    long double  var_pri = gsl_pow_2(prior->parameters->vector->data[1]);
    if (!data && likelihood->parameters){
        mu_like  = likelihood->parameters->vector->data[0];
        var_like = gsl_pow_2(likelihood->parameters->vector->data[1]);
        n        = 1;
    } else {
        n = data->matrix->size1 * data->matrix->size2;
        apop_matrix_mean_and_var(data->matrix, &mu_like, &var_like);
    }
    gsl_vector_set(outp->parameters->vector, 0, (mu_pri/var_pri + n*mu_like/var_like)/(1/var_pri + n/var_like));
    gsl_vector_set(outp->parameters->vector, 1, pow((1/var_pri + n/var_like), -.5));
    return outp;
}

/** Take in a prior and likelihood distribution, and output a posterior distribution.

This function first checks a table of conjugate distributions for the pair you
sent in. If the names match the table, then the function returns a closed-form
model with updated parameters.  If the parameters aren't in the table of conjugate
priors/likelihoods, then it uses Markov Chain Monte Carlo to sample from the posterior
distribution, and then outputs a histogram model for further analysis. Notably,
the histogram can be used as the input to this function, so you can chain Bayesian
updating procedures.

\li If the prior distribution has a \c p or \c log_likelihood element, then I use \ref apop_model_metropolis to generate the posterior.

\li If the prior does not have a \c p or \c log_likelihood but does have a \c draw element, then I make draws from the prior and weight them by the \c p given by the likelihood distribution. This is not a rejection sampling method, so the burnin is ignored.

Here are the conjugate distributions currently defined:

<table>
<tr>
<td> Prior <td></td> Likelihood  <td></td>  Notes 
</td> </tr> <tr>
<td> \ref apop_beta "Beta" <td></td> \ref apop_binomial "Binomial"  <td></td>  
</td> </tr> <tr>
<td> \ref apop_beta "Beta" <td></td> \ref apop_bernoulli "Bernoulli"  <td></td> 
</td> </tr> <tr>
<td> \ref apop_exponential "Exponential" <td></td> \ref apop_gamma "Gamma"  <td></td>  Gamma likelihood represents the distribution of \f$\lambda^{-1}\f$, not plain \f$\lambda\f$
</td> </tr> <tr>
<td> \ref apop_normal "Normal" <td></td> \ref apop_normal "Normal" <td></td>  Assumes prior with fixed \f$\sigma\f$; updates distribution for \f$\mu\f$
</td></tr> <tr>
<td> \ref apop_gamma "Gamma" <td></td> \ref apop_poisson "Poisson" <td></td> Uses sum and size of the data  
</td></tr>
</table>

\li The conjugate table is stored using a vtable; see \ref vtables for details. The typedef new functions must conform to and the hash used for lookups are:

\code
typedef apop_model *(*apop_update_type)(apop_data *, apop_model , apop_model);
#define apop_update_hash(m1, m2) ((size_t)(m1).draw + (size_t)((m2).log_likelihood ? (m2).log_likelihood : (m2).p)*33)
\endcode

\param data     The input data, that will be used by the likelihood function (default = \c NULL.)
\param  prior   The prior \ref apop_model. If the system needs to
estimate the posterior via MCMC, this needs to have a \c log_likelihood or \c p method.  (No default, must not be \c NULL.)
\param likelihood The likelihood \ref apop_model. If the system needs to
estimate the posterior via MCMC, this needs to have a \c log_likelihood or \c p method (ll preferred). (No default, must not be \c NULL.)
\param rng      A \c gsl_rng, already initialized (e.g., via \ref apop_rng_alloc). (default: an RNG from \ref apop_rng_get_thread)
\return an \ref apop_model struct representing the posterior, with updated parameters. 

Here is a test function that compares the output via conjugate table and via
Metropolis-Hastings sampling: 
\include test_updating.c

\li This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD apop_model * apop_update(apop_data *data, apop_model *prior, apop_model *likelihood, gsl_rng *rng){
    apop_data *apop_varad_var(data, NULL);
    apop_model *apop_varad_var(prior, NULL);
    apop_model *apop_varad_var(likelihood, NULL);
    gsl_rng *apop_varad_var(rng, apop_rng_get_thread());
APOP_VAR_END_HEAD
    static int setup=0; if (!(setup++)){
        apop_update_vtable_add(betabinom, apop_beta, apop_binomial);
        apop_update_vtable_add(betabernie, apop_beta, apop_bernoulli);
        apop_update_vtable_add(gammaexpo, apop_gamma, apop_exponential);
        apop_update_vtable_add(gammapoisson, apop_gamma, apop_poisson);
        apop_update_vtable_add(normnorm, apop_normal, apop_normal);
    }
    apop_update_type conj = apop_update_vtable_get(prior, likelihood);
    if (conj) return conj(data, prior, likelihood);

    apop_mcmc_settings *s = apop_settings_get_group(prior, apop_mcmc);

    bool ll_is_a_copy = false;
    likelihood = maybe_prep(data, likelihood, &ll_is_a_copy); //in apop_mcmc.c
    gsl_vector *pack = apop_data_pack(likelihood->parameters);
    int tsize = pack->size;
    gsl_vector_free(pack);
    Apop_stopif(prior->dsize != tsize, 
                return apop_model_copy(&(apop_model){.error='d'}),
                0, "Size of a draw from the prior does not match "
                   "the size of the likelihood's parameters (%i != %i).%s",
                   prior->dsize, tsize, 
                   (tsize > prior->dsize) ?  
                        " Perhaps use apop_model_fix_params to reduce the "
                        "likelihood's parameter count?" : "");
    if (prior->p || prior->log_likelihood){
        apop_model *p = apop_model_copy(product);
        //pending revision, a memory leak:
        p->more = malloc(sizeof(apop_model*)*2);
        ((apop_model**)p->more)[0] = apop_model_copy(prior);
        ((apop_model**)p->more)[1] = apop_model_copy(likelihood);
        p->more_size = sizeof(apop_model*) * 2;
        p->parameters = apop_data_alloc(prior->dsize);
        p->data = data;
        if (s) apop_settings_copy_group(p, prior, "apop_mcmc");
        apop_model *out = apop_model_metropolis(data, rng, p); 
        if (ll_is_a_copy) apop_model_free(likelihood);
        return out;
    }

    Apop_stopif(!prior->draw, return NULL, 0, "prior does not have a .p, .log_likelihood, or .draw element. I am stumped. Returning NULL.");

    if (!s) s = Apop_model_add_group(prior, apop_mcmc);

    double *draw = malloc(sizeof(double)* (tsize));
    apop_data *out = apop_data_alloc(s->periods, tsize);
    out->weights = gsl_vector_alloc(s->periods);

    apop_draw(draw, rng, prior); //set starting point.
    apop_data_fill_base(likelihood->parameters, draw);

    for (int i=0; i< s->periods; i++){
        newdraw:
        apop_draw(draw, rng, prior);
        apop_data_fill_base(likelihood->parameters, draw);
        long double p = apop_p(data, likelihood);

        Apop_notify(3, "p=%Lg for parameters:\t", p);
        if (apop_opts.verbose >=3) apop_data_print(likelihood->parameters);

        Apop_stopif(gsl_isnan(p), goto newdraw,
                1, "Trouble evaluating the "
                "likelihood function at vector beginning with %g. "
                "Throwing it out and trying again.\n"
                , likelihood->parameters->vector->data[0]);
        apop_data_pack(likelihood->parameters, Apop_rv(out, i));
        gsl_vector_set(out->weights, i, p);
    }
    apop_model *outp = apop_estimate(out, apop_pmf);
    free(draw);
    if (ll_is_a_copy) apop_model_free(likelihood);
    return outp;
}
