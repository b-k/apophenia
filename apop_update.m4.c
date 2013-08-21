/** \file 
  The \ref apop_update function.  The header is in asst.h. */ 
/* Copyright (c) 2006--2009 by Ben Klemens. Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
#include "vtables.h"

Apop_settings_init(apop_update,
   Apop_varad_set(periods, 6e3);
   Apop_varad_set(burnin, 0.05);
   Apop_varad_set(method, 'd'); //default
   //all else defaults to zero/NULL
)

Apop_settings_copy(apop_update, )
Apop_settings_free(apop_update, )

static apop_model *betabinom(apop_data *data, apop_model prior, apop_model likelihood){
    apop_model *outp = apop_model_copy(prior);
    if (!data && likelihood.parameters){
        double n = likelihood.parameters->vector->data[0];
        double p = likelihood.parameters->vector->data[1];
        *gsl_vector_ptr(outp->parameters->vector, 0) += n*p;
        *gsl_vector_ptr(outp->parameters->vector, 1) += n*(1-p);
    } else {
        double y = apop_matrix_sum(data->matrix);
        *gsl_vector_ptr(outp->parameters->vector, 0) += y;
        *gsl_vector_ptr(outp->parameters->vector, 1) += data->matrix->size1*data->matrix->size2 - y;
    }
    return outp;
}

double countup(double in){return in!=0;}

static apop_model *betabernie(apop_data *data, apop_model prior, apop_model likelihood){
    apop_model *outp = apop_model_copy(prior);
    Get_vmsizes(data);//tsize
    double sum = apop_map_sum(data, .fn_d=countup, .part='a');
    *gsl_vector_ptr(outp->parameters->vector, 0) += sum;
    *gsl_vector_ptr(outp->parameters->vector, 1) += tsize - sum;
    return outp;
}

static apop_model *gammaexpo(apop_data *data, apop_model prior, apop_model likelihood){
    apop_model *outp = apop_model_copy(prior);
    *gsl_vector_ptr(outp->parameters->vector, 0) += data->matrix->size1*data->matrix->size2;
    apop_data_set(outp->parameters, 1, -1, 1/(1/apop_data_get(outp->parameters, 1, -1) + apop_matrix_sum(data->matrix)));
    return outp;
}

static apop_model *gammapoisson(apop_data *data, apop_model prior, apop_model likelihood){
    /* Posterior alpha = alpha_0 + sum x; posterior beta = beta_0/(beta_0*n + 1) */
    apop_model *outp = apop_model_copy(prior);
    Get_vmsizes(data); //vsize, msize1
    double sum = 0;
    if (vsize)  sum = apop_sum(data->vector);
    if (msize1) sum += apop_matrix_sum(data->matrix);
    *gsl_vector_ptr(outp->parameters->vector, 0) += sum;

    double *beta = gsl_vector_ptr(outp->parameters->vector, 1);
    *beta = *beta/(*beta * tsize + 1);
    return outp;
}

static apop_model *normnorm(apop_data *data, apop_model prior, apop_model likelihood){
/*
output \f$(\mu, \sigma) = (\frac{\mu_0}{\sigma_0^2} + \frac{\sum_{i=1}^n x_i}{\sigma^2})/(\frac{1}{\sigma_0^2} + \frac{n}{\sigma^2}), (\frac{1}{\sigma_0^2} + \frac{n}{\sigma^2})^{-1}\f$

That is, the output is weighted by the number of data points for the
likelihood. If you give me a parametrized normal, with no data, then I'll take the weight to be \f$n=1\f$. 
*/
    double mu_like, var_like;
    long int n;
    apop_model *outp = apop_model_copy(prior);
    long double  mu_pri    = prior.parameters->vector->data[0];
    long double  var_pri = gsl_pow_2(prior.parameters->vector->data[1]);
    if (!data && likelihood.parameters){
        mu_like  = likelihood.parameters->vector->data[0];
        var_like = gsl_pow_2(likelihood.parameters->vector->data[1]);
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

To change the default settings (periods, burnin...),
add an \ref apop_update_settings struct to the prior.

\li If the likelihood model no parameters, I will allocate them. That means you can use
one of the stock models that ship with Apophenia. If I need to run the model's prep
routine to get the size of the parameters, then I'll make a copy of the likelihood
model, run prep, and then allocate parameters for that copy of a model.

\li Consider the state of the \c parameters element of your likelihood model to be
undefined when this exits. This may be settled at a later date.

\li If you set <tt>apop_opts.verbose=2</tt>, I will report the accept rate of the Gibbs sampler. It is a common rule of thumb to select a prior so that this is between 20% and 50%. Set <tt>apop_opts.verbose=3</tt> to see the proposal points, their likelihoods, and the acceptance odds.

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

\param data     The input data, that will be used by the likelihood function (default = \c NULL.)
\param  prior   The prior \ref apop_model. If the system needs to
estimate the posterior via MCMC, this needs to have a \c draw method.  (No default, must not be \c NULL.)
\param likelihood The likelihood \ref apop_model. If the system needs to
estimate the posterior via MCMC, this needs to have a \c log_likelihood or \c p method (ll preferred). (No default, must not be \c NULL.)
\param rng      A \c gsl_rng, already initialized (e.g., via \ref apop_rng_alloc). (default: see \ref autorng)
\return an \ref apop_model struct representing the posterior, with updated parameters. 
\todo The table of conjugate prior/posteriors (in its static \c check_conjugacy subfuction), is a little short, and can always be longer.

Here is a test function that compares the output via conjugate table and via
Gibbs sampling: 
\include test_updating.c

\li This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD apop_model * apop_update(apop_data *data, apop_model *prior, apop_model *likelihood, gsl_rng *rng){
    static gsl_rng *spare_rng = NULL;
    apop_data *apop_varad_var(data, NULL);
    apop_model *apop_varad_var(prior, NULL);
    apop_model *apop_varad_var(likelihood, NULL);
    gsl_rng *apop_varad_var(rng, NULL);
    if (!rng){
        if (!spare_rng) spare_rng = apop_rng_alloc(++apop_opts.rng_seed);
        rng = spare_rng;
    }
APOP_VAR_END_HEAD
    static int setup=0; if (!(setup++)){
        apop_update_insert(betabinom, apop_beta, apop_binomial);
        apop_update_insert(betabernie, apop_beta, apop_bernoulli);
        apop_update_insert(gammaexpo, apop_gamma, apop_exponential);
        apop_update_insert(gammapoisson, apop_gamma, apop_poisson);
        apop_update_insert(normnorm, apop_normal, apop_normal);
    }
    apop_update_type conj = apop_update_get(*prior, *likelihood);
    if (conj) return conj(data, *prior, *likelihood);

    apop_update_settings *s = apop_settings_get_group(prior, apop_update);
    if (!s) s = Apop_model_add_group(prior, apop_update);
    int ll_is_a_copy=0;
    if (!likelihood->parameters){
        if ( likelihood->vbase  >= 0 &&     // A hackish indication that
             likelihood->m1base >= 0 &&     // there is still prep to do.
             likelihood->m2base >= 0 && likelihood->prep){
                ll_is_a_copy++;
                likelihood =apop_model_copy(*likelihood);
                apop_prep(data, likelihood);
        }
        likelihood->parameters = apop_data_alloc(likelihood->vbase, likelihood->m1base, likelihood->m2base);
    }
    Get_vmsizes(likelihood->parameters) //vsize, msize1, msize2
    double    ratio, ll, cp_ll = GSL_NEGINF;
    double    *draw          = malloc(sizeof(double)* (vsize+msize1*msize2));
    apop_data *current_param = apop_data_alloc(vsize , msize1, msize2);
    Apop_stopif(s->burnin > 1, s->burnin/=(s->periods+0.0), 
                1, "Burn-in should be a fraction of the number of periods, "
                    "not a whole number of periods. Rescaling to burnin=%g", s->burnin/=(s->periods+0.0));
    apop_data *out = apop_data_alloc(s->periods*(1-s->burnin), vsize+msize1*msize2);
    int accept_count = 0;

    apop_draw(draw, rng, prior); //set starting point.
    apop_data_fill_base(current_param, draw);

    for (int i=0; i< s->periods; i++){     //main loop
        newdraw:
        apop_draw(draw, rng, prior);
        apop_data_fill_base(likelihood->parameters, draw);
        ll = apop_log_likelihood(data,likelihood);

        Apop_notify(3, "ll=%g for parameters:\t", ll);
        if (apop_opts.verbose >=3) apop_data_print(likelihood->parameters);

        Apop_stopif(gsl_isnan(ll), goto newdraw, 
                1, "Trouble evaluating the "
                "likelihood function at vector beginning with %g. "
                "Throwing it out and trying again.\n"
                , likelihood->parameters->vector->data[0]);
        ratio = ll - cp_ll;
        if (ratio >= 0 || log(gsl_rng_uniform(rng)) < ratio){
            apop_data_memcpy(current_param, likelihood->parameters);
            cp_ll = ll;
            accept_count++;
        } else {
            Apop_notify(3, "reject, with exp(ll_now-ll_prior) = exp(%g-%g) = %g.", ll, cp_ll, exp(ratio));
        }
        if (i >= s->periods * s->burnin){
            Apop_matrix_row(out->matrix, i-(s->periods *s->burnin), v)
            apop_data_pack(current_param, v);
        }
    }
    out->weights = gsl_vector_alloc(s->periods*(1-s->burnin));
    gsl_vector_set_all(out->weights, 1);
    apop_model *outp   = apop_estimate(out, apop_pmf);
    free(draw);
    if (ll_is_a_copy) apop_model_free(likelihood);
    Apop_notify(2, "Gibbs sampling accept percent = %3.3f%%\n", 100*(0.0+accept_count)/s->periods);
    return outp;
}
