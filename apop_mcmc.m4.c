/** \file 
  Markov Chain Monte Carlo. */ 
/* Copyright (c) 2014 by Ben Klemens. Licensed under the GNU GPL v2; see COPYING. */

#include "apop_internal.h"
#include <stdbool.h>


///default step and adapt fns.

static void step_to_vector(double const *d, apop_mcmc_proposal_s *ps, apop_mcmc_settings *ms){
    apop_model *m = ps->proposal;
    memcpy(m->parameters->vector->data, d, sizeof(double)*m->parameters->vector->size);

    ps->adapt_fn(ps, ms);
}

int sigma_adapt(apop_mcmc_proposal_s *ps, apop_mcmc_settings *ms){
    apop_model *m = ps->proposal;
    //accept rate. Add 1% * target to numerator; 1% to denominator, to slow early jumps
    double ar = (ps->accept_count + .01*ms->periods *ms->target_accept_rate)
               /(ps->accept_count + ps->reject_count + .01*ms->periods);
/*    double std_dev_scale= (ar > ms->target_accept_rate) 
                        ? (2 - (1.-ar)/(1.-ms->target_accept_rate))
                        : (1/(2-((ar+0.0)/ms->target_accept_rate)));
                        */
    double scale = ar/ms->target_accept_rate;
    scale = 1+ (scale-1)/100.;
    //gsl_matrix_scale(m->parameters->matrix, scale > .1? ( scale < 10 ? scale : 10) : .1);
    gsl_matrix_scale(m->parameters->matrix, scale);
    return 0;
}

/////// apop_mcmc_settings

Apop_settings_init(apop_mcmc,
   Apop_varad_set(periods, 6e3);
   Apop_varad_set(burnin, 0.05);
   Apop_varad_set(target_accept_rate, 0.35);
   Apop_varad_set(gibbs_chunks, 'b');
   Apop_varad_set(start_at, '1');
   Apop_varad_set(base_step_fn, step_to_vector);
   Apop_varad_set(base_adapt_fn, sigma_adapt);
   //all else defaults to zero/NULL
)

Apop_settings_copy(apop_mcmc, 
    if (in->block_count){
        out->proposals = calloc(in->block_count, sizeof(apop_mcmc_proposal_s));
        for (int i=0; i< in->block_count; i++){
            out->proposals[i] = in->proposals[i];
            out->proposals[i].proposal = apop_model_copy(in->proposals[i].proposal);
        }
        out->proposal_is_cp=1;
    }
)

Apop_settings_free(apop_mcmc, 
        if (in->proposal_is_cp) {
            for (int i=0; i< in->block_count; i++)
                apop_model_free(in->proposals[i].proposal);
        free(in->proposals);
        }
)

static void setup_normal_proposals(apop_mcmc_proposal_s *s, int tsize, apop_mcmc_settings *settings){
    apop_model *mvn =  apop_model_copy(apop_multivariate_normal);
    mvn->parameters = apop_data_alloc(tsize, tsize, tsize);
    gsl_vector_set_all(mvn->parameters->vector, 1);
    gsl_matrix_set_identity(mvn->parameters->matrix);
    s->proposal = mvn;
    s->step_fn = settings->base_step_fn;
    s->adapt_fn = settings->base_adapt_fn;
}

static void set_block_count_and_block_starts(apop_data *in, 
                                  apop_mcmc_settings *s, size_t total_len){
    if (s->gibbs_chunks =='a') {
        s->block_count = 1;
        s->block_starts = calloc(2, sizeof(size_t));
        s->block_starts[1] = total_len;
    } else if (s->gibbs_chunks =='b') {
        s->block_count = 0;
        for (apop_data *d = in; d; d=d->more)
            s->block_count += !!d->vector + !!d->matrix + !!d->weights;

        s->block_starts = calloc(s->block_count+1, sizeof(size_t));
        int this=1, ctr=0;
        for (apop_data *d = in; d; d=d->more){
            #define markit(test, value) if (test)  \
                s->block_starts[this++] = ctr += value; 

            markit(d->vector, d->vector->size);
            markit(d->matrix, d->matrix->size1*d->matrix->size2);
            markit(d->weights, d->weights->size);
        }
    } else { // item-by-item
        s->block_count = total_len;
        s->block_starts = calloc(total_len+1, sizeof(size_t));
        for (int i=1; i<total_len+1; i++) s->block_starts[i] = i;
    }
}

static void one_step(apop_data *d, gsl_vector *draw, apop_model *m, apop_mcmc_settings *s, gsl_rng *rng, int *constraint_fails, apop_data *out, size_t block, int out_row){
    gsl_vector *clean_copy = apop_vector_copy(draw);
    newdraw:
    apop_draw(draw->data + s->block_starts[block], rng, s->proposals[block].proposal);
    apop_data_unpack(draw, m->parameters);
    if (m->constraint && m->constraint(d, m)){
        (*constraint_fails)++;
        goto newdraw;
    }
    double ll = apop_log_likelihood(d, m);

    Apop_notify(3, "ll=%g for parameters:\t", ll);
    if (apop_opts.verbose >=3) apop_data_print(m->parameters, .output_pipe=apop_opts.log_file);

    Apop_stopif(gsl_isnan(ll) || !isfinite(ll), goto newdraw, 
            1, "Trouble evaluating the m function at vector beginning with %g. "
            "Throwing it out and trying again.\n"
            , m->parameters->vector->data[0]);

    double ratio = ll - s->last_ll;
    if (ratio >= 0 || log(gsl_rng_uniform(rng)) < ratio){//success
        if (s->proposals[block].step_fn) 
            s->proposals[block].step_fn(draw->data + s->block_starts[block], s->proposals+block, s);
        s->last_ll = ll;
        s->proposals[block].accept_count++;
        s->accept_count++;
    } else {
        s->proposals[block].reject_count++;
        s->reject_count++;
        Apop_notify(3, "reject, with exp(ll_now-ll_proposal) = exp(%g-%g) = %g.", ll, s->last_ll, exp(ratio));
        gsl_vector_memcpy(draw, clean_copy);
        apop_data_unpack(draw, m->parameters); //keep the last success in m->parameters.
    }
    if (out_row>=0) gsl_vector_memcpy(Apop_rv(out, out_row), draw);
}


/** The draw method for models estimated via \ref apop_model_metropolis.

That method produces an \ref apop_pmf, typically with a few thousand draws from the
model in a batch. If you want to get a single next step from the Markov chain, use this.

A Markov chain works by making a new draw and then accepting or rejecting the draw. If
the draw is rejected, the last value is reported as the next step in the chain. Users
sometimes mitigate this repetition by making a batch of draws (say, ten at a time) and 
using only the last.

If you run this without first running \ref apop_model_metropolis, I will run it for
you, meaning that there will be an initial burn-in period before the first draw that
can be reported to you. That run is done using \c model->data as input.

\param out An array of \c doubles, which will hold the draw, in the style of \ref apop_draw.
\param rng A \c gsl_rng, already initialized, probably via \ref apop_rng_alloc.
\param model A model which was probably already run through \ref apop_model_metropolis.
\return On return, \c out is filled with the next step in the Markov chain. The <tt>->data</tt> element of the PMF model is extended to include the additional steps in the chain.
If a proposal failed the model constraints, then return 1; else return 0. See the notes in the documentation for \ref apop_model_metropolis.

\li After pulling the attached settings group, the parent model is ignored. One expects that \c base_model in the mcmc settings group == the parent model.

\li If your settings break the model parameters into several chunks, this function returns after stepping through all chunks.
*/
int apop_model_metropolis_draw(double *out, gsl_rng* rng, apop_model *model){
    apop_mcmc_settings *s = apop_settings_get_group(model, apop_mcmc);
    if (!s || !s->pmf) {
        apop_model_metropolis(model->data, rng, model);
        s = apop_settings_get_group(model, apop_mcmc);
    }
    int constraint_fails = 0;
OMP_critical (metro_draw)
{
    apop_model *m = s->base_model;
    gsl_vector_view vv = gsl_vector_view_array(out, s->block_starts[s->block_count]);
    apop_data_pack(m->parameters, &(vv.vector));
    apop_data *earlier_draws = s->pmf->data;

    int block = 0, done = 0;
    while (!done){
        s->proposal_count++;
        earlier_draws->matrix = apop_matrix_realloc(earlier_draws->matrix, earlier_draws->matrix->size1+1, earlier_draws->matrix->size2);
        one_step(s->base_model->data, &(vv.vector), m, s, rng, &constraint_fails, 
                            earlier_draws, block, earlier_draws->matrix->size1-1);
        block = (block+1) % s->block_count;
        done = !block; //have looped back to the start.
        s->proposals[block].adapt_fn(s->proposals+block, s);
    }

    Apop_stopif(constraint_fails, , 2, "%i proposals failed to meet your model's parameter constraints", constraint_fails);
}
    return !!constraint_fails;
}


void main_mcmc_loop(apop_data *d, apop_model *m, apop_data *out, gsl_vector *draw, 
                        apop_mcmc_settings *s, gsl_rng *rng, int *constraint_fails){
    s->accept_count = 0;
    int block = 0;
    for (s->proposal_count=1; s->proposal_count< s->periods+1; s->proposal_count++){
        one_step(d, draw, m, s, rng, constraint_fails, out, block
                               , s->proposal_count-1 - s->periods*s->burnin);
        block = (block+1) % s->block_count;
        s->proposals[block].adapt_fn(s->proposals+block, s);
        //if (constraint_fails>10000) break;
    }
}

/** Use <a href="https://en.wikipedia.org/wiki/Metropolis-Hastings">Metropolis-Hastings
Markov chain Monte Carlo</a> to make draws from the given model.

The basic storyline is that draws are made from a proposal distribution, and the
likelihood of your model given your data and the drawn parameters evaluated. At each
step, a new set of proposal parameters are drawn, and if either they are more likely
than the previous set the new proposal is accepted as the next step, else with probability (prob of new params)/(prob of old params),
they are accepted as the next step anyway. Otherwise the last accepted proposal is repeated.

The output is an \ref apop_pmf model with a data set listing the draws that were
accepted, including those repetitions. The output model is modified so that subsequent
draws are one more step from the Markov chain, via \ref apop_model_metropolis_draw.

\li If a proposal fails to meet the \c constraint element of the model you input, then
the proposal is thrown out and a new one selected. By the default proposal
distribution, this is not mathematically correct (it breaks detailed balance),
and values near the constraint will be oversampled. The output model will have
<tt>outmodel->error=='c'</tt>. It is up to you to decide whether the resulting
distribution is good enough for your purposes or whether to take the time to write a
custom proposal and step function to accommodate the constraint.

Attach an \ref apop_mcmc_settings group to your model to specify the proposal
distribution, burnin, and other details of the search. See the \ref apop_mcmc_settings
documentation for details.

\li The default proposal includes an adaptive step: you specify a target accept rate
(default: .35), and if the accept rate is currently higher the variance of the proposals
is widened to explore more of the space; if the accept rate is currently lower the
variance is narrowed to stay closer to the last accepted proposal. Technically, this
breaks ergodicity of the Markov chain, but the consensus seems to be that this is
not a serious problem. If it does concern you, you can set the \c base_adapt_fn in the \ref apop_mcmc_settings group to a do-nothing function, or one that damps its adaptation as \f$n\to\infty\f$.

\li If you have a univariate model, \ref apop_arms_draw may be a suitable simpler alternative.

\li Note the \c gibbs_chunks element of the \ref apop_mcmc_settings group. If you set \c
gibbs_chunks='a', all parameters are drawn as a set, and accepted/rejected as a set. The
variances are adapted at an identical rate. If you set \c gibbs_chunks='i',
then each scalar parameter is assigned its own proposal distribution, which is adapted
at its own pace. With \c gibbs_chunks='b' (the default), then each of the vector, matrix,
and weights of your model's parameters are drawn/accepted/adapted as a block (and so
on to additional chunks if your model has <tt>->more</tt> pages). This works well for
complex models which naturally break down into subsets of parameters.

Each chunk counts as a step in the Markov chain. Therefore, if there are several chunks,
you can expect chunks to repeat from step to step. If you want a draw after cycling through all chunks, try using \ref apop_model_metropolis_draw, which has that behavior.

\param d The \ref apop_data set used for evaluating the likelihood of a proposed parameter set.

\param rng A \c gsl_rng, probably allocated via \ref apop_rng_alloc. (Default: an RNG from \ref apop_rng_get_thread)

\param m The \ref apop_model from which parameters are being drawn. (No default; must not be \c NULL)

\li If the likelihood model has \c NULL parameters, I will allocate them. That means you can use
one of the stock models that ship with Apophenia. If I need to run the model's prep
routine to get the size of the parameters, then I will make a copy of the likelihood
model, run prep, and then allocate parameters for that copy of a model.

\li On exit, the \c parameters element of your likelihood model has the last accepted parameter proposal.

\li If you set <tt>apop_opts.verbose=2</tt> or greater, I will report the accept rate of the M-H sampler. It is a common rule of thumb to select a proposal so that this is between 20% and 50%. Set <tt>apop_opts.verbose=3</tt> to see the stream of proposal points, their likelihoods, and the acceptance odds. You may want to set <tt>apop_opts.log_file=fopen("yourlog", "w")</tt> first.

\return A modified \ref apop_pmf model representing the results of the search. It has
a specialized \c draw method that returns another step from the Markov chain with each draw.

\exception out->error='c'  Proposal was outside of a constraint; see above.

\li This function uses the \ref designated syntax for inputs.
*/
APOP_VAR_HEAD apop_model *apop_model_metropolis(apop_data *d, gsl_rng *rng, apop_model *m){
    apop_data *apop_varad_var(d, NULL);
    apop_model *apop_varad_var(m, NULL);
    Apop_stopif(!m, return NULL, 0, "NULL model input.");
    gsl_rng *apop_varad_var(rng, apop_rng_get_thread());
APOP_VAR_END_HEAD
    apop_model *outp;
    OMP_critical(metropolis)
    {
    apop_mcmc_settings *s = apop_settings_get_group(m, apop_mcmc);
    if (!s)
        s = Apop_model_add_group(m, apop_mcmc);
    apop_prep(d, m); //typically a no-op
    s->last_ll = GSL_NEGINF;
    gsl_vector * drawv = apop_data_pack(m->parameters);
    Apop_stopif(s->burnin > 1, s->burnin/=(s->periods + 0.0), 
                1, "Burn-in should be a fraction of the number of periods, "
                   "not a whole number of periods. Rescaling to burnin=%g."
                   , s->burnin/(s->periods+0.0));
    apop_data *out = apop_data_alloc(s->periods*(1-s->burnin), drawv->size);

    if (!s->proposals){
        set_block_count_and_block_starts(m->parameters, s, drawv->size);
        s->proposals = calloc(s->block_count, sizeof(apop_mcmc_proposal_s));
        s->proposal_is_cp = 1;
        for (int i=0; i< s->block_count; i++){
            apop_mcmc_proposal_s *p = s->proposals+i;
            setup_normal_proposals(p, s->block_starts[i+1]-s->block_starts[i], s);
            if (!p->proposal->parameters) {
                apop_prep(NULL, p->proposal+i);
                if(p->proposal->parameters->matrix) gsl_matrix_set_all(p->proposal->parameters->matrix, 1);
                if(p->proposal->parameters->vector) gsl_vector_set_all(p->proposal->parameters->vector, 1);
            }
        }
    }

    //if s->start_at =='p', we already have m->parameters in drawv.
    if (s->start_at == '1') gsl_vector_set_all(drawv, 1);
    int constraint_fails = 0;

    main_mcmc_loop(d, m, out, drawv, s, rng, &constraint_fails);

    Apop_notify(2, "M-H sampling accept percent = %3.3f%%", 100*(0.0+s->accept_count)/s->periods);
    Apop_stopif(constraint_fails, out->error='c', 2, "%i proposals failed to meet your model's parameter constraints", constraint_fails);

    out->weights = gsl_vector_alloc(s->periods*(1-s->burnin));
    gsl_vector_set_all(out->weights, 1);
    outp = apop_estimate(out, apop_pmf);
    s->pmf = outp;
    s->base_model = m;
    outp->draw = apop_model_metropolis_draw;
    apop_settings_copy_group(outp, m, "apop_mcmc");

    gsl_vector_free(drawv);
    }
    return outp;
}
