/** \file 
 */
/*
\amodel apop_mixture The mixture model transformation: a linear combination of multiple models.  

Generated via \ref apop_model_mixture.

Note that a kernel density is a mixture of a large number of homogeneous models, where each is typically centered around a point in your data. For such situations, \ref apop_kernel_density will be easier to use.

Use \ref apop_model_mixture to produce one of these models. In the examples below, some are generated from unparameterized input models with a form like 

\code
apop_model *mf = apop_model_mixture(apop_model_copy(apop_normal), apop_model_copy(apop_normal));
Apop_model_add_group(mf, apop_mle, .starting_pt=(double[]){50, 5, 80, 5},
                                   .step_size=3, .tolerance=1e-6);
apop_model_show(apop_estimate(dd, mf));
\endcode

One can also skip the estimation and use already-parameterized models as input to \ref apop_model_mixture, e.g.:

\code
apop_model *r_ed = apop_model_mixture(apop_model_set_parameters(apop_normal, 54.6, 5.87),
                       apop_model_set_parameters(apop_normal, 80.1, 5.87));
apop_data *wts = apop_data_falloc((2), 0.36, 0.64);
Apop_settings_add(r_ed, apop_mixture, weights, wts->vector);
printf("LL=%g\n", apop_log_likelihood(dd, r_ed));
\endcode

Notice that the weights vector has to be added after the model is set up. If none is given, then equal weights are assigned to all components of the mixture.


One can think of the estimation as a missing-data problem: each data point originated
in one distribution or the other, and if we knew with certainty which data point
came from which distribution, then the estimation problem would be trivial:
just generate the subsets and call <tt>apop_estimate(dataset1, model1)</tt>, ...,
<tt>apop_estimate(datasetn, modeln)</tt> separately.  But the assignment of which
element goes where is unknown information, which we guess at using an E-M algorithm. The
standard algorithm starts with an initial set of parameters for the models, and assigns
each data point to its most likely model. It then re-estimates the
model parameters using their subsets. The standard algorithm, see e.g. <a
href="http://www.jstatsoft.org/v32/i06/paper">this PDF</a>, repeats until it arrives
at an optimum.

Thus, the log likelihood method for this model includes a step that allocates each data
point to its most likely model, and calculates the log likelihood of each observation
using its most likely model. [It would be a valuable extension to extend this to
not-conditionally IID models. Commit \c 1ac0dd44 in the repository had some notes on
this, now removed.]  As a side-effect, it calculates the odds of drawing from each model
(the vector λ). Following the above-linked paper, the probability for a given
observation under the mixture model is its probability under the most likely model
weighted by the previously calculated λ for the given model.

Apohenia modifies this routine slightly because it uses the same maximum likelihood
back-end that most other <tt>apop_model</tt>s use for estimation. The ML search algorithm
provides model parameters, then the LL method allocates observations and reports a LL to
the search algorithm, then the search algorithm uses its usual rules to step to the next
candidate set of parameters. This provides slightly more flexibility in the search.

\li <b>Estimations of mixture distributions can be sensitive to initial conditions.</b>
You are encouraged to try a sequence random starting points for your model parameters.
Some authors recommend plotting the data and eyeballing a guess as to the model parameters.

Determining to which parts of a mixture to assign a data point is a
well-known hard problem, which is often not solvable--that information is basically lost. 

\adoc    Input_format   The same data gets sent to each of the component models of the
mixture. Each row is an observation, and the estimation routine assumes that models are
conditionally IID (i.e., having chosen what component of the mixture the observation
comes from, its likelihood can be calculated independently of all other observations).

\adoc    Settings   \ref apop_mixture_settings 

\adoc    Parameter_format The parameters are broken out in a readable form in the
    settings group, so your best bet is to use those. See the sample code for usage.<br>
    The <tt>parameter</tt> element is a single vector piling up all elements, beginning
    with the first $n-1$ weights, followed by an <tt>apop_data_pack</tt> of each model's
    parameters in sequence. Because all elements are in a single vector, one could run a
    maximum likelihood search for all components (including the weights) at once. Fortunately
    for parsing, the <tt>log_likehood</tt>, <tt>estimate</tt>, and other methods unpack
    this vector into its component parts for you.

\adoc RNG Uses the weights to select a component model, then makes a draw from that component.
The model's \c dsize (draw size) element is set when you set up the model in the
model's \c prep method (automatically called by \ref apop_estimate, or call it directly)
iff all component models have the same \c dsize.

\adoc   Examples
The first example uses a text file \c faith.data, in the \c tests directory of the distribution.
\include faithful.c

\include hills2.c
*/

#include "apop_internal.h"

Apop_settings_copy(apop_mixture,
    (*out->cmf_refct)++;
)

Apop_settings_free(apop_mixture,
    if (!(--in->cmf_refct)) {
        apop_model_free(in->cmf);
        free(in->cmf_refct);
    }
    free(in->param_sizes);
) 

Apop_settings_init(apop_mixture, 
    out->cmf_refct = calloc(1, sizeof(int));
    (*out->cmf_refct)++;
)

//see apop_model_mixture in types.h
apop_model *apop_model_mixture_base(apop_model **inlist){
    apop_model *out = apop_model_copy(apop_mixture);
    int count=0;
    for (apop_model **m = inlist; *m; m++) count++;
    apop_mixture_settings *ms =Apop_settings_add_group(out, apop_mixture, .model_list=inlist, 
            .model_count=count, .param_sizes=malloc(sizeof(int)*count));
    int dsize = inlist[0]->dsize;
    for (int i=1; i< count && dsize > -99; i++) if (inlist[i]->dsize != dsize) dsize = -100;
    if (dsize > -99) out->dsize = dsize;

    ms->weights = gsl_vector_alloc(ms->model_count);
    gsl_vector_set_all(ms->weights, 1./count);
    return out;
}

void mixture_show(apop_model *m, FILE *out){
    apop_mixture_settings *ms = Apop_settings_get_group(m, apop_mixture);
    if (m->parameters){
        fprintf(out, "Mixture of %i models, with weights:\n", ms->model_count);
        apop_vector_print(ms->weights, .output_pipe=out);
    } else fprintf(out, "Mixture of %i models, with unspecified weights\n", ms->model_count);

    for (int i=0; i< ms->model_count; i++){
        fprintf(out, "\n");
        apop_model_print(ms->model_list[i], out);
    }
}


static void mixture_prep(apop_data * data, apop_model *model){
    apop_model_print_vtable_add(mixture_show, apop_mixture);
    if (model->parameters) return;
    apop_mixture_settings *ms = Apop_settings_get_group(model, apop_mixture);
    model->parameters = apop_data_alloc();

    int i=0;
    for (apop_model **m = ms->model_list; *m; m++){
        if (!(*m)->parameters) apop_prep(data, *m);
        gsl_vector *v = apop_data_pack((*m)->parameters);
        ms->param_sizes[i++] = v->size;
        model->parameters->vector = apop_vector_stack(model->parameters->vector, v, .inplace='y');
        gsl_vector_free(v);
    }
    if (!model->dsize) model->dsize = (*ms->model_list)->dsize;
}

void unpack(apop_model *min){
    apop_mixture_settings *ms = Apop_settings_get_group(min, apop_mixture);
    int posn=0, i=0;
    if (!min->parameters) return; //Trusting user that the user has added already-esimated models.
    for (apop_model **m = ms->model_list; *m; m++){
        gsl_vector v = gsl_vector_subvector(min->parameters->vector, posn, ms->param_sizes[i]).vector;
        apop_data_unpack(&v, (*m)->parameters);
        posn+=ms->param_sizes[i++];
    }
}

#define weighted_sum(fn)                                                     \
    long double total=0;                                                     \
    long double total_weight = apop_sum(ms->weights);                        \
    size_t i=0;                                                              \
    for (apop_model **m = ms->model_list; *m; m++)                           \
        total += fn(d, *m) * gsl_vector_get(ms->weights, i++)/total_weight;

//The output is a grid of log likelihoods.
apop_data* get_lls(apop_data *d, apop_model *m){
    apop_mixture_settings *ms = Apop_settings_get_group(m, apop_mixture);
    Get_vmsizes(d); //maxsize
    apop_data *out = apop_data_alloc(maxsize, ms->model_count);

    for (int i=0; i< maxsize; i++){
        Apop_row(d, i, onepoint);
        for (int j=0; j< ms->model_count; j++){
            double this_val = apop_log_likelihood(onepoint, ms->model_list[j]);
            apop_data_set(out, i, j, this_val);
        }
    }
    return out;
}

/* The trick to summing exponents: subtract the max:

let ll_M be the max LL. then
Σexp(ll) = exp(llM)*(exp(ll₁-llM)+exp(ll₂-llM)+exp(ll₃-llM))

One of the terms in the sum is exp(0)=1. The others are all less than one, and so we
are guaranteed no overflow. If any of them underflow, then that term must not have
been very important for the sum.
*/
static long double sum_exp_vector(gsl_vector const *onerow){
    long double rowtotal = 0;
    double best = gsl_vector_max(onerow);
    for (int j=0; j<onerow->size; j++) rowtotal += exp(gsl_vector_get(onerow, j)-best);
    rowtotal *= exp(best);
    return rowtotal;
}

static long double mixture_log_likelihood(apop_data *d, apop_model *model_in){
    apop_mixture_settings *ms = Apop_settings_get_group(model_in, apop_mixture);
    Apop_stopif(!ms, model_in->error='p'; return GSL_NAN, 0, "No apop_mixture_settings group. "
                                              "Did you set this up with apop_model_mixture()?");
    if (model_in->parameters) unpack(model_in);
    apop_data *lls = get_lls(d, model_in);

    //reweight by last round's lambda 
    for (int i=0; i< lls->matrix->size2; i++){
        Apop_col_v(lls, i, onecol);
        gsl_vector_add_constant(onecol, gsl_vector_get(ms->weights, i));
    }

//OK, now we need the λs, and then the max ll for each observation

/*
Draw probabilities are p₁/Σp p₂/Σp p₃/Σp (see equation (2) of the above pdf.)
But I have logs, and want to stay in log-format for as long as possible, to prevent undeflows and loss of precision.
*/
    long double total_ll=0;
    gsl_vector *ps = gsl_vector_alloc(lls->matrix->size2);
    gsl_vector *cp = gsl_vector_alloc(lls->matrix->size2);
    for (int i=0; i< lls->matrix->size1; i++){
        Apop_row_v(lls, i, onerow);
        total_ll += gsl_vector_max(onerow);
        for (int j=0; j < onerow->size; j++){
            gsl_vector_memcpy(cp, onerow);
            gsl_vector_add_constant(cp, -gsl_vector_get(onerow, j));
            gsl_vector_set(ps, j, 1./sum_exp_vector(cp));
        }
        gsl_vector_memcpy(onerow, ps);

        Apop_stopif(fabs(apop_sum(onerow) - 1) > 1e-3, /*Warn user, but continue.*/, 0,
                "One of the probability calculations is off: the total for the odds of drawing "
                "from the %i mixtures is %g but should be 1.", 
                ms->model_count, fabs(apop_sum(onerow) - 1));
    }
    gsl_vector_free(cp);
    gsl_vector_free(ps);

    for (int i=0; i< lls->matrix->size2; i++){
        Apop_col_v(lls, i, onecol);
        gsl_vector_set(ms->weights, i, apop_sum(onecol)/lls->matrix->size1);
    }
    return total_ll;
}

static int mixture_draw (double *out, gsl_rng *r, apop_model *m){
    apop_mixture_settings *ms = Apop_settings_get_group(m, apop_mixture);
    OMP_critical (mixdraw)
    if (!ms->cmf){
        ms->cmf = apop_model_copy(apop_pmf);
        ms->cmf->data = apop_data_alloc();
        ms->cmf->data->weights = apop_vector_copy(ms->weights);
        Apop_model_add_group(ms->cmf, apop_pmf, .draw_index='y');
    }
    double index; 
    Apop_stopif(apop_draw(&index, r, ms->cmf), return 1, 
            0, "Couldn't select a mixture element using the internal PMF over mixture elements.");
    return apop_draw(out, r, ms->model_list[(int)index]);
}

static long double mixture_cdf(apop_data *d, apop_model *model_in){
    Nullcheck_m(model_in, GSL_NAN)
    Nullcheck_d(d, GSL_NAN)
    apop_mixture_settings *ms = Apop_settings_get_group(model_in, apop_mixture);
    unpack(model_in);
    weighted_sum(apop_cdf);
    return total;
}

static long double mixture_constraint(apop_data *data, apop_model *model_in){
    apop_mixture_settings *ms = Apop_settings_get_group(model_in, apop_mixture);
    long double penalty = 0;
    unpack(model_in);
    //check all component models.
    for (apop_model **m = ms->model_list; *m; m++)
        penalty += (*m)->constraint(data, *m);
    if (penalty){
        int posn=0, i=0;
        for (apop_model **m = ms->model_list; *m; m++){
            gsl_vector v = gsl_vector_subvector(model_in->parameters->vector, posn, ms->param_sizes[i]).vector;
            apop_data_pack((*m)->parameters, &v);
            posn+=ms->param_sizes[i++];
        }
    }
    //weights are all positive?
    gsl_vector v = gsl_vector_subvector(model_in->parameters->vector, 0, ms->model_count).vector;
    penalty += apop_linear_constraint(&v);
    return penalty;
}

apop_model *apop_mixture=&(apop_model){"Mixture of models", .prep=mixture_prep,
    .constraint=mixture_constraint, .log_likelihood=mixture_log_likelihood,
    .cdf=mixture_cdf, .draw=mixture_draw };
