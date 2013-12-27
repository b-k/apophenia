/* is_iid option

The algorithm to allocate observations to models does so by evaluating the log likelihood (LL) of a (model, single observation) pair, which makes sense iff observations are IID. In the typical case, when all input models are IID, leave this at the default of <tt>.is_iid='y'</tt>.

Assuming IID, there is no need for a second evaluation step: the total LL is simply the sum of the individual LLs. Setting this to <tt>.is_iid='n'</tt> asks the LL algorithm to generate a separate data set for each model (which involves copying each observation), and then calculating the LL for each (model, chosen observations) set. There will still be distortion, because individual observations are allocated via evaluations of the LL of only that observation. This is also problematic if the input model fails with a one-observation data set. The idealwould be for the algorithm to calculate the LL of all partitions of the data and select the most likely, but in naïve application this is clearly infeasible; readers are welcome to submit alternative algorithms. [Default: 'y']

 * weights element
 * faithful demo.
 */


/** \file 
 */
/*
\amodel apop_mixture The mixture model transformation: a linear combination of multiple models.  

Generated via \ref apop_model_mixture.

Note that a kernel density is a mixture of a large number of homogeneous models, where each is typically centered around a point in your data. For such situations, \ref apop_kernel_density will be easier to use.

One can think of the estimation as a missing-data problem: each data point originated
in one distribution or the other, and if we knew with certainty which data point came
from which distribution, then the estimation problem would be trivial: just generate the subsets and call 
<tt>apop_estimate(dataset1, model1)</tt>, ...,  <tt>apop_estimate(datasetn, modeln)</tt> separately.
But the assignment of which element goes where is unknown information, which we guess at using an E-M algorithm. The standard algorithm starts with an initial set of parameters for the models, and assigns each data point to its most likely model. It then re-estimates the model parameters using their subsets. The algorithm repeats until it arrives at an optimum.

Thus, the log likelihood method for this model includes a step that allocates each data point to its most likely model, and calculates the log likelihood of each observation using its most likely model. (If your model is not IID, see the notes in the \ref apop_mixture_settings struct.

Apohenia modifies this routine slightly because it uses the same maximum likelihood back-end that most other <tt>apop_model</tt>s use for estimation. The ML search algorithm provides model parameters, then the LL method allocates observations and reports a LL to the search algorithm, then the search algorithm uses its usual rules to step to the next candidate set of parameters. This provides slightly more flexibility in the ---hold on.


\li Mixture searches can be sensitive to initial conditions. You are encouraged to try a sequence random starting points for your model parameters.

\li You could potentially provide a list of unparameterized models, and then have \ref apop_estimate estimate the parameters of all of them at once:

\code
    apop_model *m_blank = apop_model_mixture(apop_model_copy(apop_normal), apop_model_copy(apop_exponential));
    apop_model *m_estimated = apop_estimate(your_data, m);
\endcode

We're still experimenting with better algorithms.

In the first maximization step, I estimate the parameters of each model by sending all
available data to each model's native <tt>estimate</tt> routine.  In the expectation step, I calculate the likelihood that the data
point is drawn from each of the distributions (given the parameters to this point),
and then assign the data point to the most likely distribution. The maximization step
then re-estimates the parameters using the new selection of data points. The expectation
and maximization steps alternate until convergence.

Note that determining to which parts of a mixture to assign a data point is a
well-known hard problem, which is often not solvable--that information is basically lost. 

\li The weights are free parameters, to be estimated like any other. If you want to
fix them, then use \ref apop_model_fix_params to do so. E.g.

\code
    apop_model *m = apop_model_mixture(apop_model_copy(apop_normal), apop_model_copy(apop_exponential));
    apop_prep(NULL, m); //for some models, NULL will need to be an \ref apop_data set.
    gsl_vector_set_all(m->parameters->vector, NAN);
    gsl_vector weightsvector = gsl_vector_subvector(m->parameters->vector, 0, 2).vector;
    gsl_vector_set_all(&weightsvector, 1);
    apop_model *fixed =apop_model_fix_params(m->parameters, m);
\endcode

\adoc    Input_format   The same data gets sent to each of the submodels.
\adoc    Settings   \ref apop_mixture_settings 
\adoc    Parameter_format The parameters are broken out in a readable form in the
                          settings group, so your best bet is to use those. See the sample code for usage.<br> The <tt>parameter</tt> element is a single vector piling up all elements, beginning with the first $n-1$ weights, followed by an <tt>apop_data_pack</tt> of each model's parameters in sequence. Because all elements are in a single vector, one could run a maximum likelihood search for all components (including the weights) at once. Fortunately for parsing, the <tt>log_likehood</tt>, <tt>estimate</tt>, and other methods unpack this vector into its component parts for you.

\adoc   Examples
\include hills2.c
*/

#include "apop_internal.h"
#include "types.h"
void mixture_show(apop_model *m, FILE *out);


Apop_settings_copy(apop_mixture,
    out->cmf= in->cmf ? apop_model_copy(in->cmf): NULL;
)

Apop_settings_free(apop_mixture,
    apop_model_free(in->cmf);
    free(in->param_sizes);
) 

Apop_settings_init(apop_mixture, 
    Apop_varad_set(rng, apop_rng_alloc(apop_opts.rng_seed++))
)

apop_model *apop_model_mixture_base(apop_model **inlist){
    apop_model *out = apop_model_copy(apop_mixture);
    int count=0;
    for (apop_model **m = inlist; *m; m++) count++;
    Apop_model_add_group(out, apop_mixture, .model_list=inlist, 
            .model_count=count, .param_sizes=malloc(sizeof(int)*count));
    int dsize = inlist[0]->dsize;
    for (int i=1; i< count && dsize > -99; i++) if (inlist[i]->dsize != dsize) dsize = -100;
    if (dsize > -99) out->dsize = dsize;
    return out;
}

apop_data* allocate_to_data_sets(apop_data *d, apop_model *m, apop_data **outsets){
    //if the data is iid, the output here is a grid of log likelihoods.
    //if is_iid!='y', then outsets!=NULL, then copy data to its best choice.

    apop_mixture_settings *ms = Apop_settings_get_group(m, apop_mixture);
    static int arbitrary_counter;

    Get_vmsizes(d); //maxsize, vsize, msize1
    apop_data *out = apop_data_alloc(maxsize, ms->model_count);

    for (int i=0; i< maxsize; i++){
        Apop_data_row(d, i, onepoint);
        int best_model = 0;
        double best_val = -INFINITY;
        for (int j=0; j< ms->model_count; j++){
            double this_val = apop_log_likelihood(onepoint, ms->model_list[j]);
            apop_data_set(out, i, j, this_val);
            if (!outsets) continue;
            if (this_val > best_val){
                best_model = j;
                best_val = this_val;
            } else if (this_val == best_val && (i+arbitrary_counter)%2==1){//split data set amongst models with the same params.
                best_model = j;
                best_val = this_val;
            }
        }
        if (outsets){
            apop_data *pick = outsets[best_model]; //alias
            if (!pick) pick = outsets[best_model] = apop_data_copy(onepoint);
            else apop_data_stack(pick, onepoint, .inplace='y');
        }
    }
    arbitrary_counter++;
    return out;
}

/*
static void mixture_estimate(apop_data *d, apop_model *m){
    //The weights are a tally, while the parameters are found via ML search.
    //So fix the weights, do the search, then calculate the weights.
    apop_prep(d, m);
    gsl_vector_set_all(m->parameters->vector, NAN);
    int mcount= Apop_settings_get(m, apop_mixture, model_count);
    gsl_vector weights = gsl_vector_subvector(m->parameters->vector, 0, mcount-1).vector;
    gsl_vector_set_all(&weights, 0);

    apop_model *mf = apop_model_fix_params(m);
    apop_prep(d, mf);
    apop_maximum_likelihood(d, mf);
    apop_model_fix_params_get_base(mf); //re-fills m.
}*/

/*
    Apop_stopif(m->error, return, 0, "Trouble estimating the initial MLEs.");
    apop_mixture_settings *ms = Apop_settings_get_group(m, apop_mixture);

    apop_data **datasets = calloc(ms->model_count, sizeof(apop_data*));
    for (int ctr=0; ctr<100; ctr++){
        allocate_to_data_sets(d, m, datasets);
        for (int i=0; i< ms->model_count; i++)
            if (datasets[i]){
                apop_model *freeme=ms->model_list[i];
                ms->model_list[i] = apop_estimate(datasets[i], ms->model_list[i]);
                apop_model_free(freeme);
            }
        for (int i=0; i< ms->model_count; i++)
            apop_data_free(datasets[i]);
    }
}
*/

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

    if (!ms->weights){
        ms->weights = gsl_vector_alloc(ms->model_count);
        gsl_vector_set_all(ms->weights, 1./ms->model_count);
    }
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

#define weighted_sum(fn)                                                            \
    long double total=0;                                                            \
    if (model_in->parameters) {                                                     \
        long double total_weight = apop_vector_sum(ms->weights);                    \
        size_t i=0;                                                                 \
        for (apop_model **m = ms->model_list; *m; m++)                              \
            total += fn(d, *m) * gsl_vector_get(ms->weights, i++)/total_weight;     \
    } else {   /*assume equal weights.*/                                            \
        for (apop_model **m = ms->model_list; *m; m++)                              \
            total += fn(d, *m)/ms->model_count;                                     \
    }

/*
    ll ll ll
        I need p₁/Σp p₂/Σp p₃/Σp
        where p₁=exp(ll₁).

        Then p₁/Σp = exp(ll₁) - exp(llM)*(exp(ll₁-llM)+exp(ll₂-llM)+exp(ll₃-llM))


        which is equvalent to
                1/(1+(p₂+p₃))
        better: the subtract-out-the-max trick.
p₁/Σp = exp(ll1)/(exp(ll1)+exp(ll2)+exp(ll3))
let L be the max
p₁/Σp = exp(ll1)/

N^6 = N^2(N^(6-2))
*/

/* The trick to summing exponents: subtract the max:

let ll_M be the max LL. then
Σexp(ll) = exp(llM)*(exp(ll₁-llM)+exp(ll₂-llM)+exp(ll₃-llM))

One of the terms in the sum is exp(0)=1. The others are all less than one, and so we
are guaranteed no overflow. If any of them underflow, then that term must not have
been very important for the sum.
*/
long double sum_exp_vector(gsl_vector const *onerow){
    long double rowtotal = 0;
    double best = gsl_vector_max(onerow);
    for (int j=0; j<onerow->size; j++) rowtotal += exp(gsl_vector_get(onerow, j)-best);
    rowtotal *= exp(best);
    return rowtotal;
}

static long double mixture_log_likelihood(apop_data *d, apop_model *model_in){
    apop_mixture_settings *ms = Apop_settings_get_group(model_in, apop_mixture);
    Apop_stopif(!ms, model_in->error='p'; return GSL_NAN, 0, "No apop_mixture_settings group. "
                                              "Did you estimate this with apop_model_mixture?");
    if (model_in->parameters) unpack(model_in);
    apop_data **datasets = ms->is_iid=='y' ? calloc(ms->model_count, sizeof(apop_data*)) : NULL;
    apop_data *lls = allocate_to_data_sets(d, model_in, datasets);


    //reweight by last round's lambda 
    for (int i=0; i< lls->matrix->size2; i++){
        Apop_matrix_col(lls->matrix, i, onecol);
        gsl_vector_add_constant(onecol, gsl_vector_get(ms->weights, i));
    }

//OK, now we need the λs, and then the max ll for each guy.

    long double total_ll=0;
    gsl_vector *ps = gsl_vector_alloc(lls->matrix->size2);
    gsl_vector *cp = gsl_vector_alloc(lls->matrix->size2);
    for (int i=0; i< lls->matrix->size1; i++){
        Apop_matrix_row(lls->matrix, i, onerow);
        total_ll += gsl_vector_max(onerow);
        for (int j=0; j < onerow->size; j++){
            gsl_vector_memcpy(cp, onerow);
            gsl_vector_add_constant(cp, -gsl_vector_get(onerow, j));
            gsl_vector_set(ps, j, 1./sum_exp_vector(cp));
        }
        gsl_vector_memcpy(onerow, ps);

    //check that the row sums to one.
        assert(fabs(apop_sum(onerow) - 1) < 1e-3);
    }
    gsl_vector_free(cp);
    gsl_vector_free(ps);

    for (int i=0; i< lls->matrix->size2; i++){
        Apop_matrix_col(lls->matrix, i, onecol);
        gsl_vector_set(ms->weights, i, apop_sum(onecol)/lls->matrix->size1);
    }
    return total_ll;



    if (ms->is_iid=='y')
        for (int i=0; i< ms->model_count; i++){
            if (datasets[i]) total_ll += apop_log_likelihood(datasets[i], ms->model_list[i]);
    printf("\n%i (%g, %g):%g\t%zu", i, 
            apop_data_get(ms->model_list[i]->parameters, 0),
            apop_data_get(ms->model_list[i]->parameters, 1),
            datasets[i] ? apop_log_likelihood(datasets[i], ms->model_list[i]): NAN,
            datasets[i] ? datasets[i]->matrix->size1: 0);
    //apop_vector_print(ms->model_list[i]->parameters->vector);
            apop_data_free(datasets[i]);
        }
        free(datasets);

printf(" Σ:%Lg\n", total_ll);
#if 0
    if (model_in->parameters) {
        long double total_weight = apop_vector_sum(ms->weights);
        size_t i=0;
        for (apop_model **m = ms->model_list; *m; m++)
            total += apop_p(d, *m) * vforsum.data[i++]/total_weight;
    } else {   /*assume equal weights.*/
        for (apop_model **m = ms->model_list; *m; m++)
            total += apop_p(d, *m)/ms->model_count;
    }
#endif
    return total_ll;
}

static void mixture_draw (double *out, gsl_rng *r, apop_model *m){
    apop_mixture_settings *ms = Apop_settings_get_group(m, apop_mixture);
    if (!ms->cmf){
        ms->cmf = apop_model_copy(apop_pmf);
        ms->cmf->data = apop_data_alloc();
        ms->cmf->data->weights = apop_vector_copy(ms->weights);
        Apop_model_add_group(ms->cmf, apop_pmf, .draw_index='y');
    }
    double index; apop_draw(&index, r, ms->cmf);
    apop_draw(out, r, ms->model_list[(int)index]);
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


void mixture_show(apop_model *m, FILE *out){
    apop_mixture_settings *ms = Apop_settings_get_group(m, apop_mixture);
    if (m->parameters){
        fprintf(out, "Mixture of %i models, with weights: ", ms->model_count);
        apop_vector_print(ms->weights, .output_pipe=out);
    } else
        fprintf(out, "Mixture of %i models, with unspecified weights\n", ms->model_count);

    for (int i=0; i< ms->model_count; i++)
        apop_model_print(ms->model_list[i], out);
}

//score
//predict

apop_model *apop_mixture=&(apop_model){"Mixture of models", .prep=mixture_prep, 
    /*.estimate=mixture_estimate,*/ .constraint=mixture_constraint,
    .log_likelihood=mixture_log_likelihood, .cdf=mixture_cdf, .draw=mixture_draw };
