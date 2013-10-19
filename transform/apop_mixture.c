/** \file 
 */
/*
\amodel apop_mixture The mixture model transformation: a linear combination of multiple models.  

Generated via \ref apop_model_mixture.

Note that a kernel density is a mixture of a large number of homogeneous models, where each is typically centered around a point in your data. For such situations, \ref apop_kernel_density will be easier to use.

\li This is still in beta. Please keep an eye out for bugs, and expect the interface to change a little. 

\li You could potentially provide a list of unparameterized models, and then have \ref apop_estimate estimate the parameters of all of them at once:

\code
    apop_model *m_blank = apop_model_mixture(apop_model_copy(apop_normal), apop_model_copy(apop_exponential));
    apop_model *m_estimated = apop_estimate(your_data, m);
\endcode

One can think of the estimation as a missing-data problem: each data point originated
in one distribution or the other, and if we knew with certainty which data point came
from which distribution, then the estimation problem would be trivial: just generate the subsets and call 
<tt>apop_estimate(dataset1, model1)</tt>, ...,  <tt>apop_estimate(datasetn, modeln)</tt> separately.
But the assignment of which element goes where is unknown information, which we guess at using an E-M algorithm. 
We're still experimenting with better algorithms.

In the first maximization step, I estimate the parameters of each model by sending all
available data to each model's native <tt>estimate</tt> routine (which may or may not
be a maximization).  In the expectation step, I calculate the likelihood that the data
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
\adoc    Parameter_format You probably just want to look at the submodel parameters. 
The parent model's parameter list is entirely in <tt>yrmodel->parameters->vector</tt>. The
first N are the weights, and then we have each model's parameters packed in, in turn.
There are pack/unpack functions that run as needed. Given the estimated mixture model \c m over three models, the weights are:
\code
    gsl_vector weightsvector = gsl_vector_subvector(m->parameters->vector, 0, 3).vector;
\endcode

\adoc   Examples
\include hills2.c
*/

#include "apop_internal.h"
#include "types.h"
void mixture_show(apop_model *m, FILE *out);

typedef struct {
    apop_model **model_list;
    int *param_sizes;
    apop_model *cmf;
    int refct, model_count;
    gsl_rng *rng;
} apop_mixture_settings;/**< All of the elements of this struct should be considered private. Also, the mixture setup is in beta and these will likely change soon. */

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

void allocate_to_data_sets(apop_data *d, apop_model *m, apop_data **outsets){
    apop_model *odds = apop_model_copy(apop_pmf); //Set up a CDF for mixing each data point.
    Apop_model_add_group(odds, apop_pmf, .draw_index='y');
    odds->data = apop_data_alloc();
    apop_mixture_settings *ms = Apop_settings_get_group(m, apop_mixture);
    odds->data->weights = gsl_vector_alloc(ms->model_count);

    gsl_vector weightings = gsl_vector_subvector(m->parameters->vector, 0, ms->model_count).vector;
    gsl_vector_set_all(&weightings, 0);

    //apop_pmf_settings *ps = Apop_settings_get_group(odds, apop_pmf);
    Get_vmsizes(d); //maxsize, vsize, msize1
    for (int i=0; i< maxsize; i++){
        Apop_data_row(d, i, onepoint);
        for (int j=0; j< odds->data->weights->size; j++)
            gsl_vector_set(odds->data->weights, j,
                    apop_p(onepoint, ms->model_list[j]));
        //gsl_vector_free(ps->cmf); //This is what apop_model_clear is supposed to do.
        double chosen_model;
        apop_draw(&chosen_model, ms->rng, odds);
        apop_data *pick = outsets[(int)chosen_model]; //alias
        (*gsl_vector_ptr(&weightings, (int)chosen_model))++;
        if (!pick) outsets[(int)chosen_model] = apop_data_copy(onepoint);
        else {
            if (pick->vector) apop_vector_realloc(pick->vector, pick->vector->size+1);
            if (pick->matrix) apop_matrix_realloc(pick->matrix, pick->matrix->size1, pick->matrix->size2+1);
            if (*pick->textsize) apop_text_alloc(pick, pick->textsize[0], pick->textsize[1]+1);
            Get_vmsizes(pick); //maxsize again.
            apop_data_set_row(pick, onepoint, maxsize-1);
        }
    }
    //set weights parameters here.
    apop_data_free(odds->data);
    apop_model_free(odds);
}

static void mixture_estimate(apop_data *d, apop_model *m){
    apop_maximum_likelihood(d, m);
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

static void mixture_prep(apop_data * data, apop_model *model){
    apop_model_print_vtable_add(mixture_show, apop_mixture);
    apop_mixture_settings *ms = Apop_settings_get_group(model, apop_mixture);
    model->parameters=apop_data_alloc(ms->model_count);

    //scaffolding, to remove later
    gsl_vector_set_all(model->parameters->vector, 1./ms->model_count);

    int i=0;
    for (apop_model **m = ms->model_list; *m; m++){
        if (!(*m)->parameters) apop_prep(data, *m);
        gsl_vector *v = apop_data_pack((*m)->parameters);
        ms->param_sizes[i++] = v->size;
        apop_vector_stack(model->parameters->vector, v, .inplace='y');
        gsl_vector_free(v);
    }
    if (!model->dsize) model->dsize = (*ms->model_list)->dsize;
}

void unpack(apop_model *min){
    apop_mixture_settings *ms = Apop_settings_get_group(min, apop_mixture);
    int posn=ms->model_count, i=0;
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
        gsl_vector vforsum = gsl_vector_subvector(model_in->parameters->vector, 0, ms->model_count).vector; \
        long double total_weight = apop_vector_sum(&vforsum);                       \
        size_t i=0;                                                                 \
        for (apop_model **m = ms->model_list; *m; m++)                              \
            total += fn(d, *m) * vforsum.data[i++]/total_weight;                    \
    } else {   /*assume equal weights.*/                                            \
        for (apop_model **m = ms->model_list; *m; m++)                              \
            total += fn(d, *m)/ms->model_count;                                     \
    }


static long double mixture_log_likelihood(apop_data *d, apop_model *model_in){
    apop_mixture_settings *ms = Apop_settings_get_group(model_in, apop_mixture);
    Apop_stopif(!ms, model_in->error='p'; return GSL_NAN, 0, "No apop_mixture_settings group. "
                                              "Did you estimate this with apop_model_mixture?");
    unpack(model_in);
    long double total=0;
    if (model_in->parameters) {
        gsl_vector vforsum = gsl_vector_subvector(model_in->parameters->vector, 0, ms->model_count).vector;
        long double total_weight = apop_vector_sum(&vforsum);
        size_t i=0;
        for (apop_model **m = ms->model_list; *m; m++)
            total += apop_p(d, *m) * vforsum.data[i++]/total_weight;
    } else {   /*assume equal weights.*/
        for (apop_model **m = ms->model_list; *m; m++)
            total += apop_p(d, *m)/ms->model_count;
    }
    return log(total);
}

static void mixture_draw (double *out, gsl_rng *r, apop_model *m){
    apop_mixture_settings *ms = Apop_settings_get_group(m, apop_mixture);
    if (!ms->cmf){
        gsl_vector v, *vv=NULL;
        if (m->parameters)
            v = gsl_vector_subvector(m->parameters->vector, 0, ms->model_count).vector;
        else { //assume equiprobable models.
            vv = gsl_vector_alloc(ms->model_count);
            gsl_vector_set_all(vv, 1);
            v = *vv;
        }

        ms->cmf = apop_model_copy(apop_pmf);
        ms->cmf->data = apop_data_alloc();
        ms->cmf->data->weights = apop_vector_copy(&v);
        Apop_model_add_group(ms->cmf, apop_pmf, .draw_index='y');
        gsl_vector_free(vv);
    }
    double index; apop_draw(&index, r, ms->cmf);
    apop_draw(out, r, ms->model_list[(int)index]);
}

static long double weights_over_zero(apop_data *data, apop_model *m){
    apop_mixture_settings *ms = Apop_settings_get_group(m, apop_mixture);
    gsl_vector v = gsl_vector_subvector(m->parameters->vector, 0, ms->model_count).vector;
    return apop_linear_constraint(&v);
}

static long double mixture_cdf(apop_data *d, apop_model *model_in){
    Nullcheck_m(model_in, GSL_NAN)
    Nullcheck_d(d, GSL_NAN)
    apop_mixture_settings *ms = Apop_settings_get_group(model_in, apop_mixture);
    unpack(model_in);
    weighted_sum(apop_cdf);
    return total;
}

void mixture_show(apop_model *m, FILE *out){
    apop_mixture_settings *ms = Apop_settings_get_group(m, apop_mixture);
    if (m->parameters){
        gsl_vector weightings = gsl_vector_subvector(m->parameters->vector, 0, ms->model_count).vector;

        fprintf(out, "Mixture of %i models, with weights: ", ms->model_count);
        apop_vector_print(&weightings, .output_pipe=out);
    } else
        fprintf(out, "Mixture of %i models, with unspecified weights\n", ms->model_count);

    for (int i=0; i< ms->model_count; i++)
        apop_model_print(ms->model_list[i], out);
}

//score
//predict

apop_model *apop_mixture=&(apop_model){"Mixture of models", .prep=mixture_prep, 
    .estimate=mixture_estimate, .constraint=weights_over_zero,
    .log_likelihood=mixture_log_likelihood, .cdf=mixture_cdf, .draw=mixture_draw };
