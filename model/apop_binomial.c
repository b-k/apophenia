/* The binomial distribution as an \c apop_model.
Copyright (c) 2006--2007, 2010 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

/* \amodel apop_binomial The multi-draw generalization of the Bernoulli; the two-bin special case of the Multinomial.

 \adoc  Name  <tt> Binomial distribution</tt>

 \adoc Input_format Zeros are failures and non-zeros successes. \f$N\f$
is the size of the matrix, vector, or both (whichever is not \c NULL).
So \f$p\f$ represents the odds of a success==1; the odds of a zero is \f$1-p\f$.

\li You may be interested in \ref apop_data_to_factors to convert real numbers or text into a
vector of categories.

\li See also \ref apop_data_rank_compress for means of dealing with one more input data format.

 \adoc  Parameter_format
        The parameters are kept in the vector element of the \c apop_model parameters element. \c parameters->vector->data[0]==n;
        \c parameters->vector->data[1...]==p_1....

The numeraire is zero, meaning that \f$p_0\f$ is not explicitly listed, but is
\f$p_0=1-\sum_{i=1}^{k-1} p_i\f$, where \f$k\f$ is the number of bins. Conveniently enough,
the zeroth element of the parameters vector holds \f$n\f$, and so a full probability vector can
easily be produced by overwriting that first element. Continuing the above example: 
\code 
int n = apop_data_get(estimated->parameters, 0, -1); 
apop_data_set(estimated->parameters, 0, 1 - (apop_sum(estimated->parameters)-n)); 
\endcode
And now the parameter vector is a proper list of probabilities.

\adoc    RNG I fill an array of length \c n, with a sequence of randomly drawn ones and zeros. 
*/

#include "model.h"
#include "mapply.h"
#include "output.h"
#include "internal.h"
#include "likelihoods.h"

static double is_over_zero(double in){return in > 0;}

static void get_hits_and_misses(apop_data *data, double *hitcount, double *misscount){
    Get_vmsizes(data); //vsize, msize1, msize2;
    *hitcount = apop_map_sum(data, .fn_d=is_over_zero, .part='a');
    *misscount = vsize + msize1 * msize2 - *hitcount;
}

static void make_covar(apop_model *est){
    int size = est->parameters->vector->size;
    //the trick where we turn the params into a p-vector
    double * pv = est->parameters->vector->data;
    int n = pv[0];
    pv[0] = 1 - (apop_sum(est->parameters->vector)-n);

    apop_data *cov = apop_data_alloc(0, size, size);
    for (int i=0; i < size; i++){
        double p = apop_data_get(est->parameters, i, -1);
        apop_data_set(cov, i, i, n * p *(1-p));
        for (int j=i+1; j < size; j++){
            double pj = apop_data_get(est->parameters, j, -1);
            apop_data_set(cov, i, j, -n*p*pj);
            apop_data_set(cov, j, i, -n*p*pj);
        }
    }
    apop_data_add_page(est->parameters, cov, "<Covariance>");
    pv[0]=n;
}

static double binomial_log_likelihood(apop_data *d, apop_model *params){
  Nullcheck_m(params) Nullcheck_p(params) Nullcheck_d(d)
  double	  n       = apop_data_get(params->parameters, 0, -1),
              p       = apop_data_get(params->parameters, 1, -1);
  double hitcount, misscount;
    get_hits_and_misses(d, &hitcount, &misscount);
    return log(gsl_ran_binomial_pdf(hitcount, p, n));
}

/* \adoc estimated_parameters  As per the parameter format. Has a <tt>\<Covariance\></tt> page with the covariance matrix for the \f$p\f$s (\f$n\f$ effectively has no variance).  */
/* \adoc estimated_info   Reports <tt>log likelihood</tt>. */
static apop_model * binomial_estimate(apop_data * data,  apop_model *est){
  Nullcheck_d(data)
  double hitcount, misscount;
    get_hits_and_misses(data, &hitcount, &misscount);   
    int n = hitcount + misscount;
    apop_name_add(est->parameters->names, "n", 'r');
    apop_name_add(est->parameters->names, "p", 'r');
    apop_data_set(est->parameters, 0, -1, n);
    apop_data_set(est->parameters, 1, -1, hitcount/(hitcount + misscount));
    est->dsize = n;
    apop_data_add_named_elmt(est->info, "log likelihood", binomial_log_likelihood(data, est));
    make_covar(est);
    return est;
}

static double binomial_cdf(apop_data *d, apop_model *est){
  Nullcheck_m(est) Nullcheck_p(est) Nullcheck_d(d)
  double hitcount, misscount, psum = 0;
    get_hits_and_misses(d, &hitcount, &misscount);   
    double n = gsl_vector_get(est->parameters->vector, 0);
    double p = gsl_vector_get(est->parameters->vector, 1);
    for (int i=0; i<= hitcount; i++)
        psum += gsl_ran_binomial_pdf(hitcount, p, n);
    return psum;
}

static double multinomial_constraint(apop_data *data, apop_model *b){
  //constraint is that 0 < all elmts 
    return apop_linear_constraint(b->parameters->vector, .margin = 1e-3);
}

static void binomial_rng(double *out, gsl_rng *r, apop_model* est){
  Nullcheck_m(est); Nullcheck_p(est);
  double n = gsl_vector_get(est->parameters->vector, 0);
  double p = gsl_vector_get(est->parameters->vector, 1);
  *out = gsl_ran_binomial_knuth(r, p, n);
  /*  for (int i=0; i < n; i++)  //naive version. Knuth first uses a beta approximation, then finishes off with this.
        out[i] = (gsl_rng_uniform(r) <= p) ? 1 : 0; //one Bernoulli draw.
   */
}

static gsl_vector * get_multinomial_hitcount(const apop_data *data){
    size_t     i, j;
    gsl_vector *out;
    out = gsl_vector_calloc(1+GSL_MAX(data->vector ? gsl_vector_max(data->vector) : 0,
                                   data->matrix ? gsl_matrix_max(data->matrix) : 0));
    if (data->vector)
        for(i=0; i < data->vector->size; i ++)
            (*gsl_vector_ptr(out, apop_data_get(data, i, -1)))++;
    if (data->matrix)
        for(i=0; i < data->matrix->size1; i ++)
            for(j=0; j < data->matrix->size2; j ++)
                (*gsl_vector_ptr(out, apop_data_get(data, i, j)))++;
    return out;
}

static double multinomial_log_likelihood(apop_data *d, apop_model *params){
    Nullcheck(params); Nullcheck_p(params);
    double *pv = params->parameters->vector->data;
    size_t size = params->parameters->vector->size;

    //The GSL wants our hit count in an int*.
    gsl_vector *hits = get_multinomial_hitcount(d);
    unsigned int *hv = malloc(hits->size * sizeof(unsigned int));
    for(size_t i=0; i < hits->size; i ++)
        hv[i] = hits->data[i];
    gsl_vector_free(hits);

    double n = pv[0]; //making the params a p-vector. Put n back at the end.
    pv[0] = 1 - (apop_sum(params->parameters->vector)-n);
    double out =  gsl_ran_multinomial_lnpdf(size, pv, hv);

    pv[0]=n;
    free(hv);
    return out;
}

/*
static apop_model *multinomial_paramdist(apop_data *d, apop_model *m){
    apop_pm_settings *settings = Apop_settings_get_group(m, apop_pm);
    if (settings->index!=-1){
        int i = settings->index;
        double mu = apop_data_get(m->parameters, i, -1);
        double sigma = sqrt(apop_data_get(m->parameters, i, i, .page="<Covariance>"));
        int df = apop_data_get(m->info, .rowname="df", .page = "info");
        return apop_model_set_parameters(apop_t_distribution, mu, sigma, df);
    }

}
*/

static void multinomial_rng(double *out, gsl_rng *r, apop_model* est){
    //After the intro, cut/pasted/modded from the GSL. Copyright them.
    Nullcheck_p(est);
    double * p = est->parameters->vector->data;
    //the trick where we turn the params into a p-vector
    int N = p[0];
    p[0] = 1 - (apop_sum(est->parameters->vector)-N);
    double sum_p = 0.0;
    int draw, ctr = 0;
    unsigned int sum_n = 0;

    for (int i = 0; sum_n < N; i++) {
        if (p[i] > 0)
            draw = gsl_ran_binomial (r, p[i] / (1 - sum_p), N - sum_n);
        else
            draw = 0;
        for (int j=0; j< draw; j++)
            out[ctr++] = i;
        sum_p += p[i];
        sum_n += draw;
    }
    p[0] = N;
}

static void multinomial_show(apop_model *est){
    double * p = est->parameters->vector->data;
    int N=p[0];
    p[0] = 1 - (apop_sum(est->parameters->vector)-N);
    fprintf(apop_opts.output_pipe, "%s, with %i draws.\nBin odds:\n", est->name, N);
    apop_vector_print(est->parameters->vector, .output_pipe=apop_opts.output_pipe);
}

apop_model apop_binomial = {"Binomial distribution", 2,0,0, .dsize=1,
	.estimate = binomial_estimate, .log_likelihood = binomial_log_likelihood, 
   .constraint = multinomial_constraint, .draw = binomial_rng, .cdf =binomial_cdf,
   .print=multinomial_show};


/* \amodel apop_multinomial The \f$n\f$--option generalization of the \ref apop_binomial "Binomial distribution".
    See also the \ref apop_binomial model. 

\adoc estimated_parameters  As per the parameter format. Has a <tt>\<Covariance\></tt> page with the covariance matrix for the \f$p\f$s (\f$n\f$ effectively has no variance).  */
/* \adoc estimated_info   Reports <tt>log likelihood</tt>. */
static apop_model * multinomial_estimate(apop_data * data,  apop_model *est){
    Get_vmsizes(data); //vsize, msize1
    Nullcheck(est);
    gsl_vector * count = get_multinomial_hitcount(data);
    int n = vsize + msize1; //size of one row
/*    int n = 0;
    for (int i=0; i< count->size; i++)
        n += gsl_vector_get(count, i); */
    apop_vector_normalize(count);
    gsl_vector_set(count, 0, n);
    est->parameters=apop_data_alloc();
    est->parameters->vector = count;
    apop_name_add(est->parameters->names, "n", 'r');
    char name[100];
    for(int i=1; i < count->size; i ++){
        sprintf(name, "p%i", i);
        apop_name_add(est->parameters->names, name, 'r');
    }
    est->dsize = n;
    make_covar(est);
    apop_data_add_named_elmt(est->info, "log likelihood", multinomial_log_likelihood(data, est));
    return est;
}


/* \adoc    Input_format The default is simply a listing of bins, without regard to whether items are in the vector or
matrix of the \ref apop_data struct, or the dimensions. Here, data like <tt>0, 1, 2, 1, 1</tt>
represents one draw of zero, three draws of 1, and one draw of 2.

\li You may be interested in \ref apop_data_to_factors to convert real numbers or text into a
vector of categories.

\li See also \ref apop_data_rank_compress for means of dealing with one more input data format.

\li Please note that the number of bins is simply the largest number found. So if there
are bins {0, 1, 2} and your data set happens to consist of <tt>0 0 1 1 0</tt>, then
I won't know to generate results with three bins where the last bin has probability zero.

\adoc    Parameter_format
        The parameters are kept in the vector element of the \c apop_model parameters element. \c parameters->vector->data[0]==n;
        \c parameters->vector->data[1...]==p_1....

The numeraire is zero, meaning that \f$p_0\f$ is not explicitly listed, but is
\f$p_0=1-\sum_{i=1}^{k-1} p_i\f$, where \f$k\f$ is the number of bins. Conveniently enough,
the zeroth element of the parameters vector holds \f$n\f$, and so a full probability vector can
easily be produced by overwriting that first element. Continuing the above example: 
\code 
int n = apop_data_get(estimated->parameters, 0, -1); 
apop_data_set(estimated->parameters, 0, 1 - (apop_sum(estimated->parameters)-n)); 
\endcode
And now the parameter vector is a proper list of probabilities.

\li Because an observation is typically a single row, the value of \f$N\f$ is set to equal the length of
the first row (counting both vector and matrix elements, as appropriate). Thus, if your
data is entirely in the vector or a one-column matrix, then the \f$p\f$s are estimated
using all data, but \f$N=1\f$. The covariances are calculated accordingly, and a random
draw would return a single bin. 

\adoc    Estimate_results  Parameters are estimated. Covariance matrix is filled.   
\adoc    RNG The result of an imaginary tossing of \f$N\f$ balls into \f$k\f$ urns, with the
            given probabilities.
            
            I fill an array of length \c N, with a sequence of draws from zero to \f$N\f$. They
            are not randomly ordered: it'll look something like \f$[0 0 1 1 3 3 3]\f$, but
            will still be an accurate representation of what happens when you throw
            \f$N\f$ balls into \f$k\f$ urns and sum the results. 
            
            If you want the sequence of draws to be random at the per-item scale,
            set \f$N=1\f$ (i.e., <tt>apop_data_set(estimated->parameters, 0, 1);</tt>),
            and use a \c for loop to make the number of draws you want. This is less efficient.
            */

apop_model apop_multinomial = {"Multinomial distribution", -1,0,0, .dsize=-1,
	.estimate = multinomial_estimate, .log_likelihood = multinomial_log_likelihood, 
   .constraint = multinomial_constraint, .draw = multinomial_rng, .print=multinomial_show};
