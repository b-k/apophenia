/* The binomial distribution as an \c apop_model.
Copyright (c) 2006--2007, 2010--11 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2. 

 \amodel apop_binomial The multi-draw generalization of the Bernoulli; the two-bin special case of the \ref apop_multinomial "Multinomial distribution".
This differs from the \ref apop_multinomial only in the input data format.

It is implemented as an alias of the \ref apop_multinomial model, except that it has
a CDF, <tt>.vsize==2</tt> and <tt>.dsize==1</tt> (i.e., we know it has two parameters
and a draw returns a scalar).

\adoc    Parameter_format   a vector, v[0]=\f$n\f$; v[1]=\f$p_1\f$. Thus, \f$p_0\f$
        isn't written down; see \ref apop_multinomial for further discussion.
        If you input $v[1]>1$ and <tt>apop_opts.verbose >=1</tt>, the log likelihood
        function will throw a warning.

\adoc    Input_format Each row of the matrix is one observation, consisting of two elements.
  The number of draws of type zero (sometimes read as `misses' or `failures') are in column zero, 
  the number of draws of type one (`hits', `successes') in column one.

\adoc    RNG The RNG returns a single number representing the success count, not a
    vector of length two giving both the failure bin and success bin. This is notable
    because it differs from the input data format, but it tends to be what people expect
    from a Binomial RNG. For draws with both dimensions, use a \ref apop_multinomial model
    with <tt>.vsize =2</tt>.
*/

#include "apop_internal.h"

/* \adoc cdf At the moment, only implemented for the Binomial.
  Let the first element of the data set (top of the vector or point (0,0) in the
  matrix, your pick) be $L$; then I return the sum of the odds of a draw from the given
  Binomial distribution returning $0, 1, \dots, L$ hits.  */
static long double binomial_cdf(apop_data *d, apop_model *est){
    Nullcheck_mpd(d, est, GSL_NAN)
    Get_vmsizes(d); //firstcol
    double hitcount = apop_data_get(d, .col=firstcol);
    double n = gsl_vector_get(est->parameters->vector, 0);
    double p = gsl_vector_get(est->parameters->vector, 1);
    return gsl_cdf_binomial_P(hitcount, p, n);
}

static void make_covar(apop_model *est){
    int size = est->parameters->vector->size;
    //the trick where we turn the params into a p-vector
    double * pv = est->parameters->vector->data;
    int n = pv[0];
    pv[0] = 1 - (apop_sum(est->parameters->vector)-n);

    apop_data *cov = apop_data_add_page(est->parameters, 
                            apop_data_alloc(size, size), "<Covariance>");
    for (int i=0; i < size; i++){
        double p = apop_data_get(est->parameters, i, -1);
        apop_data_set(cov, i, i, n * p *(1-p));
        for (int j=i+1; j < size; j++){
            double pj = apop_data_get(est->parameters, j, -1);
            double thiscell = -n*p*pj;
            apop_data_set(cov, i, j, thiscell);
            apop_data_set(cov, j, i, thiscell);
        }
    }
    pv[0]=n;
}

static long double multinomial_constraint(apop_data *data, apop_model *b){
  //constraint is that 0 < all elmts, and  1>all ps.
    int size = b->parameters->vector->size;
    static threadlocal apop_data *constr;
    if (constr && constr->matrix->size2 != size)
        apop_data_free(constr);
    if (!constr){
        constr = apop_data_calloc(size*2-1, size*2-1, size);

        //top half: 0 < [param], including param 0
        Apop_submatrix(constr->matrix, 0, 0, size, size, tophalf);
        gsl_matrix_set_identity(tophalf);

        //bottom (almost) half: 1 >= [param], excluding param 0
        for (int i=size; i < size*2-1; i++){
            apop_data_set(constr, i, -1, -1);
            apop_data_set(constr, i, i-size+1, -1);
        }
    }
    return apop_linear_constraint(b->parameters->vector, constr);
}

static double binomial_ll(gsl_vector *hits, void *paramv){
    return log(gsl_ran_binomial_pdf(hits->data[1], ((gsl_vector*)paramv)->data[1], ((gsl_vector*)paramv)->data[0]));
}

static double multinomial_ll(gsl_vector *v, void *params){
    double *pv = ((apop_model*)params)->parameters->vector->data;
    size_t size = ((apop_model*)params)->parameters->vector->size;
    unsigned int hv[v->size]; //The GSL wants our hit count in an int*.
    for (size_t i=0; i < v->size; i ++)
        hv[i] = gsl_vector_get(v, i);
    return gsl_ran_multinomial_lnpdf(size, pv, hv);
}

static long double multinomial_log_likelihood(apop_data *d, apop_model *params){
    Nullcheck_mpd(d, params, GSL_NAN);
    double *pv = params->parameters->vector->data;
    double n = pv[0]; 
    Apop_assert_c(params->parameters->vector->size>=2, GSL_NAN, 0, "I need two or more input parameters "
                    "representing [n, p_1, (...)].");
    Apop_assert_c(pv[1] <=1, GSL_NAN, 1, "The input parameters should be [n, p_1, (...)], but "
        "element 1 of the parameter vector is >1."); //mostly makes sense for the binomial.
    if (n==2) return apop_map_sum(d, .fn_vp=binomial_ll, .param=params->parameters->vector);

    pv[0] = 1 - (apop_sum(params->parameters->vector)-n);//making the params a p-vector. Put n back at the end.
    double out = apop_map_sum(d, .fn_vp=multinomial_ll, .param=params);
    pv[0]=n;
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

static int multinomial_rng(double *out, gsl_rng *r, apop_model* est){
    Nullcheck_mp(est, 1);
    double * p = est->parameters->vector->data;
    //the trick where we turn the params into a p-vector
    int N = p[0];

    if (est->parameters->vector->size == 2) {
        *out = gsl_ran_binomial_knuth(r, 1-gsl_vector_get(est->parameters->vector, 1), N);
        out[1] = N-*out;
        goto done;
    }
    //else, multinomial
    //cut/pasted/modded from the GSL. Copyright them.
    p[0] = 1 - (apop_sum(est->parameters->vector)-N);
    double sum_p = 0.0;
    int sum_n = 0;

    for (int i = 0; i < est->parameters->vector->size; i++) {
        out[i] = (p[i] > 0)
                ? gsl_ran_binomial (r, p[i] / (1 - sum_p), N - sum_n)
                : 0;
        sum_p += p[i];
        sum_n += out[i];
    }
    done:
    p[0] = N;
    return 0;
}

static void multinomial_show(apop_model *est, FILE *out){
    double * p = est->parameters->vector->data;
    int N=p[0];
    p[0] = 1 - (apop_sum(est->parameters->vector)-N);
    fprintf(out, "%s, with %i draws.\nBin odds:\n", est->name, N);
    apop_vector_print(est->parameters->vector, .output_pipe=out);
    p[0]=N;
}

double avs(gsl_vector *v){return (double) apop_vector_sum(v);}

/* \amodel apop_multinomial The \f$n\f$--option generalization of the \ref apop_binomial "Binomial distribution".

\adoc estimated_parameters  As per the parameter format. Has a <tt>\<Covariance\></tt> page with the covariance matrix for the \f$p\f$s (\f$n\f$ effectively has no variance).  */
/* \adoc estimated_info   Reports <tt>log likelihood</tt>. */
static void multinomial_estimate(apop_data * data,  apop_model *est){
    Nullcheck_mpd(data, est, );
    Get_vmsizes(data); //vsize, msize1
    est->parameters= apop_map(data, .fn_v=avs, .part='c');
    gsl_vector *v = est->parameters->vector;
    int n = apop_sum(v)/data->matrix->size1; //size of one row
    apop_vector_normalize(v);
    apop_name_add(est->parameters->names, "n", 'r');
    apop_data_set(est->parameters, .val=n); //zeroth item is now n, not p_0
    char name[100];
    for(int i=1; i < v->size; i++){
        sprintf(name, "p%i", i);
        apop_name_add(est->parameters->names, name, 'r');
    }
    est->dsize = n;
    make_covar(est);
    apop_data_add_named_elmt(est->info, "log likelihood", multinomial_log_likelihood(data, est));
}

static void multinom_prep(apop_data *data, apop_model *params){
    apop_model_print_vtable_add(multinomial_show, params);
    apop_model_clear(data, params);
}

/* \adoc    Input_format Each row of the matrix is one observation: a set of draws from a single bin.
  The number of draws of type zero are in column zero, the number of draws of type one in column one, et cetera.

\li You may have a set of several Bernoulli-type draws, which could be summed together to form a single Binomial draw.
See \ref apop_data_to_dummies to do the aggregation (using the <tt>.keep_first='y'</tt> option).

\adoc    Parameter_format
        The parameters are kept in the vector element of the \c apop_model parameters element. \c parameters->vector->data[0]==n;
        \c parameters->vector->data[1...]==p_1....

The numeraire is bin zero, meaning that \f$p_0\f$ is not explicitly listed, but is
\f$p_0=1-\sum_{i=1}^{k-1} p_i\f$, where \f$k\f$ is the number of bins. Conveniently enough,
the zeroth element of the parameters vector holds \f$n\f$, and so a full probability vector can
easily be produced by overwriting that first element. Continuing the above example: 
\code 
int n = apop_data_get(estimated->parameters, 0, -1); 
apop_data_set(estimated->parameters, 0, 1 - (apop_sum(estimated->parameters)-n)); 
\endcode
And now the parameter vector is a proper list of probabilities.

\li Because an observation is a single row, the number of bins, \f$k\f$ is set to equal
the length of the first row (counting both vector and matrix elements, as appropriate).
The covariance matrix will be \f$k \times k\f$.

\li Each row should sum to \f$N\f$, the number of draws. The estimation routine doesn't check this, but instead uses the average sum across all rows.

\adoc    Estimate_results  Parameters are estimated. Covariance matrix is filled.   
\adoc    RNG Returns a single vector of length \f$k\f$, the result of an imaginary tossing 
        of \f$N\f$ balls into \f$k\f$ urns, with the given probabilities.
            */

apop_model *apop_multinomial = &(apop_model){"Multinomial distribution", -1, .dsize=-1,
	.estimate = multinomial_estimate, .log_likelihood = multinomial_log_likelihood, 
   .constraint = multinomial_constraint, .draw = multinomial_rng, .prep=multinom_prep};

apop_model *apop_binomial = &(apop_model){"Binomial distribution", 2, .dsize=1,
	.estimate = multinomial_estimate, .log_likelihood = multinomial_log_likelihood, 
   .constraint = multinomial_constraint, .draw = multinomial_rng, .prep=multinom_prep, .cdf= binomial_cdf};
