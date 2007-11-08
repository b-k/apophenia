/** \file apop_mle.c	The MLE functions. Call them with an \ref apop_model.

This file includes a number of distributions and models whose parameters
one would estimate using maximum likelihood techniques.

Each typically includes four functions: the likelihood function, the 
derivative of the likelihood function, a function that calls both of them,
and a user-usable function which takes in data and a blank vector, fills 
the vector with the most likely parameters, and returns the likelihood
of those parameters.

At the bottom are the maximum likelihood procedures themselves. There
are two: the no-derivative version and the with-derivative version.
Use the with-derivative version wherever possible---in fact, it is at
the moment entirely unused, but is just here for future use.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
*/
#include "likelihoods.h"
#include <assert.h>
#include <setjmp.h>
#include <signal.h>
#include <gsl/gsl_deriv.h>

//in apop_regress.c:
void apop_estimate_parameter_t_tests (apop_model *est);


/** \page trace_path Plotting the path of an ML estimation.

If \c trace_path (in the \c apop_mle_settings struct) has a name of positive
length, then every time the MLE evaluates the function, then the value
will be output to a table in the database/a file with the given name
(depending on the value of \c apop_opts.output_type). You can then plot this
table to get an idea of the path the estimation routine used to arrive
at its MLE.

To write to a pipe or stdout, set \c apop_opts.output_type appropriately and set \c trace_path to the literal string \c "NULL".


Below is a sample of the sort of output one would get:<br>
\image latex "search.gif" "An ML search, tracing out the surface of the function" width=\textwidth
\image html "search.gif" "An ML search, tracing out the surface of the function" 

\ingroup mle
*/

/** Use this if you already have an \c apop_model struct, but want to set the \c method_settings element to the default MLE parameters. 
 Returns the \c apop_mle_settings pointer, but the argument now has the
 method_settings element set, so you can ignore the returned pointer if
 you prefer.

 \param parent  A pointer to an allocated \c apop_model struct.
 \return A pointer to a set-up \c apop_mle_settings struct. The parent's \c method_settings element points to this struct as well.

 */
apop_mle_settings *apop_mle_settings_set_default(apop_model *parent){
  apop_mle_settings *setme =   calloc(1,sizeof(apop_mle_settings));
    setme->starting_pt      = NULL;
    setme->tolerance        = 1e-2;
    setme->method           = 1;
    setme->verbose          = 0;
    setme->use_score        = 1;
    setme->step_size        = 0.05;
    setme->delta            = 1e-2;
    setme->want_cov         = 1;
//siman:
    //siman also uses step_size  = 1.;  
    setme->n_tries          = 200; 
    setme->iters_fixed_T    = 200; 
    setme->k                = 1.0;
    setme->t_initial        = 50;  
    setme->mu_t             = 1.002; 
    setme->t_min            = 5.0e-1;
    setme->rng              = NULL;
    setme->model            = parent;
    parent->method_settings   = setme;
    parent->method_settings_size   = sizeof(apop_mle_settings);
    strcpy(setme->model->method_name, "MLE");
    setme->trace_path[0]    = '\0';
    return setme;
}

/** Neatly allocate an \ref apop_mle_settings structure. Sets a
few defaults, so you can change just one or two values and everything
else will be predictable.
When you finally call your MLE, use the \c model element of this 

\ingroup mle
 */
apop_mle_settings *apop_mle_settings_alloc(apop_data * data, apop_model model){
  apop_mle_settings *setme = apop_mle_settings_set_default(apop_model_copy(model));
    apop_model_clear(data, setme->model);
    return setme;
}

		///////////////////////
		//MLE support functions
		///////////////////////

typedef	void 	(*apop_df_with_void)(const gsl_vector *beta, void *d, gsl_vector *gradient);
typedef	void 	(*apop_fdf_with_void)(const gsl_vector *beta, void *d, double *f, gsl_vector *df);

//Including numerical differentiation and a couple of functions to
//negate the likelihood fns without bothering the user.

typedef struct {
	gsl_vector	*beta;
	apop_data	*d;
	int		    dimension;
} grad_params;

typedef struct {
    apop_model  *model;
    apop_data   *data;
    apop_fn_with_params   *f;
    apop_df_with_void   *df;
    grad_params *gp;
    gsl_vector  *beta;
    gsl_vector  *starting_pt;
    int         use_constraint;
    gsl_matrix  ***gradient_list;
    double      **gradientp;
    gsl_vector  ***score_list;
    size_t      *gsize;
    char        *trace_path;
    FILE        **trace_file;
}   infostruct;

static apop_model * apop_annealing(infostruct*);                         //below.
static int dnegshell (const gsl_vector *, void * , gsl_vector * g); //below.

static double one_d(double b, void *in){
  infostruct    *i   =in;
  int           vsize           = (i->model->parameters->vector? i->model->parameters->vector->size:0),
                msize1          = (i->model->parameters->matrix? i->model->parameters->matrix->size1:0),
                msize2          = (i->model->parameters->matrix? i->model->parameters->matrix->size2:0);
    gsl_vector_set(i->gp->beta, i->gp->dimension, b);
    apop_data_free(i->model->parameters);
    i->model->parameters    =  apop_data_unpack(i->gp->beta, vsize, msize1, msize2);
	double out= (*(i->f))(i->gp->d, i->model);
    return out;
}


#include "apop_findzeros.c"

static void apop_internal_numerical_gradient(apop_fn_with_params ll, infostruct* info, gsl_vector *out){
  int		    j;
  gsl_function	F;
  double		result, err;
  grad_params 	gp;
  infostruct    i;
  gsl_vector    *beta   = info->model->parameters->vector;
  apop_mle_settings   *mp = info->model->method_settings;
    memcpy(&i, info, sizeof(i));
    i.f         = &ll;
	gp.beta		= gsl_vector_alloc(beta->size);
	gp.d		= info->data;
    i.gp        = &gp;
	F.function	= one_d;
	F.params	= &i;
	for (j=0; j< beta->size; j++){
		gp.dimension	= j;
		gsl_vector_memcpy(gp.beta, beta);
		gsl_deriv_central(&F, gsl_vector_get(beta,j), mp->delta, &result, &err);
		gsl_vector_set(out, j, result);
	}
}

/**The GSL provides one-dimensional numerical differentiation; here's the multidimensional extension.

 If \c m has \c model_settings of type \c apop_ml_params, then the \c delta element is used for the differential.
 
 \code
 gradient = apop_numerical_gradient(beta, data, your_model);
 \endcode

 \ingroup linear_algebra
 */
gsl_vector * apop_numerical_gradient(apop_data *data, apop_model *m){
  infostruct    i;
  apop_fn_with_params ll  = m->log_likelihood ? m->log_likelihood : m->p;
    if (!ll){
        apop_error(0, 'c', "%s: Input model has neither p nor log_likelihood method. Returning zero.\n");
        return 0;
    }
  gsl_vector        *out= gsl_vector_alloc(m->parameters->vector->size);
  apop_mle_settings *mp;
  int clean = 0;
    if(strcmp(m->method_name, "MLE")){
        mp       = apop_mle_settings_alloc(data, *m);
        i.model  = mp->model;
        clean ++;
    } else 
        i.model  = m;
    i.data  = data;
    apop_internal_numerical_gradient(ll, &i, out);
    if (clean) {
        free(mp);
        apop_model_free(i.model);
    }
    return out;
}

/* They always tell you to just negate your likelihood function to turn
a minimization routine into a maximization routine---and this is the
sort of annoying little detail that Apophenia is intended to take care
of for you. The next few functions do the negation, so you have one
less sign that you have to remember. 

The negshell and dnegshell fns also take care of checking constraints,
if any, and are otherwise the primary point of contact between the
apop_models' methods and the GSL's MLE routines.
*/

static void tracepath(const gsl_vector *beta, double out, char tp[], FILE **tf){
    if (apop_opts.output_type == 'd'){
        if (beta->size == 1){
            if(!apop_table_exists(tp, 0))
                apop_query("create table  %s (beta0, ll);", tp);
            apop_query("insert into %s values (%g, %g);", tp, gsl_vector_get(beta,0), out);
        } else {
            if(!apop_table_exists(tp, 0))
                apop_query("create table  %s (beta0, beta1, ll);", tp);
            apop_query("insert into %s values (%g, %g, %g);", tp, gsl_vector_get(beta,0), gsl_vector_get(beta,1), out);
        }
    } else if (apop_opts.output_type == 'p'){
        if (beta->size == 1)
            fprintf(apop_opts.output_pipe, "%g\t %g\n",  gsl_vector_get(beta,0), out);
        else
            fprintf(apop_opts.output_pipe, "%g\t %g\t %g\n",  gsl_vector_get(beta,0), gsl_vector_get(beta,1), out);
    } else {
        if (!*tf){
            if (apop_opts.output_type == 's' && !strcmp(tp, "NULL"))
                *tf = stdout;
            else
                if (!(*tf = fopen(tp, "a"))) {
                    apop_error(0, 'c', "%s: couldn't open %s for writing. Continuing without the path trace.\n", __func__, tp);
                    return;
                }
        }
        if (beta->size == 1)
            fprintf(*tf, "%g\t %g\n",  gsl_vector_get(beta,0), out);
        else
            fprintf(*tf, "%g\t %g\t %g\n",  gsl_vector_get(beta,0), gsl_vector_get(beta,1), out);
    }
}

static void apop_pack_in_place(apop_data *from, gsl_vector *to){
    gsl_vector *v = apop_data_pack(from);
    gsl_vector_memcpy(to, v);
    gsl_vector_free(v);
}

static double negshell (const gsl_vector *beta, void * in){
  infostruct    *i              = in;
  double		penalty         = 0,
                out             = 0; 
  double 	(*f)(const apop_data *, apop_model *);
  int           vsize           = (i->model->parameters->vector? i->model->parameters->vector->size:0),
                msize1          = (i->model->parameters->matrix? i->model->parameters->matrix->size1:0),
                msize2          = (i->model->parameters->matrix? i->model->parameters->matrix->size2:0);
    f   = i->model->log_likelihood? i->model->log_likelihood : i->model->p;
    if (!f)
        apop_error(0, 's', "The model you sent to the MLE function has neither log_likelihood element nor p element.\n");
    apop_data_free(i->model->parameters);
    i->model->parameters  = apop_data_unpack(beta, vsize, msize1, msize2);
	if (i->use_constraint && i->model->constraint)
		penalty	= i->model->constraint(i->data, i->model);
    out = penalty - f(i->data, i->model); //negative llikelihood
    if (penalty)
        apop_pack_in_place(i->model->parameters, i->beta);
    if (strlen(i->trace_path))
        tracepath(i->model->parameters->vector,-out, i->trace_path, i->trace_file);
    i->model->llikelihood = -out; //negative negative llikelihood.
    return out;
}

static int dnegshell (const gsl_vector *beta, void * in, gsl_vector * g){
/* The derivative-calculating routine.
If the constraint binds
    then: take the numerical derivative of negshell, which will be the
    numerical derivative of the penalty.
    else: just find dlog_likelihood. If the model doesn't have a
    dlog likelihood or the user asked to ignore it, then the main
    maximum likelihood fn replaced model.score with
    apop_numerical_gradient anyway.
Finally, reverse the sign, since the GSL is trying to minimize instead of maximize.
*/
  infostruct    *i              = in;
  apop_mle_settings *mp       = i->model->method_settings;
  int           vsize           = (i->model->parameters->vector? i->model->parameters->vector->size:0),
                msize1          = (i->model->parameters->matrix? i->model->parameters->matrix->size1:0),
                msize2          = (i->model->parameters->matrix? i->model->parameters->matrix->size2:0);
    apop_data_free(i->model->parameters);
    i->model->parameters  = apop_data_unpack(beta, vsize, msize1, msize2);
    if(i->model->constraint)
        i->model->constraint(i->data, i->model);
    if (mp->use_score && i->model->score)
        i->model->score(i->data, g, i->model);
    else {
        apop_fn_with_params ll  = i->model->log_likelihood ? i->model->log_likelihood : i->model->p;
        apop_internal_numerical_gradient(ll, i, g);
    }
    if (strlen(i->trace_path))
        negshell (beta,  in);
    gsl_vector_scale(g, -1);
    return GSL_SUCCESS;
}

static void record_gradient_info(double f, gsl_vector *g, infostruct *i){
//By the way, the awkward pointer-to-pointer for the gradient list and
//plist is because the simulated annealing algorithm breaks if
//i->gradient_list or i->gradientp itself changes.
    if (!gsl_finite(f) || !gsl_finite(apop_vector_sum(g))) return;
    *i->score_list    = realloc(*i->score_list, sizeof(gsl_vector*)*(*i->gsize +1));
    (*i->score_list)[*i->gsize] = apop_vector_copy(g);
    (*i->gradientp)            = realloc(*i->gradientp, sizeof(double)*(*i->gsize + 1));
    (*i->gradientp)[*i->gsize] = i->model->log_likelihood? -f : -log(f);
    (*i->gsize)               ++;
}

static void fdf_shell(const gsl_vector *beta, void *i, double *f, gsl_vector *df){

	//if (negshell_model.fdf==NULL){
		*f	= negshell(beta, i);
		dnegshell(beta, i, df);
	/*} else	{
		negshell_model.fdf(beta, d, f, df);
		(*f) 	*= -1;
		gsl_vector_scale(df, -1);
	}*/
}

/** Calculate the Hessian.

  This is a synonym for \ref apop_numerical_second_derivative, q.v.
gsl_matrix * apop_numerical_hessian(apop_model dist, gsl_vector *beta, apop_data * d, void *params){
	return apop_numerical_second_derivative(dist, beta, d, params);
}
*/

/** Feeling lazy? Rather than doing actual pencil-and-paper math to find
your variance-covariance matrix, just use the negative inverse of the Hessian.

\param dist	The model
\param est	The estimate, with the parameters already calculated. The var/covar matrix will be placed in est->covariance.
\param data	The data
\ingroup basic_stats
*/
void apop_numerical_covariance_matrix(apop_model *est, apop_data *data){
    //As you can see, this is a placeholder.
    return;
}

static void cov_cleanup(infostruct *i){
  /*int j;
    for (j=0; j< *i->gsize; j++)
        gsl_vector_free((*i->score_list)[j]);
    free(*i->score_list);
    free(i->score_list);*/
    free(*i->gradientp);
    free(i->gsize);
    free(i->gradientp);
}

/** Transform the list of probabilities to a set of weightings where all
 ps are divided by the total.

 If the list of prob.s is from i->model->p, this is trivial.
 If the list of prob.s is from i->model->log_likelihood, then we don't
 want to exponentiate---that defeats the purpose of taking logs to
 begin with.  Instead, rescale against the median of the set, then
 exponentiate. The hope is that any given value is within a factor of
 maybe 100 of the median, so exp(ln(beta_n)-ln(median)) is a reasonable
 value.
 */
static void get_weightings(infostruct *i){
  int           j;
  long double   psum            = 0;
  long double   scaling_factor  = 0;
    if (i->model->log_likelihood){
        gsl_vector  *sortme = apop_array_to_vector(*i->gradientp,  *i->gsize);
        double      *pctiles= apop_vector_percentiles(sortme, 'u');
        double      median  = pctiles[50];
        gsl_vector_free(sortme);
        free(pctiles);
        for (j=0; j< *i->gsize; j++)
            scaling_factor  +=
            (*i->gradientp)[j] = exp((*i->gradientp)[j] - median);
        for (j=0; j< *i->gsize; j++)
            (*i->gradientp)[j] /= scaling_factor;
    } else {
        for (j=0; j< *i->gsize; j++)
            psum   += (*i->gradientp)[j];
        for (j=0; j< *i->gsize; j++)
            (*i->gradientp)[j]    /= psum;
    }
}


typedef struct {
    size_t  len;
    gsl_vector *score;
    gsl_vector **score_list;
} adhocstruct;

void produce_covariance_matrix(apop_model * est, infostruct *i){
  apop_mle_settings *parm = i->model->method_settings;
  if (!parm->want_cov)
      goto cov_clean;
  int        m, p, q;
  size_t    betasize    = (*i->score_list)[0]->size;
  gsl_matrix *preinv    = gsl_matrix_calloc(betasize, betasize);
  adhocstruct *out      = malloc(sizeof(adhocstruct));
    out->len        = *i->gsize;
    out->score      = gsl_vector_calloc(betasize);
    out->score_list = *i->score_list;
    get_weightings(i);
    //inv (n E(score dot score)) = info matrix.
    for (m=0; m<  *i->gsize; m++)
        if (gsl_finite((*i->gradientp)[m])){
            gsl_vector_scale((*i->score_list)[m], (*i->gradientp)[m]);
            gsl_vector_add(out->score,(*i->score_list)[m]);
        }
    for (p=0; p < betasize; p ++)
        for (q=p; q < betasize; q ++){
            long double total   = 0;
            for (m=0; m< *i->gsize; m++)
                if (gsl_finite((*i->gradientp)[m]))
                    total += gsl_vector_get((*i->score_list)[m], p) * gsl_vector_get((*i->score_list)[m], q);
            total   *= i->data ? (i->data->matrix ? i->data->matrix->size1 
                        : i->data->vector->size) : 1;
            gsl_matrix_set(preinv, p, q, total);
            gsl_matrix_set(preinv, q, p, total);
        }
    if (apop_opts.verbose > 1){
        printf("The expected Hessian:\n");
        apop_matrix_show(preinv);
    }
    gsl_matrix *inv = apop_matrix_inverse(preinv);
    est->covariance = apop_matrix_to_data(inv);
    if (est->parameters->names->row){
        apop_name_stack(est->covariance->names, est->parameters->names, 'r');
        apop_name_cross_stack(est->covariance->names, est->parameters->names, 'c', 'r');
    }
    gsl_matrix_free(preinv);
    est->more   = out;
  cov_clean:
    cov_cleanup(i);
}

static int ctrl_c;
static void mle_sigint(){ ctrl_c ++; }

/* The maximum likelihood calculations, given a derivative of the log likelihood.

If no derivative exists, will calculate a numerical gradient.

\param data	the data matrix
\param	dist	the \ref apop_model object: waring, probit, zipf, &amp;c.
\param	starting_pt	an array of doubles suggesting a starting point. If NULL, use a vector whose elements are all 0.1 (zero has too many pathological cases).
\param step_size	the initial step size.
\param tolerance	the precision the minimizer uses. Only vaguely related to the precision of the actual var.
\param verbose		Y'know.
\return	an \ref apop_model with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.

  \todo readd names */
static apop_model *	apop_maximum_likelihood_w_d(apop_data * data, infostruct *i){
  gsl_multimin_function_fdf minme;
  gsl_multimin_fdfminimizer *s;
  gsl_vector 			    *x;
  apop_model		        *est    = i->model; //just an alias.
  int                       vsize   =(est->parameters->vector ? est->parameters->vector->size :0),
                            msize1  =(est->parameters->matrix ? est->parameters->matrix->size1 :0),
                            msize2  =(est->parameters->matrix ? est->parameters->matrix->size2:0);
  int				        iter 	= 0, 
				            status  = 0,
				            betasize= vsize+ msize1*msize2;
  apop_mle_settings           *mp     = est->method_settings;
    if (mp->method == APOP_CG_BFGS)
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, betasize);
    else if (mp->method == APOP_CG_PR)
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_pr, betasize);
    else //Default:    APOP_CG_FR      conjugate gradient (Fletcher-Reeves)
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, betasize);
	if (mp->starting_pt==NULL){
		x	= gsl_vector_alloc(betasize);
  		gsl_vector_set_all (x,  0.1);
	}
	else 	x   = apop_array_to_vector(mp->starting_pt, betasize);
	minme.f		= negshell;
	minme.df	= (apop_df_with_void) dnegshell;
	minme.fdf	= (apop_fdf_with_void) fdf_shell;
	minme.n		= betasize;
	minme.params	= i;
    i->beta     = s->x;
    ctrl_c      =
    est->status = 0;
	gsl_multimin_fdfminimizer_set (s, &minme, x, mp->step_size, mp->tolerance);
    signal(SIGINT, mle_sigint);
    do { 	
        iter++;
        status 	= gsl_multimin_fdfminimizer_iterate(s);
        record_gradient_info(gsl_multimin_fdfminimizer_minimum(s), gsl_multimin_fdfminimizer_gradient(s), i);
        if (status) 	break; 
        status = gsl_multimin_test_gradient(s->gradient,  mp->tolerance);
        if (mp->verbose)
            printf ("%5i %.5f  f()=%10.5f gradient=%.3f\n", iter, gsl_vector_get (s->x, 0),  s->f, gsl_vector_get(s->gradient,0));
        if (status == GSL_SUCCESS){
            est->status	= 1;
            if(mp->verbose)	printf ("Minimum found.\n");
        }
    } while (status == GSL_CONTINUE && iter < MAX_ITERATIONS_w_d && !ctrl_c);
    signal(SIGINT, NULL);
	if (iter==MAX_ITERATIONS_w_d) {
		est->status	= -1;
		if (mp->verbose) printf("No min!!\n");
	}
	//Clean up, copy results to output estimate.
    //est->parameters = apop_data_unpack(s->x, vsize, msize1, msize2);
	gsl_multimin_fdfminimizer_free(s);
	if (mp->starting_pt==NULL) 
		gsl_vector_free(x);
    produce_covariance_matrix(est, i);
    apop_estimate_parameter_t_tests (est);
	return est;
}

static apop_model *	apop_maximum_likelihood_no_d(apop_data * data, infostruct * i){
  apop_model		        *est        = i->model;
  apop_mle_settings           *mp         = est->method_settings;
  apop_data                 *p          = est->parameters;
  int                       vsize       = (p->vector ? p->vector->size :0),
                            msize1      = (p->matrix ? p->matrix->size1:0),
                            msize2      = (p->matrix ? p->matrix->size2:0);
  int			            status,
			                iter 		= 0,
			                betasize	= vsize + msize1 * msize2; 
  size_t 			        j;
  gsl_multimin_function 	minme;
  gsl_multimin_fminimizer   *s;
  gsl_vector 		        *x, *ss;
  double			        size;
	s	= gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, betasize);
	ss	= gsl_vector_alloc(betasize);
	x	= gsl_vector_alloc(betasize);
    ctrl_c      =
	est->status	= 0;	//assume failure until we score a success.
	if (mp->starting_pt==NULL)
  		gsl_vector_set_all (x,  0);
	else
		x   = apop_array_to_vector(mp->starting_pt, betasize);
    /*if (!msize1 && !msize2)
        i->beta = i->model->parameters->vector;
    else*/
        i->beta = gsl_vector_alloc(betasize);
  	gsl_vector_set_all (ss,  mp->step_size);
	minme.f		    = negshell;
	minme.n		    = betasize;
	minme.params	= i;
	gsl_multimin_fminimizer_set (s, &minme, x,  ss);
    signal(SIGINT, mle_sigint);
      	do { 	iter++;
		status 	= gsl_multimin_fminimizer_iterate(s);
		if (status) 	break; 
		size	= gsl_multimin_fminimizer_size(s);
	   	status 	= gsl_multimin_test_size (size, mp->tolerance); 
		if(mp->verbose){
			printf ("%5d ", iter);
			for (j = 0; j < betasize; j++) {
				printf ("%8.3e ", gsl_vector_get (s->x, j)); } 
			printf ("f()=%7.3f size=%.3f\n", s->fval, size);
       			if (status == GSL_SUCCESS) {
                printf ("Optimum found at:\n");
                printf ("%5d ", iter);
                for (j = 0; j < betasize; j++) {
                    printf ("%8.3e ", gsl_vector_get (s->x, j)); } 
                printf ("f()=%7.3f size=%.3f\n", s->fval, size);
			}
		}
      	} while (status == GSL_CONTINUE && iter < MAX_ITERATIONS && !ctrl_c);
    signal(SIGINT, NULL);
	if (iter == MAX_ITERATIONS && mp->verbose)
		apop_error(1, 'c', "Optimization reached maximum number of iterations.");
    if (status == GSL_SUCCESS) 
        est->status	= 1;
    est->parameters = apop_data_unpack(s->x, vsize, msize1, msize2);
	gsl_multimin_fminimizer_free(s);
	if (mp->want_cov) 
		apop_numerical_covariance_matrix(est, data);
	return est;
}


/** The maximum likelihood calculations

\param data	The data matrix (an \ref apop_data set).
\param	dist	The \ref apop_model object: waring, probit, zipf, &amp;c. This can be allocated using 
\c apop_mle_settings_alloc (not quite: see that page), featuring:<br>
starting_pt:	an array of doubles suggesting a starting point. If NULL, use zero.<br>
step_size:	the initial step size.<br>
tolerance:	the precision the minimizer uses. Only vaguely related to the precision of the actual var.<br>
verbose:	Y'know.<br>
method:		
\li    APOP_SIMPLEX_NM      Nelder-Mead simplex (gradient handling rule is irrelevant)
\li    APOP_CG_FR      conjugate gradient (Fletcher-Reeves) (default)
\li    APOP_CG_BFGS    conjugate gradient (BFGS: Broyden-Fletcher-Goldfarb-Shanno)
\li    APOP_CG_PR      conjugate gradient (Polak-Ribiere)
\li    APOP_SIMAN       \ref simanneal "simulated annealing"
\li    APOP_RF_NEWTON   Find a root of the derivative via Newton's method
\li    APOP_RF_BROYDEN  Find a root of the derivative via the Broyden Algorithm
\li    APOP_RF_HYBRID   Find a root of the derivative via the Hybrid method
\li    APOP_RF_HYBRID_NOSCALE   Find a root of the derivative via the Hybrid method; no internal scaling
\return	an \ref apop_model with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.

By the way, the output model will have the expected score in the \c more slot. Do enough people use this to give it its own slot in the \c apop_model struct?
 \ingroup mle */
apop_model *	apop_maximum_likelihood(apop_data * data, apop_model dist){
  apop_mle_settings   *mp;
  infostruct    info    = { .data   = data,
                            .use_constraint = 1};
    info.gradientp      = malloc(sizeof(double*)); *info.gradientp = NULL;
    info.gradient_list  = malloc(sizeof(gsl_matrix*)); *info.gradient_list = NULL;
    info.score_list  = malloc(sizeof(gsl_vector*)); *info.score_list = NULL;
    info.gsize          = malloc(sizeof(size_t)); *info.gsize = 0;
    info.trace_file     = malloc(sizeof(FILE *)); *info.trace_file = NULL;
    if(strcmp(dist.method_name, "MLE")){
        mp          = apop_mle_settings_alloc(data, dist);
        info.model  = mp->model;
    } else {
        info.model  = apop_model_copy(dist);
        mp          = info.model->method_settings;
    }
    if (mp->trace_path)
        info.trace_path = mp->trace_path;
    apop_model_clear(data, info.model);
	if (mp->method == APOP_SIMAN)
        return apop_annealing(&info);  //below.
    else if (mp->method==APOP_SIMPLEX_NM)
		return apop_maximum_likelihood_no_d(data, &info);
    else if (mp->method == APOP_RF_NEWTON    ||
                mp->method == APOP_RF_BROYDEN  ||
                mp->method == APOP_RF_HYBRID  ||
                mp->method == APOP_RF_HYBRID_NOSCALE ) 
        return  find_roots (info);
        //return  find_roots (data, &dist);
	//else, Conjugate Gradient:
	return apop_maximum_likelihood_w_d(data, &info);
}

/** Input an earlier estimate, and then I will re-start the MLE search
 where the last one ended. You can specify greater precision or a new
 search method.

If the estimate converged to an OK value, then restart the converged value; else use the
starting point from the last estimate.

Only one estimate is returned, either the one you sent in or a new
one. The loser (which may be the one you sent in) is freed. That is,
there is no memory leak when you do
\code
est = apop_estimate_restart(est, alt_model);
\endcode

 \param e   An \ref apop_model that is the output from a prior MLE estimation.
 \param alt_model  Another not-yet-parametrized model that will be re-estimated with (1) the same data and (2) a <tt>starting_pt</tt> equal
 to the parameters of <tt>e</tt>. If this is <tt>NULL</tt>, then copy off <tt>e</tt> and restart from the end of the last estimation.

\return         At the end of this procedure, we'll have two \ref
    apop_model structs: the one you sent in, and the one produced using the
    new method/scale. If the new estimate includes any NaNs/Infs, then
    the old estimate is returned (even if the old estimate included
    NaNs/Infs). Otherwise, the estimate with the largest log likelihood
    is returned.

\ingroup mle
\todo The tolerance for testing boundaries are hard coded (1e4). Will need to either add another input term or a global var.
*/ 
apop_model * apop_estimate_restart (apop_model *e, apop_model *copy){
  if (!copy)
      copy = apop_model_copy(*e);
  apop_mle_settings* prm;
            //copy off the old params; modify the starting pt, method, and scale
    if (apop_vector_bounded(e->parameters->vector, 1e4)){
        double      *start_pt2 = apop_vector_to_array(e->parameters->vector);
        if (!copy->method_settings)
            prm = apop_mle_settings_alloc(copy->data,*copy);
        ((apop_mle_settings*)copy->method_settings)->starting_pt	= start_pt2;
    }
  /*
  apop_model *copy  = apop_model_copy(*e);
  apop_mle_settings *old_params   = e->method_settings;
  apop_mle_settings *new_params   = copy->method_settings; 

            //copy off the old params; modify the starting pt, method, and scale
    if (apop_vector_bounded(e->parameters->vector, 1e4)){
        apop_vector_to_array(e->parameters->vector, &start_pt2);
	    new_params->starting_pt	= start_pt2;
    }
    else
	    new_params->starting_pt	= old_params->starting_pt;
    new_params->tolerance   = old_params->tolerance * scale;
    new_params->step_size   = old_params->step_size * scale;
    new_params->method	    = new_method;
	copy                    = e->estimate(e->data, copy);
    */
            //Now check whether the new output is better than the old
//printf("orig: 1st: %g, ll %g\n", e->parameters->vector->data[0],e->llikelihood );
//printf("copy: 1st: %g, ll %g\n", copy->parameters->vector->data[0],copy->llikelihood );
  apop_model *newcopy = apop_estimate(e->data, *copy);
    apop_model_free(copy);
    if (apop_vector_bounded(newcopy->parameters->vector, 1e4) && newcopy->llikelihood > e->llikelihood){
        apop_model_free(e);
        return newcopy;
    } //else:
    apop_model_free(newcopy);
    return e;
}


//////////////////////////
// Simulated Annealing.

/** \page simanneal Notes on simulated annealing

Simulated annealing is a controlled random walk.
As with the other methods, the system tries a new point, and if it
is better, switches. Initially, the system is allowed to make large
jumps, and then with each iteration, the jumps get smaller, eventually
converging. Also, there is some decreasing probability that if the new
point is {\em less} likely, it will still be chosen. Simulated annealing
is best for situations where there may be multiple local optima. Early
in the random walk, the system can readily jump from one to another;
later it will fine-tune its way toward the optimum. The number of points
tested is basically not dependent on the function: if you give it a
4,000 step program, that is basically how many steps it will take.
If you know your function is globally convex (as are most standard
probability functions), then this method is overkill.

The GSL's simulated annealing system doesn't actually do very much. It
basically provides a for loop that calls a half-dozen functions that we
the users get to write. So, the file \ref apop_mle.c handles all of this
for you. The likelihood function is taken from the model, the metric
is the Manhattan metric, the copy/destroy functions are just the usual
vector-handling fns., et cetera. The reader who wants further control
is welcome to override these functions.

Verbosity: if ep->verbose==1, show likelihood,  temp, &c. in a table;
if ep->verbose>1, show that plus the vector of params.

 \ingroup mle
 */

static void an_record(infostruct *i, double energy){
    gsl_vector *grad = gsl_vector_alloc(i->beta->size);
    dnegshell(i->beta, i, grad);   //make sure the current param is in place. 
    record_gradient_info(energy, grad, i);
    gsl_vector_free(grad);
}

static double annealing_energy(void *in) {
  infostruct *i      = in;
  double energy      = negshell(i->beta, i);
  apop_mle_settings *p = i->model->method_settings;
    if (p->want_cov)
        an_record(i, energy);
    return energy;
}

static double annealing_distance(void *xin, void *yin) {
/** We use the Manhattan metric to correspond to the annealing_step fn below.  */
  gsl_vector    *from   = apop_vector_copy(((infostruct*)xin)->beta);
  gsl_vector    *to     = apop_vector_copy(((infostruct*)yin)->beta);
  gsl_vector_div(from, ((infostruct*)xin)->starting_pt);
  gsl_vector_div(to  , ((infostruct*)xin)->starting_pt);//starting pts are the same.
  return apop_vector_grid_distance(from, to);
  //return apop_vector_grid_distance(((infostruct*)xin)->beta, ((infostruct*)yin)->beta);
}

static void annealing_step(const gsl_rng * r, void *in, double step_size){
/** The algorithm: 
    --randomly pick dimension
    --shift by some amount of remaining step size
    --repeat for all dims
This will give a move \f$\leq\f$ step_size on the Manhattan metric.
*/
  infostruct  *i          = in;
  int         sign, j;
  double      amt, scale;
  int           vsize       =(i->model->parameters->vector ? i->model->parameters->vector->size :0),
                msize1      =(i->model->parameters->matrix ? i->model->parameters->matrix->size1 :0),
                msize2      =(i->model->parameters->matrix ? i->model->parameters->matrix->size2:0);
    double cutpoints[i->beta->size+1];
    cutpoints[0]                = 0;
    cutpoints[i->beta->size]    = 1;
    for (j=1; j< i->beta->size; j++)
        cutpoints[j] = gsl_rng_uniform(r);

    for (j=0; j< i->beta->size; j++){
        sign    = (gsl_rng_uniform(r) > 0.5) ? 1 : -1;
        scale   = gsl_vector_get(i->starting_pt, j);
        scale   = scale ? scale : 1; //if starting pt is zero, assume 1.
        amt     = cutpoints[j+1]- cutpoints[j];
        apop_vector_increment(i->beta, j,  amt * sign * scale * step_size); 
    }
    apop_data_free(i->model->parameters);
    i->model->parameters      = apop_data_unpack(i->beta, vsize, msize1, msize2);
    if (i->model->constraint && i->model->constraint(i->data, i->model)){
        gsl_vector *cv  = apop_data_pack(i->model->parameters);
        gsl_vector_memcpy(i->beta, cv);
        gsl_vector_free(cv);
    }
}

static void annealing_print(void *xp) {
    apop_vector_show(((infostruct*)xp)->beta);
}

static void annealing_print2(void *xp) { return; }

static void annealing_memcpy(void *xp, void *yp){
  infostruct    *yi = yp;
  infostruct    *xi = xp;
    if (yi->beta && yi->beta !=xi->beta) 
        gsl_vector_free(yi->beta);
    memcpy(yp, xp, sizeof(infostruct));
    yi->beta = gsl_vector_alloc(((infostruct*)xp)->beta->size);
    gsl_vector_memcpy(yi->beta, ((infostruct*)xp)->beta);
}

static void *annealing_copy(void *xp){
    infostruct *out = malloc(sizeof(infostruct));
    memcpy(out, xp, sizeof(infostruct));
    out->beta        = gsl_vector_alloc(((infostruct*)xp)->beta->size);
    gsl_vector_memcpy(out->beta, ((infostruct*)xp)->beta);
    return out;
}

static void annealing_free(void *xp){
    gsl_vector_free(((infostruct*)xp)->beta);
    free(xp);
}

//the starting point is really the list of scaling factors. They can't be zero.
static double set_start(double in){
    return in ? in : 1; }

static gsl_siman_params_t set_params(apop_model *ep, apop_mle_settings **mp){
gsl_siman_params_t simparams;
    if (ep && !strcmp(ep->method_name, "MLE")){
        *mp  = ep->method_settings;
        simparams.n_tries       = (*mp)->n_tries;
        simparams.iters_fixed_T = (*mp)->iters_fixed_T;
        simparams.step_size     = (*mp)->step_size;
        simparams.k             = (*mp)->k;
        simparams.t_initial     = (*mp)->t_initial;
        simparams.mu_t          = (*mp)->mu_t;
        simparams.t_min         = (*mp)->t_min;
    } else{
        simparams.n_tries =200; //The number of points to try for each step. 
        simparams.iters_fixed_T = 200;  //The number of iterations at each temperature. 
        simparams.step_size  = 1.;  //The maximum step size in the random walk. 
        simparams.k          = 1.0, //cooling schedule data
        simparams.t_initial   = 50,   
        simparams.mu_t        = 1.002, 
        simparams.t_min       = 5.0e-1;
    }
    return simparams;
}

jmp_buf anneal_jump;

static void anneal_sigint(){ longjmp(anneal_jump,1); }

apop_model * apop_annealing(infostruct *i){
  apop_model            *ep = i->model;
  apop_mle_settings       *mp = NULL;
  gsl_vector    *beta;
  int           vsize       =(i->model->parameters->vector ? i->model->parameters->vector->size :0),
                msize1      =(i->model->parameters->matrix ? i->model->parameters->matrix->size1 :0),
                msize2      =(i->model->parameters->matrix ? i->model->parameters->matrix->size2:0);
  int           paramct = vsize + msize1*msize2;
  gsl_siman_params_t    simparams = set_params(ep, &mp);
    static const gsl_rng   * r    = NULL;
    if (mp && mp->rng) 
        r    =  mp->rng;
    if (!r)
        r =  apop_rng_alloc(8) ; 
    if (mp && mp->starting_pt)
        beta = apop_array_to_vector(mp->starting_pt, paramct);
    else{
        beta  = gsl_vector_alloc(paramct);
        gsl_vector_set_all(beta, 1);
    }
    i->starting_pt      = apop_vector_map(beta, set_start);
	i->beta             = beta;
    i->use_constraint   = 0; //negshell doesn't check it; annealing_step does.
    gsl_siman_print_t printing_fn   = NULL;
    if (mp && mp->verbose>1)
        printing_fn = annealing_print;
    else if (mp && mp->verbose)
        printing_fn = annealing_print2;
    if (!setjmp(anneal_jump)){
        signal(SIGINT, anneal_sigint);
        gsl_siman_solve(r,        //   const gsl_rng * r
          i,                //   void * x0_p
          annealing_energy, //   gsl_siman_Efunc_t Ef
          annealing_step,   //   gsl_siman_step_t take_step
          annealing_distance, // gsl_siman_metric_t distance
          printing_fn,      //gsl_siman_print_t print_position
          annealing_memcpy, //   gsl_siman_copy_t copyfunc
          annealing_copy,   //   gsl_siman_copy_construct_t copy_constructor
          annealing_free,   //   gsl_siman_destroy_t destructor
          paramct,           //   size_t element_size
          simparams);        //   gsl_siman_params_t params
    }
    signal(SIGINT, NULL);
    i->model->parameters   = apop_data_unpack(i->beta, vsize, msize1, msize2); 
    produce_covariance_matrix(i->model, i);
    apop_estimate_parameter_t_tests(i->model);
    return i->model;
}
