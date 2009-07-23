/** \file apop_mle.c	The MLE functions. Call them with an \ref apop_model.

This file includes a number of distributions and models whose parameters
one would estimate using maximum likelihood techniques.

Each typically includes the likelihood function, the derivative of the
likelihood function, which is used by \c apop_maximum_likelihood. The \c apop_estimate function defaults to using maximum likelihood
estimation if there is no model-specific estimation routine provided.

At the bottom are the maximum likelihood procedures themselves. There
are three: the no-derivative version, the with-derivative version,
and the simulated annealing routine.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
*/
#include "model/model.h"
#include "asst.h" //apop_rng_alloc
#include "output.h"
#include "likelihoods.h"
#include <assert.h>
#include <setjmp.h>
#include <signal.h>
#include <gsl/gsl_deriv.h>

//in apop_regress.c:
void apop_estimate_parameter_t_tests (apop_model *est);

/** \page trace_path Plotting the path of an ML estimation.

If \c trace_path (in the \ref apop_mle_settings struct) has a name of positive
length, then every time the MLE evaluates the function, then the value
will be output to a table in the database/a file with the given name
(depending on the value of \ref apop_opts_type "apop_opts.output_type"). You can then plot this
table to get an idea of the path the estimation routine used to arrive
at its MLE.

To write to a pipe or stdout, set \ref apop_opts_type "apop_opts.output_type" appropriately and set \c trace_path to the literal string \c "NULL".


Below is a sample of the sort of output one would get:<br>
\image latex "search.gif" "An ML search, tracing out the surface of the function" width=\textwidth
\image html "search.gif" "An ML search, tracing out the surface of the function" 

\ingroup mle
*/

/** Allocate settings for a maximum likelihood estimation.

 \param parent  A pointer to an allocated \ref apop_model struct.
 */
apop_mle_settings *apop_mle_settings_alloc(apop_model *parent){
  apop_mle_settings *setme =   calloc(1,sizeof(apop_mle_settings));
    setme->starting_pt      = NULL;
    setme->tolerance        = 1e-2;
    setme->method           = parent->score ? APOP_CG_PR : APOP_SIMPLEX_NM;
    setme->verbose          = 0;
    setme->use_score        = 1;
    setme->step_size        = 0.05;
    setme->delta            = 1e-2;
    setme->want_cov         = 1;
//siman:
    //siman also uses step_size  = 1.;  
    setme->n_tries          = 200;  //The number of points to try for each step. 
    setme->iters_fixed_T    = 200;   //The number of iterations at each temperature. 
    setme->k                = 1.0;  //The maximum step size in the random walk. 
    setme->t_initial        = 50;   //cooling schedule data
    setme->mu_t             = 1.002; 
    setme->t_min            = 5.0e-1;
    setme->rng              = NULL;
    return setme;
}

void *apop_mle_settings_copy(apop_mle_settings * in){
  apop_mle_settings *setme = malloc(sizeof(apop_mle_settings));
    memmove(setme, in, sizeof(apop_mle_settings));
    return setme;
}

void apop_mle_settings_free(void *in){ apop_mle_settings * fin = in; free(fin); }

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
    char        *trace_path;
    FILE        **trace_file;
}   infostruct;

static apop_model * apop_annealing(infostruct*);                         //below.
static int dnegshell (const gsl_vector *, void * , gsl_vector * g); //below.

static double one_d(double b, void *in){
  infostruct    *i   =in;
    gsl_vector_set(i->gp->beta, i->gp->dimension, b);
    apop_data_unpack(i->gp->beta, i->model->parameters);
	double out= (*(i->f))(i->gp->d, i->model);
    return out;
}


#include "apop_findzeros.c"

/* For each element of the parameter set, jiggle it to find its
 gradient. Return a vector as long as the parameter list. */
static void apop_internal_numerical_gradient(apop_fn_with_params ll, infostruct* info, gsl_vector *out){
  int		    j;
  gsl_function	F;
  double		result, err;
  grad_params 	gp;
  infostruct    i;
  gsl_vector    *beta   = apop_data_pack(info->model->parameters);
  apop_mle_settings   *mp = apop_settings_get_group(info->model, "apop_mle");
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
    gsl_vector_free(beta);
}

/**The GSL provides one-dimensional numerical differentiation; here's the multidimensional extension.

 If \c m has \c method_settings of type \c apop_ml_params, then the \c delta element is used for the differential.
 
 \code
 gsl_vector *gradient = apop_numerical_gradient(data, your_parametrized_model);
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
  apop_mle_settings *mp = apop_settings_get_group(m, "apop_mle");
    if(!mp){
        Apop_settings_add_group(m, apop_mle, m);
        mp = apop_settings_get_group(m, "apop_mle");
    }
    i.model = m;
    i.data  = data;
    apop_internal_numerical_gradient(ll, &i, out);
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

static double negshell (const gsl_vector *beta, void * in){
  infostruct    *i              = in;
  double		penalty         = 0,
                out             = 0; 
  double 	(*f)(apop_data *, apop_model *);
    f   = i->model->log_likelihood? i->model->log_likelihood : i->model->p;
    if (!f)
        apop_error(0, 's', "The model you sent to the MLE function has neither log_likelihood element nor p element.\n");
    apop_data_unpack(beta, i->model->parameters);
	if (i->use_constraint && i->model->constraint)
		penalty	= i->model->constraint(i->data, i->model);
    if(penalty){
        gsl_vector *o = apop_data_pack(i->model->parameters);
        gsl_vector_memcpy((gsl_vector *)beta, o);
        gsl_vector_free(o);
    }
    out = penalty - f(i->data, i->model); //negative llikelihood
    if (penalty)
        apop_data_unpack(i->beta, i->model->parameters);
    if (i->trace_path && strlen(i->trace_path))
        tracepath(i->model->parameters->vector,-out, i->trace_path, i->trace_file);
    //The next line is not used anywhere else here. It's just for output.
    i->model->llikelihood = i->model->log_likelihood? -out : log(-out); //negative negative llikelihood.
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
  apop_mle_settings *mp       =  apop_settings_get_group(i->model, "apop_mle");
    apop_data_unpack(beta, i->model->parameters);
    if(i->model->constraint)
        if(i->model->constraint(i->data, i->model)){
            gsl_vector *o = apop_data_pack(i->model->parameters);
            gsl_vector_memcpy((gsl_vector *)beta, o);
            gsl_vector_free(o);
        }
    if (mp->use_score && i->model->score)
        i->model->score(i->data, g, i->model);
    else {
        apop_fn_with_params ll  = i->model->log_likelihood ? i->model->log_likelihood : i->model->p;
        apop_internal_numerical_gradient(ll, i, g);
    }
    if (i->trace_path && strlen(i->trace_path))
        negshell (beta,  in);
    gsl_vector_scale(g, -1);
    return GSL_SUCCESS;
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

/* * Calculate the Hessian.

  This is a synonym for \ref apop_numerical_second_derivative, q.v.
gsl_matrix * apop_numerical_hessian(apop_model dist, gsl_vector *beta, apop_data * d, void *params){
	return apop_numerical_second_derivative(dist, beta, d, params);
}
*/

/* * Feeling lazy? Rather than doing actual pencil-and-paper math to find
your variance-covariance matrix, just use the negative inverse of the Hessian.

\param est	The model with the parameters already calculated. The var/covar matrix will be placed in est->covariance.
\param data	The data
\ingroup basic_stats
*/
void apop_numerical_covariance_matrix(apop_model *est, apop_data *data){
    //As you can see, this is a placeholder.
    return;
}

typedef struct {
    apop_model *base_model;
    int *current_index;
} apop_model_for_infomatrix_struct;

/*This model is needed to make life easier with the
 apop_numerical_gradient function.  */ 

static double apop_fn_for_infomatrix(apop_data *d, apop_model *m){
    static gsl_vector *v = NULL;
    apop_model_for_infomatrix_struct *settings = m->more;
    apop_model *mm = settings->base_model;
    if (mm->score){
        if (!v || v->size != mm->parameters->vector->size){
            if (v) gsl_vector_free(v);
            v = gsl_vector_alloc(mm->parameters->vector->size);
        }
         mm->score(d, v, mm);
        return gsl_vector_get(v, *settings->current_index);
    } //else:
        gsl_vector *vv = apop_numerical_gradient(d, mm);
        double out = gsl_vector_get(vv, *settings->current_index);
        gsl_vector_free(vv);
        return out;
}


apop_model apop_model_for_infomatrix = {"Ad hoc model for working out the information matrix.", .log_likelihood = apop_fn_for_infomatrix};

/* We produce the covariance matrix via the derivative of the score
 function at the parameter. I follow Efron and Hinkley in using the
 estimated information matrix---the value of the information matrix at
 the estimated value of the score---not the expected information matrix
 that is the integral over all possible data. See Pawitan 2001 (who
 cribbed off of Efron and Hinkley) or Klemens 2008 (who cribbed off of
 both) for further details. If you really want to use E(I), check
 out a version of this file from before the end of 2007.*/
static void produce_covariance_matrix(apop_model * est, infostruct *i){
  int    j, k;
  size_t betasize  = i->model->parameters->vector->size;
  gsl_matrix *preinv    = gsl_matrix_calloc(betasize, betasize);
  gsl_vector *dscore    = gsl_vector_alloc(betasize);
    apop_model_for_infomatrix_struct m;
    m.base_model        = i->model;
    m.current_index     = &k;
    apop_model_for_infomatrix.more       =  &m;
    apop_model_for_infomatrix.parameters = i->model->parameters;
    if (apop_settings_get_group(i->model, "apop_mle"))
        apop_settings_copy_group(&apop_model_for_infomatrix, i->model, "apop_mle");
    for (k=0; k< betasize; k++){
        dscore = apop_numerical_gradient(i->data, &apop_model_for_infomatrix);
        //We get two estimates of the (k,j)th element, which are often very close,
        //and take the mean.
        for (j=0; j< betasize; j++){
            apop_matrix_increment(preinv, k, j, gsl_vector_get(dscore, j)/2);
            apop_matrix_increment(preinv, j, k, gsl_vector_get(dscore, j)/2);
        }
        gsl_vector_free(dscore);
    }
    if (apop_opts.verbose > 1){
        printf("The estimated Hessian:\n");
        apop_matrix_show(preinv);
    }
    gsl_matrix *inv = apop_matrix_inverse(preinv);
    gsl_matrix_scale(inv, -1);
    est->covariance = apop_matrix_to_data(inv);
    if (est->parameters->names->row){
        apop_name_stack(est->covariance->names, est->parameters->names, 'r');
        apop_name_cross_stack(est->covariance->names, est->parameters->names, 'c', 'r');
    }
    gsl_matrix_free(preinv);
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
  apop_mle_settings           *mp         = apop_settings_get_group(est, "apop_mle");
  assert(mp);
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
    apop_data_unpack(s->x, est->parameters);
	gsl_multimin_fdfminimizer_free(s);
	if (mp->starting_pt==NULL) 
		gsl_vector_free(x);
    if (est->parameters->vector && !est->parameters->matrix){
        produce_covariance_matrix(est, i);
        apop_estimate_parameter_t_tests (est);
    }
	return est;
}

static apop_model *	apop_maximum_likelihood_no_d(apop_data * data, infostruct * i){
  apop_model		        *est        = i->model;
  apop_mle_settings           *mp         = apop_settings_get_group(est, "apop_mle");
  assert(mp);
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
    ctrl_c      =
	est->status	= 0;	//assume failure until we score a success.
	if (mp->starting_pt==NULL){
	    x	= gsl_vector_alloc(betasize);
  		gsl_vector_set_all (x,  1);
    } else
		x   = apop_array_to_vector(mp->starting_pt, betasize);
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
    apop_data_unpack(s->x, est->parameters);
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
\li    APOP_RF_HYBRID   Find a root of the derivative via the Hybrid method
\li    APOP_RF_HYBRID_NOSCALE   Find a root of the derivative via the Hybrid method; no internal scaling
\return	an \ref apop_model with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.

By the way, the output model will have the expected score in the \c more slot. Do enough people use this to give it its own slot in the \c apop_model struct?
 \ingroup mle */
apop_model *	apop_maximum_likelihood(apop_data * data, apop_model dist){
  apop_mle_settings   *mp = apop_settings_get_group(&dist, "apop_mle");
  infostruct    info    = { .data   = data,
                            .use_constraint = 1};
    info.trace_file     = malloc(sizeof(FILE *)); *info.trace_file = NULL;
    info.model          = apop_model_copy(dist);
    info.model->data    = data;
    if(!mp){
        Apop_settings_add_group(info.model, apop_mle, info.model);
        mp          =  apop_settings_get_group(info.model, "apop_mle");
    } 
    apop_model_prep(data, info.model);
    if (mp->trace_path)
        info.trace_path = mp->trace_path;
	if (mp->method == APOP_SIMAN)
        return apop_annealing(&info);  //below.
    else if (mp->method==APOP_SIMPLEX_NM)
		return apop_maximum_likelihood_no_d(data, &info);
    else if (mp->method == APOP_RF_NEWTON    ||
                mp->method == APOP_RF_HYBRID  ||
                mp->method == APOP_RF_HYBRID_NOSCALE ) 
        return  find_roots (info);
	//else, Conjugate Gradient:
	return apop_maximum_likelihood_w_d(data, &info);
}

/** 
  The simplest use of this function is to restart a model at the current parameter estimates.

  \code
apop_model *m = apop_estimate(data, model_using_an_MLE_search);
for (int i=0; i< 10; i++)
    m = apop_estimate_restart(m);
apop_data_show(m);
  \endcode

By adding a line to reduce the tolerance each round [e.g., <tt>Apop_settings_add(m, apop_mle, tolerance, pow(10,-i))</tt>], you can hone in on a precise optimum.
 
You may have a new estimation method, such as first doing a simulated annealing search, then honing in via a conjugate gradient method:
  \code
apop_model *m = apop_estimate(data, model_with_siman_specified);
m = apop_estimate_restart(m, model_with_cg_specified);
apop_data_show(m);
  \endcode

Only one estimate is returned, either the one you sent in or a new
one. The loser (which may be the one you sent in) is freed. That is,
there is no memory leak in the above loop.

 \param e   An \ref apop_model that is the output from a prior MLE estimation. (No default, must not be \c NULL.)
 \param copy  Another not-yet-parametrized model that will be re-estimated with (1) the same data and (2) a <tt>starting_pt</tt> as per the next setting (probably
 to the parameters of <tt>e</tt>). If this is <tt>NULL</tt>, then copy off <tt>e</tt>. (Default = \c NULL)
 \param starting_pt "ep"=last estimate of the first model (i.e., its current parameter estimates); "es"= starting point originally used by the first model; "np"=current parameters of the new (second) model; "ns"=starting point specified by the new model's MLE settings. (default = "ep")
 \param boundary I test whether the starting point you give me is outside this certain bound, so I can warn you if there's divergence in your sequence of re-estimations. (default: 1e8)

\return         At the end of this procedure, we'll have two \ref
    apop_model structs: the one you sent in, and the one produced using the
    new method/scale. If the new estimate includes any NaNs/Infs, then
    the old estimate is returned (even if the old estimate included
    NaNs/Infs). Otherwise, the estimate with the largest log likelihood
    is returned.

\ingroup mle

This function uses the \ref designated syntax for inputs.
*/ 
APOP_VAR_HEAD apop_model * apop_estimate_restart (apop_model *e, apop_model *copy, char * starting_pt, double boundary){
    apop_model * apop_varad_var(e, NULL);
    apop_assert(e, 0, 0, 's', "You gave me a NULL model to restart.");
    apop_model * apop_varad_var(copy, NULL);
    char * apop_varad_var(starting_pt, NULL);
    double apop_varad_var(boundary, 1e8);
    return apop_estimate_restart_base(e, copy, starting_pt, boundary);
APOP_VAR_ENDHEAD
    gsl_vector *v;
  if (!copy)
      copy = apop_model_copy(*e);
  apop_mle_settings* prm0 = apop_settings_get_group(e, "apop_mle");
  apop_mle_settings* prm = apop_settings_get_group(copy, "apop_mle");
            //copy off the old params; modify the starting pt, method, and scale
    if (starting_pt && !strcmp(starting_pt, "es"))
        v = apop_array_to_vector(prm0->starting_pt);
    else if (starting_pt && !strcmp(starting_pt, "ns")){
        int size =sizeof(prm->starting_pt)/sizeof(double);
        v = apop_array_to_vector(prm->starting_pt, size);
        prm0->starting_pt	= malloc(sizeof(double)*size);
        memcpy(prm0->starting_pt, prm->starting_pt, sizeof(double)*size);
    }
    else if (starting_pt && !strcmp(starting_pt, "np")){
        v                   = apop_data_pack(copy->parameters); 
        prm->starting_pt	= malloc(sizeof(double)*v->size);
        memcpy(prm->starting_pt, v->data, sizeof(double)*v->size);
    }
    else {//"ep" or default.
        v                   = apop_data_pack(e->parameters); 
        prm->starting_pt	= malloc(sizeof(double)*v->size);
        memcpy(prm->starting_pt, v->data, sizeof(double)*v->size);
    }
    apop_assert(apop_vector_bounded(v, boundary), e, 0, 'c', "Your model has diverged (element(s) > %g); returning your original model without restarting.", boundary);
    gsl_vector_free(v);
        
    apop_model *newcopy = apop_estimate(e->data, *copy);
    apop_model_free(copy);
    //Now check whether the new output is better than the old
    if (apop_vector_bounded(newcopy->parameters->vector, boundary) && newcopy->llikelihood > e->llikelihood){
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

static double annealing_energy(void *in) {
  infostruct *i      = in;
  return negshell(i->beta, i);
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
  double cutpoints[i->beta->size+1];
    cutpoints[0]                = 0;
    cutpoints[i->beta->size]    = 1;
    for (j=1; j< i->beta->size; j++)
        cutpoints[j] = gsl_rng_uniform(r);

    for (j=0; j< i->beta->size; j++){
        sign    = (gsl_rng_uniform(r) > 0.5) ? 1 : -1;
        scale   = gsl_vector_get(i->starting_pt, j);
        amt     = cutpoints[j+1]- cutpoints[j];
        apop_vector_increment(i->beta, j,  amt * sign * scale * step_size); 
    }
    apop_data_unpack(i->beta, i->model->parameters);
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

static gsl_siman_params_t set_params(apop_mle_settings *mp){
gsl_siman_params_t simparams;
    simparams.n_tries       = mp->n_tries;
    simparams.iters_fixed_T = mp->iters_fixed_T;
    simparams.step_size     = mp->step_size;
    simparams.k             = mp->k;
    simparams.t_initial     = mp->t_initial;
    simparams.mu_t          = mp->mu_t;
    simparams.t_min         = mp->t_min;
    return simparams;
}

jmp_buf anneal_jump;
static void anneal_sigint(){ longjmp(anneal_jump,1); }

apop_model * apop_annealing(infostruct *i){
  apop_model            *ep = i->model;
  apop_mle_settings     *mp = apop_settings_get_group(ep, "apop_mle");
  assert(mp);
  gsl_siman_params_t    simparams = set_params(mp);
  gsl_vector    *beta;
  int           vsize       =(i->model->parameters->vector ? i->model->parameters->vector->size :0),
                msize1      =(i->model->parameters->matrix ? i->model->parameters->matrix->size1 :0),
                msize2      =(i->model->parameters->matrix ? i->model->parameters->matrix->size2:0);
  int           paramct = vsize + msize1*msize2;
  static const gsl_rng   * r    = NULL;
    if (!r)
        r = mp->rng ? mp->rng : apop_rng_alloc(8);
    if (mp->starting_pt)
        beta = apop_array_to_vector(mp->starting_pt, paramct);
    else{
        beta  = gsl_vector_alloc(paramct);
        gsl_vector_set_all(beta, 1);
    }
    i->starting_pt                  = apop_vector_map(beta, set_start);
	i->beta                         = beta;
    i->use_constraint               = 0; //negshell doesn't check it; annealing_step does.
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
    apop_data_unpack(i->beta, i->model->parameters); 
    produce_covariance_matrix(i->model, i);
    apop_estimate_parameter_t_tests(i->model);
    if (mp->rng)
        r = NULL;
    return i->model;
}
