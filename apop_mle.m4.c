/** \file apop_mle.c	The MLE functions. Call them with an \ref apop_model.

This file includes a number of distributions and models whose parameters one would estimate using maximum likelihood techniques.

It has (more-or-less) a single public function: \ref apop_maximum_likelihood, and you don't even need to use that one, because the \c apop_estimate function defaults to using it if there is no model-specific estimation routine provided.

At the bottom are the maximum likelihood procedures themselves. There are four: Newton-type derivative methods, the no-derivative version, the with-derivative version, and the simulated annealing routine.*/

/*Copyright (c) 2006--2010 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#include "apop_internal.h"
#include <setjmp.h>
#include <signal.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>

typedef double 	(*apop_fn_with_params) (apop_data *, apop_model *);
typedef	void 	(*apop_df_with_void)(const gsl_vector *beta, void *d, gsl_vector *gradient);
typedef	void 	(*apop_fdf_with_void)(const gsl_vector *beta, void *d, double *f, gsl_vector *df);

typedef struct {
	gsl_vector	*beta;
	int		    dimension;
} grad_params;

typedef struct {
    apop_model  *model;
    apop_data   *data;
    apop_fn_with_params   *f;
    grad_params *gp; //Used only by apop_internal_numerical_gradient.
    gsl_vector  *beta, *starting_pt;
    int         use_constraint;
    char        *trace_path;
    FILE        **trace_file;
    double      best_ll;
    char        want_cov, want_predicted, want_tests, want_info;
    jmp_buf     bad_eval_jump;
}   infostruct;

static apop_model * find_roots (infostruct p); //see end of file.

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

double default_delta = 1e-3;

/* Generate support fns (esp. initializers) for apop_mle_settings and apop_parts_wanted structs. */
Apop_settings_copy(apop_parts_wanted, )
Apop_settings_free(apop_parts_wanted, )
Apop_settings_init(apop_parts_wanted,
    Apop_varad_set(covariance, 'n');
    Apop_varad_set(predicted, 'n');
    Apop_varad_set(tests, 'n');
    Apop_varad_set(info, 'n');
)

Apop_settings_copy(apop_mle, )
Apop_settings_free(apop_mle, )
Apop_settings_init(apop_mle,
    Apop_varad_set(starting_pt, NULL);
    Apop_varad_set(tolerance, 1e-5);
    Apop_varad_set(max_iterations, 5000);
    Apop_varad_set(method, APOP_UNKNOWN_ML);//default picked in apop_maximum_likelihood
    Apop_varad_set(verbose, 0);
    Apop_varad_set(use_score, 'y'); //deprecated.
    Apop_varad_set(step_size, 0.05);
    Apop_varad_set(delta, default_delta);
    Apop_varad_set(want_cov, 'y');
    if (in.want_cov == 1) out->want_cov = 'y';
    Apop_varad_set(dim_cycle_tolerance, 0);
//siman:
    //siman also uses step_size  = 1.;  
    Apop_varad_set(n_tries, 5);  //The number of points to try for each step. 
    Apop_varad_set(iters_fixed_T, 5);   //The number of iterations at each temperature. 
    Apop_varad_set(k, 1.0);  //The maximum step size in the random walk. 
    Apop_varad_set(t_initial, 50);   //cooling schedule data
    Apop_varad_set(mu_t, 1.002); 
    Apop_varad_set(t_min, 5.0e-1);
    Apop_varad_set(rng, NULL);
)

//deprecated; left to make some examples in Modeling with Data coherent.
apop_mle_settings *apop_mle_settings_alloc(apop_model *parent){
    return apop_mle_settings_init((apop_mle_settings){ }); }

//      MLE support functions
//Including numerical differentiation and a couple of functions to
//negate the likelihood fns without bothering the user.

static apop_model * apop_annealing(infostruct*); //below.

static double one_d(double b, void *in){
    infostruct *i  = in;
    double penalty = 0;
    gsl_vector_set(i->gp->beta, i->gp->dimension, b);
    apop_data_unpack(i->gp->beta, i->model->parameters);
	if (i->model->constraint)
		penalty	= i->model->constraint(i->data, i->model);
	return (*(i->f))(i->data, i->model) + penalty;
}

//Numeric first and second derivatives.

/* For each element of the parameter set, jiggle it to find its
 gradient. Return a vector as long as the parameter list. */
static void apop_internal_numerical_gradient(apop_fn_with_params ll, 
                            infostruct* info, gsl_vector *out, double delta){
    double result, err;
    gsl_vector *beta = apop_data_pack(info->model->parameters, NULL, .all_pages='y');
    infostruct i = *info;
    i.f = &ll;
    i.gp = &(grad_params){ .beta = gsl_vector_alloc(beta->size)};
    gsl_function F = { .function= one_d, 
                       .params	= &i };
	for (size_t j=0; j< beta->size; j++){
		i.gp->dimension = j;
		gsl_vector_memcpy(i.gp->beta, beta);
		gsl_deriv_central(&F, gsl_vector_get(beta,j), delta, &result, &err);
		gsl_vector_set(out, j, result);
	}
    gsl_vector_free(beta);
}

/**The GSL provides one-dimensional numerical differentiation; here's the multidimensional extension.

\param data The data set to use for all evaluations. It remains constant throughout.
\param model The model, expressing the function whose derivative is sought. The gradient is taken via small changes along the model parameters.
\param delta The size of the differential. If you explicitly give me a \c delta, I'll use it. If \c delta is not specified,
 but \c model has \c method_settings of type \c apop_ml_params, then the \c delta element is used for the differential. Else, I use 1e-3.
 
 \code
 gsl_vector *gradient = apop_numerical_gradient(data, your_parametrized_model);
 \endcode

\li This function uses the \ref designated syntax for inputs.
\ingroup linear_algebra
 */
APOP_VAR_HEAD gsl_vector * apop_numerical_gradient(apop_data *data, apop_model *model, double delta){
    apop_data * apop_varad_var(data, NULL);
    apop_model * apop_varad_var(model, NULL);
    Nullcheck(model, NULL) Nullcheck_p(model, NULL)
    double apop_varad_var(delta, 0);
    if (!delta){
        apop_mle_settings *mp = apop_settings_get_group(model, apop_mle);
        delta = mp ? mp->delta : default_delta;
    }
APOP_VAR_ENDHEAD
  Get_vmsizes(model->parameters); //tsize
  apop_fn_with_params ll  = model->log_likelihood ? model->log_likelihood : model->p;
  Apop_stopif(!ll, return 0, 0, "Input model has neither p nor log_likelihood method. Returning zero.");
  gsl_vector        *out= gsl_vector_alloc(tsize);
  infostruct    i = (infostruct) {.model = model, .data = data};
    apop_internal_numerical_gradient(ll, &i, out, delta);
    return out;
}

typedef struct {
    apop_model *base_model;
    int *current_index;
} apop_model_for_infomatrix_struct;

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

apop_model apop_model_for_infomatrix = {"Ad hoc model for working out the information matrix.", 
                                                .log_likelihood = apop_fn_for_infomatrix};

/** Numerically estimate the matrix of second derivatives of the
parameter values. The math is
 simply a series of re-evaluations at small differential steps. [Therefore, it may be expensive to do this for a very computationally-intensive model.]  

\param data The data at which the model was estimated
\param model The model, with parameters already estimated
\param delta the step size for the differentials. The current default is around 1e-3.
\return The matrix of estimated second derivatives at the given data and parameter values.
 
\li This function uses the \ref designated syntax for inputs.
 */
APOP_VAR_HEAD apop_data * apop_model_hessian(apop_data * data, apop_model *model, double delta){
    apop_data * apop_varad_var(data, NULL);
    apop_model * apop_varad_var(model, NULL);
    Nullcheck(model, NULL)
    double apop_varad_var(delta, 0);
    if (!delta){
        apop_mle_settings *mp = apop_settings_get_group(model, apop_mle);
        delta = mp ? mp->delta : default_delta;
    }
APOP_VAR_ENDHEAD
    int    k;
    Get_vmsizes(model->parameters) //tsize
    size_t betasize  = tsize;
    apop_data *out    = apop_data_calloc(0, betasize, betasize);
    gsl_vector *dscore = gsl_vector_alloc(betasize);
    apop_model_for_infomatrix_struct ms = { .base_model = model, .current_index = &k, };
    apop_model *m = apop_model_copy(apop_model_for_infomatrix);
    m->parameters = model->parameters;
    m->more = &ms;
    if (apop_settings_get_group(model, apop_mle))
        apop_settings_copy_group(m, model, "apop_mle");
    for (k=0; k< betasize; k++){
        dscore = apop_numerical_gradient(data, m, delta);
        //We get two estimates of the (k,j)th element, which are often very close,
        //and take the mean.
        for (size_t j=0; j< betasize; j++){
            apop_matrix_increment(out->matrix, k, j, gsl_vector_get(dscore, j)/2);
            apop_matrix_increment(out->matrix, j, k, gsl_vector_get(dscore, j)/2);
        }
        gsl_vector_free(dscore);
    }
    if (model->parameters->names->row){
        apop_name_stack(out->names, model->parameters->names, 'r');
        apop_name_stack(out->names, model->parameters->names, 'c', 'r');
    }
    return out;
}

/** Produce the covariance matrix for the parameters of an estimated model via the derivative of the score function at the parameter. I.e., I find the second derivative via \ref apop_model_hessian , and take the negation of the inverse.

I follow Efron and Hinkley in using the estimated information matrix---the value of the information matrix at the estimated value of the score---not the expected information matrix that is the integral over all possible data. See Pawitan 2001 (who cribbed a little off of Efron and Hinkley) or Klemens 2008 (who directly cribbed off of both) for further details. 

 \param data The data by which your model was estimated
 \param model A model whose parameters have been estimated.
 \param delta The differential by which to step for sampling changes.  (default currently = 1e-3)
 \return A covariance matrix for the data. Also, if the data does not have a
 <tt>"Covariance"</tt> page, I'll set it to the result as well [i.e., I won't overwrite an
 existing covar].  

This function uses the \ref designated syntax for inputs.
 */
APOP_VAR_HEAD apop_data * apop_model_numerical_covariance(apop_data * data, apop_model *model, double delta){
    apop_data * apop_varad_var(data, NULL);
    apop_model * apop_varad_var(model, NULL);
    Nullcheck(model, NULL)
    double apop_varad_var(delta, 0);
    if (!delta){
        apop_mle_settings *mp = apop_settings_get_group(model, apop_mle);
        delta = mp ? mp->delta : default_delta;
    }
APOP_VAR_ENDHEAD
    apop_data *hessian = apop_model_hessian(data, model, delta);
    if (apop_opts.verbose > 1){
        printf("The estimated Hessian:\n");
        apop_data_show(hessian);
    }
    apop_data *out = apop_matrix_to_data(apop_matrix_inverse(hessian->matrix));
    gsl_matrix_scale(out->matrix, -1);
    if (hessian->names->row){
        apop_name_stack(out->names, hessian->names, 'r');
        apop_name_stack(out->names, hessian->names, 'c');
    }
    apop_data_free(hessian);
    if (!apop_data_get_page(model->parameters, "<Covariance>"))
        apop_data_add_page(model->parameters, out, "<Covariance>");
    return out;
}

///On to the interfaces between the models and the methods

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
                Apop_stopif(!(*tf = fopen(tp, "a")), /*keep going.*/,
                    0, "couldn't open %s for writing. Continuing without the path trace.", tp);
        }
        if (beta->size == 1)
            fprintf(*tf, "%g\t %g\n",  gsl_vector_get(beta,0), out);
        else
            fprintf(*tf, "%g\t %g\t %g\n",  gsl_vector_get(beta,0), gsl_vector_get(beta,1), out);
    }
}

/* Every actual evaluation of the function go through the negshell and dnegshell fns,
   because there are several things that have to be done beyond just getting
   model.log_likelihood:

--Negate, because statisticians and social scientists like to maximize; physicists like to minimize.
--Work out if the model provides log_likelihood or p.
--Call \ref trace_path if needed.
--Go from a single vector to a full apop_data set and back (via apop_data_pack/unpack)
--Check the derivative function if available.
--Check constraints.
*/

static double negshell (const gsl_vector *beta, void * in){
    infostruct *i = in;
    double penalty = 0,
           out     = 0; 
    double (*f)(apop_data *, apop_model *);
    f = i->model->log_likelihood? i->model->log_likelihood : i->model->p;
    Apop_stopif(!f, longjmp(i->bad_eval_jump, -1),
                0, "The model you sent to the MLE function has neither log_likelihood element nor p element.");
    apop_data_unpack(beta, i->model->parameters);
	if (i->use_constraint && i->model->constraint)
		penalty	= i->model->constraint(i->data, i->model);
    if (penalty) apop_data_pack(i->model->parameters, (gsl_vector*) beta, .all_pages='y');
    double f_val = f(i->data, i->model);
    out = penalty - f_val; //negative llikelihood
    Apop_stopif(gsl_isnan(out), longjmp(i->bad_eval_jump, -1),
                0, "I got a NaN in evaluating the objective function.%s", 
                    !i->model->constraint ? " Maybe add a constraint to your model?" : "");
    if (i->trace_path && strlen(i->trace_path))
        tracepath(i->model->parameters->vector,-out, i->trace_path, i->trace_file);
    if (i->want_info =='y'){
        //I report the log likelihood under the assumption that the final param set 
        //matches the best ll evaluated.
        double this_ll = i->model->log_likelihood? -out : log(-out); //negative negative llikelihood.

        if(gsl_isnan(this_ll)){
            Apop_stopif(!i->model->log_likelihood && penalty > f_val, /*continue*/,
                            0, "Your model's p evaluates as %g, and your penalty is %g, for an "
                               "adjusted p of %g. Please make sure that this is positive, perhaps by "
                               "rescaling your penalty.\n", f_val, penalty, f_val-penalty);
            Apop_stopif(1, apop_data_show(i->model->parameters); longjmp(i->bad_eval_jump, -1),
                        0, "NaN resulted from the following value tried by the maximum likelihood system. "
                           "Tighten your constraint? Log of a negative p?\n");
        }
        i->best_ll = GSL_MAX(i->best_ll, this_ll);
    }
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
    infostruct *i = in;
    apop_mle_settings *mp =  apop_settings_get_group(i->model, apop_mle);
    apop_data_unpack(beta, i->model->parameters);
    /* In all cases, negshell gets called first, so the constraint is already
       checked and beta nudged accordingly.
    if(i->model->constraint && i->model->constraint(i->data, i->model))
            apop_data_pack(i->model->parameters, (gsl_vector *) beta, .all_pages='y'); */
    if (i->model->score)
        i->model->score(i->data, g, i->model);
    else {
        apop_fn_with_params ll  = i->model->log_likelihood ? i->model->log_likelihood : i->model->p;
        apop_internal_numerical_gradient(ll, i, g, mp->delta);
    }
    if (i->trace_path && strlen(i->trace_path))
        negshell (beta,  in);
    gsl_vector_scale(g, -1);
    return GSL_SUCCESS;
}

//This is just to satisfy the GSL's format.
static void fdf_shell(const gsl_vector *beta, void *i, double *f, gsl_vector *df){
    *f	= negshell(beta, i);
    dnegshell(beta, i, df);
}

static int ctrl_c;
static void mle_sigint(){ ctrl_c ++; }

static int setup_starting_point(apop_mle_settings *mp, gsl_vector *x){
    Apop_stopif(!x, return -1, 0, "The vector I'm trying to optimize over is NULL.");
	if (!mp->starting_pt) gsl_vector_set_all (x, 1);
	else for (int i=0; i< x->size; i++)
            x->data[i] = mp->starting_pt[i];
    return 0;
}

static void auxinfo(apop_data *params, infostruct *i, int status, double ll){
    Get_vmsizes(params); //tsize = total # of parameters
    apop_model	*est    = i->model; //just an alias.
    /* This catches too many near-misses
       if(est->constraint)
        apop_assert(!est->constraint(i->data, est), "the maximum likelihood search ended "
                                            "at a point that doesn't satisfy the model's constraints.");*/
    if (i->want_cov=='y' && est->parameters->vector && !est->parameters->matrix){
        apop_model_numerical_covariance(i->data, est, Apop_settings_get(est,apop_mle,delta));
        if (i->want_tests=='y')
            apop_estimate_parameter_tests (est);
    }
    int param_ct = tsize;
    if (i->want_info=='y'){
        //Did the sending function save last value of f()?
        if (!ll) ll = apop_log_likelihood(i->data, i->model);

        if (!est->info) est->info = apop_data_alloc();
        apop_data_add_named_elmt(est->info, "status", status);
        apop_data_add_named_elmt(est->info, "log likelihood", ll);
        apop_data_add_named_elmt(est->info, "AIC", 2*param_ct - 2 *ll);
        if (i->data){//some models have NULL data.
            Get_vmsizes(i->data); //vsize, msize1, tsize
            apop_data_add_named_elmt(est->info, "BIC by row", param_ct * log(msize1 ? msize1: vsize) - 2 *ll);
            apop_data_add_named_elmt(est->info, "BIC by item", param_ct * log(tsize) - 2 *ll);
        }
    }
}

static apop_model *	apop_maximum_likelihood_w_d(apop_data * data, infostruct *i){
/* The maximum likelihood calculations, given a derivative of the log likelihood.

If no derivative exists, will calculate a numerical gradient.

Inside the infostruct, you'll find these elements:

\param data	the data matrix
\param	dist	the \ref apop_model object: probit, zipf, &c.
\param	starting_pt	an array of doubles suggesting a starting point. If NULL, use a vector whose elements are all 0.1 (zero has too many pathological cases).
\param step_size	the initial step size.
\param tolerance	the precision the minimizer uses. Only vaguely related to the precision of the actual var.
\return	an \ref apop_model with the parameter estimates, &c. If returned_estimate->status == 0, then optimum parameters were found; if status != 0, then there were problems.
*/
    gsl_multimin_fdfminimizer *s;
    apop_model		*est = i->model; //just an alias.
    apop_mle_settings *mp  = apop_settings_get_group(est, apop_mle);
    int iter 	= 0, 
	    status  = 0,
	    apopstatus  = 0,
	    betasize= i->beta->size;
    Apop_assert(mp, "No apop_mle settings group in the working model");
    if (mp->method == APOP_CG_BFGS)
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, betasize);
    else if (mp->method == APOP_CG_PR)
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_pr, betasize);
    else //Default:    APOP_CG_FR      conjugate gradient (Fletcher-Reeves)
	    s	= gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, betasize);
    gsl_multimin_function_fdf minme = {
        .f		= negshell,
        .df	    = (apop_df_with_void) dnegshell,
        .fdf	= (apop_fdf_with_void) fdf_shell,
        .n		= betasize,
        .params	= i};
    ctrl_c = 0;
	gsl_multimin_fdfminimizer_set (s, &minme, i->beta, mp->step_size, mp->tolerance);
    signal(SIGINT, mle_sigint);
    do { 	
        iter++;
        if (setjmp(i->bad_eval_jump)) {
            apopstatus = -1;
            break;
        }
        status 	= gsl_multimin_fdfminimizer_iterate(s);
        if(status && status!=GSL_CONTINUE) break; //commented out error msg because too many GSL_ENOPROG false positives.
        //Apop_stopif(status && status!=GSL_CONTINUE, break, 0, "GSL error: %s", gsl_strerror(status));
        status = gsl_multimin_test_gradient(s->gradient,  mp->tolerance);
        if(status && status!=GSL_CONTINUE) break; //commented out error msg because too many GSL_ENOPROG false positives.
        //Apop_stopif(status && status!=GSL_CONTINUE, break, 0, "GSL error: %s", gsl_strerror(status));
        if (mp->verbose)
            printf ("%5i %.5f  f()=%10.5f gradient=%.3f\n", iter, gsl_vector_get (s->x, 0),  s->f, gsl_vector_get(s->gradient,0));
        Apop_stopif(status == GSL_SUCCESS, apopstatus=0, 2, "Optimum found.");
    } while (status == GSL_CONTINUE && iter < mp->max_iterations && !ctrl_c);
    signal(SIGINT, NULL);
	Apop_stopif(iter==mp->max_iterations, apopstatus = -1, 1, "Max iterations reached, implying that I did not find an optimum.");
	//Clean up, copy results to output estimate.
    apop_data_unpack(s->x, est->parameters);
	gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(i->beta);
    auxinfo(est->parameters, i, apopstatus, i->best_ll);
	return est;
}

/* See apop_maximum_likelihood_w_d for notes. */
static apop_model *	apop_maximum_likelihood_no_d(apop_data * data, infostruct * i){
    apop_model *est = i->model;
    apop_mle_settings *mp = apop_settings_get_group(est, apop_mle);
    Apop_assert(mp, "No apop_mle settings group in the working model");
    int status,
        apopstatus = 0,
        iter = 0,
        betasize= i->beta->size;
    gsl_multimin_fminimizer *s;
    gsl_vector *ss;
    double size;
    s = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, betasize);
    ss = gsl_vector_alloc(betasize);
    ctrl_c      =
    apopstatus = 0; //assume failure until we score a success.
    gsl_vector_set_all (ss,  mp->step_size);
    gsl_multimin_function  minme = {.f = negshell, .n= betasize, .params = i};
    gsl_multimin_fminimizer_set (s, &minme, i->beta,  ss);
    //i->beta = s->x;
    signal(SIGINT, mle_sigint);
    do {  
        iter++;
        if (setjmp(i->bad_eval_jump)) {
            apopstatus = -1;
            break;
        }
    status  = gsl_multimin_fminimizer_iterate(s);
    if (status)  break; 
    size = gsl_multimin_fminimizer_size(s);
    status  = gsl_multimin_test_size (size, mp->tolerance); 
    if(mp->verbose){
        printf ("%5d ", iter);
        for (size_t j = 0; j < betasize; j++) 
            printf ("%8.3e ", gsl_vector_get (s->x, j)); 
        printf ("f()=%7.3f size=%.3f\n", s->fval, size);
            if (status == GSL_SUCCESS) {
                  printf ("Optimum found at:\n");
                  printf ("%5d ", iter);
                  for (size_t j = 0; j < betasize; j++)
                      printf ("%8.3e ", gsl_vector_get (s->x, j)); 
                  printf ("f()=%7.3f size=%.3f\n", s->fval, size);
            }
        }
    } while (status == GSL_CONTINUE && iter < mp->max_iterations && !ctrl_c);
    signal(SIGINT, NULL);
	Apop_stopif(iter == mp->max_iterations && mp->verbose, /*continue*/, 
                1, "Optimization reached maximum number of iterations.");
    if (status == GSL_SUCCESS) apopstatus = 0;
    apop_data_unpack(s->x, est->parameters);
	gsl_multimin_fminimizer_free(s);
    auxinfo(est->parameters, i, apopstatus, i->best_ll);
	return est;
}

/*There is a basically standard location for the log likelihood. Search there, and if you don't
find it, then recalculate it.*/
static double get_ll(apop_data *d, apop_model *est){
    if (apop_name_find(est->info->names, "Info", 'r')){
        int index = apop_name_find(est->info->names, "log likelihood", 'r');
        if (index>-2) return apop_data_get(est->info, index);
    }
    //last resort: recalculate
    return apop_log_likelihood(d, est);
}

static apop_model * dim_cycle(apop_data *d, apop_model *est, infostruct info){
    double last_ll, this_ll = GSL_NEGINF;
    int iteration = 0;
    apop_mle_settings *mp = Apop_settings_get_group(est, apop_mle);
    double tol = mp->dim_cycle_tolerance;
    int betasize = info.beta->size;
    do {
        if (mp->verbose){
            if (!(iteration++))
                printf("Cycling toward an optimum. Listing (dim):log likelihood.\n");
            printf("Iteration %i:\n", iteration);
        }
        last_ll = this_ll;
        Apop_settings_set(est, apop_mle, dim_cycle_tolerance, 0);//so sub-estimations won't use this function.
        for (int i=0; i< betasize; i++){
            gsl_vector_set(info.beta, i, GSL_NAN);
            apop_data_unpack(info.beta, est->parameters);
            apop_model *m_onedim = apop_model_fix_params(est);
            apop_prep(d, m_onedim);
            apop_maximum_likelihood(d, m_onedim);
            gsl_vector_set(info.beta, i, m_onedim->parameters->vector->data[0]);
            apop_model *full_est = apop_model_fix_params_get_base(m_onedim);//points to est, but filled.
            this_ll = get_ll(d, full_est);//only used on the last iteration.
            if (mp->verbose) printf("(%i):%g\t", i, this_ll), fflush(NULL);
            apop_model_free(m_onedim);
        }
        if (mp->verbose) printf("\n");
    } while (fabs(this_ll - last_ll) > tol);
    Apop_settings_set(est, apop_mle, dim_cycle_tolerance, tol);
    return est;
}

void get_desires(apop_model *m, infostruct *info){
    apop_parts_wanted_settings *want = apop_settings_get_group(m, apop_parts_wanted);
    apop_mle_settings *ms = apop_settings_get_group(m, apop_mle);

    info->want_tests = (want && want->tests =='y') ? 'y' : 'n';
    info->want_cov = (info->want_tests=='y' || (want && want->covariance =='y') 
            || (!want && ms && ms->want_cov =='y'))  //mle settings' want_cov only valid w/o parts_wanted
            ? 'y' : 'n';
    info->want_info = (want && want->info =='y') ? 'y' : 'n';

    //doesn't do anything at the moment.
    info->want_predicted = (want && want->predicted =='y') ? 'y' : 'n';
}

/** The maximum likelihood calculations. All of the settings are specified by adding a
  \ref apop_mle_settings struct to your model, so see the many notes there. Notably,
  the default method is the Fletcher-Reeves conjugate gradient method, and if your model
  does not have a dlog likelihood function, then a numeric gradient will be calculated
  via \ref apop_numerical_gradient. Add a \ref apop_mle_settings group to your model
  for other methods, including the Nelder-Mead simplex and simulated annealing.

\param data	The data matrix (an \ref apop_data set).
\param	dist	The \ref apop_model object: waring, probit, zipf, &amp;c. You can add
    an \c apop_mle_settings struct to it (<tt>Apop_model_add_group(your_model, apop_mle,
    .verbose=1, .method=APOP_CG_FR, and_so_on)</tt>). So, see the \c apop_mle_settings
    documentation for the many options, such as choice of method and tuning parameters.

\return	an \ref apop_model with the parameter estimates, &c. 

\li I only look at the first page of your parameter set. If you need more, perhaps use \ref apop_data_pack .

\li Get auxiliary info via, e.g.:
\code
apop_model *est = apop_estimate(your_data, apop_probit);
int status = apop_data_get(est->info, .rowname="status");
if (status)
    //trouble
else
    //optimum found
\endcode

\li During the search for an optimum, ctrl-C (SIGINT) will halt the search, and the function will return whatever parameters the search was on at the time.
 \ingroup mle */
apop_model *apop_maximum_likelihood(apop_data * data, apop_model *dist){
    apop_mle_settings   *mp = apop_settings_get_group(dist, apop_mle);
    if (!mp) mp = Apop_model_add_group(dist, apop_mle);
    if (mp->method == APOP_UNKNOWN_ML)
        mp->method = (dist->score) ? APOP_CG_FR : APOP_SIMPLEX_NM;

    Apop_assert(dist->parameters, "Not enough information to allocate parameters over which to optimize. If this was not called from apop_estimate, did you call apop_prep first?")
    infostruct info = {.data           = data,
                       .use_constraint = 1,
                       .trace_file     = malloc(sizeof(FILE *)),
                       .model          = dist};
    get_desires(dist, &info);
    info.beta = apop_data_pack(dist->parameters, NULL, .all_pages='y');
    if (setup_starting_point(mp, info.beta)) return NULL;
    *info.trace_file = NULL;
    info.model->data = data;
    if (mp->dim_cycle_tolerance)          return dim_cycle(data, dist, info);
    if (mp->trace_path)                   info.trace_path = mp->trace_path;
	if (mp->method == APOP_SIMAN)         return apop_annealing(&info);  //below.
    else if (mp->method==APOP_SIMPLEX_NM) return apop_maximum_likelihood_no_d(data, &info);
    else if (mp->method == APOP_RF_NEWTON ||
            mp->method == APOP_RF_HYBRID_NOSCALE ||
            mp->method == APOP_RF_HYBRID) return  find_roots (info);
	//else, Conjugate Gradient:
	return apop_maximum_likelihood_w_d(data, &info);
}

/** 
  The simplest use of this function is to restart a model at the latest parameter estimates.

  \code
apop_model *m = apop_estimate(data, model_using_an_MLE_search);
for (int i=0; i< 10; i++)
    m = apop_estimate_restart(m);
apop_data_show(m);
  \endcode

By adding a line to reduce the tolerance each round [e.g., <tt>Apop_settings_set(m, apop_mle, tolerance, pow(10,-i))</tt>], you can start broad and hone in on a precise optimum.
 
You may have a new estimation method, such as first doing a coarse simulated annealing search, then a fine conjugate gradient search. When reading this example, recall that the form for adding a new settings group differs from the form for modifying existing settings:
  \code
Apop_model_add_settings(your_base_model, apop_mle, .method=APOP_SIMAN);
apop_model *m = apop_estimate(data, your_base_model);
Apop_settings_set(m, apop_mle, method, APOP_CG_PR);
m = apop_estimate_restart(m);
apop_data_show(m);
  \endcode

Only one estimate is returned, either the one you sent in or a new
one. The loser (which may be the one you sent in) is freed. That is,
there is no memory leak in the above loop.

 \param e   An \ref apop_model that is the output from a prior MLE estimation. (No default, must not be \c NULL.)
 \param copy  Another not-yet-parametrized model that will be re-estimated with (1) the same data and (2) a <tt>starting_pt</tt> as per the next setting (probably
 to the parameters of <tt>e</tt>). If this is <tt>NULL</tt>, then copy <tt>e</tt>. (Default = \c NULL)
 \param starting_pt "ep"=last estimate of the first model (i.e., its current parameter estimates); "es"= starting point originally used by the first model; "np"=current parameters of the new (second) model; "ns"=starting point specified by the new model's MLE settings. (default = "ep")
 \param boundary I test whether the starting point you give me is outside this certain bound, so I can warn you if there's divergence in your sequence of re-estimations. (default: 1e8)

\return         At the end of this procedure, we'll have two \ref
    apop_model structs: the one you sent in, and the one produced using the
    new method/scale. If the new estimate includes any NaNs/Infs, then
    the old estimate is returned (even if the old estimate included
    NaNs/Infs). Otherwise, the estimate with the largest log likelihood
    is returned.

\li This function uses the \ref designated syntax for inputs.
\ingroup mle
*/ 
APOP_VAR_HEAD apop_model * apop_estimate_restart (apop_model *e, apop_model *copy, char * starting_pt, double boundary){
    apop_model * apop_varad_var(e, NULL);
    Nullcheck_m(e, NULL);
    apop_model * apop_varad_var(copy, NULL);
    char * apop_varad_var(starting_pt, "ep");
    double apop_varad_var(boundary, 1e8);
APOP_VAR_ENDHEAD
    gsl_vector *v = NULL;
    if (!copy) copy = apop_model_copy(*e);
    apop_mle_settings* prm0 = apop_settings_get_group(e, apop_mle);
    apop_mle_settings* prm = apop_settings_get_group(copy, apop_mle);
            //copy off the old params; modify the starting pt, method, and scale
    if (!strcmp(starting_pt, "es"))
        v = apop_array_to_vector(prm0->starting_pt);
    else if (!strcmp(starting_pt, "ns")){
        int size =sizeof(prm->starting_pt)/sizeof(double);
        v = apop_array_to_vector(prm->starting_pt, size);
        prm0->starting_pt	= malloc(sizeof(double)*size);
        memcpy(prm0->starting_pt, prm->starting_pt, sizeof(double)*size);
    }
    else if (!strcmp(starting_pt, "np")){
        v = apop_data_pack(copy->parameters, NULL, .all_pages='y'); 
        prm->starting_pt = malloc(sizeof(double)*v->size);
        memcpy(prm->starting_pt, v->data, sizeof(double)*v->size);
    }
    else if (e->parameters){//"ep" or default.
        v = apop_data_pack(e->parameters, NULL, .all_pages='y'); 
        prm->starting_pt = malloc(sizeof(double)*v->size);
        memcpy(prm->starting_pt, v->data, sizeof(double)*v->size);
    }
    Apop_stopif(!apop_vector_bounded(v, boundary), return e, 
                0, "Your model has diverged (element(s) > %g);"
                   " returning your original model without restarting.", boundary);
    gsl_vector_free(v);
        
    apop_model *newcopy = apop_estimate(e->data, *copy);
    apop_model_free(copy);
    //Now check whether the new output is better than the old
    if (apop_vector_bounded(newcopy->parameters->vector, boundary) 
            && get_ll(e->data, newcopy) > get_ll(e->data, e)){
        apop_model_free(e);
        return newcopy;
    } //else:
    apop_model_free(newcopy);
    return e;
}

// Simulated Annealing.

/** \page simanneal Notes on simulated annealing

Simulated annealing is a controlled random walk.  As with the other methods, the system tries a new point, and if it is better, switches. Initially, the system is allowed to make large jumps, and then with each iteration, the jumps get smaller, eventually converging. Also, there is some decreasing probability that if the new point is {\em less} likely, it will still be chosen. Simulated annealing is best for situations where there may be multiple local optima. Early in the random walk, the system can readily jump from one to another; later it will fine-tune its way toward the optimum. The number of points tested is basically not dependent on the function: if you give it a 4,000 step program, that is basically how many steps it will take.  If you know your function is globally convex (as are most standard probability functions), then this method is overkill.

The GSL's simulated annealing system doesn't actually do very much. It basically provides a for loop that calls a half-dozen functions that we the users get to write. So, the file \ref apop_mle.c handles all of this for you. The likelihood function is taken from the model, the metric is the Manhattan metric, the copy/destroy functions are just the usual vector-handling fns., et cetera. The reader who wants further control is welcome to override these functions.

Verbosity: if ep->verbose==1, show likelihood,  temp, &c. in a table; if ep->verbose>1, show that plus the vector of params.

\ingroup mle
*/

static double annealing_energy(void *in) {
    infostruct *i = in;
    return negshell(i->beta, i);
}

static double annealing_distance(void *xin, void *yin) {
/** We use the Manhattan metric to correspond to the annealing_step fn below.  */
    gsl_vector *from = apop_vector_copy(((infostruct*)xin)->beta);
    gsl_vector *to = apop_vector_copy(((infostruct*)yin)->beta);
    gsl_vector_div(from, ((infostruct*)xin)->starting_pt);
    gsl_vector_div(to, ((infostruct*)xin)->starting_pt);//starting pts are the same.
    return apop_vector_distance(from, to, .metric='m');
}

static void annealing_step(const gsl_rng * r, void *in, double step_size){
/** The algorithm: 
    --randomly pick dimension
    --shift by some amount of remaining step size
    --repeat for all dims
This will give a move \f$\leq\f$ step_size on the Manhattan metric.
*/
    infostruct *i = in;
    int sign;
    double amt, scale;
    double cutpoints[i->beta->size+1];
    cutpoints[0]             = 0;
    cutpoints[i->beta->size] = 1;
    for (size_t j=1; j< i->beta->size; j++)
        cutpoints[j] = gsl_rng_uniform(r);

    for (size_t j=0; j< i->beta->size; j++){
        sign  = (gsl_rng_uniform(r) > 0.5) ? 1 : -1;
        scale = gsl_vector_get(i->starting_pt, j);
        amt   = cutpoints[j+1]- cutpoints[j];
        apop_vector_increment(i->beta, j,  amt * sign * scale * step_size); 
    }
    apop_data_unpack(i->beta, i->model->parameters);
    if (i->model->constraint && i->model->constraint(i->data, i->model))
        apop_data_pack(i->model->parameters, i->beta, .all_pages='y');
}

static void annealing_print(void *xp) {
    apop_vector_show(((infostruct*)xp)->beta);
}

static void annealing_print2(void *xp) { return; }

static void annealing_memcpy(void *xp, void *yp){
    infostruct *yi = yp;
    infostruct *xi = xp;
    *yi = *xi;
    yi->beta = apop_vector_copy(xi->beta);
}

static void *annealing_copy(void *xp){
    infostruct *out = malloc(sizeof(infostruct));
    annealing_memcpy(xp, out);
    return out;
}

static void annealing_free(void *xp){
    gsl_vector_free(((infostruct*)xp)->beta);
    free(xp);
}

//I abuse the starting point element to hold the list of scaling factors. They can't be zero.
static double set_start(double in){ return in ? in : 1; }

jmp_buf anneal_jump;
static void anneal_sigint(){ longjmp(anneal_jump,1); }

static apop_model * apop_annealing(infostruct *i){
    apop_model *ep = i->model;
    apop_mle_settings *mp = apop_settings_get_group(ep, apop_mle);
    Apop_assert(mp, "The model you sent to the MLE function has neither log_likelihood element nor p element.");
    gsl_siman_params_t simparams = (gsl_siman_params_t) {
                         .n_tries       = mp->n_tries, 
                         .iters_fixed_T = mp->iters_fixed_T,
                         .step_size     = mp->step_size,
                         .k             = mp->k,
                         .t_initial     = mp->t_initial,
                         .mu_t          = mp->mu_t,
                         .t_min         = mp->t_min};
    static const gsl_rng *r;
    if (!r) r = mp->rng ? mp->rng : apop_rng_alloc(apop_opts.rng_seed++);
    //these two are done at apop_maximum_likelihood:
    //i->beta = apop_data_pack(ep->parameters, NULL, .all_pages='y');
    //setup_starting_point(mp, i->beta);
    int betasize = i->beta->size;
    int apopstatus = -1;
    i->starting_pt    = apop_vector_map(i->beta, set_start);
    i->use_constraint = 0; //negshell doesn't check it; annealing_step does.
    gsl_siman_print_t printing_fn = NULL;
    if (mp && mp->verbose>1)    printing_fn = annealing_print;
    else if (mp && mp->verbose) printing_fn = annealing_print2;
    if (setjmp(i->bad_eval_jump)) {
        apopstatus = -1;
        goto done;
    }
    if (!setjmp(anneal_jump)){
        signal(SIGINT, anneal_sigint);
        gsl_siman_solve(r,    // const gsl_rng * r
          i,                  // void * x0_p
          annealing_energy,   // gsl_siman_Efunc_t Ef
          annealing_step,     // gsl_siman_step_t take_step
          annealing_distance, // gsl_siman_metric_t distance
          printing_fn,        // gsl_siman_print_t print_position
          annealing_memcpy,   // gsl_siman_copy_t copyfunc
          annealing_copy,     // gsl_siman_copy_construct_t copy_constructor
          annealing_free,     // gsl_siman_destroy_t destructor
          betasize,           // size_t element_size
          simparams);         // gsl_siman_params_t params
    }
    signal(SIGINT, NULL);
    apop_data_unpack(i->beta, i->model->parameters); 
    apop_estimate_parameter_tests(i->model);
    apopstatus = 0;
done:
    if (mp->rng) r = NULL;
    auxinfo(i->model->parameters, i, apopstatus, i->best_ll);
    return i->model;
}

/* This function calls the various GSL root-finding algorithms to find the zero of the score.
   Cut/pasted/modified from the GSL documentation.  */
static apop_model * find_roots (infostruct p) {
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    apop_model *dist = p.model;
    apop_mle_settings *mlep = apop_settings_get_group(dist, apop_mle);
    int status, betasize = p.beta->size,
              apopstatus = -1;   //assume failure until we score a success.
    size_t iter = 0;
    gsl_multiroot_function f = {dnegshell, betasize, &p};
    T =   (mlep->method == APOP_RF_NEWTON)         ? gsl_multiroot_fsolver_dnewton
        : (mlep->method == APOP_RF_HYBRID_NOSCALE) ? gsl_multiroot_fsolver_hybrids
                                                   : gsl_multiroot_fsolver_hybrid;
    s = gsl_multiroot_fsolver_alloc (T, betasize);
    gsl_multiroot_fsolver_set (s, &f, p.beta);
    do {
        iter++;
        if (setjmp(p.bad_eval_jump)) break;
        status = gsl_multiroot_fsolver_iterate (s);
        if (!mlep || mlep->verbose)
            printf ("iter = %3zu x = % .3f f(x) = % .3e\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->f, 0));
        if (status)   /* check if solver is stuck */
            break;
        status = gsl_multiroot_test_residual (s->f, mlep->tolerance);
    } while (status == GSL_CONTINUE && iter < mlep->max_iterations);
    if (GSL_SUCCESS) apopstatus = 0;
    Apop_notify(2, "status = %s\n", gsl_strerror(status));
    apop_data_unpack(s->x, dist->parameters);
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (p.beta);
    auxinfo(dist->parameters, &p, apopstatus, 0); //root-finders don't store best val.
    return dist;
}
