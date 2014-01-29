/* OLS models. Much of the real work is done in apop_regression.c.
Copyright (c) 2005--2007, 2010 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.

\amodel apop_wls The Weighed Least Squares model
This is a (deprecated) synonym for \ref apop_ols, qv.  If you use the \ref apop_ols
model and provide weights in \c your_input_data->weights, then I will use them
appropriately. That is, the \ref apop_ols model really implements Weighted Least Squares,
but in most cases <tt>weights==NULL</tt> and the math reduces to the special case of
Ordinary Least Squares.

\amodel apop_ols Ordinary least squares. Weighted least squares is also handled by this model.
You can also use it for a lot of not-entirely linear models based on the form \f$Y = f(x_1) + f(x_2) + ... + \epsilon\f$.

\adoc    Input_format  See \ref dataprep.
\adoc    Parameter_format  A vector of OLS coefficients. coeff. zero
                         refers to the constant column, if any. 
\adoc    estimated_parameter_model  For the mean, a noncentral \f$t\f$ distribution (\ref apop_t_distribution).
\adoc    Prep_routine      Focuses on the data shunting. 

\adoc    settings  \ref apop_lm_settings 
\adoc    Examples
First, you will need a file named <tt>data</tt> in comma-separated form. The first column is the dependent variable; the remaining columns are the independent. For example:
\code
Y, X_1, X_2, X_3
2,3,4,5
1,2,9,3
4,7,9,0
2,4,8,16
1,4,2,9
9,8,7,6
\endcode

The program:
\include ols.c

If you saved this code to <tt>sample.c</tt>, then you can compile it with
\code
gcc sample.c -std=gnu99 -lapophenia -lgsl -lgslcblas -lsqlite3 -o run_me
\endcode

and then run it with <tt>./run_me</tt>. Alternatively, you may prefer to compile the program using a \ref makefile .

Feeling lazy? The program above was good form and demonstrated useful features, but the code below will do the same thing in two lines:

\code
#include <apop.h>
int main(){ apop_model_show(apop_estimate(apop_text_to_data("data"), apop_ols)); }
\endcode
*/

#include "apop_internal.h"

static void ols_score(apop_data *d, gsl_vector *gradient, apop_model *p);
apop_model *ols_param_models(apop_data *d, apop_model *m);
apop_data *ols_predict(apop_data *in, apop_model *m);
void ols_print(apop_model *m, FILE *ap);

Apop_settings_copy(apop_lm,
    out->instruments = apop_data_copy(in->instruments);
    if (in->input_distribution)
        out->input_distribution = apop_model_copy(in->input_distribution);
)

Apop_settings_free(apop_lm,
    apop_model_free(in->input_distribution);
) 

Apop_settings_init(apop_lm,
    if (out->want_cov == 1 || !out->want_cov) out->want_cov = 'y';
    if (out->want_expected_value == 1 || !out->want_expected_value) out->want_expected_value = 'y';
    if (!out->input_distribution) 
       out->input_distribution = apop_model_copy(apop_improper_uniform);
)

//shift first col to depvar, rename first col "one".
static void prep_names (apop_model *e){
    apop_lm_settings *p = apop_settings_get_group(e, apop_lm);
    apop_parts_wanted_settings *pwant = apop_settings_get_group(e, apop_parts_wanted);
    apop_data *predicted = apop_data_get_page(e->info, "<Predicted>");
    if (predicted){
        apop_name_add(predicted->names, (e->data->names->colct ? e->data->names->col[0] : "Observed"), 'c');
        apop_name_add(predicted->names, "Predicted", 'c');
        apop_name_add(predicted->names, "Residual", 'c');
    }
	if (e->data->names->vector) { //this is post ols shuffle.
        if (e->parameters)
            Asprintf(&e->parameters->names->title, "Regression of %s", e->data->names->vector);
        apop_name_add(e->parameters->names, "parameters", 'v');
        for(int i=0; i< e->data->names->colct; i++)
            apop_name_add(e->parameters->names, e->data->names->col[i], 'r');
        if ((pwant && pwant->covariance) || (!pwant && p && p->want_cov== 'y')){
            apop_data *cov = apop_data_get_page(e->parameters, "<Covariance>");
            if (cov && e->data->names){
                apop_name_stack(cov->names, e->data->names, 'c');
                apop_name_stack(cov->names, e->data->names, 'r', 'c');
            }
        }
	}
}

static void ols_shuffle(apop_data *d){
    if (!d) return;
    if (!d->vector){
        Apop_col_v(d, 0, independent);
        d->vector = apop_vector_copy(independent);
        gsl_vector_set_all(independent, 1);     //affine; first column is ones.
        if (d->names->colct > 0) {		
            apop_name_add(d->names, d->names->col[0], 'v');
            sprintf(d->names->col[0], "1");
        }
    }
}

static void ols_prep(apop_data *d, apop_model *m){
    apop_score_vtable_add(ols_score, apop_ols);
    apop_parameter_model_vtable_add(ols_param_models, apop_ols);
    apop_predict_vtable_add(ols_predict, apop_ols);
    apop_model_print_vtable_add(ols_print, apop_ols);
    if (!d) return;
    ols_shuffle(d);
    void *mpt = m->prep; //also use the defaults.
    m->prep = NULL;
    apop_prep(d, m);
    m->prep = mpt;
}

/* The assumption that makes a log likelihood possible is that the
errors are normally distributed.

This function is a bit inefficient, in that it calculates the error terms,
which you may have already done in the OLS estimation.  */
static long double ols_log_likelihood (apop_data *d, apop_model *p){ 
    Nullcheck_mpd(d, p, GSL_NAN); Nullcheck(d->matrix, GSL_NAN);
  long double ll = 0; 
  long double sigma, actual, weight;
  double expected, x_prob;
  apop_lm_settings *lms = Apop_settings_get_group(p, apop_lm);
  apop_model *input_distribution = lms ? lms->input_distribution : NULL;
  gsl_matrix *data = d->matrix;
  gsl_vector *errors = gsl_vector_alloc(data->size1);
	for (size_t i=0;i< data->size1; i++){
        Apop_row_v(d, i, datarow);
        gsl_blas_ddot(p->parameters->vector, datarow, &expected);
        if (d->vector){ //then this has been prepped
            actual = apop_data_get(d,i, -1);
        } else {
            actual = gsl_matrix_get(data,i, 0);
            expected += gsl_vector_get(p->parameters->vector,0) * (1 - actual); //data isn't affine.
        }
        gsl_vector_set(errors, i, expected-actual);
    }
    sigma = sqrt(apop_vector_var(errors));
	for(size_t i=0; i< data->size1; i++){
        Apop_row(d, i, justarow);
        justarow->vector = NULL;
        x_prob = (input_distribution)
                    ? apop_p(justarow, input_distribution) //probably improper uniform, and so just 1 anyway.
                    : 1;
        weight = d->weights ? gsl_vector_get(d->weights, i) : 1; 
        ll += logl(gsl_ran_gaussian_pdf(gsl_vector_get(errors, i), sigma)* weight * x_prob);
	} 
    gsl_vector_free(errors);
    return ll;
}

/* $\partial {\cal N}(x\beta - y)/\partial \beta_i = \sum{x_i} \partial {\cal N}(K)/\partial K$ (at $K=x\beta -y$) */
static void ols_score(apop_data *d, gsl_vector *gradient, apop_model *p){ 
    Nullcheck_mpd(d, p, ); Nullcheck(d->matrix, );
  long double sigma, actual, weight;
  double expected;
  gsl_matrix *data	= d->matrix;
  gsl_vector *errors = gsl_vector_alloc(data->size1);
  gsl_vector *normscore = gsl_vector_alloc(2);
  apop_data  *subdata  = apop_data_alloc(1,1);
	for(size_t i=0;i< data->size1; i++){
        Apop_row_v(d, i, datarow);
        gsl_blas_ddot(p->parameters->vector, datarow, &expected);
        if (d->vector){ //then this has been prepped
            actual       = apop_data_get(d,i, -1);
        } else {
            actual       = gsl_matrix_get(data,i, 0);
            expected    +=  gsl_vector_get(p->parameters->vector,0) * (1 - actual); //data isn't affine.
        }
        gsl_vector_set(errors, i, expected-actual);
    }
    sigma   = sqrt(apop_vector_var(errors));
    apop_model *norm = apop_model_set_parameters(apop_normal, 0.0, sigma);
    gsl_vector_set_all(gradient, 0);
	for(size_t i=0;i< data->size1; i++){
        apop_data_set(subdata, 0, 0, gsl_vector_get(errors, i));
        apop_score(subdata, normscore, norm);
        weight = d->weights ? gsl_vector_get(d->weights, i) : 1; 
        for(size_t j=0; j< data->size2; j++)
            *gsl_vector_ptr(gradient, j) += weight * apop_data_get(d, i, j) * gsl_vector_get(normscore, 0);
	} 
    gsl_vector_free(errors);
    apop_model_free(norm);
}

//xpx may be destroyed by the HH transformation.
static void xpxinvxpy(apop_data const*data, gsl_matrix *xpx, apop_data const* xpy, apop_model *out){
    apop_lm_settings   *p =  apop_settings_get_group(out, apop_lm);
    apop_parts_wanted_settings *pwant = apop_settings_get_group(out, apop_parts_wanted);
	if ( (pwant && pwant->covariance!='y' && pwant->predicted != 'y') 
       ||(!pwant && p && p->want_cov!='y' && p->want_expected_value != 'y')){	
		//then don't calculate (X'X)^{-1}
		gsl_linalg_HH_solve (xpx, xpy->vector, out->parameters->vector);
		return;
	} //else:
    double s_sq;
    gsl_vector const *y_data = data->vector; //just an alias
    apop_data *cov = apop_data_alloc();
    double det = apop_det_and_inv(xpx, &cov->matrix, 1, 1);// not yet cov, just (X'X)^-1.
    if (det < 1e-4) Apop_notify(1, "Determinant of X'X is small (%g), so matrix is near singular. "
                        "Expect the covariance matrix [based on (X'X)^-1] to be garbage.", det);
    apop_data_free(out->parameters);
    out->parameters = apop_dot(cov, xpy);               // \beta=(X'X)^{-1}X'Y
    apop_data *error = apop_dot(data, out->parameters); // X\beta ==predicted (not yet error)
	gsl_vector_sub(error->vector, y_data);              // X'\beta - Y == error
    gsl_blas_ddot(error->vector, error->vector, &s_sq); // e'e
    s_sq /= data->matrix->size1 - data->matrix->size2;  // \sigma^2 = e'e / df
	gsl_matrix_scale(cov->matrix, s_sq);                // cov = \sigma^2 (X'X)^{-1}
	if ((pwant && pwant->predicted) || (!pwant && p && p->want_expected_value)){
        apop_data *predicted_page = apop_data_get_page(out->info, "<Predicted>");
        gsl_matrix_set_col(predicted_page->matrix, 0, y_data);
        gsl_matrix_set_col(predicted_page->matrix, 2, error->vector);
        Apop_col_v(predicted_page, 1, predicted);
        gsl_vector_memcpy(predicted, y_data);
        gsl_vector_add(predicted, error->vector); //pred = y_data + error
    }
    apop_data_free(error);
    if (apop_data_get_page(out->parameters, "<Covariance>"))
        apop_data_rm_page(out->parameters, "<Covariance>");
    apop_data_add_page(out->parameters, cov, "<Covariance>");
}

/* \adoc    RNG  Linear models are typically only partially defined probability models. For
OLS, we know that \f$P(Y|X\beta) \sim {\cal N}(X\beta, \sigma)\f$, because this is
an assumption about the error process, but we don't know much of anything about the
distribution of \f$X\f$.

The \ref apop_lm_settings group includes an \ref apop_model* element named \c
input_distribution. This is the distribution of the independent/predictor/X columns
of the data set.

The default is that <tt>input_distribution = apop_improper_uniform </tt>, meaning that
\f$P(X)=1\f$ for all \f$X\f$. So \f$P(Y, X) = P(Y|X)P(X) = P(Y|X)\f$. This seems to
be how many people use linear models: the \f$X\f$ values are taken as certain (as with
actually observed data) and the only question is the odds of the dependent variable. If
that's what you're looking for, just leave the default. This is sufficient for getting
log likelihoods under the typical assumption that the observed data has probability one.

<em>But</em> you can't draw from an improper uniform. So if you draw from a linear
model with a default <tt>input_distribution</tt>, then you'll get an error.

Alternatively, you may know something about the distribution of the input data. At
     the least, you could generate a PMF from the actual data:
     \code
    apop_settings_set(your_model, apop_lm, input_distribution, apop_estimate(inset, apop_pmf));
     \endcode
Now, random draws are taken from the input data, and the dependent variable value calculated via \f$X\beta+\epsilon\f$, where \f$X\f$ is the drawn value, \f$\beta\f$ the previously-estimated parameters and \f$\epsilon\f$ is a Normally-distributed random draw. Or change the PMF to any
other appropriate distribution, such as a \ref apop_multivariate_normal,
or an \ref apop_pmf filled in with more data, or perhaps something from
http://en.wikipedia.org/wiki/Errors-in-variables_models , as desired.  */
static int ols_rng(double *out, gsl_rng *r, apop_model *m){
    //X is drawn from the input distribution, then Y = X\beta + epsilon
    apop_lm_settings *olp =  apop_settings_get_group(m, apop_lm);
    Apop_stopif(!olp, return 1, 0, "no apop_lm settings group attached. Has this model been estimated yet?");

    gsl_vector *tempdata = gsl_vector_alloc(m->parameters->vector->size);
    Apop_stopif(apop_draw(tempdata->data, r, olp->input_distribution), return 2,
            0, "Couldn't draw from the distribution of the input data.");
    gsl_blas_ddot(tempdata, m->parameters->vector, out);

    double sigma_sq = apop_data_get(m->info, .rowname="SSE")/m->data->matrix->size1;
    out[0] += gsl_ran_gaussian(r, sqrt(sigma_sq));

    if (m->dsize > 1) memcpy(out+1, tempdata->data, sizeof(double)*tempdata->size);
    gsl_vector_free(tempdata);
    return 0;
}

/* \adoc estimated_data You can specify whether the data is modified with an \ref apop_lm_settings group. If so, see \ref dataprep for details. Else, left unchanged.

\adoc estimated_parameters
The \c parameters set will hold the coefficients; the first coefficient will be the
coefficient on the constant term, and the remaining will correspond to the independent
variables. It will therefore be of size <tt>(data->size2)</tt>.

I add a page named <tt>\<Covariance\></tt>, which gives the covariance matrix for the
estimated parameters (not the data itself).

\adoc estimated_info Reports log likelihood, and runs \ref apop_estimate_coefficient_of_determination 
to add \f$R^2\f$-type information (SSE, SSR, \&c) to the info page.

Residuals: I add a page named <tt>\<Predicted\></tt>, with three columns. If this is a model
with a single dependent and lots of independent vars, then the first column is the
actual data. Let our model be \f$ Y = \beta X + \epsilon\f$. Then the second column
is the predicted values: \f$\beta X\f$, and the third column is the residuals:
\f$\epsilon\f$. The third column is therefore always the first minus the second,
and this is probably how that column was calculated internally.

Given your estimate \c est, the zeroth element is one of <br> 
<tt> apop_data_get(est->info, .page= "Predicted", .row=0, .colname="observed"),</tt><br>
<tt> apop_data_get(est->info, .page= "Predicted", .row=0, .colname="predicted") or</tt><br>
<tt> apop_data_get(est->info, .page= "Predicted", .row=0, .colname="residual").</tt><br>
*/
static void apop_estimate_OLS(apop_data *inset, apop_model *ep){
    Nullcheck_mpd(inset, ep, );
    apop_data *set;
    apop_lm_settings *olp =  apop_settings_get_group(ep, apop_lm);
    apop_parts_wanted_settings *pwant = apop_settings_get_group(ep, apop_parts_wanted);
    if (!olp) 
        olp = Apop_model_add_group(ep, apop_lm);
    ep->data = inset;
    set = olp->destroy_data ? inset : apop_data_copy(inset); 
    
    gsl_vector *weights = olp->destroy_data      //this may be NULL.
                           ? ep->data->weights 
                           : apop_vector_copy(ep->data->weights);
    if (weights)
        for (size_t i =0; i< weights->size; i++)
            gsl_vector_set(weights, i, sqrt(gsl_vector_get(weights, i)));

    if ((pwant &&pwant->predicted) || (!pwant && olp && olp->want_expected_value=='y'))
        apop_data_add_page(ep->info, apop_data_alloc(0, set->matrix->size1, 3), "<Predicted>");
    if ((pwant &&pwant->covariance) || (!pwant && olp && olp->want_cov=='y'))
        apop_data_add_page(ep->parameters, apop_data_alloc(0, set->matrix->size2, set->matrix->size2), "<Covariance>");
    if (weights)
        for (int i = -1; i < set->matrix->size2; i++){
            Apop_col_v(set, i, v);
            gsl_vector_mul(v, weights);
        }

    apop_data *xpx_d = apop_dot(set, set, .form1='t'); //(X'X)
    apop_data *xpy_d = apop_dot(set, set, .form1='t', .form2='v'); //(X'y)
    xpxinvxpy(set, xpx_d->matrix, xpy_d, ep);
    prep_names(ep);
    apop_data_free(xpx_d);
    apop_data_free(xpy_d);

    if ((pwant &&pwant->covariance) || (!pwant && olp && olp->want_cov=='y'))
        apop_estimate_parameter_tests(ep);

    add_info_criteria(ep->data, ep, ep, apop_log_likelihood(ep->data, ep)); //in apop_mle.c

    apop_data *r_sq = apop_estimate_coefficient_of_determination(ep); //Add R^2-type info to info page.
    apop_data_stack(ep->info, r_sq, .inplace='y');

    apop_data_free(r_sq);
    if (!olp->destroy_data){
        if (weights) gsl_vector_free(weights);
        apop_data_free(set);
    }
}

/* \adoc predict This function is limited to taking in a data set with a matrix, and
filling the vector with \f$X\beta\f$. Like, the OLS estimation will shuffle a matrix around
to insert a column of ones (see \ref dataprep).
 */
apop_data *ols_predict(apop_data *in, apop_model *m){
    Nullcheck_mpd(in, m, NULL);
    if (!in->vector)  ols_shuffle(in);  

    //find x dot y
    gsl_blas_dgemv (CblasNoTrans, 1, in->matrix, m->parameters->vector, 0, in->vector);
    return in;
}

apop_model *ols_param_models(apop_data *d, apop_model *m){
    Nullcheck_mpd(d, m, NULL);
    apop_pm_settings *settings = Apop_settings_get_group(m, apop_pm);
    if (settings->index!=-1){
        int i = settings->index;
        double mu = apop_data_get(m->parameters, i, -1);
        double sigma = sqrt(apop_data_get(m->parameters, i, i, .page="<Covariance>"));
        int df = apop_data_get(m->info, .rowname="df");
        return apop_model_set_parameters(apop_t_distribution, mu, sigma, df);
    }
    //else run the default
    apop_parameter_model_vtable_drop(m);
    apop_model *out = apop_parameter_model(d, m);
    apop_parameter_model_vtable_add(ols_param_models, m);
    return out;
}

void ols_print(apop_model *m, FILE *ap){
    fprintf(ap, "Parameters:\n");
    apop_data_print(m->parameters, .output_pipe=(ap? ap : stdout));
    apop_data *predict = apop_data_rm_page(m->info, "<Predicted>", .free_p='n');
    apop_data_print(m->info, .output_pipe=(ap? ap : stdout));
    apop_data_add_page(m->info, predict, predict->names->title);
}

apop_model *apop_ols = &(apop_model){.name="Ordinary Least Squares", .vsize = -1, .dsize=-1, .estimate=apop_estimate_OLS, 
            .log_likelihood = ols_log_likelihood, .prep = ols_prep, .draw=ols_rng};


/*\amodel apop_iv Instrumental variable regression

Operates much like the \ref apop_ols model, but the input parameters also need to have
a table of substitutions (like the addition of the <tt>.instruments</tt> setting in
the example below). The vector element of the table lists the column numbers to be
substituted (the dependent var is zero; first independent col is one), and then one
column for each item to substitute.

\li If the vector of your apop_data set is \c NULL, then I will use the row names to find
the columns to substitute. This is generally more robust and/or convenient.

\li If the \c instruments data set is somehow \c NULL or empty, I'll just run OLS. 

\li Don't forget that the \ref apop_lm_settings group has a \c destroy_data setting. If
you set that to \c 'y', I will overwrite the column in place, saving the trouble of
copying the entire data set.

\adoc    Input_format  See \ref apop_ols; see \ref dataprep. 
\adoc    Parameter_format  As per \ref apop_ols 
\adoc    Estimate_results  As per \ref apop_ols 
\adoc    Prep_routine  Focuses on the data shunting. 
\adoc    settings  \ref apop_lm_settings 
\adoc Examples 
\include  iv.c
*/

static apop_data *prep_z(apop_data *x, apop_data *instruments){
    apop_data *out = apop_data_copy(x);
    if (instruments->vector)
        for (int i=0; i< instruments->vector->size; i++){
            Apop_col_v(instruments, i, inv);
            Apop_col_v(out, instruments->vector->data[i], outv);
            gsl_vector_memcpy(outv, inv);
        }
    else if (instruments->names->colct)
        for (int i=0; i< instruments->names->colct; i++){
            int colnumber = apop_name_find(x->names, instruments->names->col[i], 'c');
            Apop_assert(colnumber != -2, "You asked me to substitute instrument column %i "
                    "for the data column named %s, but I could find no such name.",  i, instruments->names->col[i]);
            Apop_col_v(instruments, i, inv);
            Apop_col_v(out, colnumber, outv);
            gsl_vector_memcpy(outv, inv);
        }
    else Apop_assert(0, "Your instrument matrix has data, but neither a vector element "
                       "nor column names indicating what columns in the original data should be replaced.");
    return out;
}

static void apop_estimate_IV(apop_data *inset, apop_model *ep){
    Nullcheck_mpd(inset, ep, );
    apop_lm_settings   *olp =  apop_settings_get_group(ep, apop_lm);
    apop_parts_wanted_settings *pwant = apop_settings_get_group(ep, apop_parts_wanted);
    if (!olp) olp = Apop_model_add_group(ep, apop_lm);
    if (!olp->instruments || !(olp->instruments->matrix || olp->instruments->vector)) 
        apop_ols->estimate(inset, ep);
    ep->data = inset;
    if (ep->parameters) apop_data_free(ep->parameters);
    ep->parameters = apop_data_alloc(inset->matrix->size2);
    apop_data *set = olp->destroy_data ? inset : apop_data_copy(inset); 
    apop_data *z = prep_z(inset, olp->instruments);
    
    gsl_vector *weights = olp->destroy_data      //the weights may be NULL.
                             ? ep->data->weights 
                             : apop_vector_copy(ep->data->weights);
    if (weights)
        for (int i =0; i< weights->size; i++)
            gsl_vector_set(weights, i, sqrt(gsl_vector_get(weights, i)));

    if ((pwant && pwant->predicted) || (!pwant && olp && olp->want_expected_value))
        apop_data_add_page(ep->info, apop_data_alloc(set->matrix->size1, 3), "<Predicted>");
    prep_names(ep);
    if (weights){
        gsl_vector_mul(set->vector, weights);
        for (int i = 0; i < set->matrix->size2; i++){
            Apop_col_v(set, i, v);
            gsl_vector_mul(v, weights);
        }
    }

    apop_data *zpx = apop_dot(z, set, .form1='t');
    apop_data *zpy = apop_dot(z, set, .form1='t', .form2='v'); //z'y

    xpxinvxpy(inset, zpx->matrix, zpy, ep);

    //covariance matrix right now is sigma (Z'X)^-1. We need
    //sigma (Z'X)^-1 (Z'Z) (X'Z)^-1

    apop_data *zpz = apop_dot(z, z, .form1='t');
    apop_data zpxinv = (apop_data) {.matrix=apop_matrix_inverse(zpx->matrix)};
    apop_data *zpz_xpzinv = apop_dot(zpz, &zpxinv, .form2='t');
    apop_data *halfcov = apop_data_get_page(ep->parameters, "<Covariance>");
    apop_data *cov = apop_dot(halfcov, zpz_xpzinv);
    apop_data_rm_page(ep->parameters, "<Covariance>");
    apop_data_add_page(ep->parameters, cov, "<Covariance>");

    gsl_matrix_free(zpxinv.matrix);
    apop_data_free(zpx);
    apop_data_free(zpy);
    apop_data_free(zpz);


/*
    apop_data *zpxinv = apop_matrix_to_data(apop_matrix_inverse(zpx->matrix));
    ep->parameters = apop_dot(zpxinv, zpy);
    //cov = sigma^2 (Z'X)^-1 Z'Z (X'Z)^-1
    */

    /*
    if ((pwant &&pwant->covariance) || (!pwant && olp && olp->want_cov=='y')){
        apop_data *zpz = apop_dot(z, z, .form1='t');
        apop_data *zpz_zpxinv = apop_dot(zpz, zpxinv, .form2='t');
        apop_data_add_page(ep->parameters, apop_dot(zpx, zpz_zpxinv)
                , "<Covariance>");
        apop_data_free(zpz); apop_data_free(zpz_zpxinv);
    }
    */

    apop_data_free(zpx);// apop_data_free(zpxinv);
    apop_data_free(zpy);

    if (!olp->destroy_data) apop_data_free(set);
}

apop_model *apop_iv = &(apop_model){.name="instrumental variables", .vsize = -1, .dsize=-1,
    .estimate =apop_estimate_IV, .prep=ols_prep,
    .log_likelihood = ols_log_likelihood};
