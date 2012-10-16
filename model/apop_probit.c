/* Probit and Logit. 
Copyright (c) 2005--2008, 2010 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2. 

\amodel apop_probit The Probit model.

  Apophenia makes no distinction between the bivariate probit and the multinomial probit. This one does both.

\adoc    Input_format  
The first column of the data matrix this model expects is zeros, ones, ..., enumerating
the factors; to get there, try \ref apop_data_to_factors; if you  forget to run it,
I'll run it on the first data column for you.  The remaining columns are values of the
independent variables. Thus, the model will return [(data columns)-1]\f$\times\f$[(option
count)-1] parameters.  Column names are options; row names are input variables.

\adoc    Parameter_format  As above 
\adoc    Prep_routine You will probably want to convert some column of your data into
factors, via \ref apop_data_to_factors. If you do, then that adds a page of factors
to your data set (and of course adjusts the data itself). If I find a factor page,
I will use that info; if not, then I will run \ref apop_data_to_factors on the first
column (the vector if there is one, else the first column of the matrix.)

Also, if there is no vector, then I will move the first column of the matrix, and
replace that matrix column with a constant column of ones, just like with OLS.

\adoc    settings   None, but see above about seeking a factor page in the input data.
\adoc    RNG  See \ref apop_ols; this one is similar but produces a category number instead of OLS's continuous draw.
*/

#include "apop_internal.h"

static apop_data *get_category_table(apop_data *d){
    int first_col = d->vector ? -1 : 0;
    apop_data *out = apop_data_get_factor_names(d, .col=first_col);
    if (!out) {
        apop_data_to_factors(d, .intype='d', .incol=first_col, .outcol=first_col);
        out = apop_data_get_factor_names(d, .col=first_col);
    }
    return out;
}

static void probit_prep(apop_data *d, apop_model *m){
    apop_ols.prep(d, m);//also runs the default apop_model_clear.
    apop_data *factor_list = get_category_table(d);
    int count = factor_list->textsize[0];
    m->parameters = apop_data_alloc(d->matrix->size2, count-1);
    apop_name_stack(m->parameters->names, d->names, 'r', 'c');
    for (int i=1; i< count; i++) 
        apop_name_add(m->parameters->names, factor_list->text[i][0], 'c');
    gsl_matrix_set_all(m->parameters->matrix, 1);
    //haven't yet implemented the score for probit && $N>2$
    if (!strcmp(m->name, apop_probit.name) && count > 2) m->score = NULL;
    char *tmp = strdup(m->name);
    snprintf(m->name, 100, "%s with %s as numeraire", tmp, factor_list->text[0][0]);
    free(tmp);
}

double biprobit_ll_row(apop_data *r){
    long double n = gsl_cdf_gaussian_P(-gsl_matrix_get(r->matrix, 0, 0),1);
    n = n ? n : 1e-10; //prevent -inf in the next step.
    n = n<1 ? n : 1-1e-10; 
    return r->vector->data[0] ?  log(1-n): log(n);
}

//The case where outcome is a single zero/one option.
static double biprobit_log_likelihood(apop_data *d, apop_model *p){
    apop_data *betadotx = apop_dot(d, p->parameters); 
    betadotx->vector = d->vector;
    double total_prob = apop_map_sum(betadotx, .fn_r=biprobit_ll_row);
    betadotx->vector = NULL;
    apop_data_free(betadotx);
	return total_prob;
}

static double val;
static double unordered(double in){ return in == val; }

// This is just a for loop that runs a probit on each column.
static double multiprobit_log_likelihood(apop_data *d, apop_model *p){
    Nullcheck_mpd(d, p, GSL_NAN)
    gsl_vector *val_vector = get_category_table(d)->vector;
    if (val_vector->size==2) return biprobit_log_likelihood(d, p);
    //else, multinomial loop
    static apop_model *spare_probit = NULL;
    if (!spare_probit){
        spare_probit = apop_model_copy(apop_probit);
        spare_probit->parameters = apop_data_alloc();
    }
    Staticdef(apop_data *, working_data, apop_data_alloc());
    working_data->matrix = d->matrix;
    gsl_vector *original_outcome = d->vector;
    double ll = 0;
    double *vals = val_vector->data;
    for(size_t i=0; i < p->parameters->matrix->size2; i++){
        Apop_col(p->parameters, i, param);
        val = vals[i];
        working_data->vector = apop_vector_map(original_outcome, unordered);
        spare_probit->parameters->matrix = apop_vector_to_matrix(param);
        ll  += apop_log_likelihood(working_data, spare_probit);
        gsl_vector_free(working_data->vector); //yup. It's inefficient.
        gsl_matrix_free(spare_probit->parameters->matrix);
    }
	return ll;
}

static void probit_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
    Nullcheck_mp(p, )
    gsl_vector *val_vector = get_category_table(p->data)->vector;
    if (val_vector->size!=2){
        gsl_vector * numeric_default = apop_numerical_gradient(d, p);
        gsl_vector_memcpy(gradient, numeric_default);
        gsl_vector_free(numeric_default);
    }
    long double	cdf, betax, deriv_base;
    apop_data *betadotx = apop_dot(d, p->parameters); 
    gsl_vector_set_all(gradient,0);
    for (size_t i=0; i< d->matrix->size1; i++){
        betax = apop_data_get(betadotx, i, 0);
        cdf   = gsl_cdf_gaussian_P(-betax, 1);
        cdf = cdf ? cdf : 1e-10; //prevent -inf in the next step.
        cdf = cdf<1 ? cdf : 1-1e-10; 
        if (apop_data_get(d, i, -1))
            deriv_base      = gsl_ran_gaussian_pdf(-betax, 1) /(1-cdf);
        else
            deriv_base      = -gsl_ran_gaussian_pdf(-betax, 1) / cdf;
        for (size_t j=0; j< d->matrix->size2; j++)
            apop_vector_increment(gradient, j, apop_data_get(d, i, j) * deriv_base);
	}
	apop_data_free(betadotx);
}

apop_model apop_probit = {"Probit", .log_likelihood = multiprobit_log_likelihood, .dsize=-1,
    .score = probit_dlog_likelihood, .prep = probit_prep};


/* \amodel apop_multinomial_probit The Multinomial Probit model.

  \deprecated  Use \ref apop_probit, which handles multiple options.*/

/////////  Multinomial Logit (plain logit is a special case)

static apop_data *multilogit_expected(apop_data *in, apop_model *m){
    Nullcheck_mpd(in, m, NULL)
    gsl_matrix *params = m->parameters->matrix;
    apop_data *out = apop_data_alloc(in->matrix->size1, in->matrix->size1, params->size2+1);
    for (size_t i=0; i < in->matrix->size1; i ++){
        Apop_row(in, i, observation);
        Apop_row(out, i, outrow);
        double oneterm;
        int bestindex = 0;
        double bestscore = 0;
        gsl_vector_set(outrow, 0, 1);
        for (size_t j=0; j < params->size2+1; j ++){
            if (j == 0){
                oneterm = 0;
                gsl_vector_set(outrow, j, 1);
            } else {
                Apop_matrix_col(params, j-1, p);
                gsl_blas_ddot(observation, p, &oneterm);
                gsl_vector_set(outrow, j, exp(oneterm));
            }
            if (oneterm > bestscore){
                bestindex = j;
                bestscore = oneterm;
            }
        }
        double total = apop_sum(outrow);
        gsl_vector_scale(outrow, 1/total);
        apop_data_set(out, i, -1, bestindex);
    }
    apop_data *factor_list = get_category_table(m->data);
    apop_name_add(out->names, factor_list->text[0][0], 'c');
    apop_name_stack(out->names, m->parameters->names, 'c');
    return out;
}

static size_t find_index(double in, double *m, size_t max){
    size_t i = 0;
    while (in !=m[i] && i<max) i++;
    return i;
}

double one_logit_row(apop_data *thisobservation, void *factor_list){
    //get the $x\beta_j$ numerator for the appropriate choice:
    size_t index   = find_index(gsl_vector_get(thisobservation->vector, 0), 
                                (double*)factor_list, thisobservation->matrix->size2);
    Apop_row(thisobservation, 0, thisrow);
    double num = (index==0) ? 0 : gsl_vector_get(thisrow, index-1);

    /* Get the denominator, ln(sum(exp(xbeta))) using the subtract-the-max trick 
     mentioned in the documentation.  Don't forget the implicit beta_0, fixed at 
     zero (so we need to add exp(0-max)). */
    double max = gsl_vector_max(thisrow);
    gsl_vector_add_constant(thisrow, -max);
    apop_vector_exp(thisrow);
    return num - (max + log(apop_vector_sum(thisrow) +exp(-max)));
}

static double multilogit_log_likelihood(apop_data *d, apop_model *p){
    Nullcheck_mpd(d, p, GSL_NAN)
    Nullcheck(d->matrix, GSL_NAN)
    //Find X\beta_i for each row of X and each column of \beta.
    apop_data  *xbeta = apop_dot(d, p->parameters);
    double* factor_list = get_category_table(p->data)->vector->data;
    xbeta->vector = d->vector; //we'll need this in one_logit_row
    long double ll = apop_map_sum(xbeta, .fn_rp = one_logit_row, .param=factor_list);
    xbeta->vector = NULL;
    apop_data_free(xbeta);
	return ll;
}

static double dlogit_foreach(apop_data *x, void *gin){
  //\beta_this = choice for the row.
  //dLL/d\beta_ij = [(\beta_i==\beta_this) ? x_j : 0] - x_j e^(x\beta_i)/\sum_k e^(x\beta_k)
  //that last term simplifies: x / \sum_k e^(x(\beta_k - \beta_i))
    apop_data *gmat = gin;
    gsl_matrix *beta = gmat->more->matrix;
    apop_data* factor_list = gmat->more->more;
    Apop_row(x, 0, xdata);
    assert(gmat->matrix->size1 == x->matrix->size2);     //the j index---input vars (incl. 1 column)
    assert(gmat->matrix->size2 == beta->size2); //the i index---choices
    assert(xdata->size == beta->size1);//cols of data=variables; rows of output=var.s (cols=choices)
    size_t choice = find_index(gsl_vector_get(x->vector, 0), factor_list->vector->data, factor_list->vector->size);
    for (int i=0; i < beta->size2; i++) { //go through choices.
        gsl_vector *denom = gsl_vector_alloc(beta->size1);
        gsl_vector_set_all(denom, 1); //see below
        for (int other=0; other < beta->size2; other++){
            if (other != i) { //this block calculates exp(x (otherbeta-thisbeta))
                Apop_matrix_col(beta, i, thisbeta); 
                Apop_matrix_col(beta, other, otherbeta);
                gsl_vector *diff = apop_vector_copy(thisbeta);
                gsl_vector_scale(diff, -1);
                gsl_vector_add(diff, otherbeta);
                gsl_vector_mul(diff, xdata);
                apop_vector_exp(diff);
                gsl_vector_add(denom, diff);
                gsl_vector_free(diff);
            } //else, other==i, and \beta_i - \beta_i = 0, and e^{0x} = 1. Thus the gsl_vector_set_all(denom, 1).
        }
        for (int j=0; j< xdata->size; j++){ //add to each coefficient of the gradient matrix 
            double pick = (choice-1 == i) ? gsl_vector_get(xdata,j) : 0; //numeraire has no betas.
            apop_matrix_increment(gmat->matrix, j, i, 
                                            pick - gsl_vector_get(xdata, j)/gsl_vector_get(denom, j));
        }
        gsl_vector_free(denom);
    }
    return 0;
}

static void logit_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
    Nullcheck_mpd(d, p, );
    apop_data *gradient_matrix = apop_data_calloc(p->parameters->matrix->size1, p->parameters->matrix->size2);
    apop_data_unpack(gradient, gradient_matrix);
    gradient_matrix->more = p->parameters;  // facilitate passing to logit_foreach.
    gradient_matrix->more->more = get_category_table(d);
    apop_data_free_base( //return from apop_map is a vector of zeros.
        apop_map(d, .fn_rp=dlogit_foreach, .param=gradient_matrix, .all_pages='n')
    );
    apop_data_pack(gradient_matrix, gradient);
    gradient_matrix->more=NULL;
    apop_data_free(gradient_matrix);
}

apop_model *logit_estimate(apop_data *d, apop_model *m){
    apop_model *out = apop_maximum_likelihood(d, m);

    //That's the whole estimation. But now we need to add a row of zeros
    //for the numeraire. This is just tedious matrix-shunting.
    apop_data *p = out->parameters;
    apop_data *newparams = apop_data_calloc(p->matrix->size1, 1);//top row of zeros...
    apop_data_stack(newparams, p, 'c', .inplace='y');            //...followed by the rest
    apop_name_stack(newparams->names, p->names, 'r');

    apop_data_free(p);
    out->parameters = newparams;
    return out;
}

/* \amodel apop_logit The Logit model.

Apophenia makes no distinction between the bivariate logit and the multinomial logit. This does both.

  The likelihood of choosing item \f$j\f$ is:
  \f$e^{x\beta_j}/ (\sum_i{e^{x\beta_i}})\f$

  so the log likelihood is 
  \f$x\beta_j  - ln(\sum_i{e^{x\beta_i}})\f$

\adoc    Input_format  The first column of the data matrix this model expects is zeros,
ones, ..., enumerating the factors; to get there, try \ref apop_data_to_factors; if
you  forget to run it, I'll run it on the first data column for you.  The remaining
columns are values of the independent variables. Thus, the model will return [(data
columns)-1]\f$\times\f$[(option count)-1] parameters.  Column names are options;
row names are input variables.

\adoc    Parameter_format  As above.    
\adoc    Prep_routine You will probably want to convert some column of your data into
factors, via \ref apop_data_to_factors. If you do, then that adds a page of factors
to your data set (and of course adjusts the data itself). If I find a factor page,
I will use that info; if not, then I will run \ref apop_data_to_factors on the first
column (the vector if there is one, else the first column of the matrix.)

Also, if there is no vector, then I will move the first column of the matrix, and
replace that matrix column with a constant column of ones, just like with OLS.

\adoc    settings   None, but see above about seeking a factor page in the input data.


\li PS: Here is a nice trick used in the implementation. let \f$y_i = x\beta_i\f$.
  Then
\f[ln(\sum_i{e^{x\beta_i}}) = max(y_i) + ln(\sum_i{e^{y_i - max(y_i)}}).\f]

The elements of the sum are all now exp(something negative), so 
overflow won't happen, and if there's underflow, then that term
must not have been very important. [This trick is attributed to Tom
Minka, who implemented it in his Lightspeed Matlab toolkit.]

*/
apop_model apop_logit = {.name="Logit", .log_likelihood = multilogit_log_likelihood, .dsize=-1,
.score = logit_dlog_likelihood, .predict=multilogit_expected, .prep = probit_prep
};
