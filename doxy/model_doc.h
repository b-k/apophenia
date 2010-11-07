/** \file */
/** \defgroup models */

/* This is a hack to produce a structured list of documentation items for every model.

   Right now, the priority is to get all models to have a more uniform interface,
   and uniform documentation is a first step. Thus, every model needs to address the
   headers enumerated here. At the bottom of this file, you'll find the null model,
   which you can use as a template.

   I achieve the structured documentation by lying to Doxygen and telling it that
   the model, like apop_OLS, is actually an enum of elements including overview,
   name, settings, &c. It'll then use the enum format, which has the appropriate look,
   giving the full documentation under each header.

   There are minor problems with this lie, which are solved with a simple sed script
   to post-process (in the Makefile).  The one this is most relevant is that every
   enum item is global, so each needs to have a number tagged on, lest Doxygen gets
   confused. Add the numbers here; I'll strip them in sed.
*/

/** The Linear Model family of models

\hideinitializer \ingroup models 
 
Models based on the form \f$Y = f(x_1) + f(x_2) + ... + \epsilon\f$.

There are many customs about input formats and outputs, meaning that these models do
things slightly differently. Also, they typically have only partially-specified
probability/likelihood functions; the segment on RNGs here explains how Apophenia turns these into full probability models.
 
 */
enum apop_lm_family{
    Overview32,
        /**< You will need to provide the weights in \c yourdata->weights. 
          Otherwise, this model will behave just like \ref apop_ols. */
    Name32,         /**< <tt>Weighted Least Squares</tt> */
    Data_format32,  /**< OLS-style; see above. */
    Parameter_format32, /**< As per OLS: a vector */
    Estimate_results32, /**< 
    \return Will return an \ref apop_model <tt>*</tt>.
    <tt>The_result->parameters</tt> will hold the coefficients; the first
    coefficient will be the coefficient on the constant term, and the
    remaining will correspond to the independent variables. It will therefore
    be of size <tt>(data->size2)</tt>. Do not pre-allocate.

    If you asked for them, the covariance matrix, confidence levels, and
    residuals will also be returned. The confidence intervals give the
    level of certainty with which we can reject the hypothesis that the given 
    coefficient is zero.

    Residuals: An \ref apop_data structure with three columns. If this is a model with a single dependent and lots of independent vars, then the first column is the actual data. Let our model be \f$ Y = \beta X + \epsilon\f$. Then the second column is the predicted values: \f$\beta X\f$, and the third column is the residuals: \f$\epsilon\f$. The third column is therefore always the first minus the second, and this is probably how that column was calculated internally.

    Given your estimate \c est, the zeroth element is one of <br> 
    <tt> apop_data_get(est->parameters, .row=0, .colname="observed", .page= "Predicted"),</tt><br>
    <tt> apop_data_get(est->parameters, .row=0, .colname="predicted", .page= "Predicted") or</tt><br>
    <tt> apop_data_get(est->parameters, .row=0, .colname="residual", .page= "Predicted").</tt><br>

    Also, I'll run \f$t\f$-tests on the hypothesis that each
    parameter is different from zero, by running \ref apop_estimate_parameter_tests, qv.
                        */
    Prep_routine32, /**< Focuses on the data shunting. */
    RNG32, /**< Linear models are typically only partially defined probability
        models. For OLS, we know that \f$P(Y|X\beta) \sim {\cal N}(X\beta, \sigma)\f$,
        because this is an assumption about the error process, but we don't know much of
        anything about the distribution of \f$X\f$.

     The \ref apop_lm_settings group includes an \ref apop_model* element named
     \c input_distribution. This is the distribution of the
     independent/predictor/X columns of the data set. 

     The default is that <tt>input_distribution = apop_improper_uniform </tt>, meaning
     that \f$P(X)=1\f$ for all \f$X\f$. So \f$P(Y, X) = P(Y|X)P(X) = P(Y|X)\f$. This
     seems to be how many people use linear models: the \f$X\f$ values are taken as
     certain (as with actually observed data) and the only question is the odds of
     the dependent variable. If that's what you're looking for, just leave the default.

     <em>But</em> you can't draw from an improper uniform. So if you draw from a linear
     model with a default <tt>input_distribution</tt>, then you'll get an error.

     Alternatively, you may know something about the distribution of the input data. At
     the least, you could generate a PMF from the actual data:
     \code
    apop_settings_set(your_model, apop_lm, input_distribution, apop_estimate(inset, apop_pmf));
     \endcode
     Now, random draws are taken from the input data. Or change this to any other
     appropriate distribution, such as a \ref apop_multivariate_normal, or an \ref
     apop_pmf filled in with more data, or perhaps something from
     http://en.wikipedia.org/wiki/Errors-in-variables_models , as desired.
            */
    settings32 /**< \ref apop_lm_settings */
} ;

/** Ordinary least squares

\param inset The first column is the dependent variable, and the remaining columns the independent. Is destroyed in the process, so make a copy beforehand if you need.

\param epin    An \ref apop_model object. It may have a \ref apop_lm_settings group attached.  I'll look at the \c destroy_data element; if this is NULL or \c destroy_data==0, then the entire data set is copied off, and then mangled. If \c destroy_data==1, then this doesn't copy off the data set, but destroys it in place. I also look at \c want_cov and \c want_expected_value, though I'll not produce the covariance matrix only if both are \c 'n'.


See also the page on \ref dataprep.

  \ingroup models
  \hideinitializer  */
enum apop_ols{
    Name1, /**< <tt>Ordinary Least Squares</tt> */
    Data_format1, /**< OLS-style; see \ref dataprep  */
    Parameter_format1, /**< A vector of OLS coefficients. coeff. zero
                         refers to the constant column, if any. */
    Prep_routine1, /**<     */
    Estimate_results1, /**< As per the \ref apop_lm_family.  */
    settings1, /**< \ref apop_lm_settings */
    example1,
    /**< 
First, you will need a file named <tt>data</tt> in comma-separated form. The first column is the dependent variable; the remaining columns are the independent. For example:
\verbatim
Y, X_1, X_2, X_3
2,3,4,5
1,2,9,3
4,7,9,0
2,4,8,16
1,4,2,9
9,8,7,6
\endverbatim

The program:
\include ols1.c

If you saved this code to <tt>sample.c</tt>, then you can compile it with
\verbatim
gcc sample.c -std=gnu99 -lapophenia -lgsl -lgslcblas -lsqlite3 -o run_me
\endverbatim

and then run it with <tt>./run_me</tt>. Alternatively, you may prefer to compile the program using a \ref makefile .

Feeling lazy? The program above was good form and demonstrated useful features, but the code below will do the same thing in two lines:

\code
#include <apop.h>
int main(){ apop_model_show(apop_estimate(apop_text_to_data("data"), apop_ols)); }
\endcode
*/
};

/** The Weighed Least Squares model
\hideinitializer \ingroup models */
enum apop_wls{
    Overview3, /**< This is a (deprecated) synonym for \ref apop_ols, qv.  If you use the \ref
        apop_ols model and provide weights in \c your_input_data->weights, then I will use them
        appropriately. That is, the \ref apop_ols model really implements Weighted Least Squares,
        but in most cases <tt>weights==NULL</tt> and the math reduces to the special case of 
        Ordinary Least Squares.  */
} ;

/** Instrumental variable regression

Operates much like the \ref apop_ols model, but the input parameters also need to have
a table of substitutions (like the addition of the <tt>.instruments</tt> setting in
the example below). The vector element of the table lists the column numbers to be
substituted (the dependent var is zero; first independent col is one), and then one
column for each item to substitute.

If the vector of your apop_data set is \c NULL, then I will use the row names to find
the columns to substitute. This is generally more robust and/or convenient.

If the \c instruments data set is somehow \c NULL or empty, I'll just run OLS. 
\hideinitializer \ingroup models */
enum apop_iv{
    Name4,         /**< <tt>instrumental variables</tt> */
    Data_format4,  /**< See the \ref apop_lm_family; see \ref dataprep. */
    Parameter_format4, /**< As per the \ref apop_lm_family */
    Estimate_results4, /**< As per the \ref apop_lm_family */
    Prep_routine4, /**< Focuses on the data shunting. */
    settings4, /**< \ref apop_lm_settings */
    example4,
    /**<
\code
apop_data *submatrix =apop_data_alloc(0, data->matrix->size1, 2);
APOP_COL(submatrix, 0, firstcol);
gsl_vector_memcpy(firstcol, your_data_vector);
APOP_COL(submatrix, 1, secondcol);
gsl_vector_memcpy(firstcol, your_other_data_vector);
apop_name_add(submatrix->names, "subme_1", 'r');
apop_name_add(submatrix->names, "subme_2", 'r');

Apop_model_add_group(&apop_iv, apop_lm, .instruments = submatrix);
apop_model *est = apop_estimate(data, apop_iv);
apop_model_show(est);
\endcode
    */
} ;

/** The Bernoulli model.

A single random draw with probability \f$p\f$.
  \hideinitializer \ingroup models
*/
enum apop_bernoulli {
    Name5,         /**< <tt>Bernoulli distribution</tt> */
    Data_format5,  /**< 
  The matrix or vector can have any size, and I just count up zeros
  and non-zeros. The Bernoulli parameter \f$p\f$ is the percentage of non-zero
  values in the matrix. Its variance is \f$p(1-p)\f$.
                    */
    Estimate_results5, /**< Parameters are estimated. */
    Parameter_format5, /**< A vector of length one */
    Prep_routine5, /**<  None. */
    RNG5, /**< Yes. Returns a single zero or one. */
    settings5 /**< None. */
} ;

/** Regression via loess smoothing

    This uses a somewhat black-box routine, first written by
    Chamberlain, Devlin, Grosse, and Shyu in 1988, to fit a smoothed
    series of quadratic curves to the input data, thus producing a
    curve more closely fitting than a simple regression would.

    The curve is basically impossible to describe using a short list of
    parameters, so the representation is in the form of the \c predicted
    vector of the \c expected data set; see below.

From the 1992 manual for the package:
``The method we will use to fit local regression models is called {\em
loess}, which is short for local regression, and was chosen as the name
since a loess is a deposit of fine clay or silt along a river valley,
and thus is a surface of sorts. The word comes from the German löss,
and is pronounced löíss.''
  \ingroup models \hideinitializer  */
enum apop_loess{
    Name2,         /**< <tt>loess smoothing</tt> */
    Data_format2,  /**< 
The data is basically OLS-like:                     
the first column of the data is the dependent variable to be explained;
subsequent variables are the independent explanatory variables.  Thus,
your input data can either have a dependent vector plus explanatory
matrix, or a matrix where the first column is the dependent variable.

Unlike with OLS, I won't move your original data, and I won't add a
<b>1</b>, because that's not really the loess custom. You can of course
set up your data that way if you like.

If your data set has a weights vector, I'll use it.

In any case, all data is copied into the model's \ref
apop_loess_settings. The code is primarily FORTRAN code from 1988
converted to C; the data thus has to be converted into a relatively
obsolete internal format.  
                    
*/
    Parameter_format2, /**< The parameter vector is unused. */
    Estimate_results2, /**< 
        The \ref apop_loess_settings is filled with results (and internal 
        processing cruft). The \c expected data set has the \c actual, 
        \c predicted, and \c residual columns, which is probably what you
        were looking for.*/
    Predict2, /**<
Fills in the zeroth column (ignoring and overwriting any data there), and at the data's <tt>->more</tt> pointer, adds an \ref
apop_data set named "Confidence" (i.e., 
\code
apop_strcmp(outdata->more->names->title, "Confidence") == 1.
\endcode */
    Prep_routine2, /**< None. */
    RNG2, /**< No. */
    settings2 /**< \ref apop_loess_settings */
} ;

/** The Beta distribution.

         The beta distribution has two parameters and is
         restricted between zero and one. You may also find \ref
         apop_beta_from_mean_var to be useful. 
\hideinitializer \ingroup models */
enum apop_beta{
    Name31,         /**< <tt>Beta distribution</tt> */
    Data_format31,  /**< Any arrangement of scalar values. */
    Parameter_format31, /**<  a vector, v[0]=\f$\alpha\f$; v[1]=\f$\beta\f$    */
    Estimate_results31, /**<  Parameter estimates   */
    Prep_routine31, /**<  None. */
    RNG31, /**< Yes, producing a scalar \f$\in[0,1]\f$. */
    settings31 /**<  None. */
} ;

/** 
The multi-draw generalization of the Bernoulli; the two-bin special case of the Multinomial.
\hideinitializer \ingroup models
*/
enum apop_binomial {
    Name6,         /**< <tt> Binomial distribution</tt>*/
    Data_format6,  /**< 
Zeros are failures and non-zeros successes. \f$N\f$
is the size of the matrix, vector, or both (whichever is not \c NULL).
So \f$p\f$ represents the odds of a success==1; the odds of a zero is \f$1-p\f$.

\li You may be interested in \ref apop_data_to_factors to convert real numbers or text into a
vector of categories.

\li See also \ref apop_data_rank_compress for means of dealing with one more input data format.
                    */
    Parameter_format6, /**< 
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
                        */
    Estimate_results6, /**<   Parameters are estimated. Covariance matrix is filled.    */
    Prep_routine6, /**<    None.     */
    RNG6, /**< Yes. 
            I fill an array of length \c n, with a sequence of randomly drawn ones and zeros. */
    CDF6, /**< yes  */
    settings6, /**<  */
    Example6 /**<      */
} ;

/** The \f$n\f$--option generalization of the Binomial distribution.
    See also the \ref apop_binomial model.

\hideinitializer
\ingroup models
*/
enum apop_multinomial {
    Name7,         /**< <tt> Multinomial distribution</tt>*/
    Data_format7,  /**< 
The default is simply a listing of bins, without regard to whether items are in the vector or
matrix of the \ref apop_data struct, or the dimensions. Here, data like <tt>0, 1, 2, 1, 1</tt>
represents one draw of zero, three draws of 1, and one draw of 2.

\li You may be interested in \ref apop_data_to_factors to convert real numbers or text into a
vector of categories.

\li See also \ref apop_data_rank_compress for means of dealing with one more input data format.

\li Please note that the number of bins is simply the largest number found. So if there
are bins \{0, 1, 2\} and your data set happens to consist of <tt>0 0 1 1 0</tt>, then
I won't know to generate results with three bins where the last bin has probability zero.
                    */
    Parameter_format7, /**< 
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
                        */
    Estimate_results7, /**<  Parameters are estimated. Covariance matrix
                         is filled.   */
    Prep_routine7, /**<   None.      */
    RNG7, /**< The result of an imaginary tossing of \f$N\f$ balls into \f$k\f$ urns, with the
            given probabilities.
            
            I fill an array of length \c N, with a sequence of draws from zero to \f$N\f$. They
            are not randomly ordered: it'll look something like \f$[0 0 1 1 3 3 3]\f$, but
            will still be an accurate representation of what happens when you throw
            \f$N\f$ balls into \f$k\f$ urns and sum the results. 
            
            If you want the sequence of draws to be random at the per-item scale,
            set \f$N=1\f$ (i.e., <tt>apop_data_set(estimated->parameters, 0, 1);</tt>),
            and use a \c for loop to make the number of draws you want. This is less efficient.
            */
    settings7, /**<  */
    Example7 /**<      */
} ;

/** 
A multivariate generalization of the \ref apop_beta "Beta distribution".
\hideinitializer \ingroup models
*/
enum apop_dirichlet {
    Name8,         /**< <tt>Dirichlet distribution</tt>*/
    Data_format8,  /**<     Each row of your data is a single observation.  */
    Parameter_format8, /**<  
The estimated parameters are in the output model's <tt>parameters->vector</tt>.
                        */
    Estimate_results8, /**<  Found via maximum likelihood.   */
    Prep_routine8, /**<   None.      */
    RNG8, /**< Yes. */
    settings8, /**<  MLE-type: \ref apop_mle_settings, \ref apop_parts_wanted_settings   */
    Example8 /**<      */
} ;

/** 
A one-parameter likelihood function.

\f$Z(\mu,k) 	= \sum_k 1/\mu e^{-k/\mu} 			\f$ <br>
\f$ln Z(\mu,k) 	= \sum_k -\ln(\mu) - k/\mu			\f$ <br>
\f$dln Z(\mu,k)/d\mu 	= \sum_k -1/\mu + k/(\mu^2)			\f$ <br>

Some write the function as:
\f$Z(C,k) dx = \ln C C^{-k}. \f$
If you prefer this form, just convert your parameter via \f$\mu = {1\over
\ln C}\f$ (and convert back from the parameters this function gives you
via \f$C=\exp(1/\mu)\f$.

\li See also \ref apop_data_rank_compress for means of dealing with one more input data format.
\hideinitializer \ingroup models */
enum apop_exponential {
    Name9,         /**< <tt>Exponential distribution</tt>*/
    Data_format9,  /**< 
Ignores the matrix structure of the input data, so send in a 1 x N, an N x 1, or an N x M.
                    */
    Parameter_format9, /**<  \f$\mu\f$ is in the zeroth element of the vector.   */
    Estimate_results9, /**<  Parameter is set.   */
    Prep_routine9, /**<  None.   */
    RNG9, /**< Yes. Produces a single number.*/
    CDF9, /**< Yes. Produces a single number.*/
    settings9, /**<  None. */
    Example9 /**<   */
} ;
/** The Gamma distribution

\f$G(x, a, b)     = 1/(\Gamma(a) b^a)  x^{a-1} e^{-x/b}\f$

\f$ln G(x, a, b)= -ln \Gamma(a) - a ln b + (a-1)ln(x) + -x/b\f$

\f$d ln G/ da    =  -\psi(a) - ln b + ln(x) \f$    (also, \f$d ln \gamma = \psi\f$)

\f$d ln G/ db    =  -a/b + x/(b^2) \f$

apop_gamma.estimate() is an MLE, so feed it appropriate \ref apop_mle_settings.
  
\li See also \ref apop_data_rank_compress for means of dealing with one more input data format.
\hideinitializer \ingroup models */
enum apop_gamma{
    Name10,         /**< <tt>Gamma distribution</tt>*/
    Data_format10,  /**<    
Location of data in the grid is not relevant; send it a 1 x N, N x 1, or N x M and it will all be the same.     */
    Parameter_format10, /**<  First two elements of the vector.   */
    Estimate_results10, /**<  Parameters are estimated, using MLE.   */
    Prep_routine10, /**<    None.     */
    RNG10, /**< Yes. */
    settings10, /**<   MLE-type: \ref apop_mle_settings, \ref apop_parts_wanted_settings  */
    Example10 /**<      */
} ;

/** The kernel density smoothing of a PMF

A Kernel density is simply a smoothing of a histogram. At each point
along the histogram, put a distribution (default: Normal(0,1)) on top
of the point. Sum all of these distributions to form the output histogram.

Elements of \ref apop_kernel_density_settings that you may want to set:

\li data    a data set, which, if  not \c NULL and \c !histobase , will be converted to a histogram.
\li histobase This is the preferred format for input data. It is the histogram to be smoothed.
\li kernelbase The kernel to use for smoothing, with all parameters set and a \c p method. Popular favorites are \ref apop_normal and \ref apop_uniform.
\li set_params A function that takes in a single number and the model, and sets
the parameters accordingly. The function will call this for every point in the data
set. Here is the default, which is used if this is \c NULL. It simply sets the first
element of the model's parameter vector to the input number; this is appropriate for a
Normal distribution, where we want to center the distribution on each data point in turn.

\code
void apop_set_first_param(double in, apop_model *m){
    m->parameters->vector->data[0] = in;
}
\endcode

For a Uniform[0,1] recentered around each point, you'd want to put this function in your code:

\code
void set_midpoint(double in, apop_model *m){
    m->parameters->vector->data[0] = in-0.5;
    m->parameters->vector->data[1] = in+0.5;
}
\endcode

\hideinitializer \ingroup models */
enum apop_kernel_density{
    Name30,         /**< <tt>kernel density estimate,</tt>*/
    Data_format30,  /**< I'll estimate a \ref apop_pmf internally, so I
                      follow that format, which is one observation (of any format) per line.        */
    Parameter_format30, /**< None    */
    Estimate_results30, /**< the estimate method basically just runs
                          <tt>apop_model_add_group(your_data, apop_kernel_density);</tt>    */
    Prep_routine30, /**<  None   */
    RNG30, /**< Uses the default. */
    settings30, /**< \ref apop_kernel_density_settings.    */
    Example30 /**<    
This example sets up and uses KDEs based on a Normal and a Uniform distribution.

\include kernel.c
               */
} ;

/** A one-dimensional histogram.

  This is an empirical distribution. If you have a data set from which you want to make random draws, this is overkill; instead just use something like \code 
  gsl_rng *r = apop_rng_alloc(27);
  gsl_vector *my_data = [gather data here.];
  gsl_vector_get(my_data, gsl_rng_uniform(r)*my_data->size);
  \endcode

  But this can be used anywhere a model is needed, such as the inputs and outputs to \c apop_update.

  The model is unlike most other models in that there are no parameters
  of any sort (beyond the data itself), so 
    all the work of producing the histogram is done in \c apop_histogram_settings_alloc.

    See also the \c apop_pmf function, which implements something very
    similar, and works for multidimensional data.

\hideinitializer \ingroup models */
enum apop_histogram{
    Name11,         /**< <tt>Histogram</tt>*/
    Data_format11,  /**< A free-format pile of scalars, in vector, matrix, or both.        */
    Parameter_format11, /**<  No parameters; all information is kept in
                         the \ref apop_histogram_settings struct.   */
    Estimate_results11, /**<  Just an alias for <tt>Apop_settings_add_group(out, apop_histogram, d, 1000);</tt>. Once this is estimated, \c p, \c log_likelihood, \c draw, &c. work.
                             */
    Prep_routine11, /**<    None.     */
    RNG11, /**< Yes. The first call produces a cumulative density tally,
            and so will take several microseconds longer than later calls.  */
    settings11, /**<   \ref apop_histogram_settings   */
    Example11 /**<      */
};

/** The improper uniform returns \f$P(x) = 1\f$ for every value of x, all the
time (and thus, log likelihood(x)=0).  It has zero parameters. It is
useful, for example, as an input to Bayesian updating, to represent a
fully neutral prior.
\hideinitializer \ingroup models */
enum apop_improper_uniform {
    Name13,         /**< <tt>Improper uniform distribution</tt>*/
    Data_format13,  /**<      None   */
    Parameter_format13, /**<  \c NULL   */
    Estimate_results13, /**<   The \c estimate routine is just a dummy that returns its input.  */
    Prep_routine13, /**<    None.     */
    RNG13, /**< The \c draw function makes no sense, and therefore returns an error. */
    CDF13, /**< Half of the distribution is less than every given point, so the CDF always
             returns 0.5. One could perhaps make an argument that this should really be
             infinity, but a half is more in the spirit of the distribution's
             use to represent a lack of information. */
    settings13, /**< None.     */
    Example13 /**<      */
} ;

/** This is the two-parameter version of the Uniform, expressing a uniform distribution over [a, b].

The MLE of this distribution is simply a = min(your data); b = max(your data).
Primarily useful for the RNG, such as when you have a Uniform prior model.
\hideinitializer \ingroup models */
enum apop_uniform {
    Name14,         /**< <tt>Uniform distribution</tt>*/
    Data_format14,  /**<  An unordered set of numbers in the data set's
                      vector, matrix, or both       */
    Parameter_format14, /**<  Zeroth vector element is \f$a\f$, the min;
                          first is \f$b\f$, the max.   */
    Estimate_results14, /**<  Parameters are set.   */
    Prep_routine14, /**<  None.      */
    RNG14, /**< Yes. */
    settings14, /**<  None.    */
    Example14 /**<      */
} ;

/** 
  
  For the Poisson distribution: 
\f$p(k) = {\mu^k \over k!} \exp(-\mu), \f$
\hideinitializer \ingroup models */
enum apop_poisson {
    Name15,         /**< <tt>poisson</tt>*/
    Data_format15,  /**< 
Location of data in the grid is not relevant; send it a 1 x N, N x 1, or N x M and it will all be the same.
                     */
    Parameter_format15, /**< One parameter, the zeroth element of the vector.    */
    Estimate_results15, /**<  Parameters are set. 

  Model's \c llikelihood element is calculated.

  Unless you set
  \code
  Apop_model_add_group(your_model, apop_lm, .want_cov='n')
  \endcode
  I will also give you the variance of the parameter, via jackknife.
                         */
    Prep_routine15, /**<   None.      */
    RNG15, /**< Yes. */
    settings15, /**<  \ref apop_lm_settings, for the \c .want_cov element    */
    Example15 /**<      */
} ;

/** This is the multivariate generalization of the Normal distribution.
\hideinitializer \ingroup models */
enum apop_multivariate_normal{
    Name16,         /**< <tt>Multivariate normal distribution</tt>*/
    Data_format16,  /**<    Each row of the matrix is an observation.     */
    Parameter_format16, /**< 
  An \c apop_data set whose vector element is the vector of
  means, and whose matrix is the covariances.

  Likelihoods are set.  */
    Estimate_results16, /**< Parameters are set in the above format. 
                      <tt> apop_data_get_page(est->parameters, "cov");</tt> gives the
                      covariance matrix of the means.   */
    Prep_routine16, /**<  None.       */
    RNG16, /**< The RNG fills an input array whose length is based on the input parameters. */
    settings16, /**<  None.    */
    Example16 /**<      */
} ;

/** The Waring distribution

\f$W(x,k, b,a) 	= (b-1) \gamma(b+a) \gamma(k+a) / [\gamma(a+1) \gamma(k+a+b)]\f$

\f$\ln W(x, b, a) = \ln(b-1) + \ln\gamma(b+a) + \ln\gamma(k+a) - \ln\gamma(a+1) - \ln\gamma(k+a+b)\f$

\f$dlnW/db	= 1/(b-1)  + \psi(b+a) - \psi(k+a+b)\f$

\f$dlnW/da	= \psi(b+a) + \psi(k+a) - \psi(a+1) - \psi(k+a+b)\f$


See also \ref apop_data_rank_compress for means of dealing with one more input data format.

\hideinitializer \ingroup models */
enum apop_waring {
    Name17,         /**< <tt>Waring distribution</tt>*/
    Data_format17,  /**<    
Ignores the matrix structure of the input data, so send in a 1 x N, an N x 1, or an N x M.
     */
    Parameter_format17, /**< Two elements in the parameter set's vector.    */
    Estimate_results17, /**< Estimated via MLE.    */
    Prep_routine17, /**<  None.       */
    RNG17, /**< Yes. */
    settings17, /**<  MLE-type: \ref apop_mle_settings, \ref apop_parts_wanted_settings    */
    Example17 /**<      */
} ;

/** 
The special case of the \ref apop_waring "Waring" where \f$ \alpha = 0.	\f$<br>

\f$ Y(x, b) 	= (b-1) \gamma(b) \gamma(k) / \gamma(k+b)			\f$

\f$ \ln Y(x, b)	= \ln(b-1) + ln\gamma(b) + \ln\gamma(k) - \ln\gamma(k+b)	\f$

\f$ d\ln Y/db	= 1/(b-1)  + \psi(b) - \psi(k+b)				\f$

apop_yule.estimate() is an MLE, so feed it appropriate \ref apop_mle_settings.

See also \ref apop_data_rank_compress for means of dealing with one more input data format.
\hideinitializer \ingroup models */
enum apop_yule {
    Name18,         /**< <tt>Yule distribution</tt>*/
    Data_format18,  /**<    
Ignores the matrix structure of the input data, so send in a 1 x N, an N x 1, or an N x M.*/
    Parameter_format18, /**< One element at the top of the parameter set's vector.*/
    Estimate_results18, /**< Estimated via MLE.    */
    Prep_routine18, /**<  None.       */
    RNG18, /**< Yes. */
    settings18, /**<  MLE-type: \ref apop_mle_settings, \ref apop_parts_wanted_settings    */
    Example18 /**<      */
} ;

/** 
Wikipedia has notes on the <a href="http://en.wikipedia.org/wiki/Zipf_distribution">Zipf distribution</a>. 
\f$Z(a)        = {1\over \zeta(a) * i^a}        \f$

\f$lnZ(a)    = -(\log(\zeta(a)) + a \log(i))    \f$

apop_zipf.estimate() is an MLE, so feed it appropriate \ref apop_mle_settings.

See also \ref apop_data_rank_compress for means of dealing with one more input data format.
\hideinitializer \ingroup models */
enum apop_zipf {
    Name19,         /**< <tt>Zipf distribution</tt>*/
    Data_format19,  /**<    Ignores the matrix structure of the input data, so send in a 1 x N, an N x 1, or an N x M.*/
    Parameter_format19, /**< One item at the top of the parameter set's vector.    */
    Estimate_results19, /**< Estimates the parameter.    */
    Prep_routine19, /**<  None.       */
    RNG19, /**< Yes. */
    settings19, /**<  \ref apop_mle_settings, \ref apop_parts_wanted_settings    */
    Example19 /**<      */
} ;


/** 
A probability mass function is commonly known as a histogram, or still more commonly,
a bar chart. It indicates that at a given coordinate, there is a given mass.

The data format for the PMF is simple: each row holds the coordinates, and the
<em>weights vector</em> holds the mass at the given point. This is in contrast to the
crosstab format, where the location is simply given by the position of the data point
in the grid.

For example, here is a typical crosstab:

<table>
<tr>            <td></td><td> col 0</td><td> col 1</td><td> col 2</td></tr>
<tr><td>row 0 </td><td>0</td><td> 8.1</td><td> 3.2</td></tr>
<tr><td>row 1 </td><td>0</td><td> 0</td><td> 2.2</td></tr>
<tr><td>row 2 </td><td>0</td><td> 7.3</td><td> 1.2</td></tr>
</table>

Here it is as a sparse listing:

<table>
<tr>        <td></td> value</td></td> dimension 1<td></td> dimension 2<td></tr>
<tr> <td>8.1</td> <td>0</td> <td>1</td> </tr>
<tr> <td>3.2</td> <td>0</td> <td>2</td> </tr>
<tr> <td>2.2</td> <td>1</td> <td>2</td> </tr>
<tr> <td>7.3</td> <td>2</td> <td>1</td> </tr>
<tr> <td>1.2</td> <td>2</td> <td>2</td> </tr>
</table>

The \c apop_pmf internally represents data in this manner. The dimensions are held
in the \c matrix element of the data set, and the cell values are held in the \c weights
element (<em>not the vector</em>).

If your data is in a crosstab (with entries in the matrix element for 2-D data or the
vector for 1-D data), then use \ref apop_crosstab_to_pmf to make the conversion.

If your data is already in the sparse listing format (which is probably the case for 3-
or more dimensional data), then just point the model to your parameter set:

\code
apop_model *my_pmf = apop_model_copy(apop_pmf);
my_pmf->parameters = in_data;
//or equivalently:
apop_model *my_pmf = apop_estimate(in_data, apop_pmf);
\endcode

\li If the \c weights element is \c NULL, then I assume that all rows of the data set are
equally probable.

\li Be careful: the weights are in the \c weights element of the \c apop_data set, not in
the \c vector element. If you put the weights in the \c vector and have \c NULL \c
weights, then draws are equiprobable. This will be difficult to debug.

\li This format is ideal for holding sparse multidimensional matrices, drawing from them,
or calculating likelihoods. However, it is not ideal for linear algebra operations
on sparse matrices. For this, the column-packed format is recommended. Integration
with the CSparse library for linear algebra on sparse matrices is forthcoming.

\todo Doing MLE on this model doesn't make sense. Don't do it.
\hideinitializer \ingroup models */ 
enum apop_pmf {
    Name20,         /**< <tt>PDF or sparse matrix</tt>*/
    Data_format20,  /**<    As above, you can input to the \c estimate
                      routine a 2-D matrix that will be converted into this form.     */
    Parameter_format20, /**< The \c parameters element is the weights matrix as above.    */
    Estimate_results20, /**< Produces the parameter array.    */
    Prep_routine20, /**<  None.       */
    RNG20, /**< Yes. The first time you draw from a PMF, I will generate a CMF
(Cumulative Mass Function). For an especially large data set this
may take a human-noticeable amount of time. The CMF will be stored
in <tt>parameters->weights[1]</tt>, and subsequent draws will have no
computational overhead. */
    settings20, /**<  None.    */
    Example20 /**<      */
} ;

/** The Probit model.

  Apophenia makes no distinction between the bivariate probit and the multinomial probit.  Just use this one.
\hideinitializer \ingroup models */
enum apop_probit {
    Name21,         /**< <tt>Probit</tt>*/
    Data_format21,  /**< 
 The first column of the data matrix this model expects is zeros, ones, ..., enumerating
 the factors; to get there, try \ref apop_data_to_factors; if you  forget to run it, I'll
 run it on the first data column for you.
 The remaining columns are values of the independent variables. Thus,
 the model will return [(data columns)-1]\f$\times\f$[(option count)-1] parameters. 
 Column names are options; row names are input variables.
      */
    Parameter_format21, /**< As above */
    Estimate_results21, /**< Via MLE.    */
    Prep_routine21, 
/**< You will probably want to convert some column of your data into factors, via
    \ref apop_data_to_factors. If you do, then that adds a page of factors to your data
    set (and of course adjusts the data itself). If I find a factor page, I will use that
    info; if not, then I will run \ref apop_data_to_factors on the first column (the
    vector if there is one, else the first column of the matrix.)  

    Also, if there is no vector, then I will move the first column of the matrix, and
    replace that matrix column with a constant column of ones, just like with OLS.
 */
    RNG21, /**< No. */
    settings21, /**<  None, but see above about seeking a factor page in the input data.*/
    Example21 /**<      */
} ;

/** The Logit model.

  Apophenia makes no distinction between the bivariate logit and the multinomial logit.  Just use this one.

  The likelihood of choosing item \f$j\f$ is:
  \f$e^{x\beta_j}/ (\sum_i{e^{x\beta_i}})\f$

  so the log likelihood is 
  \f$x\beta_j  - ln(\sum_i{e^{x\beta_i}})\f$

  A nice trick used in the implementation: let \f$y_i = x\beta_i\f$.
  Then
\f[ln(\sum_i{e^{x\beta_i}}) = max(y_i) + ln(\sum_i{e^{y_i - max(y_i)}}).\f]

The elements of the sum are all now exp(something negative), so 
overflow won't happen, and if there's underflow, then that term
must not have been very important. [This trick is attributed to Tom
Minka, who implemented it in his Lightspeed Matlab toolkit.]
\hideinitializer \ingroup models */
enum apop_logit {
    Name22,         /**< <tt>Logit</tt>*/
    Data_format22,  /**< 
 The first column of the data matrix this model expects is zeros, ones, ..., enumerating
 the factors; to get there, try \ref apop_data_to_factors; if you  forget to run it, I'll
 run it on the first data column for you.
 The remaining columns are values of the independent variables. Thus,
 the model will return [(data columns)-1]\f$\times\f$[(option count)-1] parameters. 
 Column names are options; row names are input variables.
      */
    Parameter_format22, /**< As above.    */
    Estimate_results22, /**< Via MLE.    */
    Prep_routine22, 
/**< You will probably want to convert some column of your data into factors, via
    \ref apop_data_to_factors. If you do, then that adds a page of factors to your data
    set (and of course adjusts the data itself). If I find a factor page, I will use that
    info; if not, then I will run \ref apop_data_to_factors on the first column (the
    vector if there is one, else the first column of the matrix.)  

    Also, if there is no vector, then I will move the first column of the matrix, and
    replace that matrix column with a constant column of ones, just like with OLS.
 */
    RNG22, /**< No. */
    settings22, /**<  None, but see above about seeking a factor page in the input data.*/
    Example22 /**<      */
} ;

/** The Multinomial Probit model.
  \deprecated Just use \ref apop_probit, which handles multiple options fine.
\hideinitializer \ingroup models */
enum apop_multinomial_probit {
    Name23,         /**< <tt>Probit</tt>*/
} ;

/** The Wishart distribution, which is currently somewhat untested. 

Here's the likelihood function. \f$p\f$ is the dimension of the data and covariance matrix,
\f$n\f$ is the degrees of freedom, \f$\mathbf{V}\f$ is the \f$p\times
p\f$ matrix of Wishart parameters, and \f${\mathbf{W}}\f$ is the \f$p\times
p\f$ matrix whose likelihood is being evaluated.  \f$\Gamma_p(\cdot)\f$
is the \ref apop_multivariate_gamma "multivariate gamma function".

\f$
P(\mathbf{W}) = \frac{\left|\mathbf{W}\right|^\frac{n-p-1}{2}}
                         {2^\frac{np}{2}\left|{\mathbf V}\right|^\frac{n}{2}\Gamma_p(\frac{n}{2})} \exp\left(-\frac{1}{2}{\rm Tr}({\mathbf V}^{-1}\mathbf{W})\right)\f$

See also notes in \ref tfchi.
\hideinitializer \ingroup models */
enum apop_wishart  {
    Name24,         /**< <tt>Wishart distribution</tt>*/
    Data_format24,  /**<    Each row of the input matrix is a single square matrix,
                      flattened; use \ref apop_data_pack to convert your
                      sequence of matrices into rows.     */
    Parameter_format24, /**< \f$N\f$ (the degrees of freedom) is the zeroth element of the vector. The matrix holds the matrix of parameters.*/
    Estimate_results24, /**< Via MLE.    */
    Prep_routine24, /**<  Just allocates the parameters based on the size of the input data.       */
    RNG24, /**< Yes. You can use this to generate random covariance matrices, should you need them. See example below. */
    settings24, /**<  \ref apop_mle_settings, \ref apop_parts_wanted_settings    */
    Example24 /**<  
Making some random draws:

\code
gsl_matrix *rmatrix = gsl_matrix_alloc(10, 10);
gsl_rng *r = apop_rng_alloc(8765);
for (int i=0; i< 1e8; i++){
    apop_draw(rmatrix->data, r, apop_wishart);
    do_math_with_matrix(rmatrix);
}
\endcode    */
} ;

/** The t distribution, primarily for descriptive purposes.

 If you want to test a hypothesis, you probably don't need this, and should instead use \ref apop_test.  See notes in \ref tfchi.  
\hideinitializer \ingroup models */
enum apop_t_distribution  {
    Name25,         /**< <tt>t distribution</tt>*/
    Data_format25,  /**<    Unordered list of scalars in the matrix and/or vector.     */
    Parameter_format25, /**< vector->data[0] = mu<br>
                            vector->data[1] = sigma<br>
                            vector->data[2] = df */
    Estimate_results25, /**< If you do not set an \ref apop_mle_settings
                          group beforehand, I'll just count elements and
                          set \f$df = n-1\f$. Else, via MLE.    */
    Prep_routine25, /**<  None.       */
    RNG25, /**< Yes. */
    settings25, /**<  \ref apop_mle_settings, \ref apop_parts_wanted_settings   */
    Example25 /**<   */
} ;

/** The F distribution, for descriptive purposes.

 If you want to test a hypothesis, you probably don't need this, and should instead use \ref apop_test.  See notes in \ref tfchi.  
\hideinitializer \ingroup models */
enum apop_f_distribution  {
    Name26,         /**< <tt>F distribution</tt>*/
    Data_format26,  /**<    Unordered list of scalars in the matrix and/or vector.     */
    Parameter_format26, /**< Zeroth and first elements of the vector are the \f$df\f$s. */
    Estimate_results26, /**< If you do not set an \ref apop_mle_settings
                          group beforehand, I'll just count elements and
                          set \f$df=\f$ vector count minus one, and 
                          \f$df2=\f$ matrix count minus one. Else, via MLE.    */
    Prep_routine26, /**<  None.       */
    RNG26, /**< Yes. */
    settings26, /**<  \ref apop_mle_settings    */
    Example26 /**<   */
} ;

/** The \f$\chi^2\f$ distribution, for descriptive purposes.

 If you want to test a hypothesis, you probably don't need this, and should instead use \ref apop_test.  See notes in \ref tfchi.  
\hideinitializer \ingroup models */
enum apop_chi_squared  {
    Name27,         /**< <tt>Chi squared distribution</tt>*/
    Data_format27,  /**<    Unordered list of scalars in the matrix and/or vector.     */
    Parameter_format27, /**< Zeroth element of the vector is the \f$df\f$. */
    Estimate_results27, /**< If you do not set an \ref apop_mle_settings
                          group beforehand, I'll just count elements and
                          set \f$df = n-1\f$. Else, via MLE.    */
    Prep_routine27, /**<  None.       */
    RNG27, /**< Yes. */
    settings27, /**<  \ref apop_mle_settings    */
    Example27 /**<   */
} ;

/** The Lognormal distribution

The log likelihood function for the lognormal distribution:

\f$f = exp(-(ln(x)-\mu)^2/(2\sigma^2))/ (x\sigma\sqrt{2\pi})\f$
\f$ln f = -(ln(x)-\mu)^2/(2\sigma^2) - ln(x) - ln(\sigma\sqrt{2\pi})\f$
\hideinitializer \ingroup models */
enum apop_lognormal {
    Name28,         /**< <tt>Lognormal distribution</tt> */
    Data_format28,  /**<    I use the elements of the matrix, without regard to their order. */
    Parameter_format28, /**< Zeroth vector element is the mean (after
                         logging); first is the std dev (after logging)    */
    Estimate_results28, /**< Parameters are set. Log likelihood is calculated.    */
    Prep_routine28, /**<  None.       */
    RNG28, /**< Yes. */
    settings28, /**<  None.    */
    Example28 /**<   */
} ;

/** The Normal (Gaussian) distribution

You know it, it's your attractor in the limit, it's the Gaussian distribution.

\f$N(\mu,\sigma^2) = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2)\f$

\f$\ln N(\mu,\sigma^2) = (-(x-\mu)^2 / 2\sigma^2) - \ln (2 \pi \sigma^2)/2 \f$

\f$d\ln N(\mu,\sigma^2)/d\mu = (x-\mu) / \sigma^2 \f$

\f$d\ln N(\mu,\sigma^2)/d\sigma^2 = ((x-\mu)^2 / 2(\sigma^2)^2) - 1/2\sigma^2 \f$

\hideinitializer \ingroup models */
enum apop_normal{
    Name29,         /**< <tt>Normal distribution</tt>*/
    Data_format29,  /**<    I use the elements of the matrix, without regard to their order. */
    Parameter_format29, /**< 
  As is custom, the first parameter (in the vector) is the mean, the second is the standard deviation (i.e., the square root of the variance). */
    Estimate_results29, /**< Parameters are set. Log likelihood is calculated. Covariance of the parameters is calculated unless <tt>.want_cov='n'</tt>; see below.    */
    Predict29, /**< Returns the expected value. The <tt>->more</tt>
                 element holds a \ref apop_data set with the title <tt>"Covariance"</tt>, whose matrix holds the covariance of the mean. If <tt>estimated_model->data</tt> then this is going to be \f$\mu/\sqrt{n}\f$; if the model's <tt>data == NULL</tt> then cov = 0. 
                 
                 Format subject to change. */
    Prep_routine29, /**<  None.       */
    RNG29, /**< Of course. */
    settings29, /**<  
The \c apop_lm_settings group includes a \c want_cov element, which
refers to the covariance matrix for the mean and variance. You can set
that to \c 'n' if you are estimating millions of these and need to save
time (i.e. <tt>Apop_model_add_group(your_model, apop_lm, .want_cov =
'n');</tt>).  */
    Example29 /**<   */
} ;

/** Sometimes ya need a null model. The name is blank, so <tt>strlen(m->name)</tt> will tell you
    whether this is the null model.
\hideinitializer \ingroup models */
enum apop_null{
    Name0,         /**< <tt>""</tt>, an empty string (which is distinct from \c NULL.)*/
    Data_format0,  /**<    None.     */
    Parameter_format0, /**< None.    */
    Estimate_results0, /**< None.    */
    Prep_routine0, /**<  None.       */
    RNG0, /**< No. */
    settings0, /**<  None.    */
    Example0 /**<   */
} ;
