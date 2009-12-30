/** \file */
/** \defgroup models */

/* This is a hack to produce a structured list of documentation items for every model.

   The goal is to produce uniform documentation for all models. Right
   now, the priority is to get all models to have a more uniform
   interface, and uniform documentation is a first step. Thus, the
   headers that every model needs to address. At the bottom of this
   file, you'll find the null model, which you can use as a template.
   
   I achieve this by lying to Doxygen and telling it that the model,
   like apop_OLS, is actually an enum of elements including 
   overview, name, key_settings, &c. It'll then use the enum format,
   which has the appropriate look, giving the full documentation under
   each header.

   There are minor problems with this lie, which are solved with a simple
   sed script to post-process.  The one this is most relevant is that
   every enum item is global (which is why I (BK) avoid them), so each
   needs to have a number tagged on, lest Doxygen gets confused. Add the
   numbers here; I'll strip them in sed.
*/

/** Ordinary least squares
  \ingroup models
  \hideinitializer  */
enum apop_ols{
    Overview1 
        /**< 
\param inset The first column is the dependent variable, and the remaining columns the independent. Is destroyed in the process, so make a copy beforehand if you need.

\param epin    An \ref apop_model object. It may have a \ref apop_ls_settings group attached.  I'll look at the \c destroy_data element; if this is NULL or \c destroy_data==0, then the entire data set is copied off, and then mangled. If \c destroy_data==1, then this doesn't copy off the data set, but destroys it in place. I also look at \c want_cov and \c want_expected_value, though I'll not produce the covariance matrix only if both are \c 'n'.


See also the page on \ref dataprep.


\hideinitializer
 */
    Name1, /**< <tt>"Ordinary Least Squares"</tt> */
    Data_format1, /**< OLS-style; see \ref dataprep  */
    Parameter_format1, /**< A vector of OLS coefficients. coeff. zero
                         refers to the constant column, if any. */
    Prep_routine1, /**< First, data is shunted. */
    Estimate_results1, /**<
        \return Will return an \ref apop_model <tt>*</tt>.
        <tt>The_result->parameters</tt> will hold the coefficients; the first
        coefficient will be the coefficient on the constant term, and the
        remaining will correspond to the independent variables. It will therefore
        be of size <tt>(data->size2)</tt>. Do not pre-allocate.

        If you asked for them, the covariance matrix, confidence levels, and
        residuals will also be returned. The confidence intervals give the
        level of certainty with which we can reject the hypothesis that the
        given coefficient is zero.
          */
    key_settings1, /**< \ref apop_ls_settings */
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
gcc sample.c -lapophenia -lgsl -lgslcblas -lsqlite3 -o run_me
\endverbatim

and then run it with <tt>./run_me</tt>. Alternatively, you may prefer to compile the program using a \ref makefile .

Feeling lazy? The program above was good form and demonstrated useful features, but the code below will do the same thing in two lines:

\code
#include <apop.h>
int main(){ apop_model_show(apop_estimate(apop_text_to_data("data"), apop_ols)); }
\endcode
*/
};


/** The WLS model
\hideinitializer
\ingroup models
*/
enum apop_wls{
    Overview3,
        /**< You will need to provide the weights in \c yourdata->weights. 
          Otherwise, this model will behave just like \ref apop_ols. */
    Name3,         /**< <tt>"Weighted Least Squares"</tt> */
    Data_format3,  /**< OLS-style; see above. */
    Parameter_format3, /**< As per OLS: a vector */
    Estimate_results3, /**< */
    Prep_routine3, /**< Focuses on the data shunting. */
    key_settings3 /**< \ref apop_ls_settings */
} ;

/** Instrumental variable regression
\ingroup models
\hideinitializer */
enum apop_iv{
    Overview4,
        /**< 
     Operates much like the \ref apop_ols model, but the input
     parameters also need to have a table of substitutions (like the
     addition of the \c .instruments setting in the example below). The vector
     element of the table lists the column numbers to be substituted (the
     dependent var is zero; first independent col is one), and then one
     column for each item to substitute.

    If the vector of your apop_data set is NULL, then I will use the row
    names to find the columns to substitute. This is generally more robust
    and/or convenient.

    If the \c instruments data set is somehow NULL or empty, I'll just run OLS.  */
    Name4,         /**< <tt>"instrumental variables"</tt> */
    Data_format4,  /**< OLS-style; see \ref dataprep. */
    Parameter_format4, /**< As per OLS: a vector */
    Prep_routine4, /**< Focuses on the data shunting. */
    Estimate_results4, /**< */
    key_settings4, /**< \ref apop_ls_settings */
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

Apop_model_add_group(&apop_iv, apop_ls, .instruments = submatrix);
apop_model *est = apop_estimate(data, apop_iv);
apop_model_show(est);
\endcode
    */
} ;

/** The Bernoulli model.
  \hideinitializer
\ingroup models
*/
enum apop_bernoulli {
    Overview5,     /**< The Bernoulli model, of a single random draw
                     with probability \f$p\f$.
                    */
    Name5,         /**< <tt>"Bernoulli distribution"</tt> */
    Data_format5,  /**< 
  The matrix or vector can have any size, and I just count up zeros
  and non-zeros. The bernoulli paramter \f$p\f$ is the percentage of non-zero
  values in the matrix. Its variance is \f$p(1-p)\f$.
                    */
    Estimate_results5, /**< Parameters are estimated. */
    Parameter_format5, /**< A vector of length one */
    Prep_routine5, /**< First, data is shunted. */
    RNG5, /**< Yes. Returns a single zero or one. */
    key_settings5 /**< None. */
} ;

/** 
Regression via lowess smoothing
  \ingroup models
  \hideinitializer  
*/
enum apop_loess{
    Overview2,         /**< 
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
and thus is a surface of sorts. The word comes from the German l\"oss,
and is pronounced l\"o\'\i{}ss.''

                        */
    Name2,         /**< <tt>lowess smoothing</tt> */
    Data_format2,  /**< 
The data is basically OLS-like:                     
the first column of the data is the dependent variable to be explained;
subsequent varaibles are the independent explanatory variables.  Thus,
your input data can either have a dependent vector plus explanatory
matrix, or a matrix where the first column is the dependent variable.

Unlike with OLS, I won't move your original data, and I won't add a {\bf
1}, because that's not really the loess custom. You can of course set
up your data that way if you like.

If your data set has a weights vector, I'll use it.

In any case, all data is copied into the model's \ref
apop_loess_settings. The code is primarily FORTRAN code from 1988
converted to C; the data thus has to be converted into a relatively
obsolete internal format.  */
    Parameter_format2, /**< The parameter vector is unused. */
    Estimate_results2, /**< 
        The \ref apop_loess_settings is filled with results (and internal 
        processing cruft). The \c expected data set has the \c actual, 
        \c predicted, and \c residual columns, which is probably what you
        were looking for.*/
    Prep_routine2, /**< First, data is shunted. */
    RNG2, /**< No. */
    key_settings2 /**< \ref apop_loess_settings */
} ;

/** The Beta distribution.
\hideinitializer
\ingroup models */
enum apop_beta{
    Overview6,         /**< 
         The beta distribution, which has two parameters and is
         restricted between zero and one. You may also find \ref
         apop_beta_from_mean_var to be useful.  */
    Name6,         /**< <tt>"Beta distribution"</tt> */
    Data_format6,  /**< Any arrangement of scalar values. */
    Parameter_format6, /**<  a vector, v[0]=\f$\alpha\f$; v[1]=\f$\beta\f$    */
    Estimate_results6, /**<  Parameter estimates   */
    Prep_routine6, /**<  None. */
    RNG6, /**< Yes, producing a scalar \f$\in[0,1]\f$. */
    key_settings6 /**<  None. */
} ;

/** The Binomial model.
\hideinitializer
\ingroup models
*/
enum apop_binomial {
    Overview6,         /**< The multi-draw generalization of the Bernoulli.
                        */
    Name6,         /**< <tt> "Binomial distribution"</tt>*/
    Data_format6,  /**< 
The default is to take the data to have a binary form, meaning that
the system counts zeros as failures and non-zeros as successes. \f$N\f$
is the size of the matrix, vector, or both (whichever is not \c NULL).

In rank-type format, the data is taken to be the two-column miss-hit
format: a nonzero value in column zero of the matrix represents a failure
and a nonzero value in column one represents successes. Set this using,
e.g., 
\code
apop_model *estimate_me = apop_model_copy(apop_binomial);
Apop_model_add_group(estimate_me, apop_rank);
apop_model *estimated = apop_estimate(your_data, estimate_me);
\endcode

In both cases, \f$p\f$ represents the odds of a success==1; the odds of a zero is \f$1-p\f$.
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
    RNG6, /**< Yes. */
    Key_settings6, /**<  \ref apop_rank_settings    */
    Example6 /**<      */
} ;

/** The Multinomial model.
    See also the \ref apop_binomial model.
\hideinitializer
\ingroup models
*/
enum apop_multinomial {
    Overview7,         /**< 
                         The \f$n\f$ option generalization of the Binomial distribution.
                        */
    Name7,         /**< <tt> "Multinomial distribution"</tt>*/
    Data_format7,  /**< 
The default is simply a listing of bins, without regard to whether items are in the vector or
matrix of the \ref apop_data struct, or the dimensions. Here, data like <tt>0, 1, 2, 1, 1</tt>
represents one draw of zero, three draws of 1, and one draw of 2.

In rank-type format, the bins are defined by the columns: a nonzero value in column zero of the
matrix represents a draw of zero, a nonzero value in column seven a draw of seven, et cetera.
Set this form using, e.g.,
\code
apop_model *estimate_me = apop_model_copy(apop_binomial);
Apop_model_add_group(estimate_me, apop_rank);
apop_model *estimated = apop_estimate(your_data, estimate_me);
\endcode
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
                        */
    Estimate_results7, /**<  Parameters are estimated. Covariance matrix
                         is filled.   */
    Prep_routine7, /**<   None.      */
    RNG7, /**< Yes. */
    Key_settings7, /**<  \ref apop_rank_settings    */
    Example7 /**<      */
} ;

enum apop_null{
    Overview0,         /**< Sometimes you need a null model. the name is
                        blank, so <tt>strlen(m->name)</tt> will tell you
                        whether this is the null model.*/
    Name0,         /**< An empty string. (which is distinct from \c NULL.)*/
    Data_format0,  /**<         */
    Parameter_format0, /**<     */
    Estimate_results0, /**<     */
    Prep_routine0, /**<         */
    RNG0, /**< No. */
    Key_settings0, /**<      */
    Example0 /**<      */
} ;
