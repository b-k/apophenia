/* \file apop_probit.c

Copyright (c) 2005--2008 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

//The default list. Probably don't need them all.
#include "types.h"
#include "conversions.h"
#include "likelihoods.h"
#include "model.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>

/* \section dataprep Data prep rules 

There are a lot of ways your data can come in, and we would like to run estimations on a reasonably standardized form.

First, this page will give a little rationale, which you are welcome to skip, and then will present the set of rules.

  \paragraph Dealing with the ones column
Most standard regression-type estimations require or generally expect
a constant column. That is, the 0th column of your data is a constant (one), so the first parameter
\f$\beta_1\f$ is slightly special in corresponding to a constant rather than a variable.

However, there are some estimations that do not use the constant column.

\em "Why not implicitly assume the ones column?"
Some stats packages implicitly assume a constant column, which the user never sees. This violates the principle of transparency
upon which Apophenia is based, and is generally annoying.  Given a data matrix \f$X\f$ with the estimated parameters \f$\beta\f$, 
if the model asserts that the product \f$X\beta\f$ has meaning, then you should be able to calculate that product. With a ones column, a dot product is one line \c apop_dot(x, your_est->parameters, 0, 0)); without a ones column, the problem is left as an unpleasant exercise for the reader. You want to know if the regression included a ones column or not? Just look at your data.

  \paragraph Shunting columns around.

Each regression-type estimation has one dependent variable and several
independent. In the end, we want the dependent variable to be in the
vector element. However, continuing the \em "lassies faire" tradition,
doing major surgery on the data, such as removing a column and moving
in all subsequent columns, is more invasive than an estimation should be.

\subsection The rules
So those are the two main considerations in prepping data. Here are the rules,
intended to balance those considerations:

\paragraph The automatic case
There is one clever trick we can use to resolve both the need for a
ones column and for having the dependent column in the vector: given a
data set with no vector element and the dependent variable in the first
column of the matrix, we can copy the dependent variable into the vector
and then replace the first column of the matrix with ones. The result
fits all of the above expectations.

You as a user merely have to send in a \c apop_data set with no vector and a dependent column in the first column.

\paragraph The already-prepped case
If your data has a vector element, then the prep routines won't try
to force something to be there. That is, they won't move anything,
and won't turn anything into a constant column. If you don't want to use
a constant column, or your data has already been prepped by an estimation, then this is what you want.

You as a user just have to send in a \c apop_data set with a filled vector element.

\paragraph Probit, logit, and other -obits
The dependent variable for these models is a list of categories; this
option is not relevant to continuous-valued dependent variables. 

Your data source may include a list of numeric categories, in which case
you can pick one of the above cases.

The main exception is when your data is a list of text factors, in which
case your dependent variable isn't even a part of the data matrix.
In this case, you can prep the data yourself, via a call to \c apop_text_to_factors, and then insert the column yourself (and pick a constant column or not as you prefer). Or, you can have the system do it for you, via a form like
\code
int textcol = 3; //where is the list of dependent categories
apop_model *setmodel = apop_model_copy(apop_probit);
Apop_settings_add_group(setmodel, apop_category, textcol);
apop_estimate(yourdata, setmodel);
\endcode


  You'll see that there are two questions here: should there be a constant
  column of ones, and where is the dependent column to be found?

Here are the rules for preparing the data set. 

The first item is the 


There are two methods:

\li If the data set has no vector, 


for the -obit and -ogit models

   If there is a vector in place, then I won't touch anything.

   If there is no vector in place, then:
        --If you don't tell me where to find the dependent column, I'll go with column zero, and
            --move the data to the vector
            --replace the data there with ones, creating a constant column.
        --If you do tell me where to find the dependent column, via the settings, I'll turn that into a list of factors.

 */


 /* For the -obit part of the above, we need a struct to hold the
 categories. So, here's the apop_category_settings struct's methods. */
