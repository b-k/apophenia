#include <apophenia/db.h>
#include <apophenia/model.h>
#include <apophenia/types.h>
#include <apophenia/stats.h>
#include <apophenia/output.h>
#include <apophenia/output.h>
#include <apophenia/histogram.h>
#include <apophenia/regression.h>
#include <apophenia/conversions.h>
#include <apophenia/likelihoods.h>
#include <apophenia/linear_algebra.h>

//Part of the intent of a convenience header like this is that you
//don't have to remember what else you're including. So here are 
//some other common GSL headers:
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>

//And here are other things that apophenia uses, and you probably do
//too.
#include <stdio.h>
#include <assert.h>
