#include <apophenia/db.h>
#include <apophenia/asst.h>
#include <apophenia/model.h>
#include <apophenia/types.h>
#include <apophenia/stats.h>
#include <apophenia/output.h>
#include <apophenia/mapply.h>
#include <apophenia/settings.h>
#include <apophenia/regression.h>
#include <apophenia/conversions.h>
#include <apophenia/likelihoods.h>
#include <apophenia/linear_algebra.h>

//Part of the intent of a convenience header like this is that you
//don't have to remember what else you're including. So here are 
//some other common GSL headers:
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>

//And common headers for other uses (such as seeding an RNG):
#include <time.h>
#include <unistd.h>
