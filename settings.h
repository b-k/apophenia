/** \file settings.h */ 
/* Copyright (c) 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#ifndef __apop_settings_h__
#define __apop_settings_h__

#include "types.h"
#include "asst.h"

#ifdef	__cplusplus
extern "C" {
#endif

void * apop_settings_get_group(apop_model *m, char *type);
void apop_settings_rm_group(apop_model *m, char *delme);
void apop_settings_copy_group(apop_model *outm, apop_model *inm, char *copyme);
void apop_settings_group_alloc(apop_model *model, char *type, void *free_fn, void *copy_fn, void *the_group);

/** Retrieves a settings group from a model.  See \ref Apop_settings_get
 to just pull a single item from within the settings group.*/
#define Apop_settings_get_group(m, type) apop_settings_get_group(m, #type)

/** Removes a settings group from a model's list. */
#define Apop_settings_rm_group(m, type) apop_settings_rm_group(m, #type)

/* If the group already exists, it is silently removed. */


/** Add a settings group. The first two arguments (the model you are
 attaching to and the settings group name) are mandatory, and then you can use the \ref designated syntax to specify default values (if any).
 */
#define Apop_model_add_group(model, type, ...)  \
    apop_settings_group_alloc(model, #type, type ## _settings_free, type ## _settings_copy, type ##_settings_init ((type ## _settings) {__VA_ARGS__})); 

/* A convenience for your settings group init functions. 
 Gives the output item either the defaut value if there is one, or the value you specify. */
#define apop_varad_setting(in, out, name, value) (out)->name = (in).name ? (in).name : (value);

/** Add a settings group. Deprecated; use \ref Apop_model_add_group instead.
  You will need to provide arguments for the
 specific settings group you are dealing with, such as \ref apop_mle_settings_alloc, \ref apop_ls_settings_alloc, \ref apop_histogram_settings_alloc. */
#define Apop_settings_add_group(model, type, ...)  \
    apop_settings_group_alloc(model, #type, type ## _settings_free, type ## _settings_copy, type ##_settings_alloc (__VA_ARGS__)); 

/** Retrieves a setting from a model.  See \ref Apop_settings_get_group pull the entire group.*/
#define Apop_settings_get(model, type, setting)  \
    (((type ## _settings *) apop_settings_get_group(model, #type))->setting)

/** Modifies a single element of a settings group to the given value. */
#define Apop_settings_set(model, type, setting, data)  \
    do {                                                \
    apop_assert_void(apop_settings_get_group(model, #type), 0, 's', "You're trying to modify a setting in " \
                        #model "'s setting group of type " #type " but that model doesn't have such a group."); \
    ((type ## _settings *) apop_settings_get_group(model, #type))->setting = (data);    \
    } while (0);

#define Apop_settings_add Apop_settings_set
#define APOP_SETTINGS_ADD Apop_settings_set
#define APOP_SETTINGS_GET Apop_settings_get
#define APOP_MODEL_ADD_GROUP Apop_model_add_group
#define APOP_SETTINGS_ADD_GROUP Apop_settings_add_group
#define APOP_SETTINGS_GET_GROUP Apop_settings_get_group
#define APOP_SETTINGS_RM_GROUP Apop_settings_rm_group


#define Apop_settings_declarations(ysg) \
   ysg##_settings * ysg##_settings_init(ysg##_settings); \
   void * ysg##_settings_copy(ysg##_settings *); \
   void ysg##_settings_free(ysg##_settings *);

/** Settings for least-squares type models */
typedef struct {
    int destroy_data;
    gsl_vector *weights;
    apop_data *instruments;
    char want_cov;
    char want_expected_value;
} apop_ls_settings;

Apop_settings_declarations(apop_ls)

// Find apop_category_settings routines in apop_probit.c
/** for dependent-category models, send in this settings struct to
  specify which column is the dependent variable. 

  If you don't use it, these models will assume that the vector or first
  numeric column is already a coherent set of factors, but by sending
  this in, those functions have a little more information, such as names
  to use in the output.

  See also the \ref apop_category_settings_alloc function.
 */
typedef struct {
    apop_data *factors; 
    char source_type; 
    char source_column; 
    apop_data *source_data;
} apop_category_settings;

Apop_settings_declarations(apop_category)
apop_category_settings *apop_category_settings_alloc(apop_data *d, int source_column, char source_type);

/** If this settings group is present, models that can take rank data
  will read the input data as such.  Allocation is thus very simple, e.g.
  \code
  Apop_settings_group_add(your_model, apop_rank, NULL);
  \endcode */
typedef struct {
    char rank_data;
} apop_rank_settings;

//in apop_exponential.c
Apop_settings_declarations(apop_rank)
apop_rank_settings *apop_rank_settings_alloc(void *ignoreme);

#include <gsl/gsl_histogram.h>
typedef struct{
    apop_data           *data;
    gsl_histogram       *pdf;
    gsl_histogram_pdf   *cdf;
    apop_model          *histobase;
    apop_model          *kernelbase;
} apop_histogram_settings;

#define apop_kernel_density_settings apop_histogram_settings

Apop_settings_declarations(apop_histogram)
apop_histogram_settings *apop_histogram_settings_alloc(apop_data *data, int bins);


/** Method settings for a model to be put through Bayesian updating. 
\param starting_pt      The first parameter to check in the MCMC routine
\param periods How many steps should the MCMC chain run?
\param burnin  What <em>percentage</em> of the periods should be ignored as initialization. That is, this is a number between zero and one.
\param histosegments If outputting a \ref apop_histogram, how many segments should it have?
 
 */
typedef struct{
    apop_data *data;
    apop_data *starting_pt;
    long int periods;
    double burnin;
    int histosegments;
    char method;
} apop_update_settings;

apop_update_settings *apop_update_settings_alloc(apop_data *d);
apop_update_settings *apop_update_settings_init(apop_update_settings);
#define apop_update_settings_copy NULL
#define apop_update_settings_free NULL


apop_histogram_settings *apop_kernel_density_settings_alloc(apop_data *data, 
        apop_model *histobase, apop_model *kernelbase, void (*set_params)(double, apop_model*));

#define apop_kernel_density_settings_copy apop_histogram_settings_copy
#define apop_kernel_density_settings_free apop_histogram_settings_free

#ifdef	__cplusplus
}
#endif
#endif
