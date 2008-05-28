/** \file apop_settings.c Specifying model characteristics and details of estimation methods.
 
Copyright (c) 2008 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#include "types.h"
#include "settings.h"

/** \page settings Tweaks and metadata for models.

Apophenia is really only based on two objects, the \ref apop_data
set and the \ref apop_model. Data sets come in a pretty standard form,
so the data object is basically settled. But describing a statistical,
agent-based, social, or physical model in a reasonably standardized form
is much more difficult, primarily because every model has significantly
different settings. E.g., an MLE requires a method of search (conjugate
gradient, simplex, simulated annealing), and a histogram needs the number
of slots to be filled with data.

So, the \ref apop_model includes 
a single list, whose name is simply \c settings, which can hold an arbitrary
number of groups of settings. For example, you can have a set of
closed-form variables for estimating the log likelihood, and a histogram
for making random draws.

To get/set a setting, you would need to specify the model, the settings
group, and the name of the setting itself.  For the sake of giving you
a mental model, this is much like a g_lib config file (or for Windows users,
a .ini file), which would have settings divided into sections, like

\code
[MLE]
verbose = 0
method = APOP_CG_PR

[least squares]
weights = your_weight_vector

[histogram]
bins = {0.1, 0.3, 0.2,...}
\endcode


\section settings_using  Using it

If you don't need to change settings from the default, you don't need to
care about any of this, because
\code
apop_model_show(apop_estimate(your_data_here, apop_ols));
\endcode
---still works fine.

If you do need to change settings, then the process takes another step
or two. Here's a sample:

\code
1 apop_data *data = your_data_here;
2 apop_data *w = your_weights_here;

3 apop_model *m = apop_model_copy(apop_wls);
4 Apop_settings_add_group(m, apop_ls, data);
5 Apop_settings_add(m, apop_ls, weights, w);
6 apop_model *est = apop_estimate(data, *m);
\endcode

Line three establishes the baseline form of the model. Line four adds
a settings group of type apop_ls_settings to the model. If you check the
manual for apop_ls_settings_alloc, you'll see that it takes one
argument---the input data---and the inputs to the alloc function appear
after the model and settings group name. Line five sets the weights
element in that group to w, and looks like the config-file format above.

Also, the output to any of Apophenia's estimations will have an
appropriate group of settings allocated, so you can chain estimations
pretty easily. Continuing the above example, you could re-estimate with
an alternative set of weights via:

\code
Apop_settings_add(est, apop_ls, weights, weight_set_two);
apop_model *est2 = apop_estimate(data, *est);
\endcode

[Notice that the 'add' function really just modifies.]

Some people are incredibly averse to additional lines of code, so there's
a convenience to allocate and modify a setting on the same line (joining
lines 4 and 5 above):

\code
Apop_settings_alloc_add(m, apop_ls, weights, w, data);
\endcode

If you need to read a setting, such as to check whether things have
changed after an estimation, use:

\code
Apop_settings_get(m, apop_ls, weights);
\endcode

Notice the use of a single capital to remind you that you are using a
macro, and so surprising errors can crop up. Here in the modern day, we
read things like APOP_SETTINGS_ALLOC_ADD as yelling, but if you prefer
all caps to indicate macros, those work as well.

For just using a model, that's about 100% of what you need to know.


\section settings_writng  "Writing new settings"

To store the settings for your own models, you don't necessarily
need any of this. The apop_model structure has a void pointer named
more which you can use as you see fit. If more_size is larger than zero
(i.e., you set it to sizeof(your_struct)), then it will be copied via
memcpy as necessary. Apohenia's estimation routines will never impinge
on this item, so do what you feel with it.

If you do want to set up a new model, then you will need four items. 
This is the sort of boilerplate that will be familiar to users of 
object oriented languages in the style of C++ or Java. Let
your settings group be named ysg; then you will need

--The settings struct
\code
typedef struct {
...
} ysg_settings
\endcode


--The allocate function
\code
ysg_settings *ysg_settings_alloc(...){ 
    ysg_settings *out = malloc(sizeof(ysg_settings));
    // Give default values for all elements here. 
    return out; }
\endcode


--The copy function
\code
void *ysg_settings_copy(ysg_settings *copyme) {
    ysg_settings *out = malloc(sizeof(ysg_settings));
    ...
    return out; }
\endcode


--The free function, which can be as brief as:
\code
void ysg_settings_free(ysg_settings *copyme) {
    free(copyme);
}
\endcode
(but may include a freeing of pointed-to subelements as necessary).

The names are not negotiable: when you call
\code
Apop_settings_alloc(m, ysg, ...)
\endcode
the macro will look for ysg_settings, ysg_settings_alloc, et cetera.

The lines-of-code averse will cringe at having to write such boilerplate
code (I do), but after spending a year resisting it, I have to concede
that it's the least of all evils.

You can retrieve the whole group or individual elements via:
\code
Apop_settings_get_group(m, ysg)
//as well as the one-element macro mentioned above,
Apop_settings_get(m, ysg, an_element)
\endcode

As you saw above, once the typedef/alloc/copy/free machinery is written,
you can declare, get, and set in a reasonably graceful manner.


\ingroup settings
 */

void apop_settings_copy_group(apop_model *outm, apop_model *inm, char *copyme){
  int i;
    for (i=0; i< inm->setting_ct; i++)
        if (!strcmp(inm->settings[i].name, copyme)){
            apop_settings_type c = inm->settings[i];
            outm->settings = realloc(outm->settings, sizeof(apop_settings_type)*(outm->setting_ct+1)); 
            strncpy(outm->settings[outm->setting_ct].name , copyme, 100) ; 
            if (c.copy)
                outm->settings[outm->setting_ct].setting_group = ((void *(*)(void*))c.copy)(c.setting_group);
            outm->settings[outm->setting_ct].free = c.free;
            outm->settings[outm->setting_ct].copy = c.copy;
            outm->setting_ct ++; 
            return;
        }
    apop_assert_void(0, 1, 'c', "I couldn't find %s in the input model, so nothing was copied.", copyme);
}

void apop_settings_rm_group(apop_model *m, char *delme){
  int i, j;
    for (i=0; i< m->setting_ct; i++)
        if (!strcmp(m->settings[i].name, delme)){
            ((void (*)(void*))m->settings[i].free)(m->settings[i].setting_group);
            for (j=i+1; j< m->setting_ct; j++){
                strcpy(m->settings[m->setting_ct-1].name, m->settings[m->setting_ct].name);
                m->settings[m->setting_ct-1].free = m->settings[m->setting_ct].free;
                m->settings[m->setting_ct-1].copy = m->settings[m->setting_ct].copy;
                m->settings[m->setting_ct-1].setting_group = m->settings[m->setting_ct].setting_group;
            }
            m->setting_ct --; 
            return;
        }
    apop_assert_void(0, 1, 'c', "I couldn't find %s in the input model, so nothing was removed.", delme);
}

void * apop_settings_get_group(apop_model *m, char *type){
  int   i;
    for (i=0; i< m->setting_ct; i++){
        if (!strcmp(type, m->settings[i].name))
            return m->settings[i].setting_group;
    }
    return NULL;
}
