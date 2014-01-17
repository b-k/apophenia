/** \file 
         Specifying model characteristics and details of estimation methods. */
/* Copyright (c) 2008--2009, 2011, 2013 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#include "apop_internal.h"

static size_t get_settings_ct(apop_model *model){
    int ct =0;
    if (!model->settings) return 0;
    while (model->settings[ct].name[0] !='\0') ct++;
    return ct;
}

//The Dan J Bernstein string hashing algorithm.
//Could conceivably save a lot of time under certain settings-heavy circumstances.
static unsigned long apop_settings_hash(char *str){
    unsigned long int hash = 5381;
    char c;
    while ((c = *str++)) hash = hash*33 + c;
    return hash;
}

/* Remove a settings group from a model.

Use \ref Apop_settings_rm_group. That macro uses this function internally.
*/
void apop_settings_remove_group(apop_model *m, char *delme){
    if (!m->settings) return;
    int i = 0;
    int ct = get_settings_ct(m);
    unsigned long delme_hash = apop_settings_hash(delme);
 
    while (m->settings[i].name[0] !='\0'){
        if (m->settings[i].name_hash == delme_hash){
            ((void (*)(void*))m->settings[i].free)(m->settings[i].setting_group);
            for (int j=i+1; j< ct+1; j++) //don't forget the null sentinel.
                m->settings[j-1] = m->settings[j];
            i--;
        }
        i++;
    }
   // apop_assert_void(0, 1, 'c', "I couldn't find %s in the input model, so nothing was removed.", delme);
}

/* Don't use this function. It's what the \c Apop_model_add_group macro uses internally. Use that.  */
void *apop_settings_group_alloc(apop_model *model, char *type, void *free_fn, void *copy_fn, void *the_group){
    if(apop_settings_get_grp(model, type, 'c'))  
        apop_settings_remove_group(model, type); 
    int ct = get_settings_ct(model);
    model->settings = realloc(model->settings, sizeof(apop_settings_type)*(ct+2));   
    model->settings[ct] = (apop_settings_type) {
                            .setting_group = the_group,
                            .name_hash = apop_settings_hash(type),
                            .free= free_fn, .copy = copy_fn };
    strncpy(model->settings[ct].name, type, 100);
    model->settings[ct+1] = (apop_settings_type) { };
    return model->settings[ct].setting_group;
}

//need this for the apop_settings_model_group_alloc macro.
apop_model *apop_settings_group_alloc_wm(apop_model *model, char *type, void *free_fn, void *copy_fn, void *the_group){
    apop_settings_group_alloc(model, type, free_fn, copy_fn, the_group);
    return model;
}

/* This function is used internally by the macro \ref Apop_settings_get_group. Use that.  */
void * apop_settings_get_grp(apop_model *m, char *type, char fail){
    //Used only for finding the non-blank groups.
    Apop_stopif(!m, return NULL, 0, "you gave me a NULL model as input.");
    if (!m->settings) return NULL;
    int i;
    unsigned long type_hash = apop_settings_hash(type);
    for (i=0; m->settings[i].name[0] !='\0'; i++)
       if (type_hash == m->settings[i].name_hash)
           return m->settings[i].setting_group;
    Apop_assert(fail != 'f', "I couldn't find the settings group %s in the given model.", type);
    return NULL; //else, just return NULL and let the caller sort it out.
}

/** Copy a settings group with the given name from the second model to
the first.  (i.e., the arguments are in memcpy order). 

You probably won't need this often---just use \ref apop_model_copy.

\param outm The model that will receive a copy of the settings group.
\param inm The model that will provide the original.
\param copyme The string naming the group. For example, for an \ref apop_mcmc_settings group, this would be \c "apop_mcmc".

\exception outm->error=='s'  Error copying settings group.
*/
void apop_settings_copy_group(apop_model *outm, apop_model *inm, char *copyme){
    if (!copyme || !strlen(copyme)) return; //apop_settings_group_alloc takes care of the blank sentinel.
    Apop_stopif(!inm, if (outm) outm->error = 's'; return, 0, "you asked me to copy the settings of a NULL model.");
    Apop_stopif(!inm->settings, return, 0, "The input model (i.e., the second argument to this function) has no settings.");
    void *g =  apop_settings_get_grp(inm, copyme, 'c');
    Apop_stopif(!g, outm->error='s'; return, 0, "Couldn't find the group you wanted me to copy. Not copying anything; setting outmodel->error='s'.");
    int i;
    unsigned long type_hash = apop_settings_hash(copyme);
    for (i=0; inm->settings[i].name[0] !='\0'; i++)//retrieve the index.
       if (type_hash == inm->settings[i].name_hash)
           break;
    void *gnew = (inm->settings[i].copy) 
                    ? ((void *(*)(void*))inm->settings[i].copy)(g)
                    : g;
    apop_settings_group_alloc(outm, copyme, inm->settings[i].free, inm->settings[i].copy, gnew);
}
