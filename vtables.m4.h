#ifndef VTABLES_H
#define VTABLES_H

/*
This declares the vtable macros for each procedure that uses the mechanism.

--We want to have type-checking on the functions put into the vtables. Type checking
happens only with functions, not macros, so we need a type_check function for every
vtable.

--Only once in your codebase, you'll need to #define Declare_type_checking_fns to
actually define the type checking function. Everywhere else, the function is merely
declared.

--All other uses point to having a macro, such as using __VA_ARGS__ to allow any sort
of inputs to the hash.

--We want to have such a macro for every vtable. That means that we need a macro
to write macros. We can't do that with C macros.  Thus, this file uses m4 macros to
generate C macros.

--After the m4 definition of make_vtab_fns, each new vtable requires a typedef, a hash
definition, and a call to make_vtab_fns to do the rest.

--the typedef and has are defined inside of an m4 macro so that we can also paste the
same text into documentation as needed.
*/
m4_define(make_vtab_fns, <|m4_dnl
#ifdef Declare_type_checking_fns
void $1_type_check($1_type in){ };
#else
void $1_type_check($1_type in);
#endif
#define $1_vtable_add(fn, ...) ($1_type_check(fn), apop_vtable_add("$1", fn, $1_hash(__VA_ARGS__))) 
#define $1_vtable_get(...) apop_vtable_get("$1", $1_hash(__VA_ARGS__))
#define $1_vtable_drop(...) apop_vtable_drop("$1", $1_hash(__VA_ARGS__))m4_dnl
|>)

int apop_vtable_add(char *tabname, void *fn_in, unsigned long hash);
void *apop_vtable_get(char *tabname, unsigned long hash);
int apop_vtable_drop(char const *tabname, unsigned long hash);

typedef apop_model *(*apop_update_type)(apop_data *, apop_model , apop_model);
#define apop_update_hash(m1, m2) ((size_t)(m1).draw + (size_t)((m2).log_likelihood ? (m2).log_likelihood : (m2).p)*33)
make_vtab_fns(apop_update)

typedef void (*apop_score_type)(apop_data *d, gsl_vector *gradient, apop_model *params);
#define apop_score_hash(m1) ((size_t)((m1).log_likelihood ? (m1).log_likelihood : (m1).p))
make_vtab_fns(apop_score)

typedef apop_model* (*apop_parameter_model_type)(apop_data *, apop_model *);
#define apop_parameter_model_hash(m1) ((size_t)((m1).log_likelihood ? (m1).log_likelihood : (m1).p)*33 + (m1).estimate ? (size_t)(m1).estimate: 27)
make_vtab_fns(apop_parameter_model)

typedef apop_data * (*apop_predict_type)(apop_data *d, apop_model *params);
#define apop_predict_hash(m1) ((size_t)((m1).log_likelihood ? (m1).log_likelihood : (m1).p)*33 + (m1).estimate ? (size_t)(m1).estimate: 27)
make_vtab_fns(apop_predict)

#endif
