m4_define(make_vtab_fns, <|m4_dnl
static void $1_type_check($1_type in){ };
#define $1_insert(fn, ...) ($1_type_check(fn), apop_vtable_insert("$1", fn, $1_hash(__VA_ARGS__))) 
#define $1_get(...) apop_vtable_get("$1", $1_hash(__VA_ARGS__))m4_dnl
|>)

typedef apop_model *(*apop_update_type)(apop_data *, apop_model , apop_model);
#define apop_update_hash(m1, m2) ((size_t)(m1).draw + (size_t)((m2).log_likelihood ? (m2).log_likelihood : (m2).p)*33)
make_vtab_fns(apop_update)

typedef void (*apop_score_type)(apop_data *d, gsl_vector *gradient, apop_model *params);
#define apop_score_hash(m1) ((size_t)((m1).log_likelihood ? (m1).log_likelihood : (m1).p))
make_vtab_fns(apop_score)

int apop_vtable_insert(char *tabname, void *fn_in, unsigned long hash);
void *apop_vtable_get(char *tabname, unsigned long hash);
