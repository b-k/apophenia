#include <stdlib.h>
typedef struct {
    size_t hash;
    void *fn;
} apop_vtable_elmt_s;

typedef struct {
    char const *name;
    unsigned long int hashed_name;
    int elmt_ct;
    apop_vtable_elmt_s *elmts;
} apop_vtable_s;

apop_vtable_s *vtable_list;

//The Dan J Bernstein string hashing algorithm.
static unsigned long apop_settings_hash(char const *str){
    unsigned long int hash = 5381;
    char c;
    while ((c = *str++)) hash = hash*33 + c;
    return hash;
}

int apop_vtable_insert(char const *tabname, void *fn_in, unsigned long hash){
    if (!vtable_list){vtable_list = calloc(1, sizeof(apop_vtable_s));}

    //find the table we want
    apop_vtable_s *v = vtable_list;
    unsigned long h = apop_settings_hash(tabname);
    int ctr = 0;
    for ( ; v->hashed_name; ctr++, v++) if (v->hashed_name== h) break;

    //add a table if need be.
    if (!v->hashed_name){
        vtable_list=realloc(vtable_list, (ctr+2)* sizeof(apop_vtable_s));
        vtable_list[ctr] = (apop_vtable_s){.name=tabname, .hashed_name = h, .elmts=calloc(1, sizeof(apop_vtable_elmt_s))};
        vtable_list[ctr+1] = (apop_vtable_s){ };
        v = vtable_list+ctr;
    }

    //insert
    v->elmts = realloc(v->elmts, (++(v->elmt_ct))* sizeof(apop_vtable_elmt_s));
    v->elmts[v->elmt_ct-1] = (apop_vtable_elmt_s){.hash=hash, .fn=fn_in};
    return 0;
}

void *apop_vtable_get(char const *tabname, unsigned long hash){
    if (!vtable_list) return NULL;
    unsigned long thash = apop_settings_hash(tabname);
    apop_vtable_s *v = vtable_list;
    for ( ; v->hashed_name; v++) if (v->hashed_name== thash) break;
    if (!v->hashed_name) return NULL;

    for (int i=0; i< v->elmt_ct; i++)
        if (hash == v->elmts[i].hash) return v->elmts[i].fn;
    return NULL;
}
