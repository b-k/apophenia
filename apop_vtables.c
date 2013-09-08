#include <stdlib.h>
#include <string.h>
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
int ignore_me;

//The Dan J Bernstein string hashing algorithm.
static unsigned long apop_settings_hash(char const *str){
    unsigned long int hash = 5381;
    char c;
    while ((c = *str++)) hash = hash*33 + c;
    return hash;
}

static apop_vtable_s *find_tab(unsigned long h, int *ctr){
    apop_vtable_s *v = vtable_list;
    *ctr = 0;
    for ( ; v->hashed_name; (*ctr)++, v++) if (v->hashed_name== h) break;
    return v;
}

//return 0 = found; removed
//return 1 = not found; no-op
int apop_vtable_drop(char const *tabname, unsigned long hash){
    if (!vtable_list) return 1;
    unsigned long h = apop_settings_hash(tabname);
    apop_vtable_s *v = find_tab(h, &ignore_me);

    for (int i=0; i< v->elmt_ct; i++)
        if (hash == v->elmts[i].hash) {
            memmove(v->elmts+i, v->elmts+i+1, sizeof(apop_vtable_elmt_s)*(v->elmt_ct-i));
            v->elmt_ct--;
            return 0;
        }
    return 1;
}

int apop_vtable_add(char const *tabname, void *fn_in, unsigned long hash){
    if (!vtable_list){vtable_list = calloc(1, sizeof(apop_vtable_s));}

    unsigned long h = apop_settings_hash(tabname);
    int ctr;
    apop_vtable_s *v = find_tab(h, &ctr);

    //add a table if need be.
    if (!v->hashed_name){
        vtable_list = realloc(vtable_list, (ctr+2)* sizeof(apop_vtable_s));
        vtable_list[ctr] = (apop_vtable_s){.name=tabname, .hashed_name = h, .elmts=calloc(1, sizeof(apop_vtable_elmt_s))};
        vtable_list[ctr+1] = (apop_vtable_s){ };
        v = vtable_list+ctr;
    }

    //If this hash is already present, don't re-add. 
    for (int i=0; i< v->elmt_ct; i++) if (hash == v->elmts[i].hash) return 0;


    //insert
    v->elmts = realloc(v->elmts, (++(v->elmt_ct))* sizeof(apop_vtable_elmt_s));
    v->elmts[v->elmt_ct-1] = (apop_vtable_elmt_s){.hash=hash, .fn=fn_in};
    return 0;
}

void *apop_vtable_get(char const *tabname, unsigned long hash){
    if (!vtable_list) return NULL;
    unsigned long thash = apop_settings_hash(tabname);
    apop_vtable_s *v = find_tab(thash, &ignore_me);
    if (!v->hashed_name) return NULL;

    for (int i=0; i< v->elmt_ct; i++)
        if (hash == v->elmts[i].hash) return v->elmts[i].fn;
    return NULL;
}
