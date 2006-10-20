#include "birds.h"
#include <glib.h>

  GTree *flock  = NULL;

gint compare_birds(const void *L, const void *R){
  const int *Lb  = L;
  const int *Rb  = R;
    if (*Rb - *Lb < 0)
        return -1;
    else return (*Rb > *Lb);
}

void flock_init(){
    flock   = g_tree_new(compare_birds);

}

static gboolean tree_bird_plays(void *k, void *in, void *v){
    bird_plays(in, NULL);
    return 0;
}

static gboolean tree_birth_or_death(void *k, void *in, void *v){
    birth_or_death(in, NULL);
    return 0;
}

void flock_plays(){
    g_tree_foreach(flock, tree_bird_plays, NULL);
}

void add_to_flock(bird* b){
    g_tree_insert(flock, &(b->id), b);
}

void free_bird(bird* b){
    g_tree_remove(flock, &(b->id));
    free(b);
}

bird * find_opponent(int n){
    return g_tree_lookup(flock, &n);
}

int flock_size(){
    return g_tree_nnodes(flock);
}

int hawks, doves;

void bird_count(void *in, void *v){
  bird *b   = in;
    if (b){
        if (b->type == 'h')
            hawks   ++;
        else
            doves   ++;
    }
}

static gboolean tree_bird_count(void *k, void *in, void *v){
    bird_count(in, NULL);
    return 0;
}

void count(int period){
    hawks   = doves   = 0;
    g_tree_foreach(flock, tree_birth_or_death, NULL);
    g_tree_foreach(flock, tree_bird_count, NULL);
    apop_query("insert into pop values(%i, %i, %i);", period, hawks, doves);
}

void cull_flock(){
    return; }
