#include "birds.h"
#include <glib.h>

  GList *flock  = NULL;

void flock_plays(){
    g_list_foreach(flock, bird_plays, NULL);
}

void add_to_flock(bird* b){
    flock   = g_list_prepend(flock, b);
}

void free_bird(bird* b){
    flock   = g_list_remove(flock, b);
    free(b);
}

bird * find_opponent(int n){
    return g_list_nth_data(flock, n);
}

int flock_size(){
    return g_list_length(flock);
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

void count(int period){
    hawks   = doves   = 0;
    g_list_foreach(flock, birth_or_death, NULL);
    g_list_foreach(flock, bird_count, NULL);
    apop_query("insert into pop values(%i, %i, %i);", period, hawks, doves);
}

void cull_flock(){
    flock   = g_list_remove_all(flock, NULL);
}
