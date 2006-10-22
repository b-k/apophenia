#include "birds.h"
#include <glib.h>

  bird *flock;
  int   size_of_flock;

void flock_plays(){
  int   i;
    for (i=0; i< size_of_flock; i++)
        bird_plays(&(flock[i]), NULL);
}

void add_to_flock(bird* b){
    size_of_flock       = b->id;
    flock               = realloc(flock, sizeof(bird)*(size_of_flock+1));
    memcpy(&(flock[b->id]),b,sizeof(bird));
    free(b);
}

void free_bird(bird* b){
    flock[b->id].id   = -1;
}

bird * find_opponent(int n){
    if (flock[n].id > 0)
        return &(flock[n]);
    else
        return NULL;
}

int flock_size(){
    return size_of_flock;
}

int hawks, doves;

void count(int period){
  int   i;
    hawks   = doves   = 0;
    for(i=0; i< size_of_flock-1; i++)
        if (flock[i].id>0){
            birth_or_death(&(flock[i]), NULL);
        }
    for(i=0; i< size_of_flock-1; i++)
        if (flock[i].id>0){
            if (flock[i].type == 'h')
                hawks ++;
            else
                doves ++;
        }
    apop_query("insert into pop values(%i, %i, %i);", period, hawks, doves);
}

void cull_flock(){ }

void flock_init(){
  flock         = NULL;
  size_of_flock = 0;
}
