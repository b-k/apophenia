#include "birds.h"
#include <apophenia/headers.h>

gsl_rng *r;

int     periods             = 300;
int     initial_pop         = 100;
int     id_count            = 0;

void play_pd_game(bird *row, bird *col){
  double    gain    = 2, loss = -1, 
            hawk_fight  = -1.95;
  if (row->type == 'd'){
        row->wealth += loss;
        col->wealth += gain;
  }
  if (col->type == 'd'){
        col->wealth += loss;
        row->wealth += gain;
  }
  if (col->type == 'h' && row->type == 'h'){
        col->wealth += hawk_fight;
        row->wealth += hawk_fight;
  }
}

void bird_plays(void *in, void *v){
  bird *other;
    while(!(other = find_opponent((int)(gsl_rng_uniform(r)*flock_size()))));
    if (!(in==other))
        play_pd_game(in, other);
}

bird *new_chick(bird *parent){
  bird  *out    = malloc(sizeof(*out));
    if (parent)
        out->type  = parent->type;
    else{
        if (gsl_rng_uniform(r) > 0.5)
            out->type  = 'd';
        else
            out->type  = 'h';
    }
    out->wealth = 20* gsl_rng_uniform(r);
    out->id     = id_count;
    id_count    ++;
    return out;
}

void birth_or_death(void *in, void *v){
  bird  *b  = in; //cast void to bird;
    //b->wealth   -= 0.25;
    if (b->wealth <=0)
        free_bird(b);
    else if (gsl_rng_uniform(r)*75 < b->wealth)
        add_to_flock(new_chick(b));
}

void startup(int initial_flock_size){
  int   i;
    flock_init();
    r   = apop_rng_alloc(923);
    remove("history.db");
    apop_db_open("history.db");
    apop_query("create table pop(period, hawks, doves); begin;");
    for(i=0; i< initial_flock_size; i++)
        add_to_flock(new_chick(NULL));
}

int main(){
  int   i;
    startup(initial_pop);
    for (i=0; i< periods; i++){
        flock_plays();
        count(i);
        if (!(i % 10)) {
            fprintf(stderr, "%i\n", i);
            apop_query("commit; begin;");
        }
    }
    apop_query("commit;");
    return 0;
}
