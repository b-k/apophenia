#ifdef Datadir
#define DATADIR Datadir
#else
#define DATADIR "."
#endif

#include <apop.h>
int main(){ apop_model_print(apop_estimate(apop_text_to_data( DATADIR "/" "data" ), apop_ols)); }
