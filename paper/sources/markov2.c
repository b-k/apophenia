#include <apophenia/headers.h>
int main(){
    apop_data *t = apop_text_to_data("markov_data", 0, 0);
    apop_data *out  = apop_dot(t, t, 0, 0);
    apop_data_show(out);
    return 0;
}
