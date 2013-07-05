#include <apop.h>

int main(){
    apop_query("create table data (name, city, state);"
            "insert into data values ('Mike Mills', 'Rockville', 'MD');"
            "insert into data values ('Bill Berry', 'Athens', 'GA');"
            "insert into data values ('Michael Stipe', 'Decatur', 'GA');");
    apop_data *tdata = apop_query_to_text("select name, city, state from data");
    printf("Customer #1: %s\n\n", *tdata->text[0]);

    printf("The data, via apop_data_print:\n");
    apop_data_print(tdata);

    //the text alloc can be used as a text realloc:
    apop_text_alloc(tdata, 1+tdata->textsize[0], tdata->textsize[1]);
    apop_text_add(tdata, *tdata->textsize-1, 0, "Peter Buck");
    apop_text_add(tdata, *tdata->textsize-1, 1, "Berkeley");
    apop_text_add(tdata, *tdata->textsize-1, 2, "CA");

    printf("\n\nAugmented data, printed via for loop (for demonstration purposes):\n");
    for (int i=0; i< tdata->textsize[0]; i++){
        for (int j=0; j< tdata->textsize[1]; j++)
            printf("%s\t", tdata->text[i][j]);
        printf("\n");
    }

    apop_data *states = apop_text_unique_elements(tdata, 2);
    char *states_as_list = apop_text_paste(states, .between=", ");
    printf("\n States covered: %s\n", states_as_list);
}
