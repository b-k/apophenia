#include <apop.h>

void print_table_list(char *db_file){
    apop_db_open(db_file);
    apop_data *tab_list= apop_query_to_text("select name "
                    "from sqlite_master where type=='table'");
    for(int i=0; i< tab_list->textsize[0]; i++)
        printf("%s\n", tab_list->text[i][0]);
}

int main(int argc, char **argv){
    if (argc == 1){
        printf("Give me a database name, and I will print out "
               "the list of tables contained therein.\n");
        return 0; 
    }
    print_table_list(argv[1]);
}
