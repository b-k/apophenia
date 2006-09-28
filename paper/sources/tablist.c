#include <apophenia/headers.h>

void print_table_list(char *db_file){
char            ***tab_list;
int             row_ct, i;
        apop_db_open(db_file);
        tab_list= apop_query_to_chars("select name from sqlite_master where type==\"table\";");
        row_ct  =  apop_db_get_rows();
        for(i=0; i< row_ct; i++)
                printf("%s\n", tab_list[i][0]);
}

int main(int argc, char **argv){
    if (argc==1){
        printf("Usage: give one database as an argument and I will list its tables.\n");
        return 0;
    }
    print_table_list(argv[1]);
    return 0;
}
