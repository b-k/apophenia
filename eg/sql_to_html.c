#include <apop.h>

int main(){
    apop_query("create table datatab(name, age, sex);"
                "insert into datatab values ('Alex', 23, 'm');"
                "insert into datatab values ('Alex', 32, 'f');"
                "insert into datatab values ('Michael', 41, 'f');"
                "insert into datatab values ('Michael', 14, 'm');");

    apop_data *cols = apop_text_alloc(NULL, 3, 1);
    apop_text_add(cols, 0, 0, "name");
    apop_text_add(cols, 1, 0, "age");
    apop_text_add(cols, 2, 0, "sex");
    char *query= apop_text_paste(cols, .before="select ", .between=", ");
    apop_data *d = apop_query_to_text("%s from datatab", query);
    char *html_head = apop_text_paste(cols, .before="<table><tr><td>",
                                .between="</td><td>", .after="</tr>\n<tr><td>");
    char *html_table = apop_text_paste(d, .before=html_head, .after="</td></tr></table>\n",
                                .between="</tr>\n<tr><td>", .between_cols="</td><td>");
    FILE *outfile = fopen("yourdata.html", "w");                                                                                                                              fprintf(outfile, "%s", html_table);
    fclose(outfile);
}
