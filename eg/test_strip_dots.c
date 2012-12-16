#include <apop.h>
int main(){
    /* 0: replace all dots with _
      1: everything before the last dot.
      2: everything after the first dot.
      */
    char teapot[] = "tea.pot";
    char many_dots[] = "tea.pot.csv";
    char *out;
    out = apop_strip_dots(teapot, 0);
    assert(!strcmp(out, "tea_pot"));
    out = apop_strip_dots(teapot, 1);
    assert(!strcmp(out, "tea"));
    out = apop_strip_dots(teapot, 2);
    assert(!strcmp(out, "pot"));
    out = apop_strip_dots(many_dots, 0);
    assert(!strcmp(out, "tea_pot_csv"));
    out = apop_strip_dots(many_dots, 1);
    assert(!strcmp(out, "tea.pot"));
    out = apop_strip_dots(many_dots, 2);
    assert(!strcmp(out, "pot.csv"));
}
