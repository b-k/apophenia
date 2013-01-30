#include <apop.h>
int main(){
    char string1[] = "Hello. I am a string.";
    assert(apop_regex(string1, "hell"));
    apop_data *subs;
    apop_regex(string1, "(e).*I.*(xxx)*(am)", .substrings = &subs);
    //apop_data_show(subs);
    assert(!strcmp(subs->text[0][0], "e"));
    assert(!strlen(subs->text[0][1]));
    assert(!strcmp(subs->text[0][2], "am"));
    apop_data_free(subs);

    //Split a comma-delimited list, throwing out white space.
    //Notice that the regex includes only one instance of a non-comma blob 
    //ending in a non-space followed by a comma, but the function keeps 
    //applying it until the end of string.
    char string2[] = " one, two , three ,four";
    apop_regex(string2, " *([^,]*[^ ]) *(,|$) *", &subs);
    assert(!strcmp(*subs->text[0], "one"));
    assert(!strcmp(*subs->text[1], "two"));
    assert(!strcmp(*subs->text[2], "three"));
    assert(!strcmp(*subs->text[3], "four"));
    apop_data_free(subs);

    //Get a parenthetical. For EREs, \( \) match plain parens in the text.
    char string3[] = " one (but secretly, two)";
    apop_regex(string3, "(\\([^)]*\\))", &subs);
    assert(!strcmp(*subs->text[0], "(but secretly, two)"));
    apop_data_free(subs);

    //NULL input string ==> no-op.
    int match_count = apop_regex(NULL, " *([^,]*[^ ]) *(,|$) *", &subs);
    assert(!match_count);
    assert(!subs);
}
