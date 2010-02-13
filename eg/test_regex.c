void test_regex(){
    char string1[] = "Hello. I am a string.";
    assert(apop_regex(string1, "hell"));
    apop_data *subs;
    apop_regex(string1, "(e).*I.*(xxx)*(am)", .substrings = &subs);
    //apop_data_show(subs);
    assert(apop_strcmp(subs->text[0][0], "e"));
    assert(!strlen(subs->text[1][0]));
    assert(apop_strcmp(subs->text[2][0], "am"));
    apop_data_free(subs);
}
