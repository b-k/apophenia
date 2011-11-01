void test_apop_strcmp(){
    assert(!apop_strcmp("23", NULL));
    assert(!apop_strcmp("230", "23"));
    assert(!apop_strcmp("super", " uper"));
    assert(apop_strcmp(NULL, NULL));
    assert(apop_strcmp("23", "23"));
    assert(apop_strcmp("", ""));
}
