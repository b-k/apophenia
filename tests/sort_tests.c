//for inclusion into sort_example.c

void check_sorting1(apop_data *d){
    double last_val = -INFINITY;
    double last_val2 = -INFINITY;
    for (int i=0; i < d->matrix->size1; i++){
        //Check for correct sort two columns deep.
        //Col 0 is all ones, so skip it.
        double this_val = apop_data_get(d, i, 1);
        double this_val2 = apop_data_get(d, i, 2);
        assert(this_val >= last_val);
        if (this_val == last_val)
            assert(this_val2 >= last_val2);
        last_val = this_val;
        last_val2 = this_val2;
    }
}

void check_sorting2(apop_data *d){
    double last_val = -INFINITY;
    double last_val2 = -INFINITY;
    char *last_str = "Aye";
    char *last_str2 = "Dem";
    for (int i=0; i < *d->textsize; i++){
        //Check for correct sort two columns deep.
        //Col 0 is all ones, so skip it.
        char *this_str = d->text[i][1];
        char *this_str2 = d->text[i][0];
        double this_val = apop_data_get(d, i, 1);
        double this_val2 = apop_data_get(d, i, 2);
        assert(strcasecmp(this_str, last_str) >=0);
        if (!strcasecmp(this_str, last_str)){
            assert(strcasecmp(this_str2, last_str2) >=0);
            if (!strcasecmp(this_str2,last_str2)){
                assert(this_val >= last_val);
                if(this_val == last_val)
                    assert(this_val2 >= last_val2);
            }
        }
        last_str = this_str;
        last_str2 = this_str2;
        last_val = this_val;
        last_val2 = this_val2;
    }
}
