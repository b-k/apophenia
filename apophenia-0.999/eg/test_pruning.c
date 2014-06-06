#include <apop.h>

/* This sample produces a dummy times table, gets a summary, and prunes the summary table.
If you are not a test script, uncomment the last line to display the pruned table.  */
int main(){
    int i, j;
    apop_data *d = apop_data_alloc(0, 10, 4);
    for (i=0; i< 10; i++)
        for (j=0; j< 4; j++)
            apop_data_set(d, i, j, i*j);
    apop_data *summary = apop_data_summarize(d);
    apop_data_prune_columns(summary, "mean", "median");
    assert(apop_name_find(summary->names, "mean", 'c')!=-2);
    assert(apop_name_find(summary->names, "median", 'c')!=-2);
    assert(apop_name_find(summary->names, "max", 'c')==-2); //not found
    assert(apop_name_find(summary->names, "variance", 'c')==-2); //not found
    assert(apop_data_get(summary, .row=0, .colname="mean")==0);
    assert(apop_data_get(summary, .row=1, .colname="median")==4);
    assert(apop_data_get(summary, .row=2, .colname="median")==8);
    //apop_data_show(summary);
}
