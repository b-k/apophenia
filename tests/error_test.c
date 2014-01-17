#include <apop.h>
#if _POSIX_C_SOURCE >= 200809L
#define HAVE_FMEMOPEN
#endif
FILE *fmemopen(void *buf, size_t size, const char *mode); //POSIX.

/* Set up an error, and check for the presence of:
  --the correct error flag in the data set
  --the name of the function in the error log.

I'm only checking for the function name in the error log because I don't want to rewrite this every time 
the wording in the error message changes.

If fmemopen is missing, don't even bother with the error log stuff.

*/
char errorbuff[10000];

void check_log(char*fn_to_check, char*msg){
#ifdef HAVE_FMEMOPEN
    fflush(NULL);
    Apop_stopif (!apop_regex(errorbuff, fn_to_check), abort(), 0, msg);
#endif
}

void check_data_error(apop_data *in, char should_be, char *fn_to_check, char *msg){
    Apop_stopif (in->error != should_be, abort(), 0, "Didn't set %s.", msg);
    check_log(fn_to_check, msg);
}

void reset_log(){
#ifdef HAVE_FMEMOPEN
    if (apop_opts.log_file) fclose(apop_opts.log_file);
    apop_opts.log_file = fmemopen(errorbuff, 10000, "w");
#endif
}

int main(){
    apop_opts.db_engine='s';
    printf("test error checking (some systems may print error messags here)\n");
    reset_log();
    apop_data *d = apop_data_alloc();
    apop_data *d2 = apop_data_add_page(d, apop_data_alloc(), "newp");
    d2->more = d2;
    apop_data_free(d);
    check_data_error(d, 'c', "apop_data_free_base", "circular error code");
    check_data_error(d2, 'c', "apop_data_free_base", "circular error code");

    reset_log();
    apop_data *c = apop_data_copy(d);
    check_data_error(c, 'c', "apop_data_copy", "circular error code");
    check_data_error(c->more, 'c', "apop_data_copy", "circular error code");

    reset_log();
    apop_data *d3 = apop_data_alloc(2,2);
    apop_data_memcpy(d, d3);
    check_data_error(d, 'p', "apop_data_memcpy", "missing part code");

    reset_log();
    apop_data *dbig = apop_data_alloc(2e8,2e8);
    check_data_error(dbig, 'a', "apop_data_alloc", "misallocation code "
                                "(or you are using a big mother of a computer).");

    reset_log();
    check_data_error(apop_query_to_data("stelect 8 from data"), 'q', "(apop_query_to_data|<unknown>)", "query error");
    reset_log();
    check_data_error(apop_query_to_mixed_data("dd", "stelect 8 from data"), 'q', "(apop_sqlite_multiquery|<unknown>)", "query error");

    reset_log();
    apop_data *fefail = apop_data_falloc((2,2), 0, 0, -1, -1);
    apop_data *exact = apop_test_fisher_exact(fefail);
    check_data_error(exact, 'p', "apop_test_fisher_exact", "fexact internal processing code");

    reset_log();
    apop_data *d44 = apop_data_alloc(4,4);
    check_data_error(apop_dot(d44, d3), 'd', "apop_dot", "dot product dimension error");

    apop_multivariate_normal->parameters = apop_data_calloc(2, 2, 2);
    assert(apop_log_likelihood(fefail, apop_multivariate_normal) == -INFINITY);
    check_log("apop_multinormal_ll", "Failed to not take the determinant of a zero matrix.");

    if (apop_opts.log_file) fclose(apop_opts.log_file);
    apop_opts.log_file = NULL;
}
