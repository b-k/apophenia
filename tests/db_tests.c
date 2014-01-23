/*
These are database-related tests. If you have mySQL/mariaDB set up correctly (see below) then you can run them both for SQLite and m... database engines.

SQLite is the default for Apophenia, and this works out of the box for it.
For mySQL/mariaDB, you will need to either hack this to include
apop_opts.user=[your username]
apop_opts.pass=[your password]

or set up a config file, like pasting this to ~/.my.cnf :
[client]
user=[your name]
pass=[your pass]

Then, either hack this again to include
apop_opts.db_engine='m';

or at the command line:
export APOP_DB_ENGINE=mysql
*/

#include <apop.h>
#include <unistd.h>
int verbose = 1;

#define Diff(L, R, eps) {double left=(L), right=(R); Apop_stopif(isnan(left-right) || fabs((left)-(right))>(eps), abort(), 0, "%g is too different from %g (abitrary limit=%g).", (double)(left), (double)(right), eps);}

void test_data_to_db() {
  int i, j;
    if (!apop_table_exists("snps"))
        apop_text_to_db("data-mixed", "snps");
    apop_data *d = apop_query_to_mixed_data("tvttmmmt", "select * from snps");
    apop_data_print(d, "snps2", .output_type='d');
    apop_data *d2 = apop_query_to_mixed_data("vmmmtttt", "select * from snps2");
    for (i=0; i< d2->vector->size; i++)
        assert(d->vector->data[i] == d2->vector->data[i]);
    for (i=0; i< d2->matrix->size1; i++)
        for (j=0; j< d2->matrix->size2; j++)
            assert(gsl_matrix_get(d->matrix, i, j) ==  gsl_matrix_get(d2->matrix, i, j));
    for (i=0; i< d2->textsize[0]; i++)
        for (j=0; j< d2->textsize[1]; j++)
            assert(!strcmp(d->text[i][j],d2->text[i][j]));  
    unlink("snps2");
}

void test_uniform(apop_data *d){
    Apop_col_tv(d, "ab", abcol);
    apop_data ab_d = (apop_data){.vector=abcol};
    apop_model *u = apop_estimate(&ab_d, apop_uniform);
    Diff(log(apop_p(&ab_d, u)), apop_log_likelihood(&ab_d, u), 1e-5);

    apop_data_add_names(&ab_d, 'v', "a vector");
    apop_data_set(&ab_d, .colname="a vector", .row=0, .val=-297);
    assert(apop_p(&ab_d, u) == 0);
    assert(isinf(apop_log_likelihood(&ab_d, u)));

    apop_model *iu = apop_estimate(&ab_d, apop_improper_uniform);
    assert(apop_p(&ab_d, iu) == 1);
    assert(apop_log_likelihood(&ab_d, iu)==0);

    int verbosity = apop_opts.verbose;
    apop_opts.verbose = -1;
    double draw;
    apop_draw(&draw, NULL, iu);
    apop_opts.verbose = verbosity;
    assert(isnan(draw));
}


void db_to_text(){
    apop_db_close();
    apop_db_open(NULL);
    if (!apop_table_exists("d")){
        apop_data *field_params = apop_text_alloc(NULL,2,2);
        apop_text_fill(field_params, 
                "[ab][ab]", "numeric",
                ".*",     apop_opts.db_engine =='s' ? "character": "varchar(20)"
                );
        apop_text_to_db("data-mixed", "d", 0, 1, NULL, .field_params=field_params);
    }
    apop_data *d = apop_query_to_mixed_data ("tmttmmmt", "select * from d");
    int b_allele_col = apop_name_find(d->names, "b_allele", 't');
    assert(!strcmp("T",  d->text[3][b_allele_col]));
    int rsid_col = apop_name_find(d->names, "rsid", 't');
    assert(!strcmp("rs2977656",  d->text[4][rsid_col]));
    assert(apop_data_get(d, .row=5, .colname="ab")==201);

    assert(!strcmp(d->text[3][rsid_col], "rs'11804171"));

    apop_data *dcc = apop_data_copy(d); //test apop_data_copy
    assert(!strcmp("T",  dcc->text[3][b_allele_col]));
    assert(!strcmp("rs2977656",  dcc->text[4][rsid_col]));
    assert(apop_data_get(dcc, 5, .colname="ab")==201);

    apop_data *dd = apop_query_to_text ("select * from d");
    b_allele_col = apop_name_find(dd->names, "b_allele", 't');
    assert(!strcmp("T",  dd->text[3][b_allele_col]));
    rsid_col = apop_name_find(dd->names, "rsid", 't');
    assert(!strcmp("rs2977656",  dd->text[4][rsid_col]));
    
    apop_data *dc = apop_data_copy(d);
    b_allele_col = apop_name_find(dc->names, "b_allele", 't');
    assert(!strcmp("T",  dc->text[3][b_allele_col]));
    rsid_col = apop_name_find(dc->names, "rsid", 't');
    assert(!strcmp("rs2977656",  dc->text[4][rsid_col]));
    assert(apop_data_get(dc, 5, .colname="ab")==201);

    apop_data_print(dc, "mixedtest", .output_type='d');
    apop_data *de = apop_query_to_mixed_data("mmmmtttt","select * from mixedtest");
    b_allele_col = apop_name_find(de->names, "b_allele", 't');
    assert(!strcmp("T",  de->text[3][b_allele_col]));
    rsid_col = apop_name_find(de->names, "rsid", 't');
    assert(!strcmp("rs2977656",  de->text[4][rsid_col]));
    assert(apop_data_get(de, 5, .colname="ab")==201);
    unlink("mixedtest");

    gsl_matrix *as_matrix = apop_query_to_matrix("select ab from d");
    gsl_vector *as_vector = apop_query_to_vector("select ab from d");
    Apop_matrix_col(as_matrix, 0, mv);
    gsl_vector_sub(as_vector, mv);
    Diff(apop_sum(as_vector), 0, 1e-10);

    test_uniform(d);
    apop_data_free(dc); apop_data_free(dd); 
    apop_data_free(dcc); apop_data_free(d); 
}

void test_blank_db_queries(){
    apop_db_close();
    apop_db_open(NULL);
    apop_table_exists("t", 'd');
    apop_query("create table t (a integer, b integer, c integer)");
    apop_data *d = apop_query_to_data("select * from t");
    apop_data *e = apop_query_to_text("select * from t");
    gsl_matrix *f = apop_query_to_matrix("select * from t");
    gsl_vector *g = apop_query_to_vector("select * from t");
    double h = apop_query_to_float("select * from t");
    assert(d==NULL);
    assert(e==NULL);
    assert(f==NULL);
    assert(g==NULL);
    assert(gsl_isnan(h));
}

void test_nan_data(){
    apop_table_exists("nandata", 'd');
    apop_table_exists("fw", 'd');
    apop_table_exists("fww", 'd');
    apop_text_to_db("test_data_nans", "nandata");
    strcpy(apop_opts.db_name_column, "head");
    apop_opts.nan_string = "nan";
    apop_data *d = apop_query_to_data("select * from nandata");
    apop_data_print(d, "nantest", .output_type='d');
    apop_data_free(d);
    apop_data *d2  = apop_query_to_data("select * from nantest");
    assert(gsl_isnan(apop_data_get(d2, .rowname="second", .colname="c")));
    assert(gsl_isnan(apop_data_get(d2, .rowname="third", .colname="b")));
    assert(!apop_data_get(d2, .rowname="fourth", .colname="b"));
    apop_data_free(d2);
    apop_opts.nan_string = "NaN";
    
    //while we're here, test querying just names & no data.
    apop_data *justnames = apop_query_to_data("select head from nandata");
    assert(justnames->names->rowct == 4);
    assert(!justnames->vector && !justnames->matrix);
    apop_data_free(justnames);

    //Oh, and let's test fixed-width inputs.
    apop_text_to_db("test_data_fixed_width", .tabname="fw", .has_col_names='n', .field_ends=(int[]){3,6});
    assert(apop_query_to_float("select col_2 from fw")==3.14159);
    apop_data *t=apop_query_to_text("select col_1 from fw");
    assert(!strcmp(*t->text[0], "A#C"));
    assert(!strcmp(*t->text[1], " BC"));
    apop_text_to_db("test_data_fixed_width", .tabname="fww", .field_names=(char*[]){"number", "text", "foat"}, .field_ends=(int[]){3,6});
    assert(apop_query_to_float("select number from fww where number<0")==-21);
    assert(apop_query_to_float("select foat from fww where text=' BC'")==2.71828);
    unlink("nantest");
}

#include <sys/wait.h> 
static void test_printing(){
    //This compares printed output to the printed output in the attached file. 
    char outfile[] = "print_test.out";

    if (!apop_table_exists("nandata"))
        test_nan_data();
    strcpy(apop_opts.db_name_column, "head");
    gsl_matrix *m  = apop_query_to_matrix("select * from nandata");
    apop_matrix_print(m, .output_name=outfile, .output_append='w');

apop_system("cp %s xxx", outfile);

    if (!apop_table_exists("d"))
        db_to_text();
    apop_data *d = apop_query_to_mixed_data ("tvttmmwt", "select * from d");
    FILE *f = fopen(outfile, "a");
    fprintf(f, "\nand a full vector+matrix+text+weights data set, formatted for computer reading:\n");
    strcpy(apop_opts.output_delimiter, "\t| ");
    apop_name_add(d->names, "Some SNPS", 'h');
    apop_data_print(d, .output_pipe =f);

    fprintf(f, "\nand just the names:\n");
    fclose(f);
    //need to redirect stdout.
    int status;
    if (fork() == 0){
        freopen(outfile, "a", stdout);
        apop_name_print(d->names);
        fclose(stdout);
        exit(0);
    }

    wait(&status);
    f = fopen(outfile, "a");

    fprintf(f, "\nand just the weights vector:\n");
    strcpy(apop_opts.output_delimiter, "\t");
    apop_vector_print(d->weights, .output_type='p', .output_pipe=f);
    fclose(f);
    int has_diffs = apop_system("diff -b printing_sample %s", outfile);
    assert(!has_diffs);
    //apop_system("rm %s", outfile);
    unlink(outfile);
    unlink("xxx");
}

void test_crosstabbing() {
    apop_db_close(); //gotta test it somewhere
    if (!apop_table_exists("snps"))
        apop_text_to_db("data-mixed", "snps", 0, 1);
    apop_table_exists("snp_ct", 'd');
    apop_query("create table snp_ct as "
                 " select a_allele, b_allele, count(*) as ct "
                 " from snps group by a_allele, b_allele ");
    apop_data *d = apop_db_to_crosstab("snp_ct", "a_allele", "b_allele", "ct");
    assert(apop_data_get(d, .rowname="A", "G")==5);
    assert(apop_data_get(d, .rowname="C", "G")==1);

    apop_data *ct = apop_text_alloc(apop_data_alloc(3,1),3,1);
    apop_data_set(ct, 0, 0, 1); apop_text_add(ct, 0, 0, "first");
    apop_data_set(ct, 1, 0, 2); apop_text_add(ct, 1, 0, "second");
    apop_data_set(ct, 2, 0, 3); apop_text_add(ct, 2, 0, "third");
    apop_table_exists("ct", 'd');
    apop_crosstab_to_db(ct, "ct", "r", "c", "val");
    if (apop_opts.db_engine=='s'){
        assert(!strcmp(**(apop_query_to_text("select val from ct where r='r0' and c='t0'")->text), "first"));
        assert(!strcmp(**(apop_query_to_text("select val from ct where r='r2' and c='t0'")->text), "third"));
    }
    assert(apop_query_to_float("select val from ct where r='r1' and c='c0'")==2);
}

#define do_test(text, fn) {if (verbose) printf("%s:", text); \
                          fflush(NULL);                      \
                          fn;                                \
                          if (verbose) printf("\nPASS.  ");} 

int main(){
    do_test("test data to db", test_data_to_db());
    do_test("db_to_text", db_to_text());
    do_test("test queries returning empty tables", test_blank_db_queries());
    do_test("NaN handling", test_nan_data());
    do_test("test printing", test_printing());
    do_test("test db to crosstab", test_crosstabbing());
    apop_db_close();
}
