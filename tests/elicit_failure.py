#!/usr/bin/python
"""
test.c is the main test, but is intended to run through without any failure anywhere along
the way. This script is intented to test the steps that catch failures and warn the user
before a segfault happens. 

The array_of_failures is a Nx2 array. The first item in each row is a short program
designed to bring about a failure; the second item is the error message that should print
to screen. The remainder of the program generates these small test programs and checks
that their output is as promised.
"""


from subprocess import *

array_of_failures = [
["""apop_data *d1 = apop_data_alloc(10,9);
    apop_data *d2 = apop_data_alloc(10,9);
    apop_dot(d1, d2);"""
    , "apop_check_dimensions: You sent me a matrix with 9 columns to multiply against a matrix with 10 rows. Those two need to be equal."],
["""apop_data *d1 = apop_data_alloc(10,9);
    apop_data *d2 = apop_data_alloc(9,10);
    apop_dot(d1, d2, .form2='t');"""
    , "apop_check_dimensions: You sent me a matrix with 9 columns to multiply against a matrix with 10 rows (after the transposition you requested). Those two need to be equal."],
["""apop_data *d1 = apop_data_alloc(10,9);
    apop_data *d2 = apop_data_alloc(9,10);
    apop_dot(d1, d2, .form1='t');"""
    , "apop_check_dimensions: You sent me a matrix with 10 columns (after the transposition you requested) to multiply against a matrix with 9 rows. Those two need to be equal."],
["""apop_data *d1 = apop_data_alloc(10,9);
    apop_data *d2 = apop_data_alloc(10,9);
    apop_dot(d1, d2, .form1='t', .form2='t');"""
    , "apop_check_dimensions: You sent me a matrix with 10 columns to multiply against a matrix with 9 rows (after the two transpositions you requested). Those two need to be equal."],
['apop_data *d = apop_data_alloc(10,9); apop_matrix_inverse(d->matrix);'
    , "apop_det_and_inv: You asked me to invert a 10 X 9 matrix, but inversion requires a square matrix."],
['apop_data *d = apop_data_alloc(); apop_vector_increment(d->vector, 3, 3);'
    , "variadic_apop_vector_increment: You sent me a NULL vector."],
['apop_data *d = apop_data_alloc(); apop_matrix_increment(d->matrix, 3, 3);'
    , "variadic_apop_matrix_increment: You sent me a NULL matrix."],
["""
apop_data *d = apop_data_alloc(2,2);
Apop_matrix_row(d->matrix, 0, r);
apop_vector_realloc(r, 3);
    """, "apop_vector_realloc: I can't resize subvectors or other views."],
['apop_db_to_crosstab("faketab", "r1", "r2", "d");'
    , "apop_sqlite_query_to_text: no such table: faketab"],
['apop_query("create table faketab (r.1, r2, d)");' 
    , 'apop_query: near ".1": syntax error'],
['apop_query("create table faketab (r1, r2, d)"); apop_db_to_crosstab("faketab", "r1", "r2", "d");'
    , "apop_db_to_crosstab: selecting r1, r2, d from faketab returned an empty table."],
['apop_model null = {"A null model"}; apop_maximum_likelihood(NULL, &null);'
    , "setup_starting_point: The vector I'm trying to optimize over is NULL."],
['apop_model null = {"A null model",.vbase=2}; apop_maximum_likelihood(NULL, &null);'
    , "negshell: The model you sent to the MLE function has neither log_likelihood element nor p element."],
]

def onetest(i):
    compile_me = """
#include <apop.h>
int main(){
    """ + array_of_failures[i][0]+"""
    }
    """
    compile_file = file("t.c", "w")
    compile_file.write(compile_me)
    compile_file.close()

    testout = file("t.x", "w+")
    p= Popen("gcc t.c -o t -lapophenia -lsqlite3 -lgsl -lgslcblas --std=gnu99; ./t", shell=True, stderr=PIPE)

    outlines = p.stderr.readlines()
    ok = outlines[0].strip() == array_of_failures[i][1].strip()
    if (ok):
        print i, "OK."
    else:
        print "out = ", outlines[0] 
        print "s/b = ", array_of_failures[i][1]

if (__name__ == "__main__"):
    for i in range(0, len(array_of_failures)):
        onetest(i)
Popen("rm t t.c t.x", shell=True)
