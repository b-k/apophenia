void del_one_row(gsl_matrix *in){
gsl_vector  v;
int         i;
//Allocate a matrix, get a reduced view of the original, and copy.
gsl_matrix  *reduced = gsl_matrix_alloc(in->size1 - 1, in->size2);
gsl_matrix  mv     = gsl_matrix_submatrix(in, 1,0, in->size1-1, in->size2).matrix;
    gsl_matrix_memcpy(reduced, &mv);
    
    //Do math here, or just display the matrix:
    apop_matrix_show(reduced);

    for(i=0; i< in->size1; i++){
        //Get a view of row i, and copy it to position i-1 in the
        //short matrix.
        v   = gsl_matrix_row(in, i).vector;
        gsl_matrix_set_row(reduced, i, &v);

        //Do math here, or just display the matrix:
        printf("\n");
        apop_matrix_show(reduced);
    }
}
