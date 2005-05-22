//conversions.c  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
#include "conversions.h"


int apop_convert_vector_to_array(gsl_vector *in, double **out){
int		i;	
	*out	= malloc(sizeof(double) * in->size);
	for (i=0; i < in->size; i++)
		*out[i]	= gsl_vector_get(in, i);
	return in->size;
}

void apop_convert_array_to_vector(double *in, gsl_vector **out, int size){
int		i;
	*out	= gsl_vector_alloc(size);
	for(i=0; i < size; i++)
		gsl_vector_set(*out, i, in[i]);
}

void convert_array_to_line(double **in, double **out, int rows, int cols){
	//go from in[i][j] form to the GSL's preferred out[i*cols + j] form
int		i, j;
	*out	= malloc(sizeof(double) * rows * cols);
	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			(*out)[i * cols + j]	= in[i][j];
}

void apop_convert_array_to_matrix(double **in, gsl_matrix **out, int rows, int cols){
gsl_matrix_view	m;
double		*line;
	*out	= gsl_matrix_alloc(rows, cols);
	convert_array_to_line(in, &line, rows, cols);
	m	= gsl_matrix_view_array(line, rows,cols);
	gsl_matrix_memcpy(*out,&(m.matrix));
	free(line);
}

int find_index(gsl_matrix *d, double r, int start_from){
//used for apop_db_to_crosstab.
int	i	= start_from;	//i is probably the same or i+1.
	while(i!=(start_from-1)){
		if(gsl_matrix_get(d, i,0) == r) 
			return i;
		i	++;
		if (i == d->size1)
			i = 0;
	}
	printf("something went wrong in the crosstabbing.\n");
	return 0;
}


gsl_matrix * apop_db_to_crosstab(char *tabname, char *r1, char *r2, char *datacol, gsl_vector **d1, gsl_vector **d2){
//Give the name of a table in the database, and names of three of its
//columns: the x-dimension, the y-dimension, and the data.
//the output is a 2D matrix with rows indexed by r1 and cols by
//r2. if !=NULL, d1 and d2 will list the labels on the dimensions.

gsl_matrix	*pre_d1, *pre_d2, *datatab, *out;
int		i	= 0,
		j	= 0,
		k, r,c,datum;
gsl_vector_view	v;
	pre_d1	= apop_query_to_matrix("select distinct %s from %s order by %s", r1, tabname, r1);
	if (pre_d1 == NULL) printf ("Selecting %s from %s returned an empty table.\n", r1, tabname);
	pre_d2	= apop_query_to_matrix("select distinct %s from %s order by %s", r2, tabname, r2);
	if (pre_d2 == NULL) printf ("Selecting %s from %s returned an empty table.\n", r2, tabname);
	datatab	= apop_query_to_matrix("select %s, %s, %s from %s", r1, r2, datacol, tabname);
	out	= gsl_matrix_calloc(pre_d1->size1, pre_d2->size1);
	for(k =0; k< datatab->size1; k++){
		r	= gsl_matrix_get(datatab, k, 0);
		c	= gsl_matrix_get(datatab, k, 1);
		datum	= gsl_matrix_get(datatab, k, 2);
		i	= find_index(pre_d1, r, i);
		j	= find_index(pre_d2, c, j);
		gsl_matrix_set(out, i, j, datum);
	}
	if(d1!=NULL && d2!= NULL){
		*d1	= gsl_vector_alloc(pre_d1->size1);
		*d2	= gsl_vector_alloc(pre_d2->size1);
		v	= gsl_matrix_column(pre_d1, 0);
		gsl_vector_memcpy(&(v.vector), *d1);
		v	= gsl_matrix_column(pre_d2, 0);
		gsl_vector_memcpy(&(v.vector), *d2);
	}
	free(pre_d1); free(pre_d2); free(datatab);
	return out;
}


int count_cols_in_text(char *text_file){
//Open file, find the first valid row, count columns, close file.
FILE * 		infile;
char		instr[100000], *astring;
int		ct	= 0;
	infile	= fopen(text_file,"r");
	fgets(instr, 100000, infile);
	while(instr[0]=='#')	//burn off comment files
		fgets(instr, 10000, infile);
	astring	= strtok(instr,",");
	while (astring !=NULL){
		ct++;
		astring	= strtok(NULL,",");
	}
	return ct;
}


int apop_convert_text_to_array(char *text_file, double ***tab, int has_field_names){
FILE * 		infile;
char		instr[100000], *astring;
int 		i	= 0,
		ct, colno;
	ct	= count_cols_in_text(text_file);
	*tab	= malloc(sizeof(double));
	infile	= fopen(text_file,"r");
	if (has_field_names == 1){
		fgets(instr, 100000, infile);
		while(instr[0]=='#')	//burn off comment files
			fgets(instr, 10000, infile);
	}
	while(fgets(instr,10000,infile)!=NULL){
		colno	= 0;
		if(instr[0]!='#') {
			i	++;
			(*tab)	= realloc(*tab, sizeof(double) * i);
			(*tab)[i-1]= malloc(sizeof(double)*ct);
			astring	= strtok(instr,",");
			while (astring !=NULL){
				colno++;
				(*tab)[i-1][colno-1]	= atof(astring);
				astring	= strtok(NULL,",");
			}
		}
	}
	fclose(infile);
	return i;
}


int apop_convert_text_to_db(char *text_file, char *tabname, char **field_names){
FILE * 		infile;
char		q[20000], instr[100000], **fn, *astring;
int 		ct,
		i			= 0, 
		use_names_in_file	= 0,
		rows			= 0;
	ct	= count_cols_in_text(text_file);
	if (apop_table_exists(tabname,0)){
	       	printf("%s table exists; not recreating it.\n", tabname);
		return 0; //to do: return the length of the table.
	} else{
		infile	= fopen(text_file,"r");
		if (field_names == NULL){
			use_names_in_file++;
			fgets(instr, 100000, infile);
			while(instr[0]=='#')	//burn off comment files
				fgets(instr, 10000, infile);
			fn	= malloc(ct * sizeof(char*));
			astring	= strtok(instr,",");
			while(astring !=NULL){
				fn[i]	= malloc(1000 * sizeof(char));
				strcpy(fn[i], astring);
				astring	= strtok(NULL,",");
				i++;
			}
		} else	fn	= field_names;
		strcpy(q, "begin; CREATE TABLE ");
		strcat(q, tabname);
		for (i=0; i<ct; i++){
			if (i==0) 	strcat(q, " (");
			else		strcat(q, " , ");
			sprintf(q, "%s %s", q, fn[i]);
		}
		strcat(q, "); commit; begin;");
		apop_query_db(q);
		while(fgets(instr,10000,infile)!=NULL){
			rows	++;
			if(instr[0]!='#') {
				sprintf(q, "INSERT INTO %s VALUES (%s);", tabname, instr);
				apop_query_db(q);
			}
		}
		apop_query_db("commit;");
		fclose(infile);
		if (use_names_in_file){
			free(astring); 
			free(fn);
		}
		return rows;
	}
}

