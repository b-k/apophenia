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
apop_print_matrix(*out);
	free(line);
}



int apop_convert_text_to_array(char *text_file, double ***tab, int ct, int has_field_names){
FILE * 		infile;
char		instr[10000], *astring;
int 		i	= 0,
		colno;
	*tab	= malloc(sizeof(double));
	infile	= fopen(text_file,"r");
	if (has_field_names == 1){
		fgets(instr, 10000, infile);
		while(instr[0]=='#')	//burn off comment files
			fgets(instr, 10000, infile);
	}
	while(fgets(instr,1000,infile)!=NULL){
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


int apop_convert_text_to_db(char *text_file, char *tabname, int ct, char **field_names, char *field_types){
FILE * 		infile;
char		q[20000], instr[10000], **fn, *astring;
int 		i, 
		use_names_in_file	= 0,
		rows			= 0;
	if (apop_table_exists(tabname,0)){
	       	printf("%s table exists; not recreating it.\n", tabname);
		return 0; //to do: return the length of the table.
	} else{
		infile	= fopen(text_file,"r");
		if (field_names == NULL){
			use_names_in_file++;
			fgets(instr, 10000, infile);
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
			if (field_types[i] =='I')
				strcat(q, " INTEGER" );
		}
		strcat(q, "); commit; begin;");
		apop_query_db(q);
		while(fgets(instr,1000,infile)!=NULL){
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
