#include "conversions.h"


int convert_vector_to_array(gsl_vector *in, double **out){
int		i;	
	*out	= malloc(sizeof(double) * in->size);
	for (i=0; i < in->size; i++)
		*out[i]	= gsl_vector_get(in, i);
	return in->size;
}

void convert_array_to_vector(double *in, gsl_vector **out, int size){
int		i;
	*out	= gsl_vector_alloc(size);
	for(i=0; i < size; i++)
		gsl_vector_set(*out, i, in[i]);
}



void convert_text_to_db(char *text_file, char *db, char *tabname, int ct, char **field_names, char *field_types){
FILE * 		infile;
char		q[20000], instr[10000], **fn, *astring;
int 		i, use_names_in_file	= 0;
	open_db(db);
	if (table_exists(tabname,0)) printf("%s table exists; not recreating it.\n", tabname);
	else{
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
		query_db(q);
		while(fgets(instr,1000,infile)!=NULL)
			if(instr[0]!='#') {
				sprintf(q, "INSERT INTO survey VALUES (%s);", instr);
				query_db(q);
			}
		query_db("commit;");
		fclose(infile);
		if (use_names_in_file){
			free(astring); 
			free(fn);
		}
	}
}
