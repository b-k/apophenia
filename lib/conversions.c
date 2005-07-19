//conversions.c  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
#include <apophenia/conversions.h>

#define Text_Size_Limit 1000000


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

int apop_find_index(gsl_matrix *d, double r, int start_from){
//used for apop_db_to_crosstab.
int	i	= start_from;	//i is probably the same or i+1.
	do {
		if(gsl_matrix_get(d, i,0) == r) 
			return i;
		i	++;
		i	%= d->size1;	//loop around as necessary.
		//if (i == d->size1)
			//i = 0;
	} while(i!=start_from); 
	printf(" apop %s, %i: something went wrong in the crosstabbing; couldn't find %g.\n", __FILE__, __LINE__, r);
	return 0;
}

gsl_matrix * apop_db_to_crosstab(char *tabname, char *r1, char *r2, char *datacol, gsl_vector **d1, gsl_vector **d2){
//Give the name of a table in the database, and names of three of its
//columns: the x-dimension, the y-dimension, and the data.
//the output is a 2D matrix with rows indexed by r1 and cols by
//r2. if !=NULL, d1 and d2 will list the labels on the dimensions.

gsl_matrix	*pre_d1	= NULL, 
		*pre_d2	= NULL, *datatab, *out;
int		i	= 0,
		j	= 0,
		k; 
double		datum, r, c;
gsl_vector_view	v;
	pre_d1	= apop_query_to_matrix("select distinct %s, 1 from %s order by %s", r1, tabname, r1);
	if (pre_d1 == NULL) 
		printf (" apop %s, %i: selecting %s from %s returned an empty table.\n", __FILE__, __LINE__, r1, tabname);
	pre_d2	= apop_query_to_matrix("select distinct %s from %s order by %s", r2, tabname, r2);
	if (pre_d2 == NULL) 
		printf (" apop %s, %i: selecting %s from %s returned an empty table.\n", __FILE__, __LINE__, r2, tabname);
	datatab	= apop_query_to_matrix("select %s, %s, %s from %s", r1, r2, datacol, tabname);
	out	= gsl_matrix_calloc(pre_d1->size1, pre_d2->size1);
	for (k =0; k< datatab->size1; k++){
		r	= gsl_matrix_get(datatab, k, 0);
		c	= gsl_matrix_get(datatab, k, 1);
		datum	= gsl_matrix_get(datatab, k, 2);
		i	= apop_find_index(pre_d1, r, i);
		j	= apop_find_index(pre_d2, c, j);
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
	gsl_matrix_free(pre_d1); gsl_matrix_free(pre_d2); gsl_matrix_free(datatab);
	return out;
}



int apop_count_cols_in_text(char *text_file){
//Open file, find the first valid row, count columns, close file.
FILE * 		infile;
char		instr[Text_Size_Limit], *astring;
int		ct	= 0;
	infile	= fopen(text_file,"r");
	if (infile==NULL) {printf("Error opening %s", text_file); return 1;}
	fgets(instr, Text_Size_Limit, infile);
	while(instr[0]=='#')	//burn off comment lines
		fgets(instr, Text_Size_Limit, infile);
	astring	= strtok(instr,",");
	while (astring !=NULL){
		ct++;
		astring	= strtok(NULL,",");
	}
	fclose(infile);
	return ct;
}

char * strip(char *in){
//OK, OK. C sucks.
char 		*out 	= malloc(sizeof(char) * (1+strlen(in)));
int		i	= 0, 
		done	= 0,
		first_ok= 0,
		last_ok	= 0;
	while (!done)
		if (in[first_ok]==' ' || in[first_ok]=='\n' || in[first_ok]=='\t')
			first_ok	++;
		else	done		++;
	for (i=first_ok  ; i< strlen(in); i++){
		out[i-first_ok]	= in[i];
		if (in[i]!=' ' && in[i]!='\n' && in[i]!='\t')
			last_ok	= i-first_ok;
	}
	out[last_ok+1]	= '\0';
	return out;
}

void add_a_name(int *ct, char*** list, char *addme){
char		*stripped;
	(*ct)	++;
	*list	= realloc(*list, (*ct) * sizeof(char*));
	(*list)[(*ct)-1]	=malloc(sizeof(char) * (strlen(addme)+1));
	stripped	= strip(addme);
	strcpy((*list)[(*ct)-1],stripped);
	free(stripped);
}


int apop_convert_text_to_array(char *text_file, char *delimiters, double ***tab, apop_name *names){
FILE * 		infile;
char		instr[Text_Size_Limit], *astring, *str;
int 		i	= 0,
		ct, colno;
	ct	= apop_count_cols_in_text(text_file)+1;
	*tab	= malloc(sizeof(double));
	infile	= fopen(text_file,"r");
	if (names != NULL){
		fgets(instr, Text_Size_Limit, infile);
		while(instr[0]=='#')	//burn off comment lines
			fgets(instr, Text_Size_Limit, infile);
		if (names->colnamect !=0){
			astring	= strtok(instr,delimiters);
			names->colnamect= 0;
			names->colnames	= malloc(sizeof(char*));
			while (astring !=NULL){
				add_a_name(&(names->colnamect), &(names->colnames), astring);
				astring	= strtok(NULL,delimiters);
			}
		}
		if (names->rownamect !=0){
			names->rownamect= 0;
			names->rownames	= malloc(sizeof(char*));
		} else 	names->rownames	= NULL;
	}
	while(fgets(instr,Text_Size_Limit,infile)!=NULL){
		colno	= 0;
		if(instr[0]!='#') {
			i	++;
			*tab	= realloc(*tab, sizeof(double*) * i);
			(*tab)[i-1]= malloc(sizeof(double)*ct);
			astring	= strtok(instr,delimiters);
				if (names !=NULL && names->rownames !=NULL){
					add_a_name(&(names->rownamect), &(names->rownames), astring);
					astring	= strtok(NULL,delimiters);
					if (astring==NULL){
						printf("row name with no data on line %i.\n", i);
						return 1;
					}
				}
			while (astring !=NULL){
				colno++;
				(*tab)[i-1][colno-1]	= strtod(astring, &str);
				if (!strcmp(astring, str))
					printf("trouble converting item %i on line %i; using zero.\n", colno, i);
				astring	= strtok(NULL,delimiters);
			}
		}
	}
	fclose(infile);
	return i;
}

char * prep_string_for_sqlite(char *astring){
//If the string isn't a number, it needs quotes.
char		*tmpstring, 
		*out	= NULL,
		*str	= NULL;
	strtod(astring, &str);
	if (!strcmp(astring, str)){	//then it's not a number.
		tmpstring=strip(astring);
		if (tmpstring[0]!='"'){
			out	= malloc(sizeof(char) * (strlen(tmpstring)+3));
			sprintf(out, "\"%s\"",tmpstring);
		}
		free(tmpstring);
		return out;
	} else {			//sqlite wants 0.1, not .1
		tmpstring=strip(astring);
		if (tmpstring[0]=='.'){
			out	= malloc(sizeof(char) * (strlen(tmpstring)+2));
			sprintf(out, "0%s",tmpstring);
			free(tmpstring);
			return out;
		}
		free(tmpstring);
	} //If you're here, then it's a number which needs no fixing.
	out	= malloc(sizeof(char) * (strlen(astring)+1));
	strcpy(out, astring);
	return out;
}

int apop_convert_text_to_db(char *text_file, char *tabname, char **field_names){
FILE * 		infile;
char		q[20000], instr[Text_Size_Limit], **fn, *astring;
char		delimiters[]	=",";
int 		ct, one_in,
		i			= 0, 
		use_names_in_file	= 0,
		rows			= 0;
	ct	= apop_count_cols_in_text(text_file);
	if (apop_table_exists(tabname,0)){
	       	printf("apop: %s table exists; not recreating it.\n", tabname);
		return 0; //to do: return the length of the table.
	} else{
		infile	= fopen(text_file,"r");
	       	if (infile==NULL) {
			printf("apop, %s: %i. Trouble opening %s.\n", __FILE__, __LINE__, text_file);
			return 0;
		}
		if (field_names == NULL){
			use_names_in_file++;
			fgets(instr, Text_Size_Limit, infile);
			while(instr[0]=='#')	//burn off comment files
				fgets(instr, Text_Size_Limit, infile);
			fn	= malloc(ct * sizeof(char*));
			astring	= strtok(instr,delimiters);
			while(astring !=NULL){
				fn[i]	= malloc(1000 * sizeof(char));
				strcpy(fn[i], astring);
				astring	= strtok(NULL,delimiters);
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
		while(fgets(instr,Text_Size_Limit,infile)!=NULL){
			if((instr[0]!='#') && (instr[0]!='\n')) {	//comments and blank lines.
				rows	++;
				i++;
				one_in=0;
				sprintf(q, "INSERT INTO %s VALUES (", tabname);
				astring	= strtok(instr,delimiters);
				while(astring !=NULL){
					if(one_in++) 	strcat(q, ", ");
					strcat(q, prep_string_for_sqlite(astring));
					astring	= strtok(NULL,delimiters);
				}
				apop_query_db("%s);",q);
			}
		}
		apop_query_db("commit;");
		fclose(infile);
		if (use_names_in_file){
			free(fn);
		}
		return rows;
	}
}

int apop_text_to_db(char *text_file, char *tabname, char **field_names){
	//just an alias.
	return apop_convert_text_to_db(text_file, tabname, field_names);
}

int apop_crosstab_to_db(gsl_matrix *in, apop_name n, char *tabname, char *row_col_name, 
						char *col_col_name, char *data_col_name){
int		i,j;
	apop_query_db("CREATE TABLE %s (%s , %s , %s);", tabname, row_col_name, col_col_name, data_col_name);
	for (i=0; i< n.colnamect; i++){
		apop_query_db("begin;");
		for (j=0; j< n.rownamect; j++)
			apop_query_db("INSERT INTO %s VALUES (%s, %s,%g);", tabname, 
					n.rownames[j], n.colnames[i], gsl_matrix_get(in, j, i));
		apop_query_db("commit;");
	}
	return 0;
}
