apop_data	*districts, *data_set;
districts = apop_query_to_data(
        "select distinct age, gender, race, district from survey");

query   = malloc(sizeof(char)*10000);
strcpy(query,"select round(age/10) as age_group, gender, race, ");
for(i=1;i< districts->size1; i++){
        query=realloc(query, sizeof(char)*(strlen(query)+500));
        sprintf(query, 
                 "%s case district when %i then 1 else 0 end, \n",
                 query, (int) gsl_matrix_get(districts,i,0));
        }
query=realloc(query, sizeof(char)*(strlen(query)+5000));
strcat(query, "from survey;");
data_set = apop_query_to_data(query);
