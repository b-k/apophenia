gsl_matrix	*districts, *data_set;
query_to_matrix(&districts,
        "select distinct district from survey \
        where   age!=-1 and gender !=-1 and race !=-1 and \
        year !=-1 and district !=-1 and \
        (votedverified ==1 or votedverified ==3)");

query=malloc(sizeof(char)*10000);
strcpy(query,"select \ 
        case votedverified when 3 then 1 else 0 end,  \
        round(age/10) as age_group, gender, ");
for(i=1;i< districts->size1; i++){
        query=realloc(query, sizeof(char)*(strlen(query)+500));
        sprintf(query, 
                 "%s case district when %i then 1 else 0 end, \n",
                 query, (int) gsl_matrix_get(districts,i,0));
        }
query=realloc(query, sizeof(char)*(strlen(query)+5000));
strcat(query, " year-1950 \
        from survey \
        where   age!=-1 and gender !=-1 and race !=-1 and \
        year !=-1 and district !=-1 and \
        (votedverified ==1 or votedverified ==3)");
query_to_matrix(&data_set, query);
