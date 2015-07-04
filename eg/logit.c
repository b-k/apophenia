// See http://modelingwithdata.org/arch/00000160.htm for context and analysis.

#ifdef Datadir
#define DATADIR Datadir
#else
#define DATADIR "."
#endif

#include <apop.h>

int main(){
    //read the data to db, get the desired columns,
    //prep the two categorical variables
    apop_text_to_db( DATADIR "/" "amash_vote_analysis.csv" , .tabname="amash");
    apop_data *d = apop_query_to_mixed_data("mmmtt", "select 0, ideology,log(contribs+10) as contribs, vote, party from amash");
    apop_data_to_factors(d); //0th text col -> 0th matrix col
    apop_data_to_dummies(d, .col=1, .type='t', .append='y');

    //Estimate a logit model, get covariances,
    //calculate p values under popular Normality assumptions
    Apop_model_add_group(apop_logit, apop_parts_wanted, .covariance='y');
    apop_model *out = apop_estimate(d, apop_logit);
    apop_model_show(out);
    for (int i=0; i< out->parameters->matrix->size1; i++){
        printf("%s pval:\t%g\n",out->parameters->names->row[i], 
        apop_test(apop_data_get(out->parameters, i), "normal", 0, sqrt(apop_data_get(out->parameters->more, i, i))));
    }
}
