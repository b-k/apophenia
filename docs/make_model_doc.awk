#!/usr/bin/awk

# The goal: uniform model documentation for all the models that ship with Apophenia.
# 
# Doxygen provides pretty limited support for new types of documentation. It
# has the set of documentation types it does well, and that's that. So, this hack takes
# marks in the code and converts them into one documented enum block for each model,
# which Doxygen will then make beautiful.
# 
# We'll start with something like
# /* \amodel apop_beta The Beta distribution.
# 
#          The beta distribution has two parameters and is restricted between zero and
#          one. You may also find \ref apop_beta_from_mean_var to be useful.
# 
# \adoc    Input_format  Any arrangement of scalar values. 
# \adoc    Parameter_format   a vector, v[0]=\f$\alpha\f$; v[1]=\f$\beta\f$    
# \adoc    settings None. 
# */
# [...]
# /* \adoc Name  <tt>Beta distribution</tt>  */
# 
# and write out something like:
# 
# /** 
#  The Beta distribution.
# 
#          The beta distribution has two parameters and is restricted between zero and
#          one. You may also find \ref apop_beta_from_mean_var to be useful.
# 
#  \hideinitializer \ingroup models */
# enum apop_beta {
# Name_x2x_, /**<   <tt>Beta distribution</tt> */
# Input_format_x2x_, /**<   Any arrangement of scalar values. */
# Parameter_format_x2x_, /**<    a vector, v[0]=\f$\alpha\f$; v[1]=\f$\beta\f$    */
# RNG_x2x_, /**<   Produces a scalar \f$\in[0,1]\f$. */
# };
# 
# I added the _x#_x because doxygen merges enums together otherwise. There's a sed script
# (edit_group) that turns the enum documentation back into apop_model struct documentation.
# 
# The awk script can enforce some uniformity, forcing the same order in things, and writing
# file named missing, that will list documentation bits that aren't present.

BEGIN {
    IGNORECASE=1
    doc_parts[1]="Name";
    doc_parts[2]="Input_format"
    doc_parts[3]="Prep_routine"
    doc_parts[4]="Parameter_format"
    doc_parts[5]="Estimate_results"
    doc_parts[6]="Predict"
    doc_parts[7]="RNG"
    doc_parts[8]="CDF"
    doc_parts[9]="Settings"
    doc_parts[10]="Examples"
    doc_part_count=10
    model_count = 0
       #print > "missing"
}

in_doc==1 && !/\\a[model|doc]/ {
        if (sub("\\*/",""))
            in_doc=0
        items[current_model ":" current_item] = items[current_model ":" current_item] "\n" $0 
    }

/\\amodel/ {
    sub("/\\*[ \\t]*","", $0) #cut /* (if any).
    sub(".*\\\\amodel","", $0) #cut \amodel, now that I know what it is.
    current_model = $1
    in_doc=1
    oh = $0
    if (!models[current_model]) {
        models[current_model]=current_model
        model_count ++
    }
    sub(current_model,"", oh)
    sub("\\*/","", oh)
    sub("^[ \\t]","", oh)
    current_item = "intro"
    items[current_model ":" current_item] =  oh
    }

/\\adoc/ {
    sub("/\\*[ \\t]*","", $0) #cut /* (if any).
    sub("\\\\adoc[ \\t]*","", $0) #cut \adoc, now that I know what it is.
    oh = $0
    current_item = $1
    sub($1,"", oh)
    if (!sub("\\*/","", oh))
        in_doc = 1
    items[current_model ":" current_item] = oh
    }

/\*\// { in_doc = 0 }


#the declaration [like apop_model new_model = {"new model", .p=prob, .score=deriv}; ] tells us which struct elements are actually used.
/apop_model[ \t]*[^ \t]*[ \t]*=/ { in_decl=1 }

in_decl == 1 { cp = $0;
    if (match(cp, "\"([^\"]*)\"")){
        title= substr(cp,RSTART, RLENGTH)
        gsub( "\"","",title)
        items[current_model ":" "Name"] = "<tt>" title "</tt>" 
    }

    while (match(cp, "\\.[^ \\t=]+[ \\t]*=")){
        cp = substr(cp, RSTART+1, length(cp)-RSTART)
        match(cp, "[ \\t=]")
        items[current_model ":has" substr(cp, 0, RSTART-1)] = 1
        cp = substr(cp, RSTART+1) #cut at the =, keep processing.
    }
}

/\\}/ {in_decl=0;}
/}/ {in_decl=0;}

function onedot(name, isdefault){
    dot = "<td class=\"altcolor\">M</td>"
    defaultdot = "<td class=\"memitem\">D</td>"
    if (isdefault)
        class = "altcolor"
    else 
        class = "memitem"
    print "<tr><td class=\""class"\">" name "</td>"; 
    if (isdefault) 
        print dot;
    else
        print defaultdot;
    print "</tr>"
}

END {print "/** \\file */ /**\\defgroup models */"
    model_no=0
    for (m in models){
        model_no++
        #m = sorted_models[model_no];
        print "\n/** "
        if (items[m ":intro"])
            print items[m ":intro"]
        else
            print "!!! Without an intro, " m " won't print\n" >> "/dev/stderr"
        print " \\hideinitializer \\ingroup models */"
        print "enum " m " {"

#apop_model apop_normal = {"Normal distribution", 2, 0, 0, .dsize=1, 
# .estimate = normal_estimate, .log_likelihood = normal_log_likelihood, 
# .score = normal_dlog_likelihood,    
# .constraint = beta_1_greater_than_x_constraint, .draw = normal_rng,
# .cdf = normal_cdf, .predict = normal_predict}; 

            #dot = "<td class=\"memitem\" >\\f$\\bullet\\f$</td>"
            print "model_specific_x" model_no "x_, /**< <table cellpadding=3px><tr><td><table class=\"memproto\">"

            onedot("Estimation", items[m ":hasestimate"])
            onedot("Prob.", items[m ":hasp"])
            onedot("Log likelihood", items[m ":haslog_likelihood"])
            onedot("RNG", items[m ":hasdraw"])

            print "</table></td><td><table class=\"memproto\">"

            onedot("Predict", items[m ":haspredict"])
            onedot("CDF", items[m ":hascdf"])
            onedot("Score", items[m ":hasscore"])
            onedot("Prep routine", items[m ":hasprep"])

            print "</table> </td></tr></table>*/"
        for (i=1; i<=doc_part_count;i++){
            part=doc_parts[i]
            #print "processing", m, part, "\n" >> "/dev/stderr"
            if (doc_parts[i]=="Estimate_results"){
                print part "_x" model_no "x_, /**< <table>"
                print "<tr><td style=\"vertical-align:top\"><tt>data</tt></td><td style=\"vertical-align:top\">"
                if (items[m ":estimated_data"])
                    print items[m ":estimated_data"] 
                else
                    print "Unchanged."
                print "</td></tr>"
                print "<tr><td style=\"vertical-align:top\"><tt>parameters</tt></td><td style=\"vertical-align:top\">"
                if (items[m ":estimated_parameters"])
                    print items[m ":estimated_parameters"] 
                else
                    print "See parameter format."
                print "</td></tr>"
                if (items[m ":estimated_parameter_model"]){
                    print "<tr><td style=\"vertical-align:top\"><tt>parameter models</tt></td><td style=\"vertical-align:top\">"
                    print items[m ":estimated_parameter_model"] "</td></tr>"
                }
                if (items[m ":estimated_info"]){
                    print "<tr><td style=\"vertical-align:top\"><tt>info</tt></td><td style=\"vertical-align:top\">"
                    print items[m ":estimated_info"] "</td></tr>"
                }
                if (items[m ":estimated_settings"]){
                    print "<tr><td style=\"vertical-align:top\">settings</td><td style=\"vertical-align:top\">"
                    print items[m ":estimated_settings"] "</td></tr>"
                }
                print "</table> */"
            }
            else if (items[m ":" part]) print part "_x" model_no "x_, /**< " items[m ":" part] "*/"
            else print m, part >> "missing_model_parts"  #not at the moment important.
        }
        print "};"
    }
}
