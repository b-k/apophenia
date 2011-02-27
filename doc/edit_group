s/Enumeration Type Documentation/Model Documentation/
s/<h2>Enumerations/<h2>Models/
/<h2>Models/,/<h2>Functions/{
    s/[{}]//
    s/<br\/>$//
    s/<li>enum/<li>apop_model/
    s/.*>(Overview|Name|Input_format|Estimate_results|Predict|RNG|CDF|Exampe|Settings)_x[0-9].*x_//
}
s/<b>Enumerator: <\/b>//
s/<td class="memname">enum <a class="el"/<td class="memname">apop_model <a class="el"/
s/_x[0-9]*x_//g
s/model_specific/Methods are Default or Model-specific/g
s/[eE]stimate_results/Post-estimate/g
s/[iI]nput_format/Input format/g
s/[pP]arameter_format/Parameter format/g
#delete all between the two markers, but not the second marker
/name="enum-members"/,/name="func-members"/{/name="func-members"/!d}
/<div class="summary">/,/\#func-members/d
