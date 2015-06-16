s/Enumeration Type Documentation/Model Documentation/
s/<h2>Enumerations/<h2>Models/
/<h2>Models/,/<h2>Functions/{
    s/[{}]//
    s/<br\/>$//
    s/<li>enum/<li>apop_model/
    s/.*>(Overview|Name|Input_format|Estimate_results|Predict|RNG|CDF|Exampe|Settings)[+\]*_[+\]*x[0-9].*x[+\]*_//
}
s/<b>Enumerator: <\/b>//
s/.item\[Enumerator\].par//
s/^ enum//
s/\\setlength{\\rightskip}{0pt plus 5cm}enum/\\setlength{\\rightskip}{0pt plus 5cm}/
s|<tr><th colspan="2">Enumerator</th></tr>||
s/<td class="memname">enum <a class="el"/<td class="memname">apop_model <a class="el"/
s/[\+]*_[\+]*x[0-9]*x[_\+]*[_+]\+//g
s/model_specific/Methods are (D)efault<br> or (M)odel-specific/g
s/model[+\]*_[+\]*specific/\\hbox{(D)efault\/(M)odel-specific}/g
s/[eE]stimate[+\]*_[+\]*results/Post-estimate/g
s/[iI]nput[+\]*_[+\]*format/Input format/g
s/[pP]ostestimate[+\]*_[+\]*\(data\|parameters\|settings\_info\)/Post-estimate \1/g
s/[pP]ostestimate[+\]*_[+\]*parameter_model/Post-estimate parameter model/g
s/[pP]arameter[+\]*_[+\]*format/Parameter format/g
#delete all between the two markers, but not the second marker
/name="enum-members"/,/name="func-members"/{/name="func-members"/!d}
/<div class="summary">/,/\#func-members/d

s/Log likelihood &/LL \&/
s/Prep routine &/Prep \&/
