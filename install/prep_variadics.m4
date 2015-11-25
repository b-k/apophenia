m4_divert(-1)

These are the macros to filter C files to produce the headers for the compound
literal-based variadic function headers. For usage, your best bet is to learn by 
example and compare files in the base and generated source codes.  

Most of the work is in splitting the input down into the type name, function name,
function arguments, and function body. Then, the parts get reassembled as per the
template.

See also the technical notes at http://modelingwithdata.org/arch/00000022.htm

m4_changequote(`<|', `|>')
m4_changecom()
m4_define(APOP_VAR_HEAD, <|Variadify(|>)
m4_define(APOP_VAR_ENDHEAD, <|)|>)
m4_define(APOP_VAR_END_HEAD, <|)|>)
m4_define(cutm4, <|m4_define(<|stlenForCut|>, m4_regexp(<|$1|>, <|$2|>))m4_substr(<|$1|>, 0, stlenForCut)|>)
m4_define(postCutm4, <|m4_define(<|stlenForCut|>, m4_regexp(<|$1|>, <|$2|>))m4_substr(<|$1|>, stlenForCut)|>)
m4_define(Names_Only, <|m4_patsubst(m4_patsubst($1, <| *$|>,), <|.*[ *]|>,)<||>m4_ifelse(<|$#|>, <|0|>, ,<|$#|>, <|1|>, , <|, Names_Only(m4_shift($@))|>)|>)

m4_define(Variadify, <|m4_dnl
m4_define(<|PreParen|>, cutm4(<|$*|>, <| *(|>,))m4_dnl
m4_define(<|Fn_Name|>, m4_patsubst(PreParen, <|^.*\( \|\*\)|>, ))m4_dnl
m4_define(<|Type_Name|>, cutm4(cutm4(PreParen, Fn_Name), <| *$|>))m4_dnl
m4_define(<|PreBrace|>, <|cutm4(<|$*|>,<| *{|>)|>)m4_dnl
m4_define(<|PostBrace|>, <|postCutm4(<|$*|>,<| *{|>)|>)m4_dnl
m4_define(<|Args_m4|>, postCutm4(PreBrace, <|(|>))m4_dnl
m4_define(<|AArgs_m4|>, m4_patsubst(Args_m4,Output_declares,Output_vars))m4_dnl A favor to apop_output.c
#ifdef APOP_NO_VARIADIC
Type_Name Fn_Name<||>Args_m4{
#else
apop_varad_head(Type_Name, Fn_Name)PostBrace    m4_ifelse(Type_Name, <|void|>,,return) Fn_Name<||>_base(Names_Only(m4_patsubst(AArgs_m4,<| *[()]|>,)));
}

 Type_Name Fn_Name<||>_base<||>Args_m4{
#endif|>)

m4_define(Apop_var_declare,<|m4_dnl
m4_define(<|PreParen|>, cutm4(<|$*|>, <| *(|>,))m4_dnl
m4_define(<|Fn_Name|>, m4_patsubst(PreParen, <|^.*\( \|\*\)|>, ))m4_dnl
m4_define(<|Type_Name|>, cutm4(cutm4(PreParen, Fn_Name), <| *$|>))m4_dnl
m4_define(<|Args_m4|>, postCutm4(<|$*|>, <|(|>))m4_dnl
m4_define(<|Args_m4_semicolonized|>, m4_patsubst(m4_translit(Args_m4, <|,|>, ;), <|\(^ *(\|) *$\)|>, ))m4_dnl
m4_patsubst(m4_dnl
#ifdef APOP_NO_VARIADIC
 $*;
#else
 Type_Name Fn_Name<||>_base<||>Args_m4;
 apop_varad_declare(Type_Name, Fn_Name, Args_m4_semicolonized);
#define Fn_Name<||>(...) apop_varad_link(Fn_Name, __VA_ARGS__)
#endif,
<|!|>, <|,|>)
|>)

m4_define_at_some_point(<|m4_apop_version|>, <|m4_syscmd(date +%Y%m%d)|>)
m4_define(<|m4_apop_version|>, <|1.0|>)
m4_divert(0)
