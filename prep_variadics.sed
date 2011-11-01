#!/bin/sed -f 

/First, here is a comment./{
# This script filters C files to produce the headers for the compound
# literal-based variadic function headers. For usage, your best bet is to learn by 
# example and just compare files in the base and generated source codes.  
# The workings are klunky: each line is saved in the hold-space (via h)
# and re-pulled and re-processed for each line of output.  
# See also the technical notes at http://modelingwithdata.org/arch/00000022.htm
}


/APOP_VAR_DECLARE/ {
h
g
s/APOP_VAR_DECLARE/\#ifdef APOP_NO_VARIADIC\n/
s/!/,/g
p
g
s/APOP_VAR_DECLARE/\#else\n/
s/\([^ (]\) *(/\1_base(/
s/!/,/g
p
g
s/,/;/g
#annoying detail: if you take in a function, then those commas shouldn't be
#semicolons. so: declare like this: int (*infunction)(int! double *!  void)
#and I'll replace ! with , .
s/!/,/g
s/ *(/, /
s/ \([^ ]*,\)/, \1/
s/APOP_VAR_DECLARE / apop_varad_declare(/
s/!/,/g
p
g
#This form finds the line between function type and function name:
s/\([^* ]* *(\)/START_OF_FNAME\1/ 
s/.*START_OF_FNAME\([^(]*\)(.*/#define \1(...) apop_varad_link(\1, __VA_ARGS__)\n#endif/
s/[ \t]*(/(/
}

/APOP_VAR_HEAD/ {
h
g
s/APOP_VAR_HEAD/#ifdef APOP_NO_VARIADIC \n/
s/[;{][ \t]*$/{\n#else/
p
g
s/\([^* ]* *(\)/START_OF_FNAME\1/ 
s/START_OF_FNAME/, / 
s/(.*/){/
s/APOP_VAR_HEAD/apop_varad_head(/
}

#Here we construct the tail of the header function and the beginning of the base:

# return fn_base (v1, v2, v3);
# }
# fn_base(type1 v1, type2 *v2, type3 v3){
# #endif
 
/APOP_VAR_END_*HEAD/ {
g
s/APOP_VAR_HEAD//
t clear_branches #so above substitution won't matter
:clear_branches
#The C standard says that we can't "return void_fn()", so test whether the type is void.  
#(take care that it isn't void*, which does need a return)
s/^[ \t]*void[ \t]*\([^*]\)/\t\1/
t remove_types
s/^[ \t]*[^ *]*[ *]*/\treturn /

:remove_types
#comma or paren, text, spaces and stars, text2, comma or endparen 
#==> comma or paren, text2, comma or endparen
s/\([,(]\)[ \t]*[^ *,][^ *,]*[ *][ *]*\([^,][^,]*[,)]\)/\1 \2/
t remove_types

#a favor to the output system as currently written
s/Output_declares/Output_vars/

s/\(\[^ (]\)* *( */\1_base(/
s/[{ ]*$/;/
p
g
s/APOP_VAR_HEAD//
s/\(\[^ (]\)* *(/\1_base(/
s/\(.*\)[;{]/}\n\n\1{\n#endif/
h
g
}
