#!/bin/sed -f 

/First, here is a comment./{
    s/This script filters C files to produce the headers for the compound \
    literal-based variadic function headers. For usage, your best bet is to learn by    \
    example and just compare files in the base and generated source codes.  \
    The workings are klunky: each line is saved in the hold-space (via \
    h) and re-pulled and re-processed for each line of output.  //
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
p
g
s/APOP_VAR_HEAD//
s/\(\[^ (]\)* *(/\1_base(/
s/\(.*\)[;{]/}\n\n\1{\n#endif/
h
d
}
/APOP_VAR_END_*HEAD/ {
g
}
