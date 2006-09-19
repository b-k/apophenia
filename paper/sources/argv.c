#include <stdio.h>

/* This program will simply print out the command line arguments given
 to it. */

int main(int argc, char **argv){
  int   i;
    for (i=0; i< argc; i++)
        printf("command line argument %i: %s\n", i, argv[i]);
    return 0;
}
