/** \file cmd_apop_merge_dbs.c	A little command-line utility to merge
 two databases. Try <tt>apop_merge_dbs -h</tt> for help.*/

/* Copyright (c) 2005--2007, 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "db.h"
#include <unistd.h>

int main(int argc, char **argv){
int         merge_ct = 0;
char		c, 
		*delimiter,
		msg[1000], **merges = NULL;
	sprintf(msg, "%s [-v] [-t table_name] [-t next_tabname] main_db.db db_to_merge_into_main.db\n"
			     "   -t\ttable to merge. If none, do all in the source db. Use as many as you'd like.\n"
			     "   -v\tverbose\n", argv[0]); 
	if(argc<3){
		printf(msg);
		return 0;
	}
	delimiter	= malloc(5);
	strcpy(delimiter, ",");
	while ((c = getopt (argc, argv, "vht:")) != -1){
		switch (c){
		  case 'v':
			apop_opts.verbose	++;
			break;
		  case 'h':
			printf(msg);
			return 0;
          case 't':
            merges = realloc(merges, sizeof(char*)*merge_ct++);
            merges[merge_ct-1] = optarg;
		}
	}
	apop_db_open(argv[optind]);
    if (merge_ct)
        for (int i=0; i< merge_ct; i++)
            apop_db_merge_table(argv[optind +1], merges[i]);
    else
        apop_db_merge(argv[optind +1]);
	apop_db_close('n');
}
