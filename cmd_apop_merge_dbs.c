/** \file cmd_apop_merge_dbs.c	A little command-line utility to merge
 two databases. Try <tt>apop_merge_dbs -h</tt> for help.

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "db.h"
#include <unistd.h>


int main(int argc, char **argv){
char		c, 
		*delimiter,
		msg[1000];
	sprintf(msg, "%s [opts] main_db.db db_to_merge_into_main.db\n\
			-v\tverbose\n", argv[0]); 

	if(argc<3){
		printf(msg);
		return 0;
	}
	delimiter	= malloc(5);
	strcpy(delimiter, ",");
	while ((c = getopt (argc, argv, "vh")) != -1){
		switch (c){
		  case 'v':
			apop_opts.verbose	++;
			break;
		  case 'h':
			printf(msg);
			return 0;
		}
	}
	apop_db_open(argv[optind]);
	apop_db_merge(argv[optind +1]);
	apop_db_close('n');
	return 0;
}
