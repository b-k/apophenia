#include "apophenia/db.h"
#include <unistd.h>


int main(int argc, char **argv){
char		c, 
		*delimiter,
		msg[1000],
		*outfile	=NULL;
gsl_matrix	*m;

	sprintf(msg, "%s [opts] main_db.db db_to_merge_into_main.db\n", argv[0]); 

	if(argc<3){
		printf(msg);
		return 0;
	}
	delimiter	= malloc(sizeof(char) * 5);
	strcpy(delimiter, ",");
	while ((c = getopt (argc, argv, "d:f:h")) != -1){
		switch (c){
		  case 'h':
			printf(msg);
			return 0;
		}
	}
	apop_open_db(argv[optind]);
	m	= apop_merge_db(argv[optind +1]);
	apop_close_db(0);
	return 0;
}
