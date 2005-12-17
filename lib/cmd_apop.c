/** \file cmd_apop.c	The server/client which allows users to call c_functions from the command
  line/script 
 Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
 */
#include <apophenia/headers.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>	//read,write
#include <sys/socket.h>
#include <sys/un.h>
#include <sys/types.h> 
#include <netinet/in.h>
#include <netdb.h> 

#define PORT_NO 11716

int find_matrix(char *name);
char *print_matrix(gsl_matrix *m);
typedef struct operation{
	char 	name[100]; //if your command is longer than this, it's too much.
	int	argc;
	char 	internal_return;
	int 	external_return;
	char* 	(*run_me)(char **argbuffer);
	char 	help_text[1000];
} operation;


/** \page command_line Command line utilities

<b> The command-line database utilities</b><br>
Try \ref apop_db_to_crosstab, \ref apop_text_to_db, and
\ref apop_merge_dbs from the command line (<tt>-h</tt> for help). They are
simple stand-alone wrappers to the corresponding functions, in case you
just need to quickly move data in or out of a database.

<b> The server/client </b><br>
But more importantly, there is the command-line program <tt>apop</tt>,
which facilitates stats from outside of C, such as on the command line
or via a scripting language. First,

\code
apop start data.db&
\endcode

will open a server, which will in this case open a database named
<tt>data.db</tt>. The ampersand is the standard UNIX means of running a
process in background. Run <tt>apop start</tt> with no database name to
use an in-memory database.

Once the server is listening, you may send it commands. Eventually, every
function in the Apophenia library and the more useful ones from the GSL
will have a corresponding command line. At the moment, the functions
are basically those for a proof-of-concept.
Type <tt>apop help</tt> for the current list of supported functions. Below are a sample.

To query the database:
\code
apop query "select something from table"
\endcode
[Why not just use the sqlite command line? First, since the apophenia
server side persists, the database isn't loaded into memory for every
call. Second, Apophenia includes new aggregate operators---notably
<tt>var()</tt>. Finally, apophenia's server makes it possible to use an
in-memory database.]


To pull data to a matrix, which is then stored on the server:
\code
apop query_to_matrix matrix_name "select columns from db_table"
\endcode

To print the matrix to the screen:
\code
apop matrix_print matrix_name
\endcode

To do a t-test on two columns in the database:
\code
apop db_t_test tab1 col1 tab2 col2
\endcode

To stop the server:
\code
apop stop
\endcode

\section scripting_languages How to use Apophenia from Perl/Python/&c.
The easiest way is to simply call the command-line interface. In Perl, try <tt>system("apop start")</tt>; in Python, <tt>include os</tt> and then you can issue commands like <tt>os.popen3("""apop query "select * from tab" """)</tt>.

Alternatively, the reader is encouraged to write and contribute
functions to add a layer of abstraction to the command line, so users
may type \ref apop_query <tt>("select")</tt> and the \ref apop_query function
would call the appopriate command line. Since the server does all the
math, the scripting functions will only have to handle text processing
and reading matrices into native scripting-language arrays.
*/


gsl_matrix **	matrix_list;
char **		matrix_names;
int		matrix_ct	= 0;

/*/////////////
//The functions 

How to add a fn: You'll need to add a handler function here, using the same prototype as all the other handler fns.
Then, add a line in the command_list below, including the appropriate reference to your handling function. 
That's all.
////////////// */
char* null_fn(char **argbuffer) {return NULL;}

char* do_db_merge_table(char **argbuffer) {
	     apop_db_merge_table(argbuffer[0], argbuffer[1]);
	     return NULL;
	     }

char* do_db_t_test(char **argbuffer) {
char * 		printme	= malloc(sizeof(char)*2000);
double		out;
	out	= apop_db_t_test(argbuffer[0],argbuffer[1],argbuffer[2],argbuffer[3]);
	sprintf(printme, "%g\n", out);
	return printme;
}

char* do_db_paired_t_test(char **argbuffer) {
char * 		printme	= malloc(sizeof(char)*2000);
double		out;
	out	= apop_db_paired_t_test(argbuffer[0],argbuffer[1],argbuffer[2]);
	sprintf(printme, "%g\n", out);
	return printme;
}

char* do_det_and_inv(char **argbuffer) {
char * 		printme	= malloc(sizeof(char)*2000);
size_t		matrix_no, matrix_out;
double		out;
	matrix_no	= find_matrix(argbuffer[1]);
	if (strcmp(argbuffer[0],"0")){
		matrix_out	= find_matrix(argbuffer[0]);
		//free(matrix_list[matrix_out]);//didn't need it allocated...
		out		= apop_det_and_inv(matrix_list[matrix_no], &(matrix_list[matrix_out]),1,1);
	}
	else	out		= apop_det_and_inv(matrix_list[matrix_no], NULL,1,0);
	sprintf(printme, "%g\n", out);
	return printme;
}

char* do_matrix_multiply(char **argbuffer) {
int	matrix_out, matrix_1, matrix_2;
     	matrix_out	= find_matrix(argbuffer[0]);
     	matrix_1	= find_matrix(argbuffer[1]);
     	matrix_2	= find_matrix(argbuffer[2]);
	matrix_list[matrix_out]	= gsl_matrix_calloc(matrix_list[matrix_1]->size1, matrix_list[matrix_2]->size2);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, matrix_list[matrix_1], matrix_list[matrix_2], 0, matrix_list[matrix_out]); 
	return NULL;
}

char* do_query(char **argbuffer) {
	apop_query(argbuffer[0]);
	return NULL;
}

char* do_query_to_matrix(char **argbuffer) {
size_t		matrix_no;
	matrix_no		= find_matrix(argbuffer[0]);
	matrix_list[matrix_no]	= apop_query_to_matrix(argbuffer[1]);
	return NULL;
}

char* do_db_merge(char **argbuffer) {
	apop_db_merge(argbuffer[0]);
	return 0;
}

char* print_help(char **argbuffer); //This is below the command_list, since it refers to it.

/** print a matrix via command line. */
char* do_print_matrix(char **argbuffer) {
int		i, j, len= 0;
gsl_matrix *	m;
char		*out	= malloc(sizeof(char)*3),
		bite[10000];
     	m	= matrix_list[find_matrix(argbuffer[0])];
	out[0]	='\0';
	//printf("m is a %i by %i matrix.\n", m->size1, m->size2);
	for(i=0; i< m->size1; i++){
		for(j=0; j< m->size2; j++){
			sprintf(bite, "%g", gsl_matrix_get(m,i,j));
			len	+= strlen(bite)+1;
			out	 = realloc(out, sizeof(char)*(len+2));
			strcat(out, bite);
			strcat(out, "\t");
		}
		len	++;
		out	 = realloc(out, sizeof(char)*(len+1));
		strcat(out, "\n");
	}
	return out ;
}


/*
typedef struct operation{
	char 	name[100]; //if your command is longer than this, it's too much.
	int	argc;
	char 	internal_return;
	char* 	(*run_me)(char **argbuffer);
	char 	help_text[];
} operation;
*/

operation	command_list[]=
{{"start", 		0, 'x',0,  null_fn, "db_name.db: start the server using the given db (or leave blank for in-memory db)" },
{"db_merge", 		1, 'x',0,  do_db_merge, "db_name.db: merge named db into currently open db" },
{"db_merge_table", 	2, 'x',0,  do_db_merge_table, "db_name.db table: just merge in one table."},
{"db_t_test",	 	4, 'x',1,  do_db_t_test, "tab1 col1 tab2 col2: returns confidence with which we reject mean(col1)==mean(col2)"},
{"db_paired_t_test", 	3, 'x',1,  do_db_paired_t_test, "tab1 col1 tab2 col2: returns confidence with which we reject mean(col1)==mean(col2)"},
{"det_and_inv",	 	2, 'm',1,  do_det_and_inv, "inverse_name tab_to_invert: take the input table, save its inverse. Immediately returns the determinant. if inverse_name==0, don't save the inverse."},
{"query", 		1, 'x',0, do_query, "query_text: run query, return nothing."},
{"query_to_matrix", 	2, 'm',0, do_query_to_matrix, "matrix_name query: query db, dump result to matrix_name" },
{"matrix_multiply", 	3, 'm',0, do_matrix_multiply, "matrix_out left_matrix right_matrix: out = left dot right"},
{"print_matrix", 	1, 'x',1, do_print_matrix, "matrix_name: print a matrix you've already created"},
{"stop", 		0, 'x',0, null_fn, 	": close the server"},
{"help", 		0, 'x',1, print_help, 	": prints this."},
{"xxxx", 		0, 'x',0, null_fn, 	""}
};

char* print_help(char **argbuffer) {
int		i	= 0;
	while (strcmp(command_list[i].name, "xxxx")){
		printf("%s %s\n", command_list[i].name, command_list[i].help_text);
		i	++;
	}
	return NULL;
}






///// OK, that's it for the functions themselves. The rest is the server code.



//This section checks the structures declared above to find the list
//element we want.
int find_cmd_list_elmt(char *cmd){
int		i = 0;
	while (strcmp(command_list[i].name, "xxxx")){
		if (!strcmp(command_list[i].name, cmd))
			return i;
		i++;
	}
	printf("Sorry, I DK that command.\n");
	return -1;
}

int find_matrix(char *name){
int		i;
	for(i=0;i<matrix_ct; i++)
		if (!strcmp(matrix_names[i], name)){
			//printf("found matrix %s in position %i.\n", name, i);
			return i;
		}
	matrix_ct	++;
	matrix_list	= realloc(matrix_list, sizeof(gsl_matrix*)*(matrix_ct));
	matrix_names	= realloc(matrix_names, sizeof(char*)*(matrix_ct));
	matrix_names[matrix_ct-1] = malloc(sizeof(char*)*(strlen(name)+1));
	strcpy(matrix_names[matrix_ct-1],name);
	//printf("added matrix %s in position %i.\n", name, matrix_ct);
	return matrix_ct - 1;
}



/////////////
//Below is the actual server/client code.
//

void error(char *msg) {
    perror(msg);
    exit(1);
}

int open_server() {
     int sockfd, portno;
     struct sockaddr_in serv_addr;
     sockfd = socket(AF_INET, SOCK_STREAM, 0);
     if (sockfd < 0) 
        error("ERROR opening socket");
     bzero((char *) &serv_addr, sizeof(serv_addr));
     portno = PORT_NO;
     serv_addr.sin_family = AF_INET;
     serv_addr.sin_addr.s_addr = INADDR_ANY;
     serv_addr.sin_port = htons(portno);
     if (bind(sockfd, (struct sockaddr *) &serv_addr,
              sizeof(serv_addr)) < 0) 
              error("ERROR on binding");
     return sockfd;
}


/** This function does all the real work. It listens for a command,
and when it gets one, it executes it here. */
int start_listening(int sockfd){
struct sockaddr_in 	cli_addr;
int 			n, clilen, newsockfd,  cmd_no, i, stop_pt;
char 			cmdbuffer[256], **argbuffer,  
			*printme;
     clilen = sizeof(cli_addr);
     newsockfd = accept(sockfd, 
                 (struct sockaddr *) &cli_addr, 
                 &clilen);
     if (newsockfd < 0) error("Error on accept");
	//get the cmd:
     bzero(cmdbuffer,256);
     n 		= read(newsockfd,cmdbuffer,255);
     if (n < 0) error("Error reading from socket");
     cmd_no	= find_cmd_list_elmt(cmdbuffer);
     if (cmd_no == -1)
	printf("couldn't find command %s.\n", cmdbuffer);
     else {
     	n 		= write(newsockfd,"confirmed.",strlen("confirmed."));
     	if (n < 0) error("Error confirming.");
     }
     
     stop_pt	= find_cmd_list_elmt("stop");
     if(stop_pt == cmd_no){
     	apop_close_db(0);
     	close(newsockfd);
     	return 0;
     	//momma socket closes when we return from this fn.
     }

     argbuffer	= malloc(sizeof(char*) * command_list[cmd_no].argc);
     for (i=0; i < command_list[cmd_no].argc; i++){
	     argbuffer[i]	= malloc(sizeof(char) * 1000);
     		bzero(argbuffer[i],256);
     		n = read(newsockfd,argbuffer[i],255);
     		if (n < 0) error("Error reading from socket");
     		//printf("Here is the arg for %s: %s\n",cmdbuffer, argbuffer[i]);
     		n 		= write(newsockfd,"confirmed.",strlen("confirmed."));
     		if (n < 0) error("Error confirming.");
	}

     ///////////////////////////////////////////
     //The actual command handling follows here

     printme	= command_list[cmd_no].run_me(argbuffer);
     if (printme){
	n 	= write(newsockfd, printme,strlen(printme));
    	if (n < 0) error("Error writing to socket");
     }
 
	//clean
	for (i=0; i < command_list[cmd_no].argc; i++)
	     free(argbuffer[i]);
	free(argbuffer);
	if (printme !=NULL) free(printme);
	return 1; 
}

int start_client(){
    int sockfd, portno, optval=1;
    struct sockaddr_in serv_addr;
    struct hostent *server;

    portno = PORT_NO;
    sockfd = socket(AF_INET, SOCK_STREAM, 0);

  /* setsockopt: Handy debugging trick that lets 
   * us rerun the server immediately after we kill it; 
   * otherwise we have to wait about 20 secs. 
   * Eliminates "ERROR on binding: Address already in use" error. 
   */
  setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, (const void *)&optval , sizeof(int));

    if (sockfd < 0) 
        error("ERROR opening socket");
    server = gethostbyname("localhost");
    if (server == NULL) {
        fprintf(stderr,"ERROR, no such host\n");
        exit(0);
    }
    bzero((char *) &serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    bcopy((char *)server->h_addr, 
         (char *)&serv_addr.sin_addr.s_addr,
         server->h_length);
    serv_addr.sin_port = htons(portno);
    if (connect(sockfd,(struct sockaddr *)&serv_addr,sizeof(serv_addr)) < 0) 
        error("ERROR connecting");
    return sockfd;
}


#define BUFLEN 10000
int send_cmd(char *cmd, char argc, char **args){
int		n, i, sockfd;
char		buffer[BUFLEN] = "  ";
	sockfd	= start_client();
	n 	= write(sockfd, cmd,strlen(cmd));
	if (n < 0) error("ERROR writing to socket");
	n	= read(sockfd, buffer, BUFLEN);
    	if (strcmp("confirmed.", buffer))
		printf("Darn. The server didn't confirm the command.\n");
	for (i=0;i<argc; i++){
    		n = write(sockfd, args[i],strlen(args[i]));
    		if (n < 0) error("ERROR writing to socket");
    		n	= read(sockfd, buffer, BUFLEN);
		if (strcmp("confirmed.", buffer))
			printf("The server didn't confirm argument %i.\n", i);
	}
    /*
    bzero(buffer,100000);
    n = read(sockfd,buffer,100000-1);
    if (n < 0) error("ERROR reading from socket");
    printf("%s\n",buffer);
    */
    return sockfd;
}


void dump_output(int client_sock){
char		buffer[BUFLEN];
	if (read(client_sock, buffer, BUFLEN) < 0)
		error("couldn't dump stuff.");
	printf(buffer);
}

int main (int argc, char ** argv){
int		sock, client_sock, cmd_no, arg_ct,
		keep_going	= 1;
	if (!strcmp(argv[1], "help") || argc==1){
		print_help(NULL);
		return 0;
	}

	//if it's a server, open the main socket, start listening.
	if (!strcmp(argv[1], "start")){
		sock	= open_server();
		if (sock < 0) printf("fuck.");
		if (argc==3)
			apop_open_db(argv[2]);
		else
			apop_open_db(NULL);
     		listen(sock,2);
		while(keep_going)
			keep_going	= start_listening(sock);
		close(sock);
		return 0;	//You're done. Don't look thru command lists.
	}

	//else, it's a client (or a typo). No loop here: one command per call.
	cmd_no	= find_cmd_list_elmt(argv[1]);
	if (cmd_no >= 0){
		arg_ct	= command_list[cmd_no].argc;
		if (arg_ct == 0)
			client_sock	= send_cmd(argv[1], arg_ct, NULL);
		else 	client_sock	= send_cmd(argv[1], arg_ct, (argv+2));
		if (command_list[cmd_no].external_return)
			dump_output(client_sock);
		close(client_sock);
	}
	return 0;
}
