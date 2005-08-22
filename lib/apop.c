//apop.c		  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
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

//here are the commands implemented so far.
//Each line gives: name, argument count, internal return type, external return
//argument count includes any internal names.
//Internal return types: x==NULL, m=matrix.
//No types on external returns; just a zero or a one.
//Order doesn't matter except that xxxx has to be the last entry.
//
//Once you've added a line here, you'll also need to add the cmd
//handling in the start_listening function.
char cmd_list[][5][1000]= { 	
				{"start", 		"0", "x", "0", "db_name.db: start the server using the given db (or leave blank for in-memory db)"},
				{"db_merge", 		"1", "x", "0", "db_name.db: merge named db into currently open db"},
				{"db_merge_table", 	"2", "x", "0", "db_name.db table: just merge in one table."},
				{"db_t_test",	 	"4", "x", "1", "tab1 col1 tab2 col2: returns confidence with which we reject mean(col1)==mean(col2)"},
				{"db_paired_t_test", 	"3", "x", "1", "tab1 col1 tab2 col2: returns confidence with which we reject mean(col1)==mean(col2)"},
				{"det_and_inv",	 	"2", "m", "1", "inverse_name tab_to_invert: take the input table, save its inverse. Immediately returns the determinant. if inverse_name==0, don't save the inverse."},
				{"query", 		"1", "x", "0", "query_text: run query, return nothing."},
				{"query_to_matrix", 	"2", "m", "0", "matrix_name query: query db, dump result to matrix_name" },
				{"print_matrix", 	"1", "x", "1", "matrix_name: print a matrix you've already created"},
				{"stop", 		"0", "x", "0", ": close the server"},
				{"help", 		"0", "x", "1", ": prints this."},
				{"xxxx", 		"0", "x", "0"} };

gsl_matrix **	matrix_list;
char **		matrix_names;
int		matrix_ct	= 0;

int find_cmd_list_elmt(char *cmd){
int		i = 0;
	while (strcmp(cmd_list[i][0], "xxxx")){
		if (!strcmp(cmd_list[i][0], cmd))
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
			printf("found matrix %s in position %i.\n", name, i);
			return i;
		}
	matrix_ct	++;
	matrix_list	= realloc(matrix_list, sizeof(gsl_matrix*)*(matrix_ct));
	matrix_names	= realloc(matrix_names, sizeof(char*)*(matrix_ct));
	matrix_names[matrix_ct-1] = malloc(sizeof(char*)*(strlen(name)+1));
	strcpy(matrix_names[matrix_ct-1],name);
	printf("added matrix %s in position %i.\n", name, matrix_ct);
	return matrix_ct - 1;
}

void error(char *msg)
{
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

char *print_matrix(gsl_matrix *m){
int		i, j, len= 0;
char		*out	= malloc(sizeof(char)*3),
		bite[10000];
	out[0]	='\0';
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


int start_listening(int sockfd){
struct sockaddr_in 	cli_addr;
int 			n, clilen, newsockfd,  cmd_no, i, matrix_no, matrix_out;
char 			cmdbuffer[256], argbuffer[100][10000],  
			*printme	= malloc(sizeof(char)*1000);
double			out;
//gsl_matrix		**a_matrix;
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

     for (i=0; i < atoi(cmd_list[cmd_no][1]); i++){
     		bzero(argbuffer[i],256);
     		n = read(newsockfd,argbuffer[i],255);
     		if (n < 0) error("Error reading from socket");
     		//printf("Here is the arg for %s: %s\n",cmdbuffer, argbuffer[i]);
     		n 		= write(newsockfd,"confirmed.",strlen("confirmed."));
     		if (n < 0) error("Error confirming.");
	}

     ///////////////////////////////////////////
     //The actual command handling follows here
     //If the command returns nothing interesting, then it's just
     //called.
     //If the command returns a matrix, name, et cetera, then it needs
     //to be added to the appropriate array.
     ///////////////////////////////////////////

     if (!strcmp(cmdbuffer, "db_merge"))
	     apop_db_merge(argbuffer[0]);
     else if (!strcmp(cmdbuffer, "db_merge_table"))
	     apop_db_merge_table(argbuffer[0], argbuffer[1]);
     else if (!strcmp(cmdbuffer, "db_t_test")){
		out	= apop_db_t_test(argbuffer[0],argbuffer[1],argbuffer[2],argbuffer[3]);
		sprintf(printme, "%g\n", out);
     }
     else if (!strcmp(cmdbuffer, "db_paired_t_test")){
		out	= apop_db_paired_t_test(argbuffer[0],argbuffer[1],argbuffer[2]);
		sprintf(printme, "%g\n", out);
     }
     else if (!strcmp(cmdbuffer, "det_and_inv")){
     		matrix_no	= find_matrix(argbuffer[1]);
		if (strcmp(argbuffer[0],"0")){
	     		matrix_out	= find_matrix(argbuffer[0]);
			//free(matrix_list[matrix_out]);//didn't need it allocated...
	     		out		= apop_det_and_inv(matrix_list[matrix_no], &(matrix_list[matrix_out]),1,1);
		}
		else	out		= apop_det_and_inv(matrix_list[matrix_no], NULL,1,0);
		sprintf(printme, "%g\n", out);
	     }
     else if (!strcmp(cmdbuffer, "print_matrix")){
	     	matrix_no	= find_matrix(argbuffer[0]);
	     	printme		= print_matrix(matrix_list[matrix_no]);
	}
     else if (!strcmp(cmdbuffer, "query"))
	     apop_query(argbuffer[0]);
     else if (!strcmp(cmdbuffer, "query_to_matrix")){
	     matrix_no			= find_matrix(argbuffer[0]);
	     matrix_list[matrix_no]	= apop_query_to_matrix(argbuffer[1]);
	     }
     else if (!strcmp(cmdbuffer, "stop")){
	     apop_close_db(0);
	     close(newsockfd);
	     return 0;
	     //momma socket closes when we return from this fn.
	     }


     //I may need to write the output:
     if (atoi(cmd_list[cmd_no][3])){
    		n 		= write(newsockfd, printme,strlen(printme));
    		if (n < 0) error("Error writing to socket");
     }
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
char		buffer[BUFLEN];
 
	sockfd	= start_client();
    n = write(sockfd, cmd,strlen(cmd));
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

void print_help(){
int		i	= 0;
	while (strcmp(cmd_list[i][0], "xxxx")){
		printf("%s %s\n", cmd_list[i][0], cmd_list[i][4]);
		i	++;
	}
}

int main (int argc, char ** argv){
int		sock, client_sock, cmd_no, arg_ct,
		keep_going	= 1;
	if (!strcmp(argv[1], "help")){
		print_help();
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
	}

	//else, it's a client (or a typo). No loop here: one command per
	//call.
	cmd_no	= find_cmd_list_elmt(argv[1]);
	if (cmd_no >= 0){
		arg_ct	=atoi(cmd_list[cmd_no][1]);
		if (arg_ct == 0)
			client_sock	= send_cmd(argv[1], arg_ct, NULL);
		else 	client_sock	= send_cmd(argv[1], arg_ct, (argv+2));
		if (atoi(cmd_list[cmd_no][3]))
			dump_output(client_sock);
		close(client_sock);
	}
		/*
	if (!strcmp(argv[1], "stop"))
		send_cmd(argv[1], 0, NULL);
   	if (!strcmp(argv[1], "db_merge"))
                send_cmd(argv[1], 1, (argv+1));
        if (!strcmp(argv[1], "db_merge_table"))
                send_cmd(argv[1], 2, (argv+2));
	if (!strcmp(argv[1], "query"))
		send_cmd(argv[1], 1, (argv+2));
	if (!strcmp(argv[1], "query_to_matrix"))
		send_cmd(argv[1], 2, (argv+2));
	if (!strcmp(argv[1], "print_matrix")){
		client_sock	= send_cmd(argv[1], 1, (argv+2));
		dump_output(client_sock);
		}
		*/
	return 0;
}
