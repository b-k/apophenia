#include <apophenia/headers.h>
#include <stddef.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <sys/types.h> 
#include <netinet/in.h>
#include <netdb.h> 

#define PORT_NO 11716

char cmd_list[][2][100]= { 	{"stop", "0"},
				{"query", "1"},
				{"query_to_matrix", "2"},
				{"matrix_print", "1"},
				{"xxxx", "0"}};

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
		if (!strcmp(matrix_names[i], name))
			return i;
	matrix_ct	++;
	matrix_list	= realloc(matrix_list, sizeof(gsl_matrix*)*(matrix_ct));
	matrix_names	= realloc(matrix_names, sizeof(char*)*(matrix_ct));
	strcpy(matrix_names[matrix_ct-1],name);
	return matrix_ct - 1;
}

void error(char *msg)
{
    perror(msg);
    exit(1);
}

int open_server() {
     int sockfd, portno, clilen;
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
int		i,j, len= 0;
char		*out	= NULL,
		bite[10000];
	for(i=0; i< m->size1; i++){
		for(j=0; j< m->size1; j++){
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
int 			n, clilen, newsockfd, argct, cmd_no, i, matrix_no;
char 			cmdbuffer[256], argbuffer[100][10000], return_msg[10000], *printme;
     listen(sockfd,2);
     clilen = sizeof(cli_addr);
     newsockfd = accept(sockfd, 
                 (struct sockaddr *) &cli_addr, 
                 &clilen);
     if (newsockfd < 0) error("ERROR on accept");
	//get the cmd:
     bzero(cmdbuffer,256);
     n 		= read(newsockfd,cmdbuffer,255);
     if (n < 0) error("ERROR reading from socket");
     cmd_no	= find_cmd_list_elmt(cmdbuffer);

     //get the args; new socket.
     listen(sockfd,2);
     clilen = sizeof(cli_addr);
     newsockfd = accept(sockfd, 
                 (struct sockaddr *) &cli_addr, 
                 &clilen);
     if (newsockfd < 0) error("ERROR on accept");
     sprintf(return_msg,"Here is the cmd: %s\n",cmdbuffer);
     n 		= write(newsockfd,return_msg,strlen(return_msg));
     if (n < 0) error("ERROR writing to socket");
     for (i=0; i < atoi(cmd_list[cmd_no][1]); i++){
     	bzero(argbuffer[i],256);
     	n = read(newsockfd,argbuffer[i],255);
     	if (n < 0) error("ERROR reading from socket");
     	printf("Here is the arg: %s\n",argbuffer[i]);
	}
     if (!strcmp(cmdbuffer, "stop")){
	     apop_close_db(0);
	     return 0;
	     }
     if (!strcmp(cmdbuffer, "print_matrix")){
	     	matrix_no	= find_matrix(argbuffer[0]);
	     	printme	= print_matrix(matrix_list[matrix_no]);
    		n = write(newsockfd, printme,strlen(printme));
    		if (n < 0) error("ERROR writing to socket");
	}
     if (!strcmp(cmdbuffer, "query"))
	     apop_query(argbuffer[0]);
     if (!strcmp(cmdbuffer, "query_to_matrix")){
	     matrix_no	= find_matrix(argbuffer[0]);
	     matrix_list[matrix_no]	= apop_query_to_matrix(argbuffer[1]);
	     }
     return 1; 
}

int start_client(){
    int sockfd, portno;
    struct sockaddr_in serv_addr;
    struct hostent *server;
    char buffer[256];

    portno = PORT_NO;
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
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

int send_cmd(char *cmd, char argc, char **args){
int		n, i, sockfd;
char		buffer[100000];
 
	sockfd	= start_client();
    n = write(sockfd, cmd,strlen(cmd));
    if (n < 0) error("ERROR writing to socket");

	sockfd	= start_client();
    for (i=0;i<argc; i++){
    		n = write(sockfd, args[i],strlen(args[i]));
    		if (n < 0) error("ERROR writing to socket");
	}
    /*
    bzero(buffer,100000);
    n = read(sockfd,buffer,100000-1);
    if (n < 0) error("ERROR reading from socket");
    printf("%s\n",buffer);
    */
    return 0;
}

int main (int argc, char ** argv){
char		*socket_name;
int		sock, insock,
		keep_going	= 1;
struct sockaddr	client_sock;
socklen_t	client_length;
	if (!strcmp(argv[1], "start")){
		sock	= open_server();
		if (sock < 0) printf("fuck.");
		if (argc==2)
			apop_open_db(argv[2]);
		else
			apop_open_db(NULL);
		while(keep_going)
			keep_going	= start_listening(sock);
	}
	if (!strcmp(argv[1], "stop"))
		send_cmd("stop", 0, NULL);
	if (!strcmp(argv[1], "query"))
		send_cmd("query", 1, (argv+2));
	return 0;
}
