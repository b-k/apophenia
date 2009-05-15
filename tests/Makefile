PROGNAME = test

objects =$(PROGNAME).o
#CFLAGS = -g -Wall -O3
CFLAGS = -g -Wall -pthread #-O3 
LINKFLAGS = -L/usr/local/lib -lapophenia -lgsl -lgslcblas -lsqlite3
CC = gcc
#CC = g++

c: $(objects)
	$(CC) $(CFLAGS) $(objects) $(LINKFLAGS) -o $(PROGNAME)

$(objects): %.o: %.c 
		$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(PROGNAME).o  $(PROGNAME)

