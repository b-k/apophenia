#this is a sample makefile. You may need to modify it.
#Once correct for your system, you can probably use it for all Apop/GSL projects.


PROGNAME = test

objects =$(PROGNAME).o
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

