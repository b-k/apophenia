PROGNAME = test

objects =$(PROGNAME).o
#CFLAGS = -g -Wall -O3
CFLAGS = -g -Wall -pthread #-O3 
LINKFLAGS = -L/usr/local/lib -lapophenia -lgsl -lgslcblas -lsqlite3
CC = gcc
#CC = g++

#By default, this runs latex twice on the document named above.
#Alternatively, you can run any of the targets listed below.

c: $(objects)
	$(CC) $(CFLAGS) $(objects) $(LINKFLAGS) -o $(PROGNAME)

$(objects): %.o: %.c 
		$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(PROGNAME).o  $(PROGNAME)

