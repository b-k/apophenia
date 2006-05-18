OBJECTS = $(PROGNAME).o 			#If there are more files, add nextfile.o here.
CFLAGS = -g -Wall
LINKFLAGS = -L/usr/local/lib -lapophenia -lgsl -lgslcblas -lsqlite3
COMPILE   = gcc $(CFLAGS) -c $< -o $@

$(PROGNAME): $(OBJECTS)
	gcc $(CFLAGS) $(OBJECTS) $(LINKFLAGS) -o $(PROGNAME)

$(objects): %.o: %.c 
		gcc $(CFLAGS) -c $< -o $@

clean:
	rm $(PROGNAME).o $(PROGNAME)
