
objects = stats.o db.o linear_algebra.o likelihoods.o
#CFLAGS = -g -Wall -O3
CFLAGS = -g -Wall
LINKFLAGS = -L/usr/local/lib -lgsl -lgslcblas -lsqlite

$(objects): %.o: %.c 
		gcc $(CFLAGS) -c $< -o $@
