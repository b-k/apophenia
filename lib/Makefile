
LLIBDIR = /usr/local/lib
objects = apop_stats.lo apop_db.lo apop_linear_algebra.lo apop_likelihoods.lo apop_conversions.lo
#CFLAGS = -g -Wall -O3
CFLAGS = -g -Wall
LINKFLAGS = -L/usr/local/lib -lgsl -lgslcblas -lsqlite3
LIBNAME = libapophenia.la

all: $(objects)
	libtool gcc $(CFLAGS) $(LINKFLAGS) -rpath $(LLIBDIR) -o $(LIBNAME) $(objects)

$(objects): %.lo: %.c 
	#gcc $(CFLAGS) -c $< -o $@
	libtool gcc $(CFLAGS) -c $<

install:
	libtool install -c $(LIBNAME) $(LLIBDIR)/$(LIBNAME)
	libtool -n --finish /usr/local/lib
	install ./*h /usr/local/include/apophenia
