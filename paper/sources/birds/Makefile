#PROGNAME = listflock
PROGNAME = treeflock
#PROGNAME = arrayflock

objects =$(PROGNAME).o birds.o
#CFLAGS = -g -Wall -O3
CFLAGS = -pg -g -Wall `pkg-config --cflags glib-2.0` #-O3 
LINKFLAGS = -L/usr/local/lib -lapophenia -lgsl -lgslcblas -lsqlite3 `pkg-config --libs glib-2.0`
CC = gcc
#CC = g++

#By default, this runs latex twice on the document named above.
#Alternatively, you can run any of the targets listed below.

c: $(objects)
	$(CC) $(CFLAGS) $(objects) $(LINKFLAGS) -o $(PROGNAME)


$(DOCNAME).dvi: $(DOCNAME).tex
	latex $(DOCNAME);echo "-------------------------------";latex $(DOCNAME)


$(objects): %.o: %.c 
		$(CC) $(CFLAGS) -c $< -o $@

debug:
	kdbg $(PROGNAME)&
linect:
	grep -i "void\|,\|;" *.c|wc -l
	grep "[^ 	]" *.c |wc -l 


setenv:
	LD_LIBRARY_PATH=/usr/local/lib:$$LD_LIBRARY_PATH
	export LD_LIBRARY_PATH


edit:
	bash --login -c "vi $(DOCNAME)_core.tex"

latex:			#One run.
	latex $(DOCNAME)

pdf:	$(DOCNAME).pdf
	rsync -avPe ssh $(DOCNAME).pdf  ben@131.215.22.96:diss/$(DOCNAME).pdf

bib:	$(DOCNAME).dvi $(DOCNAME).bbl
#	latex $(DOCNAME); bibtex $(DOCNAME); latex $(DOCNAME); latex $(DOCNAME)

view:
	xdvi -expert -geometry 1041x797+-17+-17 -s 6 $(DOCNAME)

viewpdf: pdf
	gv -geometry 1129x824+-95+-49 -scale 2 $(DOCNAME)

clean:
	rm -f $(PROGNAME).o  $(PROGNAME)

#
#
#
#
#

$(DOCNAME).bbl:	$(DOCNAME).dvi
	bibtex $(DOCNAME); touch $(DOCNAME)_core.tex; make 
	#bibtex $(DOCNAME); latex $(DOCNAME); latex $(DOCNAME)

$(DOCNAME).ps: $(DOCNAME).dvi
	dvips  -f < $(DOCNAME).dvi > $(DOCNAME).ps

$(DOCNAME).pdf: $(DOCNAME).ps
	ps2pdf $(DOCNAME).ps  $(DOCNAME).pdf
