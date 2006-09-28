DOCNAME = gsl_stats
#DOCNAME = c_crash

all: $(DOCNAME).dvi

index: $(DOCNAME).idx

publish: $(DOCNAME).pdf
#	scp $(DOCNAME).pdf $(SFSSH):ap/htdocs/gsl_stats.pdf
	scp $(DOCNAME).pdf $(FLUFF):www/ben/pdfs/gsl_stats.pdf

onerun:
	latex $(DOCNAME)
	@#\citeindexentry will print \indexentry{{Kmenta}}, which means
	@#that Kmenta gets alphabetized under open-bracket. So
	@#turn {{...}} into {...}.
	@perl -pi.bak -e 's/indexentry{{([^}]*)}}/indexentry{$$1}/'  .texputs/$(DOCNAME).idx;

pdf:	$(DOCNAME).pdf

bib:	$(DOCNAME).dvi $(DOCNAME).bbl

book:	$(DOCNAME).pdf
#	psbook gsl_stats.ps | mpage -2 -o -M-63 | ps2pdf - > gsl_eff.pdf
	mpage -2 -o -M-63 gsl_stats.ps | ps2pdf - > gsl_eff.pdf
	scp gsl_eff.pdf $(THEBES):www/

figures:
	dot -Tps < testing.dot >testfig.eps
	dot -Tps < models.dot >modeldotfig.eps


viewps: 
	gv -seascape  -geometry 1129x824+-95+-49 -scale 2 $(DOCNAME).ps

$(DOCNAME).ps: $(DOCNAME).dvi
	dvips  -f < $(DOCNAME).dvi > $(DOCNAME).ps

$(DOCNAME).pdf: $(DOCNAME).ps
	ps2pdf $(DOCNAME).ps  $(DOCNAME).pdf

$(DOCNAME).bbl:	$(DOCNAME).dvi
	latex $(DOCNAME); bibtex $(DOCNAME); latex $(DOCNAME); latex $(DOCNAME)
	perl -pi.bak -e 's/In , /In /'  .texputs/$(DOCNAME).bbl;
	perl -pi.bak -e 's/{zarembka}, 1973./{zarembka}./'  .texputs/$(DOCNAME).bbl;
	

$(DOCNAME).idx:
	makeindex -c .texputs/$(DOCNAME).idx;

$(DOCNAME).dvi: *.tex
	make onerun
	make index
	@echo "---------------------------------------------------------------------------------"
	@echo "---------------------------------------------------------------------------------"
	@echo "---------------------------------------------------------------------------------"
	make onerun




#I should delete these:
view:
	xdvi -expert -geometry 1041x797+-17+-17 -s 6 -paper usr $(DOCNAME)

viewpdf: 
	gv -seascape  -geometry 1129x824+-95+-49 -scale 2 $(DOCNAME).pdf

pdfslides:	
	#latex $(DOCNAME)
	#dvips -t landscape -t letter -f < $(DOCNAME).dvi > $(DOCNAME).ps\
	dvips -t landscape -f < $(DOCNAME).dvi > $(DOCNAME).ps
	ps2pdf $(DOCNAME).ps  $(DOCNAME).pdf
