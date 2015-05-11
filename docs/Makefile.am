
#The Doxygen documentation, which you'll have to call yourself (via make doc).

default:

## adhoc
html man: doc

apophenia_CSOURCES = \
	$(top_srcdir)/model/*.c \
	$(top_srcdir)/transform/*.c
	echo #$(wildcard $(top_srcdir)/model/*.c) \
	echo #$(wildcard $(top_srcdir)/transform/*.c)

apophenia_IMAGES = \
	flake.gif \
	search.gif \
	right.png \
	down.png

apophenia_IMAGES_EXTRA = \
	triangle.png \
	model.png

apophenia_JS = \
	tree.js

model_doc.h: $(apophenia_CSOURCES)
	cat $^ | awk -f make_model_doc.awk > $@

doc: documentation.h model_doc.h $(apophenia_IMAGES) $(apophenia_JS)
	cp ../apop.h abkup
	sed -i -e 's/__attribute.*;/;/' -e "s/extern //" ../apop.h
	doxygen doxygen.conf
	cp abkup ../apop.h
	cp $(apophenia_IMAGES) $(apophenia_JS) html/
	sed -i -f edit_outline.sed html/index.html html/outline.html
	sed -i -f edit_group.sed html/group__models.html
	sed -i -f edit_width.sed html/*.html
	. adjust

doc-clean:
	-rm -rf html latex man

CLEANFILES = \
	missing_model_parts

MAINTAINERCLEANFILES = \
	model_doc.h \
	doxygen.log

maintainer-clean-local: doc-clean

EXTRA_DIST = \
	make_model_doc.awk \
	$(apophenia_IMAGES) \
	$(apophenia_IMAGES_EXTRA) \
	$(apophenia_JS) \
	edit_outline.sed edit_globals.sed edit_group.sed edit_width.sed \
	apop_data_fig.html head.html foot.html \
	typical.css \
	documentation.h

