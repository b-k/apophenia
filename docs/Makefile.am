
#The Doxygen documentation, which you'll have to call yourself (via make doc).

GVZDOT = /usr/bin/dot

default:

## adhoc
html man: doc

apophenia_CSOURCES = \
	$(top_srcdir)/model/*.c \
	$(top_srcdir)/transform/*.c
	echo #$(wildcard $(top_srcdir)/model/*.c) \
	echo #$(wildcard $(top_srcdir)/transform/*.c)

apophenia_DOTS = \
	structs.dot

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

apophenia_IMAGES_GENERATED = \
	structs.png

model_doc.h: $(apophenia_CSOURCES)
	cat $^ | awk -f $(top_srcdir)/docs/make_model_doc.awk > $@

doc: documentation.h model_doc.h $(apophenia_IMAGES) $(apophenia_JS) $(apophenia_IMAGES_GENERATED)
	$(MKDIR_P) include
	sed -e 's/__attribute.*;/;/;s/extern //' $(top_builddir)/apop.h > include/apop.h
	doxygen doxygen.conf
	for f in $(apophenia_IMAGES) $(apophenia_JS) ; do \
		test -f $(top_srcdir)/docs/$$f && cp $(top_srcdir)/docs/$$f html/ ;\
	done
	cp $(apophenia_IMAGES_GENERATED) html/
	sed -i -f $(top_srcdir)/docs/edit_outline.sed html/index.html html/outline.html
	sed -i -f $(top_srcdir)/docs/edit_group.sed html/group__models.html
	sed -i -f $(top_srcdir)/docs/edit_width.sed html/*.html
	$(abs_top_srcdir)/docs/adjust

doc-clean:
	-rm -rf include html latex man

install-html-local: doc
	cp -prd html $(DESTDIR)$(docdir)

maintainer-clean-local: doc-clean

CLEANFILES = \
	missing_model_parts

MAINTAINERCLEANFILES = \
	model_doc.h \
	doxygen.log

EXTRA_DIST = \
	adjust \
	make_model_doc.awk \
	$(apophenia_DOTS) \
	$(apophenia_IMAGES) \
	$(apophenia_IMAGES_EXTRA) \
	$(apophenia_JS) \
	edit_outline.sed edit_globals.sed edit_group.sed edit_width.sed \
	apop_data_fig.html head.html foot.html \
	typical.css \
	documentation.h

%.png : %.dot
	$(GVZDOT) -Tpng < $< > $@
