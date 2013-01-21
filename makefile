#!/bin/bash

distdir=dists
pkg=apophenia
version=0.99
date=`date +%d_%b_%y`
workdir=$(distdir)/$(pkg)-$(version)

#See README.developer for notes.

#This script copies everything into a junk directory within dists/, which
#GNU's prep materials will produce a thousand temp files in, concluding 
#with a file named apophenia-N.NN.tar.gz. Once that it is ready
#for distribution, the tar.gz arcive is copied to the dists/ directory.
#If the appropriate tools are available, it also produces Debian pkgs.
#Finally, the script unzips the archive and compiles from the archive.

all:
	@echo This is the \'backstage\' version of the Apophenia repository.
	@echo It builds the .deb, .rpm, and autoconf-able .tgz versions.
	@echo Therefore, it may require package-making tools that you don\'t have.
	@echo If you just want to use the library, check the download page at
	@echo https://github.com/b-k/Apophenia/downloads
	@echo [but ignore the buttons at the top of the page].
	@echo If you want to hack the library itself, welcome---the makefile 
	@echo will show you what targets are available.

install: auto
	cd $(workdir) && make
	cd $(workdir) && sudo make install

auto:
	rm -rf $(distdir)
	mkdir $(distdir)
	mkdir $(workdir)
	cp -rf `ls -I $(distdir)` $(workdir)
	rm $(workdir)/makefile #the one you're reading now---it's not for the pkg itself.
	cd $(workdir) && cp install/* .
	sed "s/PKGNAME/apophenia-$(version)-$(date).tgz/" < install/rpm.spec > $(workdir)/rpm.spec
	cd $(workdir) && cp ChangeLog NEWS
	cd $(workdir) && sed -i "s/X\.XX/$(version)/" README apop_db.c configure.ac upload
	cd $(workdir) && sed -i -f prep_variadics.sed *.c *h model/*.c
	#generate the internal header based on the public apop.h
	sed -e 's/[<>]/"/g' -e 's/apophenia\///' -e '1i#include "variadic.h"' -e '1i#include "internal.h"\n' <apop.h > $(workdir)/apop_internal.h
	#add a directory with links to all headers, for testing.
	mkdir $(workdir)/tests/apophenia
	cp $(workdir)/*.h $(workdir)/tests/apophenia
	#let the GNU prep everything and produce a distribution pkg.
	cd $(workdir) && dot -Tpng < docs/structs.dot | pngtopnm | pnmscale .5 | pnmtopng > docs/structs.png
	cd $(workdir) && echo "Ben Klemens (fluffmail@f-m.fm)" > AUTHORS
	cd $(workdir) && libtoolize
	cd $(workdir) && aclocal
	cd $(workdir) && autoconf 
	cd $(workdir) && autoheader
	cd $(workdir) && automake -a
	cd $(workdir) && ./configure #--enable-python
	-cd $(workdir) && make distcheck
	mv $(workdir)/apophenia-$(version).tar.gz $(distdir)/apophenia-$(version)-$(date).tgz 

rpm: auto
	if [ ! -x ~/rpmbuild] ; then  \
	    mkdir ~/rpmbuild;       \
	    cd ~/rpmbuild && mkdir -p SOURCES BUILD RPMS/i686 RPMS/i386 SPECS SRPMS; \
	fi
	sed "s/PKGNAME/apophenia-$(version)-$(date).tgz/" < install/rpm.spec > ~/rpmbuild/SPECS/rpm.spec
	cp dists/apophenia-$(version)-$(date).tgz ~/rpmbuild/SOURCES
	rpmbuild -ba ~/rpmbuild/SPECS/rpm.spec
	cp ~/rpmbuild/RPMS/*/apophenia*rpm dists/.

deb: rpm
	cd $(distdir) && mkdir alien_tmp && cd alien_tmp && fakeroot alien -d ~/rpmbuild/RPMS/*/apophenia*rpm \
	&& mv apop*deb .. && cd .. && rm -rf alien_tmp

upload: deb
	sed -e "s/VERSION/$(version)/" -e "s/DATE/$(date)/" < install/upload > dists/upload
	cd dists && chmod +x upload && ./upload

deb-install: deb
	sudo dpkg -i $(distdir)/apophenia*deb


deb-raw:
    #These are notes on making debian packages using the debian control file
    #I never really got it to work 100%. As above, I now just use alien to generate 
    #the deb from the RPM (which itself only needs a simple spec file to be built).
#Prerequisites: You will need the debian build tools; try:
# apt-get install dpkg-dev debhelper devscripts fakeroot linda dh-make 
#mv apophenia-*.deb dists
#mv apophenia-*.gz dists
#mv apophenia_* dists
#	./configure --enable-python
#    yes | dh_make -l -n -c GPL  
#    cp ../../install/deb.control debian/control
#    cp ../../install/deb.copyright debian/copyright
#    debuild -us -uc 
#    rm ../$pkg-$version.tar.gz
# or: sudo apt-get install build-essential autoconf automake autotools-dev dh-make debhelper devscripts fakeroot xutils lintian pbuilder
