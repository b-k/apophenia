Vendor:       Ben Klemens
Name:         apophenia
Version:      0.23
Release:      1
License:      GPLv2 w/Affero-type addendum
Provides:     apophenia, apophenia-devel
Requires:     gsl, gsl-devel, sqlite3, sqlite3-devel
Packager:     fluffmail@f-m.fm
Summary:      A library of functions and models for scientific computing.
Source:       PKGNAME
BuildRoot:    %{buildroot}
URL:          http://apophenia.info
%description
 Apophenia is a library of functions and models for scientific computing.
It is intended to make work easier---and even fun---when handling data sets, fitting
models, and designing new models in C. Facilities for managing data include an easy
link to SQLite3 or mySQL databases. Also includes a Python interface (in Beta).

%prep
%setup

%build
%configure 
make

%install
%makeinstall
#DESTDIR={%buildroot} %makeinstall

%files
%defattr(-,root,root)
/usr/lib/libapophenia.so.0.0.0
/usr/lib/libapophenia.so.0
/usr/lib/libapophenia.so
/usr/lib/libapophenia.la
/usr/lib/libapophenia.a
/usr/lib/pkgconfig/apophenia.pc
/usr/bin/apop_text_to_db
/usr/bin/apop_db_to_crosstab
/usr/bin/apop_merge_dbs
/usr/bin/apop_plot_query
/usr/bin/apop_lookup
/usr/include/apop.h
/usr/include/apophenia
/usr/include/apophenia/conversions.h
/usr/include/apophenia/db.h
/usr/include/apophenia/likelihoods.h
/usr/include/apophenia/arms.h
/usr/include/apophenia/mapply.h
/usr/include/apophenia/stats.h
/usr/include/apophenia/types.h
/usr/include/apophenia/output.h
/usr/include/apophenia/model.h
/usr/include/apophenia/asst.h
/usr/include/apophenia/settings.h
/usr/include/apophenia/variadic.h
/usr/include/apophenia/linear_algebra.h
/usr/include/apophenia/deprecated.h
