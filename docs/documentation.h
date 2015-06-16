/* Apophenia's narrative documentation
Copyright (c) 2005--2013 by Ben Klemens.  Licensed under the GPLv2; see COPYING.  */

/**  \mainpage Welcome

Apophenia is an open statistical library for working with data sets and statistical
models. It provides functions on the same level as those of the typical stats package
(such as OLS, Probit, or singular value decomposition) but gives the user more
flexibility to be creative in model-building.  The core functions are written in C,
but experience has shown them to be easy to bind to in Python/Julia/Perl/Ruby/&c.

It is written to scale well, to comfortably work with gigabyte data sets, million-step simulations, or
computationally-intensive agent-based models. If you have tried using other open source
tools for computationally demanding work and found that those tools weren't up to the
task, then Apophenia is the library for you.

<h5>The goods</h5> 

The library has been growing and improving since 2005, and has been downloaded over 10,000 times. To date, it has over two hundred functions to facilitate statistical computing, such as:

\li OLS and family, discrete choice models like Probit and Logit, kernel density estimators, and other common models
\li database querying and maintenance utilities
\li moments, percentiles, and other basic stats utilities
\li t-tests, F-tests, et cetera
\li Several optimization methods available for your own new models
\li It does <em>not</em> re-implement basic matrix operations or build yet another database
engine. Instead, it builds upon the excellent <a href="http://www.gnu.org/software/gsl/">GNU
Scientific</a> and <a href="http://www.sqlite.org/">SQLite</a> libraries. MySQL/mariaDB is also supported.

For the full list, click the <a href="group__all__public.html">index</a> link from the header.

<h5><a href="https://github.com/b-k/apophenia/archive/pkg.zip">Download Apophenia here</a>.</h5>

Most users will want to download the latest packaged version linked from the <a
href="https://github.com/b-k/apophenia/archive/pkg.zip">Download
Apophenia here</a> header.

Those who would like to work on a cutting-edge copy of the source code
can get the latest version by cutting and pasting the following onto
the command line. If you follow this route, be sure to read the development README in the
<tt>apophenia</tt> directory this command will create.

\code
git clone https://github.com/b-k/apophenia.git
\endcode

<!--git clone git://apophenia.git.sourceforge.net/gitroot/apophenia/apophenia
cvs -z3 -d:ext:<i>(your sourceforge login)</i>@cvs.sourceforge.net:/cvsroot/apophenia co -P apophenia
cvs -z3 -d:pserver:anonymous@cvs.sf.net:/cvsroot/apophenia checkout -P apophenia
svn co https://apophenia.svn.sourceforge.net/svnroot/apophenia/trunk/apophenia --> 

<h5>The documentation</h5>

To start off, have a look at this \ref gentle "Gentle Introduction" to the library.

<a href="outline.html">The outline</a> gives a more detailed narrative.

The <a href="globals.html">index</a> lists every function in the
library, with detailed reference
information. Notice that the header to every page has a link to the outline and the index.

To really go in depth, download or pick up a copy of <a
href="http://modelingwithdata.org">Modeling with Data</a>, which discusses general
methods for doing statistics in C with the GSL and SQLite, as well as Apophenia
itself. <a href="http://www.census.gov/srd/papers/pdf/rrs2014-06.pdf"><em>A Useful
Algebraic System of Statistical Models</em></a> (PDF) discusses some of the theoretical
structures underlying the library.

There is a <a href="https://github.com/b-k/apophenia/wiki">wiki</a> with some convenience
functions, tips, and so on.

<h5>Notable features</h5> 
Much of what Apophenia does can be done in any typical statistics package. The \ref
apop_data element is much like an R data frame, for example, and there is nothing special
about being able to invert a matrix or take the product of two matrices with a single
function call (\ref apop_matrix_inverse and \ref apop_dot, respectively). 
Even more advanced features like Loess smoothing (\ref apop_loess) and the Fisher Exact
Test (\ref apop_test_fisher_exact) are not especially Apophenia-specific. But here are
some things that are noteworthy.

\li It's a C library! You can build applications using Apophenia for the data-processing
back-end of your program, and not worry about the overhead associated with scripting
languages. For example, it is currently used in production for certain aspects of
processing for the U.S. Census Bureau's American Community Survey. And the numeric
routines in your favorite scripting language typically have a back-end in plain C;
perhaps Apophenia can facilitate writing your next one.

\li As well as the \ref apop_data structure, Apophenia is built around a model object,
the \ref apop_model. This allows for consistent treatment of distributions, regressions,
simulations, machine learning models, and who knows what other sorts of models you can
dream up. By transforming and combining existing models, it is easy to build complex
models from simple sub-models.

\li For example, the \ref apop_update function does Bayesian updating on any two
well-formed models. If they are on the table of conjugates, that is correctly
handled, and if they are not, an appropriate variant of MCMC produces an empirical
distribution. The output is yet another model, from which you can make random draws,
or which you can use as a prior for another round of Bayesian updating. Outside of
Bayesian updating, the \ref apop_model_metropolis function is good for approximating
other complex models.

\li The text file parser is flexible and effective. Such data files are typically
called `CSV files', meaning <em>comma-separated values</em>, but the delimiter can be
anything (or even some mix of things), and there is no requirement that text have
"special delimiters". Missing data can be specified by a simple blank or a marker
of your choosing (e.g., <tt>apop_opts.nan_string = "N/A";</tt>). Or there can be
no delimiters, as in the case of fixed-width files. If you are a heavy SQLite user,
Apophenia may be useful to you simply for its \ref apop_text_to_db function.

\li The maximum likelihood system combines a lot of different subsystems into one
form: it will do a few flavors of conjugate gradient search, Nelder-Mead Simplex,
Newton's Method, or Simulated Annealing. You pick the method by a setting attached to
your model. If you want to use a method that requires derivatives and you don't have
a closed-form derivative, the ML subsystem will estimate a numerical gradient for
you. If you would like to do EM-style maximization (all but the first parameter are
fixed, that parameter is optimized, then all but the second parameter are fixed, that
parameter is optimized, ..., looping through dimensions until the change in objective
across cycles is less than <tt>eps</tt>), add a settings group specifying the
tolerance at which the cycle should stop: <tt>Apop_settings_add_group(your_model,
apop_mle, .dim_cycle_tolerance=eps)</tt>.

\li The Iterative Proportional Fitting algorithm, \ref apop_rake, is best-in-breed,
designed to handle large, sparse matrices.



<h5>Contribute!</h5> 

\li Develop a new model object.
\li Contribute your favorite statistical routine.
\li Package Apophenia into an RPM, apt, portage, cygwin package.
\li Report bugs or suggest features.
\li Write bindings for your preferred language. For example, here are early versions of <a
href="http://modelingwithdata.org/arch/00000173.htm"> a Julia
wrapper</a> and <a href="https://github.com/b-k/Rapophenia/">an R
wrapper</a> which you could expand upon.

If you're interested,  <a href="mailto:fluffmail@f-m.fm">write to the maintainer</a> (Ben Klemens), or join the
<a href="https://github.com/b-k/apophenia">GitHub</a> project.
*/

/** \page eg Some examples
Here are a few pieces of sample code, many gathered from elsewhere in the documentation,
for testing your installation or to give you a sense of what code with Apophenia's
tools looks like. If you'd like more context or explanation, please click through to
the page from which the example was taken.

<em> Two data streams</em>

The sample program here is intended to show how one would integrate Apophenia into an existing program. For example, say that you are running a simulation of two different treatments, or say that two sensors are posting data at regular intervals. You need to gather the data in an organized form, and then ask questions of the resulting data set.  Below, a thousand draws are made from the two processes and put into a database. Then, the data is pulled out, some simple statistics are compiled, and the data is written to a text file for inspection outside of the program.  This program will compile cleanly with the sample \ref makefile.

\include draw_to_db.c

<em> Run a regression</em>

The documentation for the \ref apop_ols model provides a program to read in data and run a regression.

<em> A sequence of t-tests</em>

In \ref mapply "The section on map/apply", a new \f$t\f$-test on every row, with all
operations acting on entire rows rather than individual data points:

\include t_test_by_rows.c

In the documentation for \ref apop_query_to_text, a program to list all the tables in an SQLite database.
\include ls_tables.c

\par Marginal distribution

A demonstration of fixing parameters to create a marginal distribution, via \ref apop_model_fix_params
\include fix_params.c

\par Several uses of the \ref apop_dot function

\include dot_products.c
*/



/** \page setup Setting up
\section cast The supporting cast 
To use Apophenia, you will need to have a working C compiler, the GSL (v1.7 or higher) and SQLite installed. mySQL/mariaDB is optional.

\li Some readers are unfamiliar with modern package managers and common methods for setting up a C development environment; see 
<a href="http://modelingwithdata.org/appendix_o.html">Appendix O</a> of <em> Modeling with Data</em> for an introduction.

\li Other pages in this documentation have a few more notes for \ref windows "Windows" users, including \ref mingw users.

\li Install the basics using your package manager. E.g., try

\code
sudo apt-get install make gcc libgsl0-dev libsqlite3-dev
\endcode

or

\code
sudo yum install make gcc gsl-devel libsqlite3x-devel
\endcode

\li <a href="https://github.com/b-k/apophenia/archive/pkg.zip">Download Apophenia here</a>. 

\li Once you have the library downloaded, compile it using

\code
tar xvzf apop*tgz && cd apophenia-0.999
./configure && make && sudo make install && make check
\endcode

If you decide not to keep the library on your system, run <tt>sudo make uninstall</tt>
from the source directory to remove it.

\li If you need to install packages in your home directory because you don't have root
permissions, see the \ref notroot page.

\li A \ref makefile will help immensely when you want to compile your program.

\li You can verify that your setup works by trying some \ref eg "sample programs".

\subpage notroot

\subpage makefile

\subpage windows

*/

/** \page windows Windows

\ref mingw users, see that page.

If you have a choice, <a href="http://www.cygwin.com">Cygwin</a> is strongly recommended. The setup program is
very self-explanatory. As a warning, it will probably take up &gt;300MB on
your system. You should install at least the following programs:
\li autoconf/automake
\li binutils
\li gcc
\li gdb
\li gnuplot -- for plotting data
\li groff -- needed for the man program, below
\li gsl -- the engine that powers Apophenia
\li less -- to read text files
\li libtool -- needed for compiling programs
\li make
\li man -- for reading help files
\li more -- not as good as less but still good to have
\li sqlite3 -- a simple database engine, a requisite for Apophenia

If you are missing anything else, the program will probably tell you.
The following are not necessary but are good to have on hand as long as you are going to be using Unix and programming.

\li svn -- to partake in the versioning system
\li emacs -- steep learning curve, but people love it
\li ghostscript (for reading .ps/.pdf files)
\li openssh -- needed for cvs
\li perl, python, ruby -- these are other languages that you might also be interested in
\li tetex -- write up your documentation using the nicest-looking formatter around
\li X11 -- a windowing system

X-Window will give you a nicer environment in which to work.  After you start Cygwin, type <tt>startx</tt> to bring up a more usable, nice-looking terminal (and the ability to do a few thousand other things which are beyond the scope of this documentation).  Once you have Cygwin installed and a good terminal running, you can follow along with the remainder of the discussion without modification.

Sqlite3 is difficult to build from scratch, but you can get a packaged version by pointing Cygwin's install program to the Cygwin Ports site: http://cygwinports.dotsrc.org/ .

Second, some older versions of Cygwin have a search.h file which doesn't include the function lsearch().  If this is the case on your system, you will have to update your Cygwin installation.

Finally, windows compilers often spit out lines like:
\code
Info: resolving _gsl_rng_taus by linking to __imp__gsl_rng_taus (auto-import)
\endcode
These lines are indeed just information, and not errors. Feel free to ignore them.

[Thanks to Andrew Felton and Derrick Higgins for their Cygwin debugging efforts.]

\subpage mingw
*/

/** \page notroot  Not root? 
If you aren't root, then the common procedure for installing a library is to create a subdirectory in your home directory in which to install packages. The key is the <tt>--prefix</tt> addition to the <tt>./configure</tt> command.
\code
export MY_LIBS = myroot   #choose a directory name to be created in your home directory.
mkdir $HOME/$MY_LIBS

# From Apophenia's package directory:
./configure --prefix $HOME/$MY_LIBS
make
make install   #Now you don't have to be root.

# Adjust your paths so the compiler and the OS can find the library.
# These are environment variables, and they are usually set in the 
# shell's startup files. I assume you are using bash here.

echo "export PATH=$HOME/$MY_LIBS/include:\$PATH" >> ~/.bashrc
echo "export CPATH=$HOME/$MY_LIBS/include:\$CPATH" >> ~/.bashrc
echo "export LIBRARY_PATH=$HOME/$MY_LIBS:\$LIBRARY_PATH" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=$HOME/$MY_LIBS:\$LD_LIBRARY_PATH" >> ~/.bashrc
\endcode

Once you have created this local root directory, you can use it to install as many new libraries as desired, and your paths will already be set up to find them.
*/


/** \page makefile Makefile
 
Instead of giving lengthy compiler commands at the command prompt, you can use a Makefile to do most of the work. How to:
\li Copy and paste the following into a file named \c makefile.
\li Change the first line to the name of your program (e.g., if you have written <tt>sample.c</tt>, then the first line will read <tt>PROGNAME=sample</tt>). 
\li If your program has multiple <tt>.c</tt> files, add a corresponding <tt>.o</tt> to the currently blank <tt>objects</tt> variable, e.g. <tt>objects=sample2.o sample3.o</tt>
\li One you have a Makefile in the directory, simply type <tt>make</tt> at the command prompt to generate the executable.

 \code
PROGNAME = your_program_name_here
objects =
CFLAGS = -g -Wall -O3
LDLIBS = -lapophenia -lgsl -lgslcblas -lsqlite3

$(PROGNAME): $(objects)
\endcode

\li If your system has \c pkg-config, then you can use it for a slightly more robust and readable makefile. Replace the above C and link flags with:
\code
CFLAGS = -g -Wall `pkg-config --cflags apophenia` -O3
LDLIBS = `pkg-config --libs apophenia`
\endcode
The \c pkg-config program will then fill in the appropriate directories and libraries. Pkg-config knows Apophenia depends on the GSL and database libraries, so you need only list the most-dependent library.

\li The -O3 flag is optional, asking the compiler to run its highest level of optimization (for speed).

\li GCC users may need the <tt>--std=gnu99</tt> or <tt>--std=gnu11</tt> flag to use post-1989 C standards.

\li Order matters in the linking list: the files a package depends on should be listed after the package. E.g., since sample.c depends on Apophenia, <tt>gcc sample.c -lapophenia</tt> will work, while <tt>gcc -lapophenia sample.c</tt> is likely to give you errors. Similarly, list <tt>-lapophenia</tt> before <tt>-lgsl</tt>, which comes before <tt>-lgslcblas</tt>.

*/

/** \page designated Designated initializers

Functions so marked in this documentation use standard C designated initializers and compound literals to allow you to omit, call by name, or change the order of inputs. The following examples are all equivalent.

The standard format:
\code
apop_text_to_db("infile.txt", "intable", 0, 1, NULL);
\endcode

Omitted arguments are left at their default vaules:
\code
apop_text_to_db("infile.txt", "intable");
\endcode

You can use the variable's name, if you forget its ordering:
\code
apop_text_to_db("infile.txt", "intable", .has_col_name=1, .has_row_name=0);
\endcode

If an un-named element follows a named element, then that value is given to the next variable in the standard ordering:
\code
apop_text_to_db("infile.txt", "intable", .has_col_name=1, NULL);
\endcode

\li There may be cases where you can not use this form (it relies on a macro, which
may not be available). You can always call the underlying function directly, by adding
\c _base to the name and giving all arguments:

\code
apop_text_to_db_base("infile.txt", "intable", 0, 1, NULL);
\endcode

\li If one of the optional elements is an RNG and you do not provide one, I use one
from \ref apop_rng_get_thread.
*/


/** \page preliminaries Getting started

If you are entirely new to Apophenia, \ref gentle "have a look at the Gentle Introduction here".

As well as the information in this outline, there is a separate page covering the details of 
 \ref setup "setting up a computing environment" and another page with \ref eg "some sample code" for your perusal.

For another concrete example of folding Apophenia into a project, have a look at this \ref sample_program "sample program".

\subpage gentle

\subpage setup

\subpage eg

\subpage refstatusext
*/

/** \page refstatusext References, status, and extensions

\section mwd The book version

Apophenia co-evolved with <em>Modeling with Data: Tools and Techniques for Statistical Computing</em>. You can read about the book, or download a free PDF copy of the full text, at <a href="http://modelingwithdata.org">modelingwithdata.org</a>.

If you are at this site, there is probably something there for you, including a tutorial on C and general computing form, SQL for data-handing, several chapters of statistics from various perspectives, and more details on working Apophenia. 

As with many computer programs, the preferred manner of citing Apophenia is to cite its related book.
Here is a BibTeX-formatted entry giving the relevant information:

\code 
@book{klemens:modeling,
    title = "Modeling with Data: Tools and Techniques for Statistical Computing",
    author="Ben Klemens",
    year=2008,
    publisher="Princeton University Press"
}
\endcode

The rationale for the \ref apop_model struct, based on an algebraic system of models, is detailed in a <a href="http://www.census.gov/srd/papers/pdf/rrs2014-06.pdf">U.S. Census Bureau research report</a>:

\code 
@techreport{klemens:algebra,
    title = "A Useful Algebraic System of Statistical Models",
    author="Ben Klemens",
    month=jul,
    year=2014,
    institution="U.S.\ Census Bureau",
    number="06"
}
\endcode
\section status What is the status of the code?

[This section last updated 11 May 2015.]

Apophenia was first posted to SourceForge in February 2005, which means that we've had
several years to develop and test the code in real-world applications. 

The test suite, including the sample code and solution set for <em>Modeling with Data</em>,
is about 5,500 lines over 135 files. gprof reports that it covers over 90% of the
7,700 lines in Apophenia's code base. A broad rule of thumb for any code base is
that the well-worn parts, in this case functions like \ref apop_data_get and \ref
apop_normal's <tt>log_likelihood</tt>, are likely to be entirely reliable, while the
out-of-the-way functions (maybe the score for the Beta distribution) will always be worth a bit
of caution. Close to all of the code has been used in production, so all of it was at
least initially tested against real-world data.

It is currently at version 0.999, which is intended to indicate that it is substantially
complete. Of course, a library for scientific computing, or even for that small subset
that is statistics, will never cover all needs and all methods. But as it stands
Apophenia's framework, based on the \ref apop_data and \ref apop_model, is basically
internally consistent, has enough tools that you can get common work done quickly,
and is reasonably fleshed out with a good number of models out of the box.

The \ref apop_data structure is set, and there are enough functions there that you could
use it as a subpackage by itself (especially in tandem with the database functions)
for nontrivial dealings with data.

The \ref apop_model structure is much more ambitious---Apophenia is intended
to be a novel system for developing models.
The promise underlying the structure is that you can provide just one item, such as
an RNG or a likelihood function, and the structure will do all of the work to fill in
computationally-intensive methods for everything else; see \ref settingswriting for
the details. Some directions aren't quite there yet (such as RNG -> most other things),
the PMF model needs an internal index for faster lookups, and so on.  Readers are invited
to contribute better methods (such as an alternate means of estimating mixture models),
or filling in more of existing models (write a dlog likelihood function for a model
that does not currently have one), or submit new standard models not yet included.


\section ext How do I write extensions?

The system is written to not require a registration or initialization step to add a new
model or other such parts.  Just write your code and <tt>include</tt> it like any other
C code.  A new \ref apop_model has to conform to some rules if it is to play well with
\ref apop_estimate, \ref apop_draw, and so forth. See the notes at \ref modeldetails.
Once your new model or function is working, please post the code or a link to the code
on the <a href="https://github.com/b-k/apophenia/wiki">Apophenia wiki</a>.

\subpage c

\section links Further references

For your convenience, here are links to some other libraries you are probably using.

\li <a href="http://www.gnu.org/software/libc/manual/html_node/index.html">The standard C library</a>
\li <a href="http://www.gnu.org/software/gsl/manual/html_node/index.html">The
GSL documentation</a>, and <a href="http://www.gnu.org/software/gsl/manual/html_node/Function-Index.html">its index</a>
\li <a href="http://sqlite.org/lang.html">SQL understood by SQLite</a>


*/

/** \page outline An outline of the library

The narrative in this section goes into greater detail on how to use the 
how to use the components of Apophenia. You are encouraged to read \ref gentle first.

This overview begins with the \ref apop_data set, which is the central data structure
used all over the system. Section \ref dbs covers the use of the database interface,
because there are a lot of things that a database will always do better than a matrix
structure like the \ref apop_data struct.

Section \ref modelsec covers statistical models, in the form of the \ref apop_model structure.
This part of the system is built upon the \ref apop_data set to hold parameters, statistics, data sets, and so on.

Section \ref Histosec covers probability mass functions, which are statistical models
built around a single data set, where the chance of drawing a given observation is
proportional to how often that observation appears in the source data. There are many
situations where one would want to treat a data set as a probability distribution,
such as using \ref apop_kl_divergence to find the information loss from an observed
data set to a theoretical model fit to that data.

Section \ref testpage covers traditional hypothesis testing, beginning with common
statistics that take an \ref apop_data set or two as input, and continuing on to
generalized hypothesis testing for any \ref apop_model.

Because estimation in the \ref apop_model relies heavily on maximum likelihood
estimation, Apophenia's optimizer subsystem is extensive.  Section \ref maxipage offers
some additional notes on optimization and how it can be used in non-statistical contexts.

\subpage dataoverview

\subpage dbs

\subpage modelsec

\subpage Histosec

\subpage testpage

\subpage maxipage

\subpage moreasst
*/

/** \page c C, SQL and coding utilities
 
\par  Learning C

<a href="http://modelingwithdata.org">Modeling with Data</a> has a full tutorial for C, oriented at users of standard stats packages. More nuts-and-bolts tutorials are <a href="http://www.google.com/search?hl=es&amp;c2coff=1&amp;q=c+tutorial">in abundance</a>.  Some people find pointers to be especially difficult; fortunately, there's a <a href="http://www.youtube.com/watch?v=6pmWojisM_E">claymation cartoon</a> which clarifies everything.

\par Header aggregation

There is only one header. Put
\code
#include <apop.h>
\endcode
at the top of your file, and you're done. Everything declared in that file starts with \c apop_ (or \c Apop_).

\par Linking

You will need to link to the Apophenia library, which involves adding the <tt>-lapophenia</tt> flag to your compiler. Apophenia depends on SQLite3 and the GNU Scientific Library (which depends on a BLAS), so you will probably need something like:

\code
gcc sample.c -lapophenia -lsqlite3 -lgsl -lgslcblas -o run_me -g -Wall -O3
\endcode

Your best bet is to encapsulate this mess in a \ref makefile "Makefile". Even if you are using an IDE and its command-line management tools, see the Makefile page for notes on useful flags.

\par Standards compliance

To the best of our abilities, Apophenia complies to the C standard (ISO/IEC
9899:2011). As well as relying on the GSL and SQLite, it uses some POSIX function calls,
such as \c strcasecmp and \c popen.

\par Easier calling syntax

Many functions allow optional named arguments to functions. For
example:
\code
apop_vector_distance(v1, v2); //assumes Euclidean distance
apop_vector_distance(v1, .metric='M'); //assumes v2=0, uses Manhattan metric.
\endcode

See the \ref designated page for details of this syntax.

\subpage designated

\section debugging  Errors, logging, debugging and stopping

<h5>The \c error element</h5> 

The \ref apop_data set and the \ref apop_model both include an element named \c error. It is normally \c 0, indicating no (known) error. 

For example, \ref apop_data_copy detects allocation errors and some circular links
(when <tt>Data->more == Data</tt>) and fails in those cases. You could thus use the
function with a form like

\code
apop_data *d = apop_text_to_data("indata");
apop_data *cp = apop_data_copy(d);
if (cp->error) {printf("Couldn't copy the input data; failing.\n"); return 1;}
\endcode

There is sometimes (but not always) benefit to handling specific error codes, which are listed in the documentation of those functions that set the \c error element. E.g.,

\code
apop_data *d = apop_text_to_data("indata");
apop_data *cp = apop_data_copy(d);
if (cp->error == 'a') {printf("Couldn't allocate space for the copy; failing.\n"); return 1;}
if (cp->error == 'c') {printf("Circular link in the data set; failing.\n"); return 2;}
\endcode

The end of <a href="http://modelingwithdata.org/appendix_o.html">Appendix O</a>
of <em>Modeling with Data</em> offers some GDB macros which can make dealing with
Apophenia from the GDB command line much more pleasant. As discussed below, it
also helps to set <tt>apop_opts.stop_on_warning='v'</tt> or <tt>'w'</tt> when running
under the debugger.


\section verbsec Verbosity level and logging

The global variable <tt>apop_opts.verbose</tt> determines how many notifications and warnings get printed by Apophenia's warning mechanism:

-1: turn off logging, print nothing (ill-advised) <br>
0: notify only of failures and clear danger <br>
1: warn of technically correct but odd situations that might indicate, e.g., numeric instability <br>
2: debugging-type information; print queries  <br>
3: give me everything, such as the state of the data at each iteration of a loop.

These levels are of course subjective, but should give you some idea of where to place the
verbosity level. The default is 1.

The messages are printed to the \c FILE* handle at <tt>apop_opts.log_file</tt>. If
this is blank (which happens at startup), then this is set to \c stderr. This is the
typical behavior for a console program. Use

\code
apop_opts.log_file = fopen("mylog", "w");
\endcode

to write to the \c mylog file instead of \c stderr.

As well as the error and warning messages, some functions can also print diagnostics,
using the \ref Apop_notify macro.  For example, \ref apop_query and friends will print the
query sent to the database engine iff <tt>apop_opts.verbose >=2</tt> (which is useful
when building complex queries). The diagnostics attempt to follow
the same verbosity scale as the warning messages.

\section Stopping

Warnings and errors never halt processing. It is up to the calling function to decide
whether to stop.

When running the program under a debugger, this is an annoyance: we want to stop as
soon as a problem turns up.

The global variable <tt>apop_opts.stop_on_warning</tt> changes when the system halts:

\c 'n': never halt. If you were using Apophenia to support a user-friendly GUI, for example, you would use this mode.<br>
The default: if the variable is <tt>'\0'</tt> (the default), halt on severe errors, continue on all warnings.<br>
\c 'v': If the verbosity level of the warning is such that the warning would print to screen, then halt;
if the warning message would be filtered out by your verbosity level, continue.<br>
\c 'w': Halt on all errors or warnings, including those below your verbosity threshold.

See the documentation for individual functions for details on how each reports errors to the caller and the level at which warnings are posted.

\section Legi Legible output

The output routines handle four sinks for your output. There is a global variable that
you can use for small projects where all data will go to the same place.

\code 
apop_opts.output_type = 'f'; //named file
apop_opts.output_type = 'p'; //a pipe or already-opened file
apop_opts.output_type = 'd'; //the database
\endcode

You can also set the output type, the name of the output file or table, and other options
via arguments to individual calls to output functions. See \ref apop_prep_output for the list of options.

C makes minimal distinction between pipes and files, so you can set a
pipe or file as output and send all output there until further notice:

\code
apop_opts.output_type = 'p';
apop_opts.output_pipe = popen("gnuplot", "w");
apop_plot_lattice(...); //see https://github.com/b-k/Apophenia/wiki/gnuplot_snippets
fclose(apop_opts.output_pipe);
apop_opts.output_pipe = fopen("newfile", "w");
apop_data_print(set1);
fprintf(apop_opts.output_pipe, "\nNow set 2:\n");
apop_data_print(set2);
\endcode

Continuing the example, you can always override the global data with
a specific request:
\code
apop_vector_print(v, "vectorfile"); //put vectors in a separate file
apop_matrix_print(m, "matrix_table", .output_type = 'd'); //write to the db
apop_matrix_print(m, .output_pipe = stdout);  //now show the same matrix on screen
\endcode

I will first look to the input file name, then the input pipe, then the
global \c output_pipe, in that order, to determine to where I should
write.  Some combinations (like output type = \c 'd' and only a pipe) don't
make sense, and I'll try to warn you about those. 

What if you have too much output and would like to use a pager, like \c less or \c more?
In C and POSIX terminology, you're asking to pipe your output to a paging program. Here is
the form:
\code
FILE *lesspipe = popen("less", "w");
assert(lesspipe);
apop_data_print(your_data_set, .output_pipe=lesspipe);
pclose(lesspipe);
\endcode
\c popen will search your usual program path for \c less, so you don't have to give a full path.

\li\ref apop_data_print  
\li\ref apop_matrix_print
\li\ref apop_vector_print
\li\ref apop_data_show : alias for \ref apop_data_print limited to \c stdout.

\section sqlsec About SQL, the syntax for querying databases

For a reference, your best bet is the <a href="http://www.sqlite.org/lang.html">Structured Query Language reference</a> for SQLite.  For a tutorial; there is an abundance of <a href="http://www.google.com/search?q=sql+tutorial">tutorials online</a>.  Here is a nice blog <a href="http://fluff.info/blog/arch/00000118.htm">entry</a> about complementaries between SQL and matrix manipulation packages.

Apophenia currently supports two database engines: SQLite and mySQL/mariaDB. SQLite is the default, because it is simpler and generally more easygoing than mySQL, and supports in-memory databases.

The global <tt>apop_opts.db_engine</tt> is initially \c NUL, indicating no preference
for a database engine. You can explicitly set it:

\code
apop_opts.db_engine='s' //use SQLite
apop_opts.db_engine='m' //use mySQL/mariaDB
\endcode

If \c apop_opts.db_engine is still \c NUL on your first database operation, then I will check
for an environment variable <tt>APOP_DB_ENGINE</tt>, and set 
<tt>apop_opts.db_engine='m'</tt> if it is found and matches (case insensitive) \c mariadb or \c mysql.

\code
export APOP_DB_ENGINE=mariadb
apop_text_to_db indata mtab db_for_maria

unset APOP_DB_ENGINE
apop_text_to_db indata stab db_for_sqlite.db
\endcode


Write \ref apop_data sets to the database using \ref apop_data_print, with <tt>.output_type='d'</tt>.

\li Column names are inserted if there are any. If there are, all dots are converted
    to underscores.  Otherwise, the columns will be named \c c1, \c c2, \c c3, &c.
\li If \ref apop_opts_type "apop_opts.db_name_column" is not blank (the default is
    <tt>"row_name"</tt>), then a so-named column is created, and the row names are placed there.
\li If there are weights, they will be the last column of the table, and the column will be named \c weights.
\li If the table exists; append to. If the table does not exist, create. So perhaps
    call \ref apop_table_exists <tt>("tabname", 'd')</tt> to ensure that the table is
    removed ahead of time.
\li Use \ref apop_data_print <tt>(data, "tabname", .output_type='d', .output_append='w')</tt>
    to overwrite a new table or with <tt>.output_append='a'</tt> to append.
\li If your data set has zero data (i.e., is just a list of column names or is entirely
    blank), \ref apop_data_print returns without creating anything in the database.
\li Especially if you are using a pre-2007 version of SQLite, there may be a speed
    gain to wrapping the call to this function in a begin/commit pair:

\code
apop_query("begin;");
apop_data_print(dataset, .output_name="dbtab", .output_type='d');
apop_query("commit;");
\endcode


Finally, Apophenia provides a few nonstandard SQL functions to facilitate math via database; see \ref db_moments.


\section threads Threading

Apophenia uses OpenMP for threading. You generally do not need to know how OpenMP works
to use Apophenia, and many points of work will thread without your doing anything.

\li All functions strive to be thread-safe. Part of how this is achieved is that static
variables are marked as thread-local or atomic, as per the C standard. There still
exist compilers that can't implement thread-local or atomic variables, in which case
your safest bet is to set OMP's thread count to one as below (or get a new compiler).

\li Some functions modify their inputs. It is up to you to use those functions in
a thread-safe manner. The \ref apop_matrix_realloc handles states and global variables
correctly in a threaded environment, but if you have two threads resizing the same \c
gsl_matrix at the same time, you're going to have problems.

\li There are few compilers that don't support OpenMP. Clang on MacOS may be the only
current mainstream example as of this writing, and they are hard at work on implementing
it. In the mean time, when compiling on such a system all work will be single-threaded.

\li Set the maximum number of threads to \c N with the environment variable

\code
export OMP_NUM_THREADS N
\endcode

or the C function

\code
#include <omp.h>
omp_set_num_threads(N);
\endcode

Use one of these methods with <tt>N=1</tt> if you want a single-threaded program.
You can return later to using all available threads via <tt>omp_set_num_threads(omp_get_num_procs())</tt>.

\li \ref apop_map and friends distribute their \c for loop over the input \ref apop_data
set across multiple threads. Therefore, be careful to send thread-unsafe functions to
it only after calling \c omp_set_num_threads(1).

\li There are a few functions, like \ref apop_model_draws, that rely on \ref apop_map, and
therefore also thread by default.

\li The function \ref apop_rng_get_thread retrieves a statically-stored RNG specific
to a given thread. Therefore, if you use that function in the place of a \c gsl_rng,
you can parallelize functions that make random draws.

\li \ref apop_rng_get_thread allocates its store of threads using <tt>apop_opts.rng_seed</tt>,
then incrementing that seed by one. You thus probably have threads with seeds 479901,
479902, 479903, .... [If you have a better way to do it, please feel free to modify the
code to implement your improvement and submit a pull request on Github.]

See <a href="http://modelingwithdata.org/arch/00000175.htm">this tutorial on C
threading</a> if you would like to know more, or are unsure about whether your functions
are thread-safe or not.
*/


/** \page dataoverview Data sets

The \ref apop_data structure represents a data set.  It joins together a \c gsl_vector,
a \c gsl_matrix, an \ref apop_name, and a table of strings. It tries to be lightweight,
so you can use it everywhere you would use a \c gsl_matrix or a \c gsl_vector.

Here is a diagram showing a sample data set with all of the elements in place. Together,
they represent a data set where each row is an observation, which includes both numeric
and text values, and where each row/column may be named.

\htmlinclude apop_data_fig.html
\latexinclude apop_data_fig.tex

In a regression, the vector would be the dependent variable, and the other columns
(after factor-izing the text) the independent variables. Or think of the \ref apop_data
set as a partitioned matrix, where the vector is column -1, and the first column of
the matrix is column zero. Here is some sample code to print the vector and matrix,
starting at column -1 (but you can use \ref apop_data_print to do this).

\code
for (int j = 0; j< data->matrix->size1; j++){
    printf("%s\t", apop_name_get(data->names, j, 'r'));
    for (int i = -1; i< data->matrix->size2; i++)
        printf("%g\t", apop_data_get(data, j, i));
    printf("\n");
}
\endcode

Most functions assume that each row represents one observation, so the data vector,
data matrix, and text have the same row count: \c data->vector->size==data->matrix->size1
and \c data->vector->size==*data->textsize. This means that the \ref apop_name structure
doesn't have separate \c vector_names, \c row_names, or \c text_row_names elements:
the \c rownames are assumed to apply for all.

See below for notes on easily managing the \c text element and the row/column names.

\section pps Pages

The \ref apop_data set includes a \c more pointer, which will typically be \c NULL,
but may point to another \ref apop_data set. This is intended for a main data set
and a second or third page with auxiliary information, such as estimated parameters
on the front page and their covariance matrix on page two, or predicted data on the
front page and a set of prediction intervals on page two.

The \c more pointer is not intended as a linked list for millions of data points. In
such situations, you can often improve efficiency by restructuring your data to use
a single table (perhaps via \ref apop_data_pack and \ref apop_data_unpack).

Most functions, such as \ref apop_data_copy and \ref apop_data_free, will handle all
the pages of information. For example, an optimization search over multi-page parameter
sets would search the space given by all pages.

But pages may also be appended as output or auxiliary information, such as
covariances, and an MLE would not search over these elements. Any page with a name in
XML-ish brackets, such as <tt>\<Covariance\></tt>, is considered information about the
data, not data itself, and therefore ignored by search routines, missing data routines,
et cetera. This is achieved by a rule in \ref apop_data_pack and \ref apop_data_unpack.

Here is a toy example that establishes a baseline data set, adds a page,
modifies it, and then later retrieves it.
\code
apop_data *d = apop_data_alloc(10, 10, 10); //the base data set, a 10-item vector + 10x10 matrix
apop_data *a_new_page = apop_data_add_page(d, apop_data_alloc(2,2), "new 2 x 2 page");
gsl_vector_set_all(a_new_page->matrix, 3);

//later:
apop_data *retrieved = apop_data_get_page(d, "new", 'r'); //'r'=search via regex, not literal match.
apop_data_show(retrieved); //print a 2x2 grid of 3s.
\endcode

\section datafns Functions for using apop_data sets

There are a great many functions to collate, copy, merge, sort, prune, and otherwise
manipulate the \ref apop_data structure and its components.

\li\ref apop_data_add_named_elmt
\li\ref apop_data_copy
\li\ref apop_data_fill
\li\ref apop_data_memcpy
\li\ref apop_data_pack
\li\ref apop_data_rm_columns
\li\ref apop_data_sort
\li\ref apop_data_split
\li\ref apop_data_stack
\li\ref apop_data_transpose : transpose matrices (square or not) and text grids
\li\ref apop_data_unpack
\li\ref apop_matrix_copy
\li\ref apop_matrix_realloc
\li\ref apop_matrix_stack
\li\ref apop_text_set
\li\ref apop_text_paste
\li\ref apop_text_to_data
\li\ref apop_vector_copy
\li\ref apop_vector_fill
\li\ref apop_vector_stack
\li\ref apop_vector_realloc

Apophenia builds upon the GSL, but it would be inappropriate to redundantly replicate
the <a href="http://www.gnu.org/software/gsl/manual/html_node/index.html">GSL's documentation</a> here.
Meanwhile, here are prototypes for a few common functions. The GSL's
naming scheme is very consistent, so a simple reminder of the function name may be
sufficient to indicate how they are used.

\li <tt>gsl_matrix_swap_rows (gsl_matrix * m, size_t i, size_t j)</tt>
\li <tt>gsl_matrix_swap_columns (gsl_matrix * m, size_t i, size_t j)</tt>
\li <tt>gsl_matrix_swap_rowcol (gsl_matrix * m, size_t i, size_t j)</tt>
\li <tt>gsl_matrix_transpose_memcpy (gsl_matrix * dest, const gsl_matrix * src)</tt>
\li <tt>gsl_matrix_transpose (gsl_matrix * m) : square matrices only</tt>
\li <tt>gsl_matrix_set_all (gsl_matrix * m, double x)</tt>
\li <tt>gsl_matrix_set_zero (gsl_matrix * m)</tt>
\li <tt>gsl_matrix_set_identity (gsl_matrix * m)</tt>
\li <tt>gsl_matrix_memcpy (gsl_matrix * dest, const gsl_matrix * src)</tt>
\li <tt>void gsl_vector_set_all (gsl_vector * v, double x)</tt>
\li <tt>void gsl_vector_set_zero (gsl_vector * v)</tt>
\li <tt>int gsl_vector_set_basis (gsl_vector * v, size_t i)</tt>: set all elements to zero, but set item \f$i\f$ to one.
\li <tt>gsl_vector_reverse (gsl_vector * v)</tt>: reverse the order of your vector's elements
\li <tt>gsl_vector_ptr</tt> and <tt>gsl_matrix_ptr</tt>. To increment an element in a vector use, e.g., <tt>*gsl_vector_ptr(v, 7) += 3;</tt> or <tt>(*gsl_vector_ptr(v, 7))++</tt>.
\li <tt>gsl_vector_memcpy (gsl_vector * dest, const gsl_vector * src)</tt>

\subsection readin Reading from text files

The \ref apop_text_to_data() function takes in the name of a text file with a grid of data in (comma|tab|pipe|whatever)-delimited format and reads it to a matrix. If there are names in the text file, they are copied in to the data set. See \ref text_format for the full range and details of what can be read in.

If you have any columns of text, then you will need to read in via the database: use
\ref apop_text_to_db() to convert your text file to a database table, 
do any database-appropriate cleaning of the input data, then use \ref
apop_query_to_data() or \ref apop_query_to_mixed_data() to pull the data to an \ref apop_data set.

\subpage text_format

\section datalloc Alloc/free

You may not need to use these functions often, given that \ref apop_query_to_data, \ref apop_text_to_data, and many transformation functions will auto-allocate \ref apop_data sets for you.

The \ref apop_data_alloc function allocates a vector, a matrix, or both. After this call, the structure will have blank names, \c NULL \c text element, and \c NULL \c weights.  See \ref names for discussion of filling the names. Use \ref apop_text_alloc to allocate the \c text grid. The \c weights are a simple \c gsl_vector, so allocate a 100-unit weights vector via <tt>allocated_data_set->weights = gsl_vector_alloc(100)</tt>.

Examples of use can be found throughout the documentation; for example, see \ref gentle.

\li\ref apop_data_alloc
\li\ref apop_data_calloc
\li\ref apop_data_free
\li\ref apop_text_alloc : allocate or resize the text part of an \ref apop_data set.
\li\ref apop_text_free

    See also:

\li <tt>gsl_matrix * gsl_matrix_alloc (size_t n1, size_t n2)</tt>
\li <tt>gsl_matrix * gsl_matrix_calloc (size_t n1, size_t n2)</tt>
\li <tt>void gsl_matrix_free (gsl_matrix * m)</tt>
\li <tt>gsl_vector * gsl_vector_alloc (size_t n)</tt>
\li <tt>gsl_vector * gsl_vector_calloc (size_t n)</tt>
\li <tt>void gsl_vector_free (gsl_vector * v)</tt>


\section gslviews 	Using views

There are several macros for the common task of viewing a single row or column of a \ref
apop_data set.

\code
apop_data *d = apop_query_to_data("select obs1, obs2, obs3 from a_table");

//Get a column using its name. Note that the generated view, ov, is the
//last item named in the call to the macro.
Apop_col_t(d, "obs1", ov);
double obs1_sum = apop_vector_sum(ov);

//Get row zero of the data set's matrix as a vector; get its sum
double row_zero_sum = apop_vector_sum(Apop_rv(d, 0));

//Get a row or rows as a standalone one-row apop_data set
apop_data_print(Apop_r(d, 0));

//ten rows starting at row 3:
apop_data *d10 = Apop_rs(d, 3, 10);
apop_data_show(d10);

//Column zero's sum
gsl_vector *cv = Apop_cv(d, 0);
double col_zero_sum = apop_vector_sum(cv);
//or one one line:
double col_zero_sum = apop_vector_sum(Apop_cv(d, 0));

//Pull a 10x5 submatrix, whose origin element is the (2,3)rd
//element of the parent data set's matrix

double sub_sum = apop_matrix_sum(Apop_subm(d, 2,3, 10,5));
\endcode

Because these macros can be used as arguments to a function, these macros have abbreviated names to save line space.

\li\ref Apop_c
\li\ref Apop_r
\li\ref Apop_cv
\li\ref Apop_rv
\li\ref Apop_cs
\li\ref Apop_rs
\li\ref Apop_subm

A second set of macros have a slightly different syntax, taking the name of the object to be declared as the last argument. These can not be used as expressions such as function arguments.

\li\ref Apop_col_t
\li\ref Apop_row_t
\li\ref Apop_col_tv
\li\ref Apop_row_tv
\li\ref Apop_matrix_row
\li\ref Apop_matrix_col

The view is an automatic variable, not a pointer, and therefore disappears at the end
of the scope in which it is declared. If you want to retain the data after the function
exits, copy it to another vector:

\code
Apop_matrix_row(d->matrix, 2, rowtwo);
return apop_vector_copy(rowtwo);
\endcode

Curly braces always delimit scope, not just at the end of a function. 
These macros work by generating a number of local variables, which you may be able to see in your
debugger. When program evaluation exits a given block, all variables in that block are
erased. Here is some sample code that won't work:

\code
apop_data *outdata;
if (get_odd){
    outdata = Apop_r(data, 1);
} else {
    outdata = Apop_r(data, 0);
}
apop_data_show(outdata); //breaks: outdata points to out-of-scope variables.
\endcode

For this if/then statement, there are two sets of local variables
generated: one for the \c if block, and one for the \c else block. By the last line,
neither exists. You can get around the problem here by making sure to not put the macro
declaring new variables in a block. E.g.:

\code
apop_data *outdata = Apop_r(data, get_odd ? 1 : 0);
apop_data_show(outdata);
\endcode


This is a general rule about how variables declared in blocks will behave, but because the
macros obscure the variable declarations, it is worth watching out for here.


\section data_set_get Set/get

First, some examples:

\code
apop_data *d = apop_data_alloc(10, 10, 10);
apop_name_add(d->names, "Zeroth row", 'r');
apop_name_add(d->names, "Zeroth col", 'c');

//set cell at row=8 col=0 to value=27
apop_data_set(d, 8, 0, .val=27);
assert(apop_data_get(d, 8, .colname="Zeroth") == 27);
double *x = apop_data_ptr(d, .col=7, .rowname="Zeroth");
*x = 270;
assert(apop_data_get(d, 0, 7) == 270);

// This is invalid---the value doesn't follow the colname. Use .val=5.
// apop_data_set(d, .row = 3, .colname="Column 8", 5);  

// OK, to set (3, 8) to 5:
apop_data_set(d, 3, 8, 5);


//apop_data set holding a scalar:
apop_data *s = apop_data_alloc(1);
apop_data_set(s, .val=12);
assert(apop_data_get(s) == 12);

//apop_data set holding a vector:
apop_data *v = apop_data_alloc(12);
for (int i=0; i< 12; i++) apop_data_set(s, i, .val=i*10);
assert(apop_data_get(s,3) == 30);

//This is a common form from pulling from a list of named scalars, 
//produced via apop_data_add_named_elmt
double AIC = apop_data_get(your_model->info, .rowname="AIC");
\endcode

Some additional notes:

\li The versions that take a column/row name use  \ref apop_name_find
    for the search; see notes there on the name matching rules.
\li For those that take a column number, column -1 is the vector element. 
\li For those that take a column name, I will search the vector last---if I don't find the name among the matrix columns, but the name matches the vector name, I return column -1.
\li If you give me both a \c .row and a \c .rowname, I go with the name; similarly for \c .col and
\c .colname.
\li You can give me the name of a page, e.g.
\code
double AIC = apop_data_get(data, .rowname="AIC", .col=-1, .page="<Info>");
\endcode

\li Numeric values default to zero, which is how the examples above that treated the \ref apop_data set as a vector or scalar could do so relatively gracefully.

\li Thus, so <tt>apop_data_get(dataset, 1)</tt> gets item (1, 0) from the matrix
element of \c dataset. But as a do-what-I-mean exception, if there is no matrix element
but there is a vector, then this form will get vector element 1. Relying on this DWIM
exception is useful iff you can guarantee that a data set will have only a vector or
a matrix but not both. Otherwise, be explicit and use <tt>apop_data_get(dataset, 1, -1)</tt>.

The \ref apop_data_ptr function follows the lead of \c gsl_vector_ptr and \c
gsl_matrix_ptr, and like those functions, returns a pointer to the appropriate \c double.

\li\ref apop_data_get
\li\ref apop_data_set
\li\ref apop_data_ptr : returns a pointer to the element.
\li\ref apop_data_get_page : retrieve a named page from a data set. If you only need a few items, you can specify a page name to \c apop_data_(get|set|ptr).

    See also:

\li <tt>double gsl_matrix_get (const gsl_matrix * m, size_t i, size_t j)</tt>
\li <tt>double gsl_vector_get (const gsl_vector * v, size_t i)</tt>
\li <tt>void gsl_matrix_set (gsl_matrix * m, size_t i, size_t j, double x)</tt>
\li <tt>void gsl_vector_set (gsl_vector * v, size_t i, double x)</tt>
\li <tt>double * gsl_matrix_ptr (gsl_matrix * m, size_t i, size_t j)</tt>
\li <tt>double * gsl_vector_ptr (gsl_vector * v, size_t i)</tt>
\li <tt>const double * gsl_matrix_const_ptr (const gsl_matrix * m, size_t i, size_t j)</tt>
\li <tt>const double * gsl_vector_const_ptr (const gsl_vector * v, size_t i)</tt>
\li <tt>gsl_matrix_get_row (gsl_vector * v, const gsl_matrix * m, size_t i)</tt>
\li <tt>gsl_matrix_get_col (gsl_vector * v, const gsl_matrix * m, size_t j)</tt>
\li <tt>gsl_matrix_set_row (gsl_matrix * m, size_t i, const gsl_vector * v)</tt>
\li <tt>gsl_matrix_set_col (gsl_matrix * m, size_t j, const gsl_vector * v)</tt>


\section mapply   Map/apply

\anchor outline_mapply
These functions allow you to send each element of a vector or matrix to a function, either producing a new matrix (map) or transforming the original (apply).  The \c ..._sum functions return the sum of the mapped output.

There are two types, which were developed at different times. The \ref apop_map and
\ref apop_map_sum functions use variadic function inputs to cover a lot of different
types of process depending on the inputs. Other functions with types in their names,
like \ref apop_matrix_map and \ref apop_vector_apply, may be easier to use in some
cases. They use the same routines internally, so use whichever type is convenient.

You can do many things quickly with these functions.

Get the sum of squares of a vector's elements:

\code
  //given apop_data *dataset and gsl_vector *v:
double sum_of_squares = apop_map_sum(dataset, gsl_pow_2);
double sum_of_sqvares = apop_vector_map_sum(v, gsl_pow_2);
\endcode

Here, we create an index vector [\f$0, 1, 2, ...\f$].

\code
double index(double in, int index){return index;}
apop_data *d = apop_map(apop_data_alloc(100), .fn_di=index, .inplace='y');
\endcode

Given your log likelihood function, which acts on a \ref apop_data set with only one
row, and a data set where each row of the matrix is an observation, find the total
log likelihood via:

\code
static double your_log_likelihood_fn(apop_data * in)
     {[your math goes here]}

double total_ll = apop_map_sum(dataset, .fn_r=your_log_likelihood_fn);
\endcode

How many missing elements are there in your data matrix? 

\code
static double nan_check(const double in){ return isnan(in);}

int missing_ct = apop_map_sum(in, nan_check, .part='m');
\endcode

Get the mean of the not-NaN elements of a data set:

\code
static double no_nan_val(const double in){ return isnan(in)? 0 : in;}
static double not_nan_check(const double in){ return !isnan(in);}

static double apop_mean_no_nans(apop_data *in){
    return apop_map_sum(in, no_nan_val)/apop_map_sum(in, not_nan_check);
}
\endcode

The following program randomly generates a data set where each row is a list of numbers with a different mean. It then finds the \f$t\f$ statistic for each row, and the confidence with which we reject the claim that the statistic is less than or equal to zero.

Notice how the older \ref apop_vector_apply uses file-global variables to pass information into the functions, while the \ref apop_map uses a pointer to send parameters to the functions.

\include t_test_by_rows.c

One more toy example, demonstrating the use of \ref apop_map and \ref apop_map_sum :

\include apop_map_row.c


\li If \c apop_opts.thread_count is greater than one, then the matrix will be broken
into chunks and each sent to a different thread. Notice that the GSL is generally
threadsafe, and SQLite is threadsafe conditional on several commonsense caveats that
you'll find in the SQLite documentation. See \ref apop_rng_get_thread() to use the GSL's RNGs in a threaded environment.

\li The \c ...sum functions are convenience functions that call \c ...map and then add up the contents. Thus, you will need to have adequate memory for the allocation of the temp matrix/vector.

\li\ref apop_map
\li\ref apop_map_sum
\li\ref apop_matrix_apply
\li\ref apop_matrix_map
\li\ref apop_matrix_map_all_sum
\li\ref apop_matrix_map_sum
\li\ref apop_vector_apply
\li\ref apop_vector_map
\li\ref apop_vector_map_sum



\section  matrixmathtwo  Basic Math

\li\ref apop_vector_exp : exponentiate every element of a vector
\li\ref apop_vector_log : take the natural log of every element of a vector
\li\ref apop_vector_log10 : take the log (base 10) of every element of a vector
\li\ref apop_vector_distance : find the distance between two vectors via various metrics
\li\ref apop_vector_normalize : scale/shift a matrix to have mean zero, sum to one, have a range of exactly $[0, 1]$, et cetera
\li\ref apop_vector_entropy : calculate the entropy of a vector of frequencies or probabilities

    See also:

\li <tt>int gsl_matrix_add (gsl_matrix * a, const gsl_matrix * b)</tt>
\li <tt>int gsl_matrix_sub (gsl_matrix * a, const gsl_matrix * b)</tt>
\li <tt>int gsl_matrix_mul_elements (gsl_matrix * a, const gsl_matrix * b)</tt>
\li <tt>int gsl_matrix_div_elements (gsl_matrix * a, const gsl_matrix * b)</tt>
\li <tt>int gsl_matrix_scale (gsl_matrix * a, const double x)</tt>
\li <tt>int gsl_matrix_add_constant (gsl_matrix * a, const double x)</tt>
\li <tt>gsl_vector_add (gsl_vector * a, const gsl_vector * b)</tt>
\li <tt>gsl_vector_sub (gsl_vector * a, const gsl_vector * b)</tt>
\li <tt>gsl_vector_mul (gsl_vector * a, const gsl_vector * b)</tt>
\li <tt>gsl_vector_div (gsl_vector * a, const gsl_vector * b)</tt>
\li <tt>gsl_vector_scale (gsl_vector * a, const double x)</tt>
\li <tt>gsl_vector_add_constant (gsl_vector * a, const double x)</tt>

            
\section  matrixmath  Matrix math

\li\ref apop_dot : matrix \f$\cdot\f$ matrix, matrix \f$\cdot\f$ vector, or vector \f$\cdot\f$ matrix
\li\ref apop_matrix_determinant
\li\ref apop_matrix_inverse
\li\ref apop_det_and_inv : find determinant and inverse at the same time

See the GSL documentation for voluminous further options.


\section  sumstats  Summary stats

\li\ref apop_data_summarize
\li\ref apop_vector_moving_average
\li\ref apop_vector_percentiles
\li\ref apop_vector_bounded

        See also:

\li <tt>double gsl_matrix_max (const gsl_matrix * m)</tt>
\li <tt>double gsl_matrix_min (const gsl_matrix * m)</tt>
\li <tt>void gsl_matrix_minmax (const gsl_matrix * m, double * min_out, double * max_out)</tt>
\li <tt>void gsl_matrix_max_index (const gsl_matrix * m, size_t * imax, size_t * jmax)</tt>
\li <tt>void gsl_matrix_min_index (const gsl_matrix * m, size_t * imin, size_t * jmin)</tt>
\li <tt>void gsl_matrix_minmax_index (const gsl_matrix * m, size_t * imin, size_t * jmin, size_t * imax, size_t * jmax)</tt>
\li <tt>gsl_vector_max (const gsl_vector * v)</tt>
\li <tt>gsl_vector_min (const gsl_vector * v)</tt>
\li <tt>gsl_vector_minmax (const gsl_vector * v, double * min_out, double * max_out)</tt>
\li <tt>gsl_vector_max_index (const gsl_vector * v)</tt>
\li <tt>gsl_vector_min_index (const gsl_vector * v)</tt>
\li <tt>gsl_vector_minmax_index (const gsl_vector * v, size_t * imin, size_t * imax)</tt>


\section  moments  Moments

For most of these, you can add a weights vector for weighted mean/var/cov/..., such as
<tt>apop_vector_mean(d->vector, .weights=d->weights)</tt>

\li\ref apop_mean : the first three with short names operate on a vector.
\li\ref apop_sum
\li\ref apop_var
\li\ref apop_matrix_sum
\li\ref apop_data_correlation
\li\ref apop_data_covariance
\li\ref apop_data_summarize
\li\ref apop_matrix_mean
\li\ref apop_matrix_mean_and_var
\li\ref apop_vector_correlation
\li\ref apop_vector_cov
\li\ref apop_vector_kurtosis
\li\ref apop_vector_kurtosis_pop 
\li\ref apop_vector_mean
\li\ref apop_vector_skew
\li\ref apop_vector_skew_pop
\li\ref apop_vector_sum
\li\ref apop_vector_var
\li\ref apop_vector_var_m 


\section convsec   Conversion among types

There are no functions provided to convert from \ref apop_data to the constituent
elements, because you don't need a function.

If you need an individual element, you can use its pointer to retrieve it:

\code
apop_data *d = apop_query_to_mixed_data("vmmw", "select result, age, "
                                     "income, replicate_weight from data");
double avg_result = apop_vector_mean(d->vector, .weights=d->weights);
\endcode

In the other direction, you can use compound literals to wrap an \ref apop_data struct
around a loose vector or matrix:

\code
//Given:
gsl_vector *v;
gsl_matrix *m;

// Then this form wraps the elements into apop_data structs. Note that
// these are not pointers: they're automatically allocated and therefore
// the extra memory use for the wrapper is cleaned up on exit from scope.

apop_data *dv = &(apop_data){.vector=v}; 
apop_data *dm = &(apop_data){.matrix=m};

apop_data *v_dot_m = apop_dot(dv, dm);

//Here is a macro to hide C's ugliness:
#define As_data(...) (&(apop_data){__VA_ARGS__})

apop_data *v_dot_m2 = apop_dot(As_data(.vector=v), As_data(.matrix=m));

//The wrapped object is an automatically-allocated structure pointing to the
//original data. If it needs to persist or be separate from the original,
//make a copy:
apop_data *dm_copy = apop_data_copy(As_data(.vector=v, .matrix=m));
\endcode

\li\ref apop_array_to_vector : <tt>double*</tt>\f$\to\f$ <tt>gsl_vector</tt>
\li\ref apop_data_fill : <tt>double*</tt>\f$\to\f$  \ref apop_data.
\li\ref apop_data_falloc : macro to allocate and fill a \ref apop_data set.
\li\ref apop_text_to_data : delimited text file\f$\to\f$ \ref apop_data
\li\ref apop_text_to_db : delimited text file\f$\to\f$ database
\li\ref apop_vector_to_matrix


\section names   Name handling

If you generate your data set via \ref apop_text_to_data or from the database via
\ref apop_query_to_data (or \ref apop_query_to_text or \ref apop_query_to_mixed_data)
then column names appear as expected.  Set <tt>apop_opts.db_name_column</tt> to the
name of a column in your query result to use that column name for row names.

Sample uses, given \ref apop_data set <tt>d</tt>:

\code
int row_name_count = d->names->rowct
int col_name_count = d->names->colct
int text_name_count = d->names->textct

//Manually add names in sequence:
apop_name_add(d->names, "the vector", 'v');
apop_name_add(d->names, "row 0", 'r');
apop_name_add(d->names, "row 1", 'r');
apop_name_add(d->names, "row 2", 'r');
apop_name_add(d->names, "numeric column 0", 'c');
apop_name_add(d->names, "text column 0", 't');
apop_name_add(d->names, "The name of the data set.", 'h');

//or append several names at once
apop_data_add_names(d, 'c', "numeric column 1", "numeric column 2", "numeric column 3");

//point to element i from the row/col/text names:

char *rowname_i = d->names->row[i];
char *colname_i = d->names->col[i];
char *textname_i = d->names->text[i];

//The vector also has a name:
char *vname = d->names->vector;
\endcode

\li\ref apop_name_add : add one name
\li\ref apop_data_add_names : add a sequence of names at once
\li\ref apop_name_stack : copy the contents of one name list to another
\li\ref apop_name_find : find the row/col number for a given name.
\li\ref apop_name_print : print the \ref apop_name struct, for diagnostic purposes.


\section textsec   Text data

The \ref apop_data set includes a grid of strings, <tt>text</tt>, for holding text data. 

Text should be encoded in UTF-8. ASCII is a subset of UTF-8, so that's OK too.

There are a few simple forms for handling the \c text element of an \c apop_data set, which handle the tedium of memory-handling for you.

\li Use \ref apop_text_alloc to allocate the block of text. It is actually a realloc function, which you can use to resize an existing block without leaks.
\li Use \ref apop_text_set to write text elements. It replaces any existing text in the given slot without memory leaks.
\li The number of rows of text data in <tt>tdata</tt> is
<tt>tdata->textsize[0]</tt>; 
the number of columns is <tt>tdata->textsize[1]</tt>.
\li Refer to individual elements using the usual 2-D array notation, <tt>tdata->text[row][col]</tt>.
\li <tt>x[0]</tt> can always be written as <tt>*x</tt>, which may save some typing. The number of rows is <tt>*tdata->textsize</tt>. If you have a single column of text data (i.e., all data is in column zero), then item \c i is <tt>*tdata->text[i]</tt>. If you know you have exactly one cell of text, then its value is <tt>**tdata->text</tt>.
\li After \ref apop_text_alloc, all elements are the empty string <tt>""</tt>, which
you can check via <tt>if (!strlen(dataset->text[i][j])) printf("<blank>")</tt> or
<tt>if (!*dataset->text[i][j]) printf("<blank>")</tt>. For the sake of efficiency
when dealing with large, sparse data sets, all blank cells point to <em>the same</em>
static empty string, meaning that freeing cells must be done with care. Your best bet
is to rely on \ref apop_text_set, \ref apop_text_alloc, and \ref apop_text_free to do
the memory management for you.

Here is a sample program that uses these forms, plus a few text-handling functions.

\include eg/text_demo.c

\li\ref apop_data_transpose() : also transposes the text data. Say that you use
<tt>dataset = apop_query_to_text("select onecolumn from data");</tt> then you have a
sequence of strings, <tt>d->text[0][0], d->text[1][0], </tt>.... After <tt>apop_data
*dt = apop_data_transpose(dataset)</tt>, you will have a single list of strings,
<tt>dt->text[0]</tt>, which is often useful as input to list-of-strings handling
functions.

\li\ref apop_query_to_text
\li\ref apop_text_alloc : allocate or resize the text part of an \ref apop_data set.
\li\ref apop_text_set : replace a single cell of the text grid with new text.
\li\ref apop_text_paste : convert a table of strings into one long string.
\li\ref apop_text_unique_elements : get a sorted list of unique elements for one column of text.
\li\ref apop_text_free : you may never need this, because \ref apop_data_free calls it.
\li\ref apop_regex : friendlier front-end for POSIX-standard regular expression
            searching and pulling matches into an \ref apop_data set.

\subsection fact   Generating factors

\em Factor is jargon for a numbered category. Number-crunching programs prefer integers over text, so we need a function to produce a one-to-one mapping from text categories into numeric factors. 

A \em dummy is a variable that is either one or zero, depending on membership in a given
group. Some methods (typically when the variable is an input or independent variable
in a regression) prefer dummies; some methods (typically for outcome or dependent
variables) prefer factors. The functions that generate factors and dummies will add
an informational page to your \ref apop_data set with a name like <tt>\<categories
for your_column\></tt> listing the conversion from the artificial numeric factor to
the original data. Use \ref apop_data_get_factor_names to get a pointer to that page.

You can use the factor table to translate from numeric categories back to text (though
you probably have the original text column in your data anyway).

Having the factor list in an auxiliary table makes it easy to ensure that multiple
\ref apop_data sets use the same single categorization scheme. Generate factors in the
first set, then copy the factor list to the second, then run \ref apop_data_to_factors
on the second:

\code
apop_data_to_factors(d1);
d2->more = apop_data_copy(apop_data_get_factor_names(d1));
apop_data_to_factors(d2);
\endcode


\li\ref apop_data_to_dummies
\li\ref apop_data_to_factors
\li\ref apop_data_get_factor_names
\li\ref apop_text_unique_elements
\li\ref apop_vector_unique_elements
*/

/** \page dbs Databases

These are convenience functions to handle interaction with SQLite or mySQL/mariaDB. They open one and only one database, and handle most of the interaction therewith for you.

You will probably first use \ref apop_text_to_db to pull data into the database, then \ref apop_query to clean the data in the database, and finally \ref apop_query_to_data to pull some subset of the data out for analysis.

\li In all cases, your query may be in <tt>printf</tt> form. For example:
\code
char tabname[] = "demographics";
char colname[] = "heights";
int min_height = 175;
apop_query("select %s from %s where %s > %i", colname, tabname, colname, min_height);
\endcode


See the \ref db_moments section below for not-SQL-standard math functions that you can
use when sending queries from Apophenia, such as \c pow, \c stddev, or \c sqrt.

\li \ref apop_text_to_db : Read a text file on disk into the database. Data analysis projects often start with a call to this.
\li \ref apop_data_print : If you include the argument <tt>.output_type='d'</tt>, this prints your \ref apop_data set to the database.
\li \ref apop_query : Manipulate the database, return nothing (e.g., insert rows or create table).
\li \ref apop_db_open : Optional, for when you want to use a database on disk.
\li \ref apop_db_close : A useful (and in some cases, optional) companion to \ref apop_db_open.
\li \ref apop_table_exists : Check to make sure you aren't reinventing or destroying data. Also, a clean way to drop a table.

\li Apophenia reserves the right to insert temp tables into the opened database. They
will all have names beginning with <tt>apop_</tt>, so the reader is advised to not
generate tables with such names, and is free to ignore or delete any such tables that
turn up.
\li If you need to deal with two databases, use SQL's <a
href="https://sqlite.org/lang_attach.html"><tt>attach database</tt></a>. By default
with SQLite, Apophenia opens an in-memory database handle. It is a sensible workflow to
use the faster in-memory database as the primary database, and then attach an on-disk database
to read in data and write final output tables.

\section edftd Extracting data from the database

\li\ref apop_db_to_crosstab : take up to three columns in the database (row, column, value) and produce a table of values.
\li\ref apop_query_to_data
\li\ref apop_query_to_float
\li\ref apop_query_to_mixed_data
\li\ref apop_query_to_text
\li\ref apop_query_to_vector

\section wdttd Writing data to the database

See the print functions at \ref Legi. E.g.

\code
apop_data_print(yourdata, .output_type='d', .output_name="dbtab");
\endcode

\section cmdline Command-line utilities

A few functions have proven to be useful enough to be worth breaking out into their own programs, for use in scripts or other data analysis from the command line:

\li The \c apop_text_to_db command line utility is a wrapper for the \ref apop_text_to_db command.
\li The \c apop_db_to_crosstab function is a wrapper for the \ref apop_db_to_crosstab function.

\section db_moments Database moments (plus pow()!)

SQLite lets users define new functions for use in queries, and Apophenia uses this facility to define a few common functions.

\li <tt>select ran() from table</tt> will produce a new random number between zero and one for every row of the input table, using \c gsl_rng_uniform. 

\li The SQL standard includes the <tt>count(x)</tt> and <tt>avg(x)</tt> aggregators,
but statisticians are usually interested in higher moments as well---at least the
variance. Therefore, SQL queries using the Apophenia library may include any of these moments:

\code
select count(x), stddev(x), avg(x), var(x), variance(x), skew(x), kurt(x), kurtosis(x),
std(x), stddev_samp(x), stddev_pop(x), var_samp(x), var_pop(x)
from table
group by whatever
\endcode

<tt>var</tt> and <tt>variance</tt>; <tt>kurt</tt> and <tt>kurtosis</tt> do the same thing. Choose the one that sounds better to you. <tt>var</tt>, <tt>var_samp</tt>, <tt>stddev</tt> and <tt>stddev_samp</tt> give sample variance/standard deviation; <tt>variance</tt>, <tt>var_pop</tt>, <tt>std</tt> and <tt>stddev_pop</tt> give population standard deviation.  The plethora of variants are for mySQL compatibility. Kurtosis is the fourth central moment by itself, not adjusted by subtracting three or dividing by variance squared.

\li The  var/skew/kurtosis functions calculate sample moments. If you want the population moment of the variance/skew, multiply the result by (n-1)/n . The equation for the unbiased sample kurtosis as calculated in <a href="http://modelingwithdata.org/pdfs/moments.pdf">Appendix M of <em>Modeling with Data</em></a> is not quite as simple.

\li Also provided: wrapper functions for standard math library
functions---<tt>sqrt(x)</tt>, <tt>pow(x,y)</tt>, <tt>exp(x)</tt>, <tt>log(x)</tt>,
and trig functions. They call the standard math library function of the same name
to calculate \f$\sqrt{x}\f$, \f$x^y\f$, \f$e^x\f$, \f$\ln(x)\f$, \f$\sin(x)\f$,
\f$\arcsin(x)\f$, et cetera. For example:

\code
select sqrt(x), pow(x,0.5), exp(x), log(x), log10(x),
    sin(x), cos(x), tan(x), asin(x), acos(x), atan(x)
from table
\endcode

\li The <tt>ran()</tt> function calls <tt>gsl_rng_uniform</tt> to produce a uniform
draw between zero and one. It uses the stock of RNGs from \ref apop_rng_get_thread.

Here is a test script using many of the above.

\include db_fns.c
*/


/** \page modelsec Models
See \ref gentle_model for an overview of the intent and basic use of the \ref apop_model struct.

This segment goes into greater detail on the use of existing \ref apop_model objects.
If you need to write a new model, see \ref modeldetails.

The \c estimate function will estimate the parameters of your model. Just prep the data, select a model, and produce an estimate:

\code
    apop_data *data = apop_query_to_data("select outcome, in1, in2, in3 from dataset");
    apop_model *the_estimate = apop_estimate(data, apop_probit);
    apop_model_print(the_estimate, NULL);
\endcode

Along the way to estimating the parameters, most models also find covariance estimates for
the parameters, calculate statistics like log likelihood, and so on, which the final print statement will show.

The <tt>apop_probit</tt> model that ships with Apophenia is unparameterized:
<tt>apop_probit->parameters==NULL</tt>. The output from the estimation,
<tt>the_estimate</tt>, has the same form as <tt>apop_probit</tt>, but
<tt>the_estimate->parameters</tt> has a meaningful value.

Apophenia ships with many well-known models for your immediate use, including
probability distributions, such as the \ref apop_normal, \ref apop_poisson, or \ref
apop_beta models. The data is assumed to have been drawn from a given distribution and
the question is only what distributional parameters best fit. For example, given that
the data is Normally distributed, find the mean and variance via:
<tt>apop_estimate(your_data, apop_normal)</tt>.

There are also linear models like \ref apop_ols, \ref apop_probit, and \ref apop_logit. As in the example, they are on equal footing with the distributions, so nothing keeps you from making random draws from an estimated linear model.

  \li If you send a data set with the \c weights vector filled, \ref apop_ols estimates Weighted OLS.
  \li If the dependent variable has more than two categories, the \ref apop_probit and
\ref apop_logit models estimate a multinomial logit or probit.
  \li There are separate \ref apop_normal and \ref apop_multivariate_normal functions
because the parameter formats are slightly different: the univariate Normal keeps both
\f$\mu\f$ and \f$\sigma\f$ in the vector element of the parameters; the multivariate version uses the
vector for the vector of means and the matrix for the \f$\Sigma\f$ matrix. The univariate is so
heavily used that it merits a special-case model.

See the \ref models page for a list of with a list of models shipped with Apophenia,
including popular favorites like \ref apop_beta, \ref apop_binomial, \ref apop_iv
(instrumental variables), \ref apop_kernel_density, \ref apop_loess, \ref apop_lognormal,
\ref apop_pmf (see \ref histosec below), and \ref apop_poisson.

Simulation models seem to not fit this form, but you will see below that if you can write an objective function for the \c p method of the model, you can use the above tools. Notably, you can estimate parameters via maximum likelihood and then give confidence intervals around those parameters.

\par More estimation output

A call to \ref apop_estimate produces more than just the estimated parameters. Most models will
produce any of a covariance matrix, some hypothesis tests, a list of expected values, log
likelihood, AIC, AIC_c, BIC, et cetera.

If you don't want all that, adding to your model an \ref apop_parts_wanted_settings
group with its default values (see below on settings groups) signals to the model
that you want only the parameters and to not waste possibly significant CPU time
on covariances, expected values, et cetera. See the \ref apop_parts_wanted_settings
documentation for examples and further refinements.

You will also find in the \ref apop_model returned by \ref apop_estimate:

\li The actual parameter estimates are in an \ref apop_data set at \c your_model->parameters.
\li A pointer to the \ref apop_data set used for estimation, named \c data.
\li Scalar statistics of the model listed in the output model's \c info group, which can
be retrieved via a form like
\code
apop_data_get(your_model->info, .rowname="log likelihood");
//or
apop_data_get(your_model->info, .rowname="AIC");
\endcode
\li Covariances of the parameters as a page appended to the parameters; retrieve via
\code
apop_data *cov = apop_data_get_page(your_model->parameters, "<Covariance>");
\endcode
\li If the model calculates it, the table of expected values (typically including
expected value, actual value, and residual) is a page stapled to the main info
page. This is mostly for regression-type models. Retrieve via:
\code
apop_data *predict = apop_data_get_page(your_model->info, "<Predicted>");
\endcode


\par Post-estimation uses

But we expect much more from a model than estimating parameters from data.  

Continuing the above example where we got an estimated Probit model named \c the_estimate, we can interrogate the estimate in various familiar ways. In each of the following examples, the model object holds enough information that the generic function being called can do its work:

\code
apop_data *expected_value = apop_predict(NULL, the_estimate);

double density_under =  apop_cdf(expected_value, the_estimate);

apop_data *draws = apop_model_draws(the_estimate, .count=1000);
\endcode

\par Additional settings

But some models have to break uniformity, like how a histogram has a list of bins that
makes no sense for a Normal distribution. Modifiable characteristics of the model
that are not treated as parameters are held in <em>settings groups</em>, which you
will occasionally need to set to modify how a model is handled or estimated. The
most common example would be for maximum likelihood, eg.

\code
//Probit uses MLE. Redo the estimation using Newton's Method
Apop_settings_add_group(the_estimate, apop_mle, .verbose='y', 
                        .tolerance=1e-4, .method="Newton");
apop_model *re_est = apop_estimate(data, the_estimate);
\endcode

See below for more details on using settings groups.

To clarify the distinction between parameters and settings, note that parameters are
estimated from the data, often via a maximum likelihood search. In an ML search,
the method of search, the number of bins in a histogram, or the number of steps in a
simulation would be held fixed as the search iterates over possible parameters (and
if these settings do change, then that is a meta-model that could be encapsulated into another
\ref apop_model). As a consequence, parameters are always numeric, while settings may
be any type.

\subpage dataones

\section modelparameterization  Parameterizing or initializing a model

The models that ship with Apophenia have the requisite procedures for estimation,
making draws, and so on, but have <tt>parameters==NULL</tt> and <tt>settings==NULL</tt>. The
model is thus, for many purposes, incomplete, and you will need to take some action to
complete the model. As per the examples to follow, there are several possibilities:

  \li Estimate it! Almost all models can be sent with a data set as an argument to the
<tt>apop_estimate</tt> function. The input model is unchanged, but the output model
has parameters and settings in place.
  \li If your model has a fixed number of numeric parameters, then you can set them with
\ref apop_model_set_parameters.
  \li If your model has a variable number of parameters, you can directly set the 
\c parameters element via \ref apop_data_falloc. For most purposes, you will also need to
set the \c msize1, \c msize2, \c vsize, and \c dsize elements to the size you want. See
the example below.
  \li Some models have disparate, non-numeric settings rather than a simple matrix of
parameters. For example, an kernel density estimate needs a model as a kernel and a
base data set, which can be set via \ref apop_model_copy_set.

Here is an example that shows the options for parameterizing a model. After each
parameterization, 20 draws are made and written to a file named draws-[modelname].

\include ../eg/parameterization.c


If you need to write a new model, see \ref modeldetails.

\section transformsec Transforming models

Given that all models take the same \ref apop_model form, we can write functions that take model objects as input and product them as outputs.

In this example, Nature generated data using a mixture of three Poisson distributions,
with \f$\lambda=2.8\f$, \f$2.0\f$, and \f$1.3\f$. The resulting model is generated
using \ref apop_model_mixture.  Not knowing the true distribution, the analyst wrote
down a model with a truncated Normal(2, 1) prior (generated by sending the stock
\ref apop_normal model to the data-space constraint function \ref apop_dconstrain)
describing the parameter of a Poisson likelihood model. 
The \ref apop_update function takes three arguments: the data set, which comes from
draws from the mixture, the prior, and the likelihood. It produces an output model
which, in this case, is a PMF describing a distribution over \f$\lambda\f$, because
a truncated Normal and a Poisson are not conjugate distributions. Knowing that it is
a PMF, the <tt>->data</tt> element holds a set of draws from the posterior.

The analyst would like to present an approximation to the posterior in a simpler form,
and so finds the parameters \f$\mu\f$ and \f$\sigma\f$ of the Normal distribution that
is closest to that posterior.

Here is a program---almost a single line of code---that builds the final posterior model from the subcomponents, including draws from Nature and the analyst's prior and likelihood:

\include ../eg/transform.c

\section mathmethods Model methods

\li\ref apop_estimate : estimate the parameters of the model with data.
\li\ref apop_predict : the expected value function.
\li\ref apop_draw : random draws from an estimated model.
\li\ref apop_p : the probability of a given data set given the model.
\li\ref apop_log_likelihood : the log of \ref apop_p
\li\ref apop_score : the derivative of \ref apop_log_likelihood
\li\ref apop_model_print : write to screen, file, or database
\li\ref apop_model_copy : duplicate a model
\li\ref apop_model_set_parameters : Models ship with no parameters set. Use this to convert a Normal(μ, σ) with unknown μ and σ into a Normal(0, 1), for example.
\li\ref apop_model_free
\li\ref apop_model_clear , \ref apop_prep : remove the parameters from a parameterized model. Used infrequently.
\li\ref apop_model_draws : many random draws from an estimated model.


\section Update Filtering & updating

The model structure makes it
easy to generate new models that are variants of prior models. Bayesian updating,
for example, takes in one \ref apop_model that we call the prior, one \ref apop_model
that we call a likelihood, and outputs an \ref apop_model that we call the
posterior. One can produce complex models using simpler transformations as well. For example, to generate
a one-parameter Normal(μ, 1) given the code for for a Normal(μ, σ):

\code
apop_model *N_sigma1 = apop_model_fix_params(apop_model_set_parameters(apop_normal, NAN, 1));
\endcode

This can be used anywhere the original Normal distribution can be. If we need to truncate the distribution in the data space:

\code
//The constraint function.
double over_zero(apop_data *in, apop_model *m){
    return apop_data_get(in) > 0;
}

apop_model *trunc = apop_model_dconstrain(.base_model=N_sigma1,
                                          .constraint=over_zero);
\endcode

Chaining together simpler transformations is an easy method to produce 
models of arbitrary detail.

\li\ref apop_update : Bayesian updating
\li\ref apop_model_coordinate_transform : apply an invertible transformation to the data space
\li\ref apop_model_dconstrain : constrain the data space of a model to a subspace. E.g., truncate a Normal distribution so \f$x>0\f$.
\li\ref apop_model_fix_params : hold some parameters constant
\li\ref apop_model_mixture : a linear combination of models
\li\ref apop_model_cross : If \f$(p_1, p_2)\f$ has a Normal distribution and \f$p_3\f$ has an independent Poisson distribution, then \f$(p_1, p_2, p_3)\f$ has an <tt>apop_model_stack(apop_normal, apop_poisson)</tt> distribution.
\li\ref apop_model_dcompose : use the output of one model as a data set for another


\section modelsettings Settings groups

[For info on specific settings groups and their contents and use, see the \ref settings page.]


Describing a statistical, agent-based, social, or physical model in a standardized form is difficult because every model has significantly different settings. E.g., an MLE requires a method of search (conjugate gradient, simplex, simulated annealing), and a histogram needs the number of bins to be filled with data.

So, the \ref apop_model includes a single list which can hold an arbitrary number of groups of settings, like the search specifications for finding the maximum likelihood, a histogram for making random draws, and options about the model type.

Settings groups are automatically initialized with default values when
needed. If the defaults do no harm, then you don't need to think about
these settings groups at all.

Here is an example where a settings group is worth tweaking: the \ref apop_parts_wanted_settings group indicates which parts
of the auxiliary data you want. 


\code
1 apop_model *m = apop_model_copy(apop_ols);
2 Apop_settings_add_group(m, apop_parts_wanted, .covariance='y');
3 apop_model *est = apop_estimate(data, m);
\endcode


Line one establishes the baseline form of the model. Line two adds a settings group
of type \ref apop_parts_wanted_settings to the model. By default other auxiliary items, like the expected values, are set to \c 'n' when using this group, so this specifies that we want covariance and only covariance. Having stated our preferences, line three does the estimation we want.

Notice that the \c _settings ending to the settings group's name isn't written---macros
make it happen.  The remaining arguments to \c Apop_settings_add_group (if any) follow
the \ref designated syntax of the form <tt>.setting=value</tt>.

There is an \ref apop_model_copy_set macro that adds a settings group when it is first copied, joining up lines one and two above:

\code
apop_model *m = apop_model_copy_set(apop_ols, apop_parts_wanted, .covariance='y');
\endcode

Settings groups are copied with the model, which facilitates chaining
estimations. Continuing the above example, you could re-estimate to get the predicted
values and covariance via:

\code
Apop_settings_set(est, apop_parts_wanted, predicted, 'y');
apop_model *est2 = apop_estimate(data, est);
\endcode

\li \ref Apop_settings_set, for modifying a single setting, doesn't use the designated initializers format.
\li  The \ref settings page lists the settings structures included in Apophenia and their use.
\li Because the settings groups are buried within the model, debugging them can be a
pain. Here is a documented macro for \c gdb that will help you pull a settings group out of a 
model for your inspection, to cut and paste into your \c .gdbinit. It shouldn't be too difficult to modify this macro for other debuggers.

\code
define get_group
    set $group = ($arg1_settings *) apop_settings_get_grp( $arg0, "$arg1", 0 )
    p *$group
end
document get_group 
Gets a settings group from a model.
Give the model name and the name of the group, like
get_group my_model apop_mle 
and I will set a gdb variable named $group that points to that model, 
which you can use like any other pointer. For example, print the contents with
p *$group
The contents of $group are printed to the screen as visible output to this macro.
end 
\endcode

For using a model, that's all of what you need to know. For details on writing a new settings group, see \ref settingswriting .

\li\ref Apop_settings_add_group
\li\ref Apop_settings_set
\li\ref Apop_settings_get : get a single element from a settings group.
\li\ref Apop_settings_get_group : get the whole settings group.
*/

/** \page dataones Data format for regression-type models

Regression-type estimations typically require a constant column. That is, the 0th
column of the data is a constant (one), so the parameter \f$\beta_0\f$ is
slightly special in corresponding to a constant rather than a variable.

On the other hand, there are some estimations that do not use the constant column.

Some stats packages implicitly assume a constant column, which the user never
sees. This violates the principle of transparency upon which Apophenia is based,
and is generally annoying.  Given a data matrix \f$X\f$ with the estimated parameters
\f$\beta\f$, if the model asserts that the product \f$X\beta\f$ has meaning, then you
should be able to easily calculate that product. With a ones column, a dot product is one line:
<tt>apop_dot(x, your_est->parameters)</tt>; without a ones column, one would basically
have to construct one (using \c gsl_matrix_set_all and \c apop_data_stack).

Each regression-type estimation has one dependent variable and several independent. In
the end, we want the dependent variable to be in the vector element. Removing
a column from a <tt>gsl_matrix</tt> and adjusting all subsequent columns is relatively
difficult, because (like most structs built with the aim of very efficient processing) the
struct depends on an equal spacing in memory between each element.

<em> The automatic case</em>

We can resolve both the need for a ones column and for having the dependent column in
the vector at the same time. Given a data set with no vector element and the dependent
variable in the first column of the matrix, we can copy the dependent variable into
the vector and then replace the first column of the matrix with ones. The result fits
all of the above expectations.

You as a user merely have to send in a \c apop_data set with \c NULL vector and a dependent
column in the first column. If the data is coming from the database, then the query
is natural:

\code
apop_data *regression_data = apop_query_to_data("select depvar, indyvar1, indyvar2, indyvar3 from dataset");
apop_model_print(apop_estimate(regression_data, apop_ols), NULL);
\endcode

<em> The already-prepped case</em>

If your data has a vector element, then the prep routines won't change anything.
If you don't want to use a constant column, or your data has already been prepped by
an estimation, then this is what you want.

\code
apop_data *regression_data = apop_query_to_mixed_data("vmmm", "select depvar, indyvar1, indvar2, indvar3 from dataset");
apop_model_print(apop_estimate(regression_data, apop_logit), NULL);
\endcode
*/


/** \page testpage Tests & diagnostics

Apophenia provides a few functions for doing the more common hypothesis tests, and
enough tools associated with the \ref apop_model struct that confidence-interval type
tests can be constructed for any arbitrary distribution.

If you are producing a statistic that you know has a common form, like a central limit
theorem tells you that your statistic is Normally distributed, then the convenience
function \ref apop_test will do the final lookup step of checking where your statistic
lies on your chosen distribution.

In especially common cases, like the parameters from an OLS regression,
the commonly-associated \f$t\f$ test is included as part of the estimation
output, typically as a row in the \c info element of the output \ref apop_model.

Some tests, like ANOVA, produce a statistic using a specialized procedure, so
Apophenia includes some functions, like \ref apop_test_anova_independence and \ref
apop_test_kolmogorov, to produce the statistic and look up its significance level.

\li\ref apop_test
\li\ref apop_paired_t_test
\li\ref apop_f_test
\li\ref apop_t_test
\li\ref apop_test_anova_independence
\li\ref apop_test_fisher_exact
\li\ref apop_test_kolmogorov
\li\ref apop_estimate_coefficient_of_determination
\li\ref apop_estimate_r_squared

These are all situations where the statistic in question is known to have a textbook
distribution. \ref gentle_testing gives discussion and an example for the case where
the distribution of the statistic must be derived via Monte Carlo methods.

See also these Monte Carlo methods:

\li\ref apop_bootstrap_cov
\li\ref apop_jackknife_cov

To give another example of testing, here is a function that used to be a part of
Apophenia, but seemed a bit out of place. Here it is as a sample:

\code
// Input: any vector. Output: 1 - the p-value for a chi-squared
// test to answer the question, "with what confidence can I reject the
// hypothesis that the variance of my data is zero?"

double apop_test_chi_squared_var_not_zero(const gsl_vector *in){
    Apop_stopif(!in, return NAN, 0, "input vector is NULL. Doing nothing.");
    double sum=0;
    gsl_vector	*normed = apop_vector_normalize((gsl_vector *)in, .normalization_type='s');
    gsl_vector_mul(normed, normed);
    for(size_t i=0;i< normed->size; 
        sum +=gsl_vector_get(normed,i++));
    gsl_vector_free(normed);
    return gsl_cdf_chisq_P(sum, in->size); 
}
\endcode

Or, consider the Rao statistic, 
\f${\partial\over \partial\beta}\log L(\beta)'I^{-1}(\beta){\partial\over \partial\beta}\log L(\beta)\f$
where \f$L\f$ is your model's likelihood function and \f$I\f$ its information matrix. In code:

\code
apop_data * infoinv = apop_model_numerical_covariance(data, your_model);
apop_data * score = &(apop_data*){.vector=apop_numerical_gradient(data, your_model)};
apop_data * stat = apop_dot(apop_dot(score, infoinv), score);
\endcode

Given the correct assumptions, this is \f$\sim \chi^2_m\f$, where \f$m\f$ is the dimension of \f$\beta\f$, so the odds of a Type I error given the model is:

\code
double p_value = apop_test(stat, "chi squared", beta->size);
\endcode
*/

/** \page Histosec Empirical distributions and PMFs (probability mass functions)

The \ref apop_pmf model wraps a \ref apop_data set so it can be read as an empirical
model, with a likelihood function (equal to the associated weight for observed
values and zero for unobserved values), a random number generator (which 
simply makes weighted random draws from the data), and so on.  Setting it up is a
model estimation from data like any other, done via \ref apop_estimate(\c your_data,
\ref apop_pmf).

You have the option of cleaning up the data before turning it into a PMF. For example...

\code
apop_data_pmf_compress(your_data);          //remove duplicates
apop_data_sort(your_data);
apop_vector_normalize(your_data->weights);  //weights sum to one
apop_model *a_pmf = apop_estimate(your_data, apop_pmf);
\endcode

These are largely optional.

\li The CDF is calculated based on the percent of the weights between the zeroth row of the PMF
and the row specified. This generally makes more sense after \ref apop_data_sort.
\li Compression produces a corresponding improvement in efficiency when calculating
CDFs, but is otherwise not necessary.
\li Sorting or normalizing is not necessary for making draws or getting a likelihood or log likelihood.

It is the \c weights vector that holds the density represented by each row; the rest of the row represents the coordinates of that density. If the input data set has no \c weights segment, then I assume that all rows have equal weight.

Most models have a \c parameters \ref apop_data set that is filled when you call \ref
apop_estimate. For a PMF model, the \c parameters are \c NULL, and the \c data itself
is used for calculation. Therefore, modifying the data post-estimation can break some
internal settings set during estimation. If you modify the data, throw away any existing
PMFs (via \ref apop_model_free) and re-estimate a new one.


Using \ref apop_data_pmf_compress puts the data into one bin for each unique value in
the data set. 
You may instead want bins of fixed with, in the style of a histogram, which you can get via 
\ref apop_data_to_bins. It requires a bin specification. If you send a \c NULL binspec,
then the offset is zero and the bin size is big enough to ensure that there are
\f$\sqrt{N}\f$ bins from minimum to maximum. There are other preferred formulæ for
bin widths that minimize MSE, which might be added as a future extension.  The binspec
will be added as a page to the data set, named <tt>"<binspec>"</tt>. See the \ref
apop_data_to_bins documentation on how to write a custom bin spec.

\section histocompare Comparing histograms

There are a few ways of testing the claim that one distribution equals another, typically an empirical PMF versus a smooth theoretical distribution. In both cases, you will need two distributions based on the same binspec. 

For example, if you do not have a prior binspec in mind, then you can use the one generated by the first call to the histogram binning function to make sure that the second data set is in sync:

\code
apop_data_to_bins(first_set, NULL);
apop_data_to_bins(second_set, apop_data_get_page(first_set, "<binspec>"));
\endcode

You can use \ref apop_test_kolmogorov or \ref apop_histograms_test_goodness_of_fit to generate the appropriate statistics from the pairs of bins.

Kernel density estimation will produce a smoothed PDF. See \ref apop_kernel_density for details.
Or, use \ref apop_vector_moving_average for a simpler smoothing method.


\li\ref apop_data_pmf_compress() : merge together redundant rows in a data set before calling 
                \ref apop_estimate(\c your_data, \ref apop_pmf); optional.
\li\ref apop_vector_moving_average() : smooth a vector (e.g., your_pmf->data->weights) via moving average.
\li\ref apop_histograms_test_goodness_of_fit() : goodness-of-fit via \f$\chi^2\f$ statistic
\li\ref apop_test_kolmogorov() : goodness-of-fit via Kolmogorov-Smirnov statistic
*/

/** \page maxipage Optimization

This section includes some notes on the maximum likelihood routine. As in the section
on writing models above, if a model has a \c p or \c log_likelihood method but no \c
estimate method, then calling \c apop_estimate(your_data, your_model) executes the
default estimation routine of maximum likelihood.

If you are a not a statistician, then there are a few things you will need to keep in
mind:

\li Physicists, pure mathematicians, and the GSL minimize; economists, statisticians, and
Apophenia maximize. If you are doing a minimization, be sure that your function returns minus the objective
function's value.

\li The overall setup is about estimating the parameters of a model with data. The user
provides a data set and an unparameterized model, and the system tries parameterized
models until one of them is found to be optimal. The data is fixed.
The optimization tries a series of parameterized models,
searching for the one that is most likely.
In a non-stats setting, you will probably have \c NULL data, and will be optimizing the
parameters internal to the model. 

\li Because the unit of analysis is a parameterized model,
not just parameters, you need to have an \ref apop_model
wrapping your objective function. Here is a typical sort of function that one would
maximize; it is Rosenbrock's banana function, \f$(1-x)^2+ s(y - x^2)^2\f$, where the
scaling factor \f$s\f$ is fixed ahead of time, say at 100.

\code
    typedef struct {
            double scaling;
    } coeff_struct;

    double banana (double *params, coeff_struct *in){
            return (gsl_pow_2(1-params[0]) +
                    in->scaling*gsl_pow_2(params[1]-gsl_pow_2(params[0])));
    }
\endcode

The function returns a single number to be minimized.  You will need to write an
\ref apop_model to send to the optimizer, which is a two step process: write a log
likelihood function wrapping the real objective function, and a model that uses that
log likelihood. Here they are:

\code
double ll (apop_data *d, apop_model *in){
    return - banana(in->parameters->vector->data, in->more);
}

int main(){
    coeff_struct co = {.scaling=100};
    apop_model b = {"Bananas!", .log_likelihood= ll, .vsize=2,
                                .more = &co, .more_size=sizeof(coeff_struct)};
\endcode

The <tt>vsize=2</tt> specified that your parameters are a vector of size two.
That is, the list of <tt>double</tt>s to send to \c banana is set in
<tt>in->parameters->vector->data</tt>.
The \c more element of the \ref apop_model structure is designed to hold any arbitrary
structure of size \c more_size, which is useful for models that require additional
constants or other settings. See \ref settingswriting for more on handling model
settings.

\li Statisticians want the covariance and basic tests about the parameters. If you only
want the optimal value, then adding this line will shut off all auxiliary calculations:
\code
Apop_settings_add_group(your_model, apop_parts_wanted);
\endcode
See the documentation for \ref apop_parts_wanted_settings for details about how this
works.  It can also offer quite the speedup: especially for high-dimensional problems,
finding the covariance matrix without any information can take dozens of evaluations
of the objective function for each evaluation that is part of the search itself.
  \li MLEs have an especially large number of parameter tweaks that could be made;
see the \ref apop_mle_settings page.
  \li As a useful diagnostic, you can add a \c NULL \ref apop_data set to the MLE
settings in the <tt>.path</tt> slot, and it will be allocated and filled with the
sequence of points tried by the optimizer.
  \li Putting it all together, here is a full program to minimize Rosenbrock's banana
function. There are some extras: it uses two methods (notice how easy it is to re-run an
estimation with an alternate method, but the syntax for modifying a setting differs from
the initialization syntax) and checks that the results are accurate.

\include ../eg/banana.c

\code
Apop_settings_add_group(your_model, apop_mle, .want_path='y');
apop_model *out = apop_estimate(your_data, your_model);
apop_data_show(Apop_settings_get(out, apop_mle, path));
\endcode

\section constr Setting Constraints

The problem is that the parameters of a function must not take on certain values, either because the function is undefined for those values or because parameters with certain values would not fit the real-world problem.

The solution is to rewrite the function being maximized such that the function is continuous at the constraint boundary but takes a steep downward slope. The unconstrained maximization routines will be able to search a continuous function but will never return a solution that falls beyond the parameter limits.

If you give it an unconstrained likelihood function plus a separate constraint function, 
\ref apop_maximum_likelihood will combine them to a function that fits the above description and search accordingly.

A constraint function must do three things:
\li It must check the constraint, and if the constraint does not bind (i.e. the parameter values are OK), then it must return zero.
\li If the constraint does bind, it must return a penalty, that indicates how far off the parameter is from meeting the constraint.
\li If the constraint does bind, it must set a return vector that the likelihood function can take as a valid input. The penalty at this returned value must be zero.

The idea is that if the constraint returns zero, the log likelihood function will
return the log likelihood as usual, and if not, it will return the log likelihood at
the constraint's return vector minus the penalty. To give a concrete example, here
is a constraint function that will ensure that both parameters of a two-dimensional
input are both greater than zero, and that their sum is greater than two. As with the
constraints for many of the models that ship with Apophenia, it is a wrapper for \ref
apop_linear_constraint.

\code
static long double greater_than_zero_constraint(apop_data *data, apop_model *v){
    static apop_data *constraint = NULL;
    if (!constraint) constraint= apop_data_falloc((3,3,2), 0,  1, 0,
                                                           0,  0, 1,
                                                           2,  1, 1);
    return apop_linear_constraint(v->parameters->vector, constraint, 1e-3);
}
\endcode

\li\ref apop_linear_constraint()

\section simanneal Notes on simulated annealing

Simulated annealing is a controlled random walk.  As with the other methods, the
system tries a new point, and if it is better, switches. Initially, the system is
allowed to make large jumps, and then with each iteration, the jumps get smaller,
eventually converging. Also, there is some decreasing probability that if the new
point is less likely, it will still be chosen. Simulated annealing is best for
situations where there may be multiple local optima. Early in the random walk, the
system can readily jump from one to another; later it will fine-tune its way toward the
optimum. The number of points tested is determined by the parameters of the simulated
colling program, not the values returned by the likelihood function.  If you know your
function is globally convex (as are most standard probability functions), then this
method is overkill.


\section mlfns Useful functions

\li\ref apop_estimate_restart : Restarting an MLE with different settings can improve results.
\li\ref apop_maximum_likelihood : Rarely used. If a model has no \c estimate element, call \ref apop_estimate to run an MLE.
\li\ref apop_model_numerical_covariance
\li\ref apop_numerical_gradient
*/


/** \page moreasst Assorted

Some functions for missing data:

\li\ref apop_data_listwise_delete
\li\ref apop_ml_impute

A few more descriptive methods:

\li\ref apop_matrix_pca : Principal component analysis
\li\ref apop_anova : One-way or two-way ANOVA tables
\li\ref apop_rake : Iterative proportional fitting on large, sparse tables


General utilities:

\li\ref Apop_stopif : Apophenia's error-handling and warning-printing macro. 
\li\ref apop_opts : the global options
\li\ref apop_system : a printf-style wrapper around the standard \c system function.

Math utilities:

\li\ref apop_matrix_is_positive_semidefinite
\li\ref apop_matrix_to_positive_semidefinite
\li\ref apop_generalized_harmonic
\li\ref apop_multivariate_gamma
\li\ref apop_multivariate_lngamma
\li\ref apop_rng_alloc

*/

/** \page mingw MinGW

Minimalist GNU for Windows is indeed minimalist: it is not a full POSIX subsystem, and provides no package manager. Therefore, you will have to make some adjustments and install the dependencies yourself.

Matt P. Dziubinski successfully used Apophenia via MinGW; here are his instructions (with edits by BK):

\li get libregex (the ZIP file) from:
http://sourceforge.net/project/showfiles.php?group_id=204414&package_id=306189
\li get libintl (three ZIP files) from:
 http://gnuwin32.sourceforge.net/packages/libintl.htm .
 download "Binaries", "Dependencies", "Developer files"
\li follow "libintl" steps from:
http://kayalang.org/download/compiling/windows

\li Modify \c Makefile, adding -lpthread to AM_CFLAGS (removing -pthread) and -lregex to AM_CFLAGS and LIBS

\li Now compile the main library:
\code
make
\endcode

\li Finally, put one more expected directory in place and install:
\code
mkdir -p -- "/usr/local/Lib/site-packages"
make install
\endcode

\li You will get the usual warning about library paths, and may have to take the specified action:
\code
----------------------------------------------------------------------
Libraries have been installed in:
  /usr/local/lib

If you ever happen to want to link against installed libraries
in a given directory, LIBDIR, you must either use libtool, and
specify the full pathname of the library, or use the `-LLIBDIR'
flag during linking and do at least one of the following:
  - add LIBDIR to the `PATH' environment variable
    during execution
  - add LIBDIR to the `LD_RUN_PATH' environment variable
    during linking
  - use the `-LLIBDIR' linker flag

See any operating system documentation about shared libraries for
more information, such as the ld(1) and ld.so(8) manual pages.
----------------------------------------------------------------------
\endcode
*/


/* optionaldetails Implementation of optional arguments  [this section ignored by doxygen]
Optional and named arguments are among the most commonly commented-on features of Apophenia, so this page goes into full detail about the implementation. 

To use these features, see the all-you-really-need summary at the \ref designated
page. For a background and rationale, see the blog entry at http://modelingwithdata.org/arch/00000022.htm . 

I'll assume you've read both links before continuing.

OK, now that you've read the how-to-use and the discussion of how optional and named arguments can be constructed in C, this page will show how they are done in Apophenia. The level of details should be sufficient to implement them in your own code if you so desire.

There are three components to the process of generating optional arguments as implemented here:
\li Produce a \c struct whose elements match the arguments to the function.
\li Write a wrapper function that takes in the struct, unpacks it, and calls the original function.
\li Write a macro that makes the user think the wrapper function is the real thing.

None of these steps are really rocket science, but there is a huge amount of redundancy. 
Apophenia includes some macros that reduce the boilerplate redundancy significantly. There are two layers: the C-standard code, and the script that produces the C-standard code.

We'll begin with the C-standard header file:
\code 
#ifdef APOP_NO_VARIADIC
 void apop_vector_increment(gsl_vector * v, int i, double amt);
#else
 void apop_vector_increment_base(gsl_vector * v, int i, double amt);
 apop_varad_declare(void, apop_vector_increment, gsl_vector * v; int i; double amt);
#define apop_vector_increment(...) apop_varad_link(apop_vector_increment, __VA_ARGS__)
#endif
\endcode

First, there is an if/else that allows the system to degrade gracefully
if you are sending C code to a parser like swig, whose goals differ
too much from straight C compilation for this to work. Set \c
APOP_NO_VARIADIC to produce a plain function with no variadic support.

Else, we begin the above steps. The \c apop_varad_declare line expands to the following:

\code
typedef struct { 
    gsl_vector * v; int i; double amt ; 
} variadic_type_apop_vector_increment; 

void variadic_apop_vector_increment(variadic_type_apop_vector_increment varad_in);
  \endcode

So there's the ad-hoc struct and the declaration for the wrapper
function. Notice how the arguments to the macro had semicolons, like a
struct declaration, rather than commas, because the macro does indeed
wrap the arguments into a struct.

  Here is what the \c apop_varad_link would expand to:
  \code
#define apop_vector_increment(...) variadic_apop_increment_base((variadic_type_apop_vector_increment) {__VA_ARGS__})
  \endcode
That gives us part three: a macro that lets the user think that they are
making a typical function call with a set of arguments, but wraps what
they type into a struct.

Now for the code file where the function is declared. Again, there is is an \c APOP_NO_VARIADIC wrapper. Inside the interesting part, we find the wrapper function to unpack the struct that comes in.

\code
\#ifdef APOP_NO_VARIADIC 
 void apop_vector_increment(gsl_vector * v, int i, double amt){
\#else
apop_varad_head( void , apop_vector_increment){
    gsl_vector * apop_varad_var(v, NULL);
    Apop_stopif(!v, return, 0, "You sent me a NULL vector.");
    int apop_varad_var(i, 0);
    double apop_varad_var(amt, 1);
    apop_vector_increment_base(v, i, amt);
}

 void apop_vector_increment_base(gsl_vector * v, int i, double amt){
#endif
	v->data[i * v->stride]	+= amt;
}
\endcode

The 
\c apop_varad_head macro reduces redundancy, and will expand to
\code
void variadic_apop_vector_increment (variadic_type_variadic_apop_vector_increment varad_in)
\endcode

The function with this header thus takes in a single struct, and for every variable, there is a line like
\code
    double apop_varad_var(amt, 1);
\endcode
which simply expands to:
\code
    double amt = varad_in.amt ? varad_in.amt : 1;
\endcode
Thus, the macro declares each not-in-struct variable, and so there will need to be
one such declaration line for each argument. Apart from requiring declarations, you
can be creative: include sanity checks, post-vary the variables of the inputs, unpack
without the macro, and so on. That is, this parent function does all of the bookkeeping,
checking, and introductory shunting, so the base function can do the math. Finally,
the introductory section will call the base function.

The setup goes out of its way to leave the \c _base function in the public namespace,
so that those who would prefer speed to bounds-checking can simply call that function
directly, using standard notation. You could eliminate this feature by merging
the two functions.


<b>The m4 script</b>

The above is all you need to make this work: the varad.h file, and the above structures. But there is still a lot of redundancy, which can't be eliminated by the plain C preprocessor.

Thus, in Apophenia's code base (the one you'll get from checking out the git repository, not the gzipped distribution that has already been post-processed) you will find a pre-preprocessing script that converts a few markers to the above form. Here is the code that will expand to the above C-standard code:

\code
//header file
APOP_VAR_DECLARE void apop_vector_increment(gsl_vector * v, int i, double amt);

//code file
APOP_VAR_HEAD void apop_vector_increment(gsl_vector * v, int i, double amt){
    gsl_vector * apop_varad_var(v, NULL);
    Apop_stopif(!v, return, 0, "You sent me a NULL vector.");
    int apop_varad_var(i, 0);
    double apop_varad_var(amt, 1);
APOP_VAR_END_HEAD
	v->data[i * v->stride]	+= amt;
}
\endcode

It is obviously much shorter. The declaration line is actually a C-standard declaration with the \c APOP_VAR_DECLARE preface, so you don't have to remember when to use semicolons. The function itself looks like a single function, but there is again a marker before the declaration line, and the introductory material is separated from the main matter by the \c APOP_VAR_END_HEAD line. Done right, drawing a line between the introductory checks or initializations and the main function can improve readability.

The m4 script inserts a <tt>return function_base(...)</tt> at the end of the header
function, so you don't have to. If you want to call the function before the last line, you
can do so explicitly, as in the expansion above, and add a bare <tt>return;</tt> to
guarantee that the call to the base function that the m4 script will insert won't ever be
reached.

One final detail: it is valid to have types with commas in them---function arguments. Because commas get turned to semicolons, and m4 isn't a real parser, there is an exception built in: you will have to replace commas with exclamation marks in the header file (only). E.g.,

\code
APOP_VAR_DECLARE apop_data * f_of_f(apop_data *in, void *param, int n, double (*fn_d)(double ! void * !int));
\endcode

m4 is POSIX standard, so even if you can't read the script, you have the program needed to run it. For example, if you name it \c prep_variadics.m4, then run
\code
m4 prep_variadics.m4 myfile.m4.c > myfile.c
\endcode
*/


/**
\page gentle A quick overview

This is a "gentle introduction" to the Apophenia library. It is intended 
to give you some initial bearings on the typical workflow and the concepts and tricks that
the manual pages assume you have met.

This introduction assumes you already have some familiarity with C, how to compile a program, and how to use a debugger. If you want to install Apophenia now so you can try the samples on this page, see the \ref setup page.

An outline of this overview:

\li Apophenia fills a space between traditional C libraries and stats packages. 
\li The \ref apop_data structure represents a data set (of course). Data sets are inherently complex,
but there are many functions that act on \ref apop_data sets to make life easier.
\li The \ref apop_model encapsulates the sort of actions one would take with a model, like estimating model parameters or predicting values based on new inputs.
\li Databases are great, and a perfect fit for the sort of paradigm here. Apophenia
provides functions to make it easy to jump between database tables and \ref apop_data sets.

<em> The opening example</em>

Setting aside the more advanced applications and model-building tasks, let us begin with
the workflow of a typical fitting-a-model project using Apophenia's tools:

\li Read the raw data into the database using \ref apop_text_to_db.
\li Use SQL queries handled by \ref apop_query to massage the data as needed.
\li Use \ref apop_query_to_data to pull some of the data into an in-memory \ref apop_data set.
\li Call a model estimation such as \code apop_estimate (data_set, apop_ols)\endcode  or \code apop_estimate (data_set, apop_probit)\endcode to fit parameters to the data. This will return an \ref apop_model with parameter estimates.
\li Interrogate the returned estimate, by dumping it to the screen with \ref apop_model_print, sending its parameters and variance-covariance matrices to additional tests (the \c estimate step runs a few for you), or send the model's output to be input to another model.

Here is a concrete example of most of the above steps, which you can compile and run, to be discussed in detail below.

The program:

\include ols.c

To run this, you will need a file named <tt>data</tt> in comma-separated form, so
here is a set sufficient for the demonstration. The first column is the dependent variable; the
remaining columns are the independent:

\code
Y, X_1, X_2, X_3
2,3,4,5
1,2,9,3
4,7,9,0
2,4,8,16
1,4,2,9
9,8,7,6
\endcode

If you saved the code to <tt>sample.c</tt> and don't have a \ref makefile or other
build system, then you can compile it with

\code
gcc sample.c -std=gnu99 -lapophenia -lgsl -lgslcblas -lsqlite3 -o run_me
\endcode

or 

\code
clang sample.c -lapophenia -lgsl -lgslcblas -lsqlite3 -o run_me
\endcode

and then run it with <tt>./run_me</tt>. This compile line will work on any system with all the requisite tools,
but for full-time work with this or any other C library, you will probably want to write a \ref makefile .

The results are unremarkable---this is just an ordinary regression on too few data
points---but it does give us some lines of code to dissect. 

The first two lines in \c main() make use of a database.  
I'll discuss the value of the database step more at the end of this page, but for
now, note that there are several functions, \ref apop_query, and \ref
apop_query_to_data being the ones you will most frequently be using, that will allow you to talk to and pull data from either an SQLite or mySQL/mariaDB database. 

<em> Designated initializers</em>

Like this line,

\code
apop_text_to_db(.text_file="data", .tabname="d");
\endcode

many Apophenia functions accept named, optional arguments.  To give another example,
the \ref apop_data set has the usual row and column numbers, but also row and column
names. So you should be able to refer to a cell by any combination of name or number;
for the data set you read in above, which has column names, all of the following work:

\code
x = apop_data_get(data, 2, 3); //x == 0
x = apop_data_get(data, .row=2, .colname="X_3"); // x == 0
apop_data_set(data, 2, 3, 18);
apop_data_set(data, .colname="X_3", .row=2, .val= 18);
\endcode


Default values mean that the \ref apop_data_get, \ref apop_data_set, and \ref apop_data_ptr functions handle matrices, vectors, and scalars sensibly:
\code
//Let v be a hundred-element vector:
apop_data *v = apop_data_alloc(100);
double x1 = apop_data_get(v, 1);
apop_data_set(v, 2, .val=x1);

//A 100x1 matrix behaves like a vector
apop_data *m = apop_data_alloc(100, 1);
double m1 = apop_data_get(v, 1);

//let s be a scalar stored in a 1x1 apop_data set:
apop_data *v = apop_data_alloc(1);
double *scalar = apop_data_ptr(s);
\endcode

This form may be new to users of less user-friendly C libraries, but it it fully
conforms to the C standard (ISO/IEC 9899:2011). See the \ref designated page for details.


\section apop_data

A lot of real-world data processing is about quotidian annoyances about text versus
numeric data or dealing with missing values, and the \ref apop_data set and its
many support functions are intended to make data processing in C easy. Some users of
Apophenia use the library only for its \ref apop_data set and associated functions. See
\ref dataoverview for extensive notes on using the structure.

The structure includes seven parts:

\li a vector
\li a matrix
\li a grid of text elements
\li a vector of weights
\li names for everything: row names, a vector name, matrix column names, text names.
\li a link to a second page of data
\li an error marker

This is not a generic and abstract ideal, but is the sort of mess that real-world data sets look like. For
example, here is some data for a weighted OLS regression. It includes an outcome
variable in the vector, dependent variables in the matrix and text grid,
replicate weights, and column names in bold labeling the variables:

\htmlinclude apop_data_fig.html
\latexinclude apop_data_fig.tex

Apophenia will generally assume that one row across all of these elements
describes a single observation or data point.

See above for some examples of getting and setting individual elements.

Also, \ref apop_data_get, \ref apop_data_set, and \ref apop_data_ptr consider the vector to be the -1st column,
so using the data set in the figure, <tt>apop_data_get(sample_set, .row=0, .col=-1) == 1</tt>.

<em> Reading in data</em>

As per the example above, use \ref apop_text_to_data or \ref apop_text_to_db and then \ref apop_query_to_data.

<em> Subsets</em>

There are many macros to get subsets of the data. Each generates what is
considered to be a disposable view: once the variable goes out of scope (by the usual C
rules of scoping), it is no longer valid. However, these structures are all wrappers for pointers
to the base data, so all operations on the data view affect the base data. 

\include simple_subsets.c

All of these slicing routines are macros, because they generate several
background variables in the current scope (something a function can't do). Traditional
custom is to put macro names in all caps, like \c APOP_DATA_ROWS, which to modern
sensibilities looks like yelling. The custom has a logic: there are ways to hang
yourself with macros, so it is worth distinguishing them typographically. 
The documentation always uses a single capital.

<em> Basic manipulations</em>

\ref dataoverview lists a number of other manipulations of data sets, such as 
\ref apop_data_listwise_delete for quick-and-dirty removal of observations with <tt>NaN</tt>s,
\ref apop_data_split / \ref apop_data_stack,
or \ref apop_data_sort to sort all elements by a single column.

<em> Apply and map</em>

If you have an operation of the form <em>for each element of my data set, call this
function</em>, then you can use \ref apop_map to do it. You could basically do everything you
can do with an apply/map function via a \c for loop, but the apply/map approach is clearer
and more fun. Also, if you set OpenMP's <tt>omp_set_num_threads(N)</tt> for any \c N
greater than 1 (the default on most systems is the number of CPU cores), then the work
of mapping will be split across multiple CPU threads.  See \ref mapply for a number
of examples.

<em> Text</em>

Text in C is annoying. C already treats strings as pointer-to-characters, so a grid
of text data is a pointer-to-pointer-to-pointer-to-character. The text grid in the
\ref apop_data structure actually takes this form, but functions are provided to do
most or all the pointer work for you.  The \ref apop_text_alloc function
is really a realloc function: you can use it to resize the text grid as necessary. The
\ref apop_text_set function will do the pointer work in copying a single string to the
grid. Functions that act on entire data sets, like \ref apop_data_rm_rows, handle the
text part as well.

You have <tt>your_data->textsize[0]</tt> rows and <tt>your_data->textsize[1]</tt> columns. If you are using only the functions to this point, then empty elements are a blank string (<tt>""</tt>), not \c NULL.
For reading individual elements, refer to the \f$(i,j)\f$th text element via <tt>your_data->text[i][j]</tt>.

<em> Errors</em>

Many functions will set the <tt>error</tt> element of the \ref apop_data structure being operated on if anything goes wrong. You can use this to halt the program or take corrective action:

\code 
apop_data *the_data = apop_query_to_data("select * from d");
if (the_data->error) exit(1);
\endcode 

<em> The whole structure</em>

Here is a diagram of all of Apophenia's structures and how they
relate. It is taken from this
<a href="http://modelingwithdata.org/pdfs/cheatsheet.pdf">cheat sheet</a> (2 page PDF),
which will be useful to you if only because it lists some of the functions that act on
GSL vectors and matrices that are useful (in fact, essential) but out of the scope of the Apophenia documentation.

\image html http://apophenia.info/structs.png
\image latex ../structs.png width=18cm

All of the elements of the \ref apop_data structure are laid out at middle-left. You have
already met the vector, matrix, and weights, which are all a \c gsl_vector or \c gsl_matrix, and the text grid.

The diagram shows the \ref apop_name structure, which has received little mention so far because names
basically take care of themselves. A query will bring in column names (and row names if you set <tt>apop_opts.db_name_column</tt>), or use \ref apop_data_add_names to add names to your data
set and \ref apop_name_stack to copy from one data set to another.

The \ref apop_data structure has a \c more element, for when your data is best expressed
in more than one page of data. Use \ref apop_data_add_page, \ref apop_data_rm_page,
and \ref apop_data_get_page. Output routines will sometimes append an extra page of
auxiliary information to a data set, such as pages named <tt>\<Covariance\></tt> or
<tt>\<Factors\></tt>. The angle-brackets indicate a page that describes the data set
but is not a part of it (so an MLE search would ignore that page, for example).


Now let us move up the structure diagram to the \ref apop_model structure. 

\section gentle_model apop_model

Even restricting ourselves to the most basic operations, there are a lot of things that 
we want to do with our models: use a data set to estimate the parameters of a model (like the mean and
variance of a Normal distribution), or draw random numbers, or show the
expected value, or show the expected value of one part of the data given fixed values
for the rest of it. The \ref apop_model is intended to encapsulate most of these desires
into one object, so that models can easily be swapped around, modified to create new models, 
compared, and so on.

From the figure above, you can see that the \ref apop_model structure 
includes a number of informational items, key being the \c parameters, \c data, and
\c info elements; a list of settings to be discussed below; and a set of procedures for
many operations.  Its contents are not (entirely) arbitrary: the theoretical basis for
what is and is not included in an \ref apop_model, as well as its overall intent, are
described in this <a href="http://www.census.gov/srd/papers/pdf/rrs2014-06.pdf">U.S.
Census Bureau research report</a>.

There are helper functions that will allow you to avoid dealing with the model
internals. For example, the \ref apop_estimate helper function means you never have
to look at the model's \c estimate method (if it even has one), and you will simply
pass the model to a function, as with the above form:

\code
    apop_model *est = apop_estimate(data, apop_ols);
\endcode

\li Apophenia ships with a broad set of models, like \ref apop_ols, \ref apop_dirichlet,
    \ref apop_loess, and \ref apop_pmf (probability mass function); see the full list on <a href="http://apophenia.info/group__models.html">the models documentation page</a>. You would estimate the
parameters of any of them using \ref apop_estimate call, with the appropriate model as the second input.
\li The models that ship with Apophenia, like \ref apop_ols, include the procedures and some metadata, but are of course not yet estimated using a data set (i.e., <tt>data == NULL</tt>, <tt>parameters == NULL</tt>). The line above generated a new
model, \c est, which is identical to the base OLS model but has estimated parameters
(and covariances, and basic hypothesis tests, a log likelihood, \f$AIC_c\f$, \f$BIC\f$, et cetera), and a \c data pointer to the \ref apop_data set used for estimation. 
\li You will mostly use the models by passing them as inputs to 
functions like \ref apop_estimate, \ref apop_draw, or \ref apop_predict; more examples below.
Other than \ref apop_estimate, most require a parameterized model like \c est. After all, it doesn't make sense to
draw from a Normal distribution until its mean and standard deviation are specified.
\li If you know what the parameters should be, use \ref apop_model_set_parameters. E.g.

\code
apop_model *std_normal = apop_model_set_parameters(apop_normal, 0, 1);
apop_data *a_thousand_normals = apop_model_draws(std_normal, 1000);

apop_model *poisson = apop_model_set_parameters(apop_poisson, 1.5);
apop_data *a_thousand_waits = apop_model_draws(poisson, 1000);
\endcode

\li You can use \ref apop_model_print to print the various elements to screen.
\li You can combine and transform models with functions such as \ref apop_model_fix_params, \ref apop_model_coordinate_transform, or \ref apop_model_mixture. Each of these functions produce a new model, which can be estimated, re-combined, or otherwise used like any other model.
\li Writing your own models won't be covered in this introduction, but it can be easy to
copy and modify the procedures of an existing model to fit your needs. When in doubt, delete a procedure, because any procedures that are missing will have
defaults filled when used by functions like \ref apop_estimate (which uses \ref
apop_maximum_likelihood) or \ref apop_cdf (which uses integration via random draws). See \ref modeldetails for details.
\li There's a simple rule of thumb for remembering the order of the arguments to most of
Apophenia's functions, including \ref apop_estimate : the data always comes first.

<em> Settings</em>

How many bins are in a histogram? At what tolerance does the maximum likelihood
search end? What are the models being combined in an \ref apop_mixture distribution?

Apophenia organizes settings in <em>settings groups</em>, which are then attached
to models.  In the following snippet demonstrating Bayesian updating, we specify a Beta distribution prior.  If the
likelihood function were a Binomial distribution, \ref apop_update knows the closed-form
posterior for a Beta-Binomial pair, but in this case, with a PMF as a likelihood,
it will have to run Markov chain Monte Carlo. The \ref apop_mcmc_settings group attached
to the prior specifies details of how the run should work.

For a likelihood, we generate an empirical distribution---a PMF---from an input 
data set, via <tt>apop_estimate(your_data, apop_pmf)</tt>.
When we call \ref apop_update on the last line, it already has all of the above info
on hand.

\code
apop_model *beta = apop_model_set_parameters(apop_beta, 0.5, 0.25);
Apop_settings_add_group(beta, apop_mcmc, .burnin = 0.2, .periods =1e5);
apop_model *my_pmf = apop_estimate(your_data, apop_pmf);
apop_model *posterior = apop_update(.prior= beta, .likelihood = my_pmf);
\endcode

You will encounter model settings often when doing nontrivial work with models. All
can be set using a form like above.

<em> Databases and models</em>

Returning to the introductory example, you saw that (1) the
library expects you to keep your data in a database, pulling out the
data as needed, and (2) that the workflow is built around
\ref apop_model structures.

Starting with (2), 
if a stats package has something called a <em>model</em>, then it is
probably of the form Y = [an additive function of <b>X</b>], such as \f$y = x_1 +
\log(x_2) + x_3^2\f$. Trying new models means trying different
functional forms for the right-hand side, such as including \f$x_1\f$ in
some cases and excluding it in others. Conversely, Apophenia is designed 
to facilitate trying new models in the broader sense of switching out a 
linear model for a hierarchical, or a Bayesian model for a simulation. 
A formula syntax makes little sense over such a broad range of models.

As a result, the right-hand side is not part of 
the \ref apop_model. Instead, the data is assumed to be correctly formatted, scaled, or logged
before being passed to the model. This is where part (1), the database,
comes in, because it provides a proxy for the sort of formula specification language above:
 \code
apop_data *testme= apop_query_to_data("select y, x1, log(x2), pow(x3, 2) from data");
apop_model *est = apop_estimate(testme, apop_ols);
\endcode

Generating factors and dummies is also considered data prep, not model
internals. See \ref apop_data_to_dummies and \ref apop_data_to_factors.

Now that you have \c est, an estimated model, you can interrogate it. This is where Apophenia and its encapsulated
model objects shine, because you can do more than just admire the parameter estimates on
the screen: you can take your estimated data set and fill in or generate new data, use it
as an input to the parent distribution of a hierarchical model, et cetera. Some simple
examples:

 \code
 //If you have a new data set with missing elements (represented by NaN), you can fill in predicted values:
apop_predict(new_data_set, est);
apop_data_show(new_data_set);

 //Fill a matrix with random draws.
apop_data *d = apop_model_draws(est, .count=1000);

 //How does the AIC_c for this model compare to that of est2?
printf("ΔAIC_c=%g\n", apop_data_get(est->info, .rowname="AIC_c") 
                       - apop_data_get(est2->info, .rowname="AIC_c"));
\endcode

\section gentle_testing Testing

Here is the model for all hypothesis testing within Apophenia:

\li Calculate a statistic.
\li Describe the distribution of that statistic.
\li Work out how much of the distribution is (above|below|closer to zero than) the statistic.

There are a handful of named tests that produce a known statistic and then compare to a
known distribution, like \ref apop_test_kolmogorov or \ref apop_test_fisher_exact. For
traditional distributions (Normal, \f$t\f$, \f$\chi^2\f$), use the \ref apop_test convenience
function.

But if your model is not from the textbook, then you have the tools to apply the
above three-step process directly. First I'll give an overview of the three steps,
then another working example.

\li Model parameters are a statistic, and you know that  <tt>apop_estimate(your_data,
        your_model)</tt> will output a model with a <tt>parameters</tt> element.
\li The distribution of a parameter is also a model, so 
\ref apop_parameter_model will also return an \ref apop_model.
\li \ref apop_cdf takes in a model and a data point, and returns the area under the data
point.

Defaults for the parameter models are filled in via bootstrapping or resampling, meaning
that if your model's parameters are decidedly off the Normal path, you can still test
claims about the parameters.

The introductory example ran a standard OLS regression, whose output includes some
standard hypothesis tests; to conclude, let us go the long way and replicate those results
via the general \ref apop_parameter_model mechanism. The results here will of course be
identical, but the more general mechanism can be used in situations where the standard
models don't apply.

The first part of this program is identical to the program above. The second
half executes the three steps uses many of the above features: one of the inputs to
\ref apop_parameter_model (which row of the parameter set to use) is sent by adding a
settings group, we pull that row into a separate data set using \ref Apop_r, and we
set its vector value by referring to it as the -1st element.

\include ols2.c

Note that the procedure did not assume the model parameters had a certain form. It
queried the model for the distribution of parameter \c x_1, and if the model didn't have
a closed-form answer then a distribution via bootstrap would be provided. Then that model
was queried for its CDF. [The procedure does assume a symmetric distribution. Fixing this
is left as an exercise for the reader.] For a model like OLS, this is entirely overkill, 
which is why OLS provides the basic hypothesis tests automatically. But for models
where the distribution of parameters is unknown or has no closed-form solution, this
may be the only recourse.


This introduction has shown you the \ref apop_data set and some of the functions
associated, which might be useful even if you aren't formally doing statistical work but do have to deal with data with real-world elements like column names and mixed
numeric/text values. You've seen how Apophenia encapsulates many of a model's
characteristics into a single \ref apop_model object, which you can send with
data to functions like \ref apop_estimate, \ref apop_predict, or \ref apop_draw. Once
you've got your data in the right form, you can use this to simply estimate model
parameters, or as an input to later analysis.
*/

/** \page modeldetails Writing new models

The \ref apop_model is intended to provide a consistent expression of <em>any</em>
model that (implicitly or explicitly) expresses a likelihood of data given parameters,
including traditional linear models, textbook distributions, Bayesian hierarchies,
microsimulations, and any combination of the above.  The unifying feature is that
all of the models act over some data space and some parameter space (in some cases
one or both is the empty set), and can assign a likelihood for a fixed pair of
parameters and data given the model. This is a very broad requirement, often used
in the statistical literature.  For discussion of the theoretical structures, see <a
href="http://www.census.gov/srd/papers/pdf/rrs2014-06.pdf"><em>A Useful Algebraic System
of Statistical Models</em></a> (PDF).

This page is about writing new models from scratch, beginning with basic models and on
up to models with arbitrary internal settings, specific methods of Bayesian updating
using your model as a prior or likelihood, and so on. I assume you have already read
\ref modelsec on using models and have tried a few things with the
canned models that come with Apophenia, so you already know how a user handles basic
estimation, adding a settings group, and so on.

This page includes:

\li \ref write_likelihoods, giving a quick overview of how to write a new model from scratch.
\li \ref settingswriting, covering the writing of <em>ad hoc</em> structures to hold model- or method-specific details, like the number of periods for burning in an MCMC run or the number of bins in a histogram.
\li \ref vtables, covering the means of writing special-case routines for functions that are not part of the \ref apop_model itself, including the score or conjugate prior/likelihood pairs for \ref apop_update.
\li \ref modeldataparts, a detailed list of the requirements for the non-function elements of an \ref apop_model.
\li \ref methodsection, a detailed list of requirements for the function elements of an \ref apop_model.

\section write_likelihoods A walkthrough

Users are encouraged to always use models via the helper functions, like
\ref apop_estimate or \ref apop_cdf.  The helper functions do some boilerplate error
checking, and are where the defaults are called. For example, if your model has a \c log_likelihood
method but no \c p method, then \ref apop_p will use exp(\c log_likelihood). If you don't
give an \c estimate method, then \c apop_estimate will call \ref apop_maximum_likelihood.

So the game in writing a new model is to write just enough internal methods to give the helper functions what they need.
In the not-uncommon best case, all you need to do is write a log likelihood function or an RNG.

Here is how one would set up a model that could be estimated using maximum likelihood:

\li Write a likelihood function. Its header will look like this:

\code
long double new_log_likelihood(apop_data *data, apop_model *m);
\endcode

where \c data is the input data, and \c
m is the parametrized model (i.e. your model with a \c parameters element already filled in by the caller). 
This function will return the value of the log likelihood function at the given parameters.

\li Is this a constrained optimization? See \ref constr on how to set them. Otherwise, no constraints will be assumed.
\li Write the object:

\code
apop_model *your_new_model = &(apop_model){"The Me distribution", 
            .vsize=n0, .msize1=n1, .msize2=n2, .dsize=nd,
            .log_likelihood = new_log_likelihood };
\endcode

\li The first element is the <tt>.name</tt>, a human-language name for your model.
\li the \c vsize, \c msize1, and \c msize2 elements specify the shape of the parameter
set. For example, if there are three numbers in the vector, then set <tt>.vsize=3</tt>
and omit the matrix sizes. The default model prep routine will call

<tt>new_est->parameters = apop_data_alloc(vsize, msize1, msize2)</tt>. 
\li The \c dsize element is the size of one random draw from your model.
\li It's common to have [the number of columns in your data set] parameters; this
count will be filled in if you specify \c -1 for \c vsize, <tt>msize(1|2)</tt>, or
<tt>dsize</tt>. If the allocation is exceptional in a different way, then you will
need to allocate parameters by writing a custom \c prep method for the model.
\li If there are constraints, add a <tt>.constraint</tt> element for those too.

You already have more than enough that something like this will work (the \c dsize is used for random draws):
\code
apop_model *estimated = apop_estimate(your_data, your_new_model);
\endcode

Once that baseline works, you can fill in other elements of the \ref apop_model as needed.

For example, if you are using a maximum likelihood method to estimate parameters, you can get much faster estimates and better covariance estimates by specifying the dlog likelihood function (aka the score):

\code
void apop_new_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
    //do algebra here to find df/dp0, df/dp1, df/dp2....
    gsl_vector_set(gradient, 0, d_0);
    gsl_vector_set(gradient, 1, d_1);
}
\endcode
The score is not part of the model object, but is registered (see below) using 
\code
apop_score_insert(apop_new_dlog_likelihood, your_new_model);
\endcode

\subsection On Threading
Many procedures in Apophenia use OpenMP to thread operations. If a method can not be threaded, be sure to wrap it in an OpenMP critical region. E.g.,


\code

void apop_new_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
    #pragma omp critical (newdlog)
    {
        //un-threadable algebra here
    }
    gsl_vector_set(gradient, 0, d_0);
    gsl_vector_set(gradient, 1, d_1);
}
\endcode

\section settingswriting  Writing new settings groups

Your model may need additional settings or auxiliary information to function, which
would require associating a model-specific struct with the model.
The methods associated with such a struct usually begins with calls like
\code
long double ysg_ll(apop_data *d, apop_model *m){
    ysg_settings *sets = apop_settings_get(m, ysg);

    ...
}
\endcode

These model-specific structs are handled as expected by \ref apop_model_copy and \ref
apop_model_free, and many functions that modify or transform \ref apop_models try to
handle settings groups as expected. This section describes how to build a settings
group so all these automatic steps happen as expected, and your methods can reliably retrieve settings as needed.

But before getting into the detail of how to make model-specific groups of settings
work, note that there's a lightweight method of storing sundry settings, so in many
cases you can bypass all of the following.

The \ref apop_model structure has a \c void pointer named \c more which you can use to
point to a model-specific struct. If \c more_size is larger than zero (i.e. you set
it to <tt>your_model.more_size=sizeof(your_struct)</tt>), then it will be copied via \c
memcpy by \ref apop_model_copy, and freed by \ref apop_model_free. Apophenia's
estimation routines will never impinge on this item, so do what you wish with it.

The remainder of this subsection describes the information you'll have to provide to make
use of the conveniences described to this point: initialization of defaults, smarter
copying and freeing, and adding to an arbitrarily long list of settings groups attached
to a model.  You will need four items: a typedef for the structure itself, plus init, copy, and
free functions.  This is the sort of boilerplate that will be familiar to users of
object oriented languages in the style of C++ or Java, but it's really a list of
arbitrarily-typed elements, which makes this feel more like LISP. [And being a
reimplementation of an existing feature of LISP, this section will be macro-heavy.]

\li The settings struct will likely go into a header file, so 
here is a sample header for a new settings group named \c ysg_settings, with a dataset, its two sizes, and an owner-of-data marker. <tt>ysg</tt> stands for Your Settings Group; replace that substring with your preferred name in every instance to follow.

\code
typedef struct {
    int size1, size2;
    char *refs;
    apop_data *dataset;
} ysg_settings;

Apop_settings_declarations(ysg)
\endcode

The first item is a familiar structure definition. The last line is a macro that declares the
three functions below. This is everything you would
need in a header file, should you need one. These are just declarations; we'll write
the actual init/copy/free functions below.

The structure itself gets the full name, \c ysg_settings. Everything else is a macro
keyed on \c ysg, without the \c _settings part. Because of these macros, your struct
name must end in \c _settings.

If you have an especially simple structure, then you can generate the three functions with these three macros in your <tt>.c</tt> file:

\code
Apop_settings_init(ysg, )
Apop_settings_copy(ysg, )
Apop_settings_free(ysg, )
\endcode

These macros generate appropriate functions to do what you'd expect: allocating the
main structure, copying one struct to another, freeing the main structure.  
The spaces after the commas indicate that in these cases no special code gets added to
the functions that these macros expand into.

You'll never call the generated functions directly; they are called by \ref Apop_settings_add_group,
\ref apop_model_free, and other model or settings-group handling functions.

Now that initializing/copying/freeing of
the structure itself is handled, the remainder of this section will be about how to
add instructions for the structure internals, like data that is pointed to by the structure elements.

\li For the allocate function, use the above form if everything in your code defaults to zero/\c NULL.  
In most cases, though, you will need a new line declaring a default for every element in your structure. There is a macro to help with this too. 
These macros will define for your use a structure named \c in, and an output pointer-to-struct named \c out.
Continuing the above example:

\code
Apop_settings_init (ysg, 
      Apop_stopif(!in.size1, return NULL, 0, "I need you to give me a value for size1.");
      Apop_varad_set(size2, 10);
      Apop_varad_set(dataset, apop_data_alloc(out->size1, out->size2));
      Apop_varad_set(refs, malloc(sizeof(int)));
      *refs=1;
)
\endcode

Now, <tt>Apop_settings_add(a_model, ysg, .size1=100)</tt> would set up a group with a 100-by-10 data set, and set the owner bit to one. 

\li Some functions do extensive internal copying, so you will need a copy function even
if your code has no explicit calls to \ref apop_model_copy. The default above simply
copies every element in the structure. Pointers are copied, giving you two pointers
pointing to the same data. We have to be careful to prevent double-freeing later.

\code
//The elements of the set to copy are all copied by the function's boilerplate,
//and then make one additional modification:
Apop_settings_copy (ysg,
        (*refs)++;
)
\endcode

\li The struct itself is freed by boilerplate code, but add code in the free function
to free data pointed to by pointers in the main structure. The macro defines a
pointer-to-struct named \c in for your use. Continuing the example:

\code
Apop_settings_free (ysg,
    if (!(--in->refs)) {
        free(in->dataset);
        free(in->refs);
    }
)
\endcode

With those three macros in place and the header as above, Apophenia will treat your
settings group like any other, and users can use \ref Apop_settings_add_group to
populate it and attach it to any model.

\section vtables Registering new methods in vtables

The settings groups are for adding arbitrary model-specific nouns; vtables are for
adding arbitrary model-specific verbs.

For any given function (e.g., entropy, the dlog likelihood, Bayesian updating), there is
probably a special case for well-known models like the Normal distribution. 
Any function may maintain a registry of models and associated special-case procedures, aka a vtable.

This subsection will discuss how to add
a function to an existing vtable.

\li See \ref apop_update, \ref apop_score, \ref apop_predict, \ref apop_model_print, and \ref
apop_parameter_model for examples and procedure-specific details.
\li Write a function following the given type definition.
\li Use the associated <tt>_vtable_add</tt> function to add the function and associate it
with the given model. For example, to add a Beta-binomial routine named \c betabinom
to the registry of Bayesian updating routines, use <tt>apop_update_vtable_add(betabinom,
apop_beta, apop_binomial)</tt>.
\li Lookups happen based on a hash that takes into account the elements of the model
that will be used in the calculation. For example, the \c apop_update_hash takes in two
models and calculates the hash based on the address of the prior's \c draw method and
the likelihood's \c log_likelihood or \c p method. Thus, a vtable lookup for new models
that re-use the same methods (at the same addresses in memory) will still find the
same special-case function.
\li If you need to deregister the function, use the associated deregister function,
e.g. <tt>apop_update_vtable_drop(apop_beta, apop_binomial)</tt>. You can guarantee that a method will not be re-added by following up the <tt>_drop</tt> with, e.g., <tt>apop_update_vtable_add(NULL, apop_beta, apop_binomial)</tt>.
\li Place a call to <tt>..._vtable_add</tt> in the \c prep method of the given model, thus ensuring that the auxiliary functions are registered after the first time the model is sent to \ref apop_estimate.

The easiest way to set up a new vtable is to copy/paste/modify an existing one. Briefly:

\li See the existing setups in the vtables portion of <tt>apop.h</tt>. 
\li Cut/paste one and do a search and replace to change the name to match your desired use.
\li Set the typedef to describe the functions that get added to the vtable.
\li Rewrite the hash function to check the part of the inputs that interest you. For
example, the update vtable associates functions with the \c draw, \c log_likelihood,
and \p methods of the model. A model where these elements are identical but the name
is changed will still match.

\section modeldataparts The data elements

The remainder of this section covers the detailed expectations regarding the elements
of the \ref apop_model structure. I begin with the data (non-function) elements,
and then cover the method (function) elements. Some of the following will be
requirements for all models and some will be advice to authors; I use the accepted
definitions of <a href="http://tools.ietf.org/html/rfc2119">"must", "shall", "may"</a>
and related words.

\subsection datasubsec Data

\li Each row should be a single observation. 
For example, \ref apop_bootstrap_cov depends on each row being an iid observation to function correctly.
Calculating the Bayesian Information Criterion (BIC) requires knowing
the number of observations in the data, and assumes that row count==observation count.
For complex data, the \ref apop_data_pack and \ref  apop_data_unpack functions can help with this. 

\li Some functions (bootstrap again, or many uses of \ref apop_kl_divergence) use \ref
apop_draw to use your model's RNG (or a default) to draw a
  value, write it to the matrix element of the data set, and then move on to an
  estimation or other step. In this case, the data sent in will be entirely in the \c
  ->matrix element of the \ref apop_data set sent to model methods. Your \c likelihood, \c p, \c cdf, and \c estimate routines
  must accept data as a single row of the matrix of the \ref apop_data set for such functions to work.
  They may accept other formats. Tip: you can use \ref apop_data_pack and \ref apop_data_unpack to convert a structured set to a single row and back again.

\li Your routines may accept other data formats, as per contract with the user.
    For example, regression-type functions use a function named \c ols_shuffle
    to convert a matrix where the first column is the dependent variable to a data
    set with dependent variable in the vector and a column of ones in the first
    matrix column; see \ref dataprep.

\subsection paramsubsec Parameters, vsize, msize1,  msize2

\li The sizes will be used by the \c prep method of the model; see below. Given the model \c m and its elements \c m.vsize, \c m.msize1, \c m.msize2,
    functions that need to allocate a parameter set will do so via <tt>apop_data_alloc(m.vsize, m.msize1, m.msize2)</tt>. 


\subsection infosubsec Info

\li The first page, named \c &lt;info&gt; is typically a list of scalars. Nothing is guaranteed, but the elements may include:

\li AIC: <a href="https://en.wikipedia.org/wiki/Akaike's_Information_Criterion">Aikake Information Criterion</a>
\li AIC_c: AIC with a finite sample correction. ``<em>Generally, we advocate the use of AIC_c when the ratio \f$n/K\f$ is small (say \f$< 40\f$)</em>'' [Kenneth P. Burnham, David R. Anderson: <em>Model Selection and Multi-Model Inference</em>, p 66, emphasis in original.]
\li BIC: <a href="https://en.wikipedia.org/wiki/Bayesian_information_criterion">Bayesian Information Criterion</a>
\li R squared
\li R squared adj
\li log likelihood
\li status.

For those elements that require a count of input data, the calculations assume each row in the input \ref apop_data set is a single datum.

Get these via, e.g., <tt>apop_data_get(your_model->info, .rowname="log likelihood")</tt>.
When writing for any arbitrary function, be prepared to handle \c NaN, indicating that the element is not calculated or saved in the info page by the given model.

For OLS-type estimations, each row corresponds to the row in the original data. For
filling in of missing data, the elements may appear anywhere, so the row/col indices are
essential.

\subsection settingsgroupmention settings, more

In object-oriented jargon, settings groups are the private elements of the data set,
to be pulled out in certain contexts, and ignored in all others. Therefore, there are
no rules about internal use. The \c more element of the \ref apop_model provides a lightweight
means of attaching an arbitrary struct to a model. See \ref settingswriting above for details.

\li As many settings groups of different types as desired can be added to a single \ref apop_model.
\li One \ref apop_model can not hold two settings groups of the same type. Re-additions cause the removal of the previous version of the group.


\section methodsection Methods

\subsection psubsection p, log_likelihood

\li Function headers look like  <tt>long double your_p_or_ll(apop_data *d, apop_model *params)</tt>.
\li The inputs are an \ref apop_data set and an \ref apop_model, which should include the elements needed to fully estimate the probability/likelihood (probably a filled <tt>->parameters</tt> element, possibly a settings group added by the user).
\li We assume that the parameters have been set, by users via \ref apop_estimate or \ref apop_model_set_parameters, or by \ref apop_maximum_likelihood by its search algorithms. If the parameters are necessary, the function shall check that the parameters are not \c NULL and set the model's \c error element to \c 'p' if they are missing.
\li Return \c NaN on errors. If an error in the input model is found, the function may set the input model's \c error element to an appropriate \c char value.
\li If your model includes both \c log_likelihood and \c p methods, it must be the case that <tt>log(p(d, m))</tt> equals <tt>log_likelihood(d, m)</tt> for all \c d and \c m.
\li If observations are assumed to be iid, you can probably use \ref apop_map_sum to write the core of the log likelihood function.

\subsection prepsubsection prep

\li Function header looks like <tt>void your_prep(apop_data *data, apop_model *params)</tt>.
\li Re-prepping a model after it has already been prepped shall have no effect. Where there is ambiguity with the other requirements, this takes precedence.
\li The model's <tt>data</tt> pointer shall be set to point to the input data.
\li The \c info element shall be allocated and its title set to "<Info>".
\li If \c vsize, \c msize1, or \c msize2 are -1, then the prep function shall set them to the width of the input data.
\li If \c dsize is -1, then the prep function shall set it to the width of the input data.
\li If the \c parameters element is not allocated, the function shall allocate it via <tt>apop_data_alloc(vsize, msize1, msize2)</tt> (or equivalent).
\li The default is \ref apop_model_clear. It does all of the above.
\li The input data may be modified by the prep routine. For example, the \ref apop_ols prep routine shuffles a single input matrix as described above under \c data, and the \ref apop_pmf prep routine calls \ref apop_data_pmf_compress on the input data.
\li The prep routine may initialize any desired settings groups. Unless otherwise
stated, these should not be removed if they are already there, so that users can override defaults by adding a settings group before starting an estimation.
\li If any functions associated with the model need to be added to 
a vtable (see above), the registration shall happen here. Registration may also happen elsewhere.

\subsection estimatesubsection estimate

\li Function header looks like  <tt> void your_estimate(apop_data *data, apop_model *params)</tt>.
\li Assume that the prep routine has already been run. Notably, this means that parameters have been allocated.
\li Assume that the \c parameters hold garbage (as in a \c malloc without a subsequent assignment to the <tt>malloc</tt>-ed space).
\li The function modifies the input model, and returns nothing. Note that this is different from the wrapper function, \ref apop_estimate, which makes a copy of its input model, preps it, and then calls the \c estimate function with the prepeped copy.
\li The function shall set the \c parameters of the input model. For consistency with other models, the estimate should be the maximum likelihood estimate, unless otherwise documented.
\li Additional settings may be set.
\li The model's \c &lt;Info&gt; page may be filled with statistics. For scalars like log likelihood and AIC, use \ref apop_data_add_named_elmt.
\li Data should not be modified by the \c estimate routine; any changes to the data made by \c estimate must be documented.
\li The default called by \ref apop_estimate is \ref apop_maximum_likelihood.
\li If errors occur during processing, set the model's \c error element to a single character. Documentation should include the list of error characters and their meaning.

\subsection drawsubsection draw

\li Function header looks like <tt>void your_draw(double *out, gsl_rng* r, apop_model *params)</tt>
\li Assume that model \c paramters are set, via \ref apop_estimate or \ref apop_model_set_parameters. The author of the draw method should check that \c parameters are not \c NULL if needed and fill the output with NaNs if necessary parameters are not set.
\li Caller inputs a pointer-to-<tt>double</tt> of length \c dsize; user is expected to make sure that there is adequate space. Caller also inputs a \c gsl_rng, already allocated (probably via \ref apop_rng_alloc, possibly from \ref apop_rng_get_thread).
\li The function shall fill the space pointed to by the input pointer with a random draw from the data space, where the likelihood of any given observation is proportional to its likelihood as given by the \c p method. Data shall be reduced to a single vector via \ref apop_data_pack if it is not already a single vector.

\subsection cdfsubsection cdf

\li Function header looks like <tt>long double your_cdf(apop_data *d, apop_model *params)</tt>.
\li Assume that \c parameters are set, via \ref apop_estimate or \ref apop_model_set_parameters. The author of the CDF method should check that \c parameters are not \c NULL and return NaN if necessary parameters are not set.
\li The CDF method must accept data as a single row of data in the \c matrix of the input \ref apop_data set (as per a draw produced using the \c draw method). May accept other formats.
\li Returns the percentage of the likelihood function \f$\leq\f$ the first row of the input data. The definition of \f$\leq\f$ is chosen by the model author.
\li If one is not already present, an \c apop_cdf_settings group may be added to the model to store temp data. See the \ref apop_cdf function for details.

\subsection constraintsubsection constraint

\li Function header looks like <tt>long double your_constraint(apop_data *data, apop_model *params)</tt>.
\li Assume that \c parameters are set, via \ref apop_estimate, \ref apop_model_set_parameters, or the internals of an MLE search. The author of the constraint method should check that \c parameters are not \c NULL and return NaN if necessary parameters are not set.
\li See \ref apop_linear_constraint for a useful basis and/or example. Many constraints can be written as wrappers for this function.
\li If the constraint is met, then return zero.
\li If the constraint fails, then (1) move the \c parameters in the input model to a
constraint-satisfying value, and (2) return the distance between the input parameters and
what you've moved the parameters to. The choice of within-bounds parameters and distance function is left to the author of the constraint function.
*/
