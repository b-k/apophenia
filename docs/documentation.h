/* Apophenia's narrative documentation
Copyright (c) 2005--2013 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

/**  \mainpage  Apophenia--the intro

Apophenia is an open statistical library for working with data sets and statistical
models. It provides functions on the same level as those of the typical stats package
(such as OLS, probit, or singular value decomposition) but gives the user more
flexibility to be creative in model-building.  The core functions are written in C,
but experience has shown them to be easy to bind to in Python/Julia/Perl/Ruby/&c.

It is written to scale well, to comfortably work with gigabyte data sets, million-step simulations, or
computationally-intensive agent-based models. If you have tried using other open source
tools for computationally demanding work and found that those tools weren't up to the
task, then Apophenia is the library for you.

<h5>The goods</h5> 

The library has been growing and improving since 2005, and has been downloaded over 10,000 times. To date, it has over two hundred functions to facilitate statistical computing, such as:

\li OLS and family, discrete choice models like probit and logit, kernel density estimators, and other common models
\li database querying and maintenance utilities
\li moments, percentiles, and other basic stats utilities
\li t-tests, F-tests, et cetera
\li Several optimization methods available for your own new models
\li It does <em>not</em> re-implement basic matrix operations or build yet another database
engine. Instead, it builds upon the excellent <a href="http://sources.redhat.com/gsl/">GNU
Scientific</a> and <a href="http://www.sqlite.org/">SQLite</a> libraries. MySQL/mariaDB is also supported.

For the full list, click the <a href="globals.html">index</a> link from the header.

<h5><a href="https://github.com/b-k/Apophenia/archive/pkg.zip">Download Apophenia here</a>.</h5>

Most users will just want to download the latest packaged version linked from the <a
href="https://github.com/b-k/Apophenia/archive/pkg.zip">Download
Apophenia here</a> header.

Those who would like to work on a cutting-edge copy of the source code
can get the latest version by cutting and pasting the following onto
the command line. If you follow this route, be sure to read the development README in the
<tt>Apophenia</tt> directory this command will create.

\code
git clone https://github.com/b-k/Apophenia.git
\endcode

<!--git clone git://apophenia.git.sourceforge.net/gitroot/apophenia/apophenia
cvs -z3 -d:ext:<i>(your sourceforge login)</i>@cvs.sourceforge.net:/cvsroot/apophenia co -P apophenia
cvs -z3 -d:pserver:anonymous@cvs.sf.net:/cvsroot/apophenia checkout -P apophenia
svn co https://apophenia.svn.sourceforge.net/svnroot/apophenia/trunk/apophenia --> 

<h5>The documentation</h5>

To start off, have a look at this <a href="gentle.html">Gentle Introduction</a> to the
library.

<a href="outline.html">The outline</a> gives a more detailed narrative.

The <a href="globals.html">index</a> lists every function in the
library, with detailed reference
information. Notice that the header to every page has a link to the outline and the index.

To really go in depth, download or pick up a copy of
<a href="http://modelingwithdata.org">Modeling with Data</a>,
which discusses general methods for doing statistics in C with the GSL 
and SQLite, as well as Apophenia itself. <a
href="http://ben.klemens.org/pdfs/klemens-modelcats.pdf"><em>A Cross-paradigm Modeling
Framework</em></a> (PDF) discusses some of the theoretical structures underlying
the library.

There is a <a href="https://github.com/b-k/Apophenia/wiki">wiki</a> with some convenience
functions, tips, and so on.

<h5>Notable features</h5> 
Much of what Apophenia does can be done in any typical statistics package. The \ref
apop_data element is much like an R data frame, for example, and there is nothing special
about being able to invert a matrix or take the product of two matrices with a single
function call (\ref apop_matrix_inverse and \ref apop_dot, respectively). 
Even more advanced features like Loess smoothing (\ref apop_loess) and the Fisher Exact
Test (\ref apop_test_fisher_exact) are not especially Apophenia-specific. But here are
some things that are noteworthy.

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
across cycles is less than <tt>eps</tt>), just add a settings group specifying the
tolerance at which the cycle should stop: <tt>Apop_settings_add_group(your_model,
apop_mle, .dim_cycle_tolerance=eps)</tt>.

\li The Iterative Proportional Fitting algorithm, \ref apop_rake, is best-in-breed,
designed to handle large, sparse matrices.

\li As well as the \ref apop_data structure, Apophenia is built around a model object,
the \ref apop_model. This allows for consistent treatment of distributions, regressions,
simulations, machine learning models, and who knows what other sorts of models you can
dream up. By transforming and combining existing models, it is easy to build complex
models from simple sub-models.

\li For example, the \ref apop_update function does Bayesian updating on any two
well-formed models. If they are on the table of conjugates, that is correctly handled,
and if they are not, an appropriate variant of MCMC 
produces an empirical distribution. The output is yet another model, from which you
can make random draws, or which you can use as a prior for another round of Bayesian
updating.

\li Of course, it's a C library, meaning that you can build applications using Apophenia
for the data-processing back-end of your program. For example, it is currently used
in production for certain aspects of processing for the U.S. Census Bureau's American
Community Survey.


<h5>Contribute!</h5> 

\li Develop a new model object.
\li Contribute your favorite statistical routine.
\li Package Apophenia into an RPM, apt, portage, cygwin package.
\li Report bugs or suggest features.
\li Write bindings for your preferred language. For example, here are early versions of <a
href="http://modelingwithdata.org/arch/00000173.htm"> a Julia
wrapper</a> and <a href="https://r-forge.r-project.org/projects/rapophenia/">an R
wrapper</a> which you could expand upon.

If you're interested,  <a href="mailto:fluffmail@f-m.fm">write to the maintainer</a> (Ben Klemens), or join the
<a href="https://github.com/b-k/Apophenia">GitHub</a> project.
*/

/** \page eg Some examples
 Here are a few pieces of sample code, gathered from elsewhere in the documentation, for testing your installation or to give you a sense of what code with Apophenia's tools looks like. If you'd like more context or explanation, please click through to the page from which the example was taken.

In the documentation for the \ref apop_ols model, a program to read in data and run a regression. You'll need to go to that page for the sample data and further discussion.

From the \ref setup page, an example of gathering data from two processes, saving the input to a database, then doing a later analysis:

\include draw_to_db.c

A demonstration of \ref apop_plot_line_and_scatter . You'll need a
database from the {\em Modeling with Data} sample code, at
http://modelingwithdata.org/appendices.html.

\include scatter.c

In the \ref outline section on map/apply, a new \f$t\f$-test on every row, with all operations acting on entire rows rather than individual data points:

\include t_test_by_rows.c

In the documentation for \ref apop_query_to_text, a program to list all the tables in an SQLite database.
\include ls_tables.c

A demonstration of fixing parameters to create a marginal distribution, via \ref apop_model_fix_params
\include fix_params.c

Several uses of the \ref apop_dot function
\include dot_products.c
 */



/** \page setup Setting up
\section cast The supporting cast 
To use Apophenia, you will need to have a working C compiler, the GSL (v1.7 or higher) and SQLite installed. 

\li Some readers are unfamiliar with modern package managers and common methods for setting up a C development environment; see 
<a href="http://modelingwithdata.org/appendix_o.html">Appendix O</a> of <em> Modeling with Data</em> for an introduction.

\li Other pages on this site have a few more notes for \ref windows "Windows" users or \ref mingw users.

\li Install the basics using your package manager. E.g., try

\code
sudo apt-get install make gcc libgsl0-dev libsqlite3-dev
\endcode

or 

\code
sudo yum install make gcc gsl-devel libsqlite3x-devel
\endcode

\li <a href="https://github.com/b-k/Apophenia/archive/pkg.zip">Download Apophenia here</a>. 

\li Once you have the library downloaded, compile it using 

\code
tar xvzf apop*tgz && cd apophenia-0.999
./configure && make && sudo make install && make check
\endcode

If you decide not to keep the library on your system, run <tt>sudo make uninstall</tt>
from the source directory to remove it.

\li A \ref makefile will help immensely when you want to compile your program.


\subsection sample_program Sample programs
Here is a sample program so you can test your setup.  There is another short,
complete program in the \ref apop_ols entry which runs a simple OLS regression on a
data file. Follow the instructions there to compile and run.

The sample program here is intended to show how one would integrate Apophenia into an existing program. For example, say that you are running a simulation of two different treatments, or say that two sensors are posting data at regular intervals. You need to gather the data in an organized form, and then ask questions of the resulting data set.  Below, a thousand draws are made from the two processes and put into a database. Then, the data is pulled out, some simple statistics are compiled, and the data is written to a text file for inspection outside of the program.  This program will compile cleanly with the sample \ref makefile.

\include draw_to_db.c

*/

/** \page windows The Windows page

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
If you are missing anything else, the program will probably tell you.
The following are not necessary but are good to have on hand as long as you are going to be using Unix and programming.
\li svn -- to partake in the versioning system
\li emacs -- steep learning curve, but people love it
\li ghostscript (for reading .ps/.pdf files)
\li openssh -- needed for cvs
\li perl, python, ruby -- these are other languages that you might also be interested in
\li tetex -- write up your documentation using the nicest-looking formatter around
\li X11 -- a windowing system
X-Window will give you a nicer environment in which to work.  After you start Cygwin, just type <tt>startx</tt> to bring up a more usable, nice-looking terminal (and the ability to do a few thousand other things which are beyond the scope of this documentation).  Once you have Cygwin installed and a good terminal running, you can follow along with the remainder of the discussion without modification.

Sqlite3 is difficult to build from scratch, but you can get a packaged version by pointing Cygwin's install program to the Cygwin Ports site: http://cygwinports.dotsrc.org/ .

Second, some older (but still pretty recent) versions of Cygwin have a search.h file which doesn't include the function lsearch().  If this is the case on your system, you will have to update your Cygwin installation.

Finally, windows compilers often spit out lines like:
\code
Info: resolving _gsl_rng_taus by linking to __imp__gsl_rng_taus (auto-import)
\endcode
These lines are indeed just information, and not errors. Feel free to ignore them.

[Thanks to Andrew Felton and Derrick Higgins for their Cygwin debugging efforts.]
*/

/** \page notroot  Not root? 
If you aren't root, then you will need to create a subdirectory in your home directory in which to install packages. The GSL and SQLite installations will go like this. The key is the <tt>--prefix</tt> addition to the <tt>./configure</tt> command.
\code
export MY_LIBS = src   #choose a directory name to be created in your home directory.
tar xvzf pkg.tgz       #change pkg.tgz to the appropriate name
cd package_dir         #same here.
mkdir $HOME/$MY_LIBS
./configure --prefix $HOME/$MY_LIBS
make
make install   #Now you don't have to be root.
echo "export LD_LIBRARY_PATH=$HOME/$MY_LIBS:\$LD_LIBRARY_PATH" >> ~/.bashrc
\endcode
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
CFLAGS = -g -Wall 
LDLIBS = -lapophenia -lgsl -lgslcblas -lsqlite3

$(PROGNAME): $(objects)
\endcode

By the way, if your system has \c pkg-config, then you can use it for a slightly more robust and readable makefile. Replace the above C and link flags with:
\code
CFLAGS = -g -Wall `pkg-config --cflags apophenia`
LDLIBS = `pkg-config --libs apophenia`
\endcode
The \c pkg-config program will then fill in the appropriate directories and libraries. Pkg-config knows Apophenia depends on the GSL and database libraries, so you need only list the most-dependent library.
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

\li There is some overhead to checking for defaults. If this is noticeable for your situation, you can pass on this convenient form and call the underlying function directly, by adding \c _base to the name and giving all arguments:

\code
apop_text_to_db_base("infile.txt", "intable", 0, 1, NULL);
\endcode

\li If one of the optional elements is an RNG and you do not provide one, I use one from \ref apop_rng_get_thread.
\li For exhaustive details on implementation of the above (should you wish to write new functions that behave like this) see the \ref optionaldetails page.

*/

/** \defgroup global_vars The global variables */
/** \defgroup mle Maximum likelihood estimation */
/** \defgroup command_line "Command line programs" */


/** \page admin The admin page

Just a few page links:
\li <a href=todo.html>The to-do list</a>
\li <a href=bug.html>The known bug list</a>
\li <a href="modules.html">The documentation page list</a>
*/

/** \page outline An outline of the library
  \anchor outline

ALLBUTTON

Outlineheader preliminaries Getting started

If you are entirely new to Apophenia, \ref gentle "have a look at the Gentle Introduction here".

As well as the information in this outline, there is a separate page covering the details of 
 \ref setup "setting up a computing environment" and another page with \ref eg "some sample code" for your perusal.

For another concrete example of folding Apophenia into a project, have a look at this \ref sample_program "sample program".



Outlineheader c Some notes on C and Apophenia's use of C utilities.
 
Outlineheader learning  Learning C

<a href="http://modelingwithdata.org">Modeling with Data</a> has a full tutorial for C, oriented at users of standard stats packages. More nuts-and-bolts tutorials are <a href="http://www.google.com/search?hl=es&amp;c2coff=1&amp;q=c+tutorial">in abundance</a>.  Some people find pointers to be especially difficult; fortunately, there's a <a href="http://www.youtube.com/watch?v=6pmWojisM_E">claymation cartoon</a> which clarifies everything.

Coding often relies on gathering together many libraries; there is a section at the bottom of this outline linking to references for some libraries upon which Apophenia builds.

endofdiv

Outlineheader usagenotes  Usage notes

Here are some notes about the technical details of using the Apophenia library in your development environment.

<b> Header aggregation </b>

If you put 
\code
#include <apop.h>
\endcode
at the top of your file, then it will call virtually every header file you could need: gsl_matrix.h, gsl_blas.h, sqlite3.h, stdio.h, string.h, math.h, apophenia_all_of_them.h, et cetera. Of course, if you get `implicit declaration of...' then you will need to manually include something else.

<b> Easier calling syntax</b>

Several functions (but nowhere near all) allow optional named arguments to functions. For
example:
\code
apop_vector_distance(v1, v2); //assumes Euclidean distance
apop_vector_distance(v1, .metric='M'); //assumes v2=0, uses Manhattan metric.
\endcode

See the \ref designated page for details of this syntax.

endofdiv

Outlineheader liblinking Libraries 

Your best bet is to write yourself a \ref makefile "Makefile".  If you don't want to use the sample \ref makefile "Makefile", then here are some notes for the command line.  When compiling, you will need to tell the compiler to use the Apophenia library. That is, you will need to call GCC with <tt>gcc -lapophenia</tt> (as well as the other usual flags). For example,
\code
gcc sample.c -lapophenia -lsqlite3 -lgsl -lgslcblas -o run_me -g -Wall -O3
\endcode
will produce an output file named <tt>run_me</tt> from the input source code <tt>sample.c</tt>. It will include symbols for debugging (<tt>-g</tt>), print all warnings (<tt>-Wall</tt>), try every optimization GCC has (<tt>-O3</tt>, where <tt>-O</tt> stands for optimization) and will correctly find the functions in the Apophenia, GSL, and SQLite libraries (<tt>-lapophenia -lgsl ...</tt>). 
Order matters in the linking list: the files a package depends on should be listed after the package. E.g., since sample.c depends on Apophenia, <tt>gcc sample.c -lapophenia</tt> will work, while <tt>gcc -lapophenia sample.c</tt> is likely to give you errors. Similarly, list <tt>-lapophenia</tt> before <tt>-lgsl</tt>, which comes before <tt>-lgslcblas</tt>.

endofdiv


endofdiv

Outlineheader debugging  Errors, logging, debugging and stopping

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


<h5>Verbosity level and logging</h5> 

The global variable <tt>apop_opts.verbose</tt> determines how many notifications and warnings get printed:

-1: turn off logging, print nothing (ill-advised) <br>
0: notify only of failures and clear danger <br>
1: warn of technically correct but odd situations that might indicate, e.g., numeric instability <br>
2: debugging-type information; print queries  <br>
3: give me everything, such as the state of the data at each iteration of a loop.

These levels are of course subjective, but should give you some idea of where to place the
verbosity level. The default is 1.

The messages are printed to the \c FILE handle at <tt>apop_opts.log_file</tt>. If
this is blank (which happens at startup), then this is set to \c stderr. This is the
typical behavior for a console program. Use

\code
apop_opts.log_file = fopen("mylog", "w");
\endcode

to write to the \c mylog file instead of \c stderr.

As well as the error and warning messages, some functions can also print diagnostics,
using the \ref Apop_notify macro.  For example, \ref apop_query and friends will print the
query sent to the database engine iff <tt>apop_opts.verbose >=2</tt> (which is useful
if you are substituting in many <tt>\%s</tt>es). The diagnostics attempt to follow
the same verbosity scale as the warning messages.

<h5>Stopping</h5> 

Warnings and errors never halt processing. It is up to the calling function to decide
whether to stop.

When running the program under a debugger, this is an annoyance: we want to stop as
soon as a problem turns up.

The global variable <tt>apop_opts.stop_on_warning</tt> changes when the system halts:

\c 'n': never halt. If you were using Apophenia to support a user-friendly GUI, for example, you would use this mode.<br>
The default: if the variable is not set, halt on severe errors, continue on all warnings.<br>
\c 'v': If the verbosity level of the warning is such that the warning would print to screen, then halt;
if the warning message would be filtered out by your verbosity level, continue.<br>
\c 'w': Halt on all errors or warnings.

See the documentation page for individual functions for details on how errors are reported to the caller.

The end of <a href="http://modelingwithdata.org/appendix_o.html">Appendix O</a> of <em>Modeling with Data</em> offers some GDB 
macros which can make dealing with Apophenia from the GDB command line much more pleasant.

endofdiv

Outlineheader About SQL, the syntax for querying databases

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

Finally, Apophenia provides a few nonstandard SQL functions to facilitate math via database; see \ref db_moments.
endofdiv

Outlineheader threads Threading

Apophenia uses OpenMP for threading. You generally do not need to know how OpenMP works
to use Apophenia, and many points of work will thread without your doing anything.

\li All functions strive to be thread-safe. Part of how this is achieved is that static
variables are marked as thread-local or atomic, as per the C standard. There still
exist compilers that can't implement thread-local or atomic variables, in which case
your safest bet is to set OMP's thread count to one as below (or get a new compiler).

\li Some functions modify their inputs. It is up to you to use those functions in
a thread-safe manner. The \ref apop_matrix_realloc operates correctly in a threaded
environment, but if you have two threads resizing the same \c gsl_matrix at the same
time, you're going to have problems.

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

\li \ref apop_map and friends distribute their \c for loop over the input \ref apop_data
set across multiple threads. Therefore, be careful to send thread-unsafe functions to
it only after calling \c omp_set_num_threads(1).

\li There are a few functions, like \ref apop_model_draws, that rely on \ref apop_map, and
therefore also thread by default.

\li The function \ref apop_rng_get_thread retrieves a statically-stored RNG specific
to a given thread. Therefore, if you use that function in the place of a \c gsl_rng,
you can parallelize functions that make random draws.

\li \ref apop_rng_get_thread allocates its store of threads using apop_opts.rng_seed,
then incrementing that seed by one. You thus probably have threads with seeds 479901,
479902, 479903, .... [If you have a better way to do it, please feel free to modify the
code to implement your improvement and submit a pull request on Github.]

See <a href="http://modelingwithdata.org/arch/00000175.htm">this tutorial on C
threading</a> if you would like to know more, or are unsure about whether your functions
are thread-safe or not.

endofdiv

Outlineheader mwd The book version

Apophenia co-evolved with <em>Modeling with Data: Tools and Techniques for Statistical Computing</em>. You can read about the book, or download a free PDF copy of the full text, at <a href="http://modelingwithdata.org">modelingwithdata.org</a>.

If you are at this site, there is probably something there for you, including a tutorial on C and general computing form, SQL for data-handing, several chapters of statistics from various perspectives, and more details on working Apophenia. 

As with many computer programs, the preferred manner of citing Apophenia is to cite its related book.
Here is a BibTeX-formatted entry, which should be be easy to re-shape for other environments:

\code 
@book{klemens:modeling,
    title = "Modeling with Data: Tools and Techniques for Statistical Computing",
    author="Ben Klemens",
    year=2008,
    publisher="Princeton University Press"
}
\endcode

endofdiv

Outlineheader status What is the status of the code?

[This section last updated 3 January 2014.]

Apophenia was first posted to SourceForge in February 2005, which means that we've had
several years to develop and test the code in real-world applications. 

The test suite, including the sample code and solution set for _Modeling with Data_,
is about 5,500 lines over 135 files. gprof reports that it covers about 90% of the
7,700 lines in Apophenia's code base. A broad rule of thumb for any code base is
that the well-worn parts, in this case functions like \ref apop_data_get and \ref
apop_normal's <tt>log_likelihood</tt>, are likely to be entirely reliable, while the
out-of-the-way functions (maybe the score for the Beta distribution) are worth a bit
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

The \ref apop_model structure is much more ambitious---Apophenia is really intended
to be a novel system for developing models---and its internals can still be improved.
The promise underlying the structure is that you can provide just one item, such as
an RNG or a likelihood function, and the structure will do all of the work to fill in
computationally-intensive methods for everything else; see \ref settingswriting for
the gory details. Some directions aren't quite there yet (such as RNG -> most other
things); the default likelihood -> RNG function only works for models with a domain of
unidimensional reals; the PMF model (the essential bridge from data sets to empirical
models) needs an internal index for faster lookups.  Readers are invited perhaps
contribute better methods (such as an alternate means of estimating mixture models),
or filling in more of existing models (write a dlog likelihood function for a model
that does not currently have one), or submit new standard models not yet included.


endofdiv

Outlineheader ext How do I write extensions?

It's not a package, so you don't need an API---write your code and <tt>include</tt> it like any
other C code. The system is written to not require a registration or initialization step
to add a new model or other such parts.  You can read the notes below on generating
new models, which have to conform to some rules if they are to play well with \ref
apop_estimate, \ref apop_draw, and so forth.

Once your new model or function is working, please post the code or a link to the code on the 
<a href="https://github.com/b-k/Apophenia/wiki">Apophenia wiki</a>.

endofdiv

endofdiv

Outlineheader dataoverview Data sets

The \ref apop_data structure represents a data set.  It joins together a \c gsl_vector, a \c gsl_matrix, an \ref apop_name, and a table of strings. It tries to be lightweight, so you can use it everywhere you would use a \c gsl_matrix or a \c gsl_vector.

If you are viewing the HTML documentation, here is a diagram showing a sample data set with all of the elements in place. Together, they represet a data set where each row is an observation, which includes both numeric and text values, and where each row/column may be named.

\htmlinclude apop_data_fig.html

For example, let us say that you are running a regression: there is a vector for the one dependent variable, and a matrix for the several independent variables. Think of them as a partitioned matrix, where the vector is column -1, and the first column of the matrix is column zero. Here is some code to print the entire matrix. Notice that the column counter \c i starts counting at -1.

\code
for (int j = 0; j< data->matrix->size1; j++){
    printf("%s\t", apop_name_get(data->names, j, 'r'));
    for (int i = -1; i< data->matrix->size2; i++)
        printf("%g\t", apop_data_get(data, j, i));
    printf("\n");
}
\endcode

Most functions assume that the data vector, data matrix, and text have the same row count: \c data->vector->size==data->matrix->size1  and \c data->vector->size==*data->textsize. This means that the \ref apop_name structure doesn't have separate vector_names, row_names, and text_row_names elements: the rownames are assumed to apply for all.

See below for notes on easily managing the \c text element and the row/column names.

The \ref apop_data set includes a \c more pointer, which will typically
be \c NULL, but may point to another \ref apop_data set. This is
intended for a main data set and a second or third page with auxiliary
information: estimated parameters on the front page and their
covariance matrix on page two, or predicted data on the front page and
a set of prediction intervals on page two. \c apop_data_copy and \c apop_data_free
will handle all the pages of information. The \c more pointer is not
intended as a linked list for millions of data points---you can probably
find a way to restructure your data to use a single table (perhaps via
\ref apop_data_pack and \ref apop_data_unpack).

Easy data manipulation is essential for enjoying life as a researcher.  Thus, there are a great many functions to collate, copy, merge, sort, prune, and otherwise manipulate the \ref apop_data structure and its components.

\li\ref apop_data_add_named_elmt()
\li\ref apop_data_copy()
\li\ref apop_data_fill()
\li\ref apop_data_memcpy()
\li\ref apop_data_pack()
\li\ref apop_data_rm_columns()
\li\ref apop_data_sort()
\li\ref apop_data_split()
\li\ref apop_data_stack()
\li\ref apop_data_transpose()
\li\ref apop_data_unpack()
\li\ref apop_matrix_copy()
\li\ref apop_matrix_realloc()
\li\ref apop_matrix_rm_columns()
\li\ref apop_matrix_stack()
\li\ref apop_text_add()
\li\ref apop_text_paste()
\li\ref apop_text_to_data()
\li\ref apop_vector_bounded()
\li\ref apop_vector_copy()
\li\ref apop_vector_fill()
\li\ref apop_vector_stack()
\li\ref apop_vector_realloc()

Note: Apophenia builds upon the GSL, but it would be inappropriate to redundantly replicate the GSL's documentation here. You will find a link to the full GSL documentation at the end of this outline, and the text of the outline will list the names a few common functions.  The GSL's naming scheme is very consistent, so a simple reminder of the function name may be sufficient for your needs.

\li <tt>gsl_matrix_swap_rows (gsl_matrix * m, size_t i, size_t j)</tt>
\li <tt>gsl_matrix_swap_columns (gsl_matrix * m, size_t i, size_t j)</tt>
\li <tt>gsl_matrix_swap_rowcol (gsl_matrix * m, size_t i, size_t j)</tt>
\li <tt>gsl_matrix_transpose_memcpy (gsl_matrix * dest, const gsl_matrix * src)</tt>
\li <tt>gsl_matrix_transpose (gsl_matrix * m) : square matrices only</tt>
\li <tt>gsl_matrix_set_all (gsl_matrix * m, double x)</tt>
\li <tt>gsl_matrix_set_zero (gsl_matrix * m)</tt>
\li <tt>gsl_matrix_set_identity (gsl_matrix * m)</tt>
\li <tt>void gsl_vector_set_all (gsl_vector * v, double x)</tt>
\li <tt>void gsl_vector_set_zero (gsl_vector * v)</tt>
\li <tt>int gsl_vector_set_basis (gsl_vector * v, size_t i)</tt>: set all elements to zero, but set item \f$i\f$ to one.
\li <tt>gsl_vector_reverse (gsl_vector * v)</tt>: reverse the order of your vector's elements
\li <tt>gsl_vector_ptr</tt> and <tt>gsl_matrix_ptr</tt>. To increment an element in a vector use, e.g., <tt>*gsl_vector_ptr(v, 7) += 3;</tt> or <tt>(*gsl_vector_ptr(v, 7))++</tt>.

Outlineheader readin Reading from text files

The \ref apop_text_to_data() function takes in the name of a text file with a grid of data in (comma|tab|pipe|whatever)-delimited format and reads it to a matrix. If there are names in the text file, they are copied in to the data set. See \ref text_format for the full range and details of what can be read in.

If you have any columns of text, then you will need to go through the database: use
\ref apop_text_to_db() to convert your text file to a database table, 
do any database-appropriate cleaning of the input data, then use \ref
apop_query_to_data() or \ref apop_query_to_mixed_data() to pull the data to a \ref apop_data set.

endofdiv

Outlineheader datalloc Alloc/free

\li\ref apop_data_alloc()
\li\ref apop_data_calloc()
\li\ref apop_data_free()
\li\ref apop_text_alloc() : allocate or resize the text part of an \ref apop_data set.
\li\ref apop_text_free()

    See also:

\li <tt>gsl_vector * gsl_vector_alloc (size_t n)</tt>
\li <tt>gsl_vector * gsl_vector_calloc (size_t n)</tt>
\li <tt>void gsl_vector_free (gsl_vector * v)</tt>
\li <tt>gsl_matrix_memcpy (gsl_matrix * dest, const gsl_matrix * src)</tt>

endofdiv

Outlineheader gslviews 	Using views

You may need to pull a row of a matrix as a distinct vector. The GSL makes it easy to do so, and Apophenia provides several convenience functions to save any remaining hassle. The data in the new vector will point to the base data, not copy it. Some examples:

\code
apop_data *d = apop_query_to_data("select obs1, obs2, obs3 from a_table");

//Get a column using its name
Apop_col_t(d, "obs1", ov);
double obs1_sum = apop_vector_sum(ov);

//Get a row using its index as a vector
Apop_matrix_row(d->matrix, 0, v);
double first_row_sum = apop_vector_sum(v);

//Get a row or rows as a standalone one-row apop_data set
Apop_row(d, 0, d1);
apop_data_print(d1);

//ten rows starting at row 3:
Apop_data_rows(d, 3, 10, d10);
apop_data_show(d10);

//First column's sum
Apop_col_v(d, 0, colone);
double first_col_sum = apop_vector_sum(colone);

//Pull a 10x5 submatrix, whose first element is the (2,3)rd
//element of the parent data set's matrix
Apop_submatrix(d, 2,3, 10,5, subm);
double first_col_sum = apop_matrix_sum(subm);
\endcode

\li\ref Apop_col
\li\ref Apop_row
\li\ref Apop_cols
\li\ref Apop_rows
\li\ref Apop_col_v
\li\ref Apop_row_v
\li\ref Apop_col_t
\li\ref Apop_row_t
\li\ref Apop_col_tv
\li\ref Apop_row_tv
\li\ref Apop_matrix_row
\li\ref Apop_matrix_col
\li\ref Apop_submatrix

These macros make use of a set of GSL matrices that produce output of type <tt>gsl_matrix_view</tt> and <tt>gsl_vector_view</tt>. These types point to the source data, and add metadata to turn the data into a coherent matrix/vector. Apophenia's macros generate these views, then pull the matrix/vector from them, so you never have to deal with the <tt>_view</tt> structs directly.


The view is an automatic variable, not a pointer, and therefore disappears at the end of the scope in which it is declared. If you want to retain the data after the function exits, copy it to another vector: 

\code
Apop_matrix_row(d->matrix, 2, rowtwo);
gsl_vector *outvector = apop_vector_copy(rowtwo);
\endcode

Curly braces always delimit scope, not just at the end of a function. 
These macros work by generating a number of local variables, which you may be able to see in your
debugger. When program evaluation exits a given block, all variables in that block are
erased. Here is some sample code that won't work:

\code
if (get_odd){
    Apop_row(data, 1, outdata);
} else {
    Apop_row(data, 0, outdata);
}
apop_data_show(outdata); //breaks: no outdata in scope.
\endcode


For this if/then statement, there are two sets of local variables
generated: one for the \c if block, and one for the \c else block. By the last line,
neither exists. You can get around the problem here by making sure to not put the macro
declaring new variables in a block. E.g.:

\code
Apop_row(data, get_odd ? 1 : 0, outdata);
apop_data_show(outdata);
\endcode


This is a general rule about how variables declared in blocks will behave, but because the
macros obscure the variable declarations, it is especially worth watching out for here.

endofdiv

Outlineheader setgetsec Set/get

The set/get functions can act on both names or indices. Sample usages:


\code
double twothree = apop_data_get(data, 2, 3); //just indices
apop_data_set(data, .rowname="A row", .colname="this column", .val=13);
double AIC = apop_data_get(data, .rowname="AIC", .col=-1, .page="Info");
\endcode

\li\ref apop_data_get()
\li\ref apop_data_set()
\li\ref apop_data_ptr() : returns a pointer to the element.

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

endofdiv

Outlineheader mapplysec   Map/apply

\anchor outline_mapply
These functions allow you to send each element of a vector or matrix to a function, either producing a new matrix (map) or transforming the original (apply).  The \c ..._sum functions return the sum of the mapped output.

There is an older and a newer set of functions. The older versions, which act on <tt>gsl_matrix</tt>es or <tt>gsl_vector</tt>s have more verbose names; the newer versions, which act on the elements of an \ref apop_data set, use the \ref designated syntax to ad a few options and a more brief syntax.

You can do many things quickly with these functions.

Get the sum of squares of a vector's elements:

\code
double sum_of_squares = apop_map_sum(dataset, gsl_pow_2); //given <tt> apop_data *dataset</tt>
double sum_of_squares = apop_vector_map_sum(v, gsl_pow_2); //given <tt> gsl_vector *v</tt>
\endcode

Here, we create an index vector [\f$0, 1, 2, ...\f$].

\code
double index(double in, int index){return index;}
apop_data *d = apop_data_alloc(100);
apop_map(d, .fn_di=index, .inplace='y');
\endcode

Given your log likelihood function and a data set where each row of the matrix is an observation, find the total log likelihood:

\code
static double your_log_likelihood_fn(const gsl_vector * in){[your math goes here]}

double total_log_likelihood = apop_matrix_map_sum(dataset, your_log_likelihood_fn);
\endcode

How many missing elements are there in your data matrix? 

\code
static double nan_check(const double in){ return isnan(in);}

int nan_count_data = apop_map_sum(in, nan_check, .part='m');
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

Notice how the older \ref apop_vector_apply uses file-global variables to pass information into the functions, while the \ref apop_map uses a pointer to the constant parameters to input to the functions.

\include t_test_by_rows.c

One more toy example, demonstrating the use of the \ref Apop_row :

\include apop_map_row.c

\li\ref apop_map()
\li\ref apop_map_sum()
\li\ref apop_matrix_apply()
\li\ref apop_matrix_map()
\li\ref apop_matrix_map_all_sum()
\li\ref apop_matrix_map_sum()
\li\ref apop_vector_apply()
\li\ref apop_vector_map()
\li\ref apop_vector_map_sum()

endofdiv

Outlineheader  matrixmathtwo  Basic Math

\li\ref apop_vector_exp : exponentiate every element of a vector
\li\ref apop_vector_log : take the log of every element of a vector
\li\ref apop_vector_log10 : take the log (base 10) of every element of a vector
\li\ref apop_vector_distance : find the Euclidean distance between two vectors
\li\ref apop_vector_normalize : scale/shift a matrix to have mean zero, sum to one, et cetera
\li\ref apop_matrix_normalize : apply apop_vector_normalize to every column or row of a matrix

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

endofdiv
            
Outlineheader  matrixmath  Matrix math

\li\ref apop_dot : matrix \f$\cdot\f$ matrix, matrix \f$\cdot\f$ vector, or vector \f$\cdot\f$ matrix
\li\ref apop_matrix_determinant : returns the determinant of the input matrix
\li\ref apop_matrix_inverse : returns the inverse of the input matrix
\li\ref apop_det_and_inv : find determinant and inverse at the same time

Again, see the GSL documentation for voluminous further options

endofdiv

Outlineheader  sumstats  Summary stats

\li\ref apop_data_summarize ()
\li\ref apop_vector_moving_average()
\li\ref apop_vector_percentiles()

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

endofdiv

Outlineheader  moments  Moments

For most of these, you can add a weights vector for weighted mean/var/cov/....

\li\ref apop_data_correlation ()
\li\ref apop_data_covariance ()
\li\ref apop_data_summarize ()
\li\ref apop_matrix_mean ()
\li\ref apop_matrix_mean_and_var ()
\li\ref apop_matrix_sum ()
\li\ref apop_mean()
\li\ref apop_sum()
\li\ref apop_var()
\li\ref apop_vector_correlation ()
\li\ref apop_vector_cov ()
\li\ref apop_vector_kurtosis ()
\li\ref apop_vector_kurtosis_pop ()
\li\ref apop_vector_mean()
\li\ref apop_vector_skew()
\li\ref apop_vector_skew_pop()
\li\ref apop_vector_sum()
\li\ref apop_vector_var()
\li\ref apop_vector_var_m ()

endofdiv

Outlineheader convsec   Conversion among types

An \ref apop_data set binds together names, text, weights, a vector, and a matrix. If you need an individual element, you can naturally just point to it:

\code
apop_data *d = apop_text_to_mixed_data("vmw", "select result, age, income, sampleweight from data");
double avg_result = apop_vector_mean(d->vector, .weights=d->weights);
\endcode

In the other direction, you can use compound literals to wrap an \ref apop_data struct around a loose vector or matrix:

\code
//Given:
gsl_vector *v;
gsl_matrix *m;

// Then this form wraps the elements into \ref apop_data structs. Note that
// these are not pointers: they're automatically allocated and therefore
// the extra memory use for the wrapper is cleaned up on exit from scope.
apop_data dv = (apop_data){.vector=v}; 
apop_data dm = (apop_data){.matrix=m};

apop_data *v_dot_m = apop_dot(dv, dm);

//all at once:
apop_data *v_dot_m2 = apop_dot((apop_data){.vector=v}, (apop_data){.matrix=m});
\endcode

\li\ref apop_array_to_vector() : <tt>double*</tt>\f$\to\f$ <tt>gsl_vector</tt>
\li\ref apop_data_fill() : <tt>double*</tt>\f$\to\f$  \ref apop_data. See also \ref apop_data_falloc
\li\ref apop_matrix_to_data()
\li\ref apop_text_to_data() : delimited text file\f$\to\f$ \ref apop_data
\li\ref apop_text_to_db() : delimited text file\f$\to\f$ database
\li\ref apop_vector_to_data()
\li\ref apop_vector_to_matrix()

endofdiv

Outlineheader names   Name handling

If you generate your data set from the database via \ref apop_query_to_data (or
\ref apop_query_to_text or \ref apop_query_to_mixed_data) then column names appear
as expected.  Set <tt>apop_opts.db_name_column</tt> to the name of a column in your
query result to use that column name for row names.

Sample uses, given \ref apop_data set <tt>d</tt>:

\code
int row_name_count = d->names->rowct
int col_name_count = d->names->colct
int text_name_count = d->names->textct

//Manually add a name:
apop_name_add(d->names, "the vector", 'v');
apop_name_add(d->names, "row 0", 'r');
apop_name_add(d->names, "row 1", 'r');
apop_name_add(d->names, "row 2", 'r');
apop_name_add(d->names, "numeric column 0", 'c');
apop_name_add(d->names, "text column 0", 't');
apop_name_add(d->names, "The name of the data set.", 'h');

//or append several names at once
apop_data_add_names(d, 'c', "numeric column 1", "numeric column 2", "numeric column 3");

//point to element i from:

char *rowname_i = d->names->row[i];
char *colname_i = d->names->col[i];
char *textname_i = d->names->text[i];

//The vector also has a name:
d->names->vector;
\endcode

\li\ref apop_name_add() : add one name
\li\ref apop_data_add_names() : add a sequence of names at once
\li\ref apop_name_stack() : copy the contents of one name list to another
\li\ref apop_name_find() : find the row/col number for a given name.
\li\ref apop_name_print() : print the \ref apop_name struct, for diagnostic purposes.

endofdiv

Outlineheader textsec   Text data

The \ref apop_data set includes a grid of strings, <tt>text</tt>, for holding text data. 

Text should be encoded in UTF-8. US ASCII is a subset of UTF-8, so that's OK too.

There are a few simple forms for handling the \c text element of an \c apop_data set, which handle the tedium of memory-handling for you.

\li Use \ref apop_text_alloc to allocate the block of text. It is actually a realloc function, which you can use to resize an existing block without leaks.
\li Use \ref apop_text_add to add text elements. It replaces any existing text in the given slot without memory leaks.
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
is to rely on \ref apop_text_add, \ref apop_text_alloc, and \ref apop_text_free to do
the memory management for you.

Here is a sample program that uses these forms, plus a few text-handling functions.

\include eg/text_demo.c

\li\ref apop_data_transpose() : also transposes the text data. Say that you use
<tt>dataset = apop_query_to_text("select onecolumn from data");</tt> then you have a
sequence of strings, <tt>d->text[0][0], d->text[1][0], </tt>.... After <tt>apop_data
*dt = apop_data_transpose(dataset)</tt>, you will have a single list of strings,
<tt>dt->text[0]</tt>, which is often useful as input to list-of-strings handling
functions.

\li\ref apop_query_to_text()
\li\ref apop_text_alloc() : allocate or resize the text part of an \ref apop_data set.
\li\ref apop_text_add()
\li\ref apop_text_paste() : convert a table of little strings into one long string.
\li\ref apop_text_unique_elements() : geta sorted list of unique elements for one column of text.
\li\ref apop_text_free() : you may never need this, because \ref apop_data_free calls it.

endofdiv

Outlineheader fact   Generating factors

\em Factor is jargon for a numbered category. Number-crunching programs prefer integers over text, so we need a function to produce a one-to-one mapping from text categories into numeric factors. 

A \em dummy is a variable that is either one or zero, depending on membership in a given group. Some methods (typically when the variable is an input or independent variable) prefer dummies; some methods (typically for outcome or dependent variables) prefer factors. The functions that generate factors and dummies will add an informational page to your \ref apop_data set with a name like <tt>\<categories for your_column\></tt> listing the conversion from the artificial numeric factor to the original data. Use \ref apop_data_get_factor_names to get a pointer to that page.

\li\ref apop_data_to_dummies()
\li\ref apop_data_to_factors()
\li\ref apop_data_get_factor_names()
\li\ref apop_text_unique_elements()
\li\ref apop_vector_unique_elements()

endofdiv

endofdiv

Outlineheader dbs Databases

These are convenience functions to handle interaction with SQLite or mySQL/mariaDB. They open one and only one database, and handle most of the interaction therewith for you.

You will probably first use \ref apop_text_to_db to pull data into the database, then \ref apop_query to clean the data in the database, and finally \ref apop_query_to_data to pull some subset of the data out for analysis.

Also see the \ref db_moments page for not-SQL-standard math functions that you can
use when sending queries from Apophenia, such as \c pow, \c stddev, or \c sqrt.

\li \ref apop_text_to_db : Read a text file on disk into the database. Most data analysis projects start with a call to this.
\li \ref apop_data_print : If you include the argument <tt>.output_type='d'</tt>, this prints your \ref apop_data set to the database.
\li \ref apop_query : Manipulate the database, return nothing (e.g., insert rows or create table).
\li \ref apop_db_open : Optional, for when you want to use a database on disk.
\li \ref apop_db_close : If you used \ref apop_db_open, you will need to use this too.
\li \ref apop_table_exists : Check to make sure you aren't reinventing or destroying data. Also, a clean way to drop a table.

\li Apophenia reserves the right to insert temp tables into the opened database. They will all have names beginning with "apop_", so the reader is advised to not use tables with such names, and is free to ignore or delete any such tables that turn up.

Outlineheader dbout Out

\li\ref apop_db_to_crosstab()
\li\ref apop_query_to_data()
\li\ref apop_query_to_float()
\li\ref apop_query_to_mixed_data()
\li\ref apop_query_to_text()
\li\ref apop_query_to_vector()

endofdiv

Outlineheader dbin In

        See the print functions below. By setting <tt>apop_opts.output_type = 'd'</tt>, \ref apop_data sets (or \c gsl_matrixes and \c gsl_vectors) are `printed' to the database.

endofdiv

Outlineheader dbmath Math in the db  
    
\li\ref apop_db_paired_t_test()
\li\ref apop_db_t_test()
\li\ref db_moments 

endofdiv

endofdiv

Outlineheader Modesec Models

This segment discusses the use of existing \ref apop_model objects.
If you need to write a new model, see \ref modeldetails.

Outlineheader introtomodels Introduction

Begin with the most common use:
the \c estimate function will estimate the parameters of your model. Just prep the data, select a model, and produce an estimate:

\code
    apop_data *data = apop_query_to_data("select outcome, in1, in2, in3 from dataset");
    apop_model *the_estimate = apop_estimate(data, apop_probit);
    apop_model_print(the_estimate, NULL);
\endcode

Along the way to estimating the parameters, most models also find covariance estimates for
the parameters, calculate statistics like log likelihood, and so on, which the final print statement will show.

The <tt>apop_probit</tt> model that ships with Apophenia is unparameterized:
<tt>apop_probit.parameters==NULL</tt>. The output from the estimation,
<tt>the_estimate</tt>, has the same form as <tt>apop_probit</tt>, but
<tt>the_estimate->parameters</tt> has a meaningful value.

Outlineheader covandstuff More estimation output

A call to \ref apop_estimate produces more than just the estimated parameters. Most will
produce any of a covariance matrix, some hypothesis tests, a list of expected values, log
likelihood, AIC, AIC_c, BIC, et cetera.

First, note that if you don't want all that, 
adding to your model an \ref apop_parts_wanted_settings group with its default values (see below on settings groups) signals to
the model that you want only the parameters and to not waste CPU time on covariances,
expected values, et cetera. See the \ref apop_parts_wanted_settings documentation for examples and
further refinements.

\li The actual parameter estimates are in an \ref apop_data set at \c your_model->parameters.

\li Scalar statistics of the model are listed in the output model's \c info group, and can
be retrieved via a form like

\code
apop_data_get(your_model->info, .rowname="log likelihood");
//or
apop_data_get(your_model->info, .rowname="AIC");
\endcode

\li Covariances of the parameters are a page appended to the parameters; retrieve via

\code
apop_data *cov = apop_data_get_page(your_model->parameters, "<Covariance>");
\endcode

\li The table of expected values (typically including expected value, actual value, and
residual) is a page stapled to the main info page. Retrieve via:

\code
apop_data *predict = apop_data_get_page(your_model->info, "<Predicted>");
\endcode

endofdiv

But we expect much more from a model than just estimating parameters from data.  

Continuing the above example where we got an estimated Probit model named \c the_estimate, we can interrogate the estimate in various familiar ways. In each of the following examples, the model object holds enough information that the generic function being called can do its work:

\code
apop_data *expected_value = apop_predict(NULL, the_estimate);

double density_under =  apop_cdf(expected_value, the_estimate);

gsl_rng *rng = apop_rng_alloc(234211);
apop_data *draws = apop_model_draws(the_estimate, .count=1000, .rng=rng);
\endcode

Apophenia ships with many well-known models for your immediate use, including
probability distributions, such as the \ref apop_normal, \ref apop_poisson, or \ref apop_beta models. The data is assumed to have been drawn from a given distribution and the question is only what distributional parameters best fit; e.g., assume the data is Normally distributed and find the mean and variance: \ref apop_estimate(\c your_data, \ref apop_normal).

There are also linear models like \ref apop_ols, \ref apop_probit, and \ref apop_logit (where those last two are multinomial or binomial as the input data dictates). As in the example, they are on equal footing with the distributions, so nothing keeps you from making random draws from an estimated linear model.

Simulation models seem to not fit this form, but you will see below that if you can write an objective function for the \c p method of the model, you can use the above tools. Notably, you can estimate parameters via maximum likelihood and then give confidence intervals around those parameters.

But some models have to break uniformity, like how a histogram has a list of bins that makes no sense for a Normal distribution. These are held in <em>settings groups</em>, which you will occasionally need to tweak to modify how a model is handled or estimated. The most common example would be for maximum likelihood, eg.

\code
//Probit uses MLE. Redo the estimation using Newton's Method
Apop_settings_add_group(the_estimate, apop_mle, .verbose='y', 
                        .tolerance=1e-4, .method=APOP_RF_NEWTON);
apop_model *re_est = apop_estimate(data, the_estimate);
\endcode

See below for the details of using settings groups.

Outlineheader modelparameterization  Parameterizing or initializing a model

The models that ship with Apophenia have the requisite procedures for estimation,
making draws, and so on, but have <tt>params==NULL</tt> and <tt>settings==NULL</tt>. The
model is thus, for many purposes, incomplete, and you will need to take some action to
complete the model. There are several possibilities:

\li Estimate it! Almost all models can be sent with a data set as an argument to the
<tt>apop_estimate</tt> function. The input model is unchanged, but the output model
has parameters and settings in place. [When settings must be set beforehand, the
model documentation will tell you.]

\li If your model has a fixed number of numeric parameters, then you can set them with
\ref apop_model_set_parameters.

\li If your model has a variable number of parameters, you can directly set the \c
parameters element via \c apop_data_falloc.  For most purposes, you will also need to
set the \c msize1, \c msize2, \c vsize, and \c dsize elements to the size you want. See
the example below.

\li Some models have disparate, non-numeric settings rather than a simple matrix of
parameters. For example, an kernel density estimate needs a model as a kernel and a
base data set, which can be set via \ref apop_model_copy_set.

Here is an example that shows the options for parameterizing a model. After each
parameterization, 20 draws are made and written to a file named draws-[modelname].

\include ../eg/parameterization.c

endofdiv

Where to from here? See the \ref models page for a list of the canned models,
along with a list of basic functions that make use of them. 

If you need to write a new model, see \ref modeldetails.

endofdiv

Outlineheader modellist  Models that ship with Apophenia

This is a partial list---the full list is on the \ref models page.

Outlineheader Dist Distributions

\li\ref apop_bernoulli
\li\ref apop_beta
\li\ref apop_beta_from_mean_var()
\li\ref apop_binomial
<!-- \li\ref apop_chi_squared -->
\li\ref apop_exponential 
<!-- \li\ref apop_f_distribution -->
\li\ref apop_gamma
\li\ref apop_gaussian
\li\ref apop_improper_uniform
\li\ref apop_lognormal
\li\ref apop_multinomial
\li\ref apop_multivariate_normal
\li\ref apop_normal
\li\ref apop_poisson
\li\ref apop_t_distribution
\li\ref apop_uniform
\li\ref apop_wishart
\li\ref apop_yule
\li\ref apop_zipf

endofdiv

Outlineheader GLM  GLM family

\li\ref apop_iv
\li\ref apop_logit (incl. multinomial logit)
\li\ref apop_multinomial_probit
\li\ref apop_ols
\li\ref apop_probit
\li\ref apop_wls

endofdiv

Outlineheader moremodels More

\li\ref apop_pmf  (see the histogram section below)
\li\ref apop_kernel_density

endofdiv

endofdiv

Outlineheader mathmethods Model methods

\li\ref apop_estimate() : estimate the parameters of the model with data.
\li\ref apop_predict() : the expected value function.
\li\ref apop_draw() : random draws from an estimated model.
\li\ref apop_p() : the probability of a given data set given the model.
\li\ref apop_log_likelihood() : the log of \ref apop_p
\li\ref apop_score() : the derivative of \ref apop_log_likelihood
\li\ref apop_model_print() : display to screen

\li\ref apop_model_copy() : duplicate a model
\li\ref apop_model_set_parameters() : Models ship with no parameters set. Use this to convert a Normal(μ, σ) with unknown μ and σ into a Normal(0, 1), for example.
\li\ref apop_model_free()
\li\ref apop_model_clear(), apop_prep() : remove the parameters from a parameterized model. Used infrequently.
\li\ref apop_model_draws() : many random draws from an estimated model.

endofdiv

Outlineheader Update Filtering & updating

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

\li\ref apop_update() : Bayesian updating
\li\ref apop_model_coordinate_transform() : apply an invertible transformation to the data space
\li\ref apop_model_dconstrain() : constrain the data space of a model to a subspace. E.g., truncate a Normal distribution so \f$x>0\f$.
\li\ref apop_model_fix_params() : hold some parameters constant
\li\ref apop_model_mixture() : a linear combination of models
\li\ref apop_model_stack() : If \f$(p_1, p_2)\f$ has a Normal distribution and \f$p_3\f$ has an independent Poisson distribution, then \f$(p_1, p_2, p_3)\f$ has an <tt>apop_model_stack(apop_normal, apop_poisson)</tt> distribution.
\li\ref apop_model_dcompose() : use the output of one model as a data set for another

endofdiv

Outlineheader modelsettings Settings groups

[For info on specific settings groups and their contents and use, see the \ref settings page.]


Describing a statistical, agent-based, social, or physical model in a standardized form is difficult because every model has significantly different settings. E.g., an MLE requires a method of search (conjugate gradient, simplex, simulated annealing), and a histogram needs the number of slots to be filled with data.

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

Notice that we don't need the \c _settings ending to the settings group's name---macros
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


\li \ref Apop_settings_set, for modifying a single setting, doesn't use the nifty designated initializers format.

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

For just using a model, that's all of what you need to know. For details on writing a new settings group, see \ref settingswriting .

\li\ref Apop_settings_add_group
\li\ref Apop_settings_set
\li\ref Apop_settings_get  get a single element from a settings group.
\li\ref Apop_settings_get_group get the whole settings group.

endofdiv

endofdiv

endofdiv


Outlineheader Test Tests & diagnostics

Just about any hypothesis test consists of a few common steps:

\li  specify a statistic
\li  State the statistic's hypothesized distribution
\li  Find the odds that the statistic would lie within some given range, like <em>greater than zero</em> or <em>near 1.1</em>.

If the statistic is from a common form, like the parameters from an OLS regression, then the commonly-associated \f$t\f$ test is probably included as part of the estimation output, typically as a page appended to the \c parameters \ref apop_data set.

Some tests, like ANOVA, produce a statistic using a specialized procedure, so Apophenia includes some functions, like \ref apop_test_anova_independence and \ref apop_test_kolmogorov, to produce the statistic and look up its significance level.

If you are producing a statistic that you know has a common form, like a central limit theorem tells you that your statistic is Normally distributed, then the convenience function \ref apop_test will do the final lookup step of checking where your statistic lies on your chosen distribution.

\li\ref apop_test()
\li\ref apop_paired_t_test()
\li\ref apop_f_test()
\li\ref apop_t_test()
\li\ref apop_test_anova_independence()
\li\ref apop_test_fisher_exact()
\li\ref apop_test_kolmogorov()
\li\ref apop_plot_qq()
\li\ref apop_plot_triangle()
\li\ref apop_estimate_coefficient_of_determination()
\li\ref apop_estimate_r_squared()
\li\ref apop_estimate_parameter_tests()

See also the in-database tests above.

Outlineheader Mont Monte Carlo methods

\li\ref apop_bootstrap_cov()
\li\ref apop_jackknife_cov()

endofdiv

To give another example of testing, here is a function that used to be a part of Apophenia, but seemed a bit out of place. Here it is as a sample:


\code
// Input: any old vector. Output: 1 - the p-value for a chi-squared
// test to answer the question, "with what confidence can I reject the
// hypothesis that the variance of my data is zero?"

double apop_test_chi_squared_var_not_zero(const gsl_vector *in){
  Apop_assert_c(in, 0, 0, "input vector is NULL. Doing nothing.");
  gsl_vector	*normed;
  double 		sum=0;
    apop_vector_normalize((gsl_vector *)in, &normed, 1);
    gsl_vector_mul(normed, normed);
    for(size_t i=0;i< normed->size; 
            sum +=gsl_vector_get(normed,i++));
    gsl_vector_free(normed);
    return gsl_cdf_chisq_P(sum,in->size); 
}
\endcode

    Or, consider the Rao statistic, 
    \f${\partial\over \partial\beta}\log L(\beta)'I^{-1}(\beta){\partial\over \partial\beta}\log L(\beta)\f$
    where \f$L\f$ is your model's likelihood function and \f$I\f$ its information matrix. In code:

\code
apop_data * infoinv = apop_model_numerical_covariance(data, your_model);
apop_data * score;
apop_score(data, score->vector, your_model);
apop_data * stat    = apop_dot(apop_dot(score, infoinv), score);
\endcode

Given the correct assumptions, this is \f$\sim \chi^2_m\f$, where \f$m\f$ is the dimension of \f$\beta\f$, so the odds of a Type I error given the model is:

\code
double p_value = apop_test(stat, "chi squared", beta->size);
\endcode

endofdiv

Outlineheader Histosec Empirical distributions and PMFs (probability mass functions)

The \ref apop_pmf model wraps a \ref apop_data set so it can be read as an empirical
model, with a likelihoood function (equal to the associated weight for observed
values and zero for unobserved values), a random number generator (which consists
simply of making weighted random draws from the data), and so on.  Setting it up counts as a
model estimation from data like any other, done via \ref apop_estimate(\c your_data,
\ref apop_pmf).

You have the option of cleaning up the data before turning it into a PMF. For example...

\code
apop_data_pmf_compress(your_data);
apop_data_sort(your_data);
apop_vector_normalize(your_data->weights);
apop_model *a_pmf = apop_estimate(your_data, apop_pmf);
\endcode

...would remove all duplicates and set the weights column to reflect the original
number of copies for each observation, sort the data, normalize the weights to sum
to one, and build a PMF model around the data. Compression produces a corresponding
improvement in efficiency when calculating CDFs or making draws, but is otherwise
not necessary. Sorting or normalizing is not necessary for making draws or getting a
likelihood or log likelihood.

It is the weights vector that holds the density represented by each row; the rest of the row represents the coordinates of that density. If the input data set has no \c weights segment, then I assume that all rows have equal weight.

Most models have a \c parameters \ref apop_data set that is filled when you call \ref
apop_estimate. For a PMF model, the \c parameters are \c NULL, and the \c data itself is
used for calculation. Therefore, if you modify the data, the model will be confused. If
you modify the data, throw away any existing PMFs (\c apop_model_free) and re-estimate
a new one.


Using \ref apop_data_pmf_compress puts the data into one bin for each unique value in the data set. You may instead want bins of fixed with, in the style of a histogram. A binspec has as many columns as the data set being binned, and has one or two rows. The first row gives a bin width for the given column, so each column may have a different bin width. The second row, if present, gives the offset from zero for the bins. All bins, both above and below zero, are shifted accordingly. If there is no second row, the offset is zero.

Send this binspec to \ref apop_data_to_bins. If you send a \c NULL binspec, then the offset is zero and the bin size is big enough to ensure that there are \f$\sqrt(N)\f$ bins from minimum to maximum. There are other preferred formulæ for bin widths that minimize MSE, which might be added as a future extension.  The binspec will be added as a page to the data set, named <tt>"<binspec>"</tt>

There are a few ways of testing the claim that one distribution equals another, typically an empirical PMF versus a smooth theoretical distribution. In both cases, you will need two distributions based on the same binspec. 

For example, if you do not have a prior binspec in mind, then you can use the one generated by the first call to the histogram binning function to make sure that the second data set is in sync:

\code
apop_data_to_bins(first_set, NULL);
apop_data_to_bins(second_set, apop_data_get_page(first_set, "<binspec>"));
\endcode

You can use \ref apop_test_kolmogorov or \ref apop_histograms_test_goodness_of_fit to generate the appropriate statistics from the pairs of bins.

Kernel density estimation will produce a smoothed PDF. See \ref apop_kernel_density for details.
Or, use \ref apop_vector_moving_average for a simpler smoothing method.


Here are some functions to deal with the \ref apop_pmf.

\li\ref apop_data_pmf_compress() : merge together redundant rows in a data set before calling 
                \ref apop_estimate(\c your_data, \ref apop_pmf); optional.
\li\ref apop_vector_moving_average() : smooth a vector (e.g., your_pmf->data->weights) via moving average.
\li\ref apop_histograms_test_goodness_of_fit() : goodness-of-fit via \f$\chi^2\f$ statistic
\li\ref apop_test_kolmogorov() : goodness-of-fit via Kolmogorov-Smirnov statistic

endofdiv

Outlineheader Maxi Maximum likelihood methods

This section includes some notes on the maximum likelihood routine. As in the section
on writing models above, if a model has a \c p or \c log_likelihood function, then
calling \c apop_estimate(your_data, your_model) executes the default estimation routine
of maximum likelihood.

If you are a not a statistician, then there are a few things you will need to keep in
mind:

\li Physicists, pure mathematicians, and the GSL minimize; economists, statisticians, and
Apophenia maximize. If you are doing a minimization, be sure that your function returns minus the objective
function's value.

\li The overall setup is about estimating the parameters of a model with data. The user
provides a data set and an unparameterized model, and the system tries parameterized
models until one of them is found to be optimal. The data is fixed. The parameters make sense only in the
context of a model, so the optimization sends around a series of parameterized models,
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

The <tt>.vsize=2</tt> specified that your parameters are a vector of size two, which
means that <tt>in->parameters->vector->data</tt> is the list of <tt>double</tt>s that you
should send to \c banana. The \c more element of the structure is designed to hold any
arbitrary structure; if you use it, you will also need to use the \c more_size element, as
above. Notice that \c main in this example makes heavy use of standard C designated initializers, because
they are wonderful; you may need to tell your compiler to make use of them, such as using
the <tt>--std=gnu99</tt> flag on the \c gcc command line.

\li The function is named <em>maximum likelihood</em>, and for computational reasons,
it maximizes the log of the likelihood function. You can define your objective function
to be either the <tt>p</tt> or the <tt>log_likelihood</tt> element of your model, and
there is really only one difference: if you give me a <tt>p</tt> function, I may take
its log somewhere along the way, so the return value from <tt>p</tt> must always be non-negative.

\li Statisticians want the covariance and basic tests about the parameters. If you just
want the optimal value, then adding this line will shut off all auxiliary calculations:
\code
Apop_settings_add_group(your_model, apop_parts_wanted);
\endcode
See the documentation for \ref apop_parts_wanted_settings for details about how this
works.  It can also offer quite the speedup: especially for high-dimensional problems,
finding the covariance matrix without any information can take dozens of evaluations
of the objective funtion for each evaluation that is part of the search itself.

\li MLEs have an especially large number of parameter tweaks that could be made; see the section on MLE settings above or the \ref apop_mle_settings page.

\li Putting it all together, here is a full program to minimize Rosenbrock's banana
function. There are some extras: it uses two methods (notice how easy it is to re-run an
estimation with an alternate method, but the syntax for modifying a setting differs from
the initialization syntax) and checks that the results are accurate (because this, like
most of the code in the documentation, is part of Apophenia's test suite).

\include ../eg/banana.c

\li If you would like to see what the optimizer did, add <tt>.want_path='y'</tt> to the settings group, then get the <tt>path</tt> element from the settings group:

\code
Apop_settings_add_group(your_model, apop_mle, .want_path='y');
apop_model *out = apop_estimate(your_data, your_model);
apop_data_show(Apop_settings_get(out, apop_mle, path));
\endcode

Outlineheader constr Setting Constraints

The problem is that the parameters of a function must not take on certain values, either because the function is undefined for those values or because parameters with certain values would not fit the real-world problem.

The solution is to rewrite the function being maximized such that the function is continuous at the constraint boundary but takes a steep downward slope. The unconstrained maximization routines will be able to search a continuous function but will never return a solution that falls beyond the parameter limits.

If you give it a likelihood function with no regard to constraints plus an array of constraints, 
\ref apop_maximum_likelihood will combine them to a function that fits the above description and search accordingly.

A constraint function must do three things:
\li It must check the constraint, and if the constraint does not bind (i.e. the parameter values are OK), then it must return zero.
\li If the constraint does bind, it must return a penalty, that indicates how far off the parameter is from meeting the constraint.
\li if the constraint does bind, it must set a return vector that the likelihood function can take as a valid input. The penalty at this returned value must be zero.

The idea is that if the constraint returns zero, the log likelihood function will
return the log likelihood as usual, and if not, it will return the log likelihood at
the constraint's return vector minus the penalty. To give a concrete example, here
is a constraint function that will ensure that both parameters of a two-dimensional
input are both greater than zero. As with the constraints for many of the models that
ship with Apophenia, it is a wrapper for \ref apop_linear_constraint.

\code
static long double beta_zero_greater_than_x_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_2
    static apop_data *constraint = NULL;
    if (!constraint) constraint= apop_data_falloc((1,1,1), 0, 1);
    return apop_linear_constraint(v->parameters->vector, constraint, 1e-3);
}
\endcode

\li\ref apop_linear_constraint 

endofdiv

\li\ref apop_estimate_restart : Restarting an MLE with different settings can improve results.
\li\ref apop_maximum_likelihood()
\li\ref apop_model_numerical_covariance()
\li\ref apop_numerical_gradient()

endofdiv


Outlineheader Miss Missing data

\li\ref apop_data_listwise_delete()
\li\ref apop_ml_impute()

endofdiv

Outlineheader Legi Legible output

Most of the printing options are intended to operate in different modes
(as appropriate). The idea is that you will probably working in one
specific mode for a while---writing to screen, a file, or the
database---and so can just set one global option reflecting this:

\code 
apop_opts.output_type = 's'; //Stdout
apop_opts.output_type = 'f'; //named file
apop_opts.output_type = 'p'; //a pipe or already-opened file
apop_opts.output_type = 'd'; //the database
\endcode

C makes minimal distinction between pipes and files, so you can set a
pipe or file as output and send all output there until further notice:

\code
apop_opts.output_type = 'p';
apop_opts.output_pipe = popen("gnuplot", "w");
apop_plot_lattice(...);
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

\li\ref apop_data_print()  
\li\ref apop_matrix_print()
\li\ref apop_vector_print()
\li\ref apop_data_show() : alias for \ref apop_data_print limited to \c stdout.

The plot functions produce output for Gnuplot (so output type = \c 'd' again does not
make sense). As above, you can pipe directly to Gnuplot or write to a file. Please
consider these to be deprecated, as there is better graphics support in the works.

\li\ref apop_plot_histogram()
\li\ref apop_plot_line_and_scatter()
\li\ref apop_plot_lattice()

endofdiv

Outlineheader moreasst Assorted

A few more descriptive methods:

\li\ref apop_matrix_pca : Principal component analysis
\li\ref apop_anova : One-way or two-way ANOVA tables


General utilities:

\li\ref Apop_stopif : Apophenia's error-handling and warning-printing macro. 
\li\ref apop_opts : the global options
\li\ref apop_regex() : friendlier front-end for POSIX-standard regular expression searching and pulling matches into a \ref apop_data set.
\li\ref apop_text_paste()
\li\ref apop_system()

Math utilities:

\li\ref apop_matrix_is_positive_semidefinite()
\li\ref apop_matrix_to_positive_semidefinite()
\li\ref apop_generalized_harmonic()
\li\ref apop_multivariate_gamma()
\li\ref apop_multivariate_lngamma()
\li\ref apop_rng_alloc()

Outlineheader Prob Deprecated

These functions will probably disappear or be replaced soon.

\li\ref apop_data_to_db() Use \ref apop_data_print with <tt>.output_type='d'</tt>

endofdiv

endofdiv

Outlineheader links Further references

For your convenience, here are links to some other libraries you are probably using.

\li <a href="http://www.gnu.org/software/libc/manual/html_node/index.html">The standard C library</a>
\li <a href="http://modelingwithdata.org/c_precedence.html">The C operator precedence table</a>
\li <a href="http://www.gnu.org/software/gsl/manual/html_node/index.html">The
GSL documentation</a>, and <a href="http://www.gnu.org/software/gsl/manual/html_node/Function-Index.html">its index</a>
\li <a href="http://sqlite.org/lang.html">SQL understood by SQLite</a>

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

\li Comment out "alloc.h" in apophenia-0.22/vasprintf/vasnprintf.c
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


/** \page optionaldetails Implementation of optional arguments 
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
too much from straight C compilation for this to work. Just set \c
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
    Apop_assert(v, "You sent me a NULL vector.");
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
\c apop_varad_head macro just reduces redundancy, and will expand to
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
Thus, the macro declares each not-in-struct variable, and so there will need to be one such declaration line for each argument. Apart from requiring declarations, you can be creative: include sanity checks, post-vary the variables of the inputs, unpack without the macro, and so on. That is, this parent function does all of the bookkeeping, checking, and introductory shunting, so the base function can just do the math. Finally, the introductory section will call the base function.

The setup goes out of its way to leave the \c _base function in the public namespace, so that those who would prefer speed to bounds-checking can simply call that function directly, using standard notation. You could eliminate this feature by just merging the two functions.


<b>The m4 script</b>

The above is all you need to make this work: the varad.h file, and the above structures. But there is still a lot of redundancy, which can't be eliminated by the plain C preprocessor.

Thus, in Apophenia's code base (the one you'll get from checking out the git repository, not the gzipped distribution that has already been post-processed) you will find a pre-preprocessing script that converts a few markers to the above form. Here is the code that will expand to the above C-standard code:

\code
//header file
APOP_VAR_DECLARE void apop_vector_increment(gsl_vector * v, int i, double amt);

//code file
APOP_VAR_HEAD void apop_vector_increment(gsl_vector * v, int i, double amt){
    gsl_vector * apop_varad_var(v, NULL);
    Apop_assert(v, "You sent me a NULL vector.");
    int apop_varad_var(i, 0);
    double apop_varad_var(amt, 1);
APOP_VAR_END_HEAD
	v->data[i * v->stride]	+= amt;
}
\endcode

It is obviously much shorter. The declaration line is actually a C-standard declaration with the \c APOP_VAR_DECLARE preface, so you don't have to remember when to use semicolons. The function itself looks like a single function, but there is again a marker before the declaration line, and the introductory material is separated from the main matter by the \c APOP_VAR_END_HEAD line. Done right, drawing a line between the introductory checks or initializations and the main function can really improve readability.

The m4 script inserts a <tt>return function_base(...)</tt> at the end of the header
function, so you don't have to. If you want to call the funtion before the last line, you
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

/** \page dataprep Data prep rules 

There are a lot of ways your data can come in, and we would like to run estimations on a reasonably standardized form.

First, this page will give a little rationale, which you are welcome to skip, and then will present the set of rules.

  \section dataones Dealing with the ones column
Most standard regression-type estimations require or generally expect
a constant column. That is, the 0th column of your data is a constant (one), so the first parameter
\f$\beta_1\f$ is slightly special in corresponding to a constant rather than a variable.

However, there are some estimations that do not use the constant column.

\em "Why not implicitly assume the ones column?"
Some stats packages implicitly assume a constant column, which the user never sees. This violates the principle of transparency
upon which Apophenia is based, and is generally annoying.  Given a data matrix \f$X\f$ with the estimated parameters \f$\beta\f$, 
if the model asserts that the product \f$X\beta\f$ has meaning, then you should be able to calculate that product. With a ones column, a dot product is one line \c apop_dot(x, your_est->parameters, 0, 0)); without a ones column, the problem is left as an unpleasant exercise for the reader.

  \subsection datashunt Shunting columns around.

Each regression-type estimation has one dependent variable and several independent. In the end, we want the dependent variable to be in the vector element. However, continuing the \em "lassies faire" tradition, doing major surgery on the data, such as removing a column and moving in all subsequent columns, is more invasive than an estimation should be.

\subsection datarules The rules
So those are the two main considerations in prepping data. Here are the rules, intended to balance those considerations:

\subsubsection datamatic The automatic case
There is one clever trick we can use to resolve both the need for a ones column and for having the dependent column in the vector: given a data set with no vector element and the dependent variable in the first column of the matrix, we can copy the dependent variable into the vector and then replace the first column of the matrix with ones. The result
fits all of the above expectations.

You as a user merely have to send in a \c apop_data set with no vector and a dependent column in the first column.

\subsubsection dataprepped The already-prepped case
If your data has a vector element, then the prep routines won't try to force something to be there. That is, they won't move anything, and won't turn anything into a constant column. If you don't want to use a constant column, or your data has already been prepped by an estimation, then this is what you want.

You as a user just have to send in a \c apop_data set with a filled vector element.
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

\par The opening example

To my knowledge, Apophenia is the only open source system for doing statistical
analysis that is not a walled-garden stats package. Its data structures, syntax, and use 
reflect its unique position.

Setting aside the more advanced applications and model-building tasks, let us begin with
the workflow of a typical fitting-a-model project using Apophenia's tools:

\li Read the raw data into the database using \ref apop_text_to_db.
\li Use SQL queries handled by \ref apop_query to massage the data as needed.
\li Use \ref apop_query_to_data to pull some of the data into an in-memory \ref apop_data set.
\li Call a model estimation such as \code apop_estimate (data_set, apop_ols)\endcode  or \code apop_estimate (data_set, apop_probit)\endcode to fit parameters to the data. This will return an \ref apop_model with parameter estimates.
\li Interrogate the returned estimate, by dumping it to the screen with \ref apop_model_print, sending its parameters and variance-covariance matrices to additional tests (the \c estimate step runs a few for you), or send the model's output to be input to another model.

Here is a concrete example of most of the above steps, which you can compile and run. By the time you get to the end
of this introduction, you will have a good idea of what every line of code is doing and why.

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

If you saved the code to <tt>sample.c</tt>, then you can compile it with
\code
gcc sample.c -std=gnu99 -lapophenia -lgsl -lgslcblas -lsqlite3 -o run_me
\endcode

and then run it with <tt>./run_me</tt>. This compile line will work on any system with all the requisite tools,
but for full-time work with this or any other C library, you will probably want to write a \ref makefile .

The results should be unremarkable---this is just an ordinary regression on too few data
points---but it does give us some lines of code to dissect. 

The first two lines in \c main() make use of a database.  
I'll discuss the value of the database step more at the end of this page, but for
now, note that there are several functions, \ref apop_query, and \ref
apop_query_to_data being the ones you will most frequently be using, that will allow you to talk to and pull data from either an SQLite or mySQL/mariaDB database. 

\par Designated initializers

If anything, this line in the above sample program---

\code
apop_text_to_db(.text_file="data", .tabname="d");
\endcode

---demonstrates Apophenia's intent to balance the traditional stats package with the C library. Most C code doesn't implement variable-length argument lists or named arguments, though its only real drawback is that it bucks tradition. But this form of function appears often in the Apophenia library. It makes coding easier, less error-prone, and more pleasant.

To give another example, the \ref apop_data set has the usual row
and column numbers, but also row and column names. So you should be able to refer to a
cell by any combination of name or number; for the data set you read in above, which
doesn't have row names but does have column names, all of the following work:

\code
x = apop_data_get(data, 2, 3); //x == 0
x = apop_data_get(data, .row=2, .colname="X_3"); // x == 0
apop_data_set(data, 2, 3, 18);
apop_data_set(data, .colname="X_3", .row=2, .val= 18);
\endcode


Default values mean that the \ref apop_data_get, \ref apop_data_set, and \ref apop_data_ptr functions handle matrices, vectors, and scalars sensibly:
\code
//Let v be a vector or a one-column matrix:
double x1 = apop_data_get(v, 1);
apop_data_set(v, 2, .val=x1);
\endcode

\code
//let s be a scalar stored in a 1x1 apop_data set:
double *scalar = apop_data_ptr(s);
\endcode

Those of you coming from stats packages will enjoy working in C without having to give
up the convenience of not having to think about arguments that will take their default
value.  C traditionalists can rest assured that the Apophenia internals build this with
entirely standards-compliant C. Everybody can see the \ref designated page for details.


\section apop_data

It turns out that a lot of real-world data processing is about quotidian annoyances about text versus numeric data or dealing with missing values, and the the \ref apop_data set and its many support functions are intended to make data processing in C easy. Some users of Apophenia use the library only for its \ref apop_data set and associated functions. See the "data sets" section of the outline page (linked from the header of this page) for extensive notes on using the structure.

The structure includes seven parts:

\li a vector
\li a matrix
\li a grid of text elements
\li a vector of weights
\li names for everything: row names, a vector name, matrix column names, text names.
\li a link to a second page of data
\li an error marker

This is not a generic and abstract ideal, but is really the sort of mess that data sets look like. For
example, here is some data for a weighted OLS regression. It includes an outcome
variable in the vector, dependent variables in the matrix and text grid,
replicate weights, and column names in bold labeling the variables:

\htmlinclude apop_data_fig.html

As per the example above, Apophenia will generally assume that one row across all of these elements
describes a single observation or data point.

See above for some examples of getting and setting individual elements.

Also, \ref apop_data_get, \ref apop_data_set, and \ref apop_data_ptr consider the vector to be the -1st column,
so using the data set in the figure, \ref apop_data_get<tt>(sample_set, .row=0, .col=-1) == 1</tt>.

\par Reading in data

As per the example above, use \ref apop_text_to_data or \ref apop_text_to_db and then \ref apop_query_to_data.

\par Subsets

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

Notice that all of the slicing macros return nothing, so there is nothing to do with one
but put it on a line by itself. This limits the number of misuses.

\par Basic manipulations

The outline page (which you can get to via a link next to the snowflake header at the top of every page on this site) lists a number of other manipulations of data sets, such as 
\ref apop_data_listwise_delete for quick-and-dirty removal of observations with <tt>NaN</tt>s,
\ref apop_data_split / \ref apop_data_stack,
or \ref apop_data_sort to sort all elements by a single column.

\par Apply and map

If you have an operation of the form <em>for each element of my data set, call this
function</em>, then you can use \ref apop_map to do it. You could basically do everything you
can do with an apply/map function via a \c for loop, but the apply/map approach is clearer
and more fun. Also, if you set OpenMP's <tt>omp_set_num_threads(N)</tt> for any \c N
greater than 1 (and it is set accordingly by default on most systems), then the work
of mapping will be split across multiple CPU threads.  See the outline \> data sets \>
map/apply section for a number of examples.

\par Text

Text in C is annoying. C already treats strings as pointer-to-characters, so a grid of text data is a pointer-to-pointer-to-pointer-to-character, which is as confusing as can be.

The text grid in the \ref apop_data structure actually takes this form, but functions are provided so that most or all the pointer work is handled for you.  The \ref apop_text_alloc function is really a realloc function: you can use it to resize the text grid as necessary. The \ref apop_text_add function will do the pointer work in copying a single string to the grid. Functions that act on entire data sets, like \ref apop_data_rm_rows, handle the text part as well.

You have <tt>your_data->textsize[0]</tt> rows and <tt>your_data->textsize[1]</tt> columns. If you are using only the functions to this point, then empty elements are a blank string (<tt>""</tt>), not \c NULL.
For reading individual elements, refer to the \f$(i,j)\f$th text element via <tt>your_data->text[i][j]</tt>.

\par Errors

Many functions will set the <tt>error</tt> element of the \ref apop_data structure being operated on if anything goes wrong. You can use this to halt the program or take corrective action:

\code 
apop_data *the_data = apop_query_to_data("select * from d");
if (the_data->error) exit(0);
\endcode 

\par The whole structure

Here is a diagram of all of Apophenia's structures and how they
relate. It is taken from this
<a href="http://modelingwithdata.org/pdfs/cheatsheet.pdf">cheat sheet</a> (2 page PDF),
which will be useful to you if only because it lists some of the functions that act on
GSL vectors and matrices that are useful (in fact, essential) but out of the scope of the Apophenia documentation.

\image html http://apophenia.info/structs.png
\image latex structs.png


All of the elements of the \ref apop_data structure are laid out at middle-left. You have
already met the vector, matrix, and weights, which are all a \c gsl_vector or \c gsl_matrix.

The diagram shows the \ref apop_name structure, which has received little mention so far because names
basically take care of themselves. Just use \ref apop_data_add_names to add names to your data
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
we want to do with our models: estimating the parameters of a model (like the mean and
variance of a Normal distribution) given data, or drawing random numbers, or showing the
expected value, or showing the expected value of one part of the data given fixed values
for the rest of it. The \ref apop_model is intended to encapsulate most of these desires
into one object, so that models can easily be swapped around, modified to create new models, 
compared, and so on.

From the figure above, you can see that the \ref apop_model structure is pretty big,
including a number of informational items, key being the \c parameters, \c data, and \c
info elements, a list of settings to be discussed below, and a set of procedures for many
operations.  Its contents are not (entirely) arbitrary: the theoretical basis for what is and is not included in an \ref apop_model, as
well as its overall intent, are described in this <a
href="http://ben.klemens.org/klemens-model_objects.pdf">white paper</a>.

There are helper functions that will allow you to avoid deailing with the model internals. For example, the \ref apop_estimate helper function means you never have to look at the model's \c estimate method (if it even has one), and you will simply pass the model to a function, as with the above form:

\code
    apop_model *est = apop_estimate(data, apop_ols);
\endcode

\li Apophenia ships with a broad set of models, like \ref apop_ols, \ref apop_dirichlet,
    \ref apop_loess, and \ref apop_pmf (probability mass function); see the full list on <a href="http://apophenia.info/group__models.html">the models documentation page</a>. You would estimate the
parameters of any of them using the form above, with the appropriate model in the second
slot of the \ref apop_estimate call.
\li The models that ship with Apophenia, like \ref apop_ols, include the procedures and some metadata, but are of course not yet estimated using a one data set (i.e., <tt>parameters == NULL</tt>). The line above generated a new
model, \c est, which is identical to the base OLS model but has estimated parameters
(and covariances, and basic hypothesis tests, a log likelihood, AIC, BIC, et cetera). 
\li You will mostly use the models by passing them as inputs to 
functions like \ref apop_estimate, \ref apop_draw, or \ref apop_predict; more examples below.
After \ref apop_estimate, most require a parameterized model like \c est. After all, it doesn't make sense to
draw from a Normal distribution until you've specified its mean and standard deviation.
\li You can use \ref apop_model_print to print the various elements to screen.
\li You can combine and transform models with functions such as \ref apop_model_fix_params, \ref apop_model_coordinate_transform, or \ref apop_model_mixture. Each of these functions produce a new model, which can be estimated, re-combined, or otherwise used like any other model.
\li Writing your own models won't be covered in this introduction, but it can be pretty easy to
copy and modify the procedures of an existing model to fit your needs. When in doubt, delete a procedure, because any procedures that are missing will have
defaults filled when used by functions like \ref apop_estimate (which uses \ref
apop_maximum_likelihood) or \ref apop_cdf (which uses integration via random draws). See \ref modeldetails for details.
\li There's a simple rule of thumb for remembering the order of the arguments to most of
Apophenia's functions, including \ref apop_estimate : the data always comes first.

\par Settings

Every model, and every method one would apply to a model, is prone to have a list of
settings: how many bins in the histogram, at what tolerance does the maximum likelihood
search end, what are the models being combined in the mixture?

Apophenia organizes settings in <em>settings groups</em>, which are then attached to models. 
In the following snippet, we specify a Beta distribution prior, and attach a settings group to it with
information about how Bayesian updating should be done with this specific copy of the
Beta. For a likelihood, we generate an empirical distribution---a PMF---from an input 
data set; you will often use the \ref apop_pmf to turn a data set into an distribution.
When we call \ref apop_update on the last line, it already has all of the above info
on hand.

\code
apop_model *beta = apop_model_set_parameters(apop_beta, 0.5, 0.25);
Apop_settings_add_group(beta, apop_update, .burnin = 0.2, .periods =1e5);
apop_model *my_pmf = apop_estimate(your_data, apop_pmf);
apop_model *posterior = apop_update(.prior= beta, .likelihood = my_pmf);
\endcode

You will encounter model settings often when doing nontrivial work with models. All
can be set using a form like above. See the <a href="http://apophenia.info/group__settings.html">settings documentation page</a>
for the full list of options.
There is a full discussion of using and writing settings groups in the <a href="http://apophenia.info/outline.html">outline page</a> under the Models heading.


\par Databases and models

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

Now that you have \c est, an estimated model, you can interrogate it. This is really where Apophenia and its encapsulated
model objects shine, because you can do more than just admire the parameter estimates on
the screen: you can take your estimated data set and fill in or generate new data, use it
as an input to the parent distribution of a hierarchical model, et cetera. Some simple
examples:

 \code
 //If you have a new data set with missing elements (represented by NaN), you can fill in predicted values:
apop_predict(new_data_set, est);
apop_data_show(new_data_set)

 //Fill a matrix with random draws. The draw function needs an RNG from
 //the GNU Scientific Library, and a pointer-to-double.
gsl_rng *r = apop_rng_alloc(218);
apop_data *d = apop_model_draws(est, .count=1000, .rng=r);
\endcode

\par Testing

Here is the model for all testing within Apophenia:

\li Calculate a statistic.
\li Describe the distribution of that statistic.
\li Work out how much of the distribution is (above|below|closer to zero than) the statistic.

There are a handful of named tests that produce a known statistic and then compare to a
known distribution, like \ref apop_test_kolmogorov or \ref apop_test_fisher_exact. For
traditional distributions (Normal, \f$t\f$, \f$\chi^2\f$), use the \ref apop_test convenience
function.

But if you've gotten this far with your models, then you'll want to apply the above
three-step process to your model parameters. First I'll give an overview of the three steps, then
another working example.

\li Model parameters are a statistic, and you know that  \ref apop_estimate<tt>(your_data,
        your_model)</tt> will output a model with a <tt>parameters</tt> element.
\li Now, the distribution of a parameter is also a model, so 
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

Here is the extended program; the first part is identical to the program above. The second
half uses many of the above tricks: one of the inputs to \ref apop_parameter_model (which
row of the parameter set to use) is sent by adding a settings group, we pull that row
into a separate data set using \ref Apop_row, and we set its vector value by
referring to it as the -1st element.

The story of the code is as above: we estimate a model, query that model for the distribution associated with the
parameters, and then use that subsidiary model to find the area under the curve up to the parameter itself and zero. This is
more than sufficient to test the hypothesis.

\include ols2.c

Note that the procedure did not assume the model parameters had a certain form. It
queried the model for the distribution of parameter \c x_1, and if the model didn't have
a closed-form answer then a distribution via bootstrap would be provided. Then that model
was queried for its CDF. [The procedure does assume a symmetric distribution. Fixing this
is left as an exercise for the reader.] For a model like OLS, this is entirely overkill, 
which is why OLS provides the basic hypothesis tests automatically. But for more advanced models where the distribution of parameters is unknown or has no closed-form solution, this may be the only recourse.


This introduction has shown you the \ref apop_data set and some of the functions
associated, which might be useful even if you aren't formally doing statistical work but do have to deal with data with real-world elements like column names and mixed
numeric/text values. You've seen how Apophenia encapsulates as many of a model's
characteristics as possible into a single \ref apop_model object, which you can send with
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
href="http://ben.klemens.org/pdfs/klemens-modelcats.pdf"><em>A Cross-paradigm Modeling
Framework</em></a> (PDF).

This page includes:

\li \ref write_likelihoods, giving a quick overview of how to write a new model from scratch.
\li \ref settingswriting, covering the writing of <em>ad hoc</em> structures to hold model- or method-specific details, like the number of periods for burning in an MCMC run or the number of bins in a histogram.
\li \ref vtables, covering the means of writing special-case routines for functions that are not part of the \ref apop_model itself, including the score.
\li \ref modeldataparts, a detailed list of the requirements for the data (non-function) elements of an \ref apop_model.
\li \ref methodsection, a detailed list of requirements for the method (function) elements of an \ref apop_model.

\section write_likelihoods A walkthrough

Users are encouraged to always use models via the helper functions, like
\ref apop_estimate or \ref apop_cdf.  The helper functions do some boilerplate error
checking, and are where the defaults are called: if your model has a \c log_likelihood
method but no \c p method, then \ref apop_p will use exp(\c log_likelihood). If you don't
give an \c estimate method, then \c apop_estimate will call \ref apop_maximum_likelihood.

So the game in writing a new model is to write just enough internal methods to give the helper functions what they need.
In the not-uncommon best case, all you need to do is write a log likelihood function.

Here is how one would set up a model that could be estimated using maximum likelihood:

\li Write a likelihood function. Its header will look like this:

\code
long double new_log_likelihood(apop_data *data, apop_model *m);
\endcode

where \c data is the input data, and \c
m is the parametrized model (i.e. your model with a \c parameters element set by the caller). 
This function will return the value of the log likelihood function at the given parameters.

\li Is this a constrained optimization? See the outline page under maximum likelihood methods \f$->\f$ Setting constraints on how to set them. Otherwise, no constraints will be assumed.
\li Write the object:

\code
apop_model *your_new_model = &(apop_model){"The Me distribution", 
            .vsize=n0, .msize1=n1, .msize2=n2, .dsize=nd,
            .log_likelihood = new_log_likelihood };
\endcode

\li The first element is the human-language name for your model.
\li the \c vsize, \c msize1, and \c msize2 elements specify the shape of the parameter set. For example, if it's three numbers in the vector, then set <tt>.vsize=3</tt> and omit the matrix sizes. The default model prep routine will call 
<tt>new_est->parameters = apop_data_alloc(vsize, msize1, msize2)</tt>. 
\li The \c dsize element is the size of one random draw from your model.
\li It's common to have [the number of columns in your data set] parameters; this
count will be filled in if you specify \c -1 for \c vsize, <tt>msize(1|2)</tt>, or
<tt>dsize</tt>. If the allocation is exceptional in a different way, then you will
need to allocate parameters by writing a custom \c prep method for the model.
\li If there are constraints, add an element for those too.

You already have more than enough that something like this will work (the \c dsize is used for random draws):
\code
apop_model *estimated = apop_estimate(your_data, your_new_model);
\endcode

Once that baseline works, you can fill in other elements of the \ref apop_model as needed.

For example, if you are using a maximum likelihood method to estimate parameters, you can get much faster estimates and better covariance estimates by specifying the dlog likelihood function (aka the score):

\code
void apop_new_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
    //some algebra here to find df/dp0, df/dp1, df/dp2....
    gsl_vector_set(gradient, 0, d_0);
    gsl_vector_set(gradient, 1, d_1);
}
\endcode
The score has to be registered (see below) using 
\code
apop_score_insert(apop_new_dlog_likelihood, your_new_model);
\endcode

The random number generator for a
unidimensional case is filled in via adaptive-rejection sampling, there is a default
print and prep, et cetera. You can improve upon the defaults as your time allows.

\section settingswriting  Writing new settings groups

Your model may need additional settings or auxiliary information to function, which would require associating a model-specific struct with the model.

Before getting into the detail of how to make model-specific groups of settings work, note that there's a lightweight method of storing sundry settings, so in many cases you can bypass all of the following.

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

The structure itself gets the full name, \c ysg_settings. Everything else is a macro, and so you need only specify \c ysg, and the \c _settings part is filled in. Because of these macros, your \c struct name must end in \c _settings.

If you have an especially simple structure, then you can generate the three functions with these three macros in your <tt>.c</tt> file:

\code
Apop_settings_init(ysg, )
Apop_settings_copy(ysg, )
Apop_settings_free(ysg, )
\endcode

These macros generate appropriate functions to do what you'd expect: allocating the
main structure, copying one struct to another, freeing the main structure.  
The spaces after the commas indicate that no special code gets added to
the functions that these macros generate.

You'll never call these funtions directly; they are called by \ref Apop_settings_add_group,
\ref apop_model_free, and other model or settings-group handling functions.

Now that initializing/copying/freeing of
the structure itself is handled, the remainder of this section will be about how to
add instructions for the struture internals, like data that is pointed to by the structure elements.

\li For the allocate function, use the above form if everything in your code defaults to zero/\c NULL.  
In most cases, though, you will need a new line declaring a default for every element in your structure. There is a macro to help with this too. 
These macros will define for your use a structure named \c in, and an output pointer-to-struct named \c out.
Continuing the above example:

\code
Apop_settings_init (ysg, 
      Apop_assert(in.size1, "I need you to give me a value for size1. Stopping.");
      Apop_varad_set(size2, 10);
      Apop_varad_set(dataset, apop_data_alloc(out->size1, out->size2));
      Apop_varad_set(refs, malloc(sizeof(int)));
      *refs=1;
)
\endcode

Now, <tt>Apop_settings_add(a_model, ysg, .size1=100)</tt> would set up a group with a 100-by-10 data set, and set the owner bit to one. 

\li Some functions do extensive internal copying, so you will need a copy function even if you don't do any explicit calls to \ref apop_model_copy. The default above simply copies every element in the structure. Pointers are copied, giving you two pointers pointing to the same data. We have to be careful to prevent double-freeing later.

\code
//The elements of the set to copy are all copied, and then make one additional modification:
Apop_settings_copy (ysg,
        (*refs)++;
)
\endcode

\li The struct itself is freed by boilerplate code, but add code in the free function
to free data pointed to by pointers in the main struture. The macro defines a
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

For any given function (e.g., entropy, the dlog likelihood, Bayesian updating), there is
probably a special case for well-known models like the Normal distribution. Rather than
any procedure that could have a special-case calculation to the \c apop_model struct,
functions may maintain a registry of models and associated special-case procedures.

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
\li Calls to <tt>..._vtable_add</tt> are typically placed in the \c prep method of the given model, thus ensuring that the auxiliary functions are registered after the first time the model is sent to \ref apop_estimate.

This overview will not go into detail about setting up a new vtable. Briefly:

\li See the existing setups in vtables.h. 
\li Cut/paste one and do a search and replace to change the name to match your desired use.
\li Set the typedef to describe the functions that get added to the vtable.
\li Rewrite the hash function to check the part of the inputs that interest you. For
example, the update vtable associates functions with the \c draw, \c log_likelihood,
and \p methods of the model. A model where these elements are identical but the name
is changed will still match.

\section modeldataparts The data elements

The remainder of this page covers the detailed expectations regarding the elements
of the \ref apop_model structure. I begin with the data (non-function) elements,
and then cover the method (function) elements. Some of the following will be
requirements for all models and some will be advice to authors; I use the accepted
definitions of <a href="http://tools.ietf.org/html/rfc2119">"must", "shall", "may"</a>
and related words.

\subsection datasubsec data

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

\li Your routines may accept other data formats, as per contract with the user.
    For example, regression-type functions use a function named \c ols_shuffle
    to convert a matrix where the first column is the dependent variable to a data
    set with dependent variable in the vector and a column of ones in the first
    matrix column. By checking for a vector, the prep function knows whether to do
    the shuffling or not. Most univariate distributions take each scalar element as
    a separate data point; having one data point per row is a special case.

\subsection paramsubsec Parameters, vsize, msize1,  msize2

\li The sizes will be used by the \c prep method of the model; see below. Given the model \c m and its elements \c m.vsize, \c m.msize1, \c m.msize2,
    functions that need to allocate a parameter set will do so via <tt>apop_data_alloc(m.vsize, m.msize1, m.msize2)</tt>. 

\li As a special case, if you set any of \c .vsize, \c .msize1, or \c .msize2
to \c -1, then the default prep method will set that size to the number of columns in
the input data. This is what you want for regression methods, where there is one parameter per independent variable.


\subsection infosubsec Info

\li The first page, named \c &lt;info&gt; is typically a list of scalars. Nothing is guaranteed, but the elements may include:

\li AIC: <a href="https://en.wikipedia.org/wiki/Akaike's_Information_Criterion">Aikake Information Criterion</a>
\li AIC_c: AIC with a finite sample correction. "<b>Generally, we advocate the use of AIC_c when the ratio \f$n/K\f$ is small (say \f$< 40\f$)</b>" [Kenneth P. Burnham, David R. Anderson: <em>Model Selection and Multi-Model Inference</em>, p 66, emphasis in original.]
\li BIC: <a href="https://en.wikipedia.org/wiki/Bayesian_information_criterion">Bayesian Information Criterion</a>
\li R squared
\li R squared adj
\li log likelihood
\li status.

For those elements that require a count of input data, the calculations assume each row in the input \ref apop_data set is a single datum.

Get these via, e.g., <tt>apop_data_get(your_model->info, .rowname="log likelihood")</tt>.
When writing for any arbitrary function, be prepared to handle \c NaN, indicating that the element is not calculated or saved in the info page by the given model.

\li Several routines will include a \c predict table. The table has these rows:
    \li row (optional)
    \li col (optional)
    \li observed
    \li predicted
    \li residual

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
\li The inputs are an \ref apop_data set and an \ref apop_model, which should include a filled <tt>->parameters</tt> element.
\li We assume that the parameters have been set, by users via \ref apop_estimate or \ref apop_model_set_parameters, or by \ref apop_maximum_likelihood by its search algorithms. if the parameters are necessary, the function shall check that the parameters are not \c NULL and set the model's \c error element to \c 'p' if they are.
\li Return \c NaN on errors. If an error in the input model is found, the function may set the input model's \c error element to an appropriate \c char value.
\li If observations are assumed to be iid, you can probably use \ref apop_map_sum to write the core of the log likelihood function.
\li If your model includes both \ log_likelihood and \c p methods, it must be the case that <tt>log(p(d, m))</tt> equals <tt>log_likelihood(d, m)</tt> for all \c d and \c m.

\subsection prepsubsection prep

\li Function header looks like <tt>void your_prep(apop_data *data, apop_model *params)</tt>.
\li If \c vsize, \c msize1, or \c msize2 are -1, then the prep function will set them to the width of the input data.
\li If \c dsize is -1, then the prep function shall set it to the width of the input data.
\li If the \c parameters element is not allocated, the function shall allocate it via <tt>apop_data_alloc(vsize, msize1, msize2)</tt> (or equivalent).
\li The model's <tt>data</tt> pointer shall be set to point to the input data.
\li The input data may be modified by the prep routine. For example, the OLS prep routine shuffles a single input matrix as described above under \c data.
\li The \c info element shall be allocated and its title set to "<Info>".
\li The default is \ref apop_model_clear. It does all of the above.
\li The prep routine may initialize any desired settings groups. Unless otherwise
stated, these should not be removed if they are already there, so that users can override defaults by adding a settings group before starting an estimation.
\li If any functions associated with the model need to be added to 
a vtable (see above), the registration shall happen here. Registration may also happen elsewhere.

\subsection estimatesubsection estimate

\li Function header looks like  <tt> void your_estimate(apop_data * data, apop_model *params)</tt>.
\li Assume that the prep routine has already been run. Notably, this means that parameters have been allocated.
\li Assume that the \c parmaeters hold garbage (as in a \c malloc without a subsequent assignment to the <tt>malloc</tt>-ed space).
\li The function modifies the input model, and returns nothing. Note that this is different from the wrapper function, \ref apop_estimate, which makes a copy of its input model, preps it, and then calls the \c estimate function with the prepeped copy.
\li The function shall set the \c parameters of the input model. For consistency with other models, the estimate should be the maximum likelihood estimate, unless otherwise documented.
\li Additional settings may be set.
\li The model's \c &lt;Info&gt; page may be filled with data. For scalars like log likelihood and AIC, use \ref apop_data_add_named_elmt.
\li Data should not be modified by the \c estimate routine; any changes to the data made by \c estimate must be documented.
\li The default called by \ref apop_estimate is \ref apop_maximum_likelihood.
\li If errors occur during processing, set the model's \c error element to a single character. Documentation should include the list of error characters and their meaning.

\subsection drawsubsection draw

\li Function header looks like <tt>void your_draw(double *out, gsl_rng* r, apop_model *params)</tt>
\li Assume that model \c paramters are set, via \ref apop_estimate or \ref apop_model_set_parameters. The author of the draw method should check that \c parameters are not \c NULL and fill the output with NaNs if necessary parameters are not set.
\li User inputs a pointer-to-<tt>double</tt> of length \c dsize; user is expected to make sure that there is adequate space.  User also inputs a \c gsl_rng, already allocated (probably via \ref apop_rng_alloc).
\li The function shall fill the space pointed to by the input pointer with a random draw from the data space, where the likelihood of any given observation is proportional to its likelihood as given by the \c p method. Data shall be reduced to a single vector via \ref apop_data_pack if it is not already a single vector.

\subsection cdfsubsection cdf

\li Function header looks like <tt>long double your_cdf(apop_data *d, apop_model *params)</tt>.
\li Assume that \c paramters are set, via \ref apop_estimate or \ref apop_model_set_parameters. The author of the CDF method should check that \c parameters are not \c NULL and return NaN if necessary parameters are not set.
\li The CDF method must accept data as a single row of data in the \c matrix of the input \ref apop_data set (as per a draw produced using the \c draw method). May accept other formats.
\li Returns the percentage of the likelihood function \f$\leq\f$ the first row of the input data. The definition of \f$\leq\f$ is chosen by the model author.
\li If one is not already present, an \c apop_cdf_settings group may be added to the model. See the \ref apop_cdf function for details of its use.

\subsection constraintsubsection constraint

\li Function header looks like <tt>long double your_constraint(apop_data *data, apop_model *params)</tt>.
\li Assume that \c parameters are set, via \ref apop_estimate, \ref apop_model_set_parameters, or the internals of an MLE search. The author of the constraint method should check that \c parameters are not \c NULL and return NaN if necessary parameters are not set.
\li See \ref apop_linear_constraint for a useful basis and/or example. Many constraints can be written as wrappers for this function.
\li If the constraint is met, then return zero.
\li If the constraint fails, then (1) move the \c parameters in the input model to a
constraint-satisfying value, and (2) return the distance between the input parameters and
what you've moved the parameters to. The choice of within-bounds parameters and distance function is left to the author of the constraint function.
*/
