/* Apophenia's documentation
Copyright (c) 2005--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

/** \page eg Some examples
 Here are a few pieces of sample code, gathered from elsewhere in the documentation, for testing your installation or to give you a sense of what code with Apophenia's tools looks like. If you'd like more context or explanation, please click through to the page from which the example was taken.

In the documentation for the \ref apop_ols model, a program to read in data and run a regression. You'll need to go to that page for the sample data and further discussion.

\include ols1.c

In the \ref apop_text_to_db page, a variant, demonstrating the use of optional, named arguments:

\include ols.c

Now, from the Python section of the outline page, the same thing via Python:

\include ols.py

The same page illustrates some other calls of Apophenia functions via Python, such as the \ref apop_test_fisher_exact function.

\include fisher.py

From the \ref setup page, an example of gathering data from two processes, saving the input to a database, then doing a later analysis:

\include draw_to_db.c

A demonstration of \ref apop_plot_line_and_scatter . You'll need a
database from the {\em Modeling with Data} sample code, at
http://modelingwithdata.org/appendices.html.

\include scatter.c

In the \ref outline section on map/apply, a new \f$t\f$-test on every row, with all operations acting on entire rows rather than individual data points:

\include t_test_by_rows.c

In the documentation for \ref apop_query_to_text, a program to list all the tables in an Sqlite database.
\include ls_tables.c

Finally, a demonstration of fixing parameters to create a marginal distribution, via \ref apop_model_fix_params
\include fix_params.c

 */



/** \page setup Setting up
\section cast The supporting cast 
To use Apophenia, you will need to have a working C compiler, the GSL (v1.7 or higher) and SQLite installed. 

We've moved the setup documentation to <a href="http://modelingwithdata.org/appendix_o.html">Appendix O</a> of <em> Modeling with Data</em>. Please see that page. This site has a few more notes for \ref windows "Windows" users or \ref mingw users.

\subsection sample_program Sample programs
First, there is a short, complete program in the \ref apop_ols entry which runs a simple OLS regression on a data file. Follow
the instructions there to compile and run. 

Here is another sample program, intended to show how one would integrate Apophenia into an existing program. For example, say that you are running a simulation of two different treatments, or say that two sensors are posting data at regular intervals. You need to gather the data in an organized form, and then ask questions of the resulting data set.  Below, a thousand draws are made from the two processes and put into a database. Then, the data is pulled out, some simple statistics are compiled, and the data is written to a text file for inspection outside of the program.  This program will compile cleanly with the sample \ref makefile.

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
 \li gsl -- the engine that powers apophenia
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
\verbatim
export MY_LIBS = src   #choose a directory name to be created in your home directory.
tar xvzf pkg.tgz       #change pkg.tgz to the appropriate name
cd package_dir         #same here.
mkdir $HOME/$MY_LIBS
./configure --prefix $HOME/$MY_LIBS
make
make install   #Now you don't have to be root.
echo "export LD_LIBRARY_PATH=$HOME/$MY_LIBS:\$LD_LIBRARY_PATH" >> ~/.bashrc
\endverbatim
*/


/** \page makefile Makefile
 
Instead of giving lengthy GCC commands at the command prompt, you can use a Makefile to do most of the work. How to:
\li Copy and paste the following into a file named \c makefile.
\li Make sure the two indented lines begin with a single Tab, instead of seven or eight spaces. 
\li Change the first line to the name of your program (e.g., if you have written <tt>sample.c</tt>, then the first line will read <tt>PROGNAME=sample</tt>). 
\li If your program has multiple <tt>.c</tt> files, just add a corresponding <tt>.o</tt> on the <tt>objects</tt> line, e.g. <tt>sample2.o</tt> <tt>sample3.o</tt>
\li One you have a Makefile in the directory, simply type <tt>make</tt> at the command prompt to generate the executable.

 \code
PROGNAME = your_program_name_here
objects =$(PROGNAME).o
CFLAGS = -g -Wall 
LINKFLAGS = -lapophenia -lgsl -lgslcblas -lsqlite3
c: $(objects)
	gcc $(CFLAGS) $(objects) $(LINKFLAGS) -o $(PROGNAME)
$(objects): %.o: %.c 
	gcc $(CFLAGS) -c $< -o $@
\endcode

By the way, if your system has \c pkg-config, then you can use it for a slightly more robust and readable makefile. Replace the above C and link flags with:
\code
CFLAGS = -g -Wall `pkg-config --cflags apophenia`
LINKFLAGS = `pkg-config --libs apophenia`
\endcode
The \c pkg-config program will then fill in the appropriate directories and libraries. Pkg-config knows Apophenia depends on the GSL and (if applicable) sqlite3, so you need only list the most-dependent library.
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

\li There is some overhead to checking for defaults, which could slow the code by up to 25\%. If this is noticeable for your situation, you can pass on this convenient form and call the underlying function directly, by adding \c _base to the name and giving all arguments:

\code
apop_text_to_db_base("infile.txt", "intable", 0, 1, NULL);
\endcode

\li If one of the optional elements is an RNG, see \ref autorng on what happens when you don't provide an RNG.
\li For exhaustive details on implementation of the above (should you wish to write new functions that behave like this) see the \ref optionaldetails page.

  */

/** \page autorng Auto-allocated RNGs

Functions that use the \ref designated syntax for reading inputs and assigning default values use the following rules for handling RNGs.

- The first time a function is called with no \c gsl_rng as
input, a new \c gsl_rng is produced. The call will effectively look like this
\code  
static gsl_rng *internal_rng = gsl_rng_alloc(++apop_opts.rng_seed);
\endcode

- Because \c internal_rng is declared \c static, it will remember its state as you repeatedly call the function, so you will get appropriate random numbers.

- \c apop_opts.rng_seed is incremented at each use, so you can write down the seed used for later reference. 

- Because it increments, the next function to auto-allocate an RNG will produce different random numbers. That is, every function that uses this setup will have a different, independent RNG.

- If you would like a different outcome every time the program runs, set the seed to the time before running:
\code  
#include <time.h>
apop_opts.rng_seed = time(NULL);
\endcode  

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


/**  \mainpage An outline of the library
  \anchor outline

ALLBUTTON

Outlineheader preliminaries Overview and Preliminaries

As well as the information in this outline, there is a separate page covering the details of 
 \ref setup "setting up a computing environment" and another page with \ref eg "some sample code" for your perusal.

Outlineheader  intro A quick overview

  A typical data analysis project using Apophenia's tools would go something like this:

 \li Read the raw data into the database using \ref apop_text_to_db.
 \li Use SQL queries handled by \ref apop_query to massage the data as needed.
 \li Use \ref apop_query_to_data to pull the data into an in-memory apop_data set.
 \li Call a model estimation such as \code apop_estimate (data_set, apop_OLS)\endcode  or \code apop_estimate (data_set, apop_probit)\endcode to fit parameters to the data. This will return an \ref apop_model like the original, but with parameter estimates.
 \li Interrogate the returned estimate, by dumping it to the screen with \ref apop_model_show, sending its parameters and variance-covariance matrices to a test, et cetera.

Over the course of this,
our researcher had to parse text, use a database, manipulate matrices,
fit statistical models, and then test the results by looking up a
summary statistic in a table of distributions. Having hooks to all of
these things in one library makes the flow much more pleasant.

For a concrete example of folding Apophenia into a project, have a look at this \ref sample_program "sample program".

Of course, you may not need all this: you might just want a neat interface to SQLite, or an easy-to-use structure that holds your data 
in a matrix and allows you to label the rows and columns. Have a look at the \ref outline "library outline" for the various categories of relatively standalone functions that you can fold into your own work.

If you would like a more in-depth discussion of Apophenia's raison d'etre and logic, have a look at these
<a href="http://apophenia.info/apop_notes.html">design notes</a>.

endofdiv

Outlineheader c Some notes on C and Apophenia's use of C utilities.
 
Outlineheader learning  Learning C

<a href="http://modelingwithdata.org">Modeling with Data</a> has a full tutorial for C, oriented at users of standard stats packages.  More nuts-and-bolts tutorials are <a href="http://www.google.com/search?hl=es&amp;c2coff=1&amp;q=c+tutorial">in abundance</a>.  Some people find pointers to be especially difficult; fortunately, there's a <a href="http://cslibrary.stanford.edu/104/">claymation cartoon</a> which clarifies everything.

Coding often relies on gathering together many libraries; there is a section at the bottom of this outline linking to references for some libraries upon which Apohenia builds.

endofdiv

Outlineheader usagenotes  Usage notes

Here are some notes about the technical details of using the Apophenia library in your development enviornment.

<b> Header aggregation </b>

If you put 
\verbatim
#include <apop.h>
\endverbatim
at the top of your file, then it will call virtually every header file you could need: gsl_matrix.h, gsl_blas.h, sqlite3.h, stdio.h, string.h, math.h, apophenia_all_of_them.h, et cetera. Of course, if you get `implicit declaration of...' then you will need to manually include something else.
Bear in mind that every book on C will tell you this is bad form and you shouldn't do it.

<b> Easier calling syntax</b>

Several functions (but nowhere near all) use a script-like syntax for calling functions. For example, the \ref apop_vector_distance function finds the distance between two vectors using various metrics. The full function call, using an \f$L_3\f$ norm for vectors \f$v_1\f$ and \f$v_2\f$ would be 
\code
apop_vector_distance(v1, v2, 'L', 3);
\endcode
But some metrics don't need a number to be fully described, so you can leave that off. Here's the standard (Euclidian) distance:
\code
apop_vector_distance(v1, v2, 'E');
\endcode
Because Euclidian distance is so standard, it is assumed as the default if you don't specify a norm, so you can leave out the third input:
\code
apop_vector_distance(v1, v2);
\endcode
If you give only one vector, I'll assume the other is the zero vector, and thus return the length of the first vector. Here is the length of a vector using the Manhattan metric:
\code
apop_vector_distance(v1, .metric='M');
\endcode

See the \ref designated page for details of this syntax.

endofdiv

Outlineheader liblinking Libraries 

Your best bet is to write yourself a \ref makefile "Makefile".  If you don't want to use the sample \ref makefile "Makefile", then here are some notes for the command line.  When compiling, you will need to tell the compiler to use the Apophenia library. That is, you will need to call GCC with <tt>gcc -lapophenia</tt> (as well as the other usual flags). For example,
<tt>
gcc sample.c -lapophenia -lsqlite3 -lgsl -lgslcblas -o run_me -g 
</tt>
will produce an output file named <tt>run_me</tt> from the input source code <tt>sample.c</tt>. It will include symbols for debugging (<tt>-g</tt>) and will correctly find the functions in the Apophenia, GSL, and SQLite libraries (<tt>-lapophenia -lgsl ...</tt>). 
Order matters in the linking list: the files a package depends on should be listed after the package. E.g., since sample.c depends on Apophenia, <tt>gcc sample.c -lapophenia</tt> will work, while <tt>gcc -lapophenia sample.c</tt> is likely to give you errors. Similarly, list <tt>-lapophenia</tt> before <tt>-lgsl</tt>, which comes before <tt>-lgslcblas</tt>.

endofdiv

Outlineheader debugging  Debugging

The global variable <tt>apop_opts.verbose</tt> turns on some diagnostics, such as printing the query sent to the databse engine (which is useful if you are substituting in many <tt>\%s</tt>es). Just set <tt>apop_opts.verbose =1</tt> when you want feedback and <tt>apop_opts.verbose=0</tt> when you don't.

If you use \c gdb, you can define macros to use the pretty-printing functions on your data, which can be a significant help. Add these to your \c .gdbinit:
\code

define pv
    p apop_vector_show($arg0)
end

define pm
    p apop_matrix_show($arg0)
end

define pd
    p apop_data_show($arg0)
end

define pa
    p *($arg0)@$arg1
end
\endcode

Then just <tt>pd mydata</tt> to see all components of your data set neatly displayed, or <tt>pa myarray 5</tt> to see the first five elements of your array. 

endofdiv

Outlineheader vim Syntax highlighting 

If your text editor supports syntax highlighting, there are a few types defined in the Apophenia and GSL headers which may be worth coloring.
E.g., for <tt>vim</tt>, add the following two lines to <tt>/usr/share/vim/syntax/c.vim</tt>:
\verbatim
syn keyword     cType           gsl_matrix gsl_rng gsl_vector apop_data
syn keyword     cType           apop_name apop_model
\endverbatim
Other text editors have similar files to which you can add the above types.

endofdiv

endofdiv

 Outlineheader python The Python interface

The distribution includes a Python interface via the SWIG tool for bridging across languages. 

Installation: You will need to have the SWIG and Python-development packages installed on your computer when compiling. From there, 
\code
./configure --enable-python
make
sudo make install
\endcode 
will produce the Python
library and install it in Python's site-packages directory.

Sample script:

On the page on the \ref apop_ols model, you will find a few lines of toy data and a sample program to run an OLS regression. Once you have set up the \c data file, you can use this Python rewrite of the OLS program:

\include ols.py

\li SWIG sets all C-side global variables into a category named avar. Thus the \c apop_ols variable had to be called \c avar.apop_ols.

\li The verbose package-object-action scheme for naming functions is mostly unnecessary in Python, where objects can more closely be attached to functions. Thus, instead of calling \c apop_vector_skew(my_v), you would call \c my_v.skew(). If you want a list of methods, just use <tt>help(my_v)</tt> or any of the other familiar means of introspection.

Here is another simple example, that copies a Python-side list into a matrix using \c apop_pylist_to_data, and then runs a Fisher Exact test on the matrix. If the list were one-dimensional (flat), then the data would be copied into the vector element of the returned \ref apop_data structure.

\include fisher.py

\li The focus of the work is still in C, so there will likely always be things that you can do in C that can't be done in Python, and strage Python-side errors that will only be explicable if you understand the C-side.  That said, you can still access all of the functions from Python (including those that make little sense from Python).

endofdiv

Outlineheader About SQL, the syntax for querying databases

 For a reference, your best bet is the <a href="http://www.sqlite.org/lang.html">Structured Query Language reference</a> for SQLite.  For a tutorial; there is an abundance of <a href="http://www.google.com/search?q=sql+tutorial">tutorials online</a>.  The blog of Apophenia's author includes an <a href="http://fluff.info/blog/arch/00000118.htm">entry</a> about complementarities between SQL and matrix manipulation packages.

Apophenia currently supports two database engines: SQLite and mySQL. SQLite is the default, because it is simpler and generally more easygoing than mySQL, and supports in-memory databases.

You can switch to mySQL two ways: set <tt>apop_opts.db_engine = 'm'</tt>, or set the environment variable <tt>APOP_DB_ENGINE=mysql</tt>. Otherwise, the system will use SQLite. Ideally, after you make this switch, you need make no other changes--- \ref apop_query, \ref apop_query_to_data, \ref apop_table_exists, et cetera, will work as before. 

Finally, Apophenia provides a few nonstandard SQL functions to facilitate math via database; see \ref db_moments.
endofdiv


Outlineheader mwd The book version

Apophenia co-evolved with <em>Modeling with Data: Tools and Techniques for Statistical Computing</em>. You can read about the book, or download a free PDF copy of the full text, at <a href="http://modelingwithdata.org">modelingwithdata.org</a>.

If you are at this site, there is probably something there for you, including a tutorial on C and general computing form, SQL for data-handing, several chapters of statistics from various perspectives, and more details on working Apophenia. 

As with many computer programs, the preferred manner of citing Apophenia is to cite its related book.
Here is a BibTeX-formatted entry, which should be be easy to re-shape for other environments:

\@book{klemens:modeling,<br>
title = "Modeling with Data: Tools and Techniques for Statistical Computing",<br>
author="Ben Klemens",<br>
year=2008,<br>
publisher="Princeton University Press"<br>
}

endofdiv


endofdiv

Outlineheader dataoverview Data sets

The \ref apop_data structure represents a data set.  It joins together a \c gsl_vector, a \c gsl_matrix, an \ref apop_name, and a table of strings. It tries to be lightweight, so you can use it everywhere you would use a \c gsl_matrix or a \c gsl_vector.

For example, let us say that you are running a regression: there is a vector for the one dependent variable, and a matrix for the several independent variables. Think of them as a partitioned matrix, where the vector is column -1, and the first column of the matrix is column zero. Here is some code to print the entire matrix. Notice that the column counter \c i starts counting at -1.

  \code
  for (j = 0; j< data->matrix->size1; j++){
    printf("%s\t", apop_name_get(data->names, j, 'r'));
    for (i = -1; i< data->matrix->size2; i++)
        printf("%g\t", apop_data_get(data, j, i));
    printf("\n");
    }
    \endcode

We're generally assuming that the data vector and data matrix have the same row count: \c data->vector->size==data->matrix->size1 . This means that the \ref apop_name structure doesn't have separate vector_names and row_names elements: the rownames are assumed to apply for both.

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
    \li\ref apop_matrix_fill()
    \li\ref apop_matrix_increment()
    \li\ref apop_matrix_realloc()
    \li\ref apop_matrix_rm_columns()
    \li\ref apop_matrix_stack()
    \li\ref apop_text_add()
    \li\ref apop_vector_bounded()
    \li\ref apop_vector_copy()
    \li\ref apop_vector_fill()
    \li\ref apop_vector_increment()
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

Outlineheader datalloc Alloc/free

    \li\ref apop_data_alloc()
    \li\ref apop_data_calloc()
    \li\ref apop_data_free()
    \li\ref apop_text_alloc()
    \li\ref apop_text_free()

    See also:

\li <tt>gsl_vector * gsl_vector_alloc (size_t n)</tt>
\li <tt>gsl_vector * gsl_vector_calloc (size_t n)</tt>
\li <tt>void gsl_vector_free (gsl_vector * v)</tt>
\li <tt>gsl_matrix_memcpy (gsl_matrix * dest, const gsl_matrix * src)</tt>

endofdiv

Outlineheader gslviews 	Using GSL Views

You may need to pull a row of a matrix as a distinct vector. The GSL makes it easy to do so, and Apophenia provides several convenience functions to save any remaining hassle. The data in the new vector will point to the base data, not copy it. Some examples:

\code
apop_data *d = apop_query_to_data("select obs1, obs2, obs3 from a_table");
gsl_matrix *m = d->matrix;

//Get a column using its name
Apop_col_t(d, "obs1", ov);
double obs1_sum = apop_vector_sum(ov);

//Get a row using its index
Apop_row(d, 0, v);
double first_row_sum = apop_vector_sum(v);

//First column's sum, matrix format
Apop_matrix_col(m, 0, mv);
double first_col_sum = apop_vector_sum(mv);

//Pull a ten by five submatrix, whose first element is the (2,3)rd
//element of the parent data set's matrix
Apop_submatrix(d, 2,3, 10,5, subm);
double first_col_sum = apop_matrix_sum(subm);
\endcode

\li\ref Apop_col
\li\ref Apop_row
\li\ref Apop_col_t
\li\ref Apop_row_t
\li\ref Apop_matrix_col
\li\ref Apop_matrix_row
\li\ref Apop_submatrix


These macros make use of a set of GSL matrices that produce output of type <tt>gsl_matrix_view</tt> and <tt>gsl_vector_view</tt>. These types point to the source data, and add metadata to turn the data into a coherent matrix/vector. Apophenia's macros generate these views, then pull the matrix/vector from them, so you never have to deal with the <tt>_view</tt> structs directly.

If the macros don't work for you, you can use the GSL's macros directly.  Here's how to get the fifth row of <tt>a_matrix</tt> into a vector view:

\code
gsl_vector_view v;
v = gsl_matrix_col(a_matrix, 4);
\endcode

For rows, use <tt>gsl_matrix_row(a_matrix, n)</tt>. The vector view is a data structure which includes an element of type <tt>gsl_vector</tt> named <tt>vector</tt>; this is the only element you will be interested in. The expression <tt>&(v.vector)</tt> is of type <tt>gsl_vector *</tt>, and therefore can be used as you would any other pointer to a <tt>gsl_vector</tt>. For example, try \ref apop_mean<tt>(&(v.vector));</tt>.

The view is intended to be a common variable, not a pointer. If you want to retain the data after the function exits, copy it to another vector: 
\code
gsl_vector_view v;
gsl_vector *a_new_vector = gsl_vector_alloc(a_matrix->size1);
v = gsl_matrix_col(a_matrix, 4);
gsl_vector_memcpy(a_new_vector, &(v.vector));
\endcode

endofdiv

Outlineheader setgetsec Set/get

Let \f$t\f$ be text, and \f$i\f$ be an index. Then you can get or set an element by row title and column index, row index and column title, et cetera. The functions that have no suffix find an element by row and column index.

The \c apop_data_ptr_... functions return a pointer to the element.

        \li\ref apop_data_get()
        \li\ref apop_data_get_it()
        \li\ref apop_data_get_ti()
        \li\ref apop_data_get_tt()
        \li\ref apop_data_ptr()
        \li\ref apop_data_ptr_it()
        \li\ref apop_data_ptr_ti()
        \li\ref apop_data_ptr_tt()
        \li\ref apop_data_set()
        \li\ref apop_data_set_it()
        \li\ref apop_data_set_ti()
        \li\ref apop_data_set_tt()

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

\anchor outline_mapply Map/apply
These functions allow you to send each element of a vector or matrix to a function, either producing a new matrix (map) or transforming the original (apply).  The \c ..._sum functions return the sum of the mapped output.

There is an older and a newer set of functions. The older versions, which act on <tt>gsl_matrix</tt>es or <tt>gsl_vector</tt>s have more verbose names; the newer versions, which act on the elements of an \ref apop_data set, use the \ref designated syntax to ad a few options and a more brief syntax.

You can do many things quickly with these functions.

Get the sum of squares of a vector's elements:

\code
double sum_of_squares = apop_map_sum(dataset, gsl_pow_2); //given <tt> apop_data *dataset</tt>

double sum_of_squares = apop_vector_map_sum(v, gsl_pow_2); //given <tt> gsl_vector *v</tt>
\endcode

Here, we create an index vector [0, 1, 2, \f$\dots\f$].

\code
double index(double in, int index){return index;}
apop_data *d = apop_data_alloc(100, 0, 0);
apop_map(d, .fn_di=index, .inplace=1);
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
static double nan_check(const double in){ return isnan(in);}

static double apop_mean_no_nans(apop_data *in){
    return apop_map_sum(in, no_nan_val)/apop_map_sum(in, nan_check);
}
\endcode

This program randomly generates a data set where each row is a list of numbers with a different mean. It then finds the \f$t\f$ statistic for each row, and the confidence with which we reject the claim that the statistic is less than or equal to zero.

Notice how the older \ref apop_vector_apply uses file-global variables to pass information into the functions, while the \ref apop_map uses a pointer to the constant parameters to input to the functions.

\include t_test_by_rows.c

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
        \li\ref apop_vector_distance : find the Euclidian distance between two vectors
        \li\ref apop_vector_grid_distance : find the distance via the Manhattan metric between two vectors
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

        \li\ref apop_det_and_inv : find determinant and inverse at the same time
        \li\ref apop_dot : matrix \f$\cdot\f$ matrix, matrix \f$\cdot\f$ vector, or vector \f$\cdot\f$ matrix
        \li\ref apop_matrix_determinant : returns the determinant of the input matrix
        \li\ref apop_matrix_inverse : returns the inverse of the input matrix

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

            \li\ref apop_data_correlation ()
            \li\ref apop_data_covariance ()
            \li\ref apop_data_summarize ()
            \li\ref apop_matrix_correlation ()
            \li\ref apop_matrix_covariance ()
            \li\ref apop_matrix_mean ()
            \li\ref apop_matrix_mean_and_var ()
            \li\ref apop_matrix_sum ()
            \li\ref apop_matrix_var_m ()
            \li\ref apop_mean()
            \li\ref apop_sum()
            \li\ref apop_var()
            \li\ref apop_vector_correlation ()
            \li\ref apop_vector_cov ()
            \li\ref apop_vector_kurt ()
            \li\ref apop_vector_kurtosis ()
            \li\ref apop_vector_kurtosis_pop ()
            \li\ref apop_vector_mean()
            \li\ref apop_vector_skew()
            \li\ref apop_vector_skew_pop()
            \li\ref apop_vector_sum()
            \li\ref apop_vector_var()
            \li\ref apop_vector_var_m ()
            \li\ref apop_vector_weighted_cov ()
            \li\ref apop_vector_weighted_kurt ()
            \li\ref apop_vector_weighted_mean ()
            \li\ref apop_vector_weighted_skew ()
            \li\ref apop_vector_weighted_var ()

endofdiv

Outlineheader convsec   Conversion among types

            \li\ref apop_array_to_data()
            \li\ref apop_array_to_matrix()
            \li\ref apop_array_to_vector()
            \li\ref apop_line_to_data()
            \li\ref apop_line_to_matrix()
            \li\ref apop_matrix_to_data()
            \li\ref apop_matrix_to_db()
            \li\ref apop_text_to_data()
            \li\ref apop_text_to_db()
            \li\ref apop_vector_to_array()
            \li\ref apop_vector_to_data()
            \li\ref apop_vector_to_matrix()

endofdiv

Outlineheader names   Name handling

                \li\ref apop_name_add()
            \li\ref apop_name_alloc()
            \li\ref apop_name_copy()
            \li\ref apop_name_find()
            \li\ref apop_name_free()
            \li\ref apop_name_print()
            \li\ref apop_name_stack()

endofdiv

Outlineheader fact   Generating factors

\em Factor is jargon for a numbered category. Number-crunching programs work best on numbers, so we need a function to produce a one-to-one mapping from text categories into numeric factors. 

A \em dummy is a variable that is either one or zero, depending on membership in a given group. Some methods (typically when the variable is an input or independent variable) prefer dummies; some methods (typically for outcome or dependent variables) pefer factors.

            \li\ref apop_data_to_dummies()
            \li\ref apop_text_to_factors()
            \li\ref apop_text_unique_elements()
            \li\ref apop_vector_unique_elements()

endofdiv

endofdiv

Outlineheader dbs Databases

These are convenience functions to handle interaction with SQLite or mySQL. They open one and only one database, and handle most of the interaction therewith for you.

You will probably first use \ref apop_text_to_db to pull data into the database, then \ref apop_query to clean the data in the database, and finally \ref apop_query_to_data to pull some subset of the data out for analysis.

See also the \ref conversions, including \ref apop_text_to_db and \ref apop_matrix_to_db.

\li \ref apop_query : Manipulate the database, return nothing (e.g., input data).
\li \ref apop_db_open : Optional, for when you want to use a database on disk.
\li \ref apop_db_close : If you used \ref apop_db_open, you will need to use this too.
\li \ref apop_table_exists : Check to make sure you aren't reinventing or destroying data. Also, the clean way to drop a table.
\li \ref apop_db_merge : Import or merge the whole of another database into the currently open db.
\li \ref apop_db_merge_table : Import/merge just one table.
   \li \ref apop_crosstab_to_db : Convert between two common data layouts
    \li \ref apop_db_rng_init : Apophenia maintains a database-side RNG; set its seed it with this.

\par P.S.
Apophenia reserves the right to insert temp tables into the opened database. They will all have names beginning with "apop_", so the reader is advised to not use tables with such names, and is free to ignore or delete any such tables that turn up.

    Outlineheader dbout Out

        \li\ref apop_db_to_crosstab()
        \li\ref apop_query_to_data()
        \li\ref apop_query_to_float()
        \li\ref apop_query_to_matrix()
        \li\ref apop_query_to_mixed_data()
        \li\ref apop_query_to_text()
        \li\ref apop_query_to_vector()

endofdiv

    Outlineheader dbin In

        See the print functions below. By setting <tt>apop_opts.output_type = 'd'</tt>, \ref apop_data sets (or \c gsl_matrixes and \c gsl_vectors) are `printed' to the database.

endofdiv

    Outlineheader dbmath Math in the db  \anchor dbttest
    
        \li\ref apop_db_paired_t_test()
        \li\ref apop_db_t_test()

endofdiv

endofdiv

Outlineheader Modesec Models \anchor models

Outlineheader introtomodels Introduction

A model intermediates between data, parameters, and likelihoods.  [This definition readily includes "non-parametric" models.] When the parameters are fully specified, the model will give you a likelihood that an observation will occur.  Much of statistical analysis consists of writing down a model, estimating its parameters, and running hypothesis tests to determine the confidence with which we can make statements about those parameters.

Given data, the model can estimate the parameters that best fit the data. This is the typical model estimation.

Given parameters, the model can generate the expected value of the data, or randomly generate a full artificial data set.

Apophenia facilitates all this via its \ref apop_model objects. The model includes slots for \c log_likelihood, \c estimate, \c draw, \c expected_value, and other operations we expect from our models.

For example, a model may be a probability distribution, such as the \ref apop_normal, \ref apop_poisson, or \ref apop_beta models. The data is assumed to have been drawn from a given distribution and the question is only what distributional parameters best fit; e.g., assume the data is Normally distributed and find the mean and variance: \c apop_estimate(data,apop_normal).

The design of the objects hopes to make it as easy as possible for you, dear reader, to write new models. For the most part, all you need to do is write a log likelihood function, and \ref apop_maximum_likelihood does the rest; see below.

The most commonly-used function in the systems below is the \ref apop_estimate function. It takes in a model and data, and outputs an apop_estimate, that includes the parameter estimates and the various auxiliary data that one may need to test the estimates, such as the variance-covariance matrix. For most users, the \c estimate function will be all one needs from a model. Just prep the data, select a model, and produce an estimate:

\code
apop_data 	    *data 		    = read_in_data();
apop_model 	*the_estimate 	= apop_estimate(data, apop_probit);
apop_model_show(the_estimate);
\endcode

Outlineheader internals The internals

\image html http://apophenia.sourceforge.net/doc/model.png
\image latex model.png

The \ref apop_model struct breaks down into three parts:

\li Info like names, the pointer to the input data, and the parameters, are for the most part self-descriptive.

\li There is a broad class of functions that cover most of what you would do with a model. You can see that there is a bit of a bias toward maximum likelihood estimation. There are helper functions for most of them, like the \ref apop_estimate function above, that meant that you never had to directly handle the model's \c estimate method.  In the list at the top of this page, all of the functions whose last argument is a model are helper functions of this type. The helper functions do some boilerplate error checking, and mean that you don't have to fill in every blank in your model: if you have a \c log_likelihood method but no \c p method, then \ref apop_p will use exp(\c log_likelihood). If you don't give an \c estimate method, then \c apop_estimate will use maximum likelihood.

\li I refer to the values output from an estimation as \em parameters, and the details of how the model's machinery work as \em settings. The parameters are pretty standardized, and a \ref apop_data set is sufficient to handle the great majority of models. However, the settings of an OLS regression are drastically different from those of a histogram or an MLE.  The solution is a list of settings structures. This is probably the least pleasant structural element in the system, to the point that you will not want to handle the settings structures directly. Instead, you should use the various \ref settings "helper functions for model settings".

endofdiv

Outlineheader write_likelihoods Writing your own 

Writing \c apop_model objects is easy, because (almost) everything has a default of some sort, so you can fill in those segments you know, and let the system use computationally-intensive methods for the rest.

For example, here is how one would set up a model for an MLE:

\li Write a likelihood function. Its header will look like this:
\code
double apop_new_log_likelihood(apop_data *data, apop_model *m)
\endcode 
where \c data is the input data, and \c
m is the parametrized model (i.e., your model with a \c parameters element set by the caller). 
This function will return the value of the log likelihood function at the given parameters.

\li Is this a constrained optimization? See the \ref constraints "Constraints page" on how to set them. Otherwise, no constraints will be assumed.

\li Write the object. In your header file, include 
\code
apop_model your_new_model = {"The Me distribution", number_of_parameters, 
            .estimate = new_estimate, .log_likelihood = new_log_likelihood };
\endcode
If there are constraints, add an element for those too.

\c number_of_parameters is probably a positive integer like \c 2, but
it is often (the number of columns in your data set) -1, in which case,
set \c number_of_parameters to \c -1.

You already have enough that something like
\code
apop_model *estimated = apop_mle(your_data, your_new_model);
\endcode
will work. 

For many methods, the first thing you will write is the random number generator. Other elements, such as the 
gradient---
\code
void apop_new_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
    //some algebra here to find df/dp0, df/dp1, df/dp2....
    gsl_vector_set(gradient,0, d_0);
    gsl_vector_set(gradient,1, d_1);
}
\endcode 
---or the expected value are optional. 

endofdiv

endofdiv

    Outlineheader modellist  Models that ship with Apophenia

        Outlineheader Dist Distributions

            \li\ref apop_bernoulli
            \li\ref apop_beta
            \li\ref apop_beta_from_mean_var()
            \li\ref apop_binomial
            \li\ref apop_chi_squared
            \li\ref apop_exponential
            \li\ref apop_f_distribution
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
            \li\ref apop_waring
            \li\ref apop_wishart
            \li\ref apop_yule
            \li\ref apop_zipf

        endofdiv

        Outlineheader GLM  GLM family

            \li\ref apop_iv
            \li\ref apop_logit
            \li\ref apop_multinomial_probit
            \li\ref apop_ols
            \li\ref apop_probit
            \li\ref apop_wls

        endofdiv

        Outlineheader moremodels More

            \li\ref apop_histogram  (see the histogram section below)
            \li\ref apop_kernel_density

        endofdiv

    endofdiv

    Outlineheader Basicmodelmethods Basic object methods

        \li\ref apop_model_clear()
        \li\ref apop_model_copy()
        \li\ref apop_model_free()
        \li\ref apop_model_set_parameters()
        \li\ref apop_model_show()
        \li\ref apop_model_prep()

    endofdiv

    Outlineheader mathmethods Methods for computation

        \li\ref apop_draw()
        \li\ref apop_estimate()
        \li\ref apop_expected_value()
        \li\ref apop_log_likelihood()
        \li\ref apop_score()
        \li\ref apop_model_print()
        \li\ref apop_p()

    endofdiv

    Outlineheader modelsettings Model settings \anchor settings

Apophenia is really only based on two objects, the \ref apop_data set and the \ref apop_model. Data sets come in a pretty standard form, so the data object is basically settled. But describing a statistical, agent-based, social, or physical model in a reasonably standardized form is much more difficult, primarily because every model has significantly different settings. E.g., an MLE requires a method of search (conjugate gradient, simplex, simulated annealing), and a histogram needs the number of slots to be filled with data.

So, the \ref apop_model includes a single list, whose name is simply \c settings, which can hold an arbitrary number of groups of settings. For example, you can have a set of closed-form variables for estimating the log likelihood, and a histogram for making random draws.

To get/set a setting, you would need to specify the model, the settings group, and the name of the setting itself.  For the sake of giving you a mental model, this is much like a g_lib config file (or for Windows users, a .ini file), which would have settings divided into sections, like

\code
[MLE]
verbose = 0
method = APOP_CG_PR

[least squares]
weights = your_weight_vector

[histogram]
bins = {0.1, 0.3, 0.2,...}
\endcode


\b Using it

If you don't need to change settings from the default, you don't need to care about any of this, because
\code
apop_model_show(apop_estimate(your_data_here, apop_ols));
\endcode
---still works fine.

If you do need to change settings, then the process takes another step or two. Here's a sample:

\code
1 apop_data *data = your_data_here;
2 apop_data *w = your_weights_here;

3 apop_model *m = apop_model_copy(apop_wls);
4 Apop_model_add_group(m, apop_ls, .weights=w);
5 apop_model *est = apop_estimate(data, *m);
\endcode

Line three establishes the baseline form of the model. Line four adds a settings group of type \ref apop_ls_settings to the model, and specifies that we want it initialized with the \c weights element set to the data set \c w. The list of elements of the settings structures included in Apophenia are below; the means of specifying the elements follow the \ref designated syntax.

Having set the settings, line 5 does the weighted OLS.

Also, the output to any of Apophenia's estimations will have an appropriate group of settings allocated, so you can chain estimations pretty easily. Continuing the above example, you could re-estimate with an alternative set of weights via:

\code
Apop_settings_add(est, apop_ls, weights, weight_set_two);
apop_model *est2 = apop_estimate(data, *est);
\endcode

This uses the above-mentioned address of model/settings group/item, plus the actual value to be set.  Notice that the \c Apop_settings_add macro really just modifies the existing structure; you can think of the defaults as a blank slate, and the only values being those you'd set.

If you need to read a setting, such as to check whether things have changed after an estimation, use, e.g.:

\code
apop_data *weights_now = Apop_settings_get(m, apop_ls, weights);
\endcode

\li Notice the use of a single capital to remind you that you are using a macro, so you should beware of the sort of surprising errors associated with macros. Here in the modern day, we read things like APOP_SETTINGS_ADD as yelling, but if you prefer all caps to indicate macros, those work as well.

\li There are two additional macros which are now deprecated: \ref
Apop_settings_add_group and \c Apop_settings_alloc_add. They made sense at the time.

For just using a model, that's about 100% of what you need to know.

Outlineheader settingswritng  Writing new settings

To store the settings for your own models, you don't necessarily need any of this. The \ref apop_model structure has a \c void pointer named \c more which you can use as you see fit. If \c more_size is larger than zero (i.e., you set it to \c sizeof(your_struct)), then it will be copied via \c memcpy as necessary. Apohenia's estimation routines will never impinge on this item, so do what you feel with it.

If you do want to set up a new model, then you will need four items.  This is the sort of boilerplate that will be familiar to users of object oriented languages in the style of C++ or Java. Let your settings group be named \c ysg; then you will need

\li The settings struct
\code
typedef struct {
...
} ysg_settings
\endcode


\li The allocate function
\code
ysg_settings *ysg_settings_init(ysg_settings in){ 
    ysg_settings *out = malloc(sizeof(ysg_settings));
    // Give default values for all elements here. 
    // There is a macro, e.g.,
    // apop_varad_setting(in, out, want_cov,  'y');
    // that will help you with this. It checks in
    // for the given element (here, want_cov), and 
    // if that element is not found, sets it to the 
    // specified value (here, 1).
    return out; }
\endcode

\li The copy function
\code
void *ysg_settings_copy(ysg_settings *copyme) {
    ysg_settings *out = malloc(sizeof(ysg_settings));
    //copy elements here. If you have no pointers or other trickery, just use:
    *out = *copyme
    return out; }
\endcode

\li The free function, which can be as brief as:
\code
void ysg_settings_free(ysg_settings *freeme) {
    free(freeme);
}
\endcode
but may include a freeing of pointed-to subelements as necessary.

The names are not negotiable: when you call
\code
Apop_model_add_group(m, ysg)
\endcode
the macro will look for \c ysg_settings, \c ysg_settings_init, et cetera.

The lines-of-code averse will cringe at having to write such boilerplate code (I do), but after spending a year resisting it, I have to concede that it's the least of all evils.

You can retrieve the whole group or individual elements via:
\code
Apop_settings_get_group(m, ysg)
//as well as the one-element macro mentioned above,
Apop_settings_get(m, ysg, an_element)
\endcode

As you saw above, once the typedef/alloc/copy/free machinery is written, you can declare, get, and set in a reasonably graceful manner.

\li If for some reason you want to work with the now-deprecated \ref Apop_settings_add_group, you'll need a \c ysg_settings_alloc(...) function instead of \c ysg_settings_init. They do the same thing, but the \c alloc version takes in a variadic list of items, wile the \c init takes a single structure.

endofdiv

        \li\ref Apop_settings_set
        \li\ref Apop_settings_add_group
        \li\ref Apop_settings_alloc
        \li\ref apop_settings_copy_group
        \li\ref Apop_settings_get
        \li\ref apop_settings_get_group
        \li\ref Apop_settings_get_group
        \li\ref apop_settings_group_alloc

        Outlineheader Spec Specific settings

            \li\ref apop_category_settings_init()  
            \li\ref apop_histogram_settings_init()  
            \li\ref apop_ls_settings_init()
            \li\ref apop_mle_settings_init() \anchor mlesettingsa
            \li\ref apop_rank_settings_init()
            \li\ref apop_update_settings_init()

        endofdiv

    endofdiv

    Outlineheader Esti Estimation aids

        \li\ref apop_model_fix_params : hold some parameters constant

    endofdiv

endofdiv

Outlineheader Test Tests & diagnostics

   Just about any hypothesis test consists of a few common steps:
 
\li  specify a statistic
\li  State the statistic's (hypothesized) distribution
\li  Find the odds that the statistic would lie within some given range, like `greater than zero' or `near 1.1'

If the statistic is from a common form, like the parameters from an OLS regression, then the commonly-associated \f$t\f$ test is probably thrown in.

Some tests, like ANOVA, produce a statistic using a specialized prodecure, so Apophenia includes some functions, like \ref apop_test_anova_independence and \ref apop_test_kolmogorov, to produce the statistic and look up its significance level.

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
    \li\ref apop_estimate_parameter_t_tests()
    
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
  apop_assert(in, 0, 0, 'c', "input vector is NULL. Doing nothing.\n");
  gsl_vector	*normed;
  double 		sum=0;
	apop_vector_normalize((gsl_vector *)in,&normed, 1);
	gsl_vector_mul(normed,normed);
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

Outlineheader Histosec Histograms \anchor histograms

The GSL provides a <tt>gsl_histogram</tt> structure, that produces a PDF by accumulating data into bins. Apophenia wraps this into its \ref apop_model struct, so that the model-family machinery can be applied to Bayesian updating, kernel smoothing, or other methods that output a PDF.

To produce a PMF from \c your_data, use \ref apop_histogram "apop_estimate(your_data, apop_histogram)". If you would like to compare this histogram to other data (observed or theoretical),then you'll need to produce a second synce histogram using \ref apop_histogram_vector_reset or \ref apop_histogram_model_reset.  Then you can send both histograms to, say, \ref apop_test_kolmogorov.

The second structure from the GSL incrementally sums up the PMF's bins to produce a CMF. The CMF can be used to map from a draw from a Uniform[0,1] to a draw from the PMF.  Because it can be used to draw from the PMF, the GSL calls this the <tt>gsl_histogram_pdf</tt> structure. That's right: the the data in the <tt>gsl_histogram_pdf</tt> structure is a cumulative sum---a CMF.

Anyway, here are some functions to deal with these various histograms and such; see also the GSL documentation, linked at the bottom of this outline.

    \li\ref apop_histogram_model_reset()
    \li\ref apop_histogram_moving_average()
    \li\ref apop_histogram_normalize()
    \li\ref apop_histogram_vector_reset()
    \li\ref apop_histograms_test_goodness_of_fit()

endofdiv

Outlineheader Maxi Maximum likelihood methods

If you read the section on writing models, then you already know how to do maxium likelihood on exotic setups. Just write a model that has a \c p or \c log_likelihood function, and call \c apop_estimate(your_data,your_model).  The default estimation routine is maximum likelihood.

MLEs have an especially large number of parameter tweaks that could be made; see the section on MLE settings above.

Outlineheader constr Setting Constraints \anchor constraints

The problem is that the parameters of a function must not take on certain values, either because the function is undefined for those values or because parameters with certain values would not fit the real-world problem.

The solution is to rewrite the function being maximized such that the function is continuous at the constraint boundary but takes a steep downward slope. The unconstrained maximization routines will be able to search a continuous function but will never return a solution that falls beyond the parameter limits.

If you give it a likelihood function with no regard to constraints plus an array of constraints, 
\ref apop_maximum_likelihood will combine them to a function that fits the above description and search accordingly.

A constraint function must do three things:
\li It must check the constraint, and if the constraint does not bind (i.e., the parameter values are OK), then it must return zero.
\li If the constraint does bind, it must return a penalty, that indicates how far off the parameter is from meeting the constraint.
\li if the constraint does bind, it must set a return vector that the likelihood function can take as a valid input. The penalty at this returned value must be zero.

The idea is that if the constraint returns zero, the log likelihood function will return the log likelihood as usual, and if not, it will return the log likelihood at the constraint's return vector minus the penalty. To give a concrete example, here is a constraint function that will ensure that both parameters of a two-dimensional input are both greater than zero:

\code
static double beta_zero_greater_than_x_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_2
  static apop_data *constraint = NULL;
    if (!constraint){
        constraint= apop_data_calloc(1,1,1);
        apop_data_set(constraint, 0, 0, 1);
    }
    return apop_linear_constraint(v->parameters->vector, constraint, 1e-3);
}
\endcode

    \li\ref apop_linear_constraint 

endofdiv

    \li\ref apop_estimate_restart : Restarting an MLE with different settings can improve resuls.
    \li\ref apop_maximum_likelihood()
    \li\ref apop_model_numerical_covariance()
    \li\ref apop_numerical_gradient()

endofdiv

Outlineheader Asst More descriptive methods

    \li\ref apop_matrix_pca : Principal component analysis
    \li\ref apop_anova : One-way or two-way ANOVA tables
    \li\ref apop_update : Bayesian updating

endofdiv

Outlineheader Miss Missing data

    \li\ref apop_data_listwise_delete()
    \li\ref apop_ml_imputation()

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
apop_matrix_print(m, .outpipe = stdout);  //now show the same matrix on screen
\endcode

I will first look to the input file name, then the input pipe, then the
global \c output_pipe, in that order, to determine to where I should
write.  Some combinations (like output type = \c 'd' and only a pipe) don't
make sense, and I'll try to warn you about those. 

The plot functions produce output for Gnuplot (so output type = \c 'd'
again does not make sense). As above, you can pipe directly to Gnuplot or write to a file.

    \li\ref apop_plot_histogram()
    \li\ref apop_plot_line_and_scatter()
    \li\ref apop_plot_lattice()
    \li\ref apop_data_print()  
    \li\ref apop_data_show()
    \li\ref apop_histogram_print()
    \li\ref apop_matrix_print()
    \li\ref apop_matrix_show()
    \li\ref apop_vector_print()
    \li\ref apop_vector_show()

endofdiv

Outlineheader moreasst Assorted

Outlineheader Gene General utilities

    \li\ref Apop_assert
    \li\ref Apop_assert_void
    \li\ref apop_error()
    \li\ref apop_opts
    \li\ref apop_strip_dots()
    \li\ref apop_system()

endofdiv

Outlineheader Math Math utilities

    \li\ref apop_matrix_is_positive_semidefinite()
    \li\ref apop_matrix_to_positive_semidefinite()
    \li\ref apop_generalized_harmonic()
    \li\ref apop_multivariate_gamma()
    \li\ref apop_multivariate_lngamma()
    \li\ref apop_rng_alloc()
    \li\ref apop_random_double ()
    \li\ref apop_random_int ()

endofdiv

Outlineheader Prob Deprecated

These functions will probably disappear or be replaced soon.

    \li\ref apop_estimate_fixed_effects_OLS()
    \li\ref apop_name_cross_stack()
    \li\ref apop_kernel_density_settings_alloc()
    \li\ref apop_data_to_db()

endofdiv

endofdiv

Outlineheader links Further references

For your convenience, here are links to some other libraries you are probably using.

    \li <a href="http://www.gnu.org/software/libc/manual/html_node/index.html">The standard C library</a>
    \li <a href="c_precedence.html">The C operator precedence table</a>
    \li <a href="http://www.gnu.org/software/gsl/manual/html_node/index.html">The
    GSL documentation</a>, and <a href="http://www.gnu.org/software/gsl/manual/html_node/Function-Index.html">its index</a>
    \li <a href="http://sqlite.org/lang.html">SQL understood by SQLite</a>

*/


/** \page mingw MinGW

Minimalist GNU for Windows is indeed mimimalist: it is not a full POSIX subsystem, and provides no package manager. Therefore, you will have to make some adjustments and install the dependencies yourself.

Matt P. Dziubinski successfully used Apophenia via MinGW; here are his instructions (with edits by BK):

\li For the Python interface, get Python Windows installer from: http://python.org/download/

\li If installed to "c:\bin\prog\Python2", for example, you'll need to set:
\code
export PYTHON=/C/bin/prog/Python2/python
\endcode

\li Get SWIG at: http://www.swig.org/download.html . Being on Windows, I've opted for "swigwin".

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

\li Compile the Python interface:
\code
make apop.py
\endcode
should now run OK. Ignore warnings that
<em>'print' is a python keyword...</em>.

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
    apop_assert_void(v, 0, 's', "You sent me a NULL vector.");
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
Thus, the macro declares each not-in-struct variable, and so there will need to be one such declaration line for each argument. Apart from requiring delcarations, you can be creative: include sanity checks, post-vary the variables of the inputs, unpack without the macro, and so on. That is, this parent function does all of the bookkeeping, checking, and introductory shunting, so the base function can just do the math. Finally, the introductory section has to call the base function.

The setup goes out of its way to leave the \c _base function in the public namespace, so that those who would prefer speed to bounds-checking can simply call that function directly, using standard notation. You could eliminate this feature by just merging the two functions.


<b>The sed script</b>

The above is all you need to make this work: the varad.h file, and the above structures. But there is still a lot of redundancy, which can't be eliminated by the  plain C preprocessor.

Thus, below is another preprocessor that converts a few markers to the above form. Here is the code that will expand to the above C-standard code:

\code
//header file
APOP_VAR_DECLARE void apop_vector_increment(gsl_vector * v, int i, double amt);

//code file
APOP_VAR_HEAD void apop_vector_increment(gsl_vector * v, int i, double amt){
    gsl_vector * apop_varad_var(v, NULL);
    apop_assert_void(v, 0, 's', "You sent me a NULL vector.");
    int apop_varad_var(i, 0);
    double apop_varad_var(amt, 1);
    apop_vector_increment_base(v, i, amt);
APOP_VAR_END_HEAD
	v->data[i * v->stride]	+= amt;
}
\endcode

It is obviously much shorter. The declaration line is actually a C-standard declaration with the \c APOP_VAR_DECLARE preface, so you don't have to remember when to use semicolons. The function itself looks like a single function, but there is again a marker before the declaration line, and the introductory material is separated from the main matter by the \c APOP_VAR_END_HEAD line. Done right, drawing a line between the introductory checks or initializations and the main function can really improve readability.

One final detail: it is valid to have types with commas in them---function arguments. Because commas get turned to semicolons, and sed isn't a real parser, there is an exception built in: you will have to replace commas with exclamation marks in the header file (only). E.g.,

\code
APOP_VAR_DECLARE apop_data * f_of_f(apop_data *in, void *param, int n, double (*fn_d)(double ! void * !int));
\endcode

Sed is POSIX standard, so even if you can't read the script (included in the subversion directory), you have the program needed to run it. For example, if you name it \c prep_variadics.sed, then run
\code
./prep_variadics.sed < myfile.pre.c > myfile.c
\endcode

*/

 

/* I should do something with this:
 *
 *
 *
 *
 * This function can be used to temporarily modify the global options,
 to facilitate better encapsulation of code. Usage:

  \code
  apop_opts_type tmp_opts;
  apop_opts_memcpy(&tmp_opts, &apop_opts, sizeof(apop_opts_type));
  strcpy(apop_opts.output_name, "ad_hoc_temp_file");
  [do things here]
  apop_opts_memcpy(&apop_opts, &tmp_opts, sizeof(apop_opts_type));
  \endcode

If you just need a little more verbosity for a procedure, you probably
don't need to use this function. Just try: 
  \code
  apop_opts.verbose ++;
  [do things here]
  apop_opts.verbose --;
  \endcode

The philosophy is that the global variables are generally not going
to change over the course of a program: either you are working on the
screen, in the database, or piping out of STDOUT, and you likely won't
change mid-stream. Thus, it is easier to set these globally at the top of
the program but less convenient to switch frequently throughout the code.
\ingroup global_vars
 */


/** \page dataprep Data prep rules 

There are a lot of ways your data can come in, and we would like to run estimations on a reasonably standardized form.

First, this page will give a little rationale, which you are welcome to skip, and then will present the set of rules.

  \section Dealing with the ones column
Most standard regression-type estimations require or generally expect
a constant column. That is, the 0th column of your data is a constant (one), so the first parameter
\f$\beta_1\f$ is slightly special in corresponding to a constant rather than a variable.

However, there are some estimations that do not use the constant column.

\em "Why not implicitly assume the ones column?"
Some stats packages implicitly assume a constant column, which the user never sees. This violates the principle of transparency
upon which Apophenia is based, and is generally annoying.  Given a data matrix \f$X\f$ with the estimated parameters \f$\beta\f$, 
if the model asserts that the product \f$X\beta\f$ has meaning, then you should be able to calculate that product. With a ones column, a dot product is one line \c apop_dot(x, your_est->parameters, 0, 0)); without a ones column, the problem is left as an unpleasant exercise for the reader.

  \subsection Shunting columns around.

Each regression-type estimation has one dependent variable and several independent. In the end, we want the dependent variable to be in the vector element. However, continuing the \em "lassies faire" tradition, doing major surgery on the data, such as removing a column and moving in all subsequent columns, is more invasive than an estimation should be.

\subsection The rules
So those are the two main considerations in prepping data. Here are the rules, intended to balance those considerations:

\subsubsection The automatic case
There is one clever trick we can use to resolve both the need for a ones column and for having the dependent column in the vector: given a data set with no vector element and the dependent variable in the first column of the matrix, we can copy the dependent variable into the vector and then replace the first column of the matrix with ones. The result
fits all of the above expectations.

You as a user merely have to send in a \c apop_data set with no vector and a dependent column in the first column.

\subsubsection The already-prepped case
If your data has a vector element, then the prep routines won't try to force something to be there. That is, they won't move anything, and won't turn anything into a constant column. If you don't want to use a constant column, or your data has already been prepped by an estimation, then this is what you want.

You as a user just have to send in a \c apop_data set with a filled vector element.
 */

/* To finish and add to the dataprep section:
\paragraph Probit, logit, and other -obits
The dependent variable for these models is a list of categories; this option is not relevant to continuous-valued dependent variables. 

Your data source may include a list of numeric categories, in which case you can pick one of the above cases.

The main exception is when your data is a list of text factors, in which case your dependent variable isn't even a part of the data matrix.  In this case, you can prep the data yourself, via a call to \c apop_text_to_factors, and then insert the column yourself (and pick a constant column or not as you prefer). Or, you can have the system do it for you, via a form like \code
int textcol = 3; //where is the list of dependent categories
apop_model *setmodel = apop_model_copy(apop_probit);
Apop_settings_add_group(setmodel, apop_category, textcol);
apop_estimate(yourdata, setmodel);
\endcode


You'll see that there are two questions here: should there be a constant column of ones, and where is the dependent column to be found?

Here are the rules for preparing the data set. 

The first item is the 


There are two methods:

\li If the data set has no vector, 


for the -obit and -ogit models

   If there is a vector in place, then I won't touch anything.

   If there is no vector in place, then:
        --If you don't tell me where to find the dependent column, I'll go with column zero, and
            --move the data to the vector
            --replace the data there with ones, creating a constant column.
        --If you do tell me where to find the dependent column, via the settings, I'll turn that into a list of factors.

 */
