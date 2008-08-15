/* Apophenia's documentation
Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

/** \mainpage The Apophenia documentation

Apophenia is a library of a few hundred functions (in C) intended to make life easier when doing statistical computing.

This welcome page provides some introductory information to give you an
overview. But most of the documentation is folded into the \ref outline "outline",
linked above. This <a href="http://modelingwithdata.org">textbook</a>
on scientific computing also includes coverage of Apophenia.

 \li \ref intro "Intro": The motivation for the package.
 \li \ref setup "Setup": Installing GCC, the GSL, SQLite, and Apophenia itself.
 \li \ref c "C": Some notes on C and Apophenia's use of C utilities.
 \li \ref sql About the syntax for querying databases.

 */

/* Old, no longer used:
\section dostats Doing statistics 
 \li \ref basic_stats "Basic statistics": Mean, variance, percentiles, &c.
 \li \ref regression  
 \li \ref mle "Maximum likelihood estimation": estimators requiring a search for the maximum of a likelihood function.
 \li \ref linear_algebra "Linear Algebra": determinants, projections, numerical gradients, &c. Some convenience functions to display matrices and vectors.
 \li \ref ttest "Some tests"
 \li \ref histograms 
 \li \ref boot

\section shuntdata Shunting data 
 \li \ref conversions
 \li \ref output "Output functions": Summarize data and print tables to the screen or a file.
 \li \ref convenience_fns "Convenience functions": a few utilities to make life with the GSL a little easier.
 */

/** \page intro Intro


The key goal of Apophenia is to estimate models using data. As such, Apophenia provides 
two interlocking structures to smooth the process: the \ref apop_data and \ref apop_model. 

On the data side, the intent is to provide 
the usual tools for manipulating data, such as getting it in and out of a text file or database, dealing with metadata like column labels, stacking it, splitting it, and otherwise collating it. 

Every \ref
apop_model includes an \c estimate function, that takes in a data set
and outputs a vector (or matrix) of parameters. Notice that this broad description
includes "non-parameteric" methods, the process of fitting a distribution
to a data set, and about anything else that a statistician could want
to do.

Thus, the typical analysis using Apophenia would take the following steps:
 \li Read the data into the database using \ref apop_text_to_db.
 \li Use SQL queries handled by \ref apop_query to massage the data as needed.
 \li Use \ref apop_query_to_data to pull the data into an in-memory apop_data set.
 \li Call a model estimation such as \code apop_estimate (data_set, apop_OLS)\endcode  or \code apop_estimate (data_set, apop_probit)\endcode to fit parameters to the data. This will return an \ref apop_model like the original, but with parameter estimates.
 \li Interrogate the returned estimate, by dumping it to the screen with \ref apop_model_show, sending its parameters and variance-covariance matrices to a test, et cetera.

If this seems a bit vague, have a look at this \ref sample_program.

Of course, you may not need all this: you might just want a neat interface to SQLite, or an easy-to-use structure that holds your data 
in a matrix and allows you to label the rows and columns. Have a look at the \ref outline "library outline" for the various categories of relatively standalone functions that you can fold into your own work.

*/

/** \page setup Setting up
\section cast The supporting cast 
To use Apophenia, you will need to have a working C compiler, the GSL (v1.7 or higher) and SQLite installed. 

We've moved the setup documentation to <a href="http://modelingwithdata.org/appendix_o.html">Appendix O</a> of <em> Modeling with Data</em>. Please see that page.

\subsection testing Testing
There is a short, complete program in the \ref apop_ols entry which runs a simple OLS regression on a data file. Follow
the instructions there to compile and run. See also the 
\ref sample_program below.

\subsection using Using 
Now that it's installed, do some stats! For those who are new to C, there are some notes on the \ref c "C" page; if you are familiar with the process, see the \ref usagenotes "usage notes" for Apophenia itself.

\subsection sample_program sample program
The sample program below is intended to show how one would integrate
Apophenia into an existing program. For example, say that you are running
a simulation of two different treatments, or say that two sensors are
posting data at regular intervals. You need to gather the data in an
organized form, and then ask questions of the resulting data set.
Below, a thousand draws are made from the two processes and put into
a database. Then, the data is pulled out, some simple statistics are
compiled, and the data is written to a text file for inspection outside
of the program.
This program will compile cleanly with the sample \ref makefile.
\code


#include <apop.h>

//Your processes are probably a bit more complex.
double process_one(gsl_rng *r){
        return gsl_rng_uniform(r) * gsl_rng_uniform(r) ;
}

double process_two(gsl_rng *r){
        return gsl_rng_uniform(r);
}

int main(){
apop_data      *m;
gsl_vector      v1, v2;
double          p1, p2;
int             i;
gsl_rng *       r = apop_rng_alloc(15);

        //create the database and the data table.
        apop_db_open("runs.db");
        apop_table_exists("samples",1); //If the table already exists, delete it.
        apop_query("create table samples(iteration, process, value); begin;");

        //populate the data table with runs.
        for (i=0; i<1000; i++){
                p1      = process_one(r);
                p2      = process_two(r);
                apop_query("insert into samples values(%i, %i, %g);", i, 1, p1);
                apop_query("insert into samples values(%i, %i, %g);", i, 2, p2);
        }
        apop_query("commit;"); //the begin-commit wrapper saves writes to the drive.

        //pull the data from the database. Use the GSL's vector views to minimize copying.
        m  = apop_db_to_crosstab("samples", "iteration","process", "value");
        v1      = gsl_matrix_column(m->matrix, 0).vector;
        v2      = gsl_matrix_column(m->matrix, 1).vector;

        //print info.
        printf("\t   mean\t\t   var\n");
        printf("process 1: %f\t%f\n", apop_mean(&v1), apop_var(&v1));
        printf("process 2: %f\t%f\n\n", apop_mean(&v2), apop_var(&v2));
        printf("t test\n");
        apop_data_show(apop_t_test(&v1,&v2));
        apop_data_print(m, "the_data.txt"); //does not overwrite; appends.
        return 0;
}
\endcode
*/

/** \page windows The Windows page
Get <a href="http://www.cygwin.com">Cygwin</a>. The setup program is
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
X-Window will give you a nicer environment in which to work.  After you
start Cygwin, just type <tt>startx</tt> to bring up a more usable,
nice-looking terminal (and the ability to do a few thousand other things
which are beyond the scope of this documentation).
Once you have Cygwin installed and a good terminal running, you can
follow along with the remainder of the discussion without modification.

Sqlite3 is difficult to build from scratch, but you can get a 
packaged version by pointing Cygwin's install program to the
Cygwin Ports site: http://cygwinports.dotsrc.org/ .

Second, some older (but still pretty recent) versions of Cygwin
have a search.h file which doesn't include the function lsearch().
If this is the case on your system, you will have to update your Cygwin
installation.

Finally, windows compilers often spit out lines like:
\code
Info: resolving _gsl_rng_taus by linking to __imp__gsl_rng_taus (auto-import)
\endcode
These lines are indeed just information, and not errors. Feel free to ignore them.

[Thanks to Andrew Felton and Derrick Higgins for their Cygwin debugging efforts.]
*/

/** \page notroot  Not root? 
If you aren't root, then you will need to create a subdirectory in
your home directory in which to install packages. The GSL and SQLite
installations will go like this. The key is the <tt>--prefix</tt> addition
to the <tt>./configure</tt> command.
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

/** \page c C
 
\section learning  Learning C
<a href="http://avocado.econ.jhu.edu/modeling">Modeling with Data</a> has a full tutorial for C, oriented at users of standard stats packages.
More nuts-and-bolts tutorials are 
<a href="http://www.google.com/search?hl=es&amp;c2coff=1&amp;q=c+tutorial">in abundance</a>.
Some people find pointers to be especially difficult. Fortunately, there's a claymation cartoon which clarifies everything
<a href=http://cslibrary.stanford.edu/104/>here</a>.

\section reference References
For your convenience, here are links to the documentation for the <a href="http://www.gnu.org/software/libc/manual/html_node/">Standard library</a> and the <a href="http://www.gnu.org/software/gsl/manual/gsl-ref_toc.html">GSL</a>.

\section usagenotes  Usage notes
Here are some notes about the technical details of using the Apophenia library.

\subsection headertrick Header aggregation 
If you put 
\verbatim
#include <apop.h>
\endverbatim
at the top of your file, then it will call virtually every header file you could need: gsl_matrix.h, gsl_blas.h, sqlite3.h, stdio.h, string.h, math.h, apophenia_all_of_them.h, et cetera. Of course, if you get `implicit declaration of...' then you will need to manually include something else.
Bear in mind that every book on C will tell you this is bad form and you shouldn't do it.

\subsection liblinking Libraries 
Your best bet is to write yourself a \ref makefile "Makefile".
If you don't want to use the sample \ref makefile "Makefile", then here are some notes for the command line.
When compiling, you will need to tell the compiler to use the Apophenia library. That is, you will need to call GCC with <tt>gcc -lapophenia</tt> (as well as the other usual flags). For example,
<tt>
gcc sample.c -lapophenia -lsqlite3 -lgsl -lgslcblas -o run_me -g 
</tt>
will produce an output file named <tt>run_me</tt> from the input source code <tt>sample.c</tt>. It will include symbols for debugging (<tt>-g</tt>) and will correctly find the functions in the Apophenia, GSL, and SQLite libraries (<tt>-lapophenia -lgsl ...</tt>). 
Order matters in the linking list: the files a package depends on should be listed after the package. E.g., since sample.c depends on Apophenia, <tt>gcc sample.c -lapophenia</tt> will work, while <tt>gcc -lapophenia sample.c</tt> is likely to give you errors. Similarly, list <tt>-lapophenia</tt> before <tt>-lgsl</tt>, which comes before <tt>-lgslcblas</tt>.

\subsection debugging  Debugging
The global variable <tt>apop_opts.verbose</tt> turns on some diagnostics, such as printing the query sent to the databse engine (which is useful if you are substituting in many <tt>\%s</tt>es). Just set <tt>apop_opts.verbose =1</tt> when you want feedback and <tt>apop_opts.verbose=0</tt> when you don't.

\subsection vim Syntax highlighting 
If your text editor supports syntax highlighting, there are a few types defined in the Apophenia and GSL headers which may be worth coloring.
E.g., for <tt>vim</tt>, add the following two lines to <tt>/usr/share/vim/syntax/c.vim</tt>:
\verbatim
syn keyword     cType           gsl_matrix gsl_rng gsl_vector apop_data
syn keyword     cType           apop_name apop_model
\endverbatim
Other text editors have similar files to which you can add the above types.
*/

/** \page makefile Makefile
 
Instead of giving lengthy GCC commands at the command prompt, you can use a Makefile to do most of the work. How to:
 * Copy and paste the following into a file named Makefile (capital M).
 * Make sure the two indented lines begin with a single Tab, instead of seven or eight spaces. 
 * Change the first line to the name of your program (e.g., if you have written <tt>sample.c</tt>, then the first line will read <tt>PROGNAME=sample</tt>). 
 * If your program has multiple <tt>.c</tt> files, just add a corresponding <tt>.o</tt> on the <tt>objects</tt> line, e.g. <tt>sample2.o</tt> <tt>sample3.o</tt>
 * One you have a Makefile in the directory, simply type <tt>make</tt> at the command prompt to generate the executable.
 \verbatim
PROGNAME = your_program_name_here
objects =$(PROGNAME).o
CFLAGS = -g -Wall 
LINKFLAGS = -lapophenia -lgsl -lgslcblas -lsqlite3
c: $(objects)
	gcc $(CFLAGS) $(objects) $(LINKFLAGS) -o $(PROGNAME)
$(objects): %.o: %.c 
	gcc $(CFLAGS) -c $< -o $@
\endverbatim
*/

/** \page sql SQL
 For a reference, your best bet is the <a href="http://www.sqlite.org/lang.html">Structured Query Language reference</a> for SQLite. 
For a tutorial; there is an abundance of <a
href="http://www.google.com/search?q=sql+tutorial">tutorials online</a>.
The blog of Apophenia's author includes an <a
href="http://fluff.info/blog/arch/00000118.htm">entry</a> about
complementarities between SQL and matrix manipulation packages.

Apophenia currently supports two database engines: SQLite and
mySQL. SQLite is the default, because it is simpler and generally more
easygoing than mySQL, and supports in-memory databases.

You can switch to mySQL two ways: set <tt>apop_opts.db_engine = 'm'</tt>,
or set the environment variable <tt>APOP_DB_ENGINE=mysql</tt>. Otherwise,
the system will use SQLite. Ideally, after you make this switch, you need make no other changes---
\ref apop_query, \ref apop_query_to_data, \ref apop_table_exists, et cetera, will work
as before. 

Finally, Apophenia provides a few nonstandard SQL functions to facilitate
math via database; see \ref db_moments.
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


/*      This is debris. Not parsed by doxygen.


\li <b>The world is not linear</b>, so why use a package to fit
linear models? Apophenia facilitates writing \ref mle "likelihood functions" which
can be as crazy as the world you are modeling, and then fits them with
one function call to \ref apop_maximum_likelihood. Comparing competing
models is easy as well.
\li <b>Passing queries directly to an SQL engine.</b> if you have a gigabyte of data,
you do not want it in an in-memory database, until the last possible
minute---you're probably only using three variables out of six hundred
anyway. If your data comes from multiple sources, you will need a means
of merging the sources. In short, data management is best done via a
database. \ref apop_query allows you to dynamically write queries
that are passed directly to SQLite.
\li <b>Converting data is a pain</b>. It's funny and a bit wrong that
college and even grad school statistics classes present the student with
perfectly-cleaned data sets and spend 90% of their time on calculating
statistics and testing hypotheses; while in the real world, we spend
about 90% of our time getting, cleaning, and managing our data, and
once it's all perfect we do the statistics with a single function call. 
Apophenia provides a good number of \ref conversions to make
data wrangling as painless as it can be in the real world. The
\ref convertfromtext alone may be worth the time it takes you to install
the library.
\li <b>The GSL is just a step shy.</b> Don't get me wrong: it's a great
library, which drives much of Apophenia (see below). But the authors did not have
statistical analysis in mind, and often chose generality over ease of
use. For example, doing a singular value decomposition theoretically
just requires a single call to \c gsl_linalg_SV_decomp, but in reality
takes many steps of data massaging to get things in place. But you can
feed your data directly to \ref apop_sv_decomposition and get a principal
component space back immediately. Think of Apophenia as the lazy
statistician's supplement to the GSL.

\section datamanagement Data management and data analysis

That said, there are two simple goals of Apophenia:
\li Minimize annoyances.
\li Estimate model parameters and test hypotheses regarding those parameters.

With regard to the first goal, the hope is that every time you are writing
a program and think "darn it, now I have to write a tedious function to
massage my data" that function is already in the library. 
[If it isn't, you can write that function and contribute it.] The
library thus includes a number of functions to reformat, convert, and
do simple-but-tedious calculations on data.

With regard to the goal of estimation, Apophenia provides three
*/



/**  \page outline An outline of the library


Outlineheader dataoverview Data sets

The \ref apop_data structure represents a data set.  It joins together
a \c gsl_vector, a \c gsl_matrix, an \ref apop_name, and a table of strings. It
tries to be lightweight, so you can use it everywhere you would use a \c gsl_matrix or a \c gsl_vector.

For example, let us say that you are running a regression: there is
a vector for the dependent variable, and a matrix for the dependent
variables. Think of them as a partitioned matrix, where the vector is column -1, and the first column of the matrix is column zero. Here is some code to print the entire matrix. Notice that the column counter \c i starts counting at -1.

  \code
  for (j = 0; j< data->matrix->size1; j++){
    printf("%s\t", apop_name_get(data->names, j, 'r'));
    for (i = -1; i< data->matrix->size2; i++)
        printf("%g\t", apop_data_get(data, j, i));
    printf("\n");
    }
    \endcode

We're generally assuming that the data vector and data matrix have the
same row count: \c data->vector->size==data->matrix->size1 . This means
that the \ref apop_name structure doesn't have separate vector_names
and row_names elements: the rownames are assumed to apply for both.

Easy data manipulation is essential for enjoying life as a researcher.
Thus, there are a great many functions to collate, copy, merge, sort, prune, and otherwise manipulate the \ref apop_data structure and its components.


    \li\ref apop_data_add_named_elmt()
    \li\ref apop_data_alloc()
    \li\ref apop_data_calloc()
    \li\ref apop_data_copy()
    \li\ref apop_data_fill()
    \li\ref apop_data_free()
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
    \li\ref apop_matrix_rm_columns()
    \li\ref apop_matrix_stack()
    \li\ref apop_text_add()
    \li\ref apop_text_alloc()
    \li\ref apop_text_free()
    \li\ref apop_vector_bounded()
    \li\ref apop_vector_copy()
    \li\ref apop_vector_fill()
    \li\ref apop_vector_increment()
    \li\ref apop_vector_stack()
    \li\ref apop_vector_vfill()


Outlineheader setgetsec Set/get

Let \f$t\f$ be text, and \f$i\f$ be an index. Then you can get or set
an element by row title and column index, row index and column title,
et cetera. The functions that have no suffix find an element by row and column index.

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



endofdiv

Outlineheader  matrixmathtwo  Basic Math

        \li\ref apop_vector_exp()
        \li\ref apop_vector_log()
        \li\ref apop_vector_log10()
        \li\ref apop_vector_distance ()
        \li\ref apop_vector_grid_distance ()
        \li\ref apop_vector_normalize ()
        \li\ref apop_matrix_normalize ()

endofdiv
            
Outlineheader  matrixmath  Matrix math

        \li\ref apop_det_and_inv()
        \li\ref apop_dot()
        \li\ref apop_matrix_determinant()
        \li\ref apop_matrix_inverse()

endofdiv

Outlineheader  sumstats  Summary stats

        \li\ref apop_data_summarize ()
        \li\ref apop_vector_moving_average()
        \li\ref apop_vector_percentiles()

endofdiv


Outlineheader  moments  Moments


            \li\ref apop_data_correlation ()
            \li\ref apop_data_covariance ()
            \li\ref apop_matrix_correlation ()
            \li\ref apop_matrix_covariance ()
            \li\ref apop_matrix_mean ()
            \li\ref apop_matrix_mean_and_var ()
            \li\ref apop_matrix_sum ()
            \li\ref apop_matrix_summarize ()
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

Outlineheader mapplysec   Map/apply

These functions allow you to send each element of a vector or matrix to
a function, either producing a new matrix (map) or transforming a new one (apply). 
The \c ..._sum functions return the sum of the mapped output, which is nice for finding the count of missing elements (by mapping \c isnan to your data) or getting the total log likelihood of a set of observations.


            \li\ref apop_matrix_apply()
            \li\ref apop_matrix_map()
            \li\ref apop_matrix_map_all_sum()
            \li\ref apop_matrix_map_sum()
            \li\ref apop_vector_apply()
            \li\ref apop_vector_map()
            \li\ref apop_vector_map_sum()

endofdiv

Outlineheader names   Name handling

                \li\ref apop_name_add()
            \li\ref apop_name_alloc()
            \li\ref apop_name_copy()
            \li\ref apop_name_cross_stack()
            \li\ref apop_name_find()
            \li\ref apop_name_free()
            \li\ref apop_name_print()
            \li\ref apop_name_stack()


endofdiv

Outlineheader fact   Generating factors

            \li\ref apop_data_to_dummies()
            \li\ref apop_text_to_factors()
            \li\ref apop_text_unique_elements()
            \li\ref apop_vector_unique_elements()

endofdiv

endofdiv

Outlineheader dbs Databases


These are convenience functions to handle interaction with SQLite. They
open one and only one database, and handle most of the interaction
therewith for you.

You will probably first use \ref apop_text_to_db to pull data into
the database, then \ref apop_query to clean the data in the database,
and finally \ref apop_query_to_data to pull some subset of the data
out for analysis.


See also the \ref conversions, including \ref apop_text_to_db and \ref apop_matrix_to_db.

\li \ref apop_query : Manipulate the database, return nothing (e.g., input data).
\li \ref apop_db_open : Optional, for when you want to use a database on disk.
\li \ref apop_db_close : If you used \ref apop_db_open, you will need to use this too.
\li \ref apop_table_exists : Check to make sure you aren't reinventing or destroying data. Also, the clean way to drop a table.
\li \ref apop_count_cols : Count the columns in a table.
\li \ref apop_db_merge : Import or merge the whole of another database into the currently open db.
\li \ref apop_db_merge_table : Import/merge just one table.
   \li \ref apop_crosstab_to_db: Convert between two common data layouts
    \li \ref apop_db_rng_init: Apophenia maintains a database-side RNG; set its seed it with this.

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

        See the print functions below.

endofdiv

    Outlineheader dbmath Math in the db  \dbmath
    
        \li\ref apop_db_paired_t_test()
        \li\ref apop_db_t_test()

endofdiv

endofdiv

Outlineheader Mode Models

    Outlineheader modellist  Models that ship with Apophenia


        Outlineheader Dist Distributions

            \li\ref apop_bernoulli
            \li\ref apop_beta
            \li\ref apop_beta_from_mean_var()
            \li\ref apop_binomial
            \li\ref apop_exponential
            \li\ref apop_gamma
            \li\ref apop_gaussian
            \li\ref apop_improper_uniform
            \li\ref apop_lognormal
            \li\ref apop_multivariate_normal
            \li\ref apop_normal
            \li\ref apop_poisson
            \li\ref apop_uniform
            \li\ref apop_waring
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

A model is an equation (or system of equations) that relies on data and
has unknown parameters to be determined. [Notice that this definition
readily includes "non-parametric" models.] Much of statistical analysis
consists of writing down a model, estimating its parameters, and running
hypothesis tests to determine the confidence with which we can make
statements about those parameters.
<!--\footnote{Many statistics packages include a model structure that describes only
linear models.  "Linear" models can include a wide range of nonlinear
features, but they are still a subset of measure zero within the class of
models as described above. Currently, Apophenia has no plans to include
a summary syntax for describing linear models; the reader who has a linear
model to be estimated via OLS and friends is advised to instead manipulate
the data set to the appropriate form and then call \ref apop_ols, \ref apop_iv, et cetera.}-->

Apophenia facilitates this via its \ref apop_model objects. Each object is
a model as described above, and includes a method named <tt>estimate</tt>
which takes in data and returns an \ref apop_model which includes
the parameter estimates and the characteristics one would need for
hypothesis testing.

For example, a model may be a probability distribution. The data is assumed
to have been drawn from a given distribution and the question is
only what distributional parameters best fit; e.g., assume the data
is Normally distributed and find the mean and variance.

The design of the objects hopes to make it as easy as possible for you,
dear reader, to write new models. For the most part, all you need to
do is write a log likelihood function, and \ref apop_maximum_likelihood
does the rest; see below.

The main function in the systems below is the \ref apop_estimate
function. It takes in a model and data, and outputs an apop_estimate,
that includes the parameter estimates and the various auxiliary data
that one may need to test the estimates, such as the variance-covariance
matrix. For most users, the apop_model.estimate function will be all
one needs from a model. Just prep the data, select a model, and produce
an estimate:

\code
apop_data 	    *data 		    = read_in_data();
apop_model 	*the_estimate 	= apop_estimate(data, apop_probit);
apop_estimate_print(the_estimate);
\endcode

Outlineheader internals The internals

\image html http://apophenia.sourceforge.net/doc/model.png
\image latex model.png

The \ref apop_model struct breaks down into three parts:

\li Info like names, the pointer to the input data, and the parameters, are for the most part self-descriptive.

\li There is a broad class of functions that cover most of what you
would do with a model. You can see that there is a bit of a bias toward
maximum likelihood estimation. There are helper functions for most of
them, like the \ref apop_estimate function above, that meant that you
never had to directly handle the model's \c estimate method. 
In the list at the top of this page, 
all of the functions whose last argument is a model are helper functions
of this type. The helper functions do some boilerplate error checking,
and mean that you don't have to fill in every blank in your model: if
you have a \c log_likelihood method but no \c p method, then \ref apop_p
will use exp(\c log_likelihood). If you don't give an \c estimate method,
then \c apop_estimate will use maximum likelihood.

\li I refer to the values output from an estimation as \em parameters,
and the details of how the model's machinery work as \em settings. The
parameters are pretty standardized, and a \ref apop_data set is sufficient
to handle the great majority of models. However, the settings of an OLS
regression are drastically different from those of a histogram or an MLE.
The solution is a list of settings structures. This is probably the
least pleasant structural element in the system, to the point that you will not want to handle the settings structures directly. Instead, you should use the various \ref settings "helper functions for model settings".

endofdiv

Outlineheader write_likelihoods Writing your own 

Writing apop_model objects is easy. Here is the procedure for an MLE, possibly the most involved
type of model estimation in common use:


\li Write a likelihood function. Its header will look like this:
\code
double apop_new_log_likelihood(apop_data *data, apop_model *m)
\endcode 
where \c data is the input data, and \c
m is the parametrized model (i.e., your model with a set \c parameters element). 
This function will return the value of the log likelihood function at the given parameters.
\li Is this a constrained optimization? See the \ref constraints "Constraints page" on how to set them.
\li Write the object. In your header file, include 
\code
apop_model apop_new_likelihood = {"The Me distribution", number_of_parameters, 
            .estimate = new_estimate, .log_likelihood = new_log_likelihood };
\endcode
If there are constraints, then replace the appropriate <tt>NULL</tt> with the right constraint function: beta_zero_and_one_greater_than_x_constraint
\c number_of_parameters is probably a positive integer like \c 2, but
it is often (the number of columns in your data set) -1, in which case,
set \c number_of_parameters to \c -1.
\li Test. Debug. Retest.
\li (optional) Write a gradient for the log likelihood function. This
typically involves calculating a derivative by hand, which is an easy
problem in high-school calculus. The function's header will look like: 
\code
void apop_new_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m)
\endcode 
where \c m and \c d are as above, and \c gradient is a \c gsl_vector to be filled with the gradient of the parameters. 
At the end of this function, you will have to assign the appropriate derivative to every element of the gradient vector:
\code
gsl_vector_set(gradient,0, d_a);
gsl_vector_set(gradient,1, d_b);
\endcode 
Now add the resulting dlog likelihood function to your object, by adding a \c .dlog_likelihood=your_function to the declaration.
\li Send the code to the maintainer for inclusion in future versions of Apophenia.

endofdiv


    Outlineheader Basicmodelmethods Basic object methods

        \li\ref apop_model_clear()
        \li\ref apop_model_copy()
        \li\ref apop_model_free()
        \li\ref apop_model_set_parameters()
        \li\ref apop_model_show()


    endofdiv

    Outlineheader mathmethods Methods for computation

        \li\ref apop_draw()
        \li\ref apop_estimate()
        \li\ref apop_expected_value()
        \li\ref apop_log_likelihood()
        \li\ref apop_score()
        \li\ref apop_model_prep()
        \li\ref apop_model_print()
        \li\ref apop_p()

    endofdiv

    Outlineheader modelsettings Model settings


Apophenia is really only based on two objects, the \ref apop_data
set and the \ref apop_model. Data sets come in a pretty standard form,
so the data object is basically settled. But describing a statistical,
agent-based, social, or physical model in a reasonably standardized form
is much more difficult, primarily because every model has significantly
different settings. E.g., an MLE requires a method of search (conjugate
gradient, simplex, simulated annealing), and a histogram needs the number
of slots to be filled with data.

So, the \ref apop_model includes 
a single list, whose name is simply \c settings, which can hold an arbitrary
number of groups of settings. For example, you can have a set of
closed-form variables for estimating the log likelihood, and a histogram
for making random draws.

To get/set a setting, you would need to specify the model, the settings
group, and the name of the setting itself.  For the sake of giving you
a mental model, this is much like a g_lib config file (or for Windows users,
a .ini file), which would have settings divided into sections, like

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

If you don't need to change settings from the default, you don't need to
care about any of this, because
\code
apop_model_show(apop_estimate(your_data_here, apop_ols));
\endcode
---still works fine.

If you do need to change settings, then the process takes another step
or two. Here's a sample:

\code
1 apop_data *data = your_data_here;
2 apop_data *w = your_weights_here;

3 apop_model *m = apop_model_copy(apop_wls);
4 Apop_settings_add_group(m, apop_ls, data);
5 Apop_settings_add(m, apop_ls, weights, w);
6 apop_model *est = apop_estimate(data, *m);
\endcode

Line three establishes the baseline form of the model. Line four adds
a settings group of type \ref apop_ls_settings to the model. If you check the
manual for \ref apop_ls_settings_alloc, you'll see that it takes one
argument---the input data---and the inputs to the alloc function appear
after the model and settings group name. Line five sets the weights
element in that group to w, and looks like the config-file format above.

Also, the output to any of Apophenia's estimations will have an
appropriate group of settings allocated, so you can chain estimations
pretty easily. Continuing the above example, you could re-estimate with
an alternative set of weights via:

\code
Apop_settings_add(est, apop_ls, weights, weight_set_two);
apop_model *est2 = apop_estimate(data, *est);
\endcode

[Notice that the 'add' function really just modifies.]

Some people are incredibly averse to additional lines of code, so there's
a convenience to allocate and modify a setting on the same line (joining
lines 4 and 5 above):

\code
Apop_settings_alloc_add(m, apop_ls, weights, w, data);
\endcode

If you need to read a setting, such as to check whether things have
changed after an estimation, use:

\code
Apop_settings_get(m, apop_ls, weights);
\endcode

Notice the use of a single capital to remind you that you are using a
macro, and so surprising errors can crop up. Here in the modern day, we
read things like APOP_SETTINGS_ALLOC_ADD as yelling, but if you prefer
all caps to indicate macros, those work as well.

For just using a model, that's about 100% of what you need to know.


Outlineheader settingswritng  Writing new settings

To store the settings for your own models, you don't necessarily
need any of this. The \ref apop_model structure has a \c void pointer named
\c more which you can use as you see fit. If \c more_size is larger than zero
(i.e., you set it to \c sizeof(your_struct)), then it will be copied via
\c memcpy as necessary. Apohenia's estimation routines will never impinge
on this item, so do what you feel with it.

If you do want to set up a new model, then you will need four items. 
This is the sort of boilerplate that will be familiar to users of 
object oriented languages in the style of C++ or Java. Let
your settings group be named \c ysg; then you will need

\li The settings struct
\code
typedef struct {
...
} ysg_settings
\endcode


\li The allocate function
\code
ysg_settings *ysg_settings_alloc(...){ 
    ysg_settings *out = malloc(sizeof(ysg_settings));
    // Give default values for all elements here. 
    return out; }
\endcode


\li The copy function
\code
void *ysg_settings_copy(ysg_settings *copyme) {
    ysg_settings *out = malloc(sizeof(ysg_settings));
    //copy elements (or perhaps memcpy the whole struct) here.
    return out; }
\endcode


\li The free function, which can be as brief as:
\code
void ysg_settings_free(ysg_settings *copyme) {
    free(copyme);
}
\endcode
but may include a freeing of pointed-to subelements as necessary.

The names are not negotiable: when you call
\code
Apop_settings_alloc(m, ysg, ...)
\endcode
the macro will look for \c ysg_settings, \c ysg_settings_alloc, et cetera.

The lines-of-code averse will cringe at having to write such boilerplate
code (I do), but after spending a year resisting it, I have to concede
that it's the least of all evils.

You can retrieve the whole group or individual elements via:
\code
Apop_settings_get_group(m, ysg)
//as well as the one-element macro mentioned above,
Apop_settings_get(m, ysg, an_element)
\endcode

As you saw above, once the typedef/alloc/copy/free machinery is written,
you can declare, get, and set in a reasonably graceful manner.


endofdiv


        \li\ref Apop_settings_add
        \li\ref Apop_settings_add_group()
        \li\ref Apop_settings_alloc
        \li\ref Apop_settings_alloc_add()
        \li\ref apop_settings_copy_group()
        \li\ref Apop_settings_get()
        \li\ref apop_settings_get_group()
        \li\ref Apop_settings_get_group()
        \li\ref apop_settings_group_alloc()


        Outlineheader Spec Specific settings

            \li\ref apop_histogram_settings_alloc()  
            \li\ref apop_ls_settings_alloc()
            \li\ref apop_mle_settings_alloc() \anchor mlesettingsa
            \li\ref apop_update_settings_alloc()

        endofdiv

    endofdiv

    Outlineheader Esti Estimation aids

        \li\ref apop_model_fix_params: hold some parameters constant



    endofdiv

endofdiv

Outlineheader Test Tests & diagnostics

    \li\ref apop_paired_t_test()
    \li\ref apop_f_test()
    \li\ref apop_t_test()
    \li\ref apop_test_anova_independence()
    \li\ref apop_test_fisher_exact()
    \li\ref apop_test_kolmogorov()
    \li\ref apop_plot_qq()
    \li\ref apop_estimate_coefficient_of_determination()
    \li\ref apop_estimate_r_squared()
    \li\ref apop_estimate_parameter_t_tests()
    
    See also the in-database tests above.

    Outlineheader Mont Monte Carlo methods
    
        \li\ref apop_bootstrap_cov()
        \li\ref apop_jackknife_cov()

    endofdiv

endofdiv

Outlineheader Histosec Histograms


The GSL provides a few structures that basically accumulate data into
bins. The first is the <tt>gsl_histogram</tt> structure, that produces a PMF.

To produce a PMF from \c your_data, use \ref apop_histogram "apop_estimate(your_data, apop_histogram)", then produce a synced histogram of other data (observed or theoretical) using 
\ref apop_histogram_vector_reset or 
\ref apop_histogram_model_reset.


The second structure from the GSL incrementally sums up the PMF's bins to
produce a CMF. The CMF can be used to map from a draw from a Uniform[0,1]
to a draw from the PMF.  Because it can be used to draw from the PMF,
the GSL calls this the <tt>gsl_histogram_pdf</tt> structure. That's right:
the the data in the <tt>gsl_histogram_pdf</tt> structure is a cumulative
sum---a CMF.

Anyway, here are some functions to deal with these various histograms and such.


    \li\ref apop_histogram_model_reset()
    \li\ref apop_histogram_moving_average()
    \li\ref apop_histogram_normalize()
    \li\ref apop_histogram_vector_reset()
    \li\ref apop_histograms_test_goodness_of_fit()


endofdiv

Outlineheader Maxi Maximum likelihood methods

If you read the section on writing models, then you already know how to
do maxium likelihood on exotic setups. Just write a model that has a
\tt p or \tt log_likelihood function, and call \c apop_estimate(your_data,your_model).
The default estimation routine is maximum likelihood.

MLEs have an especially large number of parameter tweaks that could be made; see the section on MLE settings above.

Outlineheader constr Setting Constraints

The problem is that the parameters of a function must not take on certain values, either because the function is undefined for those values or because parameters with certain values would not fit the real-world problem.

The solution is to rewrite the function being maximized such that the function is continuous at the constraint boundary but takes a steep downward slope. The unconstrained maximization routines will be able to search a continuous function but will never return a solution that falls beyond the parameter limits.

If you give it a likelihood function with no regard to constraints plus an array of constraints, 
\ref apop_maximum_likelihood will combine them to a function that fits the above description and search accordingly.

A constraint function must do three things:
\li It must check the constraint, and if the constraint does not bind (i.e., the parameter values are OK), then it must return zero.
\li If the constraint does bind, it must return a penalty, that indicates how far off the parameter is from meeting the constraint.
\li if the constraint does bind, it must set a return vector that the likelihood function can take as a valid input. The penalty at this returned value must be zero.

The idea is that if the constraint returns zero, the log likelihood
function will return the log likelihood as usual, and if not, it will
return the log likelihood at the constraint's return vector minus the
penalty. To give a concrete example, here is a constraint function that
will ensure that both parameters of a two-dimensional input are both
greater than zero:

\code
static double beta_zero_greater_than_x_constraint(apop_data *returned_beta, apop_model *v){
    //constraint is 0 < beta_2
  static apop_data *constraint = NULL;
    if (!constraint){
        constraint= apop_data_calloc(1,1,1);
        apop_data_set(constraint, 0, 0, 1);
    }
    return apop_linear_constraint(v->parameters->vector, constraint, 1e-3);
}
\endcode

endofdiv

    \li\ref apop_estimate_restart()
    \li\ref apop_linear_constraint()
    \li\ref apop_maximum_likelihood()
    \li\ref apop_numerical_covariance_matrix()
    \li\ref apop_numerical_gradient()


endofdiv

Outlineheader Asst More descriptive methods

    \li\ref apop_matrix_pca: Principal component analysis
    \li\ref apop_anova: One-way or two-way ANOVA tables
    \li\ref apop_update : Bayesian updating


endofdiv

Outlineheader Miss Missing data

    \li\ref apop_data_listwise_delete()
    \li\ref apop_ml_imputation()

endofdiv

Outlineheader Legi Legible output

    \li\ref apop_plot_histogram()
    \li\ref apop_plot_line_and_scatter()
    \li\ref apop_plot_lattice()
    \li\ref apop_data_print()  \legout
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
    \li\ref apop_opts_memcpy()
    \li\ref apop_strip_dots()
    \li\ref apop_system()

endofdiv

Outlineheader Math Math utilities

    \li\ref apop_generalized_harmonic()
    \li\ref apop_rng_alloc()
    \li\ref apop_rng_GHgB3()
    \li\ref apop_random_double ()
    \li\ref apop_random_int ()

endofdiv

Outlineheader Prob Deprecated

These functions will probably disappear or be replaced soon.

    \li\ref apop_estimate_fixed_effects_OLS()
    \li\ref apop_kernel_density_settings_alloc()
    \li\ref apop_test_chi_squared_var_not_zero ()
    \li\ref apop_two_tailify()
    \li\ref apop_count_cols()
    \li\ref apop_data_to_db()

endofdiv

endofdiv

*/
