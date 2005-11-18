/** \mainpage The Apophenia documentation

This is the documentation for <a href="http://apophenia.sf.net">Apophenia</a>. 

\section Prerequisites 

 \li \ref intro "Intro": The motivation for the package.
 \li \ref setup "Setup":  Installing GCC, the GSL, SQLite, and Apophenia itself.

\section Doing statistics 

 \li \ref basic_stats "Basic statistics": Mean, variance, percentiles, &c.
 \li \ref regression  
 \li \ref ttest 
 \li \ref likelihood_fns "Maximum likelihood estimation": estimators requiring a search for the maximum of a likelihood function.
 \li \ref linear_algebra "Linear Algebra": determinants, projections, numerical gradients, &c. Some convenience functions to display matrices and vectors.

\section Shunting data 

 \li \ref db "Database utilities": an easy front end to SQLite
 \li \ref conversions
 \li \ref output "Output functions": Summarize data and print tables to the screen or a file.
 \li \ref command_line "Command-line utilities": a still easier front end to SQLite, and doing stats via Perl/Python/shell script.


\section Speaking the languages

 \li \ref c "C": For statisticians who don't know/are rusty with the C
 programming language. Those who know their C may want to check the
 Apophenia-specific notes in the \ref usagenotes "usage notes" section.
 \li \ref sql About the syntax for querying databases.
 \li \ref types "Types": New structures defined by Apophenia
 */

/** \page intro Why?

\section datamanagement Data management and data analysis
We require two things from a good statistics package: easy management of
large data sets, and the ability to crunch numbers on a large scale. [The
third place entry, \ref graphing "graphing", is discussed on another page.]

Apophenia facilitates data management by including a database
interface. By reading your data into a database instead of an in-memory
matrix, you effectively have no limits on the size of your data set.

Apophenia also makes heavy use of the GNU Scientific Library, which is
a well-optimized system for processing large matrices of numbers. Once
you have pulled out of the main data set that subset which you wish to
analyze, the GSL and Apophenia's functions will have no problem estimating
models to fit the data.

Thus, the average analysis using Apophenia would take the following steps:

 \li read the data into the database using \ref apop_convert_text_to_db 
 \li use SQL queries handled by \ref apop_query to massage the data as needed
 \li use \ref apop_query_to_matrix to pull the data into an in-memory matrix
 \li call a regression function such as \ref apop_OLS or a maximum likelihood estimator such as a \ref apop_probit "probit" to fit parameters to the data.
 \li use the results for further analysis, or just dump them to the screen with \ref apop_estimate_print.

If this seems a bit vague, have a look at this \ref sample_program.

\section sell Some complaints alleviated

\li <b>The world is not linear</b>, so why are you using a package to fit
linear models? Apophenia facilitates writing \ref likelihood_fns which
can be as crazy as the world you are modeling, and then fits them with
one function call to \ref apop_maximum_likelihood. Comparing competing
models is easy as well.

\li <b>Passing queries directly to SQL.</b> if you have a gigabyte of data,
you do not want it in an in-memory database until the last possible
minute---you're probably only using three variables out of six hundred
anyway. If your data comes from multiple sources, you will need a means
of merging the sources. In short, data management is best done via a
database. \ref apop_query will allow you to dynamically write queries
that are passed directly to SQLite.

\li <b>Converting data is a pain</b>. It's funny and a bit wrong that
college and even grad school statistics classes present  the student with
perfectly-cleaned data sets and spend 90% of their time on calculating
statistics and testing hypotheses; while in the real world, we spend
about 90% of our time getting, cleaning, and managing our data, and
once it's all perfect we do the statistics with a single function call. 
Apophenia provides a good number of \ref conversions to make
data wrangling as painless as it can be in the real world. The
\ref convertfromtext alone may be worth the time it takes you to install
apophenia.

\li <b>The GSL is just a step shy.</b> Don't get me wrong: it's a great
library, which drives much of Apophenia. But the authors did not have
statistical analysis in mind, and often chose generality over ease of
use. For example doing a singular value decomposition theoretically
just requires a single call to \c gsl_linalg_SV_decomp, but in reality
takes many steps of data massaging to get things in place. But you can
feed your data directly to \ref apop_sv_decomposition and get a principal
component space back immediately. Think of Apophenia as the lazy
statistician's supplement to the GSL.

<!--
\section librant Stats libraries vs Stats packages 
The reason I (BK) started up this library is that I was sick of learning
new languages. There are a few dozen statistics packages to choose from,
all of which are best for one task or another, but none of which was
portable and versatile enough to work for every project I had.

In most of the programming world, the evolution has been that there
a few core libraries emerged to handle the specifics of the task at hand, and
then front-ends developed to suit various tastes and specialties. In
the stats world, the process has been entirely backward: a few dozen packages
have started from scratch whith entirely new libraries of essentialy
the same functions. Each has its own specialty, each is tied to only
one front end, and each requires learning a brand new language.
Most are designed around students, and so are easy to learn but virtually
useless for one who knows stats well and wants to use the package on
an industrial-strength data set---several gigabytes of data from the
US Census or USGS or what-have-you.

Meanwhile, C has absolutely no built-in memory limitations, and using a
database makes handling large data sets still more manageable. Also, C is
the most supported programming language on Earth. Type "C tutorial" into
your favorite search engine and you will get hundreds of choices. Even
though it requires understanding one concept which is missing from more
idiot-proof lanugages (pointers), I found that it is worth the marginal
extra time spent learning because it means not having to waste time
learning yet another language when the next project comes around.

To use the database half of Apophenia of course requires learning SQL
(Structured query language), but SQL is not quite a complete language;
it is more a means of describing the rows and columns of a table. The
reader can have the basics of SQL down in a few minutes (since <tt>select
row from table</tt> is almost English), and the return on that investment
 is again huge, since so many programs support SQL.

Thus, Apophenia does not impose its own language for you to learn, but
instead relies on the existing grammar and structure of C and SQL. It
therefore ties in seamlessly with the thousands of existing C libraries,
is documented in part by the hundreds of books and web pages which
teach C and SQL, and can be used as the base for further extensions. 
The command-line utility makes it easy to call Apophenia's functions from
any other language with which you are familiar, such as Perl, Python,
or even bash.

Finally, there are also a few of the more typical reasons to write a
new library. This libarary focuses on MLE methods over linear projection
methods that most packages tend to rely on. Because it focuses on
likelihoods and distributions, the package provides certain procedures,
such as random draws from a Waring distribution, which to my knowledge
don't yet exist in Matlab, R, STATA, &c. 
-->*/

/** \page setup Setting up

\section cast The supporting cast 

To use Apophenia, you will need to have a working C compiler, the GSL and SQLite installed. 

\subsection C C
The <a href="http://gcc.gnu.org/">Gnu Compiler Collection</a> (GCC) is certainly available for your system. If you are using a unix-type system, it is probably already installed. Windows users, see \ref windows .

\subsection gsl The GSL 
This is often available as a package: try <tt>apt-get gsl-devel</tt>,
<tt>urpmi gsl-devel</tt>, or whatever means you use to install a
package. If all else fails, download the source code from the 
<a href="http://sources.redhat.com/gsl">GSL home page</a>. Compilation and
installation is simple:

\verbatim
tar xvzf pkg.tgz  #change pkg.tgz to the appropriate name
cd package_dir    #same here.
./configure
make
su -c "make install"   #see below if you don't have root priveleges.
\endverbatim


If you don't have root priveleges, then see the \ref notroot "not root" section.


\subsection sqlite SQLite
SQLite is gaining popularity, and so may or may not be an available
package for your system. If not, you will need to download it from
<a href="http://www.sqlite.org">the SQLite home page</a> and compile it. The
instructions are identical to compiling the GSL above. You will need to
set your <tt>LD_LIBRARY_PATH</tt>; see below.

\section apop_install The package 

If you're reading this on a web browser, you can download Apophenia from <a href="https://sourceforge.net/project/showfiles.php?group_id=130901">here</a>; else, go to apophenia.info and follow the download link. The library installation is just like the above:

\verbatim
tar xvzf apophenia-0.01.tgz
cd apophenia_install
./configure
make
su -c "make install"
\endverbatim

Notice that, after making the object files, we again set the <tt>LD_LIBRARY_PATH</tt> environment variable to include this directory.

\subsection libpath LD_LIBRARY_PATH 

During one or many of the installation steps above, you probably got a warning that looks like this:

\verbatim
----------------------------------------------------------------------
Libraries have been installed in:
   /usr/local/lib
 
If you ever happen to want to link against installed libraries
in a given directory, LIBDIR, you must either use libtool, and
specify the full pathname of the library, or use the `-LLIBDIR'
flag during linking and do at least one of the following:
[...]
----------------------------------------------------------------------
\endverbatim

If you did, then that means you will need to tell the linker where to look for your newly installed libraries, by adding a line in your shell's configuration file to search for the libraries you'd just installed. Cutting and pasting the following to the command prompt should do the trick. You will only need to do it once.

If your operating system's name ends in the letter X, use this:
\verbatim 
echo "export LD_LIBRARY_PATH=/usr/local/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
\endverbatim 

If you are using Cygwin:
\verbatim 
echo "export PATH=/usr/local/lib:\$PATH" >> ~/.bashrc
echo "export LIBRARY_PATH=/usr/local/lib:\$LIBRARY_PATH" >> ~/.bashrc
\endverbatim 

These commands add a line to your <tt>.bashrc</tt> file; having modified it, reload it either by restarting your shell or the command <tt>source ~/.bashrc</tt> .

\subsection testing Testing

There is a short, complete program in the \ref apop_OLS "apop_OLS" entry which runs a simple OLS regression on a data file. Follow the instructions there to compile and run. There is also a slightly longer \ref sample_program on a separate page.

\subsection using Using 

Now that it's installed, do some stats! For those who are new to C, there are some notes on the \ref C page; if you are familiar with the process, see the \ref usagenotes "usage notes" for Apophenia itself.

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
#include <apophenia/headers.h>

//Your processes are probably a bit more complex.
double process_one(gsl_rng *r){
        return gsl_rng_uniform(r) * gsl_rng_uniform(r) ;
}

double process_two(gsl_rng *r){
        return gsl_rng_uniform(r);
}

int main(){
gsl_matrix      *m;
gsl_vector_view view1, view2;
gsl_vector      *v1, *v2;
double          p1, p2;
int             i;
gsl_rng *       r;

        //set up the GSL's random number generator.
        gsl_rng_env_setup();
        r=gsl_rng_alloc(gsl_rng_default);

        //create the database and the data table.
        apop_open_db("runs.db");
        apop_table_exists("samples",1); //If the table already exists, delete it.
        apop_query_db("create table samples(iteration, process, value); begin;");

        //populate the data table with runs.
        for (i=0; i<1000; i++){
                p1      = process_one(r);
                p2      = process_two(r);
                apop_query_db("insert into samples values(%i, %i, %g);", i, 1, p1);
                apop_query_db("insert into samples values(%i, %i, %g);", i, 2, p2);
        }
        apop_query_db("commit;"); //the begin-commit wrapper saves writes to the drive.

        //pull the data from the database. Use the GSL's vector views to minimize copying.
        m  = apop_db_to_crosstab("samples", "iteration","process", "value", NULL, NULL);
        view1      = gsl_matrix_column(m, 0);
        view2      = gsl_matrix_column(m, 1);
        v1      = &(view1.vector);
        v2      = &(view2.vector);


        //print info.
        printf("\t   mean\t\t   var\n");
        printf("process 1: %f\t%f\n", apop_mean(v1), apop_var(v1));
        printf("process 2: %f\t%f\n\n", apop_mean(v2), apop_var(v2));
        printf("t test p-value: %f\n", apop_t_test(v1,v2));
        apop_matrix_print(m,"\t", "the_data.txt"); //does not overwrite; appends.
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
 \li cvs -- to partake in the versioning system
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
The Apophenia manual has a full tutorial for C; for those reading this online, it is available as a separate document <a href="http://www.fluff.info/klemens/c_crash.pdf">here</a>.

Or gosh, just search Google for <a href="http://www.google.com/search?hl=es&amp;c2coff=1&amp;q=c+tutorial">C tutorial</a>.

Some people find pointers to be especially difficult. Fortunately, there's a claymation cartoon which clarifies everything
<a href=http://cslibrary.stanford.edu/104/>here</a>.

\section reference References
For your convenience, here are links to the documentation for the [http://www.gnu.org/software/libc/manual/html_node/ Standard library] and the [http://www.gnu.org/software/gsl/manual/gsl-ref_toc.html GSL].

\section usagenotes  Usage notes
Here are some notes about the technical details of using the Apophenia library.

\subsection headertrick Header aggregation 
If you put 
\verbatim
#include <apophenia/headers.h>
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
The global variable <tt>apop_verbose</tt> turns on some diagnostics, such as printing the query sent to the databse engine (which is useful if you are substituting in many <tt>\%s</tt>es). Just set <tt>apop_verbose =1</tt> when you want feedback and <tt>apop_verbose=0</tt> when you don't.

\subsection vim Syntax highlighting 
If your text editor supports syntax highlighting, there are seven types defined in the Apophenia and GSL headers which may be worth coloring.

For <tt>vim</tt>, for example. add the following two lines to <tt>/usr/share/vim/syntax/c.vim</tt>:
\verbatim
syn keyword     cType           gsl_matrix gsl_vector apop_inventory apop_estimate 
syn keyword     cType           apop_name gsl_matrix_view gsl_vector_view
\endverbatim
Other text editors have similar files to which you can add the above seven types.
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
Until some notes show up here, your best bet is the <a href="http://www.sqlite.org/lang.html">Structured Query Language reference</a> for SQLite. 

This is a reference, not a tutorial; there is an abundance of <a
href="http://www.google.com/search?q=sql+tutorial">tutorials online</a>.
Also, the  <a href="http://apophenia.sourceforge.net/gsl_stats.pdf">PDF
documentation</a> for Apophenia includes a chapter which discusses SQL
for statisticians.

The blog of Apophenia's author includes an <a
href="http://fluff.info/blog/arch/00000118.htm">entry</a> about
complementarities between SQL and matrix manipulation packages.

*/

/** \page graphing What about graphing?

Portable graphing tools are supremely difficult to implement. The closest
thing to a truly portable setup is <a href="http://www.gnuplot.info">Gnuplot</a>,
but you may have something on your system which you prefer.

The \ref apop_output.c file includes a few functions to interface with gnuplot directly.

But every system worth its silicon will take input from a text
file. Therefore, it is easy to achieve harmony with Apophenia: do
the data management and data crunching in C, then use \ref apop_print
<tt>apop_print</tt> to dump your output to a text file, and graph away.

 In the future, Apophenia may include functions which call Gnuplot
directly. [It already includes one, which is so rudimentary it's not
even documented (see the <tt>linear_algebra.h</tt> file in the source code
if you're curious).]

Oh, but have a look at \ref apop_plot_line_and_scatter.
*/

/** \defgroup global_vars The global variables. */
