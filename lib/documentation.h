/** \mainpage The Apophenia documentation

\section Prerequisites 
 \li \ref intro "Intro": The motivation for the package.
 \li \ref setup "Setup":  Installing GCC, the GSL, SQLite, and Apophenia itself.

\section dostats Doing statistics 
 \li \ref basic_stats "Basic statistics": Mean, variance, percentiles, &c.
 \li \ref regression  
 \li \ref mle "Maximum likelihood estimation": estimators requiring a search for the maximum of a likelihood function.
 \li \ref linear_algebra "Linear Algebra": determinants, projections, numerical gradients, &c. Some convenience functions to display matrices and vectors.
 \li \ref models : How to write down a model and estimate its parameters.
 \li \ref ttest "Some tests"

\section shuntdata Shunting data 
 \li \ref db "Database utilities": an easy front end to SQLite
 \li \ref conversions
 \li \ref output "Output functions": Summarize data and print tables to the screen or a file.
 \li \ref convenience_fns "Convenience functions": a few utilities to make life with the GSL a little easier.
 \li \ref command_line "Command-line utilities": a still easier front end to SQLite, and doing stats via Perl/Python/shell script.

\section speak Speaking the languages
 \li \ref c "C": For statisticians who don't know/are rusty with the C
 programming language. Those who know their C may want to check the
 Apophenia-specific notes in the \ref usagenotes "usage notes" section.
 \li \ref sql About the syntax for querying databases.
 \li \ref types "Types": New structures defined by Apophenia

 \ref admin
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
linear models? Apophenia facilitates writing \ref mle "likelihood functions" which
can be as crazy as the world you are modeling, and then fits them with
one function call to \ref apop_maximum_likelihood. Comparing competing
models is easy as well.
\li <b>Passing queries directly to an SQL engine.</b> if you have a gigabyte of data,
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
The Apophenia manual has a full tutorial for C, available as a separate document <a href="http://ben.klemens.org/pdfs/c_crash.pdf">here</a>.
The tutorial is more conceptual and aimed at clearing up the most common failures of understanding about C. More nuts-and-bolts
tutorials are 
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
If your text editor supports syntax highlighting, there are eight types defined in the Apophenia and GSL headers which may be worth coloring.
For <tt>vim</tt>, for example. add the following two lines to <tt>/usr/share/vim/syntax/c.vim</tt>:
\verbatim
syn keyword     cType           gsl_matrix gsl_vector gsl_matrix_view gsl_vector_view
syn keyword     cType           apop_name apop_model apop_inventory apop_estimate 
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
Oh, but have a look at \ref apop_plot_line_and_scatter and \ref apop_plot_histogram.
*/

/** \defgroup global_vars The global variables. */

/** \defgroup models Models

A model is an equation (or system of equations) which rely on data and
have unknown parameters to be determined. [Notice that this definition
readily includes "non-parametric" models.] Much of statistical analysis
consists of writing down a model, estimating its parameters, and running
hypothesis tests to determine the confidence with which we can make
statements about those parameters.

Apophenia facilitates this via its \ref apop_model objects. Each object is
a model as described above, and includes a method named <tt>estimate</tt>
which takes in data and returns an \ref apop_estimate which includes
the parameter estimates and the characteristics one would need for
hypothesis testing.

The design of the objects hopes to make it as easy as possible for you,
dear reader, to write new models. For the most part, all you need to
do is write a log likelihood function, and \ref apop_maximum_likelihood
does the rest; see below.

[Many statistics packages include a model structure that describes only
linear models.  "Linear" models can include a wide range of nonlinear
features, but they are still a subset of measure zero within the class of
models as described above. Currently, Apophenia has no plans to include
a summary syntax for describing linear models; the reader who has a linear
model to be estimated via OLS or GLS is advised to instead manipulate
the data set to the appropriate form and then call \ref apop_OLS or
\ref apop_GLS.]

Frequently, a model is a probability distribution. The data is assumed
to have been drawn from a given distribution and the question is
only what distributional parameters best fit; e.g., assume the data
is Normally distributed and find the mean and variance.

In fact, one could even use the \ref apop_model to describe \ref nonstats_modeling
"non-statistical models" involving optimization subject to constraints.

The main function in the systems below is the apop_model.estimate
function. It takes in a model and data, and outputs an apop_estimate,
that includes the parameter estimates and the various auxiliary data
that one may need to test the estimates, such as the variance-covariance
matrix. For most users, the apop_model.estimate function will be all
one needs from a model. Just prep the data, select a model, and produce
an estimate:

\code
gsl_matrix 	*data 		= read_in_data();
apop_estimate 	*the_estimate 	= apop_probit(data, NULL, NULL);
apop_estimate_print(the_estimate);
\endcode

Because models are often distributions, and because it is not our place
to dictate what you will do with a model, the apop_model also includes a
number of additional functions that may be useful for additional analyis,
such as a likelihood function that could be used for ML estimation or
for estimating the Hessian, and a random number generator. Some effort
has been made to ensure that the prepackaged models include as many of
these auxiliary functions as possible; if you are writing your own,
there is no requirement that you provide all functions, and the \ref
apop_maximum_likelihood and \ref apop_numerical_hessian functions do a
good job of filling in blanks.


\section write_likelihoods Writing your own
Writing apop_model objects is easy:


\li Write a likelihood function. Its header will look like this:
\code
double apop_new_log_likelihood(const gsl_vector *beta, void *d)
\endcode 
where \c beta will be the parameters to be maximized, and \c
d is the fixed parameters---the data. In every case currently included
with Apophenia, \c d is a \c gsl_matrix, but you do not have to conform
to that. This function will return the value of the log likelihood function at the given parameters.
\li Is this a constrained optimization? See the \ref constraints "Constraints page" on how to set them.
\li Write the object. In your header file, include 
\code
apop_model apop_new_likelihood = {"The Me distribution", number_of_parameters, 
{       //what will apop_new_likelihood.estimate return?
        1,      //parameters 
        1,      //covariance
        1,      //confidence
        0,      //predicted
        0,      //residuals
        1,      //log_likelihood
        1       //names;
},          new_estimate,
            new_log_likelihood, 
            NULL,   //place dlog likelihood here.
            NULL,   //place constraint fn here.
            NULL    //place RNG here.
            };
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
void apop_new_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient)
\endcode 
where \c beta and \c d are fixed as above, and \c gradient is a \c gsl_vector with dimension matching \c beta. 
At the end of this function, you will have to assign the appropriate derivative to every element of the gradient vector:
\code
gsl_vector_set(gradient,0, d_a);
gsl_vector_set(gradient,1, d_b);
\endcode 
Now add the resulting dlog likelihood function to your object, by replacing the \c NULL labeled "place dlog likelihood here" with the name of your dlog likelihood function.
\li Send the code to the maintainer for inclusion in future versions of Apophenia.

\todo This page needs a sample model for the reader to cut 'n' paste.
*/

/** \page constraints Setting constraints

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
static double beta_zero_and_one_greater_than_x_constraint(gsl_vector *beta, void * d, gsl_vector *returned_beta){
double          limit0          = 0,
                limit1          = 0,
                tolerance       = 1e-3; // GSL_EPSILON_DOUBLE is also a popular choice, but sometimes fails.
double          beta0   = gsl_vector_get(beta, 0),
                beta1   = gsl_vector_get(beta, 1);
        if (beta0 > limit0 && beta1 > limit1)
                return 0;
        //else:
        gsl_vector_set(returned_beta, 0, GSL_MAX(limit0 + tolerance, beta0));   //create a valid return vector.
        gsl_vector_set(returned_beta, 1, GSL_MAX(limit1 + tolerance, beta1));
        return GSL_MAX(limit0 - beta0, 0) + GSL_MAX(limit1 - beta1, 0);         //return a penalty.
}
\endcode

  */

/** \page admin The admin page

Just a few page links:
\li <a href=todo.html>The to-do list</a>
\li <a href=bug.html>The known bug list</a>
\li <a href=modules.html>The documentation page list</a>
*/


/** \page nonstats_modeling Non-statistical modeling

The MLE functions of Apophenia are designed to maximize a function
subject to constraints---which sounds a lot like any of a variety of
other problems, especially in Economics. With little abuse of the package,
one could use it to solve non-statistical models involving maximization
subject to constraints. For example, here is how one could numerically
solve a utility maximization problem from Econ 101. 


\code
/*  Here is a complete program to maximize a Cobb-Douglas utility
  function subject to a budget constraint.

  This requires a few digressions from the MLE nomenclature:
  --The log likelihood is actually the utility function. Don't take logs,
  since that would mess up the marginal utilities.

  --The data is the fixed input to the MLE, which means the parameters:
  prices, budget info, utility parameters.

  --The betas in the MLE framework are the free parameters to be
  maximized, which in this case means the goods the consumer is choosing.

 Utility is U = x_1^\alpha * x_2^\beta. 
 The budget constraint dictates that P_1 x_1 + P_2 x_2 <= B.

The data vector looks like this:
0:  price0
1:  price1
2:  budget
3:  alpha
4:  beta

Most of the work is in writing down all the constraints, since the
function itself is trivial. Having written down the model, the estimation
is one function call, and calculating the marginal values is one more.
 */

#include <apophenia/headers.h>

apop_model econ_101;

static apop_estimate * econ101_estimate(gsl_matrix * data, apop_inventory *uses, void *parameters){
apop_estimate           *est;
apop_estimation_params  mle_params;
    mle_params.method       = 000;
    mle_params.starting_pt  = NULL;
    mle_params.step_size    = 1e-1;
    mle_params.tolerance    = 1e-8;
    mle_params.verbose      = 0;
    apop_inventory_filter(uses, econ_101.inventory_filter);
    est = apop_maximum_likelihood(parameters, uses, econ_101, mle_params);
    est->uses.covariance    = 0;    //MLE finds them, but it's meaningless.
    est->uses.confidence    = 0;    //MLE finds them, but it's meaningless.
    est->uses.log_likelihood  = 0;
    return est;
}

/* The constraint function, including three constraints: x0>0, x1>0, and the bundle is under budget.
 First, we check wether anything binds, and if not we return zero immediately.
 Both sets of constraints are handled in the same way: derive new
 values that are within bounds, and report how far you had to move.
*/
static double budget_constraint(gsl_vector *beta, void * d, gsl_vector *returned_beta){
gsl_vector  *budget = d;
double  price0      = gsl_vector_get(budget, 0),
        price1      = gsl_vector_get(budget, 1),
        cash        = gsl_vector_get(budget, 2),
        x0          = gsl_vector_get(beta, 0),
        x1          = gsl_vector_get(beta, 1),
        tolerance   = 1e-3,
        new_x0, new_x1, penalty;
    if ((x0 * price0 + x1 * price1<= cash) && (x0 > 0) && (x1 > 0))
        return 0;
    //else:
    if (x0 <= 0 || x1 <= 0){
        new_x0  = (new_x0 > 0)? x0 : tolerance;
        new_x1  = (new_x1 > 0)? x1 : tolerance;
        penalty = (fabs(new_x0 - x0) + fabs(new_x1 + x1));
    } else {
        new_x0  = GSL_MAX(0, cash - x1* price1);
        new_x1  = GSL_MIN(x1, cash - new_x0* price0);
        penalty = (GSL_MAX(0,x0 * price0 + x1 * price1- cash));
    }
    gsl_vector_set(returned_beta, 0, new_x0);
    gsl_vector_set(returned_beta, 1, new_x1);
    return penalty;
}

static double econ101_log_likelihood(const gsl_vector *beta, void *d){
gsl_vector  *params = d;
double      bb0     = gsl_vector_get(params, 3),
            bb1     = gsl_vector_get(params, 4),
            qty0    = gsl_vector_get(beta, 0),
            qty1    = gsl_vector_get(beta, 1);
    return pow(qty0, bb0) * pow(qty1, bb1);
}    

apop_model econ_101 = {"Max Cobb-Douglass subject to a budget constraint", 2,  {
    1,    //parameters
    0,    //covariance
    0,    //confidence
    0,    //predicted
    0,    //residuals
    0,    //log_likelihood
    0    //names;
},         
    econ101_estimate, econ101_log_likelihood, NULL, NULL, budget_constraint, NULL};

int main(){
double          param_array[]   =  {1, 3, 38.4, 0.4, 0.6};//see header.
gsl_vector      *params         = gsl_vector_alloc(5);
gsl_vector      *marginals      = gsl_vector_alloc(2);
apop_estimate   *e;
double          x1, x2;
    apop_convert_array_to_vector(param_array,&params,5);
    e   = econ_101.estimate(NULL, NULL, params);
    printf("The optimal quantities:\n");
    apop_estimate_print(e);

    x2  = param_array[2]/(param_array[3]/param_array[4]+1)/param_array[1];
    x1  = (param_array[2] - param_array[1] * x2)/param_array[0];
    printf("\nAnalytically, these should be:\n %g\n %g\n\n", x1, x2);

    apop_fn_for_derivative  = econ_101.log_likelihood;
    apop_numerical_gradient(e->parameters, params, marginals);
    printf("The marginal values:\n");
    apop_vector_print(marginals, "\n", NULL);
    printf("\nAnalytically, these should be:\n %g\n %g\n", param_array[3]*pow(x1,param_array[3]-1)*pow(x2,param_array[4])
                                                        , param_array[4]*pow(x2,param_array[4]-1)*pow(x1,param_array[3]));

    return 0;
}
*/
