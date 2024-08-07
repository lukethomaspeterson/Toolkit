README

These routines build off of publicly-available code written by Angel Jorba
for computing normal forms of autonomous Hamiltonian systems around an 
equilibrium point. Codes have been modified to compute normal form around
L1 and L2 in the Earth-Moon circular restricted 3-body problem (CR3BP). 
Most of those routines are not included in this version (v1) of the code,
so as to simplify the applications for the user. The results from these scripts
are series expansions included as tables in the Toolkit/data directory. For our
purposes, a MATLAB wrapper has been written to show a multitude of useful examples
using local action-angle orbital elements in cislunar space--these are action-angle
variables in the normal form around L1 or L2. In order to run these examples, we still
must take advatange of the speed of coordinate transformation functions written in C/C++.
The code here is modified from the original routines from Angel Jorba, found at:

http://www.maia.ub.es/~angel/soft.html

Contents
--------

When the Toolkit.zip file is opened, it creates the following
subdirectories:

Toolkit
Toolkit/bin
Toolkit/data
Toolkit/src

Directory Toolkit/bin will contain, after compilation, the binary file for
the program that transforms coordinates to/from local orbital elements, called
'tcnf' for 'transform coordinates normal form.' The MATLAB wrapper calls this 
executable from the command line.

Directory Toolkit/data contains the input data files as well as the output
of the programs.

Directory Toolkit/L1 and Toolkit/L2 contain data files needed for the examples,
and is used as a place to store data files generated from running the MATLAB wrapper.

Directory Toolkit/src contains the source code as well as the makefile.


Building the examples
---------------------

To compile the examples, go to the directory Toolkit/src and edit the
makefile. In the first lines, you can set the name of your C and C++
compilers, as well as the corresponding flags to use during
compilation and linking. The ones I've left there correspond to the
GNU C/C++ compiler (they are called gcc and g++). For instance, these
options should work smoothly if you are working on a Linux PC (all the
linux distributions I know come with those compilers installed), though
I have set it up to run on a Windows machine. If you are on a system 
without gcc/g++ installed, you should put the corresponding values at 
the beginning of the makefile. Remember that some compilers need a 
special flag to deal with ANSI C source code.

In the makefile there are a few parameters to tell the programs some
default directories for several purposes. More concretely, there is
the definition of the BIN directory (where the binaries of the
programs will be left) and the DATA directory, that is the default
place used for the programs to look for and write the several
input/output files needed. Note that I've used relative pathnames, but
it is probably better to use absolute ones. I use the default values
I've left there, starting the programs either from the bin directory
(typing, for instance, './nf') or from the data directory (typing then
'../bin/nf').

Once the makefile has been set up, you can compile the programs. You
can type 'make all'. The binaries are left in the directory Toolkit/bin
(unless you have changed this in the makefile).

The name of the different programs is the following:

 1.- nf: computation of the normal form around Lagrange point.

 2.- cvnf: computation of the expansions of the changes of variables
     for the normal form computed by nf.

 3.- tcnf: transformation of a set of points from normal form
     coordinates to synodical ones, and viceversa.

Note first that (1) and (2) have been omitted from this version, v1, of
the package; however, the outputs for Earth-Moon L1 & L2 have been included
in the Toolkit/data directory. Also note that you can also type 'make clean', 
that will erase all the object files (.o) in the Toolkit/src directory and
all the files in the Toolkit/bin directory. The output files left in the data
directory have to be erased manually.

In the next section there are some hints about how to run these
programs.


Running the examples
--------------------

Here we will explain the input taken by these programs, as well as the
output produced. The output of the program to the screen is only to
inform you about what it is doing at each moment. Since this is not
necessary for understanding neither the inner working of the programs
nor the results (that are stored in files), we will not give any
explanation. For details, look into the source code. We only recall
that the expansions computed by these programs are asymptotic series,
and that they are only accurate in a small neighbourhood of the
corresponding equilibrium point.

 1.- nf. it has two parameters. the first one is the mass parameter of
    the rtbp and it is given in the source code (it is the value
    assigned to the variable called 'mu' at the beginning of the file
    'main-nf.cc'). The second parameter is the degree ('n') of the
    power expansion of the Hamiltonian that we put in normal form
    (hence, the final normal form will contain terms up to degree
    n/2 in the final action variables). The program asks for this
    value at the begining of the execution.

    The program writes to the screen some information about what it is
    doing. The important part of the output is left in several files,
    in the directory Toolkit/data:

    nf.ctl: ascii file with the two parameters used in the present
       run. In the first line there is the degree used for the
       expansions and, in the second line, the value of the mass
       parameter.

   nf.cvl: ascii file with the (real) change of variables that puts
      the linearized vectorfield around L123 in real normal form. This
      file is needed by program tcnf.

   nf.gen: binary file with the several generating functions used to
      put the Hamiltonian in normal form. This file is used by program
      cvnf.

   nf.res: ascii file with the final normal form. This file is used by
      program ninf.

 2.- cvnf. It looks for the file Toolkit/data/nf.gen (see above). It asks
    (interactively) for the degree up to which the changes are desired
    (of course, this value must be less than or equal to the maximum
    degree of the generating function stored in nf.gen. Moreover, the
    program also asks for which change we want to compute, the direct
    one (from normal form to initial coordinates), the inverse one, or
    both).

    The output is stored in several files:

    cvnf.[1-6]: binary files with the power expansion of the direct
       change of variables. The number between 1 and 6 in the
       extension field of the name of the file refers to the
       coordinate (1 is x, 2 is px, and so on).

    cvnfi.[1-6]: binary files with the power expnasion of the inverse
       change of variables. See the remarks for files cvnf.[1-6].

 3.- tcnf. It needs three values in the command line. The first one
    contains phase space points, in the same format as the file
    ninf.res (see below). The second parameter is the name of the
    output file. It will contain the transformed points in the same
    format as the input file. The third parameter is a flag with two
    possible values: if it is 1, the transformation goes from normal
    form coordinates to initial ones; if it is -1 performs the reverse 
    transformation (we note that, in both cases, the first file is the 
    input file and the second one is the output file).

    The program needs the files nf.ctl, nf.cvl and cvnf.[1-6] (see
    above). All of them are assumed to be in the directory Toolkit/data.
    If the flag is 1, the program also writes a file called 'rtbp.mu',
    that contain the mass parameter corresponding to this data. This
    is done in order to pass this value to program rtbp, that performs
    a numerical integration of the restricted three body problem (see
    below).

Boulder, CO, April 30th, 2024.

E-mail: Luke.Peterson@colorado.edu
Web: https://sites.google.com/view/lukethomaspeterson/
