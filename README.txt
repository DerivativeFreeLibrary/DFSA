-----------------------------------------------------------
 How to use the derivative-free simulated annealing global
 optimizer DFSA for constrained global optimization 
 problems
-----------------------------------------------------------
 The package provides a FORTRAN90 version of the code.

0- Gunzip and untar the archive in a folder on your computer by
   issuing in a directory of your choice (ex. curdir) the command

   $> tar -xvf DFSA.tar.gz

2- Edit file curdir/prob.f90 to define your own objective function,
   constraint functions and bounds on the variables.
   In particular, modify the subroutines 
   SETDIM    : which sets problem dimensions
   functinit : which sets upper and lower bounds on the variables plus
               a name for the problem and a know global value, if any
   funob     : which defines the objective function
   fconstreq : which defines the equality consraint function values
   fconstrin : which defines the inequality consraint function values

   N.B. equality constraints MUST BE of the form ce(x)  = 0
      inequality constraints MUST BE of the form ci(x) <= 0

   N.B. the code is a probabilistic one. The main program in file MainGlob.f90
        calls the algorithm subroutine giving a seed for the random sequence generator.
        Different runs can be obtained ba varying the seed.
 
2- At command prompt in curdir execute 

     $> make
 
   which will create the executables 'ddfsa'

4- execute

     $> ./dfsa_exe

   The program will output and create a file named fort.3 which contains:

problem name    &  n &  mi &  me & time &  bestobj &  feasibility &  obj &  bestfeas &   nf &  nloc

	where:

	- n is the number of variables
	- mi is the number of ineq. constraints
	- me is the number of eq. constraints
	- time is the cpu time in seconds
	- bestobj = f(xstar) is the best obj. function value obtained on a point xstar with violfeas(xstar) <= 1.e-4
	- feasibility = violfeas(xstar)
	- obj = f(xbar) is the obj. function value obtained on the point with lowest feasibility violation
	- bestfeas = violfeas(xbar)
 
