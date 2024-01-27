RM=rm -f
FC=gfortran

EXE=dfsa_exe

OBJS=main.o DFSA.o DFL-dagonale-passi_strong.o halton.o prob.o

all: dfsa
	       
dfsa: $(OBJS)
	$(FC) -o $(EXE) $(OBJS)

.SUFFIXES : .f90 .o

.f90.o: $* ; $(FC) -c $*.f90

clean:
	$(RM) *.o
	$(RM) *.mod
	$(RM) fort.1
	$(RM) fort.3
	$(RM) $(EXE)
