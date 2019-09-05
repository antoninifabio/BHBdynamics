CMPLR =gfortran 
FFLAGS =-O2 
LFLAGS =

.f.o:
	$(CMPLR) -c $(FFLAGS) $<

SRCE= mod.f BHBdynamics.f cluster.f quadpack.f read.f functions.f 

OBJT = $(SRCE:.f=.o)

binaryBH: $(OBJT) $(LFLAGS)
	$(CMPLR)   $(FFLAGS) $(OBJT) -o BHBdynamics



clean:
	rm BHBdynamics *o *mod
