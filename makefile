LCLPTH = ./
#
# Define compiler and library flags.
CC =icc
CFLGS = -g -debug all -static -restrict  -O2 -ipo -xSSE4.2 -m64 -fp-model precise -fp-model source -I/apps/mgmt/appadmin/gsl-1.16/BUILD/include -I/apps/mgmt/appadmin/gsl-1.16/BUILD/include/gsl 
LDFLGS=   -vec-report2 -O2 -m64 -xSSE4.2 -fp-model precise -fp-model source -L/apps/mgmt/appadmin/gsl-1.16/BUILD/lib/ -lgsl -lgslcblas   -lm

#
# Define the object files.
#

OBJ =SHH.o RHF-F.o rk4.o


#
default: SHH
clean: ;rm -f $(OBJ); rm -f SHH; rm -f *.o; rm -f *.mod;
SHH: $(OBJ); $(CC) $(CFLGS) -o SHH $(OBJ) $(LDFLGS)
#
# Make the object files. Start with concrete ones.
SHH.o :$(LCLPATH)SHH.c;$(CC) $(CFLGS) -c $(LCLPTH)SHH.c
rk4.o :$(LCLPATH)rk4.c;$(CC) $(CFLGS) -c $(LCLPTH)rk4.c
RHF-F.o :$(LCLPTH)RHF-F.c; $(CC) $(CFLGS) -c $(LCLPTH)RHF-F.c
