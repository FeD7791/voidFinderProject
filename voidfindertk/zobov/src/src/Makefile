# Compiler choice:

# Gcc
CC = gcc
CFLAGS = -O2

# Mark's Sun
#CC = /usr/local/forte/SUNWspro/bin/cc
#CFLAGS  =   -xO4 -fsimple=2 -fsingle 

# IBM p690
#CC = xlc_r
#CFLAGS = -O2 -q64 -qsmp=omp -qreport

MLIBS	=   -lm

QLIB = -L../qhull/src/
QINC = -I../qhull/src/

###############

EXES = vozinit voz1b1 voztie jovoz boz jozov

all: $(EXES)

jovoz: jovoz.o findrtop.o
	$(CC) $(CFLAGS) -o jovoz jovoz.o findrtop.o $(MLIBS)

jozov: jozov.o findrtop.o
	$(CC) $(CFLAGS) -o jozov jozov.o findrtop.o $(MLIBS)

boz: boz.o readfiles.o findrtop.o
	$(CC) $(CFLAGS) -o boz boz.o readfiles.o findrtop.o $(MLIBS)

voz1b1: voz1b1.o readfiles.o vozutil.o voz.h
	$(CC)  -o voz1b1 $(CFLAGS) voz1b1.o readfiles.o vozutil.o -L. $(QLIB) -lqhull $(MLIBS) 

vozinit: vozinit.o readfiles.o voz.h
	$(CC)  -o vozinit $(CFLAGS) vozinit.o readfiles.o -L. $(MLIBS) 

voztie: voztie.o readfiles.o
	$(CC)  -o voztie $(CFLAGS) voztie.o readfiles.o

install: $(EXES)
	cp $(EXES) ../bin

.c.o:
	$(CC) -c -o $@ $(CFLAGS) $(QINC) $*.c

clean:
	rm -f $(EXES) *.o *~
