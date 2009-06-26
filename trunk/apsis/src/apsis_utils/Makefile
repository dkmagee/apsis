
# Mac OSX with scisoft
CFITSIO_LIB=/Applications/scisoft/i386/Packages/cfitsio/lib
CFITSIO_INC=/Applications/scisoft/i386/Packages/cfitsio/include

# Linux
#CFITSIO_LIB=/usr/lib
#CFITSIO_INC=/usr/include/cfitsio


CC=gcc
#CFLAGS=-DLINUX -DBSWAP -fPIC -DPC_LINUX -g
CFLAGS=
INSTALL_DIR=../../bin

all: simplematch superalign medianfilter maskreg findlines findpedest fillnoise whmedian

# .o : arrays.c matchup.c sort.c error.c mt19937ar.c
#	$(CC) $(CFLAGS) -c $< -o $@
# .o.c:
# 	(CC) $(CFLAGS) $@.c -o $@

clean:
	rm -rf *.o
	rm -rf *.dSYM
	
distclean:
	rm -rf simplematch superalign medianfilter maskreg findlines findpedest fillnoise whmedian
	rm -rf *.o
	rm -rf *.dSYM

simplematch: arrays.o matchup.o sort.o error.o simplealign.o mt19937ar.o simplematch.c
	$(CC) -g -o simplematch simplematch.c arrays.o matchup.o sort.o error.o simplealign.o mt19937ar.o -lm

superalign: arrays.o matchup.o sort.o error.o  simplealign.o mt19937ar.o superalign.c
	$(CC) -g -o superalign superalign.c arrays.o matchup.o sort.o error.o  simplealign.o mt19937ar.o -lm

medianfilter: arrays.o error.o  medianfilter.c
	$(CC) -g -o medianfilter medianfilter.c arrays.o error.o  -I${CFITSIO_INC} -lm -L${CFITSIO_LIB} -lcfitsio

maskreg: arrays.o error.o  maskreg.c
	$(CC) -g -o maskreg maskreg.c arrays.o error.o  -I${CFITSIO_INC} -lm -L${CFITSIO_LIB} -lcfitsio

findlines: arrays.o error.o  findlines.c
	$(CC) -g -o findlines findlines.c arrays.o error.o  mt19937ar.o sort.o -I${CFITSIO_INC} -L${CFITSIO_LIB} -lm -lcfitsio

findpedest: arrays.o error.o  findpedest.c
	$(CC) -g -o findpedest findpedest.c arrays.o error.o sort.o  -I${CFITSIO_INC} -lm -L${CFITSIO_LIB} -lcfitsio

fillnoise: arrays.o error.o  fillnoise.c
	$(CC) -g -o fillnoise fillnoise.c arrays.o error.o  mt19937ar.o -I${CFITSIO_INC} -lm -L${CFITSIO_LIB} -lcfitsio

whmedian: arrays.o whmedian.c
	$(CC) -g -o whmedian whmedian.c arrays.o error.o  mt19937ar.o -I${CFITSIO_INC} -lm -L${CFITSIO_LIB} -lcfitsio

distrmatch:
	tar -cvf simplematch.tar simplematch.c arrays.c matchup.c sort.c error.c nr.h arrays.h error.h makefile mt19937ar.c cholsl.c choldc.c nrutil.c simplealign.c superalign.h superalign.c

install:
	cp simplematch superalign medianfilter maskreg findlines findpedest fillnoise whmedian ${INSTALL_DIR}
