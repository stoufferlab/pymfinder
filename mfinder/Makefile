# 
CC=gcc
# 
CFLAGS=-fpic -c -Wall -O3 -D UNIX -I/usr/include/python2.7/
#
LDFLAGS=-lm -O3 -shared
#
LD=ld
#

#
SWIG=swig
#
SWIGFLAGS=-python

MF_FILES=clustering globals grassberger hash list mat metropolis motif_ids \
         output permutation prob random results role stubs switches
WRAPPER_FILES=wrapper
SWIG_FILES=mfinder

MF_OBJECTS=$(MF_FILES:=.o)
SWIG_OBJECTS=$(SWIG_FILES:=.o) $(WRAPPER_FILES:=.o)

LIBRARIES=_mfinder.so

all: $(LIBRARIES)

_mfinder.so: $(MF_OBJECTS) $(SWIG_OBJECTS)
	$(LD) $(LDFLAGS) -o _mfinder.so $(MF_OBJECTS) $(SWIG_OBJECTS)

# clear out all suffixes
.SUFFIXES:
# clear out all suffixes
.SUFFIXES: .o .c .i

.c.o:
	$(CC) $(CFLAGS) $<

.i.c:
	$(SWIG) $(SWIGFLAGS) -o $@ $<

clean:
	rm -rf mfinder.py $(MF_OBJECTS) $(SWIG_OBJECTS) $(LIBRARIES)
