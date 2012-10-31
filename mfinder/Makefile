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

MF_FILES=clustering grassberger hash list mat metropolis motif_ids \
         output permutation prob random results role stubs switches
SWIG_FILES=mfinder wrapper

MF_OBJECTS=$(MF_FILES:=.o)
SWIG_OBJECTS=$(SWIG_FILES:=.o) 

EXECUTABLES=mfinder
LIBRARIES=_mfinder.so

all: $(EXECUTABLES) $(LIBRARIES)

_mfinder.so: $(MF_OBJECTS) $(SWIG_OBJECTS)
	$(LD) $(LDFLAGS) -o _mfinder.so $(MF_OBJECTS) $(SWIG_OBJECTS)

$(EXECUTABLES): $(MF_OBJECTS) main.o
	$(CC) $(MF_OBJECTS) main.o -o $(EXECUTABLES) $(LDFLAGS)

# clear out all suffixes
.SUFFIXES:
# clear out all suffixes
.SUFFIXES: .o .c .i

.c.o:
	$(CC) $(CFLAGS) $<

.i.c:
	$(SWIG) $(SWIGFLAGS) -o $@ $<

clean:
	rm -rf main.o mfinder.py $(MF_OBJECTS) $(SWIG_OBJECTS) $(EXECUTABLES) $(LIBRARIES)
