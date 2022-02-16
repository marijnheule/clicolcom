MTL       = ../mtl
CORE      = ../core
CHDRS     = $(wildcard *.hpp) $(wildcard *.hh) $(wildcard $(MTL)/*.h)
EXEC      = minicsp-xcsp
CFLAGS    = -I../.. -I/usr/include/libxml2 -I/usr/local/include -I/usr/local/opt/libxml2/include/libxml2/ -Wall -ffloat-store -std=c++1z
LFLAGS    = -lz -lxml2
LIBSTDCXXDEBUG = 


COREOBJS  = $(CORE)/solver.o $(CORE)/cons.o $(CORE)/setcons.o \
	    $(CORE)/cmdline.o $(CORE)/utils.o
CSRCS     = $(wildcard *.cpp)
COBJS     = $(addsuffix .o, $(basename $(CSRCS))) $(COREOBJS)

include ../mtl/template.mk
