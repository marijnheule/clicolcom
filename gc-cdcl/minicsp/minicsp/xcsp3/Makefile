MTL       = ../mtl
CORE      = ../core
CHDRS     = $(wildcard *.hpp) $(wildcard *.hh) $(wildcard $(MTL)/*.h)
EXEC      = minicsp-xcsp3
CFLAGS    = -I../../../XCSP3-CPP-Parser/include -I../.. -I/usr/include/libxml2 -I/usr/local/include -I/usr/local/opt/libxml2/include/libxml2/  -Wall  -std=c++1z #-ffloat-store
LFLAGS    = -L../../../XCSP3-CPP-Parser/lib -lz -lxml2 $(LIBS) -lxcsp3parser -llzma -lz -licuuc -pthread -ldl -licudata
LIBSTDCXXDEBUG = 


COREOBJS  = $(CORE)/solver.o $(CORE)/cons.o $(CORE)/setcons.o \
	    $(CORE)/cmdline.o $(CORE)/utils.o $(CORE)/nvalue.o
CSRCS     = $(wildcard *.cpp)
COBJS     = $(addsuffix .o, $(basename $(CSRCS))) $(COREOBJS)

include ../mtl/template.mk
