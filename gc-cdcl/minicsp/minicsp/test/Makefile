MTL       = ../mtl
CORE      = ../core
CHDRS     = $(wildcard *.hpp) $(wildcard $(MTL)/*.h)
EXEC      = test
CFLAGS    = -I../.. -Wall -ffloat-store -std=c++17
LFLAGS    = -lz

COREOBJS  = $(CORE)/solver.o $(CORE)/cons.o $(CORE)/setcons.o $(CORE)/nvalue.o
CSRCS     = $(wildcard *.cpp)
COBJS     = $(addsuffix .o, $(basename $(CSRCS))) $(COREOBJS)

include ../mtl/template.mk
