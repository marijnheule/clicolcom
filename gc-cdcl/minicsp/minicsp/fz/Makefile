MTL       = ../mtl
CORE      = ../core
CHDRS     = $(wildcard *.hpp) $(wildcard $(MTL)/*.h)
EXEC      = minicsp-fz
CFLAGS    = -I../.. -Wall -ffloat-store -std=c++17
LFLAGS    = 

COREOBJS  = $(CORE)/solver.o $(CORE)/cons.o $(CORE)/setcons.o \
	    $(CORE)/cmdline.o $(CORE)/utils.o
CSRCS     = $(wildcard *.cpp) lexer.yy.cpp parser.tab.cpp
COBJS     = $(addsuffix .o, $(basename $(CSRCS))) $(COREOBJS)

include ../mtl/template.mk

lexer.yy.cpp: lexer.lxx
	flex -o lexer.yy.cpp lexer.lxx

parser.tab.cpp parser.tab.hh: parser.yxx
	bison -t -o parser.tab.cpp -d $<
