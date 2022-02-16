###################################################################################################################################
MAINNAME	=	main
MAINCPP		=	graphe.o \
			dsatur_algo.o \
			bornes.o \
			main.o
MAINC		=	cliquer/reorder.o \
			cliquer/graph.o \
			cliquer/cliquer.o

MAINCPP2		=	$(MAINCPP:.o=.cpp)
MAINCC2		=	$(MAINCC:.o=.cc)
MAINC2		=	$(MAINC:.o=.c)

## Compiler:
CC=g++ -O3 -g
#CC=g++ -O3 

## Compiler flags:
#CFLAGS=-Wall -O3 -fomit-frame-pointer -funroll-loops
#CFLAGS= -O3 -fomit-frame-pointer -funroll-loops 
#LONGOPTS = -DENABLE_LONG_OPTIONS

all:main

cliquer/%.o:	cliquer/%.c
	@echo "-> compiling $@"
	$(CC) $(CFLAGS) $(LONGOPTS) -o $@ -c $<

%.o:	%.cpp
	@echo "-> compiling $@"
	$(CC) $(CFLAGS) $(LONGOPTS) -o $@ -c $<

main : $(MAINCPP) $(MAINCC) $(MAINC)
	@echo "-> compiling $@"
	$(CC) $(CCFLAGS) $(CCLNDIRS) $(LDFLAGS) $(INCDIR) -o main $(MAINCPP) $(MAINCC) $(MAINC) $(CCLNFLAGS)

clean :
	rm -f *~* *.o cliquer/*.o main

###################################################################################################################################
