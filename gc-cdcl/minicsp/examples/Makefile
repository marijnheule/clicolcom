MTL       = ../mtl
CORE      = ../core
CHDRS     = $(wildcard *.hpp) $(wildcard $(MTL)/*.h)
CFLAGS    = -I.. -Wall -ffloat-store
LFLAGS    = -lz

COREOBJS  = $(CORE)/solver.o $(CORE)/cons.o $(CORE)/cmdline.o \
	    $(CORE)/utils.o
CSRCS     = $(wildcard *.cpp)
COBJS     = $(COREOBJS)

EXECS     = $(basename $(CSRCS))
EXECSp    = $(addsuffix _profile, $(basename $(CSRCS)))
EXECSpg   = $(addsuffix _pgen, $(basename $(CSRCS)))
EXECSpu   = $(addsuffix _puse, $(basename $(CSRCS)))
EXECSd    = $(addsuffix _debug, $(basename $(CSRCS)))
EXECSr    = $(addsuffix _release, $(basename $(CSRCS)))
EXECSrs   = $(addsuffix _static, $(basename $(CSRCS)))
EXECSc    = $(addsuffix _cov, $(basename $(CSRCS)))

PCOBJS     = $(addsuffix p,  $(COBJS))
PGCOBJS    = $(addsuffix pg,  $(COBJS))
PUCOBJS    = $(addsuffix pu,  $(COBJS))
DCOBJS     = $(addsuffix d,  $(COBJS))
RCOBJS     = $(addsuffix r,  $(COBJS))
CCOBJS     = $(addsuffix c,  $(COBJS))

CXX       ?= g++
CFLAGS    ?= -Wall
LFLAGS    ?= -Wall

COPTIMIZE ?= -O3 --param inline-unit-growth=60

#GOOGLE_PROFILER ?= -lprofiler
GOOGLE_PROFILER ?= 

.PHONY : s p d r rs c clean 

s:         $(EXECS)
p:         $(EXECSp)
pg:        $(EXECSpg)
pu:        $(EXECSpu)
d:         $(EXECSd)
r:         $(EXECSr)
rs:        $(EXECSrs)
c:         $(EXECSc)

## Compile options
%.o:			CFLAGS +=$(COPTIMIZE) -ggdb -D DEBUG
%.opg:			CFLAGS +=$(COPTIMIZE) -fprofile-generate -ggdb -D DEBUG
%.opu:			CFLAGS +=$(COPTIMIZE) -fprofile-use -ggdb -D DEBUG
%.op:			CFLAGS +=$(COPTIMIZE) -pg -ggdb -D NDEBUG
%.od:			CFLAGS +=-O0 -ggdb -D DEBUG -D INVARIANTS -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
%.or:			CFLAGS +=$(COPTIMIZE) -D NDEBUG
%.oc:                   CFLAGS +=-O0 -fprofile-arcs -ftest-coverage -ggdb -D DEBUG

## Link options
$(EXECS):		LFLAGS := -ggdb $(LFLAGS) $(GOOGLE_PROFILER)
$(EXECSpg):		LFLAGS := -ggdb $(LFLAGS) $(GOOGLE_PROFILER) -fprofile-generate
$(EXECSpu):		LFLAGS := -ggdb $(LFLAGS) $(GOOGLE_PROFILER)
$(EXECSp):	        LFLAGS := -ggdb -pg $(LFLAGS)
$(EXECSd):		LFLAGS := -ggdb $(LFLAGS)
$(EXECSr):	        LFLAGS := $(LFLAGS)
$(EXECSrs):		LFLAGS := --static $(LFLAGS)
$(EXECSc):              LFLAGS := -ggdb -fprofile-arcs -ftest-coverage $(LFLAGS)

# I know nothing about makefiles, i have no idea why these are needed
$(EXECS): $(COBJS)
$(EXECSp): $(PCOBJS)
$(EXECSpg): $(PGCOBJS)
$(EXECSpu): $(PUCOBJS)
$(EXECSd): $(DCOBJS)
$(EXECSr): $(RCOBJS)
$(EXECSrs): $(RCOBJS)
$(EXECSc): $(CCOBJS)

## Dependencies
%:              %.o $(COBJS)
	@echo Linking: "$@ ( $^ )"
	@$(CXX) $^ $(LFLAGS) -o $@

%_profile:	%.op $(PCOBJS)
	@echo Linking: "$@ ( $^ )"
	@$(CXX) $^ $(LFLAGS) -o $@

%_pgen:              %.o $(PGCOBJS)
	@echo Linking: "$@ ( $^ )"
	@$(CXX) $^ $(LFLAGS) -o $@

%_puse:              %.o $(PUCOBJS)
	@echo Linking: "$@ ( $^ )"
	@$(CXX) $^ $(LFLAGS) -o $@

%_debug:	%.od $(DCOBJS)
	@echo Linking: "$@ ( $^ )"
	@$(CXX) $^ $(LFLAGS) -o $@

%_release:	%.or $(RCOBJS)
	@echo Linking: "$@ ( $^ )"
	@$(CXX) $^ $(LFLAGS) -o $@

%_static:	%.or $(RCOBJS)
	@echo Linking: "$@ ( $^ )"
	@$(CXX) $^ $(LFLAGS) -o $@

%_cov:          %.oc $(CCOBJS)
	@echo Linking: "$@ ( $^ )"
	@$(CXX) $^ $(LFLAGS) -o $@

## Build rule
%.o %.op %.opg %.opu %.od %.or %.oc:	%.cpp
	@echo Compiling: "$@ ( $< )"
	@$(CXX) $(CFLAGS) -c -o $@ $<


clean:
	@rm -f $(EXECS) $(EXECSp) $(EXECSd) $(EXECSr) $(EXECSrs) $(EXECSc) \
	$(COBJS) $(PCOBJS) $(DCOBJS) $(RCOBJS) $(CCOBJS) \
	*gcno *.gcda *.gcov

## Make dependencies
depend.mk: $(CSRCS)
	@echo Making dependencies ...
	@$(CXX) $(CFLAGS) -MM $(CSRCS) > depend.mk
	@cp depend.mk /tmp/depend.mk.tmp
	@sed "s/o:/op:/" /tmp/depend.mk.tmp >> depend.mk
	@sed "s/o:/od:/" /tmp/depend.mk.tmp >> depend.mk
	@sed "s/o:/or:/" /tmp/depend.mk.tmp >> depend.mk
	@rm /tmp/depend.mk.tmp

-include depend.mk
