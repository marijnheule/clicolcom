CCC = g++ -O3
# CCC = g++ -pg -g

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------

CCOPT = -m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD

CCLNFLAGS = -lm -lpthread -ldl
CLNFLAGS  = -lm -lpthread -ldl

# --- Local directories ---

SRC = src
INC = include
OUT = obj


# --- TARGETS ---
color: $(OUT)/color.o $(OUT)/readfile.o $(INC)/readfile.h
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o color $(OUT)/color.o $(OUT)/readfile.o $(CCLNFLAGS)

$(OUT)/color.o: $(SRC)/color.cpp $(INC)/readfile.h
	$(CCC) -c $(CCFLAGS) $(SRC)/color.cpp -o $(OUT)/color.o

$(OUT)/readfile.o: $(SRC)/readfile.cpp $(INC)/readfile.h
	$(CCC) -c $(CCFLAGS) $(SRC)/readfile.cpp -o $(OUT)/readfile.o

clean:
	rm -f $(OUT)/*.o
	rm -f color
