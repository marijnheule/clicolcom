PRG =  coloring
DIR = Coloring
LIBS = 
SRC = $(wildcard *.C)
OBJ = $(SRC:.C=.o)



FLAGS   = -O2 -g -Wall
CC = g++

%.o :	%.C
	$(CC) $(FLAGS) -c $<

all:	$(OBJ)
	$(CC) $(FLAGS) -o $(PRG) $(OBJ) $(LIBS)

tar:	clean
	cd ..;tar -czvf $(PRG).tar.gz $(DIR)

floppy:	tar
	mount /floppy; cp ../$(PRG).tar.gz /floppy/.; umount /floppy

cpfromfloppy:
	mount /floppy; cp /floppy/$(PRG).tar.gz ../.; umount /floppy

unfloppy:	cpfromfloppy untar

untar:	
	cd ..;tar -xzvf $(PRG).tar.gz

clean:
	rm -f *~ *.o $(PRG) *.bak core
depend:
	makedepend -- $(FLAGS) -- $(SRC)

# DO NOT DELETE
