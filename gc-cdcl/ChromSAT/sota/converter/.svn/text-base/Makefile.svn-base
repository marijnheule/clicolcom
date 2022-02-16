SOURCES =	converter.c main.c
HEADERS = 	converter.h

CC =		gcc
DEBUG_FLAGS = -g -ansi -pedantic -Wall
CFLAGS =	$(DEBUG_FLAGS)

TARGET =	converter
OBJS =		$(SOURCES:%.c=%.o)


all: $(OBJS)
	$(CC)  $(CFLAGS) -o  $(TARGET) $(OBJS)

clean:
	rm -f $(OBJS) $(TARGET) *.cnf
	
zip::
	@tar cfz generator.tgz Makefile $(SOURCES) $(HEADERS)