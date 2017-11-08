CC = gcc
CFLAGS = -O3 -Wall
LIBS = -lm -lgsl -lgslcblas

HEADERS = initfcc.h alpotential.h
OBJECTS = initfcc.o alpotential.o main.o
PROGRAM = go

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

