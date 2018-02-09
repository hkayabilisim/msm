CC=gcc
CFLAGS=-ansi -Wall  -O3
HEADERS=msmLibrary.h ewaldHernquist.h ewaldClassical.h
OBJS=msmLibrary.o ewaldHernquist.o ewaldClassical.o
EXE=ewaldHernquist ewaldClassical
LIB=-lm

compile: $(EXE)

all: compile

ewaldHernquist: ewaldHernquist.o msmLibrary.o
	$(CC) -o ewaldHernquist msmLibrary.o ewaldHernquist.o $(LIB)

ewaldClassical: ewaldClassical.o msmLibrary.o
	$(CC) -o ewaldClassical msmLibrary.o ewaldClassical.o $(LIB)

msmLibrary.o: msmLibrary.c
	$(CC) $(CFLAGS) -I.  -c -o $@ $<

ewaldHernquist.o: ewaldHernquist.c
	$(CC) $(CFLAGS) -I.  -c -o $@ $<

ewaldClassical.o: ewaldClassical.c
	$(CC) $(CFLAGS) -I.  -c -o $@ $<

test: compile
	bash test.sh

clean:
	rm -f $(OBJS) $(EXE) *.gch a.out
