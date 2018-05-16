CC=gcc
CFLAGS=-ansi -Wall -O3 -std=c99
HEADERS=msmLibrary.h ewaldHernquist.h ewaldClassical.h
OBJS=msmLibrary.o msmPeriodic.o ewaldHernquist.o ewaldClassical.o msmExpected.o msmSimple.o
EXE=msmPeriodic msmExpected ewaldHernquist ewaldClassical msmSimple
LIB=-lm

compile: $(EXE)

all: compile

ewaldHernquist: ewaldHernquist.o msmLibrary.o
	$(CC) -o ewaldHernquist msmLibrary.o ewaldHernquist.o $(LIB)

msmPeriodic: msmPeriodic.o msmLibrary.o
	$(CC) -o msmPeriodic msmPeriodic.o msmLibrary.o $(LIB)

msmExpected: msmExpected.o msmLibrary.o
	$(CC) -o msmExpected msmExpected.o msmLibrary.o $(LIB)

msmSimple: msmSimple.o
	$(CC) -o msmSimple msmSimple.o $(LIB)

ewaldClassical: ewaldClassical.o msmLibrary.o
	$(CC) -o ewaldClassical msmLibrary.o ewaldClassical.o $(LIB)

msmLibrary.o: msmLibrary.c
	$(CC) $(CFLAGS) -I.  -c -o $@ $<

ewaldHernquist.o: ewaldHernquist.c
	$(CC) $(CFLAGS) -I.  -c -o $@ $<

ewaldClassical.o: ewaldClassical.c
	$(CC) $(CFLAGS) -I.  -c -o $@ $<

msmPeriodic.o: msmPeriodic.c
	$(CC) $(CFLAGS) -I.  -c -o $@ $<

msmExpected.o: msmExpected.c
	$(CC) $(CFLAGS) -I.  -c -o $@ $<

test: compile
	bash test.sh

forcetest: compile
	bash forcetest.sh

clean:
	rm -f $(OBJS) $(EXE) *.gch a.out
