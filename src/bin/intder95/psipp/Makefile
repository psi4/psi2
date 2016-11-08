# Makefile for standalone version of intder95's personal psipp

CODE = psipp

YACC = bison -y

CC = cc
CFLAGS = -O -I./

CSRC = \
parse.c error.c if.c include.c main.c name.c operator.c print.c scan.c 

BINOBJ = $(CSRC:%.c=%.o)

.PHONY: clean

default: $(CODE)

$(CODE): $(BINOBJ)
	$(CC) $(CFLAGS) $^ -o $(CODE)
	touch stamp-psipp

%.o: %.c
	$(CC) $(CFLAGS) -c $<

parse.c: parse.y
	$(YACC) -v -d $<
	mv y.tab.c parse.c
	mv y.tab.h parse.h

parse.h: parse.c

clean:
	/bin/rm -f $(CODE) *.o core parse.c parse.h y.output

