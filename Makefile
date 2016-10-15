CFLAGS=-Iinclude -fPIC -fpic -g -Wall -Werror -pedantic -std=c99 -O0 -Wextra
LDFLAGS=-L. -lroadef

OBJS=src/roadef.o 
HDRS=include/roadef.h include/prelude.h

HDRS+=include/machine.h
OBJS+=src/machine.o

HDRS+=include/resource.h
OBJS+=src/resource.o

HDRS+=include/service.h
OBJS+=src/service.o

HDRS+=include/process.h
OBJS+=src/process.o

HDRS+=include/balance.h
OBJS+=src/balance.o

HDRS+=include/instance.h
OBJS+=src/instance.o



.SUFFIXES: .c .o .h

.c.o: ${HDRS}
	${CC} -o $@ -c $< ${CFLAGS}

all: libroadef.a roadeftest

libroadef.a: ${OBJS}
	ar rcs $@ ${OBJS}

roadeftest: src/roadef_test.c libroadef.a ${HDRS}
	${CC} -o $@ src/roadef_test.c ${CFLAGS} ${LDFLAGS}

.PHONY: clean test

test: roadeftest
	./roadeftest

clean:
	rm -f src/*.o
	rm -f libroadef.a roadeftest
