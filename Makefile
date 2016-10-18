CXXFLAGS=-Iinclude -Iinclude/checker -fPIC -fpic -g -Wall -O0 -std=c++11
CFLAGS=-Iinclude -fPIC -fpic -g -Wall -Werror -pedantic -std=c99 -O0 -Wextra
LDFLAGS=-L. -lroadef

OBJS=src/roadef.o
HDRS=include/roadef.h include/prelude.h 

HDRS+=include/rss.h
OBJS+=src/rss.o

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

HDRS+=include/model.h
OBJS+=src/model.o

HDRS+=include/state.h
OBJS+=src/state.o

CHK_HDRS = include/checker/solution_checker.h
CHK_SRCS = src/checker/solution_checker.cc src/checker/solution_checker_run.cc

.SUFFIXES: .c .o .h

.c.o: ${HDRS}
	${CC} -o $@ -c $< ${CFLAGS}

all: libroadef.a roadeftest roadef

libroadef.a: ${OBJS}
	ar rcs $@ ${OBJS}

roadeftest: src/roadef_test.c libroadef.a ${HDRS}
	${CC} -o $@ src/roadef_test.c ${CFLAGS} ${LDFLAGS}

roadef: src/roadef.c libroadef.a ${HDRS}
	${CC} -o $@ src/roadef.c ${CFLAGS} ${LDFLAGS}

checker: ${CHK_SRCS} ${CHK_HDRS}
	${CXX} -o $@ ${CXXFLAGS} ${CHK_SRCS}

.PHONY: clean test

test: roadeftest
	./roadeftest

clean:
	rm -f src/*.o
	rm -f libroadef.a roadeftest checker roadef
