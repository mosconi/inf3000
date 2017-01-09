CHK_CXXFLAGS=  -Isrc/checker -fPIC -fpic -g -O0 -std=c++11
CFLAGS=-Iinclude -fPIC -fpic -g -Wall -Werror -pedantic -std=c99 -O0 -Wextra

CHK_HDRS = src/checker/solution_checker.h
CHK_SRCS = src/checker/solution_checker.cc src/checker/solution_checker_run.cc

.SUFFIXES: .c .o .h

.c.o: ${HDRS}
	${CC} -o $@ -c $< ${CFLAGS}

checker: ${CHK_SRCS} ${CHK_HDRS}
	${CXX} -o $@ ${CHK_CXXFLAGS} ${CHK_SRCS}

.PHONY: clean test

test: roadeftest
	./roadeftest

clean:
	rm -f src/*.o
	rm -f checker
