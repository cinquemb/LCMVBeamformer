
CFLAGS= -Wall -Os -O3 -std=gnu++11 -pedantic  /usr/local/lib/libarmadillo.dylib

all:
	ccache clang++ forward_solution.cpp -o forward_solution $(CFLAGS)

debug:
	ccache clang++ -g forward_solution.cpp -o forward_solution $(CFLAGS)

clean:
	rm forward_solution