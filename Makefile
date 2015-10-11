
CFLAGS= -Wall -Os -O3 -std=gnu++11 -pedantic  /usr/local/lib/libarmadillo.dylib

all:
	ccache clang++ forward_solution.cpp -o forward_solution $(CFLAGS)
	ccache clang++ generate_covariance_matrix.cpp -o generate_covariance_matrix $(CFLAGS)

filter:
	ccache clang++ generate_covariance_matrix.cpp -o generate_covariance_matrix $(CFLAGS)

debug:
	ccache clang++ -g forward_solution.cpp -o forward_solution $(CFLAGS)
	ccache clang++ -g generate_covariance_matrix.cpp -o generate_covariance_matrix $(CFLAGS)

clean:
	rm forward_solution
	rm generate_covariance_matrix