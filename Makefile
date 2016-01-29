CFLAGS= -Wall -Os -O2 -std=gnu++11 -pedantic  /usr/local/lib/libarmadillo.dylib
CPP= ccache clang++

all:
	$(CPP) forward_solution.cpp -o forward_solution $(CFLAGS)
	$(CPP) generate_covariance_matrix.cpp -o generate_covariance_matrix $(CFLAGS)
	$(CPP) generate_observation_matrix_from_running_moments.cpp -o gen_obs $(CFLAGS)

forward:
	$(CPP) forward_solution.cpp -o forward_solution $(CFLAGS)

filter:
	$(CPP) generate_covariance_matrix.cpp -o generate_covariance_matrix $(CFLAGS)

obs:
	$(CPP) generate_observation_matrix_from_running_moments.cpp -o gen_obs $
(CFLAGS)

debug:
	$(CPP) -g forward_solution.cpp -o forward_solution $(CFLAGS)
	$(CPP) -g generate_covariance_matrix.cpp -o generate_covariance_matrix $(CFLAGS)
	$(CPP) generate_observation_matrix_from_running_moments.cpp -o gen_obs $
(CFLAGS)

clean:
	rm forward_solution
	rm generate_covariance_matrix
	rm gen_obs
