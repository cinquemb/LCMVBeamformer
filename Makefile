#CFLAGS= -Wall -Os -O2 -std=gnu++11 -pedantic -I/usr/local/include \
#$(HOME)/gcc/gcc-5.2.0/objdir/x86_64-apple-darwin13.4.0/libgcc/addtf3.o \
#$(HOME)/gcc/gcc-5.2.0/objdir/x86_64-apple-darwin13.4.0/libgcc/divtf3.o \
#$(HOME)/gcc/gcc-5.2.0/objdir/x86_64-apple-darwin13.4.0/libgcc/eqtf2.o \
#$(HOME)/gcc/gcc-5.2.0/objdir/x86_64-apple-darwin13.4.0/libgcc/floatditf.o \
#$(HOME)/gcc/gcc-5.2.0/objdir/x86_64-apple-darwin13.4.0/libgcc/floatunditf.o \
#$(HOME)/gcc/gcc-5.2.0/objdir/x86_64-apple-darwin13.4.0/libgcc/getf2.o \
#$(HOME)/gcc/gcc-5.2.0/objdir/x86_64-apple-darwin13.4.0/libgcc/multf3.o \
#$(HOME)/gcc/gcc-5.2.0/objdir/x86_64-apple-darwin13.4.0/libgcc/subtf3.o \
#$(HOME)/gcc/gcc-5.2.0/objdir/x86_64-apple-darwin13.4.0/libgcc/sfp-exceptions.o \
#$(HOME)/gcc/gcc-5.2.0/objdir/x86_64-apple-darwin13.4.0/libgcc/letf2.o \
#/usr/local/lib/libquadmath.a \
#/usr/local/lib/libgfortran.a \
#$(HOME)/arpack/ARPACK/libarpack_OSX.a \
#/usr/local/lib/libsuperlu.a \
#/usr/local/lib/libopenblas.a \
#/usr/local/lib/libarmadillo.a \

#/usr/local/lib/libjsoncpp.a

CFLAGS= -Wall -Os -O3 -std=gnu++11 -pedantic  /usr/local/lib/libarmadillo.dylib
JSONCPP=`pkg-config --cflags --libs jsoncpp`
CPP= ccache clang++

all:
	$(CPP) forward_solution.cpp -o forward_solution $(CFLAGS)
	$(CPP) generate_covariance_matrix.cpp -o generate_covariance_matrix $(CFLAGS)
	$(CPP) generate_observation_matrix_from_running_moments.cpp -o gen_obs $(CFLAGS) $(JSONCPP)

forward:
	$(CPP) forward_solution.cpp -o forward_solution $(CFLAGS)

filter:
	$(CPP) generate_covariance_matrix.cpp -o generate_covariance_matrix $(CFLAGS)

obs:
	$(CPP) generate_observation_matrix_from_running_moments.cpp -o gen_obs $(CFLAGS) $(JSONCPP)

debug:
	$(CPP) -g forward_solution.cpp -o forward_solution $(CFLAGS)
	$(CPP) -g generate_covariance_matrix.cpp -o generate_covariance_matrix $(CFLAGS)
	$(CPP) generate_observation_matrix_from_running_moments.cpp -o gen_obs $(CFLAGS) $(JSONCPP)

clean:
	rm forward_solution
	rm generate_covariance_matrix
	rm gen_obs
