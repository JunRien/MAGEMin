# CC=gcc
CC=clang 

UNAME_S := $(shell uname -s)
EXE_NAME:= MAGEMin

ifeq ($(UNAME_S),Darwin)
	CCFLAGS = -Wall -O3 -g -fPIC -Wno-unused-variable -Wno-unused-but-set-variable
	
	# These are the flags to be used if you followed the instructions on the webpage and installed MPICH and NLopt all through homebrew 
	# Note that we use the apple Accelerate framework to call Lapack routines
    #LIBS    = -lm -framework Accelerate /opt/homebrew/lib/libnlopt.dylib /opt/homebrew/lib/libmpi.dylib
    #INC     = -I/opt/homebrew/include 

	LIBS   = -lm -framework Accelerate /usr/local/opt/lapack/lib/liblapacke.dylib /usr/local/lib/libnlopt.dylib /usr/local/lib/libmpich.dylib
	INC    = -I/usr/local/opt/lapack/include -I/usr/local/include

	# This is for a Mac, where we employ the build-in CBLAS and LAPACKE from MacPorts
	# LIBS    = -lm -framework Accelerate /opt/local/lib/lapack/liblapacke.dylib /usr/local/lib/libnlopt.dylib ~/Software/mpich-3.3.2/mpich-install/lib/libmpich.dylib
	# INC     = -I/opt/local/include/lapack -I/usr/local/include -I/Users/kausb/Software/mpich-3.3.2/mpich-install/include
endif
ifeq ($(UNAME_S),Linux)
	 # This is for a Linux:	
	 # -g, -O3, normal vs optimized compilation (~ 2/3 times faster with -O3)
	 # -Wno-unused-variable -Wno-unused-but-set-variable -Wno-maybe-uninitialized -Wno-unused-result
	 #  -lllalloc
	CCFLAGS = -Wall -O3 -g -fPIC -Wno-unused-variable -Wno-unused-result -Wno-unused-function
	LIBS   += -lm -llapacke -lnlopt -g -L/usr/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi
	INC     = -I/usr/lib/x86_64-linux-gnu/openmpi/include/
	
	## RUN MAGEMIN ON PLUTON
	#  CCFLAGS = -Wall -O3 -g -shared -fPIC -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-result -Wno-unused-function
	#  LIBS   += -lm -L/local/home/nriel/lapack-3.10.0 -llapacke -llapack -lrefblas -lgfortran -L/local/home/nriel/nlopt_install/install/lib -lnlopt -g -L/opt/mpich3/lib -lmpi
	#  INC     = -I/opt/mpich3/include -I/local/home/nriel/lapack-3.10.0/LAPACKE/include -I/local/home/nriel/nlopt_install/install/include
endif
	EXE_NAME = MAGEMin

SOURCES=src/MAGEMin.c 					\
		src/toolkit.c					\
		src/io_function.c				\
		src/gem_function.c 				\
		src/gss_init_function.c			\
		src/gss_function.c				\
		src/NLopt_opt_function.c 		\
		src/objective_functions.c		\
		src/pp_min_function.c 			\
		src/ss_min_function.c 			\
		src/simplex_levelling.c 		\
		src/PGE_function.c 				\
		src/phase_update_function.c		\
		src/dump_function.c

OBJECTS=$(SOURCES:.c=.o)
 
.c.o:
	$(CC) $(CCFLAGS) -c $< -o $@ $(INC)
 
all: $(OBJECTS)
	$(CC)  -o $(EXE_NAME) $(OBJECTS) $(INC) $(LIBS) 
	rm src/*.o

lib: $(OBJECTS)
	$(CC) -shared -fPIC  -o libMAGEMin.dylib $(OBJECTS) $(INC) $(LIBS)
 
clean:
	rm -f src/*.o *.dylib MAGEMin
