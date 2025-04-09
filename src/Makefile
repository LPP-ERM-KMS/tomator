FLAGS = -fopenmp -DEIGEN_NO_DEBUG -I../Data/Fortran -J../Data/Fortran -lgfortran -lpthread -Wall
#FLAGS = -fopenmp -DEIGEN_NO_DEBUG -IAUG4 -JAUG4 -lgfortran -lpthread -Wall ! old flags
#DEBUG_FLAGS = -fopenmp -g -DDEBUG -ggdb3 -DEIGEN_NO_DEBUG -IAUG4 -JAUG4 -lgfortran -lpthread -Wall

# Specify directories
SRC_DIRS = . Funcs Vars

# Using wildcard and patsubst functions to find source files and create object file names
CPP_FILES = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.cpp))
CPP_OBJ = $(patsubst %.cpp,%.o,$(CPP_FILES))

C_FILES = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.c))
C_OBJ = $(patsubst %.c,%.o,$(C_FILES))

F90_FILES = $(wildcard *.f90)
F90_OBJ = $(patsubst %.f90,%.o,$(F90_FILES))

HEADERS = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.h))

all: FLAGS += -O3
all: $(F90_OBJ) $(C_OBJ) $(CPP_OBJ) $(HEADERS)
	g++ -std=c++17 -o Tomator1D $(F90_OBJ) $(C_OBJ) $(CPP_OBJ) $(FLAGS)

debug: CFLAGS += -g -Ddebug -ggdb3
debug: $(F90_OBJ) $(C_OBJ) $(CPP_OBJ) $(HEADERS)
	g++ -std=c++17 -o Tomator1D $(F90_OBJ) $(C_OBJ) $(CPP_OBJ) $(FLAGS) $(CFLAGS)

%.o: %.cpp $(HEADERS)
	g++ -std=c++17 -c -o $@ $< $(FLAGS) $(CFLAGS)

%.o: %.c $(HEADERS)
	gcc -c -o $@ $< $(FLAGS) $(CFLAGS)

%.o: %.f90
	gfortran -c -o $@ $< -I./Fortran -J./Fortran -fopenmp

clean:
	rm -f Tomator1D
	rm -f $(CPP_OBJ) $(C_OBJ) $(F90_OBJ)
	@echo "Cleanup complete!"

flush:
	rm -f Tomator1D
	rm -f $(CPP_OBJ) $(C_OBJ) $(F90_OBJ)
	rm -f *.mod
	rm -f *.log
