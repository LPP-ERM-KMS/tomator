# Specify directories
SRC_DIRS = . Funcs Vars

# Using wildcard and patsubst functions to find source files and create object file names
CPP_FILES = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.cpp))
CPP_OBJ = $(patsubst %.cpp,%.o,$(CPP_FILES))

C_FILES = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.c))
C_OBJ = $(patsubst %.c,%.o,$(C_FILES))

HEADERS = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.h))

all: FLAGS = -fopenmp -O3 -lpthread -Wall 

all: $(C_OBJ) $(CPP_OBJ) $(HEADERS)
	g++ -std=c++23 -o Tomator1D $(C_OBJ) $(CPP_OBJ) $(FLAGS)

debug: CFLAGS += -g  -ggdb3 -Ddebug

debug: $(C_OBJ) $(CPP_OBJ) $(HEADERS)
	g++ -std=c++23 -o Tomator1D $(C_OBJ) $(CPP_OBJ) $(FLAGS) $(CFLAGS)

%.o: %.cpp $(HEADERS)
	g++ -std=c++23 -c -o $@ $< $(FLAGS) $(CFLAGS)

%.o: %.c $(HEADERS)
	gcc -c -o $@ $< $(FLAGS) $(CFLAGS)

clean:
	rm -f Tomator1D
	rm -f $(CPP_OBJ) $(C_OBJ) 
	@echo "Cleanup complete!"

flush:
	rm -f Tomator1D
	rm -f $(CPP_OBJ) $(C_OBJ)
	rm -f *.mod
	rm -f *.log
