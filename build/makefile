CXX = g++
CFLAGS = -c -Wall -pedantic
LIB_PATH = ../lib
INIT_PATH = ../init
SRC_PATH = ../src
VPATH = ./:$(LIB_PATH):$(INIT_PATH):$(SRC_PATH)
INCLUDE_DIRS = -I. -I$(LIB_PATH) -I$(INIT_PATH)
LDFLAGS = -lm

MAIN_OBJ = main.o integrators.o utils.o init.o
CONV_OBJ = convergence_study.o integrators.o utils.o

main: $(MAIN_OBJ)
	$(CXX) $(MAIN_OBJ) $(LDFLAGS) -o $@

convergence: $(CONV_OBJ)
	$(CXX) $(CONV_OBJ) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $<

clean:
	rm -f *.o
	@echo make clean: done
