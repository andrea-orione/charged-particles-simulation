CXX = g++
CFLAGS = -c -Wall -pedantic
LIB_PATH = ../lib
VPATH = ./:$(LIB_PATH)
INCLUDE_DIRS = -I. -I$(LIB_PATH)
LDFLAGS = -lm

CONV_OBJ = convergence_study.o integrators.o

convergence: $(CONV_OBJ)
	$(CXX) $(CONV_OBJ) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $<

clean:
	rm -f *.o
	@echo make clean: done
