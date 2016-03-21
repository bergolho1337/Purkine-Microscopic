CC = g++
IDIR = ../include
CPPFLAGS = -std=c++11 -O2 -lm -Wall -I$(IDIR)
CPPSOURCES = $(wildcard *.cpp)
ODIR = obj
_OBJ = $(CPPSOURCES:.cpp=.o)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: microscopic

microscopic: $(OBJ)
	$(CC) -o $@ $^

$(ODIR)/%.o: %.cpp
	$(CC) -c $< $(CPPFLAGS) -o $@

clean:
	-rm -r ../Runs
	-rm -rf $(ODIR)/*.o microscopic
	
remade:
	$(MAKE) clean
	$(MAKE)
