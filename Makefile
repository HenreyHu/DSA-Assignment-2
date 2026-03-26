CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

all: simplify

simplify: simplify.cpp structures.h
	$(CXX) $(CXXFLAGS) -o simplify simplify.cpp

clean:
	rm -f simplify

.PHONY: all clean
