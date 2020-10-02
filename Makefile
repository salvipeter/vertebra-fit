all: test-fit

CXXFLAGS=-std=c++17 -Wall -pedantic -I/usr/include/eigen3

test-fit: test-fit.o ccfit.o ccurve.o
	$(CXX) -o $@ $^
