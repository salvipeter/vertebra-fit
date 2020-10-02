all: test-fit

CXXFLAGS=-std=c++17 -Wall -pedantic

test-fit: test-fit.o ccfit.o ccurve.o
	$(CXX) -o $@ $^
