CC          := g++

CXXFLAGS        := -O3 -Wall -std=c++11

../bin/main: main.cpp
	$(CC) $(CXXFLAGS) -pthread -o $@ $? -fopenmp
clean:
	rm -f ../bin/main
