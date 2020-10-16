all: simulation.cpp compartment.cpp main.cpp
	gcc -shared -o simulate.so -Wall -fPIC simulation.cpp compartment.cpp main.cpp
