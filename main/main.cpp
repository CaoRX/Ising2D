#ifndef MAIN_CPP
#define MAIN_CPP
#include <lattice/lattice.hpp>
#include <config/Parameters.hpp>
#include <iostream>
ISING_LATTICE IL(10, 10);
int main()
{
	std::cout << "Hello world!" << std::endl;
	std::cout << Config::J_NN << std::endl;
	return 0;
}

#endif