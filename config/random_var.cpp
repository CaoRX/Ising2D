#ifndef RANDOM_VAR_CPP
#define RANDOM_VAR_CPP

#include <random>
#include <ctime>
using namespace std;

unsigned random_seed = (unsigned)time(NULL);
//unsigned random_seed = 812;
default_random_engine e(random_seed);
uniform_int_distribution<int> random_4(0,3);
uniform_real_distribution<double> random_Double(0,1.0);

double random_double()
{
	return random_Double(e);
}

#endif